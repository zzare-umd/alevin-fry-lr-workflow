# alevin-fry-lr: 10x 3' long reads (ONT or PacBio) -> per-cell transcript (or gene) counts.
#
#   snakemake --cores 32          (from the repository directory, after filling in config/config.yaml)
#
# Steps per library: BLAZE (cell barcodes and UMIs) -> minimap2 against the transcriptome -> CB/UMI into
# BAM tags -> alevin-fry convert / generate-permit-list / collate / quant. See README.md.
from os.path import join

# config/config.yaml next to this file is read unless another file is given with --configfile (the
# given file replaces it, it is not merged)
if not workflow.overwrite_configfiles:
    configfile: join(workflow.basedir, "config", "config.yaml")

SAMPLES = config["samples"]
RES = config["results_dir"]
# the alevin-fry binary; when the file does not exist, rule install_alevin_fry builds it there
AF = config["alevin_fry"]
AF_GIT = config.get("alevin_fry_git", "https://github.com/COMBINE-lab/alevin-fry")
AF_BRANCH = config.get("alevin_fry_branch", "long-read")
THREADS = int(config.get("threads", 16))
SCRIPT_DIR = join(workflow.basedir, "scripts")
PRESET = {"ONT": "map-ont", "PacBio": "map-hifi"}


def platform_of(sample):
    p = SAMPLES[sample]["platform"]
    if p not in PRESET:
        raise ValueError(f"sample {sample}: platform must be ONT or PacBio, not {p!r}")
    return p


# the default target: the count matrix of every library
rule all:
    input:
        expand(join(RES, "{sample}", "quant", "alevin", "quants_mat.mtx"), sample=SAMPLES),


# stage targets, to run or re-run the pipeline up to one step for every library:
#   snakemake --configfile ... --cores 32 all_blaze
rule all_blaze:
    input:
        expand(join(RES, "{sample}", "blaze", "matched_reads.fastq.gz"), sample=SAMPLES),


rule all_align:
    input:
        expand(join(RES, "{sample}", "align", "transcripts.tagged.bam"), sample=SAMPLES),


rule all_rad:
    input:
        expand(join(RES, "{sample}", "rad", "map.collated.rad"), sample=SAMPLES),


rule all_quant:
    input:
        rules.all.input,


# the long-read alevin-fry, built from its branch with cargo (rust, cmake and a C compiler are in the
# conda environment); runs only when the configured binary does not exist
rule install_alevin_fry:
    output:
        AF,
    params:
        git=AF_GIT,
        branch=AF_BRANCH,
    log:
        AF + ".install.log",
    shell:
        """
        tmp=$(mktemp -d)
        cargo install --locked --git {params.git} --branch {params.branch} --root $tmp alevin-fry > {log} 2>&1
        install -D $tmp/bin/alevin-fry {output}
        rm -rf $tmp
        {output} --version >> {log}
        """


# transcript -> gene when a tx2gene file is configured, else transcript -> transcript (per-transcript counts)
rule transcript_map:
    input:
        fa=config["transcriptome"],
    output:
        join(RES, "reference", "txp_map.tsv"),
    params:
        t2g=config.get("tx2gene") or "",
        script=join(SCRIPT_DIR, "make_t2t_map.py"),
    shell:
        """
        if [ -n "{params.t2g}" ]; then cut -f1,2 "{params.t2g}" > {output}
        else python {params.script} -i {input.fa} -o {output}; fi
        """


rule minimap2_index:
    input:
        config["transcriptome"],
    output:
        join(RES, "reference", "transcripts.{platform}.mmi"),
    params:
        preset=lambda w: PRESET[w.platform],
    threads: THREADS
    log:
        join(RES, "reference", "minimap2_index.{platform}.log"),
    shell:
        "minimap2 -t {threads} -x {params.preset} -d {output} {input} > {log} 2>&1"


# cell barcode and UMI of every read; matched reads come out in the transcript orientation with
# CB_UMI#readid names, which tag_bam reads
rule blaze:
    input:
        fastq=lambda w: SAMPLES[w.sample]["fastq"],
        whitelist=config["whitelist"],
    output:
        fastq=join(RES, "{sample}", "blaze", "matched_reads.fastq.gz"),
        whitelist=join(RES, "{sample}", "blaze", "whitelist.csv"),
    params:
        prefix=join(RES, "{sample}", "blaze") + "/",
        cells=lambda w: SAMPLES[w.sample]["expect_cells"],
    threads: THREADS
    log:
        join(RES, "{sample}", "logs", "blaze.log"),
    shell:
        """
        mkdir -p {params.prefix}
        blaze --full-bc-whitelist {input.whitelist} --expect-cells {params.cells} --output-prefix {params.prefix} \
            --output-fastq matched_reads.fastq.gz --threads {threads} --minimal_stdout {input.fastq} > {log} 2>&1
        """


# all alignments within the score range are kept (-N 100); unmapped and supplementary records are dropped
rule align:
    input:
        fastq=rules.blaze.output.fastq,
        index=lambda w: join(RES, "reference", f"transcripts.{platform_of(w.sample)}.mmi"),
    output:
        bam=join(RES, "{sample}", "align", "transcripts.bam"),
    threads: THREADS
    log:
        join(RES, "{sample}", "logs", "minimap2.log"),
    shell:
        "minimap2 -t {threads} -a -N 100 {input.index} {input.fastq} 2> {log} "
        "| samtools view -@ {threads} -h -F 2052 -b -o {output.bam} -"


rule tag_bam:
    input:
        bam=rules.align.output.bam,
    output:
        bam=join(RES, "{sample}", "align", "transcripts.tagged.bam"),
    params:
        script=join(SCRIPT_DIR, "tagged_bam_CR_UR.py"),
    log:
        join(RES, "{sample}", "logs", "tag_bam.log"),
    shell:
        "python {params.script} --input {input.bam} --output {output.bam} > {log} 2>&1"


rule convert:
    input:
        bam=rules.tag_bam.output.bam,
        af=ancient(AF),
    output:
        rad=join(RES, "{sample}", "rad", "map.rad"),
    params:
        af=AF,
        th=config.get("score_threshold", 0.95),
    threads: THREADS
    log:
        join(RES, "{sample}", "logs", "convert.log"),
    shell:
        "{params.af} convert --threads {threads} --filter_best --score-threshold {params.th} "
        "--bam {input.bam} --output {output.rad} > {log} 2>&1"


rule permit_list:
    input:
        rad=rules.convert.output.rad,
        whitelist=config["whitelist"],
        af=ancient(AF),
    output:
        json=join(RES, "{sample}", "rad", "generate_permit_list.json"),
        unmapped=join(RES, "{sample}", "rad", "unmapped_bc_count.bin"),
    params:
        af=AF,
        raddir=join(RES, "{sample}", "rad"),
        ori=config.get("expected_ori", "fw"),
    threads: THREADS
    log:
        join(RES, "{sample}", "logs", "generate_permit_list.log"),
    shell:
        """
        {params.af} generate-permit-list --threads {threads} --input {params.raddir} --expected-ori {params.ori} \
            --output-dir {params.raddir} --unfiltered-pl {input.whitelist} > {log} 2>&1
        touch {output.unmapped}
        """


rule collate:
    input:
        rules.permit_list.output.json,
        rules.permit_list.output.unmapped,
        ancient(AF),
    output:
        join(RES, "{sample}", "rad", "map.collated.rad"),
    params:
        af=AF,
        raddir=join(RES, "{sample}", "rad"),
    threads: THREADS
    log:
        join(RES, "{sample}", "logs", "collate.log"),
    shell:
        "{params.af} collate -i {params.raddir} -r {params.raddir} -t {threads} > {log} 2>&1"


rule quant:
    input:
        collated=rules.collate.output,
        map=rules.transcript_map.output,
        af=ancient(AF),
    output:
        mtx=join(RES, "{sample}", "quant", "alevin", "quants_mat.mtx"),
    params:
        af=AF,
        raddir=join(RES, "{sample}", "rad"),
        outdir=join(RES, "{sample}", "quant"),
        prior=config["quant"]["cell_prior"],
        power=config["quant"]["positional_power"],
    threads: THREADS
    log:
        join(RES, "{sample}", "logs", "quant.log"),
    shell:
        """
        mkdir -p {params.outdir}
        {params.af} quant -i {params.raddir} -m {input.map} -t {threads} -o {params.outdir} -r parsimony-em \
            --lambda-size=0.0 --tau-delta=0.0 \
            --dedup-model bf --ed-distance-model general --ed-upper-bound 2 --dedup-edges bayes --ed-far-window 10 \
            --model-truncation --truncation-interval-mode all --truncation-end-kernels absolute \
            --truncation-cell-prior {params.prior} --truncation-positional-power {params.power} > {log} 2>&1
        """
