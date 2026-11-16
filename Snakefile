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


def as_list(v):
    return v if isinstance(v, list) else [v]


SAMPLES = config["samples"]
RES = config["results_dir"]
THREADS = int(config.get("threads", 16))
SCRIPT_DIR = join(workflow.basedir, "scripts")
PRESET = {"ONT": "map-ont", "PacBio": "map-hifi"}
# one value or a list: every combination is run, each in its own directory
TH = [str(x) for x in as_list(config.get("score_threshold", 0.95))]
POWER = [str(x) for x in as_list(config["quant"]["positional_power"])]
CELL_PRIOR = config["quant"]["cell_prior"]
DEDUP_EDGES = config["quant"].get("dedup_edges", "bayes")   # bayes (calibrated per cell) or rule (edit distance only)
ORI = config.get("expected_ori", "fw")
# the alevin-fry binary; when the file does not exist, rule install_alevin_fry builds it there
AF = config["alevin_fry"]
AF_GIT = config.get("alevin_fry_git", "https://github.com/COMBINE-lab/alevin-fry")
AF_BRANCH = config.get("alevin_fry_branch", "long-read")


wildcard_constraints:
    th=r"[0-9.]+",
    power=r"[0-9.]+",


def platform_of(sample):
    p = SAMPLES[sample]["platform"]
    if p not in PRESET:
        raise ValueError(f"sample {sample}: platform must be ONT or PacBio, not {p!r}")
    return p


def whitelist_of(sample):
    # a library's own whitelist (another chemistry) overrides the common one
    return SAMPLES[sample].get("whitelist", config["whitelist"])


# the default target: the count matrix of every library, threshold and power
rule all:
    input:
        expand(join(RES, "{sample}", "quant_th{th}_p{power}", "alevin", "quants_mat.mtx"), sample=SAMPLES, th=TH, power=POWER),


# stage targets, to run or re-run the pipeline up to one step for every library:
#   snakemake --cores 32 all_blaze
rule all_blaze:
    input:
        expand(join(RES, "{sample}", "blaze", "matched_reads.fastq.gz"), sample=SAMPLES),


rule all_align:
    input:
        expand(join(RES, "{sample}", "align", "transcripts.tagged.bam"), sample=SAMPLES),


rule all_rad:
    input:
        expand(join(RES, "{sample}", "rad_th{th}", "map.collated.rad"), sample=SAMPLES, th=TH),


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
    benchmark:
        join(RES, "reference", "benchmarks", "minimap2_index.{platform}.tsv")
    shell:
        "minimap2 -t {threads} -x {params.preset} -d {output} {input} > {log} 2>&1"


# cell barcode and UMI of every read; matched reads come out in the transcript orientation with
# CB_UMI#readid names, which tag_bam reads
rule blaze:
    input:
        fastq=lambda w: SAMPLES[w.sample]["fastq"],
        whitelist=lambda w: whitelist_of(w.sample),
    output:
        fastq=join(RES, "{sample}", "blaze", "matched_reads.fastq.gz"),
        whitelist=join(RES, "{sample}", "blaze", "whitelist.csv"),
    params:
        prefix=join(RES, "{sample}", "blaze") + "/",
        cells=lambda w: SAMPLES[w.sample]["expect_cells"],
    threads: THREADS
    log:
        join(RES, "{sample}", "logs", "blaze.log"),
    benchmark:
        join(RES, "{sample}", "benchmarks", "blaze.tsv")
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
    benchmark:
        join(RES, "{sample}", "benchmarks", "minimap2.tsv")
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
    benchmark:
        join(RES, "{sample}", "benchmarks", "tag_bam.tsv")
    shell:
        "python {params.script} --input {input.bam} --output {output.bam} > {log} 2>&1"


# each read keeps the alignments within the threshold's fraction of its best alignment score
rule convert:
    input:
        bam=rules.tag_bam.output.bam,
        af=ancient(AF),
    output:
        rad=join(RES, "{sample}", "rad_th{th}", "map.rad"),
    params:
        af=AF,
    threads: THREADS
    log:
        join(RES, "{sample}", "logs", "convert_th{th}.log"),
    benchmark:
        join(RES, "{sample}", "benchmarks", "convert_th{th}.tsv")
    shell:
        "{params.af} convert --threads {threads} --filter_best --score-threshold {wildcards.th} "
        "--bam {input.bam} --output {output.rad} > {log} 2>&1"


rule permit_list:
    input:
        rad=rules.convert.output.rad,
        whitelist=lambda w: whitelist_of(w.sample),
        af=ancient(AF),
    output:
        json=join(RES, "{sample}", "rad_th{th}", "generate_permit_list.json"),
        unmapped=join(RES, "{sample}", "rad_th{th}", "unmapped_bc_count.bin"),
    params:
        af=AF,
        raddir=join(RES, "{sample}", "rad_th{th}"),
        ori=ORI,
    threads: THREADS
    log:
        join(RES, "{sample}", "logs", "generate_permit_list_th{th}.log"),
    benchmark:
        join(RES, "{sample}", "benchmarks", "generate_permit_list_th{th}.tsv")
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
        join(RES, "{sample}", "rad_th{th}", "map.collated.rad"),
    params:
        af=AF,
        raddir=join(RES, "{sample}", "rad_th{th}"),
    threads: THREADS
    log:
        join(RES, "{sample}", "logs", "collate_th{th}.log"),
    benchmark:
        join(RES, "{sample}", "benchmarks", "collate_th{th}.tsv")
    shell:
        "{params.af} collate -i {params.raddir} -r {params.raddir} -t {threads} > {log} 2>&1"


rule quant:
    input:
        collated=rules.collate.output,
        map=rules.transcript_map.output,
        af=ancient(AF),
    output:
        mtx=join(RES, "{sample}", "quant_th{th}_p{power}", "alevin", "quants_mat.mtx"),
    params:
        af=AF,
        raddir=join(RES, "{sample}", "rad_th{th}"),
        outdir=join(RES, "{sample}", "quant_th{th}_p{power}"),
        prior=CELL_PRIOR,
        edges=DEDUP_EDGES,
    threads: THREADS
    log:
        join(RES, "{sample}", "logs", "quant_th{th}_p{power}.log"),
    benchmark:
        join(RES, "{sample}", "benchmarks", "quant_th{th}_p{power}.tsv")
    shell:
        """
        mkdir -p {params.outdir}
        {params.af} quant -i {params.raddir} -m {input.map} -t {threads} -o {params.outdir} -r parsimony-em \
            --lambda-size=0.0 --tau-delta=0.0 \
            --dedup-model bf --ed-distance-model general --ed-upper-bound 2 --dedup-edges {params.edges} --ed-far-window 10 \
            --model-truncation --truncation-interval-mode all --truncation-end-kernels absolute \
            --truncation-cell-prior {params.prior} --truncation-positional-power {wildcards.power} > {log} 2>&1
        """
