# alevin-fry-lr

Isoform-level quantification of 10x single-cell long-read libraries (Oxford Nanopore or PacBio) with
`alevin-fry`. This repository is the current test setup: a user guide and the Snakemake workflow we run,
which also builds the long-read version of `alevin-fry` for you.

The long-read code lives on the `long-read` branch of [COMBINE-lab/alevin-fry](https://github.com/COMBINE-lab/alevin-fry/tree/long-read).
It is research software under active development: it works end to end on our data (10x 3' v3 libraries
on ONT and PacBio, simulated and real), but the installation and the first three steps of the pipeline
are not yet packaged the way the rest of `alevin-fry` is.

## What the pipeline does

| Step | Tool | What it does |
|---|---|---|
| 1 | [BLAZE](https://github.com/shimlab/BLAZE) | Finds the cell barcode and UMI in each read, keeps the reads of real cells, writes them in transcript orientation with `CB_UMI#read` names |
| 2 | minimap2 | Aligns the reads to the transcriptome, keeping up to 100 alignments per read |
| 3 | `tagged_bam_CR_UR.py` | Moves the barcode and UMI from the read name into BAM tags (`CB`/`CR`, `UB`/`UR`) |
| 4 | `alevin-fry convert` | Turns the BAM into a RAD file, keeping each read's alignments within 95% of its best alignment score |
| 5 | `alevin-fry generate-permit-list` | Decides which barcodes are cells (the 10x whitelist, at least 10 reads) and corrects the rest |
| 6 | `alevin-fry collate` | Groups the records by cell |
| 7 | `alevin-fry quant` | Per cell: collapses reads into molecules (UMI deduplication that tolerates indels and uses read positions), then assigns each molecule to isoforms with a positional model of 3'-anchored, truncated long reads |

## Requirements

- Linux, x86-64. Memory is set by minimap2 on your transcriptome (about 7 GB for a human
  transcriptome) and by `quant`, which uses about 0.25 GB per thread plus the matrix.
- One conda environment, `envs/alevin-fry-lr.yaml`: [BLAZE](https://github.com/shimlab/BLAZE) 2.5,
  minimap2 ≥ 2.28, samtools, Python with pysam, Snakemake ≥ 8, and the Rust toolchain, `cmake` and C
  compiler that build `alevin-fry` (cmake builds the compression library it links).

## Installation

```bash
git clone git@github.com:zzare-umd/alevin-fry-lr-workflow.git && cd alevin-fry-lr-workflow
conda env create -f envs/alevin-fry-lr.yaml
conda activate alevin-fry-lr
```

That is all. The first `snakemake` run builds `alevin-fry` from its `long-read` branch into `bin/alevin-fry`
(the `alevin_fry` path in the config), a few minutes once; it fetches the `libradicl` release the branch
depends on by itself. If you already have the binary, set `alevin_fry` to its path and nothing is built.
To build it by hand instead:

```bash
cargo install --locked --git https://github.com/COMBINE-lab/alevin-fry --branch long-read --root . alevin-fry
```

## What you need

1. The reads of the library as `fastq.gz`, as they come from basecalling (no demultiplexing; BLAZE does
   that). ONT or PacBio; both are handled, only the minimap2 preset differs.
2. The 10x barcode whitelist of your chemistry. For 3' v3 it is `3M-february-2018.txt`, shipped with
   Cell Ranger under `lib/python/cellranger/barcodes/` (gunzip it). BLAZE and `generate-permit-list`
   both use it.
3. A transcriptome FASTA (for example GENCODE `transcripts.fa`). Transcript names are the first word of
   each FASTA header. Optionally a `transcript<TAB>gene` file if you want counts per gene instead of
   per transcript.
4. The expected number of cells, roughly (BLAZE uses it to find the cell barcodes).

## Running the workflow

Fill in the paths in `config/config.yaml` and run from the repository directory:

```bash
snakemake --cores 32 --dry-run   # shows the plan
snakemake --cores 32
```

`--cores` is the total number of CPU cores Snakemake may use. The default target `all` is the count
matrix of every library; `all_blaze`, `all_align`, `all_rad` and `all_quant` run the pipeline up to
that stage for every library, for example `snakemake --cores 32 all_align` to stop after the
alignments. To keep several configurations, copy `config/config.yaml` and pass the copy with
`--configfile path/to/copy.yaml`.

The configuration has one entry per library (`fastq`, `platform`, `expect_cells`), the whitelist, the
transcriptome, the optional transcript-to-gene file, the path of the binary, the output directory and
the thread count. Everything is written under `results_dir/<sample>/`:

```
results/<sample>/blaze/      matched_reads.fastq.gz, whitelist.csv (the cells BLAZE called), summary.txt
results/<sample>/align/      transcripts.bam, transcripts.tagged.bam
results/<sample>/rad/        map.rad, generate_permit_list.json, map.collated.rad
results/<sample>/quant/      alevin/quants_mat.mtx, quants_mat_rows.txt, quants_mat_cols.txt, quant.json,
                             featureDump.txt, truncation_model.tsv
results/<sample>/logs/       one log per step
results/reference/           txp_map.tsv, transcripts.<platform>.mmi
bin/alevin-fry               the binary, when the workflow built it
```

## The output

`quants_mat.mtx` is a Matrix Market file, cells in rows (`quants_mat_rows.txt`, barcodes) and transcripts
or genes in columns (`quants_mat_cols.txt`). Entries are molecules (UMIs); they can be fractional,
because a molecule whose reads fit several isoforms is split among them by its posterior. To load it in
Python:

```python
import scipy.io, pandas as pd, anndata as ad
d = "results/pbmc_ont/quant/alevin"
X = scipy.io.mmread(f"{d}/quants_mat.mtx").tocsr()
adata = ad.AnnData(X,
    obs=pd.DataFrame(index=open(f"{d}/quants_mat_rows.txt").read().split()),
    var=pd.DataFrame(index=open(f"{d}/quants_mat_cols.txt").read().split()))
```

## Parameters that matter

Defaults are the settings we run; all of them live in the config file.

- `score_threshold` (`convert --score-threshold`, 0.95): a read keeps every alignment whose score is at
  least this fraction of its best one. Those alignments are the molecule's candidate isoforms. 1.0 keeps
  only the best-scoring alignments and is slightly stricter; both work.
- `expected_ori` (`generate-permit-list --expected-ori fw`): BLAZE writes reads in the transcript
  orientation, so alignments are expected forward. Leave it unless you skip BLAZE.
- The UMI deduplication (`quant`): `--ed-distance-model general --ed-upper-bound 2` joins UMIs within an
  edit distance of 2 (indels included, which matters for nanopore UMIs), `--dedup-model bf` resolves each
  cell's UMI graph with a likelihood model, and `--dedup-edges bayes --ed-far-window 10` decides whether
  two near-identical UMIs are one molecule from their distance and whether their reads lie at the same
  place on the transcript, with the error rates estimated on each cell. Nothing to set here.
- The isoform model (`quant --model-truncation ...`): reads of 10x 3' libraries start at the poly-A end
  and usually stop before the 5' end, so a read's position on each candidate isoform says something about
  which isoform it came from. The model learns, per cell, how reads are placed and uses it together with
  the isoform abundances. Two knobs: `positional_power` (0.5) tempers how strongly positions weigh in
  against abundances, and `cell_prior` (0.25) is a small pseudo-count toward equal isoform shares that
  keeps a cell with two molecules from declaring an isoform absent. We do not recommend changing them
  for a first run.
- `threads`: every step scales; BLAZE and minimap2 take most of the time.
  


## Questions

Open an issue in this repository, or reply to the thread that pointed you here. Feedback on real data is
exactly what we are after at this stage: what failed, what looked odd, how long it took.
