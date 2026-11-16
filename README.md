# alevin-fry-lr

Isoform-level quantification of 10x single-cell long-read libraries (Oxford Nanopore or PacBio) with
`alevin-fry`. This repository is the current test setup: a user guide and the Snakemake workflow we run,
which also builds the long-read version of `alevin-fry` for you.

The long-read code lives on the `long-read` branch of [COMBINE-lab/alevin-fry](https://github.com/COMBINE-lab/alevin-fry/tree/long-read).
It is research software under active development: it works end to end on our data (10x 3' v3 libraries
on ONT and PacBio, simulated and real), but the installation and the first three steps of the pipeline
are not yet packaged the way the rest of `alevin-fry` is. The section [What we are streamlining](#what-we-are-streamlining)
says what will change.

## What the pipeline does

| Step | Tool | What it does | Our ONT test: 72 M reads, 5,000 cells, 64 threads |
|---|---|---|---|
| 1 | [BLAZE](https://github.com/shimlab/BLAZE) | Finds the cell barcode and UMI in each read, keeps the reads of real cells, writes them in transcript orientation with `CB_UMI#read` names | 60 min, 1.3 GB |
| 2 | minimap2 | Aligns the reads to the transcriptome, keeping up to 100 alignments per read | 16 min, 6.5 GB |
| 3 | `tagged_bam_CR_UR.py` | Moves the barcode and UMI from the read name into BAM tags (`CB`/`CR`, `UB`/`UR`) | 42 min (single-threaded), 0.1 GB |
| 4 | `alevin-fry convert` | Turns the BAM into a RAD file, keeping each read's alignments within 95% of its best alignment score | 2 min, 0.1 GB |
| 5 | `alevin-fry generate-permit-list` | Decides which barcodes are cells (the 10x whitelist, at least 10 reads) and corrects the rest | 2 s |
| 6 | `alevin-fry collate` | Groups the records by cell | 3 s, 1.9 GB |
| 7 | `alevin-fry quant` | Per cell: collapses reads into molecules (UMI deduplication that tolerates indels and uses read positions), then assigns each molecule to isoforms with a positional model of 3'-anchored, truncated long reads | 20 s wall, 15 CPU-min, 5.3 GB |

The result is a cells × transcripts count matrix (or cells × genes with a transcript-to-gene file).

## Requirements

- Linux, x86-64. Memory is set by minimap2 on your transcriptome (about 7 GB for a human
  transcriptome) and by `quant`, which uses about 0.25 GB per thread plus the matrix.
- One conda environment, `workflow/envs/alevin-fry-lr.yaml`: [BLAZE](https://github.com/shimlab/BLAZE) 2.5,
  minimap2 ≥ 2.28, samtools, Python with pysam, Snakemake ≥ 8, and the Rust toolchain, `cmake` and C
  compiler that build `alevin-fry` (cmake builds the compression library it links).

## Installation

```bash
git clone <this repository> alevin-fry-lr && cd alevin-fry-lr
conda env create -f workflow/envs/alevin-fry-lr.yaml
conda activate alevin-fry-lr
```

That is all. The first `snakemake` run builds `alevin-fry` from its `long-read` branch into `bin/alevin-fry`
(the `alevin_fry` path in the config), a few minutes once; it fetches the `libradicl` release the branch
depends on by itself. If you already have the binary, set `alevin_fry` to its path and nothing is built.
To build it by hand instead:

```bash
cargo install --locked --git https://github.com/COMBINE-lab/alevin-fry --branch long-read --root . alevin-fry
```

The commit this guide was tested with is d63db0c.

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

Copy `workflow/config/config.yaml`, fill in the paths, and run from the repository directory:

```bash
snakemake -s workflow/Snakefile --configfile workflow/config/my_config.yaml --cores 32 --dry-run   # the plan
snakemake -s workflow/Snakefile --configfile workflow/config/my_config.yaml --cores 32
```

Without `--configfile`, the template `workflow/config/config.yaml` next to the Snakefile is read. The
default target `all` is the count matrix of every library; `all_blaze`, `all_align`, `all_rad` and
`all_quant` run the pipeline up to that stage for every library, for example
`snakemake -s workflow/Snakefile --configfile ... --cores 32 all_align` to stop after the alignments.

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

## Running the steps by hand

Every step is one rule in `workflow/Snakefile` with its exact command line, so the file doubles as the
list of commands to run by hand; `snakemake --printshellcmds --dry-run` prints them with your paths
filled in.

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

Two things to know about the rows:

- The matrix holds every whitelist barcode with at least 10 reads, not only real cells. Filter cells
  afterwards, for example to the barcodes in `results/<sample>/blaze/whitelist.csv` (BLAZE's call) or by
  a knee on total counts.
- `truncation_model.tsv` has one row per cell with the fitted positional model (read-length laws, the
  fractions of 3'-anchored, 5'-anchored, full-length and uniformly placed reads, and convergence flags).
  `FitStatus` tells whether a cell had enough molecules to fit anything.

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

An alternative model pools isoform abundances across cells of the same type instead of fitting each cell
alone (`quant --truncation-preset recommended --truncation-discover-groups <tx2gene>`); it is more
accurate where many cells share one isoform usage and we are still comparing the two on real data. The
per-cell setting above is the one we ship because it makes no assumption about which cells are alike.

## Troubleshooting

- **BLAZE finds almost no cells.** Check the whitelist matches the chemistry (3' v3: `3M-february-2018.txt`;
  other kits ship other lists). With some basecalling or library setups the barcodes are on the other
  strand: reverse-complement the whitelist and run again.
- **Other 10x chemistries.** BLAZE has `--kit-version` and `--umi-len`; the rest of the pipeline only needs
  the `CB_UMI#` read names and the matching whitelist. We have run 3' v3 (16-bp barcode, 12-bp UMI).
- **`collate` complains about `unmapped_bc_count.bin`.** The workflow creates that file after the permit
  list; if you run the steps by hand, `touch` it in the RAD directory.
- **`convert` finds no barcodes.** It reads the `CB`/`UB` tags; the tagging step must have run on reads
  named by BLAZE. Reads without the name pattern pass through untagged and are dropped.
- **Building `alevin-fry` fails.** The log is next to the binary's path (`bin/alevin-fry.install.log`).
  Rust, cmake and the compiler come from the conda environment; outside it, a Rust toolchain ≥ 1.85
  and cmake ≥ 3.5 are needed (`CMAKE=/path/to/cmake` if cmake is not in `PATH`).
- **Memory.** `quant` memory grows with threads (about 0.25 GB each) and with the largest cell; use
  fewer threads if it is tight. minimap2 needs the index in memory.

## What we are streamlining

These are the parts of the current setup that we know are rough, and what is being done about them:

- **Installation.** Building `alevin-fry` from a branch is temporary; the long-read code will be merged
  and released so that a bioconda package suffices.
- **Three external steps.** BLAZE, minimap2 and the tagging script are run separately today. The tagging
  step in particular is a single-threaded Python loop (40 minutes on 72 M reads) that exists only to move
  the barcode from the read name into a tag; `convert` will read BLAZE's names directly, and a single
  command will drive demultiplexing and alignment.
- **Defaults.** `quant` needs a long list of flags today. They are the settings we found best and will
  become the defaults of a long-read mode, so that `quant -r parsimony-em --long-read` is enough.
- **Model choice.** Two isoform models are implemented (the per-cell one above and a pooled one); we are
  finishing the comparison on real data with cell types of different isoform usage before fixing the
  default.
- **Documentation.** The flags are documented in `docs/source/quant.rst` on the branch; a long-read page
  on the `alevin-fry` documentation site will replace this README.

## How it works, briefly

`alevin-fry` processes RAD files, which hold each read's barcode, UMI and candidate alignments; the
long-read records (in `libradicl` 0.13) carry the alignment intervals the positional model needs and no
read names, which keeps the files small. `quant` works one cell at a time. It first builds the UMI graph of the
cell and resolves it into molecules (He et al. 2022 describe the parsimony resolution; the long-read
branch adds a general edit distance and the position-aware edge model). It then fits a mixture model of
read placement (after the one in [oarfish](https://github.com/COMBINE-lab/oarfish)): each read is 3'-anchored,
5'-anchored, full-length or uniformly placed, with learned read-length laws and end kernels. The
assignment of a molecule to an isoform is the product of its alignment score, the positional likelihood
of its read on that isoform (raised to `positional_power`) and the isoform's abundance in the cell,
estimated by EM with the `cell_prior` pseudo-count. Counts are the sums of these posteriors.

If you use this, please cite the `alevin-fry` paper (He et al., *Nature Methods* 2022) and `oarfish`
(Zare Jousheghani and Patro, 2024) for the positional model; a preprint describing the long-read method
is in preparation.

## Questions

Open an issue in this repository, or reply to the thread that pointed you here. Feedback on real data is
exactly what we are after at this stage: what failed, what looked odd, how long it took.
