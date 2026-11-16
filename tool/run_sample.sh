#!/usr/bin/env bash
# alevin-fry-lr on one library without Snakemake: the same steps and settings as workflow/Snakefile.
#
#   tool/run_sample.sh -f reads.fastq.gz -w 3M-february-2018.txt -x transcripts.fa -o out_dir \
#       [-p ONT|PacBio] [-c expected_cells] [-t threads] [-m tx2gene.tsv] [-b path/to/alevin-fry] [-s score_threshold]
#
# Needs blaze, minimap2, samtools, python with pysam in PATH (workflow/envs/alevin-fry-lr.yaml) and the
# alevin-fry-lr binary (tool/install.sh). Output: out_dir/quant/alevin/quants_mat.mtx (+ rows/cols).
set -euo pipefail
HERE=$(cd "$(dirname "$0")" && pwd)
PLATFORM=ONT; CELLS=5000; THREADS=16; T2G=""; AF="$HERE/bin/alevin-fry"; TH=0.95
PRIOR=0.25; POWER=0.5; ORI=fw
while getopts "f:w:x:o:p:c:t:m:b:s:" opt; do
  case $opt in
    f) FASTQ=$OPTARG;; w) WL=$OPTARG;; x) TX=$OPTARG;; o) OUT=$OPTARG;; p) PLATFORM=$OPTARG;; c) CELLS=$OPTARG;;
    t) THREADS=$OPTARG;; m) T2G=$OPTARG;; b) AF=$OPTARG;; s) TH=$OPTARG;; *) exit 1;;
  esac
done
: "${FASTQ:?-f reads.fastq.gz}" "${WL:?-w whitelist}" "${TX:?-x transcripts.fa}" "${OUT:?-o out_dir}"
case $PLATFORM in ONT) PRESET=map-ont;; PacBio) PRESET=map-hifi;; *) echo "-p must be ONT or PacBio" >&2; exit 1;; esac
mkdir -p "$OUT"/{reference,blaze,align,rad,quant,logs}
echo "[1/8] transcript map";      if [ -n "$T2G" ]; then cut -f1,2 "$T2G" > "$OUT/reference/txp_map.tsv"; else python "$HERE/../workflow/scripts/make_t2t_map.py" -i "$TX" -o "$OUT/reference/txp_map.tsv"; fi
echo "[2/8] minimap2 index";      minimap2 -t "$THREADS" -x $PRESET -d "$OUT/reference/transcripts.mmi" "$TX" > "$OUT/logs/minimap2_index.log" 2>&1
echo "[3/8] BLAZE";               blaze --full-bc-whitelist "$WL" --expect-cells "$CELLS" --output-prefix "$OUT/blaze/" --output-fastq matched_reads.fastq.gz --threads "$THREADS" --minimal_stdout "$FASTQ" > "$OUT/logs/blaze.log" 2>&1
echo "[4/8] minimap2";            minimap2 -t "$THREADS" -a -N 100 "$OUT/reference/transcripts.mmi" "$OUT/blaze/matched_reads.fastq.gz" 2> "$OUT/logs/minimap2.log" | samtools view -@ "$THREADS" -h -F 2052 -b -o "$OUT/align/transcripts.bam" -
echo "[5/8] CB/UMI tags";         python "$HERE/../workflow/scripts/tagged_bam_CR_UR.py" --input "$OUT/align/transcripts.bam" --output "$OUT/align/transcripts.tagged.bam" > "$OUT/logs/tag_bam.log" 2>&1
echo "[6/8] convert";             "$AF" convert --threads "$THREADS" --filter_best --score-threshold "$TH" --bam "$OUT/align/transcripts.tagged.bam" --output "$OUT/rad/map.rad" > "$OUT/logs/convert.log" 2>&1
echo "[7/8] permit list, collate"; "$AF" generate-permit-list --threads "$THREADS" --input "$OUT/rad" --expected-ori $ORI --output-dir "$OUT/rad" --unfiltered-pl "$WL" > "$OUT/logs/generate_permit_list.log" 2>&1
touch "$OUT/rad/unmapped_bc_count.bin"
"$AF" collate -i "$OUT/rad" -r "$OUT/rad" -t "$THREADS" > "$OUT/logs/collate.log" 2>&1
echo "[8/8] quant";               "$AF" quant -i "$OUT/rad" -m "$OUT/reference/txp_map.tsv" -t "$THREADS" -o "$OUT/quant" -r parsimony-em \
  --lambda-size=0.0 --tau-delta=0.0 --dedup-model bf --ed-distance-model general --ed-upper-bound 2 --dedup-edges bayes --ed-far-window 10 \
  --model-truncation --truncation-interval-mode all --truncation-end-kernels absolute --truncation-cell-prior $PRIOR --truncation-positional-power $POWER > "$OUT/logs/quant.log" 2>&1
echo "done: $OUT/quant/alevin/quants_mat.mtx"
