#!/usr/bin/env bash
# Run `matpredict detect --phylum auto` on every sample in samples.txt as a
# SLURM array. Submit from the repository root:
#     sbatch --array=1-$(grep -cv '^\s*#\|^\s*$' samples.txt) scripts/03_run_matpredict.sbatch
# or use scripts/03_submit.sh, which computes the array size.
# One task = one genome. Output: results/matpredict/<sample>/.
#SBATCH --job-name=matpredict
#SBATCH --output=logs/matpredict_%A_%a.out
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=4:00:00
set -euo pipefail
source "${PROJECT_DIR:-${SLURM_SUBMIT_DIR:-$PWD}}/config.sh"

N="${SLURM_ARRAY_TASK_ID:-1}"
abbrev=$(grep -v '^\s*#' "$SAMPLES_FILE" | grep -v '^\s*$' | sed -n "${N}p")
[ -n "$abbrev" ] || { echo "no sample at line $N of $SAMPLES_FILE" >&2; exit 1; }
R="FungiDB-${FUNGIDB_RELEASE}_${abbrev}"
GENOME_GZ="$DATA_DIR/genome/${R}_Genome.fasta.gz"
PROTEIN_GZ="$DATA_DIR/protein/${R}_AnnotatedProteins.fasta.gz"
OUT="$RUNS_DIR/$abbrev"
CPUS="${SLURM_CPUS_PER_TASK:-4}"

[ -s "$GENOME_GZ" ] || { echo "missing $GENOME_GZ (run scripts/02_fetch_fungidb_data.sh)" >&2; exit 1; }
if [ -s "$OUT/detection_report.yaml" ]; then
    echo "$abbrev: $OUT/detection_report.yaml exists, skipping (delete it to re-run)"; exit 0
fi

# Node-local scratch for the uncompressed FASTA and the work files.
# $SCRATCH is set by SLURM on UCR HPCC; fall back to a temp dir for local runs.
WORK="${SCRATCH:-$(mktemp -d)}/matpredict_$abbrev"
mkdir -p "$WORK" "$OUT" "$LOGS_DIR"
zcat "$GENOME_GZ" > "$WORK/$abbrev.genome.fasta"
PROTEIN_ARG=()
if [ -s "$PROTEIN_GZ" ]; then
    zcat "$PROTEIN_GZ" > "$WORK/$abbrev.proteins.fasta"
    PROTEIN_ARG=(--proteins "$WORK/$abbrev.proteins.fasta")
fi

echo "$abbrev: start $(date)  cpus=$CPUS  work=$WORK"
export MATPREDICT_POLISH_WORKERS="$CPUS"   # MATPredict has no --threads flag; this sets its worker pool
matpredict_run matpredict detect \
    --genome "$WORK/$abbrev.genome.fasta" \
    "${PROTEIN_ARG[@]}" \
    --phylum auto \
    --sample "$abbrev" \
    --out-dir "$OUT" \
    2>&1 | tee "$OUT/detect.log"
echo "$abbrev: done $(date)"
rm -rf "$WORK"
ls -l "$OUT"
