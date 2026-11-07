#!/usr/bin/env bash
# Submit scripts/03_run_matpredict.sbatch as an array sized to samples.txt.
set -euo pipefail
source "${PROJECT_DIR:-$PWD}/config.sh"
mkdir -p "$LOGS_DIR" "$RUNS_DIR"
n=$(grep -v '^\s*#' "$SAMPLES_FILE" | grep -cv '^\s*$')
[ "$n" -gt 0 ] || { echo "samples.txt is empty" >&2; exit 1; }
sbatch -p "$SLURM_PARTITION" --cpus-per-task="$SLURM_CPUS" --mem="$SLURM_MEM" --time="$SLURM_TIME" \
    --array="1-$n" "$PROJECT_DIR/scripts/03_run_matpredict.sbatch"
