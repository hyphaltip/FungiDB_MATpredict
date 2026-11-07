#!/usr/bin/env bash
# Submit scripts/03_run_matpredict.sbatch as an array sized to samples.txt.
set -euo pipefail
source "${PROJECT_DIR:-$PWD}/config.sh"
mkdir -p "$LOGS_DIR" "$RUNS_DIR"
n=$(grep -v '^\s*#' "$SAMPLES_FILE" | grep -cv '^\s*$')
[ "$n" -gt 0 ] || { echo "$SAMPLES_FILE is empty" >&2; exit 1; }
tasks=$(( (n + SAMPLES_PER_TASK - 1) / SAMPLES_PER_TASK ))
echo "release $FUNGIDB_RELEASE: $n samples, $SAMPLES_PER_TASK per task -> array of $tasks tasks"
sbatch -p "$SLURM_PARTITION" --job-name="matpredict-r$FUNGIDB_RELEASE" \
    --output="$LOGS_DIR/matpredict_%A_%a.out" --export=ALL,FUNGIDB_RELEASE="$FUNGIDB_RELEASE" \
    --cpus-per-task="$SLURM_CPUS" --mem="$SLURM_MEM" --time="$SLURM_TIME" \
    --array="1-$tasks" "$PROJECT_DIR/scripts/03_run_matpredict.sbatch"
