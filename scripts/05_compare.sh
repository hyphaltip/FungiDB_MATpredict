#!/usr/bin/env bash
# Compare MATPredict gene calls with the FungiDB GFF annotation for every sample.
set -euo pipefail
source "${PROJECT_DIR:-$PWD}/config.sh"
mkdir -p "$COMPARE_DIR"
matpredict_run python "$PROJECT_DIR/scripts/compare_matpredict_to_gff.py" \
    --runs "$RUNS_DIR" --gff-dir "$DATA_DIR/gff" --release "$FUNGIDB_RELEASE" \
    --samples "$SAMPLES_FILE" --out-dir "$COMPARE_DIR"
echo "--- $COMPARE_DIR:"; ls -l "$COMPARE_DIR"
echo "--- summary:"; column -t -s $'\t' "$COMPARE_DIR/summary.tsv"
