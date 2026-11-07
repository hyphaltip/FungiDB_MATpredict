#!/usr/bin/env bash
# Compare MATPredict gene calls with the FungiDB GFF annotation for every sample.
set -euo pipefail
source "${PROJECT_DIR:-$PWD}/config.sh"
mkdir -p "$COMPARE_DIR"
matpredict_run python "$PROJECT_DIR/scripts/compare_matpredict_to_gff.py" \
    --runs "$RUNS_DIR" --gff-dir "$DATA_DIR/gff" --release "$FUNGIDB_RELEASE" \
    --samples "$SAMPLES_FILE" --out-dir "$COMPARE_DIR"
# Per-genome pages: MAT locus genes above the FungiDB annotation of the same region, with a verdict per gene.
python3 "$PROJECT_DIR/scripts/annotation_view.py" \
    --runs "$RUNS_DIR" --gff-dir "$DATA_DIR/gff" --release "$FUNGIDB_RELEASE" \
    --samples "$SAMPLES_FILE" --organisms "$TABLES_DIR/sample_sheet.tsv" \
    --out-dir "$COMPARE_DIR/annotation" --report-rel "index.html"
cp "$COMPARE_DIR/annotation/annotation_verdicts.tsv" "$COMPARE_DIR/"
echo "--- $COMPARE_DIR:"; ls -l "$COMPARE_DIR"
echo "--- summary:"; column -t -s $'\t' "$COMPARE_DIR/summary.tsv"
