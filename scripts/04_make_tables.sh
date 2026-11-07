#!/usr/bin/env bash
# Combine every detect run into genomes.tsv and loci.tsv, plus a summary.
set -euo pipefail
source "${PROJECT_DIR:-$PWD}/config.sh"
mkdir -p "$TABLES_DIR"
matpredict_run matpredict report table --runs "$RUNS_DIR" --out-dir "$TABLES_DIR"
matpredict_run matpredict report summary --table-dir "$TABLES_DIR" --by phylum
echo "--- $TABLES_DIR:"; ls -l "$TABLES_DIR"
echo "--- genomes.tsv:"; column -t -s $'\t' "$TABLES_DIR/genomes.tsv" | cut -c1-200
