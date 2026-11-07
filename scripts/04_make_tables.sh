#!/usr/bin/env bash
# Combine every detect run into genomes.tsv and loci.tsv, plus a summary.
set -euo pipefail
source "${PROJECT_DIR:-$PWD}/config.sh"
mkdir -p "$TABLES_DIR"

# Sample sheet: the organism name from each genome FASTA header
# ("| organism=Genus_species_strain |"), joined into genomes.tsv by `report table --samples`.
SHEET="$TABLES_DIR/sample_sheet.tsv"
printf "ASMID\torganism\n" > "$SHEET"
grep -v '^\s*#' "$SAMPLES_FILE" | grep -v '^\s*$' | while read -r s; do
    g=$(genome_file "$s"); org=""
    # `|| true`: head closes the pipe early, so gzip exits non-zero under pipefail
    [ -z "$g" ] || org=$( (gzip -cdf "$g" 2>/dev/null || true) | head -1 | sed -n 's/.*organism=\([^ |]*\).*/\1/p' | tr '_' ' ') || true
    printf "%s\t%s\n" "$s" "$org"
done >> "$SHEET"
matpredict_run matpredict report table --runs "$RUNS_DIR" --out-dir "$TABLES_DIR" --samples "$SHEET" --samples-key ASMID
matpredict_run matpredict report summary --table-dir "$TABLES_DIR" --by phylum
echo "--- $TABLES_DIR:"; ls -l "$TABLES_DIR"
echo "--- genomes.tsv:"; column -t -s $'\t' "$TABLES_DIR/genomes.tsv" | cut -c1-200
