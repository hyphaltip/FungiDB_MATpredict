#!/usr/bin/env bash
# Copy the per-genome MATPredict HTML reports and all result tables into docs/,
# and rebuild docs/index.md. Commit and push docs/ to publish on GitHub Pages.
set -euo pipefail
source "${PROJECT_DIR:-$PWD}/config.sh"
DOCS="$PROJECT_DIR/docs"
mkdir -p "$DOCS/reports" "$DOCS/tables"

INDEX="$DOCS/index.md"
{
cat <<'HDR'
# FungiDB MATpredict

[MATPredict](https://github.com/stajichlab/MATPredict) mating-type (MAT) locus
calls on FungiDB genomes, compared with the FungiDB gene annotation.
Every step is a script in the
[repository](https://github.com/hyphaltip/FungiDB_MATpredict); see its README.

HDR
echo "Generated $(date -u +%Y-%m-%d) with MATPredict $(cat "$TOOLS_DIR/MATPredict.version" 2>/dev/null || echo '?'), $(cat "$TOOLS_DIR/gmg.version" 2>/dev/null || echo 'GmG ?'), FungiDB release $FUNGIDB_RELEASE."
echo
echo "## Per-genome reports"
echo
for run in "$RUNS_DIR"/*/; do
    s=$(basename "$run")
    if [ -s "$run/report.html" ]; then
        mkdir -p "$DOCS/reports/$s"
        cp "$run/report.html" "$DOCS/reports/$s/index.html"
        cp "$run/detected_loci.gff3" "$DOCS/reports/$s/detected_loci.gff3" 2>/dev/null || true
        cp "$run/detection_report.yaml" "$DOCS/reports/$s/detection_report.yaml" 2>/dev/null || true
        echo "- [$s](reports/$s/) ([GFF3](reports/$s/detected_loci.gff3), [YAML](reports/$s/detection_report.yaml))"
    else
        echo "- $s: no report (run failed or not run)"
    fi
done
echo
echo "## Tables"
echo
for f in "$TABLES_DIR"/*.tsv "$TABLES_DIR"/*.md "$COMPARE_DIR"/*.tsv; do
    [ -s "$f" ] || continue
    cp "$f" "$DOCS/tables/"
    echo "- [$(basename "$f")](tables/$(basename "$f"))"
done
if [ -s "$COMPARE_DIR/summary.tsv" ]; then
    echo; echo "## Comparison summary"; echo
    awk -F'\t' 'NR==1{printf "|"; for(i=1;i<=NF;i++) printf " %s |",$i; printf "\n|"; for(i=1;i<=NF;i++) printf " --- |"; printf "\n"; next}
                {printf "|"; for(i=1;i<=NF;i++) printf " %s |",$i; printf "\n"}' "$COMPARE_DIR/summary.tsv"
fi
} > "$INDEX"
echo "wrote $INDEX"; ls -R "$DOCS" | head -40
