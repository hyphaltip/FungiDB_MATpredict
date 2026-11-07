#!/usr/bin/env bash
# Publish one release's results into docs/release-<N>/ (per-genome MATPredict
# HTML reports, all tables, an index), then rebuild the top docs/index.md
# across every published release (primary release first). Commit and push
# docs/ to publish on GitHub Pages.
set -euo pipefail
source "${PROJECT_DIR:-$PWD}/config.sh"
mkdir -p "$DOCS_RELEASE_DIR/reports" "$DOCS_RELEASE_DIR/tables"

# 1. copy this release's reports and tables
for run in "$RUNS_DIR"/*/; do
    [ -d "$run" ] || continue
    s=$(basename "$run")
    [ -s "$run/report.html" ] || continue
    mkdir -p "$DOCS_RELEASE_DIR/reports/$s"
    GFF="$DATA_DIR/gff/FungiDB-${FUNGIDB_RELEASE}_${s}.gff.gz"
    if [ -s "$GFF" ]; then
        # add a FungiDB gene column (links to fungidb.org/gene/<ID>) to the report's Gene evidence tables
        python3 "$PROJECT_DIR/scripts/annotate_report.py" --report "$run/report.html" --gff "$GFF" \
            --loci-gff3 "$run/detected_loci.gff3" --annotation-link annotation.html \
            --out "$DOCS_RELEASE_DIR/reports/$s/index.html" 2>/dev/null
    else
        cp "$run/report.html" "$DOCS_RELEASE_DIR/reports/$s/index.html"
    fi
    cp "$run/detected_loci.gff3" "$DOCS_RELEASE_DIR/reports/$s/" 2>/dev/null || true
    cp "$run/detection_report.yaml" "$DOCS_RELEASE_DIR/reports/$s/" 2>/dev/null || true
    [ -s "$COMPARE_DIR/annotation/$s.html" ] && cp "$COMPARE_DIR/annotation/$s.html" "$DOCS_RELEASE_DIR/reports/$s/annotation.html"
done
for f in "$TABLES_DIR"/*.tsv "$TABLES_DIR"/*.md "$COMPARE_DIR"/*.tsv; do
    [ -s "$f" ] && cp "$f" "$DOCS_RELEASE_DIR/tables/"
done
[ -s "$PROJECT_DIR/results/release_changes.tsv" ] && cp "$PROJECT_DIR/results/release_changes.tsv" "$DOCS_DIR/"

# 2. this release's index and the top index
VERSIONS="MATPredict $(cat "$TOOLS_DIR/MATPredict.version" 2>/dev/null || echo '?'), $(cat "$TOOLS_DIR/gmg.version" 2>/dev/null || echo 'GmG ?')"
python3 "$PROJECT_DIR/scripts/build_docs_index.py" \
    --docs "$DOCS_DIR" --results "$PROJECT_DIR/results" --release "$FUNGIDB_RELEASE" \
    --releases $FUNGIDB_RELEASES --versions "$VERSIONS"
echo "wrote $DOCS_RELEASE_DIR/index.md and $DOCS_DIR/index.md"
