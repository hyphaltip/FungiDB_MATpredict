#!/usr/bin/env bash
# Copy the per-genome MATPredict HTML reports and all result tables into docs/,
# and rebuild docs/index.md. Commit and push docs/ to publish on GitHub Pages.
set -euo pipefail
source "${PROJECT_DIR:-$PWD}/config.sh"
DOCS="$PROJECT_DIR/docs"
mkdir -p "$DOCS/reports" "$DOCS/tables"

INDEX="$DOCS/index.md"
# copy per-genome reports and all tables into docs/
for run in "$RUNS_DIR"/*/; do
    s=$(basename "$run")
    [ -s "$run/report.html" ] || continue
    mkdir -p "$DOCS/reports/$s"
    cp "$run/report.html" "$DOCS/reports/$s/index.html"
    cp "$run/detected_loci.gff3" "$DOCS/reports/$s/" 2>/dev/null || true
    cp "$run/detection_report.yaml" "$DOCS/reports/$s/" 2>/dev/null || true
done
for f in "$TABLES_DIR"/*.tsv "$TABLES_DIR"/*.md "$COMPARE_DIR"/*.tsv; do
    [ -s "$f" ] && cp "$f" "$DOCS/tables/"
done

# build index.md: one results table, one row per genome, then the table links
VERSIONS="MATPredict $(cat "$TOOLS_DIR/MATPredict.version" 2>/dev/null || echo '?'), $(cat "$TOOLS_DIR/gmg.version" 2>/dev/null || echo 'GmG ?'), FungiDB release $FUNGIDB_RELEASE"
python3 - "$RUNS_DIR" "$TABLES_DIR" "$COMPARE_DIR" "$DOCS" "$VERSIONS" <<'PY' > "$INDEX"
import csv, sys, datetime
from pathlib import Path
runs, tables, compare, docs = map(Path, sys.argv[1:5]); versions = sys.argv[5]

def read(p):
    return list(csv.DictReader(open(p), delimiter="\t")) if p.is_file() else []

genomes = {r["sample"]: r for r in read(tables / "genomes.tsv")}
loci = {}
for r in read(tables / "loci.tsv"):
    loci.setdefault(r["sample"], []).append(r)
summary = {r["sample"]: r for r in read(compare / "summary.tsv")}
samples = sorted(set(genomes) | set(summary) | {p.name for p in runs.iterdir() if p.is_dir()})

print("## Results\n")
print(f"Generated {datetime.date.today()} with {versions}. One row per genome; "
      "click the sample for MATPredict's full report.\n")
print("| Genome | Report | Phylum | MAT idiomorph | Confidence | Genes found | Annotated in FungiDB | Not annotated | Locus coordinates |")
print("| --- | --- | --- | --- | --- | --- | --- | --- | --- |")
for s in samples:
    g = genomes.get(s, {}); c = summary.get(s, {}); L = loci.get(s, [])
    has_report = (docs / "reports" / s / "index.html").is_file()
    link = f"[{s}](reports/{s}/)" if has_report else s
    rep = (f"[HTML](reports/{s}/) · [GFF3](reports/{s}/detected_loci.gff3) · [YAML](reports/{s}/detection_report.yaml)"
           if has_report else g.get("status", "no run"))
    coords = "<br>".join(f"{l['contig']}:{l['start']}-{l['end']}" for l in L) or "—"
    conf = "<br>".join(l["confidence"] for l in L) or "—"
    idio = g.get("idiomorphs") or g.get("not_searched_reason") or "—"
    found = "<br>".join(l["genes_found"].replace(";", ", ") for l in L) or "—"
    phylum = g.get("phylum") or (g.get("families_attempted") or "").split(":")[0] or "—"
    print(f"| {link} | {rep} | {phylum} | {idio} | {conf} | {found} | "
          f"{c.get('n_annotated', '—')} | {c.get('n_unannotated', '—')} | {coords} |")

print("\n## Tables\n")
desc = {
    "genomes.tsv": "one row per genome: routing, idiomorph, locus count",
    "loci.tsv": "one row per called locus: family, idiomorph, confidence, genes found and missing",
    "annotated_MAT_genes.tsv": "MATPredict genes that overlap a FungiDB gene, with the FungiDB ID and product",
    "unannotated_MAT_genes.tsv": "MATPredict genes with no overlapping FungiDB gene, with nearest neighbours",
    "locus_genes.tsv": "every FungiDB gene inside a MATPredict locus span",
    "summary.tsv": "comparison counts per genome",
    "summary_by_phylum.tsv": "MATPredict recovery counts by phylum",
}
print("| File | Contents |"); print("| --- | --- |")
for p in sorted((docs / "tables").glob("*.tsv")):
    print(f"| [{p.name}](tables/{p.name}) | {desc.get(p.name, '')} |")

print("\n## Method\n")
print("Each genome is routed to a phylum by [GmG](https://github.com/stajichlab/GmG) (`matpredict detect --phylum auto`), "
      "then [MATPredict](https://github.com/stajichlab/MATPredict) locates the MAT locus and models its genes. "
      "Each modelled gene is matched to the FungiDB GFF by coordinate overlap on the same contig. "
      "Every step is a script in the [repository](https://github.com/hyphaltip/FungiDB_MATpredict); "
      "see its README to repeat the run.")
PY
echo "wrote $INDEX"; ls -R "$DOCS" | head -40
