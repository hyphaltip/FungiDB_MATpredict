#!/usr/bin/env python3
"""Write docs/release-<N>/index.md for one release and docs/index.md across releases.

Called by scripts/06_publish_docs.sh. Standard library only.
"""
import argparse
import csv
import datetime
from pathlib import Path


def read(p: Path):
    return list(csv.DictReader(open(p), delimiter="\t")) if p.is_file() else []


def results_table(res: Path, docs_rel: Path, prefix: str) -> tuple[list[str], dict]:
    """Markdown rows for one release and its summary counts."""
    genomes = {r["sample"]: r for r in read(res / "tables" / "genomes.tsv")}
    loci = {}
    for r in read(res / "tables" / "loci.tsv"):
        loci.setdefault(r["sample"], []).append(r)
    summary = {r["sample"]: r for r in read(res / "compare" / "summary.tsv")}
    runs = res / "matpredict"
    samples = sorted(set(genomes) | set(summary) | ({p.name for p in runs.iterdir() if p.is_dir()} if runs.is_dir() else set()))
    lines = ["| Genome | Report | Phylum | MAT idiomorph | Confidence | Genes found | Annotated in FungiDB | Not annotated | Locus coordinates |",
             "| --- | --- | --- | --- | --- | --- | --- | --- | --- |"]
    counts = {"genomes": len(samples), "with_locus": 0, "not_searched": 0, "failed": 0, "unannotated_genes": 0}
    for s in samples:
        g = genomes.get(s, {}); c = summary.get(s, {}); L = loci.get(s, [])
        has_report = (docs_rel / "reports" / s / "index.html").is_file()
        link = f"[{s}]({prefix}reports/{s}/)" if has_report else s
        rep = (f"[HTML]({prefix}reports/{s}/) · [GFF3]({prefix}reports/{s}/detected_loci.gff3) · [YAML]({prefix}reports/{s}/detection_report.yaml)"
               if has_report else (g.get("status") or "no run"))
        if g.get("status") == "failed":
            counts["failed"] += 1
        if L:
            counts["with_locus"] += 1
        if g.get("not_searched_reason"):
            counts["not_searched"] += 1
        try:
            counts["unannotated_genes"] += int(c.get("n_unannotated") or 0)
        except ValueError:
            pass
        phylum = g.get("phylum") or (g.get("families_attempted") or "").split(":")[0] or "—"
        coords = "<br>".join(f"{l['contig']}:{l['start']}-{l['end']}" for l in L) or "—"
        conf = "<br>".join(l["confidence"] for l in L) or "—"
        idio = g.get("idiomorphs") or g.get("not_searched_reason") or "—"
        found = "<br>".join(l["genes_found"].replace(";", ", ") for l in L) or "—"
        lines.append(f"| {link} | {rep} | {phylum} | {idio} | {conf} | {found} | "
                     f"{c.get('n_annotated', '—')} | {c.get('n_unannotated', '—')} | {coords} |")
    return lines, counts


TABLE_DESC = {
    "genomes.tsv": "one row per genome: routing, idiomorph, locus count",
    "loci.tsv": "one row per called locus: family, idiomorph, confidence, genes found and missing",
    "annotated_MAT_genes.tsv": "MATPredict genes that overlap a FungiDB gene, with the FungiDB ID and product",
    "unannotated_MAT_genes.tsv": "MATPredict genes with no overlapping FungiDB gene, with nearest neighbours",
    "locus_genes.tsv": "every FungiDB gene inside a MATPredict locus span",
    "summary.tsv": "comparison counts per genome",
    "summary_by_phylum.tsv": "MATPredict recovery counts by phylum",
}

METHOD = ("Each genome is routed to a phylum by [GmG](https://github.com/stajichlab/GmG) "
          "(`matpredict detect --phylum auto`), then [MATPredict](https://github.com/stajichlab/MATPredict) "
          "locates the MAT locus and models its genes. Each modelled gene is matched to the FungiDB GFF by "
          "coordinate overlap on the same contig. Every step is a script in the "
          "[repository](https://github.com/hyphaltip/FungiDB_MATpredict); see its README to repeat the run.")


def write_release(docs: Path, results: Path, rel: str, versions: str) -> dict:
    docs_rel = docs / f"release-{rel}"
    res = results / f"release-{rel}"
    rows, counts = results_table(res, docs_rel, "")
    out = [f"# FungiDB release {rel}", "",
           f"[All releases](../) · Generated {datetime.date.today()} with {versions}, FungiDB release {rel}.", "",
           f"{counts['genomes']} genomes run, {counts['with_locus']} with a MAT locus called, "
           f"{counts['not_searched']} not searched, {counts['failed']} failed.", "",
           "## Results", "", "Click a genome for MATPredict's full report.", "", *rows, "",
           "## Tables", "", "| File | Contents |", "| --- | --- |"]
    for p in sorted((docs_rel / "tables").glob("*.tsv")):
        out.append(f"| [{p.name}](tables/{p.name}) | {TABLE_DESC.get(p.name, '')} |")
    out += ["", "## Method", "", METHOD, ""]
    docs_rel.mkdir(parents=True, exist_ok=True)
    (docs_rel / "index.md").write_text("\n".join(out))
    return counts


def write_top(docs: Path, results: Path, releases: list[str], versions: str, counts_by_rel: dict):
    out = ["# FungiDB MATpredict", "",
           "[MATPredict](https://github.com/stajichlab/MATPredict) mating-type (MAT) locus calls on FungiDB "
           "genomes, compared with the FungiDB gene annotation. Results are kept per FungiDB release; "
           f"release {releases[0]} is the primary set.", "",
           f"Generated {datetime.date.today()} with {versions}.", "",
           "## Releases", "",
           "| Release | Genomes run | MAT locus called | Not searched | Failed | MAT genes without annotation |",
           "| --- | --- | --- | --- | --- | --- |"]
    for i, rel in enumerate(releases):
        c = counts_by_rel.get(rel)
        if c is None:
            out.append(f"| release {rel} | not run yet | | | | |")
            continue
        tag = " (primary)" if i == 0 else ""
        out.append(f"| [release {rel}](release-{rel}/){tag} | {c['genomes']} | {c['with_locus']} | "
                   f"{c['not_searched']} | {c['failed']} | {c['unannotated_genes']} |")
    changes = results / "release_changes.tsv"
    if changes.is_file():
        rows = read(changes)
        by_status = {}
        for r in rows:
            by_status.setdefault(r["status"], []).append(r["sample"])
        out += ["", "## Changes between releases", "",
                f"From [release_changes.tsv](release_changes.tsv), built by `scripts/compare_releases.sh` "
                f"(releases {', '.join(releases)}).", "",
                "| Status | Genomes | Samples |", "| --- | --- | --- |"]
        for st, ss in sorted(by_status.items()):
            shown = ", ".join(ss[:12]) + (f", … ({len(ss)} total)" if len(ss) > 12 else "")
            out.append(f"| {st} | {len(ss)} | {shown} |")
    out += ["", "## Method", "", METHOD, ""]
    (docs / "index.md").write_text("\n".join(out))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--docs", type=Path, required=True)
    ap.add_argument("--results", type=Path, required=True)
    ap.add_argument("--release", required=True, help="the release whose index to (re)build")
    ap.add_argument("--releases", nargs="+", required=True, help="all releases, primary first")
    ap.add_argument("--versions", default="")
    a = ap.parse_args()
    counts = {}
    for rel in a.releases:
        if (a.results / f"release-{rel}" / "tables" / "genomes.tsv").is_file() or rel == a.release:
            counts[rel] = write_release(a.docs, a.results, rel, a.versions)
    write_top(a.docs, a.results, a.releases, a.versions, counts)


if __name__ == "__main__":
    main()
