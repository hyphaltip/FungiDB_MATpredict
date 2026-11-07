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


def phylum_call(run_dir: Path) -> dict:
    """GmG's phylum call from detection_report.yaml: {phylum, reason, note}. Empty when unavailable."""
    y = run_dir / "detection_report.yaml"
    if not y.is_file():
        return {}
    try:
        import yaml
        doc = yaml.safe_load(open(y))
    except Exception:
        return {}

    def find(d):
        if isinstance(d, dict):
            if "phylum_guess" in d and isinstance(d["phylum_guess"], dict):
                return d["phylum_guess"]
            for v in d.values():
                r = find(v)
                if r:
                    return r
        return None
    g = find(doc) or {}
    return {"phylum": g.get("phylum") or "", "reason": g.get("reason") or "", "note": g.get("note") or ""}


def run_provenance(run_dir: Path) -> dict:
    """Tool versions recorded in detection_report.yaml: matpredict (git describe), commit, database sha256, gmg."""
    y = run_dir / "detection_report.yaml"
    if not y.is_file():
        return {}
    try:
        import yaml
        run = (yaml.safe_load(open(y)) or {}).get("run") or {}
    except Exception:
        return {}
    code = run.get("code") or {}
    db = run.get("database") or {}
    guess = (run.get("phylum_guess") or {})
    ver = code.get("git_describe") or run.get("matpredict_version") or ""
    if code.get("dirty"):
        ver += " (dirty)"
    return {"matpredict": ver, "commit": (code.get("commit") or "")[:10],
            "database": (db.get("content_sha256") or "")[:10], "gmg": guess.get("gmg_version") or "",
            "gmg_index": guess.get("index") or "", "taxonomy": run.get("taxonomy_source") or ""}


def provenance_summary(runs: Path) -> tuple[str, dict]:
    """One line naming every distinct MATPredict version (with run counts) used in a release, plus the counts."""
    from collections import Counter
    vers, extra = Counter(), {}
    if runs.is_dir():
        for d in sorted(runs.iterdir()):
            p = run_provenance(d)
            if p.get("matpredict"):
                key = f"MATPredict {p['matpredict']} ({p['commit']}), database {p['database']}"
                vers[key] += 1
                extra = p
    if not vers:
        return "no runs", {}
    parts = [f"{k}: {n} genomes" for k, n in vers.most_common()]
    line = "; ".join(parts)
    if extra.get("gmg"):
        line += f". {extra['gmg']}, index {extra['gmg_index']}. {extra['taxonomy']}"
    if len(vers) > 1:
        line = f"**Runs come from {len(vers)} different MATPredict builds.** " + line
    return line, dict(vers)


def badge(text: str, cls: str) -> str:
    return f'<span class="badge {cls}">{text}</span>' if text else "—"


def results_table(res: Path, docs_rel: Path, prefix: str) -> tuple[list[str], dict, dict]:
    """Markdown rows for one release (ordered by phylum, then genome), summary counts, and per-phylum counts."""
    genomes = {r["sample"]: r for r in read(res / "tables" / "genomes.tsv")}
    loci = {}
    for r in read(res / "tables" / "loci.tsv"):
        loci.setdefault(r["sample"], []).append(r)
    summary = {r["sample"]: r for r in read(res / "compare" / "summary.tsv")}
    runs = res / "matpredict"
    samples = set(genomes) | set(summary) | ({p.name for p in runs.iterdir() if p.is_dir()} if runs.is_dir() else set())
    calls = {s: phylum_call(runs / s) for s in samples}

    def phylum_of(s):
        g = genomes.get(s, {})
        return (calls[s].get("phylum") or g.get("phylum")
                or (g.get("families_attempted") or "").split(":")[0] or "unknown")

    ordered = sorted(samples, key=lambda s: (phylum_of(s), s.lower()))
    lines = ["| Genome | Organism | Phylum (GmG) | MAT idiomorph | Confidence | Genes found | Annotated | Not annotated | Locus | Files |",
             "| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |"]
    counts = {"genomes": len(samples), "with_locus": 0, "not_searched": 0, "failed": 0, "unannotated_genes": 0}
    by_phylum = {}
    for s in ordered:
        g = genomes.get(s, {}); c = summary.get(s, {}); L = loci.get(s, []); ph = phylum_of(s)
        bp = by_phylum.setdefault(ph, {"genomes": 0, "called": 0, "not_searched": 0})
        bp["genomes"] += 1
        has_report = (docs_rel / "reports" / s / "index.html").is_file()
        link = f"[{s}]({prefix}reports/{s}/)" if has_report else s
        ann = f"[annotation]({prefix}reports/{s}/annotation.html) " if (docs_rel / "reports" / s / "annotation.html").is_file() else ""
        files = (f'<span class="links">{ann}[GFF3]({prefix}reports/{s}/detected_loci.gff3) [YAML]({prefix}reports/{s}/detection_report.yaml)</span>'
                 if has_report else (g.get("status") or "no run"))
        organism = g.get("organism") or "—"
        if g.get("status") == "failed":
            counts["failed"] += 1
        if L:
            counts["with_locus"] += 1; bp["called"] += 1
        if g.get("not_searched_reason"):
            counts["not_searched"] += 1; bp["not_searched"] += 1
        try:
            counts["unannotated_genes"] += int(c.get("n_unannotated") or 0)
        except ValueError:
            pass
        coords = "<br>".join(f'<span class="mono">{l["contig"]}:{l["start"]}-{l["end"]}</span>' for l in L) or "—"
        conf = "<br>".join(badge(l["confidence"], "conf-" + (l["confidence"] or "none")) for l in L) or "—"
        if L:
            idio = "<br>".join(f"**{l['idiomorph'] or 'undetermined'}**" for l in L)
        elif g.get("not_searched_reason"):
            idio = badge("not searched", "warn") + f' <span class="mono">{g["not_searched_reason"]}</span>'
        else:
            idio = "none called"
        found = "<br>".join(l["genes_found"].replace(";", ", ") for l in L) or "—"
        n_ann = c.get("n_annotated", "—"); n_un = c.get("n_unannotated", "—")
        if c.get("status", "ok") != "ok" and c:
            n_ann = n_un = f'<span class="mono">{c["status"]}</span>'
        lines.append(f"| {link} | *{organism}* | {ph} | {idio} | {conf} | {found} | {n_ann} | {n_un} | {coords} | {files} |")
    return lines, counts, by_phylum


def stats_strip(items: list[tuple[str, str]]) -> str:
    return '<div class="stats">' + "".join(f'<div class="stat"><b>{v}</b><span>{k}</span></div>' for k, v in items) + "</div>"


TABLE_DESC = {
    "genomes.tsv": "one row per genome: routing, idiomorph, locus count",
    "loci.tsv": "one row per called locus: family, idiomorph, confidence, genes found and missing",
    "annotated_MAT_genes.tsv": "MATPredict genes that overlap a FungiDB gene, with the FungiDB ID and product",
    "unannotated_MAT_genes.tsv": "MATPredict genes with no overlapping FungiDB gene, with nearest neighbours",
    "locus_genes.tsv": "every FungiDB gene inside a MATPredict locus span",
    "summary.tsv": "comparison counts per genome",
    "summary_by_phylum.tsv": "MATPredict recovery counts by phylum",
    "sample_sheet.tsv": "sample to organism name, from the FungiDB FASTA headers",
    "annotation_verdicts.tsv": "per MATPredict gene: annotated, partial or novel against the FungiDB GFF, with the matched gene",
}

METHOD = ("Each genome is routed to a phylum by [GmG](https://github.com/stajichlab/GmG) "
          "(`matpredict detect --phylum auto`), then [MATPredict](https://github.com/stajichlab/MATPredict) "
          "locates the MAT locus and models its genes. Each modelled gene is matched to the FungiDB GFF by "
          "coordinate overlap on the same contig. Every step is a script in the "
          "[repository](https://github.com/hyphaltip/FungiDB_MATpredict); see its README to repeat the run.")


def write_release(docs: Path, results: Path, rel: str, versions: str) -> dict:
    docs_rel = docs / f"release-{rel}"
    res = results / f"release-{rel}"
    rows, counts, by_phylum = results_table(res, docs_rel, "")
    phyla = sorted(by_phylum.items(), key=lambda kv: -kv[1]["genomes"])
    prov_line, prov = provenance_summary(res / "matpredict")
    counts["versions"] = ", ".join(k.split(" (")[0].replace("MATPredict ", "") for k in prov) or "—"
    out = ["---", f"title: FungiDB release {rel}", "---",
           f"# FungiDB release {rel}", "",
           f"Generated {datetime.date.today()}. Software recorded in the run reports: {prov_line}", "",
           "Rows are ordered by GmG phylum, then genome; "
           "click a column header to sort, or type in the filter box (for example a phylum name).", "",
           stats_strip([("genomes", counts["genomes"]), ("MAT locus called", counts["with_locus"]),
                        ("not searched", counts["not_searched"]), ("failed", counts["failed"]),
                        ("MAT genes not annotated", counts["unannotated_genes"])]), "",
           "Download all results: [genomes.tsv](tables/genomes.tsv) (one row per genome) · "
           "[loci.tsv](tables/loci.tsv) (one row per MAT locus) · "
           "[annotation_verdicts.tsv](tables/annotation_verdicts.tsv) (one row per MAT gene vs the FungiDB GFF) · "
           "[all tables](#tables)", "",
           "## By phylum", "",
           "| Phylum (GmG call) | Genomes | MAT locus called | Not searched |", "| --- | --- | --- | --- |",
           *[f"| {ph} | {c['genomes']} | {c['called']} | {c['not_searched']} |" for ph, c in phyla], "",
           "## Genomes", "", *rows, "",
           "## Tables", "", "| File | Contents |", "| --- | --- |"]
    for p in sorted((docs_rel / "tables").glob("*.tsv")):
        out.append(f"| [{p.name}](tables/{p.name}) | {TABLE_DESC.get(p.name, '')} |")
    out += ["", "## Method", "", METHOD, ""]
    docs_rel.mkdir(parents=True, exist_ok=True)
    (docs_rel / "index.md").write_text("\n".join(out))
    return counts


def write_top(docs: Path, results: Path, releases: list[str], versions: str, counts_by_rel: dict):
    out = ["---", "title: Releases", "---", "# Releases", "",
           "[MATPredict](https://github.com/stajichlab/MATPredict) mating-type (MAT) locus calls on FungiDB "
           "genomes, compared with the FungiDB gene annotation. Results are kept per FungiDB release; "
           f"release {releases[0]} is the primary set. Generated {datetime.date.today()}; the MATPredict version is the one "
           "recorded in each run's report (see a release page for the commit and database checksum).", "",
           "## Releases", "",
           "| Release | MATPredict version | Genomes run | MAT locus called | Not searched | Failed | MAT genes without annotation |",
           "| --- | --- | --- | --- | --- | --- | --- |"]
    for i, rel in enumerate(releases):
        c = counts_by_rel.get(rel)
        if c is None:
            out.append(f"| release {rel} | | not run yet | | | | |")
            continue
        tag = " (primary)" if i == 0 else ""
        out.append(f"| [release {rel}](release-{rel}/){tag} | {c.get('versions', '—')} | {c['genomes']} | {c['with_locus']} | "
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
