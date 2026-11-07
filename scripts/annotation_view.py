#!/usr/bin/env python3
"""Per-genome page: the MATPredict MAT locus drawn above the FungiDB annotation of the same region.

For every sample with a MATPredict `detected_loci.gff3` and a FungiDB GFF, write
`<out-dir>/<sample>.html` with, per locus:
  * track 1: the genes MATPredict modelled (coloured by role),
  * track 2: every FungiDB gene in the locus span plus a flank, coloured by whether a
    MATPredict gene overlaps it on the same strand,
  * a table with a verdict per MATPredict gene: annotated (a FungiDB gene covers it on the
    same strand), partial (overlap but a different extent or strand), or novel (no FungiDB gene),
    and the FungiDB genes in the span that MATPredict did not model.
Also writes `<out-dir>/annotation_verdicts.tsv` across samples. Standard library only.
"""
from __future__ import annotations

import argparse
import csv
import gzip
import html
import io
import subprocess
import sys
from pathlib import Path
from urllib.parse import unquote

GENE_TYPES = {"gene", "protein_coding_gene", "ncRNA_gene", "pseudogene"}
FLANK = 5000
W, TRACK_H, LANE_H = 1100, 34, 16
ROLE_COLOR = {"core_MAT": "#b8471a", "flanking_conserved": "#2a6fb0"}
OTHER = "#6b7480"
FUNGIDB_GENE_URL = "https://fungidb.org/gene/"
VERDICT_CLASS = {"annotated": "conf-high", "partial": "conf-medium", "novel": "warn"}

CSS = """
:root{--bg:#fff;--surface-2:#f4f5f7;--ink:#1b1f24;--ink-2:#424a54;--muted:#5f6873;--line:#dde1e6;--accent:#1f5f99;
--high:#1d7a46;--high-bg:#e3f3ea;--medium:#7a4f00;--medium-bg:#fbf0d9;--warn:#8a4b00;--warn-bg:#fff6e5;
--mono:ui-monospace,SFMono-Regular,Menlo,Consolas,monospace;--sans:system-ui,-apple-system,"Segoe UI",Roboto,Arial,sans-serif}
@media (prefers-color-scheme: dark){:root{--bg:#111418;--surface-2:#1e232a;--ink:#e6e9ed;--ink-2:#c3c9d0;--muted:#9aa3ad;--line:#2c333b;--accent:#7fb6ea;
--high:#6fd39b;--high-bg:#173527;--medium:#f0c060;--medium-bg:#3a2e12;--warn:#f0b866;--warn-bg:#33270f}}
body{margin:0;background:var(--bg);color:var(--ink);font:13px/1.45 var(--sans)}
main{max-width:1140px;margin:0 auto;padding:10px 16px 40px}
h1{font-size:18px;margin:8px 0 2px} h2{font-size:14px;margin:18px 0 6px;text-transform:uppercase;letter-spacing:.05em;color:var(--ink-2)}
p{color:var(--ink-2);margin:4px 0 8px} a{color:var(--accent);text-decoration:none} a:hover{text-decoration:underline}
.nav{font-size:12.5px;margin-bottom:6px} .nav a{margin-right:12px}
svg{display:block;max-width:100%;height:auto;border:1px solid var(--line);border-radius:6px;background:var(--bg)}
.lbl{font-size:10.5px;fill:var(--ink-2);font-family:var(--sans)} .tk{font-size:11px;font-weight:600;fill:var(--ink-2)} .ax{font-size:10px;fill:var(--muted)}
table{border-collapse:collapse;width:100%;font-size:12.5px;margin:6px 0 14px}
th,td{padding:4px 8px;border-bottom:1px solid var(--line);text-align:left;vertical-align:top}
th{background:var(--surface-2);font-size:11.5px;text-transform:uppercase;letter-spacing:.04em;color:var(--ink-2)}
.badge{display:inline-block;padding:0 7px;border-radius:999px;font-size:11px;font-weight:600;border:1px solid currentColor;white-space:nowrap}
.conf-high{color:var(--high);background:var(--high-bg)} .conf-medium{color:var(--medium);background:var(--medium-bg)} .warn{color:var(--warn);background:var(--warn-bg)}
.mono{font-family:var(--mono);font-size:11.5px} .muted{color:var(--muted)}
.legend span{display:inline-block;margin-right:14px;font-size:12px} .legend i{display:inline-block;width:12px;height:10px;margin-right:4px;vertical-align:middle;border-radius:2px}
"""


def open_text(path: Path):
    if path.suffix == ".gz":
        return gzip.open(path, "rt")
    if path.suffix == ".zst":
        return io.TextIOWrapper(subprocess.Popen(["zstd", "-dc", str(path)], stdout=subprocess.PIPE).stdout)
    return open(path)


def parse_attrs(field: str) -> dict:
    out = {}
    for part in field.strip().split(";"):
        if "=" in part:
            k, v = part.split("=", 1)
            out[k.strip()] = unquote(v.strip())
    return out


def read_gff(path: Path, types):
    with open_text(path) as fh:
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue
            c = line.rstrip("\n").split("\t")
            if len(c) < 9 or (types and c[2] not in types):
                continue
            yield {"seqid": c[0].split()[0], "type": c[2], "start": int(c[3]), "end": int(c[4]),
                   "strand": c[6], "attrs": parse_attrs(c[8])}


def load_fungidb(path: Path) -> dict[str, list[dict]]:
    by = {}
    for f in read_gff(path, GENE_TYPES):
        a = f["attrs"]
        by.setdefault(f["seqid"], []).append({
            "id": a.get("ID", ""), "name": a.get("Name", ""), "type": f["type"], "start": f["start"],
            "end": f["end"], "strand": f["strand"], "desc": a.get("description", "")})
    for v in by.values():
        v.sort(key=lambda g: g["start"])
    return by


def load_matpredict(path: Path):
    loci, genes = [], []
    for f in read_gff(path, None):
        a = f["attrs"]
        if f["type"] == "MAT_locus":
            loci.append({"id": a.get("ID", ""), "family": a.get("family", ""), "idiomorph": a.get("idiomorph", ""),
                         "confidence": a.get("confidence", ""), "contig": f["seqid"], "start": f["start"], "end": f["end"]})
        elif f["type"] == "gene" and a.get("present") == "true":
            genes.append({"id": a.get("ID", ""), "locus": a.get("Parent", "").split(".seg")[0], "name": a.get("Name", ""),
                          "role": a.get("role", ""), "identity": a.get("identity", ""), "method": a.get("method", ""),
                          "contig": f["seqid"], "start": f["start"], "end": f["end"], "strand": f["strand"]})
    return loci, genes


def overlap(a, b) -> int:
    return max(0, min(a["end"], b["end"]) - max(a["start"], b["start"]) + 1)


def verdict(mp: dict, hits: list[dict]) -> tuple[str, str]:
    """annotated / partial / novel, with a short reason."""
    if not hits:
        return "novel", "no FungiDB gene overlaps this MATPredict gene"
    best = max(hits, key=lambda h: overlap(mp, h))
    frac = overlap(mp, best) / (mp["end"] - mp["start"] + 1)
    same = best["strand"] == mp["strand"]
    if same and frac >= 0.5:
        return "annotated", f"{best['id']} covers {frac:.0%} of it on the same strand"
    if not same:
        return "partial", f"{best['id']} overlaps {frac:.0%} but on the opposite strand"
    return "partial", f"{best['id']} overlaps only {frac:.0%} of it"


def gene_link(g: dict | None) -> str:
    if not g:
        return "—"
    label = html.escape(g["id"]) + (f" {html.escape(g['name'])}" if g["name"] else "")
    return f'<a href="{FUNGIDB_GENE_URL}{html.escape(g["id"])}">{label}</a>'


def assign_lanes(items: list[dict]) -> list[int]:
    """Greedy lane packing by genomic extent so labels do not collide."""
    lanes_end, out = [], []
    for g in sorted(items, key=lambda g: g["start"]):
        for i, e in enumerate(lanes_end):
            if g["start"] > e:
                lanes_end[i] = g["end"]; g["_lane"] = i; break
        else:
            lanes_end.append(g["end"]); g["_lane"] = len(lanes_end) - 1
    return [g["_lane"] for g in items]


def arrow(x1, x2, y, h, strand, color, title, label, cls="lbl"):
    tip = min(6, max(2, (x2 - x1) / 3))
    if strand == "-":
        pts = f"{x1 + tip},{y} {x2},{y} {x2},{y + h} {x1 + tip},{y + h} {x1},{y + h / 2}"
    else:
        pts = f"{x1},{y} {x2 - tip},{y} {x2},{y + h / 2} {x2 - tip},{y + h} {x1},{y + h}"
    t = html.escape(title, quote=True)
    return (f'<polygon points="{pts}" fill="{color}" opacity="0.9"><title>{t}</title></polygon>'
            f'<text class="{cls}" x="{x1}" y="{y - 2}">{html.escape(label)}</text>')


def locus_svg(locus: dict, mp_genes: list[dict], fdb_genes: list[dict], matched_ids: set[str]) -> str:
    s0, s1 = locus["start"] - FLANK, locus["end"] + FLANK
    span = max(1, s1 - s0)
    lm = 70
    x = lambda bp: lm + (bp - s0) / span * (W - lm - 10)
    mp_lanes = assign_lanes(mp_genes)
    fd_lanes = assign_lanes(fdb_genes)
    n1, n2 = (max(mp_lanes) + 1 if mp_lanes else 1), (max(fd_lanes) + 1 if fd_lanes else 1)
    y1 = 26
    y2 = y1 + n1 * (LANE_H + 14) + 24
    H = y2 + n2 * (LANE_H + 14) + 30
    out = [f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {W} {H}" width="{W}" height="{H}">']
    # locus span band
    out.append(f'<rect x="{x(locus["start"]):.1f}" y="{y1 - 14}" width="{x(locus["end"]) - x(locus["start"]):.1f}" '
               f'height="{H - y1 - 6}" fill="#1f5f99" opacity="0.06"/>')
    out.append(f'<text class="tk" x="4" y="{y1 + 8}">MATPredict</text>')
    out.append(f'<text class="tk" x="4" y="{y2 + 8}">FungiDB</text>')
    for g, lane in zip(mp_genes, mp_lanes):
        y = y1 + lane * (LANE_H + 14) + 12
        out.append(arrow(x(g["start"]), x(g["end"]), y, LANE_H - 4, g["strand"], ROLE_COLOR.get(g["role"], OTHER),
                         f"{g['name']} ({g['role']}) {g['contig']}:{g['start']}-{g['end']} {g['strand']} identity {g['identity']}",
                         g["name"]))
    for g, lane in zip(fdb_genes, fd_lanes):
        y = y2 + lane * (LANE_H + 14) + 12
        color = "#1d7a46" if g["id"] in matched_ids else OTHER
        label = g["id"] + (f" {g['name']}" if g["name"] else "")
        out.append(arrow(x(g["start"]), x(g["end"]), y, LANE_H - 4, g["strand"], color,
                         f"{label} {g['start']}-{g['end']} {g['strand']}: {g['desc']}", label))
    # axis
    ya = H - 12
    out.append(f'<line x1="{lm}" y1="{ya}" x2="{W - 10}" y2="{ya}" stroke="#9aa3ad"/>')
    step = 10 ** max(0, len(str(span)) - 2)
    while span / step > 12:
        step *= 2
    t = (s0 // step + 1) * step
    while t < s1:
        out.append(f'<line x1="{x(t):.1f}" y1="{ya - 3}" x2="{x(t):.1f}" y2="{ya + 3}" stroke="#9aa3ad"/>'
                   f'<text class="ax" x="{x(t):.1f}" y="{ya + 12}" text-anchor="middle">{t / 1000:g} kb</text>')
        t += step
    out.append(f'<text class="ax" x="{lm}" y="{H - 1}">{html.escape(locus["contig"])}</text>')
    out.append("</svg>")
    return "".join(out)


def build_sample(sample: str, organism: str, mp_gff: Path, fdb_gff: Path, out_html: Path, report_rel: str) -> list[dict]:
    loci, genes = load_matpredict(mp_gff)
    fdb = load_fungidb(fdb_gff)
    rows = []
    parts = [f"<h1>{html.escape(sample)} <span class='muted'>{html.escape(organism)}</span></h1>",
             f"<div class='nav'><a href='../../'>Release table</a><a href='{report_rel}'>MATPredict report</a>"
             f"<a href='detected_loci.gff3'>GFF3</a></div>",
             "<p>Top track: the genes MATPredict modelled at the MAT locus. Bottom track: every FungiDB gene in the "
             f"same region (locus span, shaded, plus {FLANK // 1000} kb each side). A FungiDB gene is green when a "
             "MATPredict gene overlaps it on the same strand. Hover a gene for its coordinates and product; FungiDB gene IDs link to fungidb.org.</p>",
             "<p class='legend'><span><i style='background:#b8471a'></i>core MAT gene</span><span><i style='background:#2a6fb0'></i>"
             "flanking conserved gene</span><span><i style='background:#6b7480'></i>other</span><span><i style='background:#1d7a46'></i>"
             "FungiDB gene matched by MATPredict</span></p>"]
    if not loci:
        parts.append("<p>MATPredict called no MAT locus in this genome.</p>")
    for L in loci:
        mp_genes = [g for g in genes if g["locus"] == L["id"]]
        region = [g for g in fdb.get(L["contig"], []) if g["end"] >= L["start"] - FLANK and g["start"] <= L["end"] + FLANK]
        matched = set()
        table = []
        for g in mp_genes:
            hits = [h for h in fdb.get(g["contig"], []) if overlap(g, h) > 0]
            v, why = verdict(g, hits)
            for h in hits:
                if h["strand"] == g["strand"]:
                    matched.add(h["id"])
            best = max(hits, key=lambda h: overlap(g, h)) if hits else None
            rows.append({"sample": sample, "organism": organism, "locus": L["id"], "idiomorph": L["idiomorph"],
                         "mp_gene": g["name"], "role": g["role"], "contig": g["contig"], "start": g["start"],
                         "end": g["end"], "strand": g["strand"], "verdict": v, "reason": why,
                         "fdb_gene": best["id"] if best else "", "fdb_name": best["name"] if best else "",
                         "fdb_description": best["desc"] if best else ""})
            table.append(f"<tr><td><b>{html.escape(g['name'])}</b></td><td>{html.escape(g['role'])}</td>"
                         f"<td class='mono'>{g['start']}-{g['end']} {g['strand']}</td>"
                         f"<td><span class='badge {VERDICT_CLASS[v]}'>{v}</span></td>"
                         f"<td>{gene_link(best)}</td>"
                         f"<td>{html.escape(best['desc'] if best else '')}</td><td class='muted'>{html.escape(why)}</td></tr>")
        extra = [g for g in region if g["id"] not in matched and L["start"] <= g["end"] and g["start"] <= L["end"]]
        parts.append(f"<h2>{html.escape(L['id'])}: {html.escape(L['idiomorph'] or 'undetermined')} "
                     f"<span class='muted'>{html.escape(L['family'])}, confidence {html.escape(L['confidence'])}, "
                     f"{html.escape(L['contig'])}:{L['start']}-{L['end']}</span></h2>")
        parts.append(locus_svg(L, mp_genes, region, matched))
        parts.append("<table><thead><tr><th>MATPredict gene</th><th>Role</th><th>Coordinates</th><th>Verdict</th>"
                     "<th>FungiDB gene</th><th>FungiDB product</th><th>Why</th></tr></thead><tbody>"
                     + "".join(table) + "</tbody></table>")
        if extra:
            parts.append("<p>FungiDB genes inside the locus span that MATPredict did not model:</p><table><thead><tr>"
                         "<th>FungiDB gene</th><th>Coordinates</th><th>Product</th></tr></thead><tbody>"
                         + "".join(f"<tr><td>{gene_link(g)}</td>"
                                   f"<td class='mono'>{g['start']}-{g['end']} {g['strand']}</td><td>{html.escape(g['desc'])}</td></tr>"
                                   for g in extra) + "</tbody></table>")
    out_html.parent.mkdir(parents=True, exist_ok=True)
    out_html.write_text(f"<!doctype html><html lang='en'><head><meta charset='utf-8'><meta name='viewport' content='width=device-width, initial-scale=1'>"
                        f"<title>{html.escape(sample)} – MAT locus vs FungiDB annotation</title><style>{CSS}</style></head>"
                        f"<body><main>{''.join(parts)}</main></body></html>")
    return rows


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--runs", required=True, type=Path)
    ap.add_argument("--gff-dir", required=True, type=Path)
    ap.add_argument("--release", default="71")
    ap.add_argument("--samples", required=True, type=Path)
    ap.add_argument("--organisms", type=Path, default=None, help="TSV with ASMID and organism columns (optional)")
    ap.add_argument("--out-dir", required=True, type=Path)
    ap.add_argument("--report-rel", default="index.html", help="relative link from the page to the MATPredict report")
    a = ap.parse_args()
    organisms = {}
    if a.organisms and a.organisms.is_file():
        organisms = {r["ASMID"]: r.get("organism", "") for r in csv.DictReader(open(a.organisms), delimiter="\t")}
    samples = [l.strip() for l in a.samples.read_text().splitlines() if l.strip() and not l.startswith("#")]
    all_rows, n = [], 0
    for s in samples:
        mp = a.runs / s / "detected_loci.gff3"
        fg = a.gff_dir / f"FungiDB-{a.release}_{s}.gff.gz"
        if not (mp.is_file() and fg.is_file()):
            continue
        all_rows += build_sample(s, organisms.get(s, ""), mp, fg, a.out_dir / f"{s}.html", a.report_rel)
        n += 1
    a.out_dir.mkdir(parents=True, exist_ok=True)
    cols = ["sample", "organism", "locus", "idiomorph", "mp_gene", "role", "contig", "start", "end", "strand",
            "verdict", "reason", "fdb_gene", "fdb_name", "fdb_description"]
    with open(a.out_dir / "annotation_verdicts.tsv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t"); w.writeheader(); w.writerows(all_rows)
    from collections import Counter
    print(f"wrote {n} annotation pages to {a.out_dir}; verdicts: {dict(Counter(r['verdict'] for r in all_rows))}", file=sys.stderr)


if __name__ == "__main__":
    main()
