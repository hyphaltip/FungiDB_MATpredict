#!/usr/bin/env python3
"""Add a "FungiDB gene" column to the Gene evidence tables of a MATPredict report.html.

For every locus card in the report, each gene row's position (absolute contig
coordinates, as MATPredict prints them) is intersected with the FungiDB GFF genes
on that locus's contig. Overlapping gene IDs are written as links to
https://fungidb.org/gene/<ID>, with the FungiDB product as a tooltip. A link to
the annotation comparison page is added under the report title.

The report is edited as text with regular expressions anchored on MATPredict's
own markup (`<section class="card locus" ...>`, `<table class="genes">`), so the
original file is left unchanged and the output is written to --out.
Standard library only; the GFF may be plain, .gz or .zst.
"""
from __future__ import annotations

import argparse
import gzip
import html
import io
import re
import subprocess
import sys
from pathlib import Path
from urllib.parse import unquote

GENE_TYPES = {"gene", "protein_coding_gene", "ncRNA_gene", "pseudogene"}
FUNGIDB_GENE_URL = "https://fungidb.org/gene/"


def open_text(path: Path):
    if path.suffix == ".gz":
        return gzip.open(path, "rt")
    if path.suffix == ".zst":
        return io.TextIOWrapper(subprocess.Popen(["zstd", "-dc", str(path)], stdout=subprocess.PIPE).stdout)
    return open(path)


def load_fungidb(path: Path) -> dict[str, list[dict]]:
    """Genes by contig. Each gene carries its transcripts' exons and CDS: gene["models"] = [{id, exons, cds}]."""
    by: dict[str, list[dict]] = {}
    genes_by_id: dict[str, dict] = {}
    tx_parent: dict[str, str] = {}      # transcript id -> gene id
    tx_parts: dict[str, dict] = {}      # transcript id -> {"exons": [...], "cds": [...]}
    with open_text(path) as fh:
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue
            c = line.rstrip("\n").split("\t")
            if len(c) < 9:
                continue
            a = dict(kv.split("=", 1) for kv in c[8].split(";") if "=" in kv)
            fid, parent = unquote(a.get("ID", "")), unquote(a.get("Parent", ""))
            start, end = int(c[3]), int(c[4])
            if c[2] in GENE_TYPES:
                g = {"id": fid, "name": unquote(a.get("Name", "")), "start": start, "end": end, "strand": c[6],
                     "desc": unquote(a.get("description", "")), "models": []}
                by.setdefault(c[0].split()[0], []).append(g)
                genes_by_id[fid] = g
            elif c[2] in ("exon", "CDS"):
                for t in parent.split(","):
                    tx_parts.setdefault(t, {"exons": [], "cds": []})["exons" if c[2] == "exon" else "cds"].append((start, end))
            elif parent:   # mRNA, transcript, pseudogenic_transcript, ncRNA ...
                tx_parent[fid] = parent
                tx_parts.setdefault(fid, {"exons": [], "cds": []})
    for tid, parts in tx_parts.items():
        g = genes_by_id.get(tx_parent.get(tid, ""))
        if g is not None:
            g["models"].append({"id": tid, "exons": sorted(parts["exons"]), "cds": sorted(parts["cds"])})
    return by


def locus_contigs(gff3: Path) -> list[str]:
    """Contigs of the MAT_locus features in detected_loci.gff3, in file order (locus-0, locus-1, ...)."""
    out = []
    if gff3.is_file():
        for line in open(gff3):
            c = line.rstrip("\n").split("\t")
            if len(c) >= 9 and c[2] == "MAT_locus":
                out.append(c[0])
    return out


POS_RE = re.compile(r'<td class="mono num">([\d,]+)\s*[–-]\s*([\d,]+)\s*<span class="strand">([+\-.])</span></td>')
SECTION_RE = re.compile(r'(<section class="card locus"[^>]*>)(.*?)(</section>)', re.S)
TABLE_RE = re.compile(r'(<table class="genes">.*?</table>)', re.S)


TICK_RE = re.compile(r'<line class="tick" x1="([\d.]+)"[^>]*/><text class="tick-label" x="[\d.]+"[^>]*>([\d.,]+)\s*(bp|kb|Mb)</text>')
AXIS_RE = re.compile(r'<line class="axis" x1="([\d.]+)" y1="\d+" x2="([\d.]+)"')
FIG_RE = re.compile(r'(<figure class="fig">.*?</figure>)', re.S)
TRACK_W = 960


def figure_range(fig: str) -> tuple[float, float, int, int] | None:
    """(x_left, x_right, bp_left, bp_right) of a MATPredict locus figure, from its axis line and two tick labels."""
    am = AXIS_RE.search(fig)
    ticks = [(float(x), float(v.replace(",", "")) * {"bp": 1, "kb": 1000, "Mb": 1_000_000}[u]) for x, v, u in TICK_RE.findall(fig)]
    if not am or len(ticks) < 2:
        return None
    (x0, b0), (x1, b1) = ticks[0], ticks[-1]
    if x1 == x0:
        return None
    bp_per_px = (b1 - b0) / (x1 - x0)
    xl, xr = float(am.group(1)), float(am.group(2))
    return xl, xr, int(round(b0 - (x0 - xl) * bp_per_px)), int(round(b0 + (xr - x0) * bp_per_px))


def fungidb_track(genes: list[dict], contig: str, xl: float, xr: float, lo: int, hi: int, matched: set[str], url: str,
                  release: str) -> str:
    """SVG of FungiDB gene models (exons as boxes, CDS filled, introns as lines) over the same x range."""
    span = max(1, hi - lo)
    x = lambda bp: xl + (max(lo, min(hi, bp)) - lo) / span * (xr - xl)
    region = sorted((g for g in genes if g["end"] >= lo and g["start"] <= hi), key=lambda g: g["start"])
    lane_h, label_h, top = 14, 14, 10
    lanes_end, lanes = [], []
    for g in region:
        label_px = len(g["id"] + (" " + g["name"] if g["name"] else "")) * 6.2
        ext0 = min(x(g["start"]), (x(g["start"]) + x(g["end"])) / 2 - label_px / 2)
        ext1 = max(x(g["end"]), (x(g["start"]) + x(g["end"])) / 2 + label_px / 2)
        for i, e in enumerate(lanes_end):
            if ext0 > e + 6:
                lanes_end[i] = ext1; lanes.append(i); break
        else:
            lanes_end.append(ext1); lanes.append(len(lanes_end) - 1)
    n = max(lanes, default=0) + 1
    height = top + n * (lane_h + label_h) + 14
    out = [f'<svg class="locus-fig fdb-track" viewBox="0 0 {TRACK_W} {height}" role="img" preserveAspectRatio="xMinYMin meet">'
           f'<title>FungiDB release {html.escape(release)} gene models on {html.escape(contig)}, {lo:,} to {hi:,}</title>']
    if not region:
        out.append(f'<text x="{xl}" y="{top + 12}" fill="#5f6873" font-size="11" font-family="Helvetica, Arial, sans-serif">'
                   f'no FungiDB gene annotated in this window</text>')
    for g, lane in zip(region, lanes):
        y = top + lane * (lane_h + label_h) + label_h
        color = "#1d7a46" if g["id"] in matched else "#6b7480"
        label = g["id"] + (f" {g['name']}" if g["name"] else "")
        title = html.escape(f"{g['id']} {g['start']:,}-{g['end']:,} ({g['strand']}): {g['desc']}", quote=True)
        gx0, gx1 = x(g["start"]), x(g["end"])
        out.append(f'<a href="{url}{html.escape(g["id"])}" target="_blank" rel="noopener"><title>{title}</title>')
        # gene line with a strand arrow head
        ym = y + lane_h / 2
        out.append(f'<line x1="{gx0:.1f}" y1="{ym:.1f}" x2="{gx1:.1f}" y2="{ym:.1f}" stroke="{color}" stroke-width="1.5"/>')
        tip = 5 if gx1 - gx0 > 8 else max(2, gx1 - gx0)
        if g["strand"] == "-":
            out.append(f'<polygon points="{gx0:.1f},{ym:.1f} {gx0 + tip:.1f},{ym - 4:.1f} {gx0 + tip:.1f},{ym + 4:.1f}" fill="{color}"/>')
        else:
            out.append(f'<polygon points="{gx1:.1f},{ym:.1f} {gx1 - tip:.1f},{ym - 4:.1f} {gx1 - tip:.1f},{ym + 4:.1f}" fill="{color}"/>')
        model = max(g["models"], key=lambda m: len(m["cds"]) + len(m["exons"]), default=None)
        if model:
            for (a, b) in model["exons"]:
                out.append(f'<rect x="{x(a):.1f}" y="{y + 3:.1f}" width="{max(1.0, x(b) - x(a)):.1f}" height="{lane_h - 6}" '
                           f'fill="{color}" opacity="0.35"/>')
            for (a, b) in model["cds"]:
                out.append(f'<rect x="{x(a):.1f}" y="{y:.1f}" width="{max(1.0, x(b) - x(a)):.1f}" height="{lane_h}" fill="{color}"/>')
        else:
            out.append(f'<rect x="{gx0:.1f}" y="{y + 3:.1f}" width="{max(1.0, gx1 - gx0):.1f}" height="{lane_h - 6}" fill="{color}" opacity="0.5"/>')
        out.append(f'<text x="{(gx0 + gx1) / 2:.1f}" y="{y - 3:.1f}" text-anchor="middle" fill="#424a54" font-size="11" '
                   f'font-family="Helvetica, Arial, sans-serif">{html.escape(label)}</text></a>')
    out.append("</svg>")
    return "".join(out)


def fungidb_cell(genes: list[dict], start: int, end: int, strand: str, url: str) -> str:
    hits = [g for g in genes if g["start"] <= end and g["end"] >= start]
    if not hits:
        return '<td class="small fdb">none</td>'
    parts = []
    for g in sorted(hits, key=lambda g: g["start"]):
        same = "" if g["strand"] == strand or strand == "." else " (opposite strand)"
        label = html.escape(g["id"]) + (f" {html.escape(g['name'])}" if g["name"] else "")
        title = html.escape(f"{g['start']}-{g['end']} {g['strand']}: {g['desc']}", quote=True)
        parts.append(f'<a href="{url}{html.escape(g["id"])}" title="{title}">{label}</a>{same}')
    return '<td class="small fdb">' + "<br>".join(parts) + "</td>"


def annotate(report: str, fdb: dict[str, list[dict]], contigs: list[str], url: str, annotation_link: str | None,
             release: str = "") -> tuple[str, int]:
    n_links = 0

    def do_section(m: re.Match, idx: list[int]) -> str:
        head, body, tail = m.group(1), m.group(2), m.group(3)
        i = idx[0]; idx[0] += 1
        # the contig: "on <contig>" in the card text, else the GFF3 locus order
        text = re.sub(r"<[^>]+>", " ", body)
        cm = re.search(r"\bon\s+(\S+?)[\s<,.]", text)
        contig = None
        if cm and cm.group(1) in fdb:
            contig = cm.group(1)
        elif i < len(contigs):
            contig = contigs[i]
        genes = fdb.get(contig or "", [])
        matched: set[str] = set()   # FungiDB genes a MATPredict gene overlaps on the same strand
        for pm in POS_RE.finditer(body):
            s_, e_ = int(pm.group(1).replace(",", "")), int(pm.group(2).replace(",", ""))
            for g in genes:
                if g["start"] <= e_ and g["end"] >= s_ and (g["strand"] == pm.group(3) or pm.group(3) == "."):
                    matched.add(g["id"])

        def do_figure(fm: re.Match) -> str:
            fig = fm.group(1)
            rng = figure_range(fig)
            if not rng or not contig:
                return fig
            xl, xr, lo, hi = rng
            track = fungidb_track(genes, contig, xl, xr, lo, hi, matched, url, release)
            return (fig + f'<figure class="fig fdb-fig"><div class="fig-scroll">{track}</div>'
                    f'<figcaption>FungiDB release {html.escape(release)} annotation of the same window: gene models from the GFF '
                    f'(filled = CDS, pale = UTR exon, line = intron, arrow = strand). Green genes overlap a MATPredict gene on the '
                    f'same strand; grey genes do not. Click a gene to open it on fungidb.org.</figcaption></figure>')

        body = FIG_RE.sub(do_figure, body)

        def do_table(tm: re.Match) -> str:
            nonlocal n_links
            t = tm.group(1)
            t = t.replace('<th scope="col">Reference record</th>',
                          '<th scope="col">Reference record</th><th scope="col">FungiDB gene</th>', 1)

            def do_row(rm: re.Match) -> str:
                nonlocal n_links
                row = rm.group(0)
                pm = POS_RE.search(row)
                if not pm:
                    return row.replace("</tr>", '<td class="small fdb">–</td></tr>')
                s, e = int(pm.group(1).replace(",", "")), int(pm.group(2).replace(",", ""))
                cell = fungidb_cell(genes, s, e, pm.group(3), url)
                n_links += cell.count("<a ")
                return row.replace("</tr>", cell + "</tr>")
            t = re.sub(r"<tr[^>]*>.*?</tr>", lambda rm: do_row(rm) if "<td" in rm.group(0) else rm.group(0), t, flags=re.S)
            return t
        body = TABLE_RE.sub(do_table, body)
        return head + body + tail

    idx = [0]
    out = SECTION_RE.sub(lambda m: do_section(m, idx), report)
    if annotation_link:
        out = re.sub(r"(</h1>)", r"\1" + f'<p class="meta"><a href="{annotation_link}">MAT locus above the FungiDB annotation '
                     f'(gene by gene)</a></p>', out, count=1)
    out = out.replace("</style>", ".fdb a{word-break:break-all}\n.fdb-fig{margin-top:-6px}\n.fdb-track a:hover text{fill:#1f5f99;text-decoration:underline}\n</style>", 1)
    return out, n_links


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--report", required=True, type=Path, help="MATPredict report.html")
    ap.add_argument("--gff", required=True, type=Path, help="FungiDB GFF3 for the same genome")
    ap.add_argument("--loci-gff3", type=Path, default=None, help="MATPredict detected_loci.gff3 (locus contig fallback)")
    ap.add_argument("--out", required=True, type=Path)
    ap.add_argument("--annotation-link", default=None, help="href of the annotation comparison page to add under the title")
    ap.add_argument("--fungidb-url", default=FUNGIDB_GENE_URL)
    ap.add_argument("--release", default="", help="FungiDB release number, for the track caption")
    a = ap.parse_args()
    fdb = load_fungidb(a.gff)
    contigs = locus_contigs(a.loci_gff3) if a.loci_gff3 else []
    out, n = annotate(a.report.read_text(), fdb, contigs, a.fungidb_url, a.annotation_link, a.release)
    a.out.parent.mkdir(parents=True, exist_ok=True)
    a.out.write_text(out)
    print(f"{a.out}: {n} FungiDB gene links", file=sys.stderr)


if __name__ == "__main__":
    main()
