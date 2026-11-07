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
    by: dict[str, list[dict]] = {}
    with open_text(path) as fh:
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue
            c = line.rstrip("\n").split("\t")
            if len(c) < 9 or c[2] not in GENE_TYPES:
                continue
            a = dict(kv.split("=", 1) for kv in c[8].split(";") if "=" in kv)
            by.setdefault(c[0].split()[0], []).append({
                "id": unquote(a.get("ID", "")), "name": unquote(a.get("Name", "")),
                "start": int(c[3]), "end": int(c[4]), "strand": c[6], "desc": unquote(a.get("description", ""))})
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


def annotate(report: str, fdb: dict[str, list[dict]], contigs: list[str], url: str, annotation_link: str | None) -> tuple[str, int]:
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
    out = out.replace("</style>", ".fdb a{word-break:break-all}\n</style>", 1)
    return out, n_links


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--report", required=True, type=Path, help="MATPredict report.html")
    ap.add_argument("--gff", required=True, type=Path, help="FungiDB GFF3 for the same genome")
    ap.add_argument("--loci-gff3", type=Path, default=None, help="MATPredict detected_loci.gff3 (locus contig fallback)")
    ap.add_argument("--out", required=True, type=Path)
    ap.add_argument("--annotation-link", default=None, help="href of the annotation comparison page to add under the title")
    ap.add_argument("--fungidb-url", default=FUNGIDB_GENE_URL)
    a = ap.parse_args()
    fdb = load_fungidb(a.gff)
    contigs = locus_contigs(a.loci_gff3) if a.loci_gff3 else []
    out, n = annotate(a.report.read_text(), fdb, contigs, a.fungidb_url, a.annotation_link)
    a.out.parent.mkdir(parents=True, exist_ok=True)
    a.out.write_text(out)
    print(f"{a.out}: {n} FungiDB gene links", file=sys.stderr)


if __name__ == "__main__":
    main()
