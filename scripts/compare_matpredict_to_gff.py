#!/usr/bin/env python3
"""Compare MATPredict MAT locus gene calls with a FungiDB GFF annotation.

For every sample, read MATPredict's `detected_loci.gff3` (one `MAT_locus`
feature per locus, one `gene` feature per MAT gene) and the FungiDB GFF3 for
the same genome. Match on genomic coordinates: a MATPredict gene and a FungiDB
gene are matched when they lie on the same contig and their intervals overlap.

Outputs (TSV, in --out-dir):
  annotated_MAT_genes.tsv    MATPredict genes that overlap a FungiDB gene
  unannotated_MAT_genes.tsv  MATPredict genes with no overlapping FungiDB gene
                             (plus the nearest annotated gene on each side)
  locus_genes.tsv            every FungiDB gene inside a MATPredict locus span
  summary.tsv                one row per sample

Only the standard library is used. Compressed GFF input (.gz, .zst) is read
directly.
"""
from __future__ import annotations

import argparse
import csv
import gzip
import io
import subprocess
import sys
from pathlib import Path
from urllib.parse import unquote

GENE_TYPES = {"gene", "protein_coding_gene", "ncRNA_gene", "pseudogene"}


def open_text(path: Path):
    """Open plain, .gz or .zst text transparently."""
    if path.suffix == ".gz":
        return gzip.open(path, "rt")
    if path.suffix == ".zst":
        proc = subprocess.Popen(["zstd", "-dc", str(path)], stdout=subprocess.PIPE)
        return io.TextIOWrapper(proc.stdout)
    return open(path)


def parse_attrs(field: str) -> dict:
    out = {}
    for part in field.strip().split(";"):
        if "=" in part:
            k, v = part.split("=", 1)
            out[k.strip()] = unquote(v.strip())
    return out


def read_gff(path: Path, types: set[str] | None):
    """Yield dicts for GFF3 rows, optionally filtered by feature type."""
    with open_text(path) as fh:
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 9:
                continue
            if types and cols[2] not in types:
                continue
            yield {
                "seqid": cols[0].split()[0], "type": cols[2],
                "start": int(cols[3]), "end": int(cols[4]), "strand": cols[6],
                "attrs": parse_attrs(cols[8]),
            }


def load_fungidb_genes(path: Path) -> dict[str, list[dict]]:
    """FungiDB genes grouped by contig, sorted by start."""
    by_contig: dict[str, list[dict]] = {}
    for f in read_gff(path, GENE_TYPES):
        a = f["attrs"]
        by_contig.setdefault(f["seqid"], []).append({
            "gene_id": a.get("ID", ""), "name": a.get("Name", ""),
            "type": f["type"], "start": f["start"], "end": f["end"], "strand": f["strand"],
            "description": a.get("description", ""),
        })
    for genes in by_contig.values():
        genes.sort(key=lambda g: g["start"])
    return by_contig


def load_matpredict(path: Path):
    loci, genes = [], []
    for f in read_gff(path, None):
        a = f["attrs"]
        if f["type"] == "MAT_locus":
            loci.append({"locus_id": a.get("ID", ""), "family": a.get("family", ""),
                         "idiomorph": a.get("idiomorph", ""), "confidence": a.get("confidence", ""),
                         "contig": f["seqid"], "start": f["start"], "end": f["end"]})
        elif f["type"] == "gene":
            genes.append({"mp_gene_id": a.get("ID", ""), "locus_id": a.get("Parent", ""),
                          "mp_gene": a.get("Name", ""), "role": a.get("role", ""),
                          "present": a.get("present", ""), "not_searchable": a.get("not_searchable", "false"),
                          "identity": a.get("identity", ""), "method": a.get("method", ""),
                          "status": a.get("status", ""),
                          "contig": f["seqid"], "start": f["start"], "end": f["end"], "strand": f["strand"]})
    return loci, genes


def overlapping(genes: list[dict], start: int, end: int) -> list[dict]:
    return [g for g in genes if g["start"] <= end and g["end"] >= start]


def nearest(genes: list[dict], start: int, end: int):
    """Nearest annotated gene upstream (ends before start) and downstream (starts after end)."""
    up = max((g for g in genes if g["end"] < start), key=lambda g: g["end"], default=None)
    down = min((g for g in genes if g["start"] > end), key=lambda g: g["start"], default=None)
    return up, down


def fmt_gene(g: dict | None, ref_start: int, ref_end: int) -> str:
    if not g:
        return ""
    dist = ref_start - g["end"] if g["end"] < ref_start else g["start"] - ref_end
    return f"{g['gene_id']}({g['start']}-{g['end']},{g['strand']},{dist}bp)"


def write_tsv(path: Path, rows: list[dict], columns: list[str]):
    with open(path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=columns, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)


def read_samples(path: Path) -> list[str]:
    out = []
    for line in path.read_text().splitlines():
        line = line.strip()
        if line and not line.startswith("#"):
            out.append(line)
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--runs", required=True, type=Path, help="directory of MATPredict detect run folders, one per sample")
    ap.add_argument("--gff-dir", required=True, type=Path, help="directory of FungiDB GFF files")
    ap.add_argument("--release", default="68", help="FungiDB release number used in file names")
    ap.add_argument("--samples", required=True, type=Path, help="samples.txt: one FungiDB abbreviation per line")
    ap.add_argument("--out-dir", required=True, type=Path)
    args = ap.parse_args()
    args.out_dir.mkdir(parents=True, exist_ok=True)

    annotated, unannotated, locus_genes, summary = [], [], [], []
    for sample in read_samples(args.samples):
        mp_gff = args.runs / sample / "detected_loci.gff3"
        fdb_gff = args.gff_dir / f"FungiDB-{args.release}_{sample}.gff.gz"
        row = {"sample": sample, "status": "", "n_loci": 0, "idiomorphs": "", "n_MAT_genes_called": 0,
               "n_annotated": 0, "n_unannotated": 0, "n_absent": 0, "n_genes_in_locus_span": 0}
        if not mp_gff.is_file():
            row["status"] = "no_matpredict_output"
            summary.append(row)
            continue
        if not fdb_gff.is_file():
            row["status"] = f"missing {fdb_gff.name}"
            summary.append(row)
            continue
        loci, genes = load_matpredict(mp_gff)
        fdb = load_fungidb_genes(fdb_gff)
        contigs_missing = {g["contig"] for g in genes} - set(fdb)
        row["status"] = "ok" if not contigs_missing else "contig_names_differ:" + ",".join(sorted(contigs_missing))
        row["n_loci"] = len(loci)
        row["idiomorphs"] = ";".join(sorted({l["idiomorph"] for l in loci if l["idiomorph"]}))
        locus_info = {l["locus_id"]: l for l in loci}

        for l in loci:
            for g in overlapping(fdb.get(l["contig"], []), l["start"], l["end"]):
                locus_genes.append({"sample": sample, **l, "fdb_gene_id": g["gene_id"], "fdb_name": g["name"],
                                    "fdb_type": g["type"], "fdb_start": g["start"], "fdb_end": g["end"],
                                    "fdb_strand": g["strand"], "fdb_description": g["description"]})
                row["n_genes_in_locus_span"] += 1

        for g in genes:
            if g["present"] != "true":
                row["n_absent"] += 1
                continue
            row["n_MAT_genes_called"] += 1
            l = locus_info.get(g["locus_id"].split(".seg")[0], {})
            base = {"sample": sample, "family": l.get("family", ""), "idiomorph": l.get("idiomorph", ""),
                    "locus_confidence": l.get("confidence", ""), **g}
            hits = overlapping(fdb.get(g["contig"], []), g["start"], g["end"])
            if hits:
                row["n_annotated"] += 1
                for h in hits:
                    ov = min(g["end"], h["end"]) - max(g["start"], h["start"]) + 1
                    annotated.append({**base, "fdb_gene_id": h["gene_id"], "fdb_name": h["name"],
                                      "fdb_type": h["type"], "fdb_start": h["start"], "fdb_end": h["end"],
                                      "fdb_strand": h["strand"], "overlap_bp": ov,
                                      "same_strand": str(h["strand"] == g["strand"]).lower(),
                                      "fdb_description": h["description"]})
            else:
                row["n_unannotated"] += 1
                up, down = nearest(fdb.get(g["contig"], []), g["start"], g["end"])
                unannotated.append({**base, "nearest_upstream": fmt_gene(up, g["start"], g["end"]),
                                    "nearest_downstream": fmt_gene(down, g["start"], g["end"])})
        summary.append(row)

    mp_cols = ["sample", "family", "idiomorph", "locus_confidence", "locus_id", "mp_gene_id", "mp_gene", "role",
               "contig", "start", "end", "strand", "identity", "method", "status"]
    write_tsv(args.out_dir / "annotated_MAT_genes.tsv", annotated,
              mp_cols + ["fdb_gene_id", "fdb_name", "fdb_type", "fdb_start", "fdb_end", "fdb_strand",
                         "overlap_bp", "same_strand", "fdb_description"])
    write_tsv(args.out_dir / "unannotated_MAT_genes.tsv", unannotated,
              mp_cols + ["nearest_upstream", "nearest_downstream"])
    write_tsv(args.out_dir / "locus_genes.tsv", locus_genes,
              ["sample", "locus_id", "family", "idiomorph", "confidence", "contig", "start", "end",
               "fdb_gene_id", "fdb_name", "fdb_type", "fdb_start", "fdb_end", "fdb_strand", "fdb_description"])
    write_tsv(args.out_dir / "summary.tsv", summary, list(summary[0].keys()) if summary else ["sample"])
    print(f"wrote 4 tables to {args.out_dir}: {len(annotated)} annotated, {len(unannotated)} unannotated "
          f"MAT gene calls, {len(locus_genes)} genes inside locus spans, {len(summary)} samples", file=sys.stderr)


if __name__ == "__main__":
    main()
