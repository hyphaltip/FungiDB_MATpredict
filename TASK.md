# Task: Compare MATPredict mating-type calls with FungiDB gene annotations

## Goal

Run MATPredict on a set of FungiDB genomes. Compare the MAT locus genes
MATPredict finds with the genes annotated in the matching FungiDB GFF files.
Report which annotated genes belong to the MAT locus and which MAT genes
MATPredict finds that have no annotation.

## Inputs

FungiDB files are kept per release under `data/release-<N>/`. Release 71 is
the primary set; release 68 is kept to see which genomes changed.

- `data/release-<N>/genome/`  : one genome assembly (FASTA) per species or strain.
- `data/release-<N>/protein/` : the FungiDB annotated proteins, used for the MATPredict DIAMOND fast path.
- `data/release-<N>/gff/`     : the FungiDB GFF3 annotation for each genome. The file
  names pair each GFF with its genome.

`scripts/02_fetch_fungidb_data.sh` fills these folders from the FungiDB
release named in `config.sh`, for every sample in `samples.txt`.

## Steps

Each step is a numbered script in `scripts/`; see README.md.

1. Install MATPredict (https://github.com/stajichlab/MATPredict) locally.
   Use `pixi` (the method the README recommends). Use `uv` only if pixi is
   not suitable. Record the installed version and commit.

2. For each genome, guess its taxonomic placement with GmG ("Guess My
   Genome", https://github.com/stajichlab/GmG). Use the MATPredict
   `detect --phylum auto` option, which runs GmG and routes the genome to
   the correct MAT gene families. Pin the GmG index version. If compute
   nodes have no internet access, fetch the index on a login node first.

3. Run `matpredict detect` on each genome. Keep the per-genome output
   (`detection_report.yaml`, `detected_loci.gff3`, `report.html`).

4. Build a combined report table with `matpredict report table`. The table
   must list, for each genome: the MAT idiomorph or allele call, each MAT
   gene found, its contig, start, end and strand, and the call confidence.

5. Compare the MATPredict calls with the FungiDB GFF annotation for the same
   genome. Match on genomic coordinates (overlap of the MATPredict gene
   interval with annotated gene or CDS features).

## Deliverables

1. **Annotated MAT genes**: a table of FungiDB genes that overlap a
   MATPredict MAT locus gene. Give the FungiDB gene ID, its product
   description, the MATPredict gene name (for example MAT1-1-1, MAT1-2-1),
   and the coordinates from both sources.

2. **Unannotated MAT genes**: a table of MAT genes MATPredict calls that
   overlap no gene in the FungiDB GFF. These are candidate missing or
   mis-annotated genes. Give the MATPredict gene name, coordinates,
   confidence and any nearby annotated genes.

3. A short summary per genome: MAT call, number of annotated MAT genes,
   number of unannotated MAT genes, and any genomes MATPredict could not
   route or search (`not_searched`).

## Publishing

Publish results on the GitHub Pages site (https://hyphaltip.github.io/FungiDB_MATpredict/),
which serves `docs/` on `main`, one folder per release:

- Copy each per-genome MATPredict `report.html` into `docs/release-<N>/reports/<genome>/`.
- Copy the combined tables (`genomes.tsv`, `loci.tsv`) and the comparison
  tables from the Deliverables section into `docs/release-<N>/tables/`.
- Link every report and table from `docs/release-<N>/index.md`.
- The top `docs/index.md` lists the releases, primary first, with summary
  counts and the genomes that changed between releases.
