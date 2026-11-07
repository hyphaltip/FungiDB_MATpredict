# FungiDB_MATpredict

Run [MATPredict](https://github.com/stajichlab/MATPredict) on FungiDB genomes
and compare the mating-type (MAT) locus genes it finds with the FungiDB gene
annotation. Every step is a script, so the whole analysis can be repeated
from a fresh clone with no other help. The task statement is in
[TASK.md](TASK.md). Results are published at
<https://hyphaltip.github.io/FungiDB_MATpredict/>.

## Requirements

- Linux x86-64. Tested on the UCR HPCC (SLURM, Rocky 8).
- [pixi](https://pixi.sh) on PATH (`module load pixi` on UCR HPCC).
- Internet access on the machine that runs steps 0 to 2 (a login node).
  The SLURM compute step runs offline.
- `curl`, `sha256sum`, `zcat`, `sbatch`.

## Steps

Run every command from the repository root.

| Step | Script | What it does | Where |
| --- | --- | --- | --- |
| 0 | `scripts/00_install_matpredict.sh` | Clone MATPredict at the pinned tag into `tools/MATPredict` and `pixi install` it. | login node |
| 1 | `scripts/01_install_gmg_and_taxonomy.sh` | Download the GmG ("Guess My Genome") binary and reference index, and build MATPredict's offline NCBI taxonomy table. | login node |
| 2 | `scripts/02_fetch_fungidb_data.sh` | For each sample in `samples.txt`, link or download the FungiDB genome, GFF and protein files into `data/`. | login node |
| 3 | `scripts/03_submit.sh` | Submit `scripts/03_run_matpredict.sbatch` as a SLURM array, one genome per task. Runs `matpredict detect --phylum auto`. Output: `results/matpredict/<sample>/`. | SLURM |
| 4 | `scripts/04_make_tables.sh` | Combine all runs into `results/tables/genomes.tsv` and `loci.tsv`, plus a summary. | login node |
| 5 | `scripts/05_compare.sh` | Compare MATPredict gene calls with the FungiDB GFF by coordinate overlap (`scripts/compare_matpredict_to_gff.py`). Output: `results/compare/`. | login node |
| 6 | `scripts/06_publish_docs.sh` | Copy per-genome HTML reports and all tables into `docs/` and rebuild `docs/index.md`. Commit and push `docs/` to publish. | login node |

```bash
bash scripts/00_install_matpredict.sh
bash scripts/01_install_gmg_and_taxonomy.sh
bash scripts/02_fetch_fungidb_data.sh
bash scripts/03_submit.sh          # wait for the array to finish (squeue -u $USER)
bash scripts/04_make_tables.sh
bash scripts/05_compare.sh
bash scripts/06_publish_docs.sh
git add docs && git commit -m "Publish results" && git push
```

## Settings

All paths, versions and SLURM resources are in `config.sh`. Pinned versions:
MATPredict, GmG, the GmG index and the NCBI taxonomy snapshot date. Change
`samples.txt` to add genomes: one FungiDB organism abbreviation per line,
as used in the FungiDB download tree
(`https://fungidb.org/common/downloads/release-68/<ABBREV>/`).

Environment variables you can set before a script:

| Variable | Default | Meaning |
| --- | --- | --- |
| `FUNGIDB_LOCAL_MIRROR` | `/bigdata/stajichlab/shared/FungiDB/release-68` | Local copy of the FungiDB download tree. Files found there are symlinked; others are downloaded. |
| `SLURM_PARTITION`, `SLURM_CPUS`, `SLURM_MEM`, `SLURM_TIME` | `epyc`, `8`, `16G`, `4:00:00` | Resources per genome. |

## How the comparison works

MATPredict writes `detected_loci.gff3` per genome: one `MAT_locus` feature
per locus and one `gene` feature per MAT gene it modelled. The comparison
script matches each called gene (`present=true`) to FungiDB genes on the same
contig whose interval overlaps it. Four tables come out:

- `annotated_MAT_genes.tsv`: MATPredict genes that overlap a FungiDB gene,
  with the FungiDB gene ID, name, product description, coordinates, overlap
  length and whether the strands agree.
- `unannotated_MAT_genes.tsv`: MATPredict genes with no overlapping FungiDB
  gene. These are candidate missing or mis-annotated genes. The nearest
  annotated gene on each side is listed with its distance.
- `locus_genes.tsv`: every FungiDB gene inside a MATPredict locus span.
- `summary.tsv`: counts per genome.

## Layout

```
config.sh        shared settings, sourced by every script
samples.txt      genomes to run
scripts/         the numbered steps above
data/            inputs (ignored by git; recreated by step 2)
tools/           installed software (ignored by git; recreated by steps 0 and 1)
results/         run output and tables (ignored by git)
docs/            GitHub Pages site, rebuilt by step 6
```

## License

GPLv3. See [LICENSE](LICENSE).
