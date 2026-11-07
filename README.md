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

## Releases

Everything is keyed by FungiDB release. `config.sh` names the primary release
(`FUNGIDB_RELEASE=71`) and the list shown on the site (`FUNGIDB_RELEASES="71 68"`).
Data, sample lists, results, logs and site pages live in per-release folders:

```
data/release-<N>/{genome,gff,protein}/   inputs
samples/release-<N>.txt                  genomes to run
results/release-<N>/{matpredict,tables,compare}/
logs/release-<N>/
docs/release-<N>/{index.md,reports,tables}/
```

Every step runs on the primary release by default. Prefix a command with
`FUNGIDB_RELEASE=68` to run it on another release.

## Steps

Run every command from the repository root.

| Step | Script | What it does | Where |
| --- | --- | --- | --- |
| 0 | `scripts/00_install_matpredict.sh` | Clone MATPredict at the pinned tag into `tools/MATPredict` and `pixi install` it. | login node |
| 1 | `scripts/01_install_gmg_and_taxonomy.sh` | Download the GmG ("Guess My Genome") binary and reference index, and build MATPredict's offline NCBI taxonomy table. | login node |
| 2 | `scripts/02_fetch_fungidb_data.sh` | For each sample in `samples/release-<N>.txt`, link the FungiDB genome, GFF and protein files from a local mirror into `data/release-<N>/`. With `FUNGIDB_DOWNLOAD=1` it downloads missing files from fungidb.org. Entries FungiDB has no annotation for get a genome only. Skip this step when the files are already in `data/`. | login node |
| 3 | `scripts/03_submit.sh` | Submit `scripts/03_run_matpredict.sbatch` as a SLURM array. Each task runs `matpredict detect --phylum auto` on `SAMPLES_PER_TASK` genomes (default 25). Output: `results/matpredict/<sample>/`. | SLURM |
| 4 | `scripts/04_make_tables.sh` | Read the organism name from each genome FASTA header into `tables/sample_sheet.tsv`, then combine all runs into `genomes.tsv` and `loci.tsv` plus a summary. | login node |
| 5 | `scripts/05_compare.sh` | Compare MATPredict gene calls with the FungiDB GFF by coordinate overlap (`compare_matpredict_to_gff.py`), and write one page per genome with the MAT locus drawn above the FungiDB genes of the same region and a verdict per gene (`annotation_view.py`). Output: `results/release-<N>/compare/`. | login node |
| 6 | `scripts/06_publish_docs.sh` | Copy each MATPredict report into `docs/release-<N>/reports/<sample>/`, adding a FungiDB gene column with links to fungidb.org to its Gene evidence tables (`annotate_report.py`); copy the annotation pages and all tables; write the release index and the top `docs/index.md` (`build_docs_index.py`). Commit and push `docs/` to publish. | login node |
| - | `scripts/make_samples.sh` | Rebuild `samples/release-<N>.txt` from the genome files in `data/release-<N>/genome/`. | login node |
| - | `scripts/compare_releases.sh 68 71` | Which genomes are new, dropped, unchanged or changed between two releases (md5 of the uncompressed sequence). Writes `results/release_changes.tsv`, shown on the site. | login node |

```bash
bash scripts/00_install_matpredict.sh
bash scripts/01_install_gmg_and_taxonomy.sh
# put the FungiDB files in data/release-71/{genome,gff,protein}/ (or run step 2 with a mirror)
bash scripts/make_samples.sh       # samples/release-71.txt from data/release-71/genome/
bash scripts/03_submit.sh          # prints the array size and job id
squeue -u $USER                    # wait until the array is gone
bash scripts/04_make_tables.sh
bash scripts/05_compare.sh
bash scripts/compare_releases.sh 68 71   # optional, once both releases are in data/
bash scripts/06_publish_docs.sh
git add docs && git commit -m "Publish release 71 results" && git push
```

The same for release 68, with `FUNGIDB_RELEASE=68` in front of each command
from `make_samples.sh` onwards.

Steps 0 to 2 only need to be repeated when a version in `config.sh` or the
sample list changes. Step 3 skips every genome that already has a
`detection_report.yaml`, so re-running `scripts/03_submit.sh` after a failure
only runs what is missing. A failed genome is marked with
`results/matpredict/<sample>/failed_reason.txt`; its log is in `logs/`.

Two ways to fill `samples/release-<N>.txt`:

- Write the FungiDB abbreviations by hand, then step 2 links or downloads the files.
- Put the genome files in `data/release-<N>/genome/` first (copies or symlinks
  into a local FungiDB mirror), then `scripts/make_samples.sh` derives the list
  from the file names and reports which samples lack a GFF or protein file.
  A genome may be `_Genome.fasta.gz` or uncompressed `_Genome.fasta`.

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
| `FUNGIDB_RELEASE` | `71` | Release every step works on. |
| `FUNGIDB_DOWNLOAD` | `0` | Set to `1` to let step 2 download missing files from fungidb.org. |
| `SAMPLES_PER_TASK` | `25` | Genomes per SLURM array task. One genome takes about 1.5 minutes on 8 CPUs. |
| `SLURM_PARTITION`, `SLURM_CPUS`, `SLURM_MEM`, `SLURM_TIME` | `epyc`, `8`, `16G`, `4:00:00` | Resources per array task. Allow about 2 minutes per genome in the task. |

## The site

`docs/` is a Jekyll site served by GitHub Pages from `main`. The page shell
(header, CSS, the filter box and column sorting) is `docs/_layouts/default.html`;
`docs/_config.yml` applies it to every page. All page content under `docs/` is
generated by step 6 and is never edited by hand; delete `docs/release-<N>/` and
run step 6 again to rebuild it.

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
- `annotation_verdicts.tsv`: one row per MATPredict gene with a verdict against
  the FungiDB GFF: `annotated` (a FungiDB gene covers at least half of it on the
  same strand), `partial` (overlap with a different extent or strand), or
  `novel` (no overlapping FungiDB gene).

## Layout

```
config.sh        shared settings, sourced by every script
samples/         genomes to run, one file per release
scripts/         the numbered steps above
data/            inputs per release (ignored by git)
tools/           installed software (ignored by git; recreated by steps 0 and 1)
results/         run output and tables per release (ignored by git)
docs/            GitHub Pages site, one folder per release, rebuilt by step 6
```

## License

GPLv3. See [LICENSE](LICENSE).
