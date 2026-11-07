# FungiDB MATpredict

[MATPredict](https://github.com/stajichlab/MATPredict) mating-type (MAT) locus calls on FungiDB genomes, compared with the FungiDB gene annotation. Results are kept per FungiDB release; release 71 is the primary set.

Generated 2026-11-06 with MATPredict v0.8.1, gmg v0.3.2, index bfd-reps-2026.10.5-s20000.gmgi, taxonomy 2026-10-01.

## Releases

| Release | Genomes run | MAT locus called | Not searched | Failed | MAT genes without annotation |
| --- | --- | --- | --- | --- | --- |
| release 71 | not run yet | | | | |
| [release 68](release-68/) | 1 | 1 | 0 | 0 | 0 |

## Method

Each genome is routed to a phylum by [GmG](https://github.com/stajichlab/GmG) (`matpredict detect --phylum auto`), then [MATPredict](https://github.com/stajichlab/MATPredict) locates the MAT locus and models its genes. Each modelled gene is matched to the FungiDB GFF by coordinate overlap on the same contig. Every step is a script in the [repository](https://github.com/hyphaltip/FungiDB_MATpredict); see its README to repeat the run.
