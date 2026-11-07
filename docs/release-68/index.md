# FungiDB release 68

[All releases](../) · Generated 2026-11-06 with MATPredict v0.8.1, gmg v0.3.2, index bfd-reps-2026.10.5-s20000.gmgi, taxonomy 2026-10-01, FungiDB release 68.

1 genomes run, 1 with a MAT locus called, 0 not searched, 0 failed.

## Results

Click a genome for MATPredict's full report.

| Genome | Report | Phylum | MAT idiomorph | Confidence | Genes found | Annotated in FungiDB | Not annotated | Locus coordinates |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| [AfumigatusA1163](reports/AfumigatusA1163/) | [HTML](reports/AfumigatusA1163/) · [GFF3](reports/AfumigatusA1163/detected_loci.gff3) · [YAML](reports/AfumigatusA1163/detection_report.yaml) | Ascomycota | MAT1-1 | high | MAT1-1-1, MAT1-2-1, APN2, SLA2, COX13 | 5 | 0 | scf_000003_A_fumigatus_A1163:2517390-2530649 |

## Tables

| File | Contents |
| --- | --- |
| [annotated_MAT_genes.tsv](tables/annotated_MAT_genes.tsv) | MATPredict genes that overlap a FungiDB gene, with the FungiDB ID and product |
| [genomes.tsv](tables/genomes.tsv) | one row per genome: routing, idiomorph, locus count |
| [loci.tsv](tables/loci.tsv) | one row per called locus: family, idiomorph, confidence, genes found and missing |
| [locus_genes.tsv](tables/locus_genes.tsv) | every FungiDB gene inside a MATPredict locus span |
| [summary.tsv](tables/summary.tsv) | comparison counts per genome |
| [summary_by_phylum.tsv](tables/summary_by_phylum.tsv) | MATPredict recovery counts by phylum |
| [unannotated_MAT_genes.tsv](tables/unannotated_MAT_genes.tsv) | MATPredict genes with no overlapping FungiDB gene, with nearest neighbours |

## Method

Each genome is routed to a phylum by [GmG](https://github.com/stajichlab/GmG) (`matpredict detect --phylum auto`), then [MATPredict](https://github.com/stajichlab/MATPredict) locates the MAT locus and models its genes. Each modelled gene is matched to the FungiDB GFF by coordinate overlap on the same contig. Every step is a script in the [repository](https://github.com/hyphaltip/FungiDB_MATpredict); see its README to repeat the run.
