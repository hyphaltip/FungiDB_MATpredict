---
title: Releases
---
# Releases

[MATPredict](https://github.com/stajichlab/MATPredict) mating-type (MAT) locus calls on FungiDB genomes, compared with the FungiDB gene annotation. Results are kept per FungiDB release; release 71 is the primary set. Generated 2026-11-06; the MATPredict version is the one recorded in each run's report (see a release page for the commit and database checksum).

## Releases

| Release | MATPredict version | Genomes run | MAT locus called | Not searched | Failed | MAT genes without annotation |
| --- | --- | --- | --- | --- | --- | --- |
| [release 71](release-71/) (primary) | v0.8.1 | 313 | 274 | 8 | 0 | 83 |
| [release 68](release-68/) | v0.8.1 | 331 | 261 | 41 | 0 | 83 |

## Changes between releases

From [release_changes.tsv](release_changes.tsv), built by `scripts/compare_releases.sh` (releases 71, 68).

| Status | Genomes | Samples |
| --- | --- | --- |
| changed_genome | 8 | AnigerCBS513-88, CneoformansH99, CneoformansH992018, CneoformansKN99, HcapsulatumH88, PaphanisDRCT72020, Spombe972h, ZtriticiIPO323 |
| dropped_after_68 | 42 | AastaciAPO3, Acandida2VRR, AeuteichesMF1, AinvadansNJM9701, AlaibachiiNc14, AmacrogynusATCC38327, BdendrobatidisJEL423, BlactucaeSF5, BsalamandrivoransBS, CgattiiCA1873, CgattiiIND107, CgattiiVGIIR265, … (42 total) |
| new_in_71 | 24 | AdauciA2016, AfumigatusA1160, BmaydisC4, BsorokinianaND90Pr, CbacillisporusCA1873, CbeticolaCb09-40, CdeneoformansB-3501A, CdeneoformansJEC21, CdeneoformansXL280, CdeuterogattiiR265, CkikuchiiMAFF305040, CtetragattiiIND107, … (24 total) |
| same_genome | 281 | AaculeatusATCC16872, AalternataSRC1lrK2f, AapisARSEF7405, AbrasiliensisCBS101740, AcampestrisIBT28561, AcarbonariusITEM5010, AclavatusNRRL1, AcristatusGZAAS20.1005, AellipticusCBS707.79, AeucalypticolaCBS122712, AfijiensisCBS313.89, AfischeriNRRL181, … (281 total) |

## Method

Each genome is routed to a phylum by [GmG](https://github.com/stajichlab/GmG) (`matpredict detect --phylum auto`), then [MATPredict](https://github.com/stajichlab/MATPredict) locates the MAT locus and models its genes. Each modelled gene is matched to the FungiDB GFF by coordinate overlap on the same contig. Every step is a script in the [repository](https://github.com/hyphaltip/FungiDB_MATpredict); see its README to repeat the run.
