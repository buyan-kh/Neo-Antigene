# Sources

Patient CU04 from Bulik-Sullivan et al., *Deep learning using tumor HLA peptide mass spectrometry datasets improves neoantigen identification*, Nature Biotechnology, published 17 December 2018.

- Article: https://doi.org/10.1038/nbt.4313
- Live article page: https://www.nature.com/articles/nbt.4313

Every variant row in `variants.tsv` and every peptide row in `labels.tsv` is copied from the files below. Nothing in those tables was filled in from memory.

## Files used

| What | Where it lives in this case | Source |
| --- | --- | --- |
| Peptide sequence, Y/N call, genomic mutation, gene, protein change, TPM, model-predicted restriction | `variants.tsv`, `expression.tsv`, `labels.tsv` | Supplementary Data 4, "Peptides tested for T-cell recognition in NSCLC patients". CSV downloaded 26 September 2026 from https://media.springernature.com/original/springer-static/esm/art%3A10.1038%2Fnbt.4313/MediaObjects/41587_2019_BFnbt4313_MOESM58_ESM.csv |
| HLA-A/B/C alleles, diagnosis | `sample.yaml` | Supplementary Data 5, "Demographics of NSCLC patients". XLSX downloaded the same day from https://media.springernature.com/original/springer-static/esm/art%3A10.1038%2Fnbt.4313/MediaObjects/41587_2019_BFnbt4313_MOESM59_ESM.xlsx |
| What the Y/N column means | `labels.tsv` notes | Supplementary Information PDF, Supplementary Data 4 caption and Supplementary Figure 10: https://media.springernature.com/original/springer-static/esm/art%3A10.1038%2Fnbt.4313/MediaObjects/41587_2019_BFnbt4313_MOESM46_ESM.pdf |
| Transcript protein sequences used to regenerate the peptides | `reference/Homo_sapiens.GRCh38.pep.all.fa.gz` | Ensembl release 116, https://ftp.ensembl.org/pub/release-116/fasta/homo_sapiens/pep/Homo_sapiens.GRCh38.pep.all.fa.gz . SHA-256 `9b43da92651b35814597af6a8b18f500b768679a49fa4678224f384917ce7668` |

The journal's 18 December 2018 correction says several supplementary-data numbers were swapped on first posting. The captions above are the ones on the live article page, which already includes that correction. Supplementary Data 4 on that page is the NSCLC peptide CSV.

## License

The article page has a Rights and permissions section and does not present a Creative Commons license. Copyright stays with the publisher. The supplementary CSV, XLSX, and PDF above downloaded without an account and without a data-use agreement. This directory stores the derived rows the pipeline needs (coordinates, sequences, calls, HLA, gene-level TPM). It does not include the original CSV or XLSX. The Ensembl proteome is the public FTP file, copied so the manifest path resolves.

## Sample

Supplementary Data 5, row `patient ID` = CU04:

- age range 61-70, Female, Hispanic or Latino, year of diagnosis 2013
- tumor stage I, primary tumor Lung, histological type Adenocarcinoma
- current anti-PD(L)-1 therapy: durvalumab plus tremelimumab
- HLA-A `A*24:26` and `A*26:01`, HLA-B `B*18:01` and `B*38:01`, HLA-C `C*12:03` and `C*12:03`

The homozygous `C*12:03` is listed once in `sample.yaml`. The same row also reports 336 expressed mutations, 511 nonsynonymous mutations, and a sample median VAF of 0.224. Those are sample-level counts. They were not copied onto individual variants, because the sheet does not give a VAF, depth, or copy number per mutation. `variants.tsv` leaves `dna_vaf`, `rna_vaf`, `tumor_depth`, `copy_number`, `population_af`, and `filter` empty. Known driver `TP53_R158G` is named in that row and is not in Supplementary Data 4, so it is not a variant here.

`cancer_type` is `lung adenocarcinoma`, from the histological type on that row.

## How a row got in

Supplementary Data 4 has 20 CU04 peptide rows. Kept a row only when `individual_peptide_response_any_timepoint`, after stripping whitespace, was `Y` or `N`. `Y` is `positive`. `N` is `negative`. The CSV stores the negatives as `N` plus a trailing space; the space was stripped and the letter was not changed.

The call is the paper's own any-timepoint flag. Supplementary Figure 10b shows CU04 individual-peptide IFN-gamma ELISpot spot-forming units at the initial visit and at 2 and 14 months, after in vitro expansion of PBMCs. This table does not re-threshold those plots. Spot counts are not in the CSV, so `effect_size` is empty. The authors' `mhcflurry_rank` and `full_ms_model_rank` columns are not copied into `predicted_rank` or `predicted_score`.

The restricting allele on each label is `full_ms_model_most_probable_restriction`. The Supplementary Data 4 caption says this column is the allele the authors' model predicted was most likely to present the peptide. It is not a tetramer or monoallelic measurement. The Y/N itself is the ELISpot call.

One row was left out: peptide `QTKPASLLY`, gene BIRC6, protein effect `G2619fs`, mutation `chr2_32487684_AG_A`, `mutation_type` `del_fs`, TPM 111.74, model allele HLA-A*26:01, individual call `N`. The peptide generator does not build frameshift peptides, so this negative is absent from both methods.

That leaves 19 peptides from 17 missense mutations. ADGRA3 C734F and CFL2 D66Y each contribute two peptides.

## Transcript chosen for each mutation

Supplementary Data 4 does not publish an Ensembl transcript id. The id in `variants.tsv` was chosen so the existing generator can rebuild the published peptide, and then checked against Ensembl VEP.

Rule, applied the same way to every mutation:

1. Parse `protein_effect` as a one-letter substitution (for example `C734F`).
2. Take every `protein_coding` translation of the published gene symbol in the Ensembl 116 proteome above. For EPRS, use HGNC symbol EPRS1, whose previous symbols include EPRS (https://rest.genenames.org/fetch/symbol/EPRS1, HGNC:3418). The `gene` column stays `EPRS`, the symbol printed in the paper.
3. Keep a translation only if the reference residue matches, `generate_for_variant` emits every published peptide for that mutation as an 8–11-mer overlapping the substituted residue, and Ensembl VEP (`https://rest.ensembl.org/vep/homo_sapiens/region/{chrom}:{pos}-{pos}:1/{alt}`, queried 26 September 2026) reports that transcript as a `missense_variant` with the same amino acids and protein position. VEP also had to show the published reference base on GRCh38 (`allele_string` starting with the published ref).
4. Among those translations, take the longest protein. Break remaining ties by the smallest transcript id.

All 17 mutations matched. `n` is how many translations passed the rule before the length tie-break.

| Gene | Change | Mutation field in Supplementary Data 4 | Transcript used | Gene TPM | n |
| --- | --- | --- | --- | --- | --- |
| STX5 | E134Q | `chr11_62827178_C_G` | ENST00000996436 | 83.43 | 15 |
| C1S | P295L | `chr12_7066530_C_T` | ENST00000328916 | 157.54 | 19 |
| CFL2 | D66Y | `chr14_34713369_C_A` | ENST00000298159 | 16.65 | 7 |
| NEK9 | D252H | `chr14_75117203_C_G` | ENST00001044055 | 20.29 | 12 |
| DHX58 | M513L | `chr17_42104792_T_A` | ENST00000251642 | 35.87 | 11 |
| EPRS | M277I | `chr1_220024376_C_G` | ENST00000366923 | 76.64 | 13 |
| INPP5B | Q606E | `chr1_37874128_G_C` | ENST00000373024 | 36.85 | 3 |
| ZFYVE9 | K845T | `chr1_52268541_A_C` | ENST00000287727 | 70.08 | 8 |
| ETAA1 | E493Q | `chr2_67404159_G_C` | ENST00000272342 | 38.47 | 3 |
| CAPG | E314K | `chr2_85395579_C_T` | ENST00000951335 | 151.69 | 28 |
| CRELD1 | Q347H | `chr3_9943508_G_C` | ENST00000383811 | 29.9 | 7 |
| ADGRA3 | C734F | `chr4_22413213_C_A` | ENST00001029712 | 20.67 | 8 |
| GCNT2 | P94L | `chr6_10556704_C_T` | ENST00000316170 | 25.19 | 4 |
| NUP205 | L691V | `chr7_135598004_C_G` | ENST00000921547 | 42.37 | 11 |
| ATP6V0A4 | P163H | `chr7_138762364_G_T` | ENST00000310018 | 47.21 | 14 |
| RNF216 | M45V | `chr7_5752914_T_C` | ENST00000389902 | 49.2 | 3 |
| ATP6AP2 | E145K | `chrX_40597563_G_A` | ENST00000901377 | 88.26 | 9 |

TPM is the gene-level `tpm` value from Supplementary Data 4, written onto the chosen `transcript_id` because the paper does not report transcript-level TPM. The expression reader looks up transcript id first.

## Every peptide

Call and allele are from Supplementary Data 4 columns `individual_peptide_response_any_timepoint` and `full_ms_model_most_probable_restriction`. Assay in `labels.tsv` is `ifng_elispot`.

| Peptide | Call | Allele on the label | Gene | Change | Mutation |
| --- | --- | --- | --- | --- | --- |
| `DEERIPVL` | negative | HLA-B*18:01 | RNF216 | M45V | `chr7_5752914_T_C` |
| `DENITTIQF` | positive | HLA-B*18:01 | ADGRA3 | C734F | `chr4_22413213_C_A` |
| `DHFETIIKY` | negative | HLA-B*18:01 | EPRS | M277I | `chr1_220024376_C_G` |
| `DTVEYPYTSF` | positive | HLA-A*26:01 | CFL2 | D66Y | `chr14_34713369_C_A` |
| `EEADFLLAY` | negative | HLA-B*18:01 | GCNT2 | P94L | `chr6_10556704_C_T` |
| `EHIPESAGF` | negative | HLA-B*38:01 | CRELD1 | Q347H | `chr3_9943508_G_C` |
| `ENITTIQFY` | negative | HLA-A*26:01 | ADGRA3 | C734F | `chr4_22413213_C_A` |
| `EVADAATLTM` | positive | HLA-A*26:01 | ZFYVE9 | K845T | `chr1_52268541_A_C` |
| `FHATNPLNL` | negative | HLA-B*38:01 | NEK9 | D252H | `chr14_75117203_C_G` |
| `IEVEVNEI` | negative | HLA-B*18:01 | NUP205 | L691V | `chr7_135598004_C_G` |
| `IQDQIQNCI` | negative | HLA-B*38:01 | ETAA1 | E493Q | `chr2_67404159_G_C` |
| `LELKAVHAY` | negative | HLA-B*18:01 | ATP6V0A4 | P163H | `chr7_138762364_G_T` |
| `MELKVESF` | negative | HLA-B*18:01 | INPP5B | Q606E | `chr1_37874128_G_C` |
| `QAVAAVQKL` | negative | HLA-C*12:03 | DHX58 | M513L | `chr17_42104792_T_A` |
| `VAKGFISRM` | negative | HLA-C*12:03 | CAPG | E314K | `chr2_85395579_C_T` |
| `VEIEQLTY` | negative | HLA-B*18:01 | STX5 | E134Q | `chr11_62827178_C_G` |
| `VEYPYTSF` | negative | HLA-B*18:01 | CFL2 | D66Y | `chr14_34713369_C_A` |
| `VFKDLSVTL` | negative | HLA-B*38:01 | ATP6AP2 | E145K | `chrX_40597563_G_A` |
| `YHGDPMPCL` | negative | HLA-B*38:01 | C1S | P295L | `chr12_7066530_C_T` |

Positives: `DENITTIQF` (ADGRA3 C734F, HLA-B*18:01), `EVADAATLTM` (ZFYVE9 K845T, HLA-A*26:01), `DTVEYPYTSF` (CFL2 D66Y, HLA-A*26:01).

## Why this paper, and the IEDB number

The label is an IFN-gamma ELISpot of in vitro expanded PBMCs, deconvolved to individual peptides (Supplementary Figure 10), not a bulk IEDB download.

A 2026 bioRxiv preprint (Deepflare; https://doi.org/10.64898/2026.03.30.710191, version 2 at https://www.biorxiv.org/content/10.64898/2026.03.30.710191v2) says that as of January 2025, 55.8% of assessable IEDB data are labeled by computational models. The results section makes the denominator explicit: an audit of mass-spectrometry eluted-ligand entries (3,971,663 records). Of 3,417,839 entries whose provenance could be determined, 1,906,376 (55.8%) carry predictor-dependent allele labels, and 1,511,463 (44.2% of the assessable set) were classed clean. Extended Data Table 1a is that split. The quotable form is eluted-ligand allele labels. It is a different denominator from IEDB evidence-code counts (about 16.4% of `mhc_search` rows tagged as inferred by motif or alleles present), and it is not a count of T-cell assay records. This case does not use IEDB. The same sentence is in `docs/BENCHMARK.md` and in `docs/CITATIONS.md`.

## Coordinates

The mutation field is `chr{chrom}_{pos}_{ref}_{alt}`. The `chr` prefix was removed. Chromosome, position, ref, and alt are otherwise the published string. Ensembl VEP on GRCh38 returned the published reference base, and the published amino-acid change, at all 17 positions.

