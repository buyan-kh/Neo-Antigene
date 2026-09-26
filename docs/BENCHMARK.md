# Benchmark

**Pooled null: two published cases, 17 scored positives, no method separation.**

The pipeline runs end to end, label-blind, on Ott et al. 2017 and on Bulik-Sullivan et al. 2018 patient CU04. On Ott a deterministic shuffle has the highest top-10 hit count, and Neo Antigene is one peptide ahead of an MHCflurry affinity sort. On CU04 the MHCflurry sort is ahead, 3 hits to 1. At k=20 the two methods tie on both cases. No within-case test rejects the hypothesis that the top 10 is a random draw from the assayed pool. That is the finding. It is why the next labels have to be generated rather than downloaded.

| | Ott 2017, patients 1–6 | Bulik-Sullivan 2018, CU04 | Sum of the two top-10s |
| --- | ---: | ---: | ---: |
| Assayed pairs the pipeline ranked | 139 | 19 | 158 |
| Positives in that pool | 14 | 3 | 17 |
| neoantigene hits in top 10 | 2 | 1 | 3 |
| `binding_only` hits in top 10 | 1 | 3 | 4 |
| `arbitrary` hits in top 10 | 3 | 1 | 4 |
| Chance expectation, hits in top 10 | 1.01 | 1.58 | 2.59 |

The last column adds two separate rankings. It is not one ranking of 158 peptides. Under chance the two top-10s together contain about 2.6 positives. Neo Antigene returned 3, binding-only 4, the shuffle 4.

## The null, case by case

Exact hypergeometric tail, P(at least this many positives in a random top 10):

| Case | Best point estimate | Hits | P |
| --- | --- | ---: | ---: |
| Ott | `arbitrary` | 3 of 14 | 0.064 |
| Ott | neoantigene | 2 of 14 | 0.265 |
| Ott | `binding_only` | 1 of 14 | 0.667 |
| CU04 | `binding_only` | 3 of 3 | 0.124 |
| CU04 | neoantigene | 1 of 3 | 0.913 |

One peptide moves Ott recall@10 by 7.1 points. On CU04, precision@10 cannot exceed 30% even if all three positives are ranked first, because k=10 and there are only three positives. k=20 does not add information on CU04: the assayed pool is 19, so every method's top 20 is the whole pool and the rate collapses to 3/19.

### Ott et al. 2017

Ott PA et al., *Nature* 2017;547:217–221. DOI [10.1038/nature22991](https://doi.org/10.1038/nature22991). Six vaccinated melanoma patients. Inputs are Supplementary Tables 2 and 4 (the somatic mutations and the HLA types). Labels are Supplementary Table 5, IFN-γ ELISPOT, opened only after the rankings were frozen. Case files and the per-patient table are in the Ott benchmark; the pooled retrieval numbers below are the ones that run reported.

| Method | recall@10 | recall@20 | precision@10 | Best positive rank |
| --- | ---: | ---: | ---: | ---: |
| neoantigene | 14.3% (2/14) | 14.3% (2/14) | 20% | 2 |
| `binding_only` | 7.1% (1/14) | 14.3% (2/14) | 10% | 9 |
| `presentation_only` | 7.1% (1/14) | 21.4% (3/14) | 10% | 8 |
| `expression_only` | 14.3% (2/14) | 21.4% (3/14) | 20% | 7 |
| `arbitrary` | 21.4% (3/14) | 28.6% (4/14) | 30% | 4 |

The shuffle beat both real methods. Neo Antigene's only advantage over `binding_only` is that one extra top-10 hit, and a best positive at rank 2 rather than rank 9. `expression_only` matching the ranker is expected: the paper's TPM is published for the vaccine peptides, so loading it would mark the assayed set, and it was left unloaded. Expression, tumor selectivity, and clonality are therefore constant. Those three weights are a large fraction of the model, and they did not run.

Four of 18 published positives were never scored: three are frameshift or neo-ORF peptides the generator does not build, and one (`CASP1` p.P172S) matched no Ensembl 116 isoform. They are absent from the 14.

This case's mutation input is the full table, 4,314 usable missense and in-frame variants out of 11,092 MAF rows. The score is still only the 139 pairs the trial assayed. Unassayed peptides are excluded, so a method cannot win by nominating something nobody tested, and it also cannot be credited for finding an epitope the trial never labeled.

### Bulik-Sullivan et al. 2018, patient CU04

Bulik-Sullivan B et al., *Nature Biotechnology*, published 17 December 2018. DOI [10.1038/nbt.4313](https://doi.org/10.1038/nbt.4313). One lung adenocarcinoma. Supplementary Data 4 is the peptide list; Supplementary Data 5 is the HLA type. Terminal output is in [`benchmarks/bulik-sullivan-2018-nsclc-cu04/RESULTS.md`](../benchmarks/bulik-sullivan-2018-nsclc-cu04/RESULTS.md). Row provenance is in [`SOURCES.md`](../benchmarks/bulik-sullivan-2018-nsclc-cu04/SOURCES.md).

`benchmark` on the default config, 26 September 2026:

```
neoantigene          k=10   validated 1/10 = 10.0%  (assayed pool 19, positives 3)  lift x0.33
binding_only         k=10   validated 3/10 = 30.0%  (assayed pool 19, positives 3)
neoantigene          k=19   validated 3/19 = 15.8%  (assayed pool 19, positives 3)  lift x1.00
binding_only         k=19   validated 3/19 = 15.8%  (assayed pool 19, positives 3)
```

`--k 20` clamped to 19. The same sort as `validation_rate_at_k` places the best positive at rank 2 for Neo Antigene (`EVADAATLTM`) and rank 3 for `binding_only` (`DENITTIQF`). Binding-only also places the other two positives at ranks 4 and 9. Neo Antigene places them at 13 and 15. `presentation_only`, `expression_only`, and `arbitrary` each returned 1/10.

Same hits, stated as recall so they can sit next to Ott:

| Method | recall@10 | recall@20 | precision@10 | Best positive rank |
| --- | ---: | ---: | ---: | ---: |
| neoantigene | 33.3% (1/3) | 100% (3/3) | 10% | 2 |
| `binding_only` | 100% (3/3) | 100% (3/3) | 30% | 3 |

Recall@20 is 100% for both because every positive sits inside a pool of 19.

This case is a re-rank of a list the authors already chose. Supplementary Data 5 records 511 nonsynonymous mutations and 336 expressed mutations in this tumor. The assayed list is two pools of 10. The restricting allele on each label is `full_ms_model_most_probable_restriction`, which the supplement's caption defines as the allele their model predicted. The Y/N is an IFN-γ ELISpot after in vitro expansion, not a tetramer and not an ex vivo call. No per-variant VAF, depth, or purity was published, so clonality is the same number on every variant. One frameshift negative was omitted because the generator cannot build it. The ranker is being scored on a shortlist whose hard selection already happened, with the clonality term switched off.

## These two runs did not share a config file

CU04 used `config/default.yaml` at main `4632e34`, SHA-256 `ba17954987d4bc45…`. Ott used the config on its own branch, SHA-256 `f36218f21aea1706…`, which sets `hydrophobicity` to 0, adds `self_dissimilarity`, `wt_dissimilarity`, and `mutation_exposure`, and the runner sets `output.include_failed = true` so a gated assayed peptide still receives a rank. Weights were priors in both runs. Neither run was refit on its labels. The hit counts above are what those two runs returned. They are not a paired experiment on one frozen weight vector.

## What was checked and not used

TESLA (Wells et al., *Cell* 2020) publishes the immunogenicity calls. The matched somatic variants are on Synapse (`syn21048999`) and in dbGaP (`phs000452`), behind a data-use agreement this project does not have. Sahin 2017, Hilf 2019, Keskin 2019, and Parkhurst 2019 do not, in the openly downloadable supplements, combine minimal 8–11-mers, a per-peptide experimental call, and genomic coordinates. IEDB was not used. The two public numbers people quote about it are different denominators, and neither is a T-cell assay label set:

- **55.8%** is eluted-ligand allele labels. Preibisch et al., bioRxiv, DOI [10.64898/2026.03.30.710191](https://doi.org/10.64898/2026.03.30.710191). Of 3,417,839 assessable mass-spectrometry eluted-ligand entries as of January 2025, 1,906,376 (55.8%) carry a predictor-dependent HLA assignment (Extended Data Table 1a). The peptides are experimental observations. The allele label is what the model supplied. Cite it as eluted-ligand allele labels.
- **16.4%** is an IEDB evidence-code count on `mhc_search`: "Inferred by motif or alleles present", 953,634 of 5,805,610 rows, measured 26 September 2026. A separate T-cell-table code, "MHC binding prediction", is 11.6% of `tcell_search`. That field does not record whether an allele was assigned by deconvolution before the row was submitted. It does not contradict 55.8%.

Both figures, and the curation rule that lets a prediction replace a coarser experimental restriction, are in [`docs/CITATIONS.md`](CITATIONS.md).

## What the owned labels have to be

A CRO pilot that could actually separate these methods:

- Allele-resolved. Tetramer, multimer, or a monoallelic measurement. A model's most probable restriction is not a label.
- Ex vivo, with the post-expansion result recorded separately if expansion is also done.
- Spot counts, or whatever the assay's native magnitude is, in `effect_size`.
- The full somatic mutanome as input, not the 20 peptides a previous model already kept.
- Per-variant VAF, depth, and purity, so clonality is computed rather than filled with a constant.
- More than three positives, and more than one patient. A doubling of precision@10 over `binding_only` at p < 0.05 takes on the order of 10² positives.

Fourteen positives in one melanoma trial and three in one lung-adenocarcinoma sample do not.
