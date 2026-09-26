# CU04 results

The pooled read of this case and Ott 2017 is in [`docs/BENCHMARK.md`](../../docs/BENCHMARK.md). Binding-only beats Neo Antigene on the only cut of this case that can separate them.

Commands were run from this repo on 26 September 2026, on the default config (`config/default.yaml`, presentation backend `mhcflurry`). Weights were not changed. Labels were not used to pick a threshold, a transcript, or a patient after seeing the scores. Patient CU04 was chosen before the run because Supplementary Data 4 gives that patient three individual-peptide positives, the most among patients whose peptides were tested one at a time.

## Commands

```
uv run neoantigene benchmark benchmarks/bulik-sullivan-2018-nsclc-cu04/sample.yaml benchmarks/bulik-sullivan-2018-nsclc-cu04/labels.tsv --k 10
uv run neoantigene benchmark benchmarks/bulik-sullivan-2018-nsclc-cu04/sample.yaml benchmarks/bulik-sullivan-2018-nsclc-cu04/labels.tsv --k 20
```

`benchmark` prints `ValidationSummary.describe()` for every built-in ranker. It does not print the rank of the best positive. That rank is computed below with the same sort `validation_rate_at_k` uses (`sorted` by descending score, then by key), on a separate pipeline run. The k=10 and k=20 lines from that run match the CLI.

## Terminal output

```
===== benchmark --k 10 =====
2026-09-26 10:16:29,506 INFO    [20260926T171551Z-baa769f2] neoantigene.pipeline: run 20260926T171551Z-baa769f2: sample CU04, 5 alleles, backend mhcflurry
2026-09-26 10:16:29,507 INFO    [20260926T171551Z-baa769f2] neoantigene.pipeline: read 17 protein-altering variants
2026-09-26 10:16:29,508 INFO    [20260926T171551Z-baa769f2] neoantigene.pipeline: somatic filtering: kept 17, dropped 0
2026-09-26 10:16:31,219 INFO    [20260926T171551Z-baa769f2] neoantigene.pipeline: loaded proteome with 382428 transcripts
2026-09-26 10:16:31,222 INFO    [20260926T171551Z-baa769f2] neoantigene.pipeline: generated 646 unique mutant peptides
2026-09-26 10:16:32,390 INFO    [20260926T171551Z-baa769f2] neoantigene.pipeline: predicted 6460 peptide-allele pairs
2026-09-26 10:16:44,401 INFO    [20260926T171551Z-baa769f2] neoantigene.pipeline: shortlisted 34 candidates
neoantigene          k=10   validated 1/10 = 10.0%  (assayed pool 19, positives 3)  lift x0.33
binding_only         k=10   validated 3/10 = 30.0%  (assayed pool 19, positives 3)
presentation_only    k=10   validated 1/10 = 10.0%  (assayed pool 19, positives 3)
expression_only      k=10   validated 1/10 = 10.0%  (assayed pool 19, positives 3)
arbitrary            k=10   validated 1/10 = 10.0%  (assayed pool 19, positives 3)
===== benchmark --k 20 =====
2026-09-26 10:16:49,873 INFO    [20260926T171645Z-22aac089] neoantigene.pipeline: run 20260926T171645Z-22aac089: sample CU04, 5 alleles, backend mhcflurry
2026-09-26 10:16:49,874 INFO    [20260926T171645Z-22aac089] neoantigene.pipeline: read 17 protein-altering variants
2026-09-26 10:16:49,875 INFO    [20260926T171645Z-22aac089] neoantigene.pipeline: somatic filtering: kept 17, dropped 0
2026-09-26 10:16:51,746 INFO    [20260926T171645Z-22aac089] neoantigene.pipeline: loaded proteome with 382428 transcripts
2026-09-26 10:16:51,753 INFO    [20260926T171645Z-22aac089] neoantigene.pipeline: generated 646 unique mutant peptides
2026-09-26 10:16:53,387 INFO    [20260926T171645Z-22aac089] neoantigene.pipeline: predicted 6460 peptide-allele pairs
2026-09-26 10:17:05,614 INFO    [20260926T171645Z-22aac089] neoantigene.pipeline: shortlisted 34 candidates
neoantigene          k=19   validated 3/19 = 15.8%  (assayed pool 19, positives 3)  lift x1.00
binding_only         k=19   validated 3/19 = 15.8%  (assayed pool 19, positives 3)
presentation_only    k=19   validated 3/19 = 15.8%  (assayed pool 19, positives 3)
expression_only      k=19   validated 3/19 = 15.8%  (assayed pool 19, positives 3)
arbitrary            k=19   validated 3/19 = 15.8%  (assayed pool 19, positives 3)
```

`--k 20` is clamped to the assayed pool. Every method then returns all 3 positives out of 19 peptides, which is the prevalence, not a ranking result.

## Numbers

| Method | Validation rate at 10 | Validation rate at 20 (clamped to 19) | Rank of the best positive, among the 19 assayed peptides |
| --- | --- | --- | --- |
| neoantigene | 1/10 = 10.0% | 3/19 = 15.8% | 2, `EVADAATLTM` / HLA-A*26:01 |
| binding_only | 3/10 = 30.0% | 3/19 = 15.8% | 3, `DENITTIQF` / HLA-B*18:01 |

`binding_only` places all three positives inside the top 10 (ranks 3, 4, and 9). Neo Antigene places one (`EVADAATLTM` at 2). The other two are at 13 (`DENITTIQF`) and 15 (`DTVEYPYTSF`).

Order among the 19 assayed peptides, same sort as `validation_rate_at_k`. Score for `binding_only` is the negated MHCflurry affinity percentile from this run, not the `mhcflurry_rank` column in the paper.

```
neoantigene
   1  neg  YHGDPMPCL  HLA-B*38:01  0.883166
   2  POS  EVADAATLTM  HLA-A*26:01  0.858262
   3  neg  VAKGFISRM  HLA-C*12:03  0.849803
   4  neg  MELKVESF  HLA-B*18:01  0.836585
   5  neg  LELKAVHAY  HLA-B*18:01  0.836118
   6  neg  EEADFLLAY  HLA-B*18:01  0.833407
   7  neg  FHATNPLNL  HLA-B*38:01  0.832419
   8  neg  DHFETIIKY  HLA-B*18:01  0.832162
   9  neg  EHIPESAGF  HLA-B*38:01  0.828968
  10  neg  VEIEQLTY  HLA-B*18:01  0.81879
  11  neg  DEERIPVL  HLA-B*18:01  0.818759
  12  neg  VEYPYTSF  HLA-B*18:01  0.806788
  13  POS  DENITTIQF  HLA-B*18:01  0.785337
  14  neg  QAVAAVQKL  HLA-C*12:03  0.774849
  15  POS  DTVEYPYTSF  HLA-A*26:01  0.767163
  16  neg  IEVEVNEI  HLA-B*18:01  0.628002
  17  neg  IQDQIQNCI  HLA-B*38:01  0.624424
  18  neg  ENITTIQFY  HLA-A*26:01  0.533891
  19  neg  VFKDLSVTL  HLA-B*38:01  0.365538

binding_only
   1  neg  EEADFLLAY  HLA-B*18:01  -0.001375
   2  neg  LELKAVHAY  HLA-B*18:01  -0.002125
   3  POS  DENITTIQF  HLA-B*18:01  -0.003875
   4  POS  EVADAATLTM  HLA-A*26:01  -0.0165
   5  neg  FHATNPLNL  HLA-B*38:01  -0.01725
   6  neg  EHIPESAGF  HLA-B*38:01  -0.02225
   7  neg  MELKVESF  HLA-B*18:01  -0.0285
   8  neg  VEIEQLTY  HLA-B*18:01  -0.0285
   9  POS  DTVEYPYTSF  HLA-A*26:01  -0.030875
  10  neg  VAKGFISRM  HLA-C*12:03  -0.04475
  11  neg  YHGDPMPCL  HLA-B*38:01  -0.045125
  12  neg  QAVAAVQKL  HLA-C*12:03  -0.058125
  13  neg  VEYPYTSF  HLA-B*18:01  -0.102375
  14  neg  DEERIPVL  HLA-B*18:01  -0.1105
  15  neg  IQDQIQNCI  HLA-B*38:01  -0.131
  16  neg  IEVEVNEI  HLA-B*18:01  -0.205625
  17  neg  ENITTIQFY  HLA-A*26:01  -0.222625
  18  neg  DHFETIIKY  HLA-B*18:01  -0.35475
  19  neg  VFKDLSVTL  HLA-B*38:01  -3.35112
```

All 19 labelled peptides were in the intersection. None were dropped before scoring.

`evaluate` reads a shortlist TSV, not `all_scored`. `rank` wrote 34 candidates (the default shortlist; 3,053 of 3,230 pairs were gated as weak binders and left out). Two assayed negatives did not make that shortlist, so the evaluate pool is 17, still with 3 positives.

```
uv run neoantigene rank benchmarks/bulik-sullivan-2018-nsclc-cu04/sample.yaml --out-dir /tmp/neo-bench/cu04-rank
uv run neoantigene evaluate benchmarks/bulik-sullivan-2018-nsclc-cu04/labels.tsv --ranked /tmp/neo-bench/cu04-rank/CU04.ranked.tsv --k 10
uv run neoantigene evaluate benchmarks/bulik-sullivan-2018-nsclc-cu04/labels.tsv --ranked /tmp/neo-bench/cu04-rank/CU04.ranked.tsv --k 20
```

```
run id             : 20260926T171802Z-e5ad7e5e
variants read      : 17
somatic filtering  : kept 17, dropped 0
peptides generated : 646
peptide-HLA pairs  : 3230
gated pairs        : {'weak_binder': 3053}
shortlisted        : 34 -> /tmp/neo-bench/cu04-rank/CU04.ranked.tsv
neoantigene          k=10   validated 1/10 = 10.0%  (assayed pool 17, positives 3)
neoantigene          k=17   validated 3/17 = 17.6%  (assayed pool 17, positives 3)
```

In that shortlist of 34, which also contains peptides the paper did not assay, `EVADAATLTM` is rank 2, `DENITTIQF` is rank 23, and `DTVEYPYTSF` is rank 26. `evaluate` has no binding-only comparator unless a second TSV is passed. The head-to-head is the `benchmark` output above.

## What this number is

This is a re-rank of a list the authors already chose. Supplementary Figure 11 says predicted neoantigens were combined into two pools of 10 peptides according to model ranking. Supplementary Data 5 says this tumor had 511 nonsynonymous mutations and 336 expressed mutations. The metric scores the 19 assayed peptides the generator can emit. It does not ask whether the model would have found these three among the 511.

The following are also true, and they all limit what the 10% vs 30% can mean:

- The allele on each label is `full_ms_model_most_probable_restriction`. The Supplementary Data 4 caption says that is the allele the authors' model predicted. The Y/N is an IFN-gamma ELISpot after in vitro expansion of PBMCs (Supplementary Figure 10), not an ex vivo or tetramer-defined restriction.
- No per-variant VAF, depth, or purity was published in Supplementary Data 4, so those columns are empty. Somatic filtering kept all 17. Clonality cannot separate peptides when every variant is missing the inputs clonality needs. The sample-level median VAF of 0.224 in Supplementary Data 5 was not copied onto the variants.
- Ensembl transcript ids are not in the paper. They were chosen by the rule in `SOURCES.md` so the published peptide is regenerated. TPM is gene-level, written onto that transcript.
- One frameshift negative (`QTKPASLLY`, BIRC6 G2619fs) was omitted because the generator does not build frameshifts. That removes one negative from both methods.
- Default weights are the priors in `config/default.yaml`. They were not refit on these labels. The presentation calls are MHCflurry 2.2.1 from this run. They are not the paper's EDGE ranks and not the paper's MHCflurry 1.2.0 ranks.
- n = 19 assayed peptides, 3 positives, 1 patient. That does not support a conclusion about the ranker.

## What a better label set would change

A useful owned set would be allele-resolved (tetramer or monoallelic), measured ex vivo as well as after expansion, with spot counts, on the full somatic mutanome rather than 20 preselected peptides, and with per-variant VAF, depth, and purity so clonality is not a constant. More than three positives, in more than one patient, is what would make a win or a loss mean something.
