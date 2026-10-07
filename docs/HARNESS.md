# The evaluation harness

A blinded, auditable evaluation of neoantigen rankings. It works on output
from any tool, not just this one.

```bash
uv run neoantigene compare labels.tsv \
  -m ours=ranked.tsv \
  -m netmhcpan=netmhcpan.tsv \
  --baseline netmhcpan --k 20 -o report.md
```

## Why this exists

Published neoantigen immunogenicity predictors cannot currently be compared,
and the reasons are specific rather than vague:

- A model given **only the HLA allele, with the peptide masked out**, outscored
  all fifteen entrants in a published comparison, because the benchmark rewards
  memorizing each allele's positive rate. Zhang et al., *Cell Genomics* 2026
  (ImmuBPI / ImmUni), DOI
  [10.1016/j.xgen.2026.101214](https://doi.org/10.1016/j.xgen.2026.101214).

  Three qualifications are mandatory whenever this is cited, because the
  unqualified version says something stronger and less true. The fifteen are
  predictor *configurations* from one comparison table, not fifteen
  independent tools, and three of them are that table's own authors' models.
  The benchmark is IEDB's **infectious-disease** immunogenicity set, not
  neoantigens and not CEDAR. And the whole comparison sits near chance — the
  strongest prior entrant scored AUROC 0.595 — so the finding is that the
  benchmark is weak, **not** that allele frequency predicts immunogenicity.
  It is also a single group with no independent replication, and the group's
  own debiasing method is the proposed fix, so cite the problem rather than
  the solution.
- **PRIME 1.0** had roughly 70% of its training peptides inside a
  CEDAR-derived evaluation set, and its performance dropped when that overlap
  was removed. Attribute this carefully: the measurement is third-party, by
  the IMPROVE authors (*Front Immunol* 2024), not a self-report, and it
  concerns version 1.0 rather than PRIME 2.0. PRIME's own authors used
  leave-one-study-out cross-validation precisely because "standard
  cross-validation results can be artificially boosted by batch effects", so
  the overlap is with a benchmark that postdates the model, not carelessness.
- **Released NetMHCpan-4.2** was fine-tuned on the whole CEDAR set, including
  the split the paper held out. All 1,486 `cedar_test` peptide–allele pairs
  are in the shipped `c00*_cedar` training files (5,172 records; Nilsson et
  al., *Front Immunol* 2025, DOI
  [10.3389/fimmu.2025.1616113](https://doi.org/10.3389/fimmu.2025.1616113)).
  CEDAR is not a clean test set for that neoepitope mode. It is the wrong
  conclusion for NetMHCpan-4.1, whose released training files have no CEDAR
  partition, and it is a weaker claim for 4.2's binding-affinity partition,
  which contains 706 of 5,027 CEDAR peptides as binding measurements.
- Contamination is **not** uniform, and assuming it is would be its own error.
  IMPROVE's in-house training set (17,520 peptide-HLA pairs, 467 positives)
  has **exactly zero** overlap with the CEDAR benchmark (2,436 pairs, 548
  positives), by peptide-plus-allele and by peptide alone. Some groups
  separate their data properly; the point of an audit is to find out which,
  not to assume.
- **IEDB's curation manual** instructs curators to record the most *precise*
  restriction rather than the most *direct* one, so a prediction can supersede
  a coarser experimental call. Prediction-derived labels are not confined to
  obviously low-confidence records.

Separately, the benchmarks in circulation are mostly too small to settle
anything. Both cases in [`BENCHMARK.md`](BENCHMARK.md) are underpowered, and
one of them — 19 assayed peptides with 3 positives — **cannot produce a
significant result at any outcome**, for any method. That fact is a property
of the pool, knowable before a single peptide is scored, and it is the first
thing this harness reports.

## What it reports, in order

1. **A verdict**, as the first line, phrased to be quotable without its
   surrounding context.
2. **Input digests.** SHA-256 of the label file and every ranking file, so a
   report is pinned to exact bytes.
3. **What the design could detect** — computed from the assayed pool alone,
   before any ranking is read. If no outcome reaches your alpha, it says so
   here.
4. **Retrieval** — recall@10, recall@20, precision@10, best positive rank.
5. **Against chance** — exact hypergeometric tail, not a normal
   approximation, which is poor in the tail.
6. **Ordering quality** — concordance over every positive-negative pair, which
   extracts far more from scarce labels than a top-k cut. Reported both pooled
   and stratified within patients.
7. **Paired comparison against the baseline** — a permutation test on the gap
   between two rankings over the peptides they both scored.
8. **Caveats**, generated from the data rather than written by hand.

## Input formats

**Labels** are the `assays.schema` TSV: `sample_id`, `peptide`, `allele`,
`assay`, `call`, with `call` one of `positive`, `negative`, `indeterminate`.
Indeterminate rows are counted and never scored — "not determined" is not a
negative. Columns `effect_size`, `replicate_count`, `run_date`, `operator`,
`notes` are optional and recommended; `notes` is the right place to record
which assay readout set the label.

**Rankings** are a TSV with `sample_id`, a peptide column (`mutant_peptide` or
`peptide`), optionally `allele`, and `score`, where higher is better. Any tool
that can emit four columns can be evaluated. Peptides are joined on
`sample_id|peptide|allele`, and `*` is accepted for an unknown allele.

## The rules it enforces

**Every method is restricted to the peptides they all scored.** A method
cannot win by having nominated an easier set.

**Peptides that were ranked but never assayed are excluded, not counted as
negatives.** Not testing a peptide is not evidence against it.

**Ties resolve against the ranking.** A tied block is ordered negatives-first,
so a method is never credited for a hit it did not actually rank above
anything. A ranking of all-equal scores reports zero hits, not luck.

**Patients are never pooled across.** `stratified auc` counts only
within-patient pairs. Concatenating patients into one ranking compares a
peptide from one tumor against a peptide from another and measures cohort
heterogeneity as ranking skill.

**Intervals come from a two-stage bootstrap** over patients and then peptides
within them. The plain cluster bootstrap covered 74–83% at three to six
patients while claiming 95%; the two-stage version measures 97–99%. Coverage
was measured, not assumed.

## Exit codes

| code | meaning |
| --- | --- |
| 0 | The report supports a claim: some method cleared chance, or separated from the baseline. |
| 2 | Nothing in the report could sustain a claim. |

Exit 2 exists so a script or CI job cannot mistake an inconclusive benchmark
for a passing one. An inconclusive result is not a failure — it is frequently
the correct scientific outcome of a small pool — but it must not be silently
consumed as a pass.

## How to read an underpowered result

If the capacity section says no outcome reaches significance, **do not quote a
hit rate from that report**, however good it looks. A perfect ranking on 19
peptides with 3 positives still yields p = 0.124 at k=10. The ordering
columns are the only ones with power at that size, and even they need
replication across patients.

Detecting a doubling of precision@20 over a presentation-sorted baseline takes
on the order of 100 positives, and — this is the counterintuitive part —
**power against chance does not improve with pool size** at fixed prevalence,
because top-k hits are approximately Binomial(k, prevalence). Assaying more
peptides from the same patients cannot rescue an underpowered benchmark. More
patients, a larger k, and a paired design can.

## Recording predictor provenance

If a ranking was produced by a predictor fine-tuned on the same labels, this
harness measures that overlap rather than the method, and it cannot detect
that for you. Record the predictor and its **weights version** alongside any
report. MHCflurry is trained on binding affinity and mass-spectrometry eluted
ligands rather than on immunogenicity labels, which makes it the cleaner
backend for evaluation against T-cell assay outcomes.
