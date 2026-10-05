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

- A model trained on the **HLA allele alone, with no peptide**, outscored all
  fifteen entrants in one published benchmark. Those benchmarks partly measure
  intra-HLA class imbalance.
- **PRIME** had roughly 70% of its training data inside the CEDAR benchmark it
  was evaluated on; performance dropped when the overlap was removed.
- **NetMHCpan-4.2** is fine-tuned on CEDAR neoepitopes, so CEDAR is not a
  clean test set for any method using recent NetMHCpan as a feature.
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
