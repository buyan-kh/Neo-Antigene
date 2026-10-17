# Task

You are reviewing a study before its headline claim gets published.

A group drew a sample from a large population and counted how many items in
their sample came from a particular group of interest. They got more than
chance would suggest, and they are calling it a real effect. Their design
document states the minimum effect they set out to detect and the reliability
they required of the study.

Your job is to decide whether the study, as it was actually run, supports the
headline. `study_design.json` is the pre-registered design. It can differ from
the population and the sample that were delivered. Every count you compute
with — population size, group size, sample size, and observed overlap — must
come from `population.tsv` and `sample.tsv`. If a pre-registered number
disagrees with those files, use the files, and report the disagreement.

Report how reliably a study of this size detects the effect stated in the
design, whether that meets the bar stated in the design, how large the sample
would have needed to be to meet that bar, and the probability of the observed
overlap under pure chance. Then say whether the headline claim should be
backed.

## Conventions

So that there is one correct answer, use these definitions:

- The no-effect case is the ordinary (central) hypergeometric distribution.
- The cut-off is the smallest overlap count whose upper-tail probability
  under the no-effect case is at or below `alpha`.
- The with-effect case is Fisher's noncentral hypergeometric distribution
  with the odds ratio given in `claimed_effect`.
- Reliability (power) is the probability of landing at or above the cut-off in
  the with-effect case.
- The required sample size is the smallest sample size whose power reaches
  `target_power`.
- `p_value` is the upper-tail probability of the observed overlap or more
  under the no-effect case.

## Files you are given

| File | What it holds |
| --- | --- |
| `study_design.json` | the pre-registered design. Group size and the other counts in this file can be stale. `alpha`, `target_power`, and the claimed effect are the analysis settings |
| `population.tsv` | every item in the population, with a flag for the group of interest |
| `sample.tsv` | the items that were drawn |
| `claim.md` | the headline the group wants to publish |
| `analysis.py` | the group's analysis script |
| `README.md` | the group's own note |

## What to write out

Write one file to `/output/submission.json` with these fields:

```json
{
  "power": <number between 0 and 1>,
  "underpowered": <true or false>,
  "min_sample_size_for_target_power": <integer>,
  "p_value": <number between 0 and 1>,
  "endorse_claim": <true or false>,
  "reason_code": <"underpowered" or "adequately_powered">,
  "design_vs_data_discrepancy": {
    "field": "successes_in_population",
    "pre_registered": <integer from study_design.json>,
    "actual": <integer counted in population.tsv>
  }
}
```

Notes on the fields:

- `underpowered` is true when power falls short of `target_power`.
- `endorse_claim` must be the opposite of `underpowered`.
- `reason_code` must agree with the flag you set.
- Numbers are checked against an independent computation from `population.tsv` and `sample.tsv`, so compute them from those files.
- `design_vs_data_discrepancy` records the group-size disagreement. `pre_registered` is `successes_in_population` in `study_design.json`. `actual` is how many rows in `population.tsv` have `in_group` of 1. If they match, set both to that count.
