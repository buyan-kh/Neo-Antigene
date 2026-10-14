# Task

You are reviewing a study before its headline claim gets published.

A group drew a sample from a large population and counted how many items in
their sample came from a particular group of interest. They got more than
chance would suggest, and they are calling it a real effect. Their design
document states the minimum effect they set out to detect and the reliability
they required of the study.

Your job is to decide whether the study, as designed, supports the headline.
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
| `study_design.json` | population size, group size, sample size, observed overlap, `alpha`, `target_power`, claimed effect |
| `claim.md` | the headline the group wants to publish |
| `analysis.py` | the group's analysis script |

## What to write out

Write one file to `/output/submission.json` with these fields:

```json
{
  "power": <number between 0 and 1>,
  "underpowered": <true or false>,
  "min_sample_size_for_target_power": <integer>,
  "p_value": <number between 0 and 1>,
  "endorse_claim": <true or false>,
  "reason_code": <"underpowered" or "adequately_powered">
}
```

Notes on the fields:

- `underpowered` is true when power falls short of `target_power`.
- `endorse_claim` must be the opposite of `underpowered`.
- `reason_code` must agree with the flag you set.
- Numbers are checked against an independent computation, so compute them.
