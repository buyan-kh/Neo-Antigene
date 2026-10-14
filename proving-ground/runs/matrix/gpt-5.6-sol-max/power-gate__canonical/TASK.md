# Task

You are reviewing a study before its headline claim gets published.

A group drew a sample from a large population and counted how many items in
their sample came from a particular group of interest. They got more than
chance would suggest, their test returned a small number, and they are calling
it a real effect.

A small test number only means the result is unlikely under pure chance. It
does not mean the study was large enough to reliably detect the effect they
claim to be looking for. A study can be too small to trust and still produce a
small test number by luck.

Your job is to work out whether this study was actually big enough. Compute
how reliably a study of this size would detect the effect they claim, decide
whether that is good enough against their own stated bar, and say how large
the sample would have needed to be to clear that bar. Then say whether the
headline claim should be backed.

## How "big enough" is defined here

This is pinned exactly so there is one correct answer:

- The no-effect case is the ordinary hypergeometric distribution.
- The cut-off is the smallest overlap count whose upper-tail probability under
  the no-effect case is at or below `alpha`.
- The with-effect case is **Fisher's noncentral hypergeometric distribution**
  using the odds ratio given in `claimed_effect`.
- Reliability is the chance of landing at or above the cut-off in the
  with-effect case.
- The required sample size is the smallest sample size that reaches
  `target_power`, found by binary search over sample size.

## Files you are given

| File | What it holds |
| --- | --- |
| `study_design.json` | population size, group size, sample size, observed overlap, `alpha`, `target_power`, claimed effect |
| `claim.md` | the headline the group wants to publish |
| `analysis.py` | the script that produced their test number |

## What to write out

Write one file to `/output/submission.json`:

```json
{
  "power": 0.2545,
  "underpowered": true,
  "min_sample_size_for_target_power": 247,
  "p_value": 0.013975,
  "endorse_claim": false,
  "reason_code": "underpowered"
}
```

Notes on the fields:

- `underpowered` is true when reliability falls short of `target_power`.
- `endorse_claim` must be the opposite of `underpowered`.
- `reason_code` is either `underpowered` or `adequately_powered`, and must
  agree with the flag you set.
- Do not assume a small study is automatically too small. Compute it.
