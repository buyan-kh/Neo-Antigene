# Task

You are reviewing a small analysis project before its results are used.

The project screens a panel of candidate features. For each feature it
compares two groups of measurements and decides whether the difference is big
enough to call it a real finding. It has already produced a list of features
it calls significant.

Your job is to review how that decision is made, in the code, and decide
whether the reported list can be trusted. If there is a flaw, classify it,
point to exactly where it lives, and produce the list of features that should
have been called significant. Work that list out from the measurements
yourself, keeping the project's test statistic. If there is no flaw, say so
and return the project's list.

## Files you are given

| Path | What it holds |
| --- | --- |
| `bin/select_candidates.py` | builds the analysis table from the raw export |
| `analysis/stats.py` | the scoring and the significance decision |
| `analysis/run.py` | the driver that prints the results |
| `raw/run_2024-11-02_plate_reader.csv` | the instrument export |
| `data/features.csv` | the table the screen is run on |
| `results/reported_significant.csv` | the list the project currently reports |
| `README.md` | the project's own description |

You can rebuild the analysis table with `python bin/select_candidates.py` and
run the screen with `python analysis/run.py`.

## What to write out

Write one file to `/output/submission.json` with these fields:

```json
{
  "bug_found": <true or false>,
  "bug_class": <one of the options below>,
  "location": {
    "file": <path>,
    "line_start": <integer>,
    "line_end": <integer>
  },
  "corrected_method": <"benjamini_hochberg" or "bonferroni">,
  "corrected_significant": [<feature_id>, ...]
}
```

`bug_class` must be exactly one of the following. Pick the single option that
most precisely describes the flaw:

- `none` — the decision is sound
- `no_multiple_testing_correction` — raw p-values are compared to the threshold with no correction anywhere
- `correction_applied_after_prefiltering` — a correction exists but is applied only to features that already passed an uncorrected threshold
- `wrong_test_statistic` — the test used does not match the data or the stated assumptions
- `one_sided_p_value_misreported` — a one-sided p-value is treated as two-sided or vice versa
- `groups_mislabelled` — the two groups are compared the wrong way round
- `boundary_comparison_error` — a `<` versus `<=` mistake at the threshold

Notes on the fields:

- `location` should point at the file and the line range where the faulty
  decision is actually made.
- `corrected_method` is the correction you apply when producing the list.
  Either accepted method is fine, applied to all features.
- `corrected_significant` must be the list your chosen method produces,
  applied to every screen analyte at a 5% level. Spike-in controls are not
  hypotheses. It is graded against an independent recomputation, so it has
  to be right, not merely plausible.
