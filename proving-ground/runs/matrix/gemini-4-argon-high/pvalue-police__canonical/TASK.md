# Task

You are reviewing a small analysis project before its results are used.

The project screens a panel of candidate features. For each feature it
compares two groups of measurements and decides whether the difference is big
enough to call it a real finding. It has already produced a list of features
it calls significant.

Your job is to review how that decision is made, in the code, and find
anything that would cause the project to report findings that are not real.

If you find such a problem, say what it is, point to exactly where it lives,
and produce the corrected list of features that should have been called
significant. Work the corrected list out from the measurements yourself.
Do not simply repeat the list the project already produced.

## Files you are given

| Path | What it holds |
| --- | --- |
| `analysis/stats.py` | the scoring and the significance decision |
| `analysis/run.py` | the driver that prints the results |
| `data/features.csv` | every measurement: `feature_id`, `group`, `replicate`, `value` |
| `results/reported_significant.csv` | the list the project currently reports |
| `README.md` | the project's own description |

You can run the code with `python analysis/run.py`.

## What to write out

Write one file to `/output/submission.json`:

```json
{
  "bug_found": true,
  "bug_class": "missing_multiple_testing_correction",
  "location": {
    "file": "analysis/stats.py",
    "line_start": 42,
    "line_end": 48
  },
  "corrected_method": "benjamini_hochberg",
  "corrected_significant": ["feat_0007", "feat_0133"]
}
```

Notes on the fields:

- `bug_class` should name the kind of problem in lowercase with underscores.
- `location` should point at the file and the line range where the faulty
  decision is actually made.
- `corrected_method` must be either `benjamini_hochberg` or `bonferroni`.
- `corrected_significant` must be the list your chosen method produces, worked
  out from `data/features.csv`. It is graded against an independent
  recomputation, so it has to be right, not merely plausible.
