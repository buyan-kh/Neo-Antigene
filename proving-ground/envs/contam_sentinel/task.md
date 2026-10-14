# Task

You are reviewing a prediction tool before it gets published.

Someone built a tool that guesses which short protein fragments will stick to
a given immune-system receptor. They tested it on a set of examples they say
was kept separate from the data used to build it. The tool scored very well on
that test, and the team wants to publish.

Your job is to decide whether that test was genuinely kept separate, and
whether the reported score can be believed.

If any of the test examples also appear in the data the tool was built from,
list exactly which ones, say how they got there, and report an honest score
computed only from test examples that were genuinely unseen. If nothing is
wrong, say so.

## Files you are given

| File | What it holds |
| --- | --- |
| `predictor_training_set/` | the records the tool was built from, as the team received them |
| `benchmark.csv` | the supposedly separate test examples, with the right answer for each |
| `predictions.csv` | the tool's score for every test example |
| `README.md` | the team's own note about the project |

## What to write out

Write one file to `/output/submission.json` with these fields:

```json
{
  "contaminated": <true or false>,
  "leaked_pair_count": <integer>,
  "leaked_pairs": [[<peptide>, <allele>], ...],
  "mechanism": <one of the options below>,
  "trusted_metric": <true or false>,
  "corrected_metric": {
    "name": "auroc",
    "value": <number, or null>,
    "n_clean_pairs": <integer>,
    "uncomputable": <true or false>
  }
}
```

`mechanism` must be exactly one of:

- `none` — the test set is genuinely separate
- `train_test_overlap` — test examples are present in the build data
- `label_leakage` — the build data encodes the test answers without containing the test examples
- `duplicate_benchmark_rows` — the test set repeats its own examples
- `benchmark_distribution_shift` — the test set is separate but unrepresentative

Notes on the fields:

- `leaked_pairs` is every test example you found that also appears in the
  build data, written as `[peptide, allele]`. If the list is long, you may
  instead write it to `/output/leaked_pairs.csv` with one `peptide,allele`
  per line and leave `leaked_pairs` out of the JSON.
- `leaked_pair_count` must match the number of pairs you actually report.
- `corrected_metric.value` is the AUROC of the tool's scores against the true
  labels, using only the test examples that do not appear in the build data.
  `n_clean_pairs` is how many such examples there are. If none are left, set
  `value` to `null`, `uncomputable` to `true`, and `n_clean_pairs` to `0`.
  Do not invent a number.
- The leaked set and the corrected score are checked against an independent
  recomputation, so they have to be right, not merely plausible.
