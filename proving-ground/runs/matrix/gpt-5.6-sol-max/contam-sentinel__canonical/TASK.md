# Task

You are reviewing a prediction tool before it gets published.

Someone built a tool that guesses which short protein fragments will stick to
a given immune-system receptor. They tested it on a set of examples they say
was kept separate from the data used to build it. The tool scored very well on
that test, and the team wants to publish.

Your job is to decide whether that test was genuinely kept separate.

Work out whether any of the test examples also appear in the data the tool was
built from. If they do, say which ones and how they got there. Then say whether
the reported score can be believed, and report an honest score computed only
from test examples that were genuinely unseen.

Be careful: the same fragment and the same receptor can be written more than
one way across the two files. Comparing the text exactly as written will miss
matches that are really the same thing.

## Files you are given

| File | What it holds |
| --- | --- |
| `predictor_training_set.csv` | the records the tool was built from |
| `benchmark.csv` | the supposedly separate test examples, with the right answer for each |
| `predictions.csv` | the tool's score for every test example |
| `README.md` | the team's own note about the project |

## What to write out

Write one file to `/output/submission.json`:

```json
{
  "contaminated": true,
  "leaked_pair_count": 1486,
  "leaked_pairs": [["SIINFEKL", "HLA-A*02:01"]],
  "mechanism": "train_test_overlap",
  "trusted_metric": false,
  "corrected_metric": {
    "name": "auroc",
    "value": null,
    "n_clean_pairs": 0,
    "uncomputable": true
  }
}
```

Notes on the fields:

- `leaked_pairs` is every test example you found that also appears in the
  build data. If the list is long, you may instead write it to
  `/output/leaked_pairs.csv` with one `peptide,allele` per line and leave
  `leaked_pairs` out of the JSON.
- `leaked_pair_count` must match the number of pairs you actually report.
- `corrected_metric` is the score recomputed using only the untouched test
  examples. If none are left, set `value` to `null`, `uncomputable` to `true`,
  and `n_clean_pairs` to `0`. Do not invent a number.
