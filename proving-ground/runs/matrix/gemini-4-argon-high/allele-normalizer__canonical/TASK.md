# Task

You are merging two tables that describe the same things.

Two groups recorded a value for each of a set of tissue types. The two tables
were written years apart by different teams, so the same tissue type is often
written differently in the two files. Some names use an older style that only
the provided alias table can resolve; others are the same name formatted a
different way.

Produce the correct merged table. Then report how many rows actually matched,
and how many of those matches a plain exact-text comparison would have missed
entirely.

Be careful. If you match the names exactly as they are written, you will
quietly lose rows that genuinely belong together, and the merged table will
look fine while being badly incomplete.

## Files you are given

| File | What it holds |
| --- | --- |
| `source_a.csv` | `allele`, `payload_a` |
| `source_b.csv` | `allele`, `payload_b` |
| `gene_aliases.tsv` | `canonical_symbol`, `prev_symbol`, `synonyms` (pipe separated) |
| `README.md` | the teams' own note |

The alias table covers the older naming style. It does **not** cover every
way a name can be formatted, so you will need to handle the rest yourself.

## What to write out

Write one file to `/output/submission.json`:

```json
{
  "merged": [
    {"canonical_allele": "HLA-A*02:01", "a_payload": 12.5, "b_payload": 88.1}
  ],
  "joined_count": 30,
  "naive_dropped_count": 24,
  "mapping": {"A*0201": "HLA-A*02:01", "A2": "HLA-A*02:01"}
}
```

Notes on the fields:

- `canonical_allele` should use the full standard form, for example
  `HLA-A*02:01`.
- `joined_count` is how many entries appear in both files once names are
  reconciled.
- `naive_dropped_count` is how many of those matches would be missed by
  comparing the raw text exactly as written.
- `mapping` must cover **every** raw name appearing in either file, pointing
  at the standard form it resolves to.
- Payloads must stay attached to the right source. Do not swap them.
