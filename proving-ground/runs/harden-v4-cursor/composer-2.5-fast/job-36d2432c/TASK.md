# Task

You are merging two tables that describe the same things.

Two groups recorded a value for each of a panel of HLA alleles. The two tables
were written years apart by different teams. An alias table from the older
team's records is included.

Produce the merged table of alleles that appear in both files, each with its
value from each source. Report how many alleles matched, and how many of those
matches are written differently in the two files.

## Files you are given

| File | What it holds |
| --- | --- |
| `source_a.csv` | `allele`, `payload_a` |
| `source_b.csv` | `allele`, `payload_b` |
| `gene_aliases.tsv` | `canonical_symbol`, `prev_symbol`, `synonyms` (pipe separated) |
| `README.md` | the teams' own note |

## What to write out

Write one file to `/output/submission.json` with these fields:

```json
{
  "merged": [
    {"canonical_allele": <standard name>, "a_payload": <number>, "b_payload": <number>}
  ],
  "joined_count": <integer>,
  "naive_dropped_count": <integer>,
  "mapping": {<raw name as written in either file>: <standard name>}
}
```

Notes on the fields:

- `canonical_allele` uses the full standard two-field form, `HLA-<gene>*<dd>:<dd>`,
  keeping any trailing suffix letter the name carries (expression or group).
- `joined_count` is how many alleles appear in both files.
- `naive_dropped_count` is how many of those matched alleles are not written
  with the identical string in the two files.
- `mapping` must cover every raw name appearing in either file.
- Payloads must stay attached to the right source. Do not swap them.
- The merged table is checked against an independent reconciliation, so it
  has to be exactly right: no missing matches and no false matches.
