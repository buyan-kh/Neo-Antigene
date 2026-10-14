# Task

You are checking a pipeline that works out what a shifted reading frame
produces.

When a small insertion or deletion knocks a gene's reading frame out of
alignment, everything after that point is read differently and produces a new
stretch of protein that does not exist in a healthy copy. This pipeline
computes that new stretch for each variant and has already written its
results out.

The annotation files and the enumeration code were not produced against the
same version of the reference. Where the versions disagree, the starting
point of the new stretch is placed wrong, and every result from that point on
is wrong with it.

Work out whether that has happened here. If it has, say what the mismatch is
and how big the shift is, list which variants are affected, and produce the
corrected results. Compute the corrected protein stretches yourself from the
sequences and the annotation. Do not reuse the numbers already in the output
file.

## Files you are given

| File | What it holds |
| --- | --- |
| `variants.vcf` | the variants, with transcript and position in the coding sequence |
| `transcripts.fasta` | coding sequence for each transcript |
| `annotation.gff3` | coordinates, with the release it came from in the header |
| `enumerate_downstream.py` | the code that computes the new stretch |
| `reported_neoorfs.csv` | the pipeline's current output |
| `README.md` | the team's own note |

## What to write out

Write one file to `/output/submission.json`:

```json
{
  "drift_detected": true,
  "mechanism": "ensembl_release_coordinate_offset",
  "magnitude": {"offset": 1, "direction": "downstream"},
  "affected_variants": ["ENST00000000001:142"],
  "corrected": [
    {
      "variant_id": "ENST00000000001:142",
      "downstream_peptide": "MKT*",
      "start_offset": 47
    }
  ]
}
```

Notes on the fields:

- `mechanism` should name the kind of mismatch in lowercase with underscores.
- `magnitude.offset` is how many codons the start is out by.
- `downstream_peptide` is the amino-acid sequence of the new stretch,
  including the stop character `*` if the stretch runs into one.
- `start_offset` is the codon index where the new stretch begins, counting
  from zero at the start of the coding sequence.
- Each corrected peptide is checked against an independent translation of the
  sequence, so it must be exactly right.
