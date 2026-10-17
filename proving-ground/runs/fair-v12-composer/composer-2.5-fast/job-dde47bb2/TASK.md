# Task

You are checking a pipeline that works out what a shifted reading frame
produces.

When a small insertion or deletion knocks a gene's reading frame out of
alignment, everything after that point is read differently and produces a new
stretch of protein that does not exist in a healthy copy. This pipeline
computes that new stretch for each variant and has already written its
results out.

Your job is to decide whether those results are correct. If they are not, say
what went wrong and how big the error is, list exactly which variants are
affected, and produce the corrected results for those variants. Compute the
corrected protein stretches yourself from the sequences and the annotation.
Do not reuse the numbers already in the output file. If the results are
correct, say so.

## Files you are given

| File | What it holds |
| --- | --- |
| `variants.vcf` | the variants. `CDSPOS` is the 1-based position within the coding sequence in `annotation.gff3` |
| `transcripts.fasta` | full transcript sequences |
| `annotation.gff3` | UTR and coding-sequence coordinates on the reference chromosome |
| `transcript_origin.tsv` | genomic origin of each transcript sequence |
| `enumerate_downstream.py` | the code that computes the new stretch |
| `cds_offsets.py` | the coding-sequence spans the code uses |
| `reported_neoorfs.csv` | the pipeline's current output |
| `README.md` | the team's own note |

## What to write out

Write one file to `/output/submission.json` with these fields:

```json
{
  "drift_detected": <true or false>,
  "mechanism": <one of the options below>,
  "magnitude": {"offset": <integer, in codons>, "direction": <"upstream" or "downstream" or "none">},
  "affected_variants": [<variant_id>, ...],
  "corrected": [
    {
      "variant_id": <variant_id>,
      "downstream_peptide": <amino-acid string>,
      "start_offset": <integer>
    }
  ]
}
```

`mechanism` must be exactly one of:

- `none` — the output is correct
- `annotation_release_mismatch` — the code and the annotation files come from different releases of the reference
- `zero_vs_one_based_position` — positions are treated as 0-based in one place and 1-based in another
- `downstream_start_skips_variant_codon` — the start of the new stretch is placed at the wrong codon
- `wrong_strand` — the sequence is read in the wrong orientation
- `wrong_codon_table` — the wrong translation table is used
- `indel_applied_at_wrong_sequence_level` — the variant is applied to the full transcript where the coding sequence was meant, or the reverse
- `translation_not_stopped_at_stop_codon` — translation runs past the first stop

Notes on the fields:

- `affected_variants` is only the variants whose reported result is wrong.
- `downstream_peptide` is the amino-acid sequence of the new stretch,
  including the stop character `*` if the stretch runs into one.
- `start_offset` is the zero-based codon index, within the coding sequence in
  `annotation.gff3`, at which the new stretch begins.
- `magnitude.offset` is how many codons separate the coding-sequence start the
  code used from the start in `annotation.gff3`. `direction` is where the
  annotation start sits relative to the start the code used. Use `none` and
  offset 0 when the two starts agree.
- Each corrected peptide is checked against an independent translation of the
  sequence, so it must be exactly right.
