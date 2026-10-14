# Frameshift downstream enumeration

This pipeline reads frameshift variants and works out the new protein tail
each one produces.

```
variants.vcf              the frameshift variants (transcript space)
transcripts.fasta         full transcript sequences, soft-masked
annotation.gff3           UTR and CDS coordinates per transcript
enumerate_downstream.py   the enumeration code
cds_offsets.py            cached CDS spans used by the code
reported_neoorfs.csv      the current output
```

Run with `python enumerate_downstream.py transcripts.fasta variants.vcf out.csv`.

The sequence and annotation files were refreshed to the current Ensembl
release last quarter.
