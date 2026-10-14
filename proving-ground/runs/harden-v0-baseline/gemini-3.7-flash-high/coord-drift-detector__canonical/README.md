# Frameshift downstream enumeration

This pipeline reads frameshift variants and works out the new protein tail
each one produces.

```
variants.vcf              the frameshift variants
transcripts.fasta         coding sequences for the affected transcripts
annotation.gff3           exon and CDS coordinates
enumerate_downstream.py   the enumeration code
reported_neoorfs.csv      the current output
```

The annotation files were recently upgraded to the latest release. The
enumeration code has been in use for a while and was not changed during the
upgrade.
