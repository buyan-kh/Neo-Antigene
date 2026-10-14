# Merge two allele tables

Two groups measured the same set of tissue types and recorded a value for
each. The two tables were produced years apart by different teams, so the
names are not written the same way in both files.

Files:

```
source_a.csv     allele, payload_a
source_b.csv     allele, payload_b
gene_aliases.tsv canonical_symbol, prev_symbol, synonyms (pipe separated)
```

The alias table covers the older naming style. It does not cover every way a
name can be written.
