# Merge two allele tables

Two groups typed the same panel of HLA alleles and recorded a value for each.
The two tables were produced years apart by different teams. An alias table
from the older team's records is included. It is not guaranteed to be current.

`payload_b` is stored as a percentage of the scale used for `payload_a`.
Divide it by 100 when merging.

Files:

```
source_a.csv     allele, payload_a
source_b.csv     allele, payload_b
gene_aliases.tsv canonical_symbol, prev_symbol, synonyms (pipe separated)
```
