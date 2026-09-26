# Benchmark provenance

Every peptide and label in this directory traces to the supplementary tables of
one paper. Regenerate with:

```bash
uv run python scripts/build_ott2017_benchmark.py
```

## Citation

Ott PA, Hu Z, Keskin DB, Shukla SA, Sun J, Bozym DJ, et al. "An immunogenic
personal neoantigen vaccine for patients with melanoma." *Nature*
2017;547(7662):217-221. doi:[10.1038/nature22991](https://doi.org/10.1038/nature22991). PMID 28678778.
PMC5577644.

## Files downloaded, with digests verified at build time

| Table | URL | SHA-256 |
| --- | --- | --- |
| Supplementary Table 2 | `https://media.springernature.com/original/springer-static/esm/art%3A10.1038%2Fnature22991/MediaObjects/41586_2017_BFnature22991_MOESM2_ESM.xlsx` | `bc537bca05393bad3c29afc04ea130cfe5b0d889dd9445985e3f3a4f6a83f00d` |
| Supplementary Table 4 | `https://media.springernature.com/original/springer-static/esm/art%3A10.1038%2Fnature22991/MediaObjects/41586_2017_BFnature22991_MOESM4_ESM.xlsx` | `c761b08f529ad21f1ff3c3a817a0187ee8913f5605687de4c38a893131e9b3b9` |
| Supplementary Table 5 | `https://media.springernature.com/original/springer-static/esm/art%3A10.1038%2Fnature22991/MediaObjects/41586_2017_BFnature22991_MOESM5_ESM.xlsx` | `9b10a6b6a430a331d88973889caed3cbe9328609774c936c7c2a851f41476b60` |

Table 2 is the complete somatic mutation list for Patients 1-10. Table 4 is
HLA typing. Table 5 is the immunizing peptides with IFN-gamma ELISPOT results.

## What went into the pipeline, and what did not

`inputs/` derives **only** from Tables 2 and 4. `validated.tsv` derives from
Table 5's ELISPOT columns and is read only after ranking is complete.

Table 5 also publishes each peptide's predicted MHC affinity and per-gene TPM.
Neither is used: the affinities are another predictor's output, and the TPM
values exist only for vaccine-selected peptides, so using them would advantage
precisely the peptides that were assayed.

## HLA typing used

| Patient | Class I alleles |
| --- | --- |
| Patient 1 | HLA-A*02:01 HLA-A*24:02 HLA-B*15:01 HLA-B*44:02 |
| Patient 2 | HLA-A*01:01 HLA-B*38:01 HLA-B*56:01 |
| Patient 3 | HLA-A*02:01 HLA-A*03:01 HLA-B*27:05 HLA-B*47:01 |
| Patient 4 | HLA-A*02:01 HLA-A*25:01 HLA-B*18:01 HLA-B*27:02 |
| Patient 5 | HLA-A*23:01 HLA-A*66:01 HLA-B*35:01 HLA-B*41:02 |
| Patient 6 | HLA-A*01:03 HLA-A*66:01 HLA-B*08:01 |

Supplementary Table 4 has HLA-A and HLA-B columns only, with no HLA-C, so
C-restricted candidates cannot be scored for this cohort. The allotypes are
stored as Excel time values and are decoded by total minutes; the decoded
values were cross-checked against the plain-text allele column in Table 5.

## Mutation accounting

| outcome | rows |
| --- | --- |
| dropped: no Ensembl isoform matches annotation | 377 |
| dropped: patient not vaccinated | 1563 |
| dropped: unparseable protein change | 5 |
| dropped: unsupported class 3'UTR | 773 |
| dropped: unsupported class 5'Flank | 49 |
| dropped: unsupported class 5'UTR | 437 |
| dropped: unsupported class Frame_Shift_Del | 29 |
| dropped: unsupported class Frame_Shift_Ins | 4 |
| dropped: unsupported class IGR | 30 |
| dropped: unsupported class Intron | 220 |
| dropped: unsupported class Non-coding_Transcript | 139 |
| dropped: unsupported class Nonsense_Mutation | 341 |
| dropped: unsupported class Silent | 2573 |
| dropped: unsupported class Splice_Site | 238 |
| kept | 4314 |
| maf rows read | 11092 |

Frameshift and nonsense variants are excluded because neo-ORF peptide
enumeration is not implemented in this package. That exclusion is visible in
the benchmark as assayed peptides the ranker never scored.

## Label accounting

| outcome | rows |
| --- | --- |
| dropped: ELISPOT not determined (n.d.) | 1 |
| dropped: no clean mutant peptide sequence | 6 |
| label: negative | 149 |
| label: positive | 18 |
| table 5 rows read | 174 |

The call comes from the **peptide-pulsed autologous APC** IFN-gamma ELISPOT.
Rows marked `n.d.` are dropped rather than treated as negative. The minigene
and autologous-tumor readouts are preserved in the `notes` column as secondary
evidence, not merged into the call.

Peptides that were never included in a vaccine pool have no ELISPOT result and
are therefore **unlabeled, not negative**. Metrics are computed only over the
assayed intersection.
