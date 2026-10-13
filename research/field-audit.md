# Field audit

This is the multi-agent research. A fifteen-agent pass read primary files,
not press summaries. Agents were read-only. When an agent's reading of a
paper disagreed with the released code or data, the file won. Every number
below was checked against a file. The line references and the exact
commands live in `docs/AUDIT.md`. This page is the record of what those
checks concluded.

Three of the conclusions correct something this repository had said before.
They are marked.

## IMPROVE kept its training set clean

*Frontiers in Immunology* 2024, DOI 10.3389/fimmu.2024.1360281.
Repository `github.com/SRHgroup/IMPROVE_paper`.

Two released label tables:

| Table | Rows | Yeses | Rate |
| --- | ---: | ---: | ---: |
| In-house, used to train | 17,520 | 467 | 2.67% |
| CEDAR benchmark | 2,436 | 548 | 22.50% |

Counted overlap: zero peptides in common, and zero peptide-plus-allele pairs
in common. Including zero of the 548 CEDAR yeses. Contamination is not
universal. This group separated its data. An audit that assumes every paper
leaked would be wrong about this one.

The tables have no patient id. Row-wise cross-validation on data like this
inflates AUC by about 0.14 even when there is no real signal. Patient-grouped
splits are impossible here. The table can show the sign of a feature. It
cannot set a weight you would spend synthesis money on.

The two sets also have very different yes rates (about eightfold). CEDAR
yeses are peptides somebody chose to publish. The in-house set is closer to
a screen. A precision number from one is not a precision number from the
other.

Their feature ranking, read from the released `results.zip`, puts a
hydrophobicity-of-the-core feature first and their agretopicity feature 13th
of 22. Their agretopicity is a ratio of eluted-ligand ranks, not a ratio of
binding affinities. That comparison has already been run, at 467 yeses, and
the rank ratio was mid-table. Reformulating our agretopicity as a rank ratio
is not supported by their result.

Their hydrophobicity feature is also not our feature. It is computed on
NetMHCpan's predicted binding core, then bulged residues are added back.
Their own code defines the core two different ways, one of which matches our
slice of the raw peptide. Dropping position 3 from our feature would not
make the two features the same. That recommendation is withdrawn.

A separate code check: several of their published p-values were passed to a
statistical function that silently ignores the multiple-testing argument.
Treat those p-values as unadjusted. The gap between their top feature and
their agretopicity feature is large enough that this does not reorder them.

The IMPROVE repository has no license file. Do not build on their code or
their tables until that is resolved with the authors.

## PRIME 1.0 overlapped a later benchmark

The repeated figure is that about 70% of PRIME's training peptides sit
inside a CEDAR-derived evaluation. That figure is real, with three limits:

- It is PRIME 1.0, not PRIME 2.0. No published overlap number for 2.0 was
  found.
- The measurement was made by the IMPROVE authors, not by PRIME's authors.
- CEDAR was published after PRIME 1.0. This is not "they trained on the
  test they later claimed." PRIME's own authors already refused ordinary
  cross-validation because peptides from the same study in both halves
  inflate the score. The same reason we group folds by patient.

## A model that ignores the peptide can win a weak benchmark

Zhang and colleagues, *Cell Genomics* 2026. DOI 10.1016/j.xgen.2026.101214.
A model given only the HLA allele, with the peptide masked, beat the other
entries in one comparison.

The limits, which have to travel with the citation:

- The other entries are configurations from one table, not fifteen separate
  labs. Three belong to the table's own authors.
- The benchmark is infectious-disease data from IEDB. It is not neoantigens
  and it is not CEDAR. On the neoepitope benchmark, the HLA-only version was
  only the best of the authors' own three variants.
- The whole comparison is near chance. The strongest prior score was an
  AUROC of 0.595. The HLA-only model is roughly 0.61 to 0.65. The benchmark
  rewards remembering how often each allele was a yes. It does not show that
  allele frequency predicts tumor immunity.

One group, no independent repeat. Cite the problem. Do not cite their
debiasing method as the fix.

## MHCflurry is clean for this use, and loosely pinned

MHCflurry is trained on binding affinity and on mass-spectrometry eluted
ligands. It is not trained on T-cell yeses or on CEDAR. The authors say so:
their work covers presentation, not T-cell recognition. Weighting its score
does not leak the immune-cell label.

Flanking sequence adds a little (mean about 2% on positive predictive value).
Combining binding with presentation adds much more (about 40% to 120%). The
flanks are worth passing in. A separate cleavage model is not where the
signal is.

The license is Apache-2.0. Cite PMID 32711842. PMID 33091335 is an erratum
about a figure axis.

Two defects in how this repository uses it, not yet fixed:

- The dependency is `mhcflurry>=2.0` with no upper bound, so a run can pick
  up newer code and newer weights than the paper. A run manifest should
  record the package version, the download release, and the weight-archive
  URL.
- The allele string we store is the string we sent in. MHCflurry will accept
  several spellings and we throw away its canonical form.

## Allele spelling will fake a clean overlap

The same allele is written four ways in files we might join:

| Source | Example |
| --- | --- |
| PRIME | `A0201` |
| IEDB / CEDAR | `HLA-A*02:01` |
| IMPROVE, and NetMHCpan 4.2 CEDAR files | `HLA-A02:01` |

An overlap check that does not normalize these will report zero shared
peptides when the peptides are the same. `normalize_allele` exists in
`src/neoantigene/hla.py`. The overlap command that would use it on two
ranking files is not built.

## NetMHCpan-4.2's cancer-peptide mode saw all of CEDAR

Nilsson, Greenbaum, Peters, and Nielsen, *Frontiers in Immunology* 2025.
DOI 10.3389/fimmu.2025.1616113.

The paper says the shipped models were fine-tuned on the entire IEDB or
CEDAR set, without holding points out for an external test. The released
training files match. `c000_cedar` through `c004_cedar` hold 5,172 records:
1,188 yes and 3,984 no, over 92 alleles. One pair is repeated, so 5,171
distinct pairs. The companion test file has 1,486 pairs, and all 1,486 are
inside those training files.

Say "cancer-peptide mode" or "neoepitope mode." Do not say the authors
cheated. Do not retrain NetMHCpan and sell it. The academic license does
not allow that.

Two limits, both from the files:

- NetMHCpan-4.1 has no CEDAR training partition. A pipeline on 4.1 is not
  hit by this. Most published neoantigen papers used 4.1.
- In 4.2's binding-affinity files, 706 of 5,027 distinct CEDAR peptides also
  appear (14.0%). Those rows are binding measurements, not the immunogenicity
  labels. It is a smaller, different exposure.

A separate open question: the 4.2 filter kept neoepitopes "with HLA allele
typing" and did not require the restriction to have been measured. No
public breakdown was found of how many CEDAR restrictions are predicted
rather than experimental. A predicted restriction would bias every
binding-model score on CEDAR, on top of the training overlap.

## Frameshift tails were being joined wrong

Read from Ensembl's VEP `DownstreamProtein` plugin. The plugin calls itself
version 2.4 on every branch, so the version number cannot tell you the
convention. The Ensembl release can.

Before release 112, the tail starts at the annotated protein position. In
releases 112 and 113, a variant that hits the third base of a codon drops
that codon from the tail. From release 114, a length-change field says where
the tail starts, but that index is not permission to copy the normal protein
into the gap.

An insertion that skips an unchanged residue should keep the reference
residue. A deletion that skips a new residue must use the amino acid the
annotation states, or refuse to emit a peptide. Copying the normal residue
into a deletion builds a protein that is neither the normal one nor the true
mutant, and the old safety check does not catch it.

That correction is in the code. Releases 112 and 113, and any VCF with no
VEP header, stay on the annotated join. When both a full frameshift sequence
and a length change are present, the full sequence is used and the length
change is ignored, because the length change describes the other field.

## What the audit did not answer

Whether each CEDAR HLA restriction was measured or predicted. No source
found in the pass breaks that fraction out.
