# Audit findings

Verified observations about published neoantigen immunogenicity work, each
traceable to a file we read or a computation we ran. Claims we could not
verify from a primary source are marked as such and should not be quoted.

This exists because the field's benchmarks are hard to compare and the
reasons are specific. It is also where our own errors get recorded: three of
the entries below are corrections to things this repository previously stated.

Method: a fifteen-agent sweep over primary sources, with every quantitative
claim required to resolve to a retrievable file. Agents were read-only. Where
an agent's finding disagreed with a paper's prose, the released code or data
was treated as authoritative.

## IMPROVE

Front Immunol 2024, DOI [10.3389/fimmu.2024.1360281](https://doi.org/10.3389/fimmu.2024.1360281).
Repository `github.com/SRHgroup/IMPROVE_paper`. All figures below were read
directly from the released files rather than from the paper.

### The two released label tables

| | rows | positives | prevalence | distinct alleles |
| --- | ---: | ---: | ---: | ---: |
| `In_house_neoepitope_for_CV.tsv` | 17,520 | 467 | 2.67% | 36 |
| `Neoepitopes_CEDAR_benchmark_data.tsv` | 2,436 | 548 | 22.50% | 81 |

Both have exactly three columns: `Mut_peptide`, `HLA_allele`, `response`.
Allele spelling is `HLA-A02:01`, with no asterisk.

**There is no patient identifier, and no feature columns.** The absence of a
patient identifier is the binding limitation: patient-grouped cross-validation
is impossible on this table, and row-wise splitting inflates AUC by about 0.14
on data with no signal at all (see `learning.train`). The table can establish
a feature's sign and rough effect size. It cannot support a weight you would
stake a synthesis budget on.

### The two tables do not overlap

Exact set intersection, computed on the released files:

```
peptide + allele : in-house 17,353   CEDAR 2,436   overlap 0  (0.0%)
peptide alone    : in-house 15,293   CEDAR 2,404   overlap 0  (0.0%)
CEDAR positives 548, of which already in the in-house set:  0  (0.0%)
```

**Contamination is not uniform.** IMPROVE separated its training data from
CEDAR completely. An audit that assumed otherwise would be committing the
error it exists to detect.

The 8.4-fold prevalence difference is its own warning: CEDAR positives are
peptides somebody chose to test and publish, the in-house set is a prospective
screen. Absolute precision and PPV are not comparable across the two.

### Feature importances reproduce exactly

From `results.zip`, `5_fold_CV/TME_excluded/Feature_importance_TME_excluded.txt`,
22 features × 250 records each (5 partitions × 50 resamples):

| rank | feature | mean importance |
| ---: | --- | ---: |
| 1 | `HydroCore` | 0.1112 |
| 2 | `Prime` | 0.0849 |
| 3 | `RankBA` | 0.0788 |
| 4 | `PropHydroAro` | 0.0777 |
| 11 | `SelfSim` | 0.0428 |
| 13 | `DAI` | 0.0416 |
| 22 | `Foreigness` | 0.0104 |

### `DAI` is already a rank ratio, and it is weak

Two independent confirmations in their own source:

- `bin/R_script/02_feature_data/02_2_1_add_features.R:42` —
  `mutate(DAI = RankEL/RankEL_wt)`, where `RankEL` is renamed from
  `Mut_MHCrank_EL` and `RankEL_wt` from `Norm_MHCrank_EL`.
- `bin/python_script/src/multimerPatientTools.py:173` and `:200` —
  `dfMerge['Agrotopicity'] = dfMerge['%Rank_EL_mut']/dfMerge['%Rank_EL_wt']`.

So IMPROVE's agretopicity is a ratio of NetMHCpan **eluted-ligand %Ranks**,
not of affinities. It ranks 13th of 22 at 467 positives.

**This settles a pending question against us.** ICERFIRE reports that
rank-based agretopicity beats affinity-based (AUC 0.589 vs 0.539) on CEDAR,
which invited reformulating our `agretopicity` as a %rank ratio. That
reformulation has already been tested, at 467 positives, on a prospective
screen rather than a literature-curated set, and it placed mid-table. The
reformulation is not supported.

### `HydroCore` differs from our `hydrophobicity` in more ways than P3

`HydroCore` is a rename of `MeanHydroph_coreNoAnc`
(`bin/python_script/feature_calculations.py:88`), computed as a plain mean of
Kyte-Doolittle values (`multimerPatientTools.py:296-307`) over a field called
`CoreNonAnchor`. The scale is the same as ours. What it is computed *over* is
not:

- It operates on **NetMHCpan-4.1's predicted 9-mer binding core**, not on the
  raw peptide. Ours slices the raw peptide and infers anchors from length.
- It then **appends bulged residues** —
  `CoreNonAnchor + PeptMut[Gp : Gp+Gl]` (`:179`, `:205`) — the residues
  NetMHCpan's alignment reports as protruding from the groove. We have no
  notion of a core, a gap position, or a bulge.
- **Their own codebase defines the non-anchor core two different ways**:
  `core[3:-1]` at `:177` (positions P4–P8 of the core) and `core[2:-1]` at
  `:203` (P3–P8). The paper's Methods describe the first.

Our `tcr_contact_positions` excludes P1, P2 and the C-terminus, giving P3–P8
for a 9-mer — which is **exactly** their `core[2:-1]` variant, modulo
operating on the raw peptide rather than a predicted core.

**This corrects a recommendation we had accepted.** The advice was to drop P3
from our feature to align with IMPROVE, on the grounds that P3 is an auxiliary
anchor where Chowell documented a sign reversal. But IMPROVE is not
self-consistent about P3, and our definition matches one of its two variants.
The substantive gap is core alignment and bulge handling, not P3. Dropping P3
alone would not align the definitions and would change our feature for a
reason the source does not support.

### The Table 1 p-values may not be adjusted

`bin/R_script/99_functions.R:157` and `:226` pass
`test.args = list(p.adjust.method = "bonferroni")` to `wilcox.test`.
`wilcox.test` has no `p.adjust.method` parameter, so the argument is absorbed
by `...` and silently ignored. Only two of roughly eight `wilcox.test` calls
in the plotting code even attempt an adjustment. Treat the published p-values
as unadjusted unless the authors clarify.

This does not change the direction of any finding here: `HydroCore` at
p ≈ 1.6e-12 and `DAI` at p ≈ 0.96 are far enough apart that a Bonferroni
factor over 22 features does not reorder them.

### Licensing

The repository contains **no LICENSE file**. There is no explicit grant to
reuse the data or code. Resolve this with the authors before building on it.

## PRIME

The widely repeated figure is that PRIME had ~70% of its training data inside
the CEDAR benchmark. Verified, with three qualifications that the short form
drops:

- It concerns **PRIME 1.0**, not PRIME 2.0.
- The measurement is **third-party**, by the IMPROVE authors, not a
  self-report.
- CEDAR **postdates** PRIME 1.0, so this is not a case of evaluating on a
  benchmark one trained against. PRIME's authors used leave-one-study-out
  cross-validation specifically because "standard cross-validation results can
  be artificially boosted by batch effects because different peptides from the
  same study are found in both the training and testing sets" — the same class
  of precaution as grouping our own folds by patient.

Training peptides are public: PRIME 1.0 in Supplementary Table S1, PRIME 2.0
in Supplementary Table S4. CEDAR is bulk-downloadable from IEDB. A mechanical
overlap computation is therefore feasible for both versions. A published
PRIME 2.0 ∩ CEDAR figure was not found.

Allele spelling in PRIME's released files is `A0201`.

## The HLA-only result

Zhang et al., *Cell Genomics* 2026, DOI
[10.1016/j.xgen.2026.101214](https://doi.org/10.1016/j.xgen.2026.101214)
(ImmuBPI, model renamed ImmUni). A model given only the HLA allele, with the
peptide masked, outscored all fifteen entrants in a published comparison.
Peer-reviewed; the claim survived review unchanged.

Three qualifications are mandatory, because the unqualified version implies
something stronger and false:

- The fifteen are predictor **configurations from a single comparison table**,
  not fifteen independent tools. Three belong to that table's own authors.
- The benchmark is IEDB's **infectious-disease** immunogenicity set — not
  neoantigens, not CEDAR. On the neoepitope benchmark the HLA-only variant was
  only best among the authors' own three input variants.
- The entire comparison sits **near chance**: the strongest prior entrant
  scored AUROC 0.595, and the HLA-only model reaches roughly 0.61–0.65. The
  finding is that the benchmark rewards memorizing each allele's positive
  rate. It is **not** that allele frequency predicts immunogenicity.

Single group, no independent replication, and the group's own debiasing method
is the proposed fix. Cite the problem, not the solution.

## MHCflurry

Our own presentation backend, audited before anyone else's.

**It is clean for immunogenicity evaluation.** Every component is trained on
peptide-MHC binding affinity and mass-spectrometry eluted ligands. No
component uses T-cell immunogenicity labels, IEDB's T-cell assay table, or
CEDAR. IEDB is used, but only its MHC-ligand table. The authors state the
scope limit themselves: their work "only addresses the steps contributing to
MHC class I ligand presentation, not T cell recognition of presented
epitopes."

That matters for us specifically: the presentation score is genuinely upstream
of the label it is evaluated against, so weighting it at 2.5 does not leak the
outcome.

Flanking context contributes a mean **+2.1% PPV** (range 0.14–3.9%), against
**+41% to +120%** for integrating binding with presentation. Feeding real
flanks is worth doing and is already done; a separate cleavage predictor is
not where the signal is.

Code and weights are Apache-2.0, commercial use permitted. Cite PMID 32711842,
not 33091335, which is a figure-axis erratum.

### Two defects in how we use it

- `pyproject.toml` pins `mhcflurry>=2.0` with **no upper bound**, so a run
  loads current 2.3.x code and the current dated weight release rather than
  the 2.0 weights the paper describes. The "clean" property is inherited from
  current documentation, not from the published methods. Three handles are
  queryable and all three belong in the run manifest:
  `mhcflurry.__version__`, `mhcflurry.downloads.get_current_release()`, and
  the dated tarball URL recorded in the weights directory's
  `DOWNLOAD_INFO.csv`.
- Our backend returns MHCflurry's `sample_name`, which is the allele string we
  passed in, so `PresentationCall.allele` is never canonicalized. MHCflurry
  accepts `A*02:01`, `HLA-A0201` and `HLA-A*02:01` and emits a canonical form
  we discard. Upstream spelling inconsistency propagates through our outputs.

## Allele spelling is not standardized, and that is a hazard

Three conventions observed in three datasets we intend to join:

| source | spelling |
| --- | --- |
| PRIME released files | `A0201` |
| IEDB / CEDAR | `HLA-A*02:01` |
| IMPROVE released tables | `HLA-A02:01` |
| NetMHCpan 4.2 `c00*_cedar` | `HLA-A02:01` |

Any overlap computation must normalize centrally. Leaving normalization to the
caller is how a leakage audit silently reports zero overlap because the keys
never matched.

## NetMHCpan-4.2 and CEDAR

Nilsson, Greenbaum, Peters and Nielsen, *Front Immunol* 2025, DOI
[10.3389/fimmu.2025.1616113](https://doi.org/10.3389/fimmu.2025.1616113).
The paper's own methods, and the training files shipped beside the method.

The paper states that the models included in NetMHCpan-4.2 "were fine-tuned
on the entire IEDB or CEDAR epitope training sets without removed data points
for the purpose of external test set construction." The CEDAR construction in
the same methods section is 5,172 data points (1,188 positive, 3,984
negative).

The released files match that count and then go further. `NetMHCpan_train`
contains `c000_cedar` through `c004_cedar`: **5,172 records, 1,188 labeled
positive and 3,984 labeled negative, over 92 alleles.** One pair is repeated
(`SLLSLPLSL` / `HLA-A02:01`), so there are 5,171 distinct pairs. The
companion `cedar_test` file has **1,486 pairs, and all 1,486 are in those
training files.** For the neoepitope mode that ships, the paper's held-out
CEDAR split is not held out.

Two limits, both checked on the released files rather than inferred:

- **NetMHCpan-4.1's training directory has no CEDAR partition.** It contains
  `c00*_ba`, `c00*_el`, `allelelist` and `MHC_pseudo.dat` only. A pipeline on
  4.1 is not contaminated by this finding. Most published neoantigen work
  used 4.1.
- **The 4.2 binding-affinity partition is a different exposure.** 706 of the
  5,027 distinct CEDAR peptides also appear in `c00*_ba` (14.0%). Those rows
  are binding measurements. That is peptide-level exposure of the default
  presentation mode, not the immunogenicity labels sitting in `c00*_cedar`.

Allele strings in the CEDAR training files are `HLA-A02:01`: locus, hyphen,
no asterisk. Joining them to an IEDB export of `HLA-A*02:01` without
normalizing reports a false zero.

## VEP `DownstreamProtein` does not always start at `protein_position`

Read from `Downstream.pm` on the Ensembl VEP_plugins branches. The plugin's
own version string is `2.4` on every one of them, so the version cannot
select the convention.

| Ensembl release | where `DownstreamProtein` starts | `ProteinLengthChange` |
| --- | --- | --- |
| 111 and earlier | `protein_position` | `min(translation start, end) + len(tail) - len(reference)` |
| 112 and 113 | one residue later, when the variant hits the third base of a codon | the same older definition |
| 114 and later, including current `main` | the same one-residue-later case | `len(mutant peptide) - len(reference)` |

The 112 change is `last_complete_codon = low_pos - (low_pos % 3)` after the
variant has already been written into the CDS. Offset `3m` is the start of
codon `m+1`, so a deletion or insertion whose first affected base is the
third base of a codon drops that codon from the tail. The header still says
the tail includes "any amino acids overlapped by the variant itself."

**This corrects a claim this repository made.** The earlier note said a
wrong join would surface as `ReferenceMismatch`, because a tail identical to
the reference is rejected. The dropped residue is not the reference residue.
`wildtype[:protein_position - 1] + tail` then builds a protein that differs
from both the reference and the true mutant, and the check does not fire.
Peptides from that protein are silently wrong.

From release 114 the length change says where the tail starts in the mutant:
`tail_start = len(reference) + ProteinLengthChange - len(tail) + 1`. That
index is not a license to copy the reference into the gap. Two cases share
`tail_start = protein_position + 1`:

- An insertion between codons skips a residue that was never changed.
  Keeping that reference residue and then the tail reconstructs the mutant.
- A deletion of the third base of a codon skips a **new** amino acid, present
  in neither the reference nor the tail. Writing the reference residue there
  is a different wrong protein. The generator uses that residue only when
  `Amino_acids` states it (`K/N`, or `K/NX`). When the annotation is `K/X`
  or empty, it raises `FrameshiftNotSupported` instead of emitting a peptide.

The same formula applied to a 112 or 113 length change would move tails that
were already aligned, so those releases are left on the annotated position. A
VCF with no `##VEP` header is left there too: the number required to choose
is the Ensembl release, and the plugin version cannot supply it.

pVACtools' `FrameshiftSequence` is the whole mutant protein from residue 1,
including the altered residue. It is preferred over `DownstreamProtein` when
a VCF carries both, sliced at `protein_start` on read, and a
`ProteinLengthChange` sitting next to it is discarded, because that number
describes `DownstreamProtein`.

## Open, pending a primary source

- **Whether a given CEDAR HLA restriction was measured or predicted.** The
  4.2 filter kept neoepitopes "with HLA allele typing." It did not require
  the restriction to be experimental. A prediction-guided restriction would
  bias every binding-predictor evaluation on CEDAR, separately from the
  training-file overlap above. No breakdown of that fraction was found.
