# Neoantigen immunogenicity: open labeled data, evaluation methodology, and leakage

Research report, compiled 30 September 2026.

Scope: what experimentally-labeled neoantigen immunogenicity data can actually be
obtained and reused; how a small-sample neoantigen ranker should be evaluated so the
result means something; and the specific, named ways immunogenicity models in this
field have been shown to be overfit, leaky, or non-replicable.

Everything below is separated into **verified** (I fetched the artifact and inspected
it, or read the primary text), **reported** (stated in a primary source I read but not
independently checked), and **unverified** (I could not confirm; flagged as such).
Byte counts, row counts, and label counts marked "verified" were computed during this
research from the actual downloaded files.

---

## 0. Executive summary

**The measurement problem is real and worse than it looks.** I computed the exact power
of your current design (`scripts/power_precision_at_k.py`). On the Ott 2017 pool
(N=139 assayed pairs, P=14 positives, k=20), an exact one-sided hypergeometric test
needs **5 of 20** top-ranked peptides to be hits before p ≤ 0.05. A ranker that is
genuinely **2x** better than chance has **37% power**; you need a **3x** ranker to
reach 79%. On Bulik-Sullivan CU04 (N=19, P=3) there is no test at all: k=20 exceeds the
pool, and even at k=10 **no attainable hit count reaches p ≤ 0.05**. Your null results
are the expected output of these designs, not evidence about your model.

**The non-obvious part:** power against chance at fixed prevalence and fixed k is
**invariant to pool size**. Going from N=139 to N=4,000 at 10% prevalence leaves power
at ~37%. Hits in the top-k are approximately Binomial(k, prevalence), so growing the
assayed pool buys nothing. Only three things buy power: larger k, more patients
(independent strata), or a paired design against the baseline.

**There is far more open, reusable, minimal-peptide labeled data than your README
assumes.** Three resources are openly downloadable today with no DUA:

| Resource | Screened rows | Positives | Minimal peptides | HLA | Variant coords |
| --- | --- | --- | --- | --- | --- |
| NCI Surgery Branch (Gartner, NIH figshare) | 1.09M pMHC / 13.3k variants | 185 variant-level; **27** with confirmed minimal epitope | yes, 8–11mers | yes (mostly inferred) | yes |
| IMPROVE in-house (GitHub TSV) | 17,520 | 467 | yes, 8–11mers | yes | no |
| IMPROVE/CEDAR benchmark (GitHub TSV) | 2,436 | 548 | yes | yes | no |

The NCI Mmps files also ship precomputed MHCflurry 1.6 presentation, MHCflurry WT:Mut
rank, T-cell-contact hydrophobicity, anchor/TCR-contact flags, expression decile and
VAF decile — i.e. **near-exact analogues of seven of your nine features**, already
joined to labels.

**Two of your features are contradicted by the two largest independent analyses.**
Agretopicity/DAI shows no significant separation between immunogenic and
non-immunogenic neo-peptides in CEDAR (n=16,602; Wilcoxon p=0.3056) and in ITSNdb
(n=199; p=0.25). Your `agretopicity` weight of 0.5 is the third-largest in the model
and is not supported by the two best-powered tests of it. Conversely `hydrophobicity`,
which you ship at 0.0, is the **top-ranked feature** in IMPROVE.

**The leakage risks that apply specifically to you.** (1) NetMHCpan-4.2 (2025) is
fine-tuned on CEDAR neoepitopes, so CEDAR is no longer a clean test set for any model
using recent NetMHCpan as a feature. (2) HLA-A\*02:01 is enriched for *negatives* in
CEDAR (OR 0.55), the opposite of the naive assumption. (3) Two independent groups
showed that per-allele positive/negative imbalance is a shortcut that models learn
instead of immunogenicity — an HLA-allele-only model beat all 15 published predictors
on one benchmark.

---

## Part A — Open, reusable, experimentally-labeled datasets

### A.1 Summary table

Access: **Open** = direct download, no registration. **Registered** = free account.
**DUA** = data use agreement with institutional signing official.

| Resource | What is labeled | Assay | Restriction | Pos / Neg | Min. peptides | Variant coords | Access | Locator |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| **NCI SB Mmps (test)** | peptide–HLA, variant-level label propagated | TIL co-culture: IFN-γ ELISPOT + 41BB/OX40 upregulation, post-expansion | 27 experimental; rest model-paired | 13,400 propagated / 1,092,359 screened-neg (2.62M rows total) | **yes** | **yes** | Open | [10.35092/yhjc.c.4792338.v2](https://doi.org/10.35092/yhjc.c.4792338.v2) |
| **NCI SB Nmers** | 25-mer variant | same | n/a (HLA-agnostic screen) | 139 train / 46 test positives; 9,543 + 3,768 screened | no (25-mers) | **yes** | Open | same |
| **IMPROVE in-house** | peptide–HLA | DNA-barcoded pMHC multimer, PBMC/TIL/infusion product | experimental (multimer built with that HLA) | 467 / 17,053 | **yes** | no | Open | [GitHub TSV](https://raw.githubusercontent.com/SRHgroup/IMPROVE_paper/main/neoepitope_tabels/In_house_neoepitope_for_CV.tsv) |
| **IMPROVE CEDAR bench** | peptide–HLA | mixed (CEDAR-derived) | 4-digit required | 548 / 1,888 | **yes** | no | Open | [GitHub TSV](https://raw.githubusercontent.com/SRHgroup/IMPROVE_paper/main/neoepitope_tabels/Neoepitopes_CEDAR_benchmark_data.tsv) |
| **CEDAR** | peptide, assay-level | ELISPOT 68%, multimer 12%, bioactivity 10%, ELISA 7%, ICS 2% | 11,460/20,601 assays have any; **6,369 (31%) 4-digit** | 2,178 / 14,424 | yes | partial | Open | [cedar.iedb.org/database_export_v3.php](https://cedar.iedb.org/database_export_v3.php) |
| **ITSNdb** | peptide–HLA | tetramer titration + IFN-γ/TNF ELISPOT, *plus* mandatory MS/transfection presentation proof | experimental | 129 / 70 | yes (9–10mers) | gene+AA only | Open | [github.com/elmerfer/ITSNdb](https://github.com/elmerfer/ITSNdb) |
| **ICERFIRE training** | peptide–HLA | CEDAR-derived | 4-digit + WT annotated | 631 / 2,402 (63 alleles) | yes | no | Open | doi [10.1093/narcan/zcae002](https://doi.org/10.1093/narcan/zcae002) |
| **NetMHCpan-4.2 CEDAR test** | peptide–HLA | CEDAR-derived, 8-mer-overlap purged | 4-digit | 315 / 1,171 | yes | no | Reported | doi [10.3389/fimmu.2025.1616113](https://doi.org/10.3389/fimmu.2025.1616113) |
| **TESLA Table S4** | peptide–HLA | **pMHC multimer** (two independent labs) | experimental | **37 / 571** (608 total) | yes | no (no WT, no gene in public ver.) | Publisher supplement | Cell doi [10.1016/j.cell.2020.09.015](https://doi.org/10.1016/j.cell.2020.09.015) |
| **TESLA Table S7** | peptide–HLA | tetramer | experimental | **4 / 306** (310 total) | yes | no | Publisher supplement | same |
| **TESLA Synapse** | raw WES/RNA FASTQ + BED **only** | — | — | **no peptide labels at all** | no | n/a | **DUA** | [syn21048999](https://www.synapse.org/Synapse:syn21048999) |
| **NEPdb** | peptide–HLA | literature-curated | mixed | **173 / 17,376** | yes | partial | Registered | [nep.whu.edu.cn](http://nep.whu.edu.cn/) |
| **TSNAdb v2.0** | peptide–HLA | literature + DB aggregation | mixed | 1,856 total (67 tier-1) / **no negatives** | yes | yes | Open | [pgx.zju.edu.cn/tsnadb](https://pgx.zju.edu.cn/tsnadb/) |
| **dbPepNeo2.0** | peptide | MS + TCR | mixed | 801 HC / no curated negatives | yes | partial | Open | [biostatistics.online/dbPepNeo2](http://www.biostatistics.online/dbPepNeo2) |
| **NDD v0.2** | peptide–HLA | literature-curated | mixed | **257 / 0** | yes | yes | Open | [github.com/NeoDiscovery/NDD](https://github.com/NeoDiscovery/NDD) |

### A.2 TESLA / Tumor Neoantigen Selection Alliance

Wells et al., *Cell* 183:818–834.e13 (2020). PMID **33038342**, PMC **PMC7652061**,
doi **10.1016/j.cell.2020.09.015**.

**What the study is.** 28 teams submitted ranked neoepitope predictions on shared
tumor sequencing from 6 subjects (3 melanoma, 3 NSCLC); 25 teams were analysed.
Submissions ranged 7–81,904 ranked pMHC per sample (median 204). **608** peptides
selected from across all teams' top-ranked lists were tested by pMHC multimer assay
(median 97/subject, range 73–144); **37 (6%)** were immunogenic. An independent
validation cohort of 3 melanoma patients contributed **310** tetramer-tested pMHC of
which **4** were immunogenic. *(Verified — read from the Cell full text.)*

**Assay detail.** Multimer binding was performed independently by two groups: the
Immunomonitoring Laboratory at Washington University in St. Louis and the Netherlands
Cancer Institute. Fluorescent-streptavidin pMHC-I multimers, two fluorochromes per
monomer, combinatorial encoding; plus a DNA-barcoded nanoparticle format ("pNP") for
multiplexed analysis. Read out on PBMC or tumor lysate with CD45/CD8 positive gating
and a CD4/CD14/CD16/CD19 dump channel. This is **multimer binding, not functional
cytokine release**, and it is **post-expansion** in the general case — a materially
different label from ELISPOT. *(Verified.)*

**What is actually public — this matters and is widely misstated.** I enumerated
Synapse `syn21048999` anonymously via the Sage REST API. It contains exactly one child
folder, `syn23446508` ("TESLA phase 1 data release"), which contains three folders:

- `syn10141118` "lung cancer" — 12 files, all `TESLA_{86,87,90,91,94,97}_{1,2}.fastq.gz`
- `syn16810855` "melanoma 1" — one `EXOME_REGIONS_*.bed.gz` plus a FASTQ folder
- `syn15672147` "melanoma 2" — one BED zip plus `patient{04,08,09}_{normal,tumor}.fastq`
  and tumor RNA FASTQs

**There are no peptide validation tables on Synapse.** The 608-pMHC and 310-pMHC label
tables exist only as Cell **Supplemental Table S4** and **Table S7** ("Contains all
relevant information and features for all 608 pMHC tested in TESLA" / "...all 310 pMHC
tested in the validation cohort"). *(Verified: folder enumeration and supplement
captions both read directly.)*

Consequences:

- The **labels are not behind the DUA** — they are in a journal supplement. The
  **sequencing is** behind the DUA (Synapse) and, for the melanoma cohort,
  additionally in dbGaP **phs000452.v3.p1** (from Liu et al. 2019).
- I could **not** download the `mmc*.xlsx` files from this environment: every request
  to `cell.com/cms/10.1016/j.cell.2020.09.015/attachment/mmcN.xlsx` returned HTTP 403
  with a 5,577-byte HTML body (bot protection). They are normally retrievable from the
  article page in a browser. PMC7652061 is an NIH author manuscript and hosts **no**
  supplementary files. *(Verified negative.)*
- Synapse access friction is real and ongoing: a thread on `syn21048999` dated
  **2 July 2026** reports `AccessDenied` on download *despite* approved green access
  status. Budget for this if you pursue the sequencing.

**Independent critique of TESLA's public table.** The ITSNdb authors excluded TESLA
from their curated database because, in the public version, the negative peptides
"do not present the wild type (WT) sequence or the gene name", and the positive
peptides' processing and presentation were never experimentally validated. They used
TESLA negatives only as a false-positive stress test and deliberately excluded them
from their scoring, because "the TSNs selection procedure is not reproducible in a
clinical setting." *(Verified — read from the Frontiers text.)* That last point is the
important one for you: **every TESLA peptide was pre-selected by a prediction pipeline**,
so the pool is range-restricted on exactly the quantity you want to evaluate.

**TESLA's own metrics are open and worth copying.** `github.com/ParkerICI/tesla`
contains `performance-metric-functions.R` (MIT-licensed, 3,661 bytes). I read it. It
defines three metrics, each computed on an `inner_join` of the ranked list and the
tested list keyed on `(hla, sequence)` — **structurally identical to your `evaluate`
restricting to the assayed intersection**, which is good independent validation of
that design choice:

- `auprc.calculation` — AUPRC via `PRROC::pr.curve` on negated ranks.
- `fr.calculation` — "fraction ranked": of all validated peptides in the tested list,
  the fraction the method placed in its **top 100**. This is recall@100.
- `ttif.calculation` — "top-twenty immunogenic fraction": validated / total among
  merged rows with `rank < 20`. This is precision@k.

Two caveats if you adopt these. First, `ttif` filters `rank < 20`, i.e. **ranks 1–19,
not the top 20** — an off-by-one in the reference implementation. Second, both `fr`
and `ttif` are computed after the join, so the denominator of `ttif` is "tested
peptides in your top 19", which varies by method and is not k. Define your own k
explicitly rather than inheriting this.

### A.3 NCI Surgery Branch / Rosenberg lab — the systematic release you were looking for

**Yes, it exists, and it is open.** This is the most actionable finding in the report.

Gartner et al., "A machine learning model for ranking candidate HLA class I
neoantigens based on known neoepitopes from multiple human tumor types," *Nature
Cancer* 2:563–574 (2021). PMID **34927080**, doi **10.1038/s43018-021-00197-6**.
Data availability statement: "Source data are available from the NIH figshare
repository at https://doi.org/10.35092/yhjc.c.4792338.v2"; sequencing in dbGaP
**phs001003.v1.p1**; models at `github.com/JaredJGartner/SB_neoantigen_Models`.

I resolved the figshare collection (id 4792338) through the figshare v2 API. Six
articles, four of which are the data:

| figshare id | File | Size |
| --- | --- | --- |
| 11400966 | `NmersTrainingSet.txt` — all long peptides screened | 13,157,762 B |
| 11400984 | `NmersTestingSet.txt` — test-set long peptides | 4,304,986 B |
| 11400969 | `MmpsTrainingSet.txt` — all mutated minimal peptides | 7,717,448,288 B |
| 11400987 | `MmpsTestingSet.txt` — test-set mutated minimal peptides | 2,364,187,325 B |

(The other two, 11400972 and 11400975, are expression-debiased subsamples.)
No registration, no DUA. *(Verified — downloaded the two Nmers files in full and
range-requested + streamed the Mmps test set.)*

**Nmers (long peptide / variant level) — verified contents.** 23 columns:
`Variant key` (e.g. `14:31085636-31085636 C>T`), `seen in RNA`, patient `ID`,
`tumor type`, 1000G frequency, dbSNP ID, COSMIC info, CGC tier, `Gene Name`,
per-transcript AA changes, `transcript ID`, `Mutation type`, `Exon`, `cDNA change`,
`AA change`, `Wt Epitope` (25-mer), `Mut Epitope` (25-mer),
`Gene Expression Decile`, `Exome VAF Decile`, `Screening Status`, comments, key.

Row and label counts I computed:

| File | Rows | `unscreened` | screened negative | positive |
| --- | --- | --- | --- | --- |
| `NmersTrainingSet.txt` | 30,887 | 21,344 | 9,404 (`-`) | **139** (`CD8`) |
| `NmersTestingSet.txt` | 8,782 | 5,014 | 3,722 (`0`) | **46** (`1`) |

139 + 46 = **185**, matching the paper's "dataset of 185 neoepitopes." Note the label
encoding is **inconsistent between the two files** (`-`/`CD8` vs `0`/`1`) — a real
trap. You must drop `unscreened` rows; treating them as negatives would add 26,358
unlabeled rows to your negative class.

**Mmps (minimal peptide / peptide–HLA level) — verified contents.** 62 columns. Beyond
the Nmers genomic fields it adds `HLA`, `Peptide Wild Type`, `Peptide Mutant`,
`core Mutant`, and precomputed features that map almost one-to-one onto yours:

| NCI column | Your feature |
| --- | --- |
| `MHCflurry1.6 Mutant Presentation score` | `presentation` |
| `MHCflurry1.6 WT:Mut Rank`, `NetSTAB Mut:WT` | `agretopicity` |
| `T-cell Contact Residues Hydrophobicity` | `hydrophobicity` |
| `is anchor`, `T-Cell contact` | `mutation_exposure` |
| `Gene expression decile` | `expression` |
| `Exome VAF Decile` | `clonality` (proxy, no purity correction) |

plus netChop C-term and 20S scores, TAP binding, NetMHCstabpan, IEDB immunogenicity
score, MixMHCpred, NetMHCpan EL and BA ranks, and HLAthena ranks — each for both
mutant and WT, with deltas. There is no `self_dissimilarity` or `tumor_selectivity`
analogue.

For `MmpsTestingSet.txt` I streamed all 2.36 GB and counted:

- **2,622,623** peptide–HLA rows, **61** distinct HLA alleles
- `Screening Status`: 13,400 `CD8`, 1,092,359 `-`, 1,516,864 `unscreened`
- `epitope status`: **27** `1.0`, 1,105,484 `0.0`, 1,517,112 blank

**This 13,400-vs-27 gap is the single most important caveat in this whole report.**
The NCI screen is done with tandem minigenes and 25-mer peptides, i.e. **at the
variant level, HLA-agnostically**. A reactive variant's `CD8` status is then
propagated to *every* 8–11mer × *every* patient allele derived from it — that is where
13,400 comes from. Only **27** rows in the test set have `epitope status = 1`, meaning
the minimal epitope *and* its restricting allele were actually mapped experimentally.

So: if you train on `Screening Status`, your positive labels have **model-inferred HLA
restriction** and a mostly-wrong minimal peptide, and your loss will be dominated by
the pairing heuristic rather than immunogenicity. If you train on `epitope status`,
your labels are trustworthy but there are 27 of them in the test set. Neither is a
free lunch, and the distinction is not flagged in the figshare metadata.

**Upstream screening papers (context, not additional open label tables).**

- Parkhurst et al., *Cancer Discov* 9:1022–1035 (2019), doi
  **10.1158/2159-8290.CD-18-1494**, PMC7138461. 75 MSS GI-cancer patients;
  **7,654 / 10,715 (71%)** variant transcripts screened; **124 (1.6%)** immunogenic;
  neoantigen-reactive TIL in 62/75 (83%); 57 CD8 / 67 CD4; all but one neoantigenic
  determinant unique across patients. Screening is TMG/long-peptide, so this is a
  variant-level table (Supplementary Tables S2/S3), not minimal epitopes.
- Malekzadeh et al., *J Clin Invest* 129:1109–1114 (2019), doi **10.1172/JCI123791**.
  TP53 hotspot TMGs, 25-mers with 12 aa WT flanks, DC electroporation or peptide
  pulsing, read out by 41BB/OX40 flow and IFN-γ ELISPOT. Useful as an assay-protocol
  reference and a shared-antigen source for your PDAC wedge.
- The 2024 NeoExpand paper (PMC11141192) states data "available upon reasonable
  request" — **not** an open release.

One more consequence worth stating plainly: the Parkhurst 1.6% and IMPROVE 2.7%
immunogenicity rates are the **realistic prevalence** for an unselected expressed-
variant pool. Your Ott benchmark's 10.1% prevalence is inflated ~4–6x because the pool
is vaccine-selected. Precision@k measured on a 10% pool does not transfer to a 2% pool.

### A.4 IEDB and CEDAR

CEDAR: Koşaloğlu-Yalçın et al., *Nucleic Acids Res* 51(D1):D845–D852 (2023), doi
**10.1093/nar/gkac1046**, PMC9825495. Blueprint paper: *Front Immunol* 12:735609
(2021). Export endpoint: `https://cedar.iedb.org/database_export_v3.php`.

**The authoritative answer to "what fraction is tumor neoantigen with experimental
restriction" is in a 2025 meta-analysis by the CEDAR team itself:** Sette, Carri,
Marrama, ..., Koşaloğlu-Yalçın, *Cancer Immunol Immunother* 74:362 (2025), doi
**10.1007/s00262-025-04209-7** (open access; I read the full text). Snapshot January
2025, restricted to human SNV-derived neo-peptides tested against T cells **from the
same patient** (healthy-donor and cell-line assays excluded):

- **16,602** neo-peptides, from 13,490 unique SNVs in 7,887 source proteins,
  **20,601** T-cell assays, **180** studies, 28 cancer types
- **2,178 (13%)** positive, **14,424 (87%)** negative under "≥1 positive assay ⇒ positive"
- Stricter "repeatedly tested" subset (≥3 assays, ≥2 positive ⇒ positive): only
  **518** peptides survive — **301 positive (58%)**, 217 negative
- Assay mix: **ELISPOT 68%**, multimer/tetramer 12%, bioactivity/degranulation 10%,
  ELISA 7%, **ICS 2%**. IFN-γ dominates; IL-2, TNF-α, IL-5, IL-17 rare.
- **Restriction: 11,460 of 20,601 assays (56%) report any HLA restriction; only
  6,369 (31%) at four-digit resolution.** 83% of assays are class I. 143 distinct
  alleles (99 class I, 44 class II). Class I median 9 peptides/allele (range 1–1,410);
  class II median 2 (range 1–27).
- Of 5,019 neo-peptide:HLA pairs with NetMHCpan predictions, 81% have EL %Rank < 2 and
  51% < 0.5. Median %Rank 0.27 for neo-epitopes vs 0.528 for negatives (Wilcoxon
  p < 2.2e-16) — binding separates, but weakly.

**Findings that directly contradict two of your feature weights:**

- **DAI / agretopicity: no significant difference** between neo-epitopes and negative
  neo-peptides (Wilcoxon **p = 0.3056**). Mutant peptides bind more strongly than WT
  in *both* classes (p < 2.2e-16 each), so the *difference* carries no signal. This is
  the best-powered test of agretopicity that exists (n=16,602), and it is null. Your
  `agretopicity` weight is 0.5, third largest.
- **HLA-A\*02:01 is enriched for NEGATIVES** (OR 0.55, p-adj 1.36e-11), as are
  B\*07:02 (OR 0.30) and B\*08:01 (OR 0.31). Enriched for *positives*:
  B\*40:01 (OR 10.63), C\*15:02 (OR 25.42), DRB1\*11:01 (OR 28.85).
- Class II presents a **higher** proportion of validated neo-epitopes than class I
  (Wilcoxon p = 0.0001) despite being far less tested.

The authors are appropriately skeptical about their own enrichment result: common
alleles are studied to saturation so additional testing yields mostly negatives, while
rare alleles are tested mainly when immunogenicity is already suspected. Treat the
allele ORs as **confounded**, not causal. But the direction still falsifies the naive
"A\*02:01 is where the epitopes are" prior.

They also flag the negative-label reliability problem crisply: unlike infectious
disease, "neo-peptides are typically tested in a single patient at a single time,"
so false negatives are structural, and they cite prior work (Carri et al., *Front
Immunol* 15:1496204, 2024) showing peptides labeled negative later elicited responses
when retested under different conditions. Exhausted/dysfunctional tumor-reactive T
cells further depress detectability.

**IEDB restriction provenance — actionable.** IEDB does not curate predictions
(inclusion criteria: "Computer derived predictions without functional experimental
data will not be included"). But restriction *within* a curated assay is whatever the
original authors asserted, and IEDB exposes an **MHC restriction evidence code** filter
in the MHC Restriction search pane. The options include:

- `assay` — restriction established by MHC binding assay and/or elution
- **`MHC binding prediction`** — restriction assigned by prediction with no in vitro
  or in vivo confirmation
- `Not determined` — pre-dates evidence codes

**Filter out `MHC binding prediction` before using IEDB/CEDAR restrictions as labels.**
If you do not, your "experimentally restricted" positives include pairs whose allele
was assigned by the very class of model you are benchmarking. IEDB separately ships
RATE (Paul et al., *J Immunol Methods* 2015, PMC4458389; `tools.iedb.org/rate`) which
infers restriction by odds-ratio association in HLA-typed cohorts and, in its
promiscuity-aware iteration, **incorporates predicted binding affinity** — so
RATE-derived restrictions are also not prediction-free.

### A.5 Curated neoantigen databases — a critical assessment

**NEPdb** — Xia et al., *Front Immunol* 12:644637 (2021), PMID **33927717**,
PMC8078594. `nep.whu.edu.cn`.

The advertised "more than 17,000 validated human immunogenic neoantigens" is
misleading. The paper's own Table 1 text states: "the immunogenic neoantigen dataset
encompassed **173 neoepitopes and 17,376 non-immunogenic peptides** of human cancers
from **41 published papers**." So the Validated Neopeptide Dataset is **99% negatives**
and the positive set (173) is smaller than ITSNdb's (129 with much stricter criteria)
and far smaller than IMPROVE's (467). Curation was by a "semi-automatic pipeline" —
NLP keyword filtering of abstracts followed by manual extraction — which is a
reasonable method but not a guarantee of per-record accuracy. The separate "Predicted
Neopeptide Dataset" (516,036 peptides scored by NetMHCpan 4.0 and HLAthena over 16,745
COSMIC mutations) is **pure prediction with no labels**; do not let it leak into a
training set. Verdict: usable as a negative source, weak as a positive source, and the
headline number should not be quoted.

**TSNAdb v2.0** — Wu et al., *Genomics Proteomics Bioinformatics* 21:259–266 (2023).
`pgx.zju.edu.cn/tsnadb`. 1,856 experimentally validated neoantigens in three tiers:
**tier 1 = 67** (immunogenic *and* presentation-confirmed), tier 2 = 1,190
(immunogenic only), tier 3 = 599 (presentation only). Critically, these were
**aggregated from dbPepNeo, NeoPeptide, NEPdb and CAPD plus literature mining** — so
TSNAdb is largely a union of the other databases and **shares their errors**. Counting
it as independent corroboration is double-counting. It ships **no negatives**. The
bulk of the site (372,273 SNV + 137,130 INDEL + 11,093 fusion "neoantigens") is
DeepHLApan/MHCflurry/NetMHCpan prediction output, not data.

**dbPepNeo2.0** — Lu et al., *Front Immunol* 13:855976 (2022), PMC9043652.
801 high-confidence neoantigens (up from ~300), 842,289 low-confidence HLA
immunopeptidome entries, 55 class II HC, 630 neoantigen-reactive TCRβ sequences. The
842,289 "low-confidence" rows are **MS-derived HLA ligands, i.e. presentation evidence,
not immunogenicity**; conflating the LC tier with negatives would be a serious error.
It also bundles its own predictor (DeepCNN-Ineo) trained on its own contents.

**NeoPeptide** — Zhou et al., *Database* 2019:baz128, PMC6901387. Literature-curated
T-cell-defined neoantigens, positives only, no systematic negatives. Largely superseded.

**Caspepdb** — I could **not** verify this resource. Searches did not surface a primary
publication or a live server I could confirm. The most likely intended resource is
**CAPED / CAPD**, the Cancer Antigenic Peptide Database at
`caped.icp.ucl.ac.be` (Boon/de Plaen lab, UCLouvain), which catalogs human tumor
antigenic peptides but is dominated by shared tumor-associated antigens (MAGE, NY-ESO-1,
gp100) rather than private neoantigens and carries no systematic negatives. Treated as
unverified; do not cite "Caspepdb" on my authority.

**NDD v0.2** — `github.com/NeoDiscovery/NDD` and `huggingface.co/NeoDiscovery`. 257
experimentally validated **positive** peptides, 25 cancer types, 46 class I alleles,
8–13mers with mutation/HLA/assay annotations plus derived features. **Positives only
by design** in the public release; patient-level data is gated behind their leaderboard.
It is recent (2025) and I found no peer-reviewed publication describing its curation —
treat provenance as unestablished.

**ITSNdb — the best-designed small benchmark, with a fatal prevalence caveat.**
Immunogenic Tumor Specific Neoantigen database, in González et al. (Fernández group),
*Front Immunol* 14:1094236 (2023), doi **10.3389/fimmu.2023.1094236**.
`github.com/elmerfer/ITSNdb`. Inclusion criteria are the strictest in the field, and
worth restating because they are the right criteria:

1. derived from non-silent somatic SNVs, with the WT sequence verified present in the
   referenced protein;
2. **experimentally validated MHC-I binding** (mass spec or competition binding assay);
3. positive **or negative** immunogenicity by tetramer titration or IFN-γ/TNF ELISPOT;
4. experimentally validated **processing and presentation** (MS, or transfection of
   mutated genes into APCs) for the immunogenic ones.

Result: **199** 9- and 10-mers — **129 immunogenic, 70 non-immunogenic** — curated from
45 publications out of >70 reviewed, with WT counterpart, restricting HLA, gene, tumor
tissue, and anchor/non-anchor annotation. Because binding and presentation are already
proven for every entry, ITSNdb isolates *immunogenicity* from *presentation*, which no
other resource does.

Two caveats, both severe:

- **HLA-A\*02:01 is 80/199 = 40.2% of entries.** This is the most concentrated allele
  imbalance of any resource here, and it is inherited from the literature.
- **Prevalence is 65% positive.** Real assayed pools run 1.6–6%. ITSNdb is therefore
  appropriate for **AUC/ranking-discrimination** questions and **completely
  inappropriate for precision@k**, which is prevalence-dependent. Do not use it to
  estimate your KPI.

### A.6 Vaccine trials with per-peptide supplementary tables

The structural problem is the same across all of them: **vaccine trials assay only the
peptides the vaccine selected**, which a prediction pipeline already chose (typically
top ~20). The negatives are "predicted-good peptides that failed," not a representative
candidate pool. This is textbook range restriction / selection on the predictor: the
pool has been conditioned on the very variable under evaluation, which compresses its
variance and biases any measured association toward null. **This is the primary reason
your Ott 2017 benchmark is uninformative, independent of sample size.** No amount of
extra vaccine-trial data fixes it; you need pools selected *without* reference to the
predictor, which is what the NCI SB, IMPROVE and TESLA designs provide (TESLA only
partially, since it pooled many pipelines' top ranks).

| Trial | Reference | Per-peptide labels | Assay | Assessment |
| --- | --- | --- | --- | --- |
| Ott 2017 NeoVax melanoma | *Nature* 547:217, doi 10.1038/nature22991 | Suppl. Table 5 | ex vivo + post-IVS IFN-γ ELISPOT, long peptides deconvoluted | Your current benchmark. Vaccine-selected only. |
| Sahin 2017 IVAC MUTANOME | *Nature* 547:222, doi 10.1038/nature23003 | per-patient response figures | IFN-γ ELISPOT / ICS | Warehouse+mutanome design; no clean candidate-pool negatives. **Access not verified.** |
| Hilf 2019 GAPVAC-101 GBM | *Nature* 565:240, doi 10.1038/s41586-018-0810-y | Suppl. Tables 1–4 "represent raw data" | ELISPOT, ICS, multimer | **Raw data for Figs 2–4 "available from the corresponding author upon reasonable request"** — not open. MS at PeptideAtlas **PASS01284**; microarray **GSE122498**. DNA seq withheld (consent). |
| Keskin 2019 NeoVax GBM | *Nature* 565:234, doi 10.1038/s41586-018-0792-9, PMC6546179 | Suppl. Table 5 (expression + class I predictions); ELISPOT mapping in Ext. Data Figs 2–3 | ex vivo IFN-γ ELISPOT, positive = **≥55 SFU/1e6 PBMC or ≥3x baseline**, pool→subpool deconvolution, up to 21 d IVS | Explicit, reusable positivity threshold. 8 evaluable patients, median 12 peptides (7–20) of 24 aa (15–30). Tiny. |
| Ott 2020 NEO-PV-01 (mel/NSCLC/bladder) | *Cell* 183:347, doi 10.1016/j.cell.2020.08.053 | per-peptide response tables | ELISPOT + ICS | Vaccine-selected. **Contents not verified.** |
| Awad 2022 NEO-PV-01 NSCLC | *Cancer Cell* 40:1010, doi 10.1016/j.ccell.2022.08.003 | interim immune analysis | ELISPOT/ICS, de novo CD4+CD8 | 38 ITT / 21 vaccinated / 16 completed; interim analysis of 8 patients found responses to **48% of vaccine peptides**. High hit rate precisely because pool is pre-selected. |
| Palmer 2022 GRANITE (Gritstone) | *Nat Med* 28:1619, doi 10.1038/s41591-022-01937-6 | per-patient ELISPOT | IFN-γ ELISPOT, SFU/1e6 | ChAd68 + samRNA, 14 phase 1 patients, CD8 responses in all. Cassette-encoded, not individually assayed minimal epitopes. |
| Rojas 2023 autogene cevumeran PDAC | *Nature* 618:144, doi 10.1038/s41586-023-06063-y | Source Data behind Ext. Data Figs 2 & 4 | **ex vivo IFN-γ ELISPOT on pools of 15-mers overlapping by 11 aa** | 16 patients, ≤20 neoantigens each, 8/16 responders. **Pool-level, not per-minimal-peptide**; normalized spot counts are given for *immunogenic* neoantigens, no full tested-and-negative table. Unusable as ranking labels. |
| Khattak 2023 / Weber 2024 mRNA-4157/V940 + pembro | *Lancet* 403:632, doi 10.1016/S0140-6736(23)02268-7 | — | — | Randomized phase 2b (KEYNOTE-942) reporting RFS. **I found no per-peptide immunogenicity supplement.** Do not assume one exists. |
| Bulik-Sullivan 2018 (EDGE) | *Nat Biotechnol* 37:55, doi 10.1038/nbt.4313 | per-patient tables incl. CU04 | mixed | Your second benchmark. Published tables are small (CU04: 3/19). |

### A.7 2024–2026 large-scale screens

- **IMPROVE** — Borch et al., *Front Immunol* 15:1360281 (2024), doi
  **10.3389/fimmu.2024.1360281**. **17,520 neopeptides, 467 (2.7%) T-cell recognized,
  70 patients** across three cohorts (melanoma TIL-ACT 3.45%, metastatic urothelial
  carcinoma under PD-L1 CPI 2.36%, a CPI basket trial 2.16%), screened with
  **DNA-barcoded pMHC multimers** on PBMC, TIL, and (melanoma) the TIL-ACT infusion
  product. 100–1,092 neopeptides per patient. 36 class I alleles covered, responses
  restricted to 27. Candidate selection: NetMHCpan RankEL < 2, tightened to < 0.5 for
  patients with many candidates, relaxed to hit a 100-peptide floor for patients with
  few; TPM > 0.1 required.

  **This is the closest published analogue to your setup and the most useful single
  reference in this report.** I downloaded both public tables and verified them:

  - `In_house_neoepitope_for_CV.tsv` — 17,520 rows, 3 columns
    (`Mut_peptide`, `HLA_allele`, `response`); **467 positive / 17,053 negative**.
    Lengths: 10,561 9-mers, 4,716 10-mers, 1,728 11-mers, 515 8-mers.
  - `Neoepitopes_CEDAR_benchmark_data.tsv` — 2,436 rows; **548 positive / 1,888 negative**.
  - **I verified zero overlap** on `(peptide, HLA)` between the two tables — their
    deduplication holds. Good hygiene, and it means you can use one to train and the
    other to test without re-deduplicating.

  Per-allele counts I computed from the in-house table show the negative-enrichment
  pattern independently of CEDAR: B\*07:02 1,778 tested / **26** positive (1.5%),
  A\*03:01 1,243 / 13 (1.0%), A\*02:01 1,560 / 44 (2.8%), versus A\*01:01 794 / **50**
  (6.3%) and B\*40:01 798 / **38** (4.8%). A\*01:01 and B\*40:01 are ~4x richer than
  B\*07:02.

  **Limitation:** the public TSVs carry **no patient identifier**, so you cannot
  reproduce their patient-grouped partitioning or do patient-level bootstrap from these
  files alone. Features and patient IDs are in `data.zip` (103,626,234 B) in the same
  repo, which I did not download.

- **NetMHCpan-4.2** — *Front Immunol* (2025), doi **10.3389/fimmu.2025.1616113**.
  Adds transfer learning and structural features, and **fine-tunes on IEDB epitopes and
  CEDAR neoepitopes**. Constructed external test sets with ≥5 positives and ≥5 negatives
  per allele after removing 8-mer overlap with all training sources: **10,621 IEDB
  epitopes (2,298 pos / 8,323 neg)** and **1,486 CEDAR neoepitopes (315 pos / 1,171
  neg)**; reduced training sets of 23,891 IEDB and 3,297 CEDAR points. Uses label
  smoothing (0.95/0.05) and focal loss (α=0.5, γ=1.0) with positive upsampling during a
  20-epoch burn-in. *(Reported — read from the accepted text, not independently
  reproduced.)*

- **Shared neoepitope–HLA discovery** — Nat Biotechnol (2023), doi
  **10.1038/s41587-023-01945-y**. 47 common cancer neoantigens × 15 prevalent HLA
  alleles; from **>24,000** possible neoepitope–HLA combinations, biochemical +
  computational assessment yielded 844 candidates, of which **86 were verified** by
  immunoprecipitation MS in **monoallelic** engineered cell lines; TCRs (via Adaptive
  MIRA) confirmed response for selected pairs. This is **presentation** evidence with
  immunogenicity on a subset — the right resource for validating a presentation
  feature, not for immunogenicity labels.

- **Myeloid leukemia driver-mutation screen** — conference abstract. 49 driver
  mutations, all possible 8–11mers, UV-mediated ligand-exchange ELISA binding across
  **8 common class I alleles**, **15,704 unique pHLA complexes** screened; **8%** bound,
  and **89% of binders were not predicted to bind** by HLA binding algorithms; then a
  pooled DNA-barcoded pHLA-tetramer library with 10x single-cell TCR readout identified
  polyclonal responses against **29** known and novel neoantigen pHLA. The "89% of
  experimental binders were missed by predictors" figure is a strong indictment of
  binding-first funnels. **Abstract only — I could not verify a peer-reviewed version
  or a data release.**

- **Miller et al.**, *Sci Transl Med* 16:eabj9905 (2024), doi
  **10.1126/scitranslmed.abj9905**. The "identify-prioritize-validate" (IPV) platform:
  an **HLA-agnostic** bioinformatic prioritizer plus short-term in vitro functional
  assay on autologous PBMCs, interrogating 50–75 expressed mutations from a single
  50 ml blood draw, at a stated **>40% positive predictive value**. 13 patients, 8
  hard-to-treat tumor types; all 13 had detectable neoantigen-specific T cells;
  reported subsequently extended to >130 patients / 25 cancer types. Validated
  neoantigens included both driver and passenger mutations and "would not have been
  otherwise detected using an in silico prediction approach." Data files S1–S4 are
  journal supplements. Two notes: the >40% PPV is at the **mutation** level with an
  HLA-agnostic prioritizer, so it is not comparable to your per-pMHC precision@k; and
  the >130-patient extension is a press/news claim, not a published table.

- **PredIG** — *Genome Med* 17 (2025), doi **10.1186/s13073-025-01569-8**. 17,448
  pHLA pairs with validated immunogenicity, aggregated from public databases and
  literature (Sept 2023). Deliberately **mixes pathogen epitopes, tumor-associated
  antigens, non-canonical cancer antigens and neoantigens**. Their label rule is worth
  noting and is defensible: for pHLAs with multiple experiments, any positive ⇒
  positive, and remaining negative instances for that pHLA are **discarded** rather
  than kept, because "a negative finding does not invalidate" a positive one. If you
  build a training set from mixed sources, adopt this rule — but see §C.3 on why the
  pathogen fraction is a transfer hazard.

- **NeoPrecis** (2025, PMC12932759) — trains on CEDAR, tests externally on the
  Parkhurst NCI dataset. Their fairness protocol is instructive: they identified pMHCs
  **absent from the training sets of every predictor compared** (n=1,147), then excluded
  pMHCs with alleles unsupported by any predictor, leaving **438 pMHCs (228 positive)**.
  They also **rebalanced by bootstrap to a 1:3 positive:negative ratio** to match the
  balance used in the PRIME and ICERFIRE benchmarks, noting this "better reflect[s]
  realistic immunogenicity rates" — an explicit acknowledgment that published benchmarks
  run at prevalences unlike reality. External validation on the NCI data: MHC-I
  n=1,089 with **36** positives; MHC-II n=1,189 with **33**.

---

## Part B — Evaluation methodology

### B.1 What the right test is when positives are in the single digits

**Use exact tests; never a normal approximation.** With P=14 in N=139 and k=20, the
expected hit count under random ranking is 20·14/139 = 2.01. Asymptotics are worthless
here.

**Test 1 — against chance: exact hypergeometric (= one-sided Fisher).** Hits in top-k
under a random permutation of the assayed pool are exactly Hypergeometric(N, P, k).
This is the test TESLA itself used for threshold-set evaluation (they report Fisher's
exact p = 9e-5, OR 116.5 on the validation cohort and p = 8e-6, OR 348 for the combined
presentation+recognition filter). It is also the test your README already applies, and
your p-values (0.265 for the pipeline, 0.064 for arbitrary) are correct.

The critical-value structure is what you should report alongside them. From
`scripts/power_precision_at_k.py`:

```
Ott 2017 (N=139, P=14, k=20): E[hits] under random = 2.01
  hits needed for p<=0.05 one-sided: 5 of 20
  power at 1.5x enrichment: 15.8%
  power at 2.0x enrichment: 37.3%
  power at 3.0x enrichment: 78.7%
  power at 4.0x enrichment: 96.5%

Bulik-Sullivan CU04 (N=19, P=3):
  k=20 exceeds the pool (N=19); top-k is the whole pool -> no test exists
  at k=5:  hits needed for p<=0.05 = 3 of 5
  at k=10: hits needed for p<=0.05 = unreachable
```

**Report the critical value, not just the p-value.** "This pool requires 5/20 hits for
significance; we observed 4" is an honest and immediately interpretable statement, and
it makes the design limit visible in a way a bare p = 0.265 does not.

**Test 2 — against a baseline on the same pool: paired.** This is the test you actually
need and the one you are not running. Both rankers order the *same* assayed pool, so the
two top-k sets overlap heavily and the arms are not independent. Treating them as
independent throws away most of your power. The correct statistic is the difference in
hits restricted to the **discordant** slots: peptides in your top-k but not the
baseline's, versus the reverse. Peptides in both cancel. With one patient this is
McNemar/exact-binomial on discordant pairs; across patients it is a paired permutation
test on per-patient precision@k deltas (§B.2).

**Test 3 — use the whole ranking, not just the top-k.** Precision@k throws away all
information below rank k and is a discontinuous function of the scores — a peptide
moving from rank 21 to 20 changes it, a move from 139 to 21 does not. The
Mann-Whitney U statistic (equivalently AUC) uses every positive–negative pair and is
therefore far better powered at fixed P. TESLA used Mann-Whitney U for its
feature-level comparisons (binding affinity p = 0.012, tumor abundance p = 0.033,
binding stability p = 0.067 in the 310-pMHC cohort with only **4** positives —
significance at P=4, which precision@20 could never have achieved on that pool).

The field's standard compromise is **partial AUC restricted to the low-false-positive
region**. IMPROVE reports **AUC01** (partial AUC at 10% FPR) alongside global AUC
throughout, precisely "to focus on the high specificity part of the ROC curve." That is
the right metric for a shortlist product: it is continuous and well-powered like AUC,
but it only credits performance in the region a lab would actually synthesize.

**Recommended primary/secondary split:** make **AUC01 (or AUPRC)** your primary
statistical endpoint and keep **precision@k** as the reported business KPI. They answer
different questions and precision@k cannot carry an inferential claim at your sample
sizes.

**Bayesian alternative.** With P in the single digits, a Beta-Binomial posterior on
precision@k, or a hierarchical Beta-Binomial pooling patients with partial shrinkage, is
more honest than a p-value: it yields a credible interval that visibly spans most of
[0,1] and a directly usable P(new > baseline). I found **no** published Bayesian
analysis of neoantigen ranker precision@k, so this is a defensible choice rather than a
convention you can cite.

### B.2 Pooling across patients

**Do not pool peptides across patients into one flat list.** Patients differ in pool
size (IMPROVE: 100–1,092 peptides/patient), in positive count, in prevalence (cohort
means 2.16%–3.45%), and in immune competence. A flat pool lets large-pool or
high-prevalence patients dominate, and it makes Simpson's paradox available: a method
can win within every patient and lose overall, or vice versa.

**Do compute the metric per patient, then aggregate over patients.** This is the
information-retrieval convention (per-query metric, then mean over queries) and it is
what this field does:

- IMPROVE reports **per-patient partial AUC** and compares models with a **paired
  Wilcoxon test** across patients (p = 0.95 for IMPROVE TME vs IMPROVE on the
  per-patient metric).
- They also report the fraction of patients with ≥1 true immunogenic neoepitope in the
  **top 20 and top 50** per patient, comparing against RankEL and random sampling — a
  patient-level hit-rate, which is the statistic your KPI should become.

**The IMPROVE TME result is the single most important methodological warning in this
literature for you.** Adding tumor-microenvironment features raised global AUC from
0.630 to 0.652 (roc.test p = 0.01) — and produced **no improvement per patient**
(paired Wilcoxon p = 0.95). The authors diagnosed it: the delta between the two models'
scores correlated with patient cytolytic activity at Spearman 0.76–0.80, so the TME
features "were unable to distinguish the immunogenicity of peptides within patients but
favored the patients with an immunocompetent TME." The global AUC gain came entirely
from **re-ordering patients**, not peptides.

Every patient-constant feature in your model has this failure mode. `tumor_selectivity`
is patient-constant when GTEx medians are used, and `expression` and `clonality` are
patient-constant in their *scale*. A pooled AUC improvement from any of them may be
pure between-patient sorting and worth nothing to a lab ranking one patient's peptides.
**Always report the per-patient metric next to the pooled one; if they disagree, the
per-patient one is the one that matters for your product.**

**Which aggregation test.** The best empirical evidence comes from information
retrieval, where this exact problem (per-topic metrics, few topics) has been studied at
scale. Urbano, Lima & Hanjalic, "Statistical Significance Testing in Information
Retrieval: An Empirical Analysis of Type I, Type II and Type III Errors," SIGIR 2019
(arXiv **1905.11096**), computed **over 500 million p-values** across tests, systems,
measures, topic-set sizes and effect sizes with known null hypotheses. Conclusions,
which transfer directly with "topic" → "patient":

- The **paired t-test** and the **permutation test** behave almost identically; the
  t-test is "remarkably close to ideal behavior even with small sample sizes" and is
  their top recommendation for hypotheses about mean effectiveness.
- The **permutation test** remains valuable because it accommodates test statistics
  other than the mean — which is your case, since you want precision@k deltas.
- The **Wilcoxon signed-rank and sign tests** make more errors than expected and are
  "overconfident... with a clear bias towards small p-values." They recommend
  discontinuing both. (Note this conflicts with IMPROVE's use of paired Wilcoxon.)
- The **bootstrap-shift test** shows "a systematic bias towards small p-values and is
  therefore more prone to Type I errors"; they recommend discontinuation.

**So: paired t-test or paired permutation test on per-patient deltas. Not Wilcoxon, not
bootstrap-shift.** Bootstrap is still the right tool for *confidence intervals* — just
resample **patients**, not peptides, so the interval reflects the variance you will
actually face on the next patient. Note the standard caution that percentile bootstrap
CIs under-cover at small n; use BCa or bootstrap-t.

**Mixed-effects models** (patient random intercept, method fixed effect, on a
per-peptide binary outcome) are the natural generalization and let you keep
patient-level covariates. At 2 patients and 17 positives they are not estimable. Revisit
at ≳15–20 patients.

### B.3 Beating your own dominant input feature

This is the hardest claim in your product and it needs the most careful design. Your
model is `sigmoid(bias + Σ wᵢxᵢ)` with `presentation` (MHCflurry) at weight 2.5, and
your comparator is sorting by `presentation` alone. The two rankings are **nested**:
your score is a monotone function of `presentation` plus perturbations. Nested
comparisons need nested tests.

**1. DeLong's test for correlated ROC curves.** Two AUCs computed on the *same* samples
are correlated; an unpaired comparison is invalid. DeLong, DeLong & Clarke-Pearson,
"Comparing the areas under two or more correlated receiver operating characteristic
curves: a nonparametric approach," *Biometrics* 44:837–845 (1988), gives the covariance
structure. This is what the field uses: IMPROVE compared every pair of models with
`roc.test` from **pROC** v1.18.0 (Robin et al., *BMC Bioinformatics* 12:77, 2011),
reporting IMPROVE vs RankEL **p = 4.3e-6**, IMPROVE vs NNAlign **p = 0.039**, and
IMPROVE+TME vs IMPROVE **p = 0.01**. ICERFIRE and NetMHCpan-4.2 use the same machinery.
**Use `pROC::roc.test` (or `scipy`/`sklearn` + a DeLong implementation) with the paired
option.** Caveat: DeLong is asymptotic and will be unreliable at 14–17 positives; fall
back to a paired permutation test on the AUC difference, permuting the method label
within each peptide.

**2. Likelihood-ratio test for the nested model.** Since `binding_only` is literally
your model with all non-presentation weights set to zero, that is a nested hypothesis:
H₀: w₂ = … = w₉ = 0. Fit both by maximum likelihood on the same labeled pairs and
compare with an LRT (χ² with 8 df), or use AIC/BIC if you prefer information criteria.
This is the cleanest possible statement of "the extra features carry information beyond
presentation" and it is strictly stronger than comparing two rank orderings, because it
tests the coefficients directly. With 17 positives across both your cohorts, an 8-df
LRT is badly underpowered — so also fit the 1-df version (presentation + one candidate
feature) per feature, which is what your evidence actually supports.

**3. Ablation conventions.** The field's convention is leave-one-feature-out with
re-training, reported as a table. Two examples worth copying:

- IMPROVE dropped features correlated above |0.7| (removing `HydroAll` for `HydroCore`,
  `VarAlFrac` for `PriorScore`), ran backward and forward selection, and **removed
  features that added no predictive power** (their one-hot mutation-consequence and
  mutation-position features were dropped for this reason — note that
  mutation position is your `mutation_exposure`).
- More importantly, IMPROVE **retrained without the PRIME feature** specifically
  because PRIME was implicated in leakage, and showed performance held up. That is the
  gold-standard move: when a feature is suspect, retrain without it and report both.

**Do the analogous thing: fit and report your model with `presentation` removed
entirely.** If performance collapses to chance, you have learned that your model is
MHCflurry with decoration, which is a legitimate and publishable finding about the
field, and it is the honest framing of your KPI. If it degrades but stays above chance,
you have evidence of orthogonal signal. Either way, the ablation is more informative
than any head-to-head against `binding_only`.

**4. Residual / orthogonalization framing.** The sharpest version: regress your score
on `presentation`, take the residual, and test whether the **residual alone** ranks
better than chance. This directly asks "is there signal in my model that MHCflurry does
not already have," with no nesting artifacts. Equivalently, stratify the assayed pool
into presentation deciles and test for enrichment *within* strata (Cochran-Mantel-
Haenszel). I found no published example of this in the neoantigen literature, so it
would be a novel analysis rather than a convention — but it is the correct one, and it
is cheap to run.

**5. Choose the right comparator and say so.** Note that IMPROVE's headline comparison
is against **RankEL** (NetMHCpan eluted-ligand %rank), which achieved AUC 0.539 —
barely above chance on their data. Your `binding_only` baseline is the same class of
comparator. Beating a 0.539 baseline is a much weaker claim than beating a strong one,
and reviewers and partners will make that point. IMPROVE itself reached only AUC 0.630
/ AUC01 0.0139. **Realistic target: AUC ~0.63–0.65. Anything above ~0.75 on tumor
neoantigens should be treated as evidence of leakage until proven otherwise** (see
Part C).

### B.4 How many positives to detect a doubling of precision@20

I found **no published power analysis for this setting.** The closest statements are
IMPROVE's "the model needs validation with more data, especially more immunogenic
neoepitopes" and ITSNdb's "present weaknesses are: the current low amount of peptides"
— acknowledgments, not calculations. Sakai's topic-set-size design work in IR is the
right formal analogue but has not been applied here.

So I computed it. `scripts/power_precision_at_k.py`, α = 0.05 one-sided, power 0.80.

**Result 1 — pool size does not help.** Holding prevalence at 10% and k = 20, a ranker
2x better than chance has essentially constant power regardless of N:

```
  N=  139  P=  14  crit=  5/20  power=37.7%
  N=  300  P=  30  crit=  5/20  power=37.3%
  N=  600  P=  60  crit=  5/20  power=37.1%
  N= 1000  P= 100  crit=  5/20  power=36.9%
  N= 2000  P= 200  crit=  5/20  power=37.1%
  N= 4000  P= 400  crit=  5/20  power=37.6%
```

Because hits in the top-k are ~Binomial(k, prevalence) once N ≫ k, power depends only
on **k** and the **prevalence ratio** — not on how many peptides you assayed. Assaying
a bigger pool from one patient does not fix this. This is counterintuitive and it
invalidates the most natural plan ("assay more peptides").

**Result 2 — unpaired two-proportion requirement.** For a doubling of precision@20,
slots needed per arm (each patient contributes k = 20 slots):

| Baseline p@20 | Target | Top-20 slots/arm | Patients/arm | Expected positives in winning arm |
| --- | --- | --- | --- | --- |
| 5% | 10% | 343 | 18 | ~34 |
| 10% | 20% | 157 | 8 | ~31 |
| 20% | 40% | 64 | 4 | ~26 |

**Result 3 — paired permutation on per-patient deltas, which is the design you should
use.** Power depends critically on how much the two top-20 lists overlap, because only
discordant slots carry signal:

| Patients | Power (50% overlap) | Power (80% overlap) |
| --- | --- | --- |
| 5 | 9.2% | — |
| 10 | 49.8% | 18.1% |
| 20 | 84.2% | 44.3% |
| 40 | 98.7% | 75.8% |
| 60 | 100% | 89.1% |
| 80 | 100% | 96.4% |
| 120 | — | 99.5% |

**The headline numbers.** To detect a doubling of precision@20 from 10% to 20% at 80%
power you need roughly:

- **~20 patients** if your ranking differs substantially from the baseline (50% top-20
  overlap), which implies **~40–60 labeled positives** across the cohort;
- **~40–45 patients** if your ranking mostly agrees with its dominant input feature
  (80% overlap), implying **~80–120 labeled positives**.

The second row is your realistic case: a re-ranker whose largest weight is MHCflurry
presentation will produce top-20 lists that overlap the MHCflurry-sorted top-20
heavily. **Overlap is the hidden driver of your power budget**, and you can measure it
today from data you already have — compute the Jaccard overlap of your top-20 and
`binding_only`'s top-20 on the Ott pool and read the required cohort size off the table.

At 1.6–2.7% realistic prevalence rather than the vaccine-inflated 10%, the requirement
grows further; the 5%→10% row (18 patients/arm unpaired, ~34 positives) is the better
guide for an unselected pool, and paired requirements scale similarly.

**Practical implication.** ~40 patients with ~100 assayed pMHC each is ~4,000 assays —
comparable to one IMPROVE cohort, and plainly beyond a bootstrapped budget. This is the
argument for using the open NCI SB and IMPROVE data for model development and reserving
your own assay spend for prospective confirmation on a shared-antigen wedge where the
same synthesized peptide amortizes across donors. Your PDAC/KRAS wedge is the right
instinct for exactly this reason.

---

## Part C — Data leakage and replication failures

### C.1 Train/test peptide overlap and homology leakage

**The field's standard solutions, with citations you can point at.**

- **Common-motif clustering** (Nielsen lab). NetMHCpan-4.1 (Reynisson et al., *Nucleic
  Acids Res* 48:W449–W454, 2020, doi **10.1093/nar/gkaa379**) and NetMHCpan-4.2 split
  training data into five cross-validation partitions "using the common motif
  approach... ensuring that peptides sharing at least an **8-mer overlap** were placed
  in the same partition." NetMHCIIpan-4.0 / NNAlign_MA does the same at **9-mer**
  subsequences, clustering BA and EL data simultaneously before separating them.
- **Hobohm-1 plus peptide-kernel similarity.** ICERFIRE (Wan, Koşaloğlu-Yalçın, Peters
  & Nielsen, *NAR Cancer* 6:zcae002, 2024, doi **10.1093/narcan/zcae002**) identified
  redundant points with **Hobohm-1**, computed pairwise similarity with the peptide
  kernel over window sizes 3–8 at **threshold 0.9**, split only the *dissimilar*
  peptides into 10 folds, then added the held-out similar peptides back to the fold
  containing their closest match. Their training set: 3,033 unique neo-epitope–HLA
  pairs (2,926 unique sequences) over 63 alleles, **631 positive / 2,402 negative**.
- **Patient-grouped partitioning** (the one most often skipped). IMPROVE clusters on
  *both* shared motifs between immunogenic neoepitopes *and* patient identity, "so that
  all neoepitopes from the same patient and all similar neoepitopes were grouped
  together in the same partition," then distributes negatives into the established
  partitions, and deselects any test peptide that also appears in training.

**MHCflurry 2.0** (O'Donnell, Rubinsteyn & Laserson, *Cell Syst* 11:42–48, 2020, doi
**10.1016/j.cels.2020.06.010**) uses a different and relevant trick: it deliberately
trains the binding-affinity predictor on *in vitro* affinity measurements "which are
largely independent of AP... one of several design choices intended to limit the BA
predictor's tendency to learn AP signals," then trains the antigen-processing predictor
on the **residual**. That is feature-level orthogonalization to prevent one sub-model
from absorbing another's signal — the same idea as §B.3.5, and a reason your
`presentation` feature is not simply "binding."

**What this means for you.** Your `self_dissimilarity` feature (BLOSUM62 distance to
the nearest peptide in the whole human proteome) and your hard gate dropping peptides
that exactly match the reference proteome are both forms of similarity control, but
neither is train/test partitioning. When you refit on accumulated labels you must
partition by **(patient, motif cluster)**, not randomly. With `refit` documented to run
at ~40 labeled peptides, random splitting will produce a wildly optimistic estimate,
because the same patient's peptides share HLA, expression context, and often source
protein.

### C.2 IEDB-derived label contamination

Two distinct problems.

**Predictor-inferred HLA restriction.** Covered in §A.4. IEDB never curates predictions
as data, but the restriction *asserted within* a curated assay can itself have been
derived from a binding prediction, and IEDB explicitly flags this with a
`MHC binding prediction` evidence code. **Filter it.** RATE-inferred restrictions also
incorporate predicted binding in the promiscuity-aware mode.

**Benchmark contamination.** The CEDAR meta-analysis documents the structural cause:
only **2%** of the 13,490 mutations have neo-peptides tested in more than one study, and
of those 255 mutations, **213 (83%) share at least one author** with another study
testing the same mutation. Across 63 distinct study groups testing the same mutations,
the median number of shared authors is 5, range 1–21. So the apparent independent
replication in the literature is mostly **the same groups following up their own
findings**. Any benchmark assembled by pooling published studies inherits this: your
"independent test set" is likely to contain the same peptides, from the same labs, as
your training set. Only **42** mutations have been examined by genuinely
non-overlapping author groups.

The same analysis quantifies the recurrent-mutation testing bias: TP53 mutation
frequency in cBioPortal correlates with the number of studies testing its neo-peptides
at **R = 0.74 (p = 1.6e-10)**. TP53 alone contributes 152 neo-peptides from 68 amino-acid
mutations, of which **53 (35%)** are positive versus 6% overall
(χ², p = 1.51e-14). Of 72 RAS neo-peptides, **37 (51%)** are positive. Hotspot-enriched
databases will therefore report inflated positivity and reward models that memorize
hotspots.

**Most consequential for you specifically: NetMHCpan-4.2 is fine-tuned on CEDAR
neoepitopes.** So is ICERFIRE trained on CEDAR, and IMPROVE benchmarked on CEDAR. If
you use a recent NetMHCpan (or any tool refined on CEDAR) as your `presentation`
feature, **CEDAR is no longer a clean external test set for your model.** MHCflurry 2.0
predates this and is trained on affinity + MS eluted-ligand data rather than T-cell
immunogenicity labels, so MHCflurry + CEDAR is currently a *cleaner* pairing than
NetMHCpan-4.2 + CEDAR. That is an argument for keeping MHCflurry as your backend,
and for pinning the version in your run manifest (which you already do).

### C.3 Positive sets dominated by viral/pathogen epitopes

**Demonstrated, not merely asserted.** Buckley et al. (Oxford Experimental Medicine),
"Evaluating performance of existing computational models in predicting CD8+ T cell
pathogenic epitopes and cancer neoantigens," *Briefings in Bioinformatics* 23:bbac141
(2022), doi **10.1093/bib/bbac141**, PMID **35471658**. Findings:

- For SARS-CoV-2 epitopes, "**none of the models perform substantially better than
  random or offer considerable improvement beyond HLA ligand prediction**."
- "Suboptimal performance for predicting cancer neoantigens."
- They "compared key parameters associated with immunogenicity between pathogenic
  peptides and cancer neoantigens and observed **evidence for differences in the
  thresholds of binding affinity and stability**, which suggested the need to modulate
  different features in identifying immunogenic pathogen versus cancer peptides."

The mechanism is not subtle. Pathogen epitopes are non-self and face a naive,
unpurged TCR repertoire. Neoantigens differ from self by one residue and face a
repertoire shaped by thymic negative selection against the WT sequence — the ITSNdb
authors make exactly this point, that cross-reactive WT/mutant TCRs are deleted "due to
their high sequence similarity, which is what makes neoantigens different from pathogen
epitopes." A model calibrated on pathogens learns the wrong threshold and the wrong
features.

**Consequences for specific resources:** the IEDB general immunogenicity predictor,
DeepImmuno and CIImm were built for *general* immunogenicity, not neoantigens. ITSNdb
found AUCs of **0.52–0.60** across 16–19 metrics from 7 software packages, concluding
that "none of the methods could be considered the best immunogenicity predictor, since
they show very low AUCs and loss of Specificity in favor of Sensitivity." PredIG (2025)
mixes pathogens, TAAs, non-canonical antigens and neoantigens in one 17,448-pair
training set — powerful, but the pathogen fraction is a transfer hazard you must
stratify out when evaluating on tumors.

ITSNdb's explicit recommendation, which I endorse: **do not** build negatives by
sampling random human proteome sequences, and **do not** use pathogen epitopes as
positives, when training a neoantigen immunogenicity model. Negatives must be peptides
that **demonstrably bind MHC** and still failed to elicit a response; otherwise the
model learns binding, which you already have as a feature.

### C.4 Allele imbalance and intra-allele label imbalance

**Inter-allele imbalance** (a few alleles dominate) is widely known. ITSNdb:
**HLA-A\*02:01 = 80/199 = 40.2%**. CEDAR class I: median 9 peptides/allele but range up
to **1,410**.

**Intra-allele imbalance** — the positive:negative ratio *within* each allele — is the
subtler and more damaging problem, and two groups found it independently.

**Finding 1 (Buckley et al. 2022, above):** "cross-HLA variation in the distribution of
immunogenic and non-immunogenic peptides in the training data of the models seems to
**substantially confound the predictions**."

**Finding 2 (ImmuBPI; Wang lab, Tsinghua; bioRxiv 2024.02.07.579420v2; code at
`github.com/WangLabTHU/ImmuBPI`)** — a mechanistic demonstration. Collecting the
training sets of four published immunogenicity predictors, they found that **over 50%
of HLA alleles show a tenfold or greater positive:negative imbalance in all four
datasets.** Then:

- In their transformer's attention map for immunogenicity, **HLA received abnormally
  high attention** compared with the same architecture trained for binding or
  presentation — biologically backwards.
- Mean prediction score per allele in the test set correlated strongly (Spearman) with
  that allele's train-set log₂ positive:negative ratio — **and this held for
  BigMHC_IM, DeepHLApan and DeepImmuno with their own training sets**, i.e. the effect
  is **model-agnostic**, not an artifact of their architecture.
- **The decisive ablation: a model given only the HLA allele (no peptide) performed
  best on both benchmarks.** The HLA-only model "greatly surpassed (AUROC: +3.18% ∼
  +15.95%, AUPRC: +5.13% ∼ +14.01%) all 15 models evaluated" in a prior benchmark.

**That last result means published immunogenicity benchmarks can be won without looking
at the peptide.** Any reported AUC on such a benchmark is partly or wholly a measure of
how well the model's allele prior matches the benchmark's allele composition.

They also note a second-order leak: mutual information between residue identity and
label peaks at **peptide positions 2 and 9** — the MHC **anchor** positions, which have
only a secondary role in TCR contact. Since anchors determine the binding motif, "the
anchor position would to some degree leak the HLA information." So even a
peptide-sequence-only model can back-door the allele prior.

**What you must do.** (1) Report per-allele performance, never a single pooled AUC.
(2) Run the HLA-only ablation as a sanity floor — if allele alone reproduces your
performance, you have measured nothing. (3) Stratify or match allele composition
between train and test. (4) Do not assume A\*02:01-heavy data generalizes; in CEDAR it
is *negatively* enriched (OR 0.55). For your PDAC wedge, note that your target alleles
C\*08:02 and A\*11:01 are exactly the under-tested ones — good for novelty, bad for
available labels, and you will have almost no prior data to refit on.

### C.5 Named replication failures and critiques

| Predictor | Finding | Source |
| --- | --- | --- |
| **PRIME / PRIME2** | **70% of PRIME's training data was inside the CEDAR benchmark** used to evaluate it, "likely resulting in a performance overestimation"; confirmed causally — "the PRIME performance dropped when removing the peptides overlapping with the PRIME training data from the CEDAR evaluation data, showing that Prime lacks performance with the independent dataset." | Borch et al. 2024 (IMPROVE), doi 10.3389/fimmu.2024.1360281 |
| **PRIME2, BigMHC_IM** | "At least 6 out of 8 test data points actually appeared in the training set of PRIME2 and BigMHC IM, **yet the two models failed to make distinctive predictions**" on immunogenicity-changing mutation variants. Memorized without generalizing. | ImmuBPI, bioRxiv 2024.02.07.579420v2 |
| **BigMHC_IM, DeepHLApan, DeepImmuno** | Per-allele prediction scores track train-set intra-allele imbalance; HLA-only model beats all 15 benchmarked models. Model-agnostic shortcut learning. | ImmuBPI |
| **INeo-Epp** | Excluded from a 2025 benchmark outright because "its training data includes a significant portion of the benchmarking dataset, and the results would therefore be heavily biased." | neoIM, *Vaccines* 13:865 (2025), doi 10.3390/vaccines13080865 |
| **DeepHLApan** | Top-ranked on ITSNdb but "**show[s] discrepancy on validation databases**"; gave the **highest false-positive rate, 100% (297)** on TESLA negatives. | ITSNdb, doi 10.3389/fimmu.2023.1094236 |
| **DeepImmuno** | Best on the validation set but **worst on ITSNdb**; 78.9% (234) FPR on TESLA negatives. Rank order inverts across benchmarks. | ITSNdb |
| **All 7 packages / 19 metrics** | AUC **0.52–0.60** on ITSNdb. "None of the methods could be considered the best immunogenicity predictor." Predictions were **not shared between methods** — they disagree on which peptides are immunogenic. | ITSNdb |
| **All models, SARS-CoV-2 + neoantigens** | "None of the models perform substantially better than random or offer considerable improvement beyond HLA ligand prediction." | Buckley et al., doi 10.1093/bib/bbac141 |
| **DAI / agretopicity** | No significant separation: **p = 0.3056** (CEDAR, n=16,602) and **p = 0.25** (ITSNdb, n=199). TESLA separately found mutation position useless for filtering. | Sette et al. 2025; ITSNdb; Wells et al. 2020 |
| **Tumor neoantigen burden (TNB) as ICB biomarker** | "The original observed association between TNB and outcome **loses significance** upon the immunogenicity predictor used," while plain TMB stayed significant in all cohorts but one. "**TMB is still a superior biomarker than TNB.**" | ITSNdb |
| **Hydrophobicity, sign of effect** | Chowell 2015 found hydrophobic TCR-contact residues enriched in immunogenic epitopes; **TESLA found immunogenic tumor pMHC significantly *less* hydrophobic (p = 0.04)**; IMPROVE found `HydroCore` and `PropHydroAro` the **most important** features (p = 1.6e-12, p < 2.22e-16). Unresolved and sign-contested. | as cited |

**The rank-order instability is the tell.** DeepHLApan is best on ITSNdb and
inconsistent elsewhere; DeepImmuno is best on one validation set and worst on ITSNdb;
PRIME is best by AUC01 and collapses once its training overlap is removed. When
benchmark choice determines the winner, the benchmarks are measuring dataset
composition, not model quality. Expect your own numbers to move substantially between
the NCI SB, IMPROVE and CEDAR sets, and **report all three** rather than picking one.

### C.6 Leakage risks specific to your pipeline

1. **`presentation` provenance.** MHCflurry 2.0 is trained on affinity + MS eluted
   ligands, not T-cell immunogenicity, so it is comparatively clean. Do **not** silently
   swap in NetMHCpan-4.2 or PRIME as the backend: both have seen CEDAR immunogenicity
   labels, and doing so would contaminate any CEDAR-based evaluation. Your pluggable
   backend registry makes this easy to do by accident; consider recording the backend's
   training-data provenance in the run manifest alongside the version.

2. **Your `docs/CITATIONS.md` priors are literature-derived, which is a mild prior
   leak.** You freeze weights beforehand and do not tune on the benchmark — good, and
   stated honestly. But the *literature* those priors come from overlaps the studies
   that generated Ott 2017 and the CEDAR/TESLA label pools. This is weak leakage and
   unavoidable; it is worth one sentence in the benchmark doc.

3. **`refit` at ~40 labeled peptides will overfit without grouped partitioning.**
   Partition by (patient, motif cluster) per §C.1. At 40 labeled peptides from one or
   two patients, expect a fitted 9-weight model to be unidentifiable; consider
   regularization toward your priors (ridge with the prior as the penalty center)
   rather than unpenalized MLE.

4. **Two features are actively contraindicated by current evidence.**
   `agretopicity` (weight 0.5) is null in the two largest tests of it.
   `hydrophobicity` (weight 0.0) is the top feature in IMPROVE. Your README's reasoning
   for zeroing hydrophobicity — that TESLA found the opposite sign to Chowell — is
   sound and well-documented, but IMPROVE is now the larger, more recent, multimer-based
   study and it points the other way. Consider demoting `agretopicity` to near zero on
   the same "contested sign" principle you already applied to `hydrophobicity`, and
   flagging hydrophobicity for re-evaluation once you have labels.

5. **Your `arbitrary` baseline beating your model is a power artifact, not a finding.**
   Your README already says this correctly. Worth adding: on that pool, `arbitrary` at
   p = 0.064 with 4/20 hits is *one hit short* of the 5/20 critical value. Two
   deterministic orderings differing by one hit is noise at this sample size.

---

## Part D — Concrete recommendations

**Measurement, in priority order.**

1. **Report the critical value with every precision@k claim.** "This pool needs 5/20
   for p ≤ 0.05; we saw 4." Cheap, honest, and makes the design limit legible.
2. **Switch your primary statistical endpoint to AUC01 or AUPRC**, keeping precision@k
   as the business KPI. Precision@k cannot support an inferential claim at P ≈ 14.
3. **Make the baseline comparison paired.** McNemar / exact binomial on discordant
   top-k slots within a patient; paired t-test or paired permutation on per-patient
   deltas across patients. Not Wilcoxon (SIGIR 2019), not bootstrap-shift.
4. **Compute the per-patient metric alongside the pooled one, always.** IMPROVE's TME
   result shows a pooled AUC gain can be entirely between-patient re-sorting. Bootstrap
   over **patients**, not peptides, for CIs.
5. **Measure your top-20 overlap with `binding_only` now.** It sets your power budget:
   ~20 patients at 50% overlap, ~40–45 at 80%.
6. **Add three ablations to `neoantigene benchmark`:** (a) model with `presentation`
   removed; (b) HLA-allele-only predictor as a sanity floor; (c) leave-one-feature-out.
   Report per-allele performance, never a single pooled number.
7. **Add the nested LRT** comparing your full model against the presentation-only model
   fitted on the same labels, plus per-feature 1-df versions.
8. **Fix the `evaluate` k semantics** if you borrow TESLA's metrics; their `ttif` uses
   `rank < 20` (top 19) and a method-dependent denominator.

**Data, in priority order.**

1. **Ingest the NCI Surgery Branch figshare release** (doi 10.35092/yhjc.c.4792338.v2).
   It is the only open resource with minimal peptides, HLA, genomic coordinates,
   negatives at scale, *and* precomputed analogues of seven of your nine features.
   Start with the small `Nmers*.txt` files (17 MB total) to validate your variant
   ingestion end-to-end against real coordinates, then take the `Mmps` test set.
   **Drop `unscreened` rows. Handle the `-`/`CD8` vs `0`/`1` encoding difference between
   train and test. Decide explicitly whether you train on `Screening Status`
   (13,400 rows, model-inferred restriction) or `epitope status` (27 rows, experimental
   restriction) — they are different labels and the metadata does not warn you.**
2. **Ingest the two IMPROVE TSVs** (~470 KB, verified zero overlap between them) as a
   peptide+HLA immunogenicity benchmark. Note the absence of patient IDs; get
   `data.zip` if you need patient-grouped CV.
3. **Pull CEDAR** from `cedar.iedb.org/database_export_v3.php`, and **filter the MHC
   restriction evidence code to exclude `MHC binding prediction`**. Only ~31% of assays
   have four-digit restriction at all, so expect heavy attrition.
4. **Use ITSNdb (`github.com/elmerfer/ITSNdb`) for discrimination only, never for
   precision@k** — its 65% prevalence and 40% A\*02:01 share make it unrepresentative
   of any real pool.
5. **Get TESLA Tables S4 and S7 from the Cell article page in a browser** (I was
   blocked programmatically). They are the only open **multimer-validated** pMHC label
   tables with two-lab independent assays. Remember the public version lacks WT
   sequence and gene name, and every peptide was prediction-preselected.
6. **Treat NEPdb / TSNAdb / dbPepNeo / NDD as positive-example sources at most.** NEPdb
   is 173 positives, not 17,000. TSNAdb is an aggregate of the others and is not
   independent. None supplies trustworthy negatives.
7. **Do not spend on the TESLA DUA or dbGaP phs000452 / phs001003** unless you
   specifically need raw sequencing to re-derive variants. The labels are not there.

**Model.** Demote `agretopicity` toward zero on the same contested-evidence principle
you applied to `hydrophobicity`; flag `hydrophobicity` for re-evaluation given
IMPROVE's contrary finding; record backend training-data provenance in the run manifest;
and regularize `refit` toward your priors rather than fitting 9 free weights on 40
labels.

**Expectation setting.** IMPROVE, trained on 467 positives across 70 patients with
multimer labels and 22 features, reaches **AUC 0.630 / AUC01 0.0139** versus RankEL's
0.539. That is the current state of the art on tumor neoantigens. A model reporting
much above ~0.65 on this task should be assumed leaky until the ablations in §C.5 say
otherwise. The commercially meaningful claim is not a large AUC; it is a **reproducible
per-patient top-20 enrichment over the MHCflurry-sorted baseline, demonstrated on a
pool that was not selected by a predictor**, with the critical value and the overlap
fraction stated.

---

## Verification log

**Verified by direct retrieval and inspection.**

- TESLA: Cell full text (608/37, 310/4, two-lab multimer protocol, Table S4/S7
  captions); Synapse `syn21048999` → `syn23446508` → three folders containing only
  FASTQ/BED, enumerated anonymously via the Sage REST API; `ParkerICI/tesla`
  `performance-metric-functions.R` read in full (AUPRC, FR, TTIF definitions; the
  `rank < 20` off-by-one); Cell `mmc*.xlsx` return HTTP 403 from this environment;
  PMC7652061 hosts no supplements; PMID/PMCID/DOI resolved via NCBI eutils.
- NCI SB: figshare collection 4792338 resolved via the figshare v2 API; four data
  articles and exact file sizes; `NmersTrainingSet.txt` and `NmersTestingSet.txt`
  downloaded in full (30,887 and 8,782 rows; 23 columns enumerated; label
  distributions 139/9,404/21,344 and 46/3,722/5,014); `MmpsTestingSet.txt` header
  range-requested (62 columns enumerated) and all 2.36 GB streamed to count
  2,622,623 rows, 61 alleles, `Screening Status` 13,400/1,092,359/1,516,864,
  `epitope status` 27/1,105,484/1,517,112.
- IMPROVE: both public TSVs downloaded; 17,520 rows (467/17,053) and 2,436 rows
  (548/1,888); peptide-length distribution; per-allele tested/positive counts;
  **zero `(peptide, HLA)` overlap between the two tables** confirmed by set
  intersection; repo contents listed via the GitHub API (`data.zip` 103,626,234 B,
  `neoepitope_tabels/`); full paper text read.
- CEDAR meta-analysis (Sette et al. 2025): full open-access text read — all counts,
  assay percentages, restriction resolution, allele ORs, DAI p = 0.3056, author-overlap
  statistics.
- ITSNdb paper: full text read — inclusion criteria, 129/70, 40.2% A\*02:01, AUC
  0.52–0.60, DAI p = 0.25, TESLA public-version critique, TNB-vs-TMB result.
- ImmuBPI preprint: full text read — intra-HLA imbalance, HLA-only ablation, PRIME2 /
  BigMHC_IM memorization, positions 2 and 9 mutual information.
- Power analysis: computed in `scripts/power_precision_at_k.py`, run in this repo; all
  numbers in §B.4 are its output.

**Reported (read in a primary source, not independently reproduced).**

- NetMHCpan-4.2 test-set sizes and training procedure; ICERFIRE counts and partitioning;
  Parkhurst 2019 screening totals; Hilf 2019 and Keskin 2019 supplementary structure and
  ELISPOT thresholds; Awad 2022, Palmer 2022, Rojas 2023 assay formats; Buckley et al.
  2022 conclusions (abstract and publisher summaries; I did not obtain the full text);
  Gartner 2021 dbGaP accession; NEPdb / TSNAdb / dbPepNeo / NDD counts.

**Unverified — treat with suspicion.**

- **"Caspepdb"** — no primary publication or live server confirmed. Probably CAPED/CAPD
  (`caped.icp.ucl.ac.be`), which is TAA-dominated with no systematic negatives. Do not
  cite on my authority.
- **Weber 2024 / Khattak 2023 (mRNA-4157/V940)** — I found **no** per-peptide
  immunogenicity supplement. Do not assume one exists.
- **Sahin 2017 and Ott 2020** supplementary contents — not inspected.
- **Myeloid leukemia 15,704-pHLA screen** — conference abstract only; no peer-reviewed
  version or data release confirmed. The "89% of experimental binders were not predicted
  to bind" figure is striking but unverified.
- **Miller et al. >130 patients / 25 cancer types** — institutional news claim, not a
  published table. The published study is 13 patients.
- **NDD v0.2 provenance** — no peer-reviewed curation description found.
- TESLA supplementary tables — captions verified, **contents not inspected** (403).
