# The frontier of tumor neoantigen immunogenicity prediction (as of Sept 2026)

A citation-backed technical review, written against the specific question: *does anything
in the current literature reliably beat "sort by predicted presentation/binding" for
ranking neoantigens, and if so, what is it doing that a MHCflurry + 8-feature logistic
model is not?*

The honest headline, stated up front so the rest can be read skeptically:

- **No method robustly and reproducibly beats presentation/binding sorting on a truly
  independent, allele-resolved, experimentally-labeled neoantigen test set.** Several
  methods beat it *on their own cross-validation or on enriched benchmarks*, by modest
  margins, and those margins repeatedly shrink or invert when the test set changes lab,
  assay, or curation source.
- The **best reported precision** on a held-out neoepitope set with experimental labels is
  **BigMHC IM's mean PPVn ≈ 0.44** (Albert et al., *Nat Mach Intell* 2023) — but that number
  is measured against an already-enriched pool whose random baseline is 0.21, is
  retrospective, and draws its positives partly from literature databases with plausible
  train/test leakage. The realistic *prospective* ceiling is closer to the **6% validation
  rate** the TESLA consortium observed among in-silico top-ranked candidates (Wells et al.,
  *Cell* 2020).
- Two 2024 analyses (Wan et al.'s ICERFIRE paper; the ImmuBPI preprint) show that a large
  fraction of the apparent skill of published immunogenicity predictors is **dataset
  bias**: viral-epitope contamination of the positive set, and intra-HLA positive/negative
  ratio imbalance that a model can exploit with the HLA allele alone.

Everything below distinguishes *demonstrated* from *claimed*, gives DOIs/PMIDs and GitHub
URLs, and flags widely-repeated claims that the underlying papers do not support.

---

## 1. PRIME 1.0 and 2.0 (Gfeller lab, UNIL/CHUV Lausanne)

**Primary sources.**
- PRIME 1.0: Schmidt J, Schmidt J, ... Gfeller D. "Prediction of neo-epitope immunogenicity
  reveals TCR recognition determinants and provides insight into immunoediting." *Cell Rep
  Med* 2021;2(2):100194. DOI [10.1016/j.xcrm.2021.100194](https://doi.org/10.1016/j.xcrm.2021.100194).
  PMC7897774.
- PRIME 2.0: Gfeller D, Schmidt J, Croce G, ... Bassani-Sternberg M. "Improved predictions
  of antigen presentation and TCR recognition with MixMHCpred2.2 and PRIME2.0 reveal potent
  SARS-CoV-2 CD8+ T-cell epitopes." *Cell Syst* 2023;14(1):72-83.e5.
  DOI [10.1016/j.cels.2022.12.002](https://doi.org/10.1016/j.cels.2022.12.002). PMID 36603583.
- Code: <https://github.com/GfellerLab/PRIME> (binding predictor
  <https://github.com/GfellerLab/MixMHCpred>). License: academic / non-commercial use
  (the repositories ship a custom academic license; **not** OSI-approved — verify the
  `LICENSE` file before any commercial use).

**The idea it adds over binding.** PRIME does not try to re-predict HLA binding. It takes
MixMHCpred's presentation score (−log(%rank)) as one input and adds a model of **TCR
recognition propensity** built from the amino-acid identities at the peptide positions that
*least* affect HLA binding — the "minimal-impact-on-affinity" (MIA) positions, roughly P4
through PΩ-1, chosen per allele from curated HLA motifs. The premise is that those residues
face the TCR, so their composition carries recognition signal that is (by construction)
nearly decorrelated from affinity. This is the one genuinely load-bearing idea in the PRIME
line, and it is the single feature class your pipeline most clearly lacks (your
`hydrophobicity`/TCR-contact term is a scalar; PRIME learns a 20-dimensional residue
propensity at TCR-facing positions).

**Architecture.**
- *PRIME 1.0*: L2-regularized **logistic regression** (glmnet, α=0) on a 21-dim vector =
  MixMHCpred2.1 −log(%rank) + frequency of each of 20 amino acids at MIA positions.
  Hard rule: peptides with MixMHCpred %rank > 5% are scored 0.
- *PRIME 2.0*: a small **neural network**. Input = MixMHCpred2.2 −log(%rank) (1 node) +
  20 MIA-position amino-acid frequencies + 7-node one-hot peptide length (8–14). The paper
  attributes the NN's gain over the logistic version to its ability to model the
  *interaction* between affinity and TCR-propensity: aromatic/hydrophobic residues at
  TCR-facing positions matter more for *weak* binders than strong ones (Fig. 3E).

**Training data.**
- 1.0: 4,958 peptides (1,282 immunogenic, 3,676 non-immunogenic) pooled from pathogen,
  cancer-testis, and cancer-mutation studies; the *cancer-mutation* subset is only
  **129 immunogenic / 3,200 non-immunogenic**. Plus 2,800 random human 9-mers (50 per each
  of 56 alleles) added as negatives to mimic the real positive:negative ratio. Positives =
  observed CD8+ reactivity (IFN-γ ELISpot or tetramer); negatives = tested-but-no-reactivity.
- 2.0: 596 experimentally verified immunogenic neo-epitopes + 6,084 non-immunogenic tested
  peptides, plus **99 random same-source-protein peptides per neo-epitope** as additional
  negatives to correct the fact that the tested peptides were themselves preselected on
  predicted binding.

**Reported performance.** On 10-fold, leave-one-study-out, and leave-one-allele-out CV,
PRIME beat MixMHCpred, NetMHCpan(EL/BA), MHCflurry, HLAthena, NetMHCstabpan, NetChop, TAP,
and the IEDB (Calis) immunogenicity model on AUC and PR-AUC for the cancer neo-epitope
subset. The authors deliberately used leave-one-study-out/allele-out precisely because
standard CV is inflated by per-study batch effects — a good-practice detail worth copying.
Absolute AUCs on neo-epitopes are modest (mid-0.6s to low-0.7s), not the 0.9s seen for
presentation.

**Documented limitations / criticisms.**
- The "tryptophan enrichment in immunogenic peptides" mechanistic story is **largely a
  viral-data artifact.** Wan et al. (ICERFIRE, §3) reproduced it and showed tryptophan
  importance rises monotonically as you add viral peptides to training and is *not*
  significantly enriched in pure neo-epitopes (W frequency 0.0176 vs 0.0161 in
  immunogenic vs non-immunogenic neo-epitopes, p=0.24). Cite the MIA-position idea; do not
  cite "add tryptophan."
- PRIME's performance drops from ~0.71 to ~0.66 AUC when moved from its own data to CEDAR
  (ICERFIRE Fig. 6/7), i.e. it overfits its training distribution like everything else.
- The cancer-specific positive set is tiny (129 in v1), so much of the model's fitted
  signal is carried by viral/CTA peptides.

---

## 2. BigMHC / BigMHC_IM (Karchin lab, Johns Hopkins)

**Primary source.** Albert BA, Yang Y, Shao XM, Singh D, Smith KN, Anagnostou V, Karchin R.
"Deep neural networks predict class I MHC epitope presentation and transfer learn neoepitope
immunogenicity." *Nat Mach Intell* 2023;5(8):861-872.
DOI [10.1038/s42256-023-00694-6](https://doi.org/10.1038/s42256-023-00694-6). PMID 37829001.
PMC10569228. Code + weights: <https://github.com/KarchinLab/bigmhc> (Zenodo
10.5281/zenodo.8023523). License: the repo ships a `LICENSE` file — it is a source-available
academic license, not confirmed OSI-standard here; **verify before commercial use.**

**Architecture.** Ensemble of **seven pan-allelic deep nets**, ~87M parameters each (~612M
total). Each net = single-headed self-attention → a "wide" bidirectional LSTM (each cell
consumes an 8-residue window of the peptide plus the full MHC pseudosequence) + an "anchor
block" (first/last 4 residues + MHC). Novel MHC pseudosequence: top-30 positions by
information content from an 18,929-sequence cross-species MHC alignment (vs NetMHCpan's
34-mer contact residues). The scalar output is a linear combination of the one-hot MHC
encoding, which gives an interpretable per-position attention overlay.

**The transfer-learning idea.** Because immunogenicity labels are scarce, they first train
**BigMHC EL** on presentation (eluted-ligand) data, then **retrain only the final two
fully-connected layers** on immunogenicity data to get **BigMHC IM**. Because presented ⊃
immunogenic, this is framed as *narrowing* the task, not switching domains.

**Training data.**
- EL: 288,032 eluted ligands + 16.7M random negatives across 149 alleles (from NetMHCpan-4.1
  and MHCflurry-2.0 single-allele sets).
- IM: the *non-random* portion of the PRIME 1.0 + 2.0 datasets = **1,580 positive / 5,293
  negative** (mix of neoepitopes, CTAs, viral). Only the last two layers are updated.

**Reported performance (the numbers that matter for you).**
- Presentation: BigMHC EL AUROC 0.9733 / AUPRC 0.8779 on 45,409 EL among 900,592 decoys,
  beating NetMHCpan-4.1, MixMHCpred-2.1, TransPHLA, MHCnuggets. (Presentation is the solved
  problem; this is not the interesting part.)
- **Neoepitope immunogenicity (198 immunogenic / 739 non-immunogenic, from MANAFEST +
  NEPdb + Neopepsee + TESLA):** BigMHC IM **mean PPVn = 0.4375** (95% CI 0.4108–0.4642),
  vs the best prior method HLAthena-rank 0.2638, and vs BigMHC EL 0.2704. The top 9
  predictions were all immunogenic. **This PPVn ≈ 0.44 is the single best clean-ish
  precision number in the literature** for held-out neoepitopes — but note the pool's
  random baseline is already 0.2113 (the set is enriched), and AUROC on the same set is only
  **0.53–0.57**, i.e. ranking skill across the whole list is weak; the gain is concentrated
  at the very top.
- Infectious-disease epitopes: BigMHC IM PPVn 0.7999 vs PRIME-2.0 0.7991 — **a tie** (the
  paper says so explicitly). So the neoepitope advantage, not a general immunogenicity
  advantage, is the claim.

**Documented limitations (mostly from the authors).**
- All evaluations are **retrospective**; the prospective plan (NSCLC/mesothelioma/GE trials)
  was future work.
- They **could not** show BigMHC IM distinguishes immunogenic from *presented-but-non-
  immunogenic* neoepitopes — no dataset with both MS and immunogenicity labels exists. This
  is the exact discrimination your pipeline needs, and BigMHC does not demonstrate it.
- New alleles require full retraining (the information-content pseudosequence is global).
- Positive neoepitope labels come partly from literature DBs (NEPdb, Neopepsee) that overlap
  the training distribution; the ImmuBPI preprint (§4) reports that ≥6/8 points in a
  mutation-variant test actually appeared in BigMHC IM's / PRIME2's training data, and both
  models then failed to give distinctive predictions across those variants.

---

## 3. ICERFIRE (Nielsen lab, DTU + Peters lab, LJI)

**Primary source.** Wan Y-tR, Koşaloğlu-Yalçın Z, Peters B, Nielsen M. "A large-scale study
of peptide features defining immunogenicity of cancer neo-epitopes." *NAR Cancer*
2024;6(1):zcae002. DOI [10.1093/narcan/zcae002](https://doi.org/10.1093/narcan/zcae002).
PMID 38288446. PMC10823584. Web server (no standalone weights published):
<https://services.healthtech.dtu.dk/services/ICERFIRE-1.0/>. Also runnable inside IEDB's PVC
tool. License: DTU Health Tech academic terms.

This is the most directly relevant paper to your situation, because it is a *feature study*
that trains the same model class you do and reports what actually generalizes.

**Architecture.** Ensemble of **90 random forests** (nested 10-fold CV, 300 trees, max depth
8). Hobohm-1 redundancy reduction between folds to prevent leakage.

**Input / features.** The decisive design choice: operate on the **ICORE** — the nested
submer (length 8–12) with the best NetMHCpan-4.1 presentation %rank — rather than the full
neopeptide. Base input = ICORE amino-acid composition (21-dim) + ICORE %rank. The final
shipped "consensus" model adds only two features on top: **BLOSUM62 mutation score** and
**wild-type antigen expression** (via IEDB pepX over TCGA). Feature importance: %rank
dominates (14.7–22.7%), BLOSUM mutation score ~4–6%.

**How it handled rank vs affinity (your explicit question).** They rebuilt agretopicity on
**presentation %rank instead of binding affinity**, defining a *scaled rank-agretopicity* =
(mutant ICORE rank / WT ICORE rank) scaled by |mutant−WT|. **Rank-agretopicity beat
affinity-agretopicity: AUC 0.589 vs 0.539.** They conclude %rank/eluted-ligand scores are
the better substrate for every binding-derived feature — directly applicable to your
`agretopicity` term, which currently uses an affinity ratio.

**What they found predictive vs not (the skeptical core).**
- Predictive and *generalizing*: ICORE %rank, BLOSUM mutation score, WT expression. These
  went into the shipped model.
- Predictive only *within one dataset* (overfit signals): foreignness score, Boman index,
  physicochemical properties, self-similarity, and various positional-weighting schemes —
  each helped on CEDAR **or** PRIME but not both, and the "optimal CEDAR" model collapsed
  from 0.742 → 0.655 AUC when moved to PRIME.
- **Foreignness (antigen.garnish) ≈ random**: mean AUC 0.516 across three datasets. If you
  are carrying a foreignness/self-dissimilarity term expecting it to rank neoepitopes, this
  is the strongest published evidence that it does not, on pure neo-epitope sets.
- **Tryptophan is a viral artifact** (see PRIME above).
- Anchor-masking (PRIME's scheme) is **dataset-specific**, not universally better.

**Reported performance.** Cross-validation and external AUCs sit ~0.65–0.74. ICERFIRE's
selling point is not a high peak AUC but the *lowest variance across CEDAR/PRIME/NEPdb* —
i.e. it is chosen to generalize rather than to win any single benchmark. It modestly beats
MHC-binding predictors, which themselves score a respectable ~0.64+ AUC and are *more robust
across datasets* than most immunogenicity predictors.

**Documented limitations.** Allele imbalance (HLA-A*02:01 = 38% of training); label
definition noise from collapsing multi-assay IEDB/CEDAR entries; and database curation
errors — they call out "neo-epitopes" in CEDAR that are actually **anchor-optimized synthetic
variants** (SLLMWITQV, ELAGIGILTV) rather than real tumor mutations.

---

## 4. The deep-learning / TCR-co-model cluster (which are real, which to skip)

**Real and worth knowing:**

- **DeepImmuno** (Griffith lab; *Brief Bioinform* 2021; PMC7781330; code
  <https://github.com/griffithlab/DeepImmuno>). (I verified the journal/year and PMC ID;
  I did not independently confirm the author list or issue number.) CNN on a PCA of 566 AAindex physicochemical
  features (to beat one-hot sparsity), trained on IEDB immunogenicity with a **beta-binomial
  "immunogenicity potential"** (responders/tested) as the label, plus dengue/TESLA
  validation. Also ships DeepImmuno-GAN to *generate* immunogenic peptides. Honest caveat
  from the authors' own repo: their GCN variant fails due to "shortcut" learning. In the
  ITSNdb benchmark below, DeepImmuno called **79% of a curated negative set as immunogenic.**

- **Seq2Neo** (Diao K, et al. *Int J Mol Sci* 2022;23(19):11624; PMC9569519; code
  <https://github.com/XSLiuLab/Seq2Neo>, license AFL v3.0). A full raw-reads→neoantigen
  pipeline (SNV/indel/fusion) whose immunogenicity head, **Seq2Neo-CNN**, consumes
  NetMHCpan-4.1 IC50 + NetCTLpan TAP + peptide. Reports beating DeepHLApan, IEDB, and
  DeepImmuno-CNN on TESLA. It is a pipeline convenience more than a modeling advance.

- **NeoaPred** (Jiang D, et al. *Bioinformatics* 2024;40(9):btae547;
  DOI [10.1093/bioinformatics/btae547](https://doi.org/10.1093/bioinformatics/btae547);
  code <https://github.com/Dulab2020/NeoaPred>, **Apache-2.0**). The only structure-based
  entry that is genuinely novel: **PepConf** builds the pHLA-I conformation (82% of test
  structures < 1 Å RMSD, beating PANDORA), then **PepFore** computes a *foreignness score*
  from the surface/structural/atom-group **difference between the mutant and WT peptide
  conformations**. On a 625-positive/2,672-negative test it reports **AUROC 0.81 / AUPRC
  0.54, vs BigMHC 0.70/0.30.** Skeptical reading: this is the strongest structural result,
  but (a) it is a single self-constructed test set, (b) performance craters on alleles
  absent from training (AUROC 0.73→0.62), and (c) HLA-C is poor. Treat 0.81 as an upper
  bound on one distribution, not a portable number. Notably, a *structural* foreignness
  works where *sequence* foreignness (ICERFIRE) does not — that is the interesting signal.

- **DeepNeo / DeepNeo-v2** (Kim K, et al. *Nucleic Acids Res* 2023;51(W1);
  PMID 37070174; web server <https://deepneo.net>). Two nets: pMHC binding + a "T-cell
  reactivity" model trained on structural properties of pMHC pairs. MHC-I and MHC-II.
  Webserver-only; no downloadable weights emphasized.

- **NetTCR / TCR-ESM / ImmuBPI / THLANet / GRAPE** are **TCR-pMHC binding** models, not
  neoantigen-immunogenicity rankers. They require the TCR sequence as input, which you do
  not have in a ranking pipeline. They matter to you only as (a) evidence that pMHC-only
  immunogenicity has a hard ceiling without the TCR, and (b) the source of the sharpest
  bias critique (ImmuBPI, §7). TCR-ESM (Sattari Vayghan? — *Comput Struct Biotechnol J*
  2024; PMC10749252) and GRAPE use ESM-2 embeddings; THLANet
  (*PLoS Comput Biol* 2025;21:e1013050) fine-tunes ESM-2 for TCR-pHLA. Useful background,
  not drop-in rankers.

**Skip / could not verify as a distinct significant neoantigen-immunogenicity method:**
"TRAP" does not correspond to a verifiable, significant neoantigen **immunogenicity**
predictor in the primary literature I could confirm (the acronym is reused across unrelated
tools); I could not verify it, so I am not reporting a phantom method.

---

## 5. NeoDisc / NeoScreen (Bassani-Sternberg & Coukos, UNIL/CHUV Lausanne)

These are **discovery + validation pipelines**, not standalone immunogenicity scores — but
they are the most clinically validated end of the field and their ML ranker is the Müller
2023 model in §6.

**NeoDisc** — Huber F, et al. "A comprehensive proteogenomic pipeline for neoantigen
discovery to advance personalized cancer immunotherapy." *Nat Biotechnol* 2024.
DOI [10.1038/s41587-024-02420-y](https://doi.org/10.1038/s41587-024-02420-y). Singularity
container; web <https://neodisc.unil.ch>. Integrates WES/WGS + RNA-seq + **MS
immunopeptidomics** for HLA-I/II neoantigens, TAAs, viral, and noncanonical antigens, plus
antigen-presentation-machinery (APPM) defect and HLA-LOH analysis. Its HLA-I SNV ranker is
the trained ML model of Müller 2023; HLA-II/TAA/viral use rule-based schemes. The immunogenicity-
relevant features it annotates that yours does not: **NetMHCstabpan stability, NetCTLpan/NetChop
processing, differential agretopicity = log(MT %rank) − log(WT %rank), MixMHCpred-derived
anchor position, ipMSDB immunopeptidome presentation evidence, and CScape/IntOGen oncogenicity.**
Reported benefit: NeoDisc's ML prioritization beat pTuneos, pVACtools, and Gartner et al. on
the NCI-test set (15 samples, 24 immunogenic peptides). Illustrative clinical hit rate
(cervical adenocarcinoma CESC-1): of 66 rule-prioritized short peptides screened by ELISpot,
**11/66 were immunogenic, 2 within the top 10**, and ML re-ranking improved the top-n count.

**NeoScreen** — Arnaud M, et al. "Sensitive identification of neoantigens and cognate TCRs
in human solid tumors." *Nat Biotechnol* 2022;40(5):656-660.
DOI [10.1038/s41587-021-01072-6](https://doi.org/10.1038/s41587-021-01072-6). The
*experimental* half: engineered CD40-activated B-cell/APC presentation to TILs enriches rare
antigen-specific CD8 T cells. Reported effect: **average 3 tumor epitopes/patient with
NeoScreen vs 1 with conventional TIL culture (p=0.02); 19 epitopes across 7 patients, 10 of
which were found *only* with NeoScreen.** The point for you: the field's clinical progress is
coming as much from more sensitive *assays* as from better *predictors* — a reminder that a
weak predictor may be partly an artifact of insensitive validation.

**NeoDiscMS** — Shapiro EI, Huber F, Michaux J, Bassani-Sternberg M. "Sensitive neoantigen
discovery by real-time mutanome-guided immunopeptidomics." *Nat Commun* 2025;16.
DOI [10.1038/s41467-025-62647-4](https://doi.org/10.1038/s41467-025-62647-4). Targeted
real-time MS acquisition guided by NGS-prioritized inclusion lists; improves detection of
mutated neoantigens and TAAs. Direct MS evidence of presentation is the strongest positive
signal there is; NeoDiscMS is about getting more of it per sample.

---

## 6. The most useful paper for you: Müller et al. 2023 (harmonized datasets)

Müller M, Huber F, ... Bassani-Sternberg M. "Machine learning methods and harmonized
datasets improve immunogenic neoantigen prediction." *Immunity* 2023;56(11):2650-2663.
DOI [10.1016/j.immuni.2023.09.002](https://doi.org/10.1016/j.immuni.2023.09.002).
PMID 37816353. Classifiers: <https://figshare.com/s/a000b0990465ab3e9d33> (neo-peptide LR/XGB).

This is the paper whose recipe most directly addresses "my linear model isn't beating a
binding sort." It is a **logistic-regression + XGBoost voting ensemble** — i.e. the same
model class as yours — and it *does* beat binding-affinity ranking, so the gap is features
and training data, not model capacity.

- **Data**: reprocessed WES/RNA-seq from 131 patients across NCI (Gartner et al.), TESLA, and
  in-house HiTIDE, uniformly re-called: 46,017 SNVs / 1,781,445 neo-peptides, of which only
  **212 mutations / 178 neo-peptides** are immunogenic. Trained on **NCI-train** specifically
  because the NCI mini-gene screen has the **least selection bias** (binding was *not* a
  screening criterion there).
- **New predictive features beyond the usual set** (their stated findings):
  1. **HLA presentation "hotspots"** — whether the WT counterpart falls in an
     immunopeptidome-dense region (ipMSDB). Your pipeline has nothing like this.
  2. **Binding promiscuity** — the neo-peptide binds *multiple* of the patient's HLA-I
     alleles (10 NCI immunogenic peptides missed NetMHCpan %rank≤0.5 but passed MixMHCpred —
     the two predictors carry complementary signal, argue for using both).
  3. **Oncogenicity / driver status** of the mutated gene (IntOGen), motivated by the
     oncogenicity–immunogenicity trade-off.
  4. RNA-seq **mutation coverage**, not just gene TPM.
- **Result**: trained on NCI, the classifiers generalized to TESLA and HiTIDE (different
  tumor types, alleles, labs, assays) and **increased immunogenic peptides ranked in the
  top 20 by up to 30%** over Gartner et al., and achieved the top rank among TESLA
  participants on all three metrics.
- **The negative-control lesson you must copy**: a model trained on **HiTIDE** learned to
  weight RNA expression, CCF, and ipMSDB — exactly the features HiTIDE used to *pre-select*
  peptides for screening — and then ranked *worse* on TESLA/NCI. They call this out as
  shortcut learning and train on the least-biased cohort to avoid it. This is the same trap
  your `docs/BENCHMARK.md` already worries about with CU04's pre-selected shortlist.

---

## 7. Protein language models and co-folding (2024–2026)

Verdict: **no PLM or co-folding method has demonstrated a robust, portable advantage for
neoantigen immunogenicity ranking yet.** What exists:

- **ImmugenX** ("A modular protein language modelling approach to immunogenicity
  prediction," *PLoS Comput Biol* 2024;20(10):e1012511,
  DOI [10.1371/journal.pcbi.1012511](https://doi.org/10.1371/journal.pcbi.1012511)). Transformer pMHC submodule trained
  sequentially on pMHC tasks, then an immunogenicity head; extendable with TCR embeddings.
  Reports state-of-the-art pMHC-only immunogenicity on a cancer holdout. **Crucially, the
  authors tried replacing their embedding layer with ESM-2 and found no benefit on current
  datasets** — a candid negative result on PLM transfer for this task, consistent with the
  "not enough labels" ceiling.
- **ESM-2-based TCR-pMHC models** (TCR-ESM, GRAPE, THLANet) improve TCR-epitope *binding*
  prediction, not pMHC-only neoantigen ranking; they need the TCR.
- A dedicated benchmark ("Do Domain-Specific Protein Language Models Outperform General
  Models on Immunology-Related Tasks?") finds that **small ESM-2 variants (8M params, early
  layers) often match or beat larger ones**, and that domain-specific PLMs do not
  reliably win — i.e. the embedding is not the bottleneck; the labels are.
- **Co-folding**: NeoaPred (§4) is the closest to a co-folding immunogenicity method
  (it folds the pMHC and scores the mutant-vs-WT surface delta) and is the most promising
  structural direction, but with the generalization caveats noted. AlphaFold3-class pMHC-TCR
  co-folding for immunogenicity is, as of this writing, aspirational (the authors of these
  papers point to it as future work; I could not verify a published, benchmarked
  AF3-immunogenicity ranker that beats sequence baselines).

---

## 8. Cross-cutting answers

### 8a. Best reported precision@20 / hit rate on a held-out, allele-resolved, labeled set

**BigMHC IM: mean PPVn ≈ 0.4375** on 198 immunogenic / 739 non-immunogenic held-out
neoepitopes (Albert 2023). Top-9 all immunogenic. **Caveats that lower the realistic
ceiling:**
- The pool's random baseline is 0.2113 — it is a pre-enriched set, so 0.44 is ~2× chance,
  not 0.44 out of a raw mutanome.
- Retrospective; positives partly from literature DBs with plausible leakage (ImmuBPI shows
  overlap with BigMHC/PRIME training).
- AUROC on the same set is only 0.53–0.57 — *ranking* skill is weak; the precision is
  concentrated in the extreme top.

The best number that reflects a **prospective, in-silico-top-ranked** reality is the TESLA
consortium's **6%** (37 immunogenic of 608 multimer-tested top candidates; median 3/51 per
team; Wells et al., *Cell* 2020;183:818-834,
DOI [10.1016/j.cell.2020.09.015](https://doi.org/10.1016/j.cell.2020.09.015), PMID 33038342).
Müller 2023's "+30% in the top 20" is a *relative* gain over another method, not an absolute
precision — the absolute number of immunogenic peptides per patient in the top 20 remains
small (single digits).

**Net**: plan your CRO pilot around a realistic top-20 precision of **~15–40% on an enriched
shortlist, ~5–10% on a raw mutanome** — and treat anything above that in a paper as
CV-inflated until shown on an independent lab's assay.

### 8b. Does anything robustly beat "sort by predicted presentation/binding"?

**On curated CV / enriched benchmarks: yes, modestly** (PRIME in leave-one-study/allele-out;
BigMHC IM PPVn 0.27→0.44 after transfer; Müller +30% top-20). **On truly independent,
allele-resolved sets: the advantage is small, fragile, and often bias-driven.** The
skeptical evidence, which you should weight heavily:
- **ICERFIRE**: MHC-binding predictors score ~0.64+ AUC and are *more robust across
  datasets* than most immunogenicity predictors, which overfit; foreignness ≈ random.
- **ImmuBPI** (bioRxiv 2024.02.07.579420,
  DOI [10.1101/2024.02.07.579420](https://doi.org/10.1101/2024.02.07.579420)): an
  immunogenicity model trained on **HLA allele alone** (no peptide) *beat all 15 models* from
  a published benchmark (AUROC +3–16%). The benchmarks are substantially measuring
  intra-HLA positive/negative imbalance, not immunogenicity.
- **ITSNdb / Frontiers** (*Front Immunol* 2023;14:1094236,
  DOI [10.3389/fimmu.2023.1094236](https://doi.org/10.3389/fimmu.2023.1094236)): on a
  clean 297-peptide negative set, MixMHCpred and NetMHCpan had the *lowest* false-positive
  rates (30–31%), while DeepHLApan (100%) and DeepImmuno (79%) were worse — i.e. plain
  binding predictors beat dedicated immunogenicity predictors.
- **Your own benchmark** (`docs/BENCHMARK.md`): `binding_only` ties or beats your model on
  Ott 2017 and CU04, and a deterministic shuffle beat both — consistent with the field.

So: the frontier's honest position is that presentation/binding is a *strong, hard-to-beat
baseline*, and the reproducible increments come from a **small set of orthogonal features**
(below) plus **rigorous de-biasing**, not from a better model family.

### 8c. Features in winning models that your list is missing

Your 8: presentation, clonality/CCF, expression TPM, agretopicity (affinity ratio),
tumor selectivity, self-dissimilarity (BLOSUM62), WT dissimilarity, mutation position,
TCR-contact hydrophobicity.

Consistently-helpful features you do **not** have (ranked by how well they replicate):

1. **TCR-facing residue amino-acid propensity** at non-anchor/MIA positions — a learned
   20-dim composition, not a scalar hydrophobicity. (PRIME's core, echoed by ImmuBPI's
   post-debias finding that hydrophobicity/polarity matter once HLA shortcut is removed.)
   *This is your single biggest omission.*
2. **Rank-based agretopicity** (ratio of MT/WT presentation %rank, not affinity) — ICERFIRE
   AUC 0.589 vs 0.539. Cheap swap to your existing `agretopicity`.
3. **Operate on the ICORE** (best-binding submer) rather than the full peptide for every
   binding-derived feature (ICERFIRE).
4. **Binding stability** (NetMHCstabpan half-life) — used by NeoDisc and Müller.
5. **Binding promiscuity** across the patient's HLA alleles (Müller).
6. **Immunopeptidome presentation hotspot** evidence for the WT source region — MS-backed,
   not predicted (Müller's ipMSDB; NeoDisc). The strongest *positive* signal in the field.
7. **Oncogenicity / driver status** of the mutated gene (IntOGen/CScape) — Müller, NeoDisc.
8. **RNA mutation coverage** (variant-supporting reads), distinct from gene TPM (Müller).
9. **Use two binding predictors** (NetMHCpan *and* MixMHCpred) — they carry complementary
   signal (Müller).
10. **Structural mutant-vs-WT surface delta** (NeoaPred) — promising but not yet portable.

Two of your features are, per the primary sources, **weak or sign-unstable on neo-epitopes**:
`self_dissimilarity`/foreignness (≈ random in ICERFIRE) and `hydrophobicity` (viral artifact;
TESLA found immunogenic pMHC *less* hydrophobic). Your `docs/CITATIONS.md` already sets
hydrophobicity to 0 — that matches the evidence.

### 8d. Why published immunogenicity predictors fail to replicate

1. **IEDB/CEDAR train–test contamination.** Literature-curated positives recur across
   "independent" datasets; ImmuBPI found ≥6/8 mutation-variant test points already in
   PRIME2/BigMHC training. Databases also contain **anchor-optimized synthetic "neo-epitopes"**
   (SLLMWITQV, ELAGIGILTV) that are not real tumor mutations (ICERFIRE).
2. **Positive-set viral bias.** Pooling viral/CTA epitopes with neo-epitopes imports
   features (tryptophan, hydrophobicity) that are viral, not neoantigenic; the IEDB/Calis
   model, trained on viral peptides, fails on neo-epitopes (ICERFIRE §Results/Discussion).
3. **Intra-HLA positive/negative imbalance → HLA shortcut.** Benchmarks with per-allele
   pos/neg ratios that mirror training let a model "predict" immunogenicity from the allele;
   HLA-only beat 15 models (ImmuBPI). This inflates AUROC/AUPRC without any peptide-level
   skill.
4. **Allele imbalance / population bias.** HLA-A*02:01 dominates every dataset (38% in CEDAR);
   "pan-allele" predictors are demonstrably worse on alleles common in
   under-represented populations (HLAEquity, *iScience* 2023, article
   [S2589-0042(23)02690-1](https://www.cell.com/iscience/fulltext/S2589-0042(23)02690-1)).
5. **Selection bias in the negatives.** Negatives were themselves pre-filtered on predicted
   binding, so the model is scored on an artificially hard discrimination and/or learns the
   selection pipeline's features (Müller's HiTIDE cautionary result). TESLA candidates were
   in-silico pre-selected, which is part of why only 6% validated.
6. **Label noise / definition drift.** Collapsing multi-assay IEDB/CEDAR entries into one
   binary label, in-vitro-only readouts that do not reflect in vivo, and immunodominance/
   stochasticity place a hard cap on achievable accuracy (ICERFIRE discussion). A
   science-studies critique (*Soc Stud Sci* 2023,
   DOI [10.1177/03063127231192857](https://doi.org/10.1177/03063127231192857)) makes the same
   point about TESLA's ground truth.
7. **Cross-dataset overfitting**, quantified: optimal-CEDAR 0.742→0.655 on PRIME; PRIME
   0.710→0.659 on CEDAR; an NNAlign model 0.691→0.439 (below random) on NEPdb (ICERFIRE).

---

## 9. Concrete implications for your pipeline

Not a rewrite plan, just what the evidence says is worth trying, ordered by evidence
strength:

- **Add a TCR-facing residue-propensity feature** (PRIME-style MIA-position composition, or
  just call PRIME2 as a feature). This is the best-supported thing you lack.
- **Switch `agretopicity` to rank-based** (MT/WT presentation %rank), and **compute features
  on the ICORE** submer. Cheap, evidence-backed (ICERFIRE).
- **Add MS presentation-hotspot evidence** for the WT source region if you have any
  immunopeptidome reference — the strongest positive signal in the field (Müller/NeoDisc).
- **Add binding stability, promiscuity across alleles, oncogenicity, and RNA coverage**
  (Müller). Use NetMHCpan *and* MixMHCpred.
- **Reconsider `self_dissimilarity`/foreignness and keep hydrophobicity at 0** — both are
  near-random or sign-unstable on pure neo-epitopes (ICERFIRE; TESLA). If you want a
  foreignness that works, the structural surface-delta (NeoaPred) is the only variant with a
  positive published signal.
- **Guard against the shortcut ImmuBPI/Müller describe**: never train or select on a cohort
  whose peptides were pre-filtered by the features you are fitting. Your instinct in
  `docs/BENCHMARK.md` (owning allele-resolved, ex-vivo, full-mutanome labels, ≥10² positives)
  is exactly right and is the actual rate limiter — not the model.
- **Expect the ceiling to be low.** Even the best method (BigMHC IM) is ~0.44 precision on an
  enriched set and ~0.55 AUROC; presentation sorting is a genuinely strong baseline. A
  realistic, defensible goal is a *small* but *reproducible* top-20 precision lift over
  `binding_only`, demonstrated on labels from an assay and cohort your model never touched.

---

### Appendix: method-by-method quick reference

| Method | Model | Immunogenicity training set (pos/neg) | Best reported immuno metric | Code | License |
| --- | --- | --- | --- | --- | --- |
| PRIME 1.0 | Logistic reg (21-dim) | 129/3,200 cancer (of 1,282/3,676 total) | beats binding in LOSO/LOAO CV (neo AUC ~0.65–0.7) | [GfellerLab/PRIME](https://github.com/GfellerLab/PRIME) | academic |
| PRIME 2.0 | Small NN (+length) | 596/6,084 (+99 random/neo) | tie w/ BigMHC IM on viral PPVn 0.80 | same | academic |
| BigMHC IM | 7× BiLSTM+anchor, ~612M | 1,580/5,293 (transfer) | **neo PPVn 0.4375** (AUROC ~0.55) | [KarchinLab/bigmhc](https://github.com/KarchinLab/bigmhc) | academic src-available |
| ICERFIRE | 90× random forest | CEDAR 631/2,402 | robust AUC ~0.65–0.74; rank-agreto 0.589 | [web server](https://services.healthtech.dtu.dk/services/ICERFIRE-1.0/) | DTU academic |
| DeepImmuno | CNN (AAindex PCA) | IEDB beta-binomial | 79% FP on ITSNdb negatives | [griffithlab/DeepImmuno](https://github.com/griffithlab/DeepImmuno) | see repo |
| Seq2Neo-CNN | CNN (IC50+TAP) | in-house + TESLA | beats DeepImmuno/DeepHLApan on TESLA (self-report) | [XSLiuLab/Seq2Neo](https://github.com/XSLiuLab/Seq2Neo) | AFL v3.0 |
| NeoaPred | PepConf+PepFore (structure) | pHLA structures + 625/2,672 | AUROC 0.81/AUPRC 0.54 (one test set) | [Dulab2020/NeoaPred](https://github.com/Dulab2020/NeoaPred) | Apache-2.0 |
| DeepNeo-v2 | 2 nets (bind + reactivity) | IEDB structural | webserver only | [deepneo.net](https://deepneo.net) | webserver |
| Müller 2023 | LR + XGBoost voting | NCI-train (178 immuno neo-pep total) | +30% top-20 vs Gartner; best in TESLA | [figshare](https://figshare.com/s/a000b0990465ab3e9d33) | open |
| NeoDisc | pipeline + Müller ML | — | beats pVACtools/pTuneos on NCI-test | [neodisc.unil.ch](https://neodisc.unil.ch) | academic |

*License entries are best-effort from repositories/services; verify the `LICENSE` file in
each repo before any non-academic use — I did not independently confirm every license text.*
