# Structure-based, TCR-based and PLM-based approaches to neoantigen immunogenicity

**Research report — state of the field as of 30 September 2026**

Scope: whether structural modelling (pMHC and TCR-pMHC), TCR-specificity prediction, protein
language models, or self-tolerance/foreignness modelling add predictive signal for tumour
neoantigen immunogenicity over MHCflurry-2.0-class presentation plus sequence proxies.

Conventions used throughout:

- **Demonstrated** = quantitative result in a primary paper with a held-out or independent test set.
- **Claimed** = reported by the method's own authors, not independently reproduced.
- **Null** = a primary paper reports no significant effect, or near-random performance.
- **Unverified** = I could not find a primary source; stated as such rather than asserted.

---

## 0. Executive summary and the direct answer on the AF3 pMHC secondary filter

**Recommendation: do not implement the "AlphaFold 3 pMHC secondary filter" as a scoring
component. Implement it as a modelling-quality gate only, or leave it unimplemented.**

Four independent reasons, each sourced below:

1. **There is no published evidence that pMHC structure or pMHC confidence metrics correlate
   with immunogenicity.** Every pMHC confidence result I could find validates confidence against
   *model accuracy* (RMSD), not against binding and not against T-cell recognition. The single
   clean head-to-head on *binding* (Motmaen et al., PNAS 2023) found out-of-the-box AlphaFold
   PAE/pLDDT discriminate binders from non-binders "considerably poorer than NetMHCpan", and
   AlphaFold happily docks non-binders into the groove. After task-specific fine-tuning with a
   logistic head on 10,340 pMHC examples, the structural model reaches *parity* with NetMHCpan —
   not superiority — and even that comparison is contaminated because the test peptides are inside
   NetMHCpan's training data.
2. **A pMHC-structure filter is largely redundant with what you already have.** The quantity a
   fine-tuned pMHC structural model recovers is peptide-MHC binding/presentation. MHCflurry 2.0
   presentation already estimates that, at four to five orders of magnitude lower cost.
3. **Licensing blocks it for any commercial use.** AlphaFold 3 source is Apache-2.0 but the
   *weights* are under the AlphaFold 3 Model Parameters Terms of Use (non-commercial, request-only,
   non-redistributable). AlphaFold Server is non-commercial and its Additional Terms of Service
   explicitly forbid using it "in any automated system that predicts the binding or interaction of
   the protein with ligands or peptides", which is exactly what a pMHC secondary filter is.
4. **Cost is prohibitive at candidate-list scale for the structural task that actually has
   evidence** (TCR-pMHC, not pMHC): ~5 days for 1,500 complexes on a single NVIDIA RTX A6000,
   extrapolating to >250 days for a full repertoire-scale dataset (Noakes et al. 2026).

**The one defensible use of the existing demote-only interface.** The demote-only contract is the
right *shape*, but for the wrong reason. The only demonstrated use of pMHC confidence is detecting
peptides whose structure you cannot model: TFold's peptide-core pLDDT correlates with Cα peptide
RMSD at Spearman ρ = 0.55 (discovery) and 0.42 (test). That is enough to flag "this pMHC is
unmodellable / register is ambiguous" and nothing more. So: keep the interface, have it demote only
on *modelling failure*, and do not let low confidence stand in for low immunogenicity.

**What does have signal, in descending order of evidential strength:**

| Approach | Best demonstrated performance | Verdict for your pipeline |
|---|---|---|
| pMHC structure accuracy (TFold, PANDORA, APE-Gen2) | median Cα-pRMSD 0.73–0.86 Å | Solved as a *geometry* problem. Useless as an immunogenicity signal. |
| pMHC structure → binding | AUROC 0.97 (class I), parity with NetMHCpan, only after fine-tuning | Redundant with MHCflurry. Skip. |
| pMHC structure → immunogenicity | No evidence. NeoaPred claims AUROC 0.81 but is unreplicated. | Not actionable. |
| TCR-pMHC (AF3) → unseen-epitope specificity | macro-AUC₀.₁ 0.601 (IMMREP25); ROC-AUC ≈ 0.6 ceiling | Real but weak signal. Needs patient TCRs. Shortlist-only. |
| Sequence TCR-specificity models → unseen epitopes | AUPRC 0.55 best of 50 models; near-random | Do not build. |
| PLM embeddings for 8–11mers | Parity-to-modest gains; generic protein pretraining can *hurt* | Not the bottleneck. Skip. |
| IEDB foreignness | p = 0.24 (IMPROVE); 1.7% feature importance (ICERFIRE) | **Do not add. Consider removing if present.** |
| Agretopicity / DAI | p = 0.96 (IMPROVE); P = 0.25, AUC 0.59 (ITSNdb) | **Do not add.** |
| Self-proteome dissimilarity | AUC 0.85 on original small set; p = 0.24 on broad T-cell-validated set | Dataset-dependent. Treat as unproven. |
| Boring features (mutant %Rank, WT %Rank separately, hydrophobicity/AA composition, TPM, stability, clonality) | ceiling AUC 0.62–0.73 on independent sets | This is where the field actually is. |

**The uncomfortable headline finding:** every rigorous, T-cell-validated neoepitope immunogenicity
benchmark published since 2023 lands in the AUC 0.60–0.73 band regardless of which exotic feature
family is bolted on, and the additional feature families that are *specifically* intended to model
self-tolerance and foreignness are the ones that test null. Section D is the part of this report
you should read most carefully, because it is the part where the literature most strongly
contradicts the prevailing narrative.

---

## A. pMHC class I structure prediction

### A.1 Can current models predict pMHC-I geometry accurately? Yes. This is essentially solved.

**TFold (AlphaFold + hybrid pMHC templates)** — Structure, 2024.
DOI [10.1016/j.str.2023.11.011](https://doi.org/10.1016/j.str.2023.11.011), PMC10872456.

Demonstrated, class I:

- Median Cα peptide RMSD (superimposing on MHC chains only): **0.73 Å** on the discovery set,
  **0.77 Å** on the test set. Median all-atom peptide RMSD 1.55 Å and 1.77 Å respectively.
- Selecting the best model per target rather than the highest-predicted-accuracy model improves
  median Cα-pRMSD only slightly, to 0.64 Å on both sets — i.e. model *selection* is not the
  limiting factor.
- Beats a MODELLER-based pipeline (PANDORA) on 79% of targets, 0.77 Å vs 1.48 Å, Wilcoxon
  p < 10⁻⁷.
- Peptide **register** prediction beats NetMHCpan 4.1 on the discovery set.

Critically for your use case, the paper's own confidence analysis:

- The error score (100 − pLDDT averaged over the peptide core) correlates with Cα-pRMSD at
  Spearman **ρ = 0.55** (discovery) and **ρ = 0.42** (test). That is a *moderate-to-weak*
  correlation, and it is with structural accuracy, not with any biological endpoint.
- Even with pLDDT-based selection, a substantial fraction of models are poor or unacceptable
  (>2.5 Å): 12.5% in the better-predicted half and 31.7% in the worse-predicted half.
- The authors explicitly tried enriching the score (anchor/non-anchor pLDDT weighting, template
  quality terms) and observed **no improvement** in model discrimination.

**PANDORA** — Front. Immunol., 2022. DOI [10.3389/fimmu.2022.878762](https://doi.org/10.3389/fimmu.2022.878762).
Anchor-restrained MODELLER homology modelling; reported ~0.70–0.86 Å median backbone L-RMSD on its
own benchmark. Usable, fast, low resource footprint. Its AF2 comparison is worth knowing: naive
AlphaFold2 multimer and linker approaches misplaced the **P2 anchor residue outside its pocket**,
producing high backbone L-RMSD, and cost up to 18 GB GPU for 20 min per complex. Note this is
*naive* AF2 without pMHC templates; TFold's whole contribution is fixing exactly this by supplying
positionally-aligned peptide templates and dropping MSAs.

**Interpretation.** pMHC-I geometry is a template-lookup problem. Sub-ångström medians are achieved
because the PDB contains close structural neighbours for most alleles. This is a *strength* for
geometry and a *weakness* for inference: a method that succeeds by finding a near-identical template
is not learning anything new about the peptide you care about.

### A.2 Does pMHC structure or confidence predict binding? Only at parity, only after fine-tuning.

**Motmaen, Dauparas, Baek, Abedi, Baker, Bradley — "Peptide-binding specificity prediction using
fine-tuned protein structure prediction networks"**, PNAS 2023.
Code: [github.com/phbradley/alphafold_finetune](https://github.com/phbradley/alphafold_finetune).

This is the single most important paper for your decision, and it contains both the positive and
the negative result.

The negative result first (verbatim paraphrase of their findings):

> PAE (computed between MHC and peptide) and peptide-averaged pLDDT "provided some discrimination of
> binders from non-binders". However, "the discrimination between binders and non-binders was
> considerably poorer than NetMHCpan and AlphaFold tended to dock non-binding peptides in the
> MHC-peptide-binding groove."

The authors attribute this to AlphaFold never having been trained on negative binding examples.
This is the direct answer to "do pLDDT/PAE from an off-the-shelf structure predictor tell you
anything useful about a pMHC pair": they tell you less than a 2020 sequence model.

The positive result: adding a logistic-regression head on the peptide-MHC PAE values, fitting the
head first and then fine-tuning all AlphaFold parameters, trained on 10,340 pMHC examples (203
crystallographic binders, 5,102 modelled binders, 5,035 non-binders; 68 class I alleles, 39 class II
allele pairs):

- Test set: 2,402 binders / 2,717 non-binders, 32 class I and 26 class II alleles, no peptide with
  fewer than two mismatches to training.
- **AUROC 0.97 (class I), 0.93 (class II).** Class I training was restricted to 9-mers; the test set
  included 8-, 9- and 10-mers.
- Binding-affinity regression on HLA-A\*02:01 9-mers: Pearson **r = 0.79**, versus **0.78** for
  NetMHCpan and **0.57** for default AlphaFold.

The authors' own framing is careful and you should adopt it: the fine-tuned model "approaches the
overall performance of the state-of-the-art NetMHCpan sequence-based method", and the direct
comparison is confounded because "the vast majority of available binding data is already contained
in [NetMHCpan's] training set". Their claimed advantage is *generalisation to systems with little
data* — they demonstrate this by transferring the pMHC-optimised model to SH3 and PDZ domains,
which is a genuine and interesting result, and entirely irrelevant to HLA class I where data is
abundant.

**Bottom line for section A.2.** Structural pMHC modelling buys you, at absolute best, parity with
sequence-based binding prediction, and only after you spend the effort to fine-tune on binding data.
You already run MHCflurry 2.0 presentation. There is no reported configuration in which a pMHC
structural score adds orthogonal binding signal on top of a modern presentation predictor for common
HLA class I alleles.

### A.3 Does pMHC structure predict IMMUNOGENICITY? No direct evidence. One unreplicated claim.

I found **no** study demonstrating that pMHC-I structural confidence metrics (pLDDT, PAE, ipTM)
correlate with measured immunogenicity. If such a study exists I could not locate it, and I am
flagging that as a gap rather than as a negative.

The closest positive claim is **NeoaPred** — Bioinformatics 40(9), 2024.
DOI [10.1093/bioinformatics/btae547](https://doi.org/10.1093/bioinformatics/btae547);
code [github.com/Dulab2020/NeoaPred](https://github.com/Dulab2020/NeoaPred) (Apache-2.0).

Claimed:

- Its structure module (PepConf) produces 82.37% of pHLA-I structures at RMSD < 1 Å.
- Its prediction module (PepFore) computes *differences* in surface, spatial-structural and
  atom-group features between mutant and wild-type peptide-HLA complexes, and outputs a learned
  "foreignness score" (note: this is **not** the IEDB foreignness of Łuksza/Richman — same word,
  different quantity).
- On a test set of 625 immunogenic and 2,672 non-immunogenic peptides: **AUROC 0.81, AUPRC 0.54**,
  versus BigMHC at 0.70 / 0.30 and SimToIEDB at 0.64 / 0.40. Ablating the surface, structural or
  atom-group blocks degrades validation performance.
- Performance was "significantly lower" for HLA-C alleles.

Skeptical read: this is a single-group result on a self-assembled split, benchmarked against tools
not designed for that split, in a field where every externally-validated benchmark (Section D.5)
caps out near AUC 0.65. AUPRC 0.54 at a 1:4.3 positive:negative ratio is a real number and would be
notable if it replicated, but it has not been independently replicated, and it was not evaluated in
any of the community benchmarks. Treat as **claimed, unreplicated**. Do not architect around it.

A useful secondary observation from NeoaPred's design: the signal it claims comes from *mutant minus
wild-type* structural deltas, not from absolute confidence. If you ever do revisit structure, that
is the more mechanistically sensible framing than "low pLDDT ⇒ demote".

### A.4 AlphaFold 3, Boltz-1/-2, Chai-1 on pMHC and TCR-pMHC

There is no pMHC-*only* head-to-head of AF3/Boltz/Chai that I could find. All the current
head-to-heads are on the harder TCR-pMHC task, so I report those.

**Pierce lab benchmark** — "Benchmarking AlphaFold and related deep learning approaches for modeling
antibody and TCR antigen recognition", PMC13370930.

- 25 TCR-pMHC test cases (expanded from 20), 20 seeds each (~100 ranked structures/complex).
- **AF3 generally outperforms Boltz-1 and Chai-1** for antibody-protein, antibody-peptide and
  TCR-pMHC. For TCR-pMHC and antibody-peptide, **Boltz-1 > Chai-1**. TCRmodel2 success was
  **comparable to AF3**.
- With a 14-class-I / 6-class-II MHC benchmark: all protocols exceeded 60% top-ranked success for
  medium-or-better accuracy; **high**-accuracy models remained a major challenge. TCRmodel2 with
  massive sampling (1,500 structures/complex) reached 100% top-5 success at medium-or-better and
  >50% top-25 success at high accuracy.
- **Confidence-metric classification accuracy is the key result for you:** for antibody-protein
  interfaces, AF3 confidence scores gave ROC-AUC ≥ 0.94 for classifying model accuracy, but
  "classification accuracy dropped considerably for antibody-peptide and TCR-pMHC interfaces".
  Interface-specific metrics (RL-ipTM, I-pLDDT) outperformed global ipTM and the AlphaFold-Multimer
  model-confidence formulation (0.8·ipTM + 0.2·pTM), because global scores are diluted by the
  peptide-MHC and TCRα-TCRβ interfaces. **If you use any AF confidence metric, use an
  interface-restricted one, never the global model confidence.**

**"A Unified Framework for TCR-pMHC Structural Model Assessment"** (2026 preprint; benchmark of 20
complexes, pre- and post-training-cutoff splits):

- AF3: 90% of complexes above acceptability in both splits; 6 high-quality pre-cutoff, 2 post-cutoff.
- TCRmodel2: comparable acceptability, lower quality (1 HQ pre-cutoff, 0 post-cutoff).
- Chai-1: similar to AF3 pre-cutoff; **collapses post-cutoff** (50% low-quality, 50% acceptable).
- **Boltz-2 and tFold-TCR underperformed in both splits, reaching only acceptable quality with no
  medium- or high-quality models.** They achieved low interface RMSD (<2 Å pre-cutoff, <5 Å
  post-cutoff) but low DockQ — i.e. reasonable local geometry, wrong interface.
- Context: only **232 TCR-pMHC class I complexes in the PDB as of January 2025**. That is the
  binding constraint on this entire subfield.

**Boltz-2** (PMC12262699) on general PDB-2024/2025 holdouts: matches or moderately improves on
Boltz-1, "performs competitively, edging the other commercially available models Chai-1 and
ProteinX, but **lagging a bit behind AlphaFold3**". Boltz-2's headline contribution is affinity
prediction (Pearson 0.66 on an FEP subset) — small-molecule, not peptide-MHC.

**Independent antigen-specificity comparison** (Noakes et al. 2026, below, 1,500 TCR-pMHC complexes,
same MSAs/templates for all models): Chai-1 ≈ AF3 for antigen-specificity prediction, both >
AlphaFold2; **Boltz-2 only comparable to AF2** despite sharing the AF3 architecture
(ANOVA p = 1.90 × 10⁻¹⁹). For structure prediction on 21 holdout complexes AF3 was best
(Kruskal-Wallis p = 2.27 × 10⁻⁵), with the highest fraction of CAPRI "High" quality models.

### A.5 AlphaFold 3 availability and licensing

Verified from primary sources:

- **Source code:** Apache License 2.0
  ([github.com/google-deepmind/alphafold3](https://github.com/google-deepmind/alphafold3)).
- **Model parameters:** *not* open. Released under the
  [AlphaFold 3 Model Parameters Terms of Use](https://github.com/google-deepmind/alphafold3/blob/main/WEIGHTS_TERMS_OF_USE.md).
  Non-commercial, obtained by request from Google DeepMind, not redistributable.
- **AlphaFold Server** (alphafoldserver.com): non-commercial use by individuals, universities,
  non-profits, research institutes, educational and government bodies, or journalism. From the
  Additional Terms of Service, you must not use the Server or its outputs:
  - "in connection with any commercial activities, including research on behalf of commercial
    organizations";
  - "**in any automated system that predicts the binding or interaction of the protein with ligands
    or peptides, such as Glide or AutoDock**";
  - "to train machine learning models or related technology for biomolecular structure prediction
    similar to AlphaFold Server".
  - Server limits: minimum chain length 4 residues (so peptides are supported), 5,000-token overall
    cap, no MSA or template customisation.
- **AlphaFold 2:** Apache-2.0 for both code and weights; usable commercially. EMBL-EBI's own training
  material states plainly that "for some types of projects, you must use AlphaFold 2".
- **Boltz-1 / Boltz-2:** released by the MIT-affiliated Boltz team with training code and pipelines
  available, explicitly positioned against AF3's and Chai-1's restrictive licensing. I did not
  independently verify the current license text in the repository; **confirm before relying on it.**
- **Chai-1:** I could not verify the current commercial status of the Chai-1 weights. **Unverified.**

**Practical consequence.** The second bullet of the Server terms is dispositive: an "AlphaFold 3
pMHC secondary filter" is, on its face, an automated system that predicts peptide binding. Running
it on the Server would breach the terms. Running it on local weights requires the non-commercial
Parameters Terms. If your pipeline is or may become commercial, AF3 is unavailable for this purpose
and AF2 (with the TFold-style template protocol, which is Apache-2.0 via
`phbradley/alphafold_finetune`) or Boltz is the only route.

### A.6 Specific pMHC tools: status

| Tool | Status (Sept 2026) | What it demonstrates | Access |
|---|---|---|---|
| **TFold** | Usable | Median Cα-pRMSD 0.73–0.77 Å class I; best pMHC-I geometry published | Structure 2024, [10.1016/j.str.2023.11.011](https://doi.org/10.1016/j.str.2023.11.011) |
| **PANDORA** | Usable, fast, low-resource | ~0.70–0.86 Å median backbone; anchor-restrained MODELLER | [10.3389/fimmu.2022.878762](https://doi.org/10.3389/fimmu.2022.878762) |
| **APE-Gen2.0** | Usable; webserver + Docker + local | Extends pMHC-I modelling to **PTMs** (phosphorylation, citrullination) and **non-canonical geometries** (N-/C-terminal extensions) — the one genuine capability gap that sequence tools cannot cover | JCIM 64(5):1730–1750, 2024; PMID 38415656; [apegen.kavrakilab.org](https://apegen.kavrakilab.rice.edu/) |
| **DockTope** | **Offline.** IEDB page carries "Docktope is currently not functional (2024-03-26)" | Only ever supported 4 alleles (HLA-A\*02:01, HLA-B\*27:05, H2-Db, H2-Kb); ~1 Å average | [10.1038/srep18413](https://doi.org/10.1038/srep18413) |
| **Rosetta FlexPepDock refinement** | Usable | Sub-ångström pHLA-I refinement, but critically dependent on template choice | [10.1021/ci500393h](https://doi.org/10.1021/ci500393h) |
| **GradDock, DINC2.0, pDOCK** | Usable, legacy | Anchor-based / incremental docking; superseded on accuracy | see APE-Gen2.0 related work |
| **"NetMHCpan-structure hybrids"** | Does not exist as a named tool | The real instantiation of this idea is Motmaen et al.'s fine-tuned AlphaFold (Section A.2), which reaches NetMHCpan parity | [github.com/phbradley/alphafold_finetune](https://github.com/phbradley/alphafold_finetune) |
| **MHCIIFold-GNN** | Class II only | AF3 pMHC-II models + pLDDT + contact-probability features → GNN. Mean FRANK 0.225 fine-tuned; combining with NetMHCIIpan 4.3 gives 0.198, an **11% boost** — the authors' point is *orthogonality*, not superiority | bioRxiv [10.1101/2025.07.10.664203](https://doi.org/10.1101/2025.07.10.664203) |

The MHCIIFold-GNN result is the strongest existing case that structure can be orthogonal to sequence
for peptide-MHC — and note that it is **class II**, where NetMHCIIpan is weaker, core-register
inference is genuinely hard, and the authors generated a new highly-multiplexed binding assay to
train on. None of those conditions hold for your class I problem.

---

## B. TCR-pMHC modelling

### B.1 The tools, with architecture, data, performance and code

**TCRdock** — Bradley, eLife 2023. DOI [10.7554/eLife.82813](https://doi.org/10.7554/eLife.82813).
Code: [github.com/phbradley/TCRdock](https://github.com/phbradley/TCRdock).

- Architecture: AlphaFold (monomer `model_2_ptm` parameters, deliberately *not* Multimer, to reduce
  complex-training bias) with **hybrid multi-chain templates**. Four templates each for pMHC, TCRα,
  TCRβ selected by sequence identity, combined with **12 representative docking geometries** obtained
  by hierarchically clustering experimental TCR-pMHC docking geometries. Three AF runs × 4 templates.
  **No MSA** — single-sequence input, which is what makes it fast.
- Structure benchmark: 130 non-redundant TCR-pMHC complexes. Beats AlphaFold-Multimer full
  (p < 10⁻⁷) and trimmed (p < 10⁻¹²); on the 20 targets unrelated to any pre-May-2018 structure,
  p < 10⁻³ for both. Beats TCRpMHCmodels. Final model beats the best of its 12 template docking
  geometries in 30% of cases and the median template in 94%.
- **Specificity benchmark (the number that matters):** 8 class I epitopes (HLA-A\*02:01 and H2-Db,
  9- and 10-mers) each with up to 50 cognate TCRs, each TCR docked against its true peptide plus
  **9 NetMHCpan-4.1-matched decoy peptides** with binding scores in the same range. Score =
  summed TCR:pMHC residue-residue PAE, then **double-normalised**: subtract each pMHC's mean PAE
  over 50 irrelevant background TCRs (to remove pMHC-intrinsic confidence), then centre each TCR's
  row (to remove TCR-intrinsic confidence).
  - Overall pooled **AUROC = 0.82**.
  - Per-epitope spread is wide: **AUROC ≥ 0.96** for A\*02:01-YLQ9 and A\*02:01-ELA10, down to
    "only slightly better than random" for others.
  - At repertoire level the correct peptide ranks first for **6 of 8** epitopes.
  - Epitopes with more sequence-diverse repertoires (A\*02:01-GLC9, A\*02:01-NLV9) are harder,
    quantified via the TCRdiv diversity measure.

**Three caveats on TCRdock that are decisive for neoantigen ranking.** (i) The normalisation
requires modelling each pMHC against 50 background TCRs, so the per-candidate cost is ~50× the naive
cost. (ii) It requires a *repertoire* of cognate TCRs per epitope — 50 of them. (iii) The task is
inverted relative to yours: it ranks candidate peptides for a known TCR set, which maps onto "I have
TILs, which neoepitope are they against?", not "I have 300 candidates and no TCRs".

**TCRmodel2** — Yin, Ribeiro-Filho, Lin, Gowthaman, Cheung, Pierce, Nucleic Acids Research
51(W1):W569–W576, 2023. DOI [10.1093/nar/gkad356](https://doi.org/10.1093/nar/gkad356).
Server: [tcrmodel.ibbr.umd.edu](https://tcrmodel.ibbr.umd.edu).

- AF2 adapted with a TCR/MHC-focused MSA database and optimised TCR-chain and pMHC template
  selection. 48 benchmark cases (32 class I, 16 class II) deposited after 30 April 2018.
- Medium-or-high CAPRI accuracy for **over 50%** of cases; higher accuracy than AF2 v2.2.
- ~15 min per complex via the server (an independent evaluation measured ~25 min for five ranked
  models). Structure only — no binding classifier.

**Sequence-based TCR specificity models:**

| Model | Architecture | Training data | Reported performance | Code |
|---|---|---|---|---|
| **NetTCR-2.2** | Dual CNN: pan-specific block (32 1-D filters, widths 1/3/5/7/9 over peptide + all 6 CDRs) + peptide-specific block (16 filters); plus loss-scaling and similarity integration | IEDB + VDJdb + 10x; negatives by swapping TCRs to peptides at Levenshtein ≥ 3 | SOTA on IMMREP22; usable accuracy down to **15 positive TCRs per peptide**; authors state challenges "persist for unseen peptides, especially those distant from training examples" | [github.com/mnielLab/NetTCR-2.2](https://github.com/mnielLab/NetTCR-2.2); eLife [10.7554/eLife.93934](https://doi.org/10.7554/eLife.93934) |
| **ERGO-II** | TCR via autoencoder or LSTM; peptide via LSTM; V/J genes, MHC, CD4/CD8 as learned embeddings; two MLP heads (with/without TCRα) | McPAS-TCR, VDJdb | Outperformed by STAPLER by +9.6% mAP on leak-controlled data; drops to near chance under UURA negative sampling (TULIP re-evaluation) | [github.com/mapo9/ERGO-II](https://github.com/mapo9/ERGO-II) |
| **TITAN** | Bimodal attention (`paccmann_predictor`); epitopes optionally as SMILES | BindingDB pretraining + TCR-epitope fine-tuning | Near-random AUPRC (~0.5) on seen epitopes in the Nature Methods 2025 assessment; collapses under strict split | [github.com/PaccMann/TITAN](https://github.com/PaccMann/TITAN) |
| **TULIP** | Encoder-decoder transformer, **unsupervised** — no negative sampling | VDJdb + IEDB + McPAS; split by epitopes with <20 vs >20 examples (1,796 AUCs) | Retains some above-chance predictability on unseen epitopes under the stringent UURA scheme where PanPep, ERGO2 and DLpTCR drop to chance. Performance decays with Levenshtein distance from nearest training epitope | [github.com/barthelemymp/TULIP-TCR](https://github.com/barthelemymp/TULIP-TCR); PNAS [10.1073/pnas.2316401121](https://doi.org/10.1073/pnas.2316401121) |
| **STAPLER** | Transformer on **full-length** TCR + peptide, masked-LM pretrained then fine-tuned | VDJdb+ and DFCI | +9.6% mAP over ERGO-II_LSTM_vj on the leak-controlled VDJdb+-ETN set. **Its own key negative** (below) | bioRxiv [10.1101/2023.04.25.538237](https://doi.org/10.1101/2023.04.25.538237) |
| **pMTnet-omni** | Hybrid sequence + structure; pan-MHC, cross-species, class I and II, handles missing fields by zero-embedding | 111,202 human/mouse class I/II TCR-pMHC pairs from 69 public datasets; train/test from different datasets | Overall AUROC 0.888 / AUPRC 0.786. **Unseen epitopes (158/253 = 62.5% of test): AUROC 0.855.** Unseen MHC 0.827; both unseen 0.782. αβ ablation: masking β drops AUROC to 0.685, masking α to 0.848 | Nat. Commun. 2026, [10.1038/s41467-026-73396-3](https://doi.org/10.1038/s41467-026-73396-3) |
| **MixTCRpred** | Peptide-specific transformers, 146 pMHC models | — | Consistently top-10 for *seen* epitopes across IMMREP23 and the Nature Methods 2025 benchmark | — |
| **TAPIR3** | ESM fine-tuned to **predict Chai-1 structural metrics** (pTM, ipTM, iPAE) — a structure-distillation model with sequence input | IEDB + IMMREP Kaggle set + proprietary (Vcreate) | macro-AUC₀.₁ **0.538** on IMMREP25 unseen peptides — best non-structural submission, and it beat native Chai-1 ranking | TAPIR: bioRxiv [10.1101/2023.09.12.557285](https://doi.org/10.1101/2023.09.12.557285) |
| **ImmSET** | "Immune Synapse Encoding Transformer" — set-based over variable-length sequences | Proprietary MIRA/pairSEQ assay data | Beats AF2 and eventually AF3 on IMMREP25 **A\*02:01** with enough A\*02:01 training data; ~10 ms/pair on a T4 (>10,000× faster than AF pipelines). **Does not transfer to B\*40:01**, where AF pipelines retain a clear advantage | PMLR v297 (noceda26a); arXiv 2603.26994 |
| **UniPMT** | Unified peptide-MHC-TCR framework | — | Claims up to 15% AUPRC improvement on pMHC-TCR, and up to 17.62% AUPRC improvement on experimentally-validated neoantigen-specific binding | Nat. Mach. Intell. 2025, [10.1038/s42256-025-01002-0](https://doi.org/10.1038/s42256-025-01002-0); Zenodo 10.5281/zenodo.14630611 |
| **TRAP** | Contrastive alignment of pMHC structure+sequence features with TCR sequence | — | Claims AUPR 0.84 / AUC 0.92 random split; **AUC 0.75 unseen epitope** | Chem. Sci. 16:9881, 2025, [10.1039/D4SC08141B](https://doi.org/10.1039/D4SC08141B) |
| **TCRBinder** | TCRα/β Roformer encoders (12 blocks, 12 heads, hidden 640) masked-LM pretrained on **5.3M** unpaired full-length TCRs, plus ESM2 encoders for antigen and HLA pseudo-sequence | 5.3M TCRs pretraining | — | PLoS Comput. Biol., [10.1371/journal.pcbi.1014396](https://doi.org/10.1371/journal.pcbi.1014396) |
| **STAG** | Graph neural network on modelled 3-D TCR-pMHC structures, spatial + physicochemical features only | — | "Comparable or better performance than existing methods" | IEEE TCBB, Jan 2025, [10.1109/TCBBIO.2024.3504235](https://doi.org/10.1109/TCBBIO.2024.3504235) |
| **NetTCR-struc** | GNN for docking-quality scoring and model selection on AF-Multimer TCR-pMHC models | — | Spearman(predicted quality, DockQ) 0.681 → **0.855**; avoids selecting failed structures entirely. Honest self-assessment: "the structural pipeline struggled to generate sufficiently accurate TCR-pMHC models for reliable binding classification" | Front. Immunol. 16:1616328, 2025, [10.3389/fimmu.2025.1616328](https://doi.org/10.3389/fimmu.2025.1616328) |
| **"ImmunoBind"** | — | — | **Could not verify.** No publication or repository found under this name. The nearest real methods are NetTCR-struc, STAG, UniPMT, TRAP, TCRBinder, ImmSET and DecoderTCR. If you have a source for ImmunoBind, it did not surface in this search. | — |

Also noted but not verified in a primary source: **OmniTCR** (described in a search synthesis as a
113M-parameter autoregressive model pretrained on 328M immune records). **Unverified.**

### B.2 Is TCR specificity prediction to unseen epitopes solved? No. Emphatically not.

This is the best-documented negative result in the field, replicated across five independent
evaluations with escalating rigour.

**(1) Grazioli et al., "On TCR binding predictors failing to generalize to unseen peptides"**,
Front. Immunol. 2022. DOI [10.3389/fimmu.2022.1014256](https://doi.org/10.3389/fimmu.2022.1014256).

Built the **TChard** dataset (IEDB + VDJdb + McPAS-TCR + MIRA positives; negatives from both
randomisation and 10x assays) and introduced the **hard split**, which guarantees no test peptide
appears in training. Result: in the peptide+CDR3β setting under HS(RN), "the predictions on the test
set barely exceed random-level performance", **AUROC ≈ 0.55**, and this held whether negatives came
from randomisation only or from assays. Conclusion as written: "robust prediction of TCR recognition
is still far for being solved."

**(2) Deng et al., "Performance comparison of TCR-pMHC prediction tools reveals a strong data
dependency"**, Front. Immunol. 2023. DOI [10.3389/fimmu.2023.1128326](https://doi.org/10.3389/fimmu.2023.1128326),
PMC10152969.

Harmonised all major public TCR-pMHC data and retrained five models (TITAN, NetTCR-2.0, ERGO,
DLpTCR, ImRex) under matched conditions. "Model performance collapsed for strict splitting …
indicating that current models do not generalize to unseen peptides." Also: performance is strongly
dependent on data balance and size, i.e. low robustness.

**(3) IMMREP22 and IMMREP23.** IMMREP23: Nielsen et al., ImmunoInformatics 2024.
DOI [10.1016/j.immuno.2024.100045](https://doi.org/10.1016/j.immuno.2024.100045).

53 teams, 398 submissions. Findings, in the organisers' words:

- Reasonable performance in the "seen" pMHC setting (IMMREP22/23 median AUC₀.₁ ≥ 0.7).
- "**Most participating methods had close to random performance on the subset of 'unseen' peptides,
  underlining that this prediction challenge remains essentially unsolved.**"
- A second, self-inflicted lesson: the **negative-set construction in IMMREP23 itself introduced
  data leakage**, identified by several participating teams, which complicated interpretation of
  the rankings. Best-performing on seen pMHC: IMW DETECT (code not available), MixTCRpred, NetTCR.

**(4) Nature Methods 2025 — "Assessment of computational methods in predicting TCR-epitope binding
recognition."** DOI [10.1038/s41592-025-02910-0](https://doi.org/10.1038/s41592-025-02910-0).

The largest systematic assessment: **50 models**, 21 datasets, **762 epitopes**, hundreds of
thousands of binding TCRs; 46 models testable as published, 31 retrainable under standardised
conditions.

Hard numbers:

- **Seen epitopes**, CDR3β-only models with "assay-swapped" negatives (S_Data1: 978 TCRs,
  3 epitopes): best AUPRC **0.70** (ATM-TCR), then TEIM 0.68, TEPCAM 0.67. PiTE-epiSplit, TITAN and
  TCRfinder were near-random (~0.5). Only ATM-TCR showed a decent precision/recall trade-off
  (F1 0.57); TEIM had recall ~0.2; epiTCR and AttnTAP-vdj had recall >0.8 at precision ~0.5.
- **Unseen epitopes** (U_Data1: 345 TCRs, 40 epitopes): best AUPRC **0.55** (ImRex), then ATM-TCR
  and others at **0.52**. **13 of 28 models (46.4%) had AUPRC ≤ 0.50.**
- With patient-source or healthy-source negatives instead of assay-swapped, unseen-epitope AUPRC was
  "near 0.5 for the majority of models", which the authors note "diminished the interpretability of
  relative model rankings".
- Structural conclusions: the **source of negative TCRs substantially determines apparent accuracy**;
  external negatives introduce uncontrolled confounders; performance on independent test sets was
  consistently lower than on intra-dataset test sets across almost all models; models using MHC
  class and αβ information beat CDR3β-only models; performance saturates around **~1,000 TCRs per
  epitope**; optimal positive:negative ratio ≈ 1:1.
- Excluded from this benchmark: TULIP and TCRdock (unsupervised, no negatives) and TCRen (needs
  experimental structures). So this benchmark is silent on the structural approaches.

**(5) IMMREP25 — the unseen-epitope-only contest.** bioRxiv
[10.64898/2026.03.30.715276](https://doi.org/10.64898/2026.03.30.715276). (DOI prefix is as
reported by the source; note it is a preprint.)

Design: 20 unseen peptides restricted by HLA-A\*02:01 or HLA-B\*40:01; 1,000 TCRs, each confirmed
*not* specific for any of the other 19 pMHCs; 1,000 positives + 9,000 in-MHC-context negatives,
downsampled to **8,938 records** so contestants could not assume a uniform TCR-per-pMHC
distribution. Metric: macro-AUC₀.₁ over the 18-pMHC private set. **126 named submissions.**

Leaderboard (top method per team, macro-AUC₀.₁):

| Rank | Team | Type | Modelling | Score used | macro-AUC₀.₁ |
|---|---|---|---|---|---|
| 1 | Bradley | Structural | TCRdock updated to **AF3 / AF3-TD** (encourages canonical TCR-pMHC geometry) + MSAs | **pLDDT** averaged over peptide and CDR loops, row/column mean-normalised | **0.601** |
| 2 | Altin | Structural | AF3 (MSAs + templates on, 5 seeds, 1 diffusion sample, 10 recycles) | mean **PAE** of peptide:TCR interface, features pre-trained on IEDB | 0.582 |
| 3 | Pierce | Structural | AF3 | average inter-chain score | 0.576 |
| 4 | Gowthaman | Structural | AF3 | XGBoost on AF3 metrics only (ipTM, pTM, ranking score, multimer confidence, mean peptide pLDDT); trained on 384 VDJdb-confidence-3 positives vs 404 random pairs. The variant adding 1,280-d ESM-2 embeddings of CDR3α/CDR3β/peptide was **worse** than AF3-metrics-only | — | ~0.57 |
| 5 | Wang | Structural | AF3 | **ipTM** | 0.569 |
| 6 | TAPIR3 | Sequence (structure-distilled) | — | ESM fine-tuned on Chai-1 metrics | 0.538 |
| — | STAPLER | Sequence (pure) | — | — | 0.518 |
| — | EPACT, MINT, TULIP, TEIM, MixTCRpred, NetTCR-2.2 | Sequence | — | — | low range, mostly ≤0.52 |

Distribution: **73.0% of submissions exceeded random (0.50), but 50.0% fell in (0.50, 0.52], and
only 10.3% exceeded 0.55.** Rankings by macro-AUC (ρ = 0.83) and PR-AUC (ρ = 0.92) were largely
concordant.

Four further findings from IMMREP25 that matter:

- **TCR-pMHC-specific scoring functions did not beat AF3's own confidence metrics.** One team ran
  AF3-TD for structure and then scored with **TCRen** (a statistical residue-contact energy
  potential over ≤5 Å contacts); this "significantly reduced performance across all metrics" versus
  AF3's native confidence. The organisers' conclusion: TCR-pMHC-specific knowledge can usefully be
  injected *during* complex prediction (AF3-TD), but "there is not evidence that existing
  TCR-pMHC-specific scoring functions can outperform AF3's confidence metrics".
- **Exploiting dataset structure did not help.** Some competitors clustered TCRs to exploit the
  assumed dataset structure and gained on the 2-pMHC public leaderboard, but across the 18-pMHC
  private set there were "no substantial improvements in macro-AUC₀.₁".
- **The only truly sequence-based method in the top table was STAPLER at 0.518.** TAPIR3 at 0.538
  is a sequence-*input* Chai-1 distillation model, so it carries implicit structural information.
- The organisers' framing: "structural methods are the current leading method for identifying TCRs
  binding to unseen pMHCs … while also falling far short of perfect prediction", and the
  computational cost "would exclude the bulk repertoire scale".

**(6) The most rigorous characterisation of what AF3 is actually doing.** Noakes et al.,
"Characterising AlphaFold 3's ability to predict T cell antigen specificity", bioRxiv
[10.64898/2026.07.08.737208](https://doi.org/10.64898/2026.07.08.737208) (preprint).

This paper is essential reading before you build anything structural, because it is simultaneously
the most optimistic and the most deflationary source available.

- Engineering: a >100× speedup by replacing AF3's MSA stage (the default takes ~51 min/complex for
  MSAs alone) with a precomputed/iterative scheme — verified on 50 complexes (ANOVA n.s. for
  specificity) and 21 holdout structures (Wilcoxon n.s. for DockQ) to be lossless.
- Model comparison on 1,500 TCR-pMHC complexes (750 binders, 750 non-binders), identical MSAs and
  templates: Chai-1 ≈ AF3 > AF2; **Boltz-2 only ≈ AF2**. ANOVA p = 1.90 × 10⁻¹⁹.
- **Ceiling: "The ROC-AUC performance of AF3 reached a maximum of ≈ 0.6. Although this is better
  than previous models, which have near random performance (ROC-AUC ≈ 0.5), it does not achieve
  good enough discrimination of binders and non-binders to remove the need for experimental
  validation or to be considered for use in clinical settings. Performance would need to be upwards
  of ROC-AUC ≈ 0.9 to match the ability of current experimental assays."**
- **Mechanism (the most important finding): AF3's predictive power "is strongly associated with the
  relative placement of the TCR over the pMHC and not chemical interactions between them."** They
  parameterised 251 experimental complexes with 12 geometric features (pitch, tilt, roll, scanning
  angle, X/Y/Z distance, shape complementarity, buried surface area, number and types of
  interactions, contact maps) and fit these against AF3's CDR:pMHC PAE. AF3 is doing geometric
  plausibility, not biophysics.
- **The PAE Aggregator is a partial negative result.** Their learned score, which weights individual
  CDR-residue × pMHC-residue PAE pairs, produced a **statistically significant** reduction in
  intrinsic TCR/pMHC confidence confounding (two-sided t-test, p = 2.15 × 10⁻⁹) — but **"the ability
  of the scores to predict T cell antigen specificity showed no significance after the same
  statistical test was performed."** Better-calibrated, not more predictive.
- **Ablation:** removing paired MSAs, unpaired MSAs and templates all together drops antigen
  specificity to random (ROC-AUC ≈ 0.5) and DockQ to ≈ 0.0. "The sequence information alone is not
  enough for the model." Templates alone give a small significant gain. So AF3's TCR-pMHC signal is
  substantially MSA-derived, i.e. evolutionary-coupling-derived, which is a strange basis for
  predicting hypervariable non-germline CDR3 contacts.
- **Where it does work well:** AF3 clusters sequence-similar TCRs according to binding mode, and
  detects disrupting point mutations accurately — ROC-AUC **0.92** and **0.83** for detecting
  disrupting mutations in two coeliac CD4 epitopes (α1 and ω1). That is mutational-scanning
  interpretation, not de-novo candidate ranking.
- **Cost, measured:** 1,500 complexes ≈ **5 days** on one NVIDIA RTX A6000 (48 GB); the full
  unsampled dataset would have taken **>250 days**.

**Synthesis for B.2.** Unseen-epitope TCR specificity prediction is *not solved*. It has moved from
"indistinguishable from random" (2022–2023) to "reproducibly, significantly better than random but
far from useful" (2025–2026), and the entire gain came from structural modelling with AF3, at ROC-AUC
≈ 0.6 / macro-AUC₀.₁ ≈ 0.60, against a stated clinical-utility bar of ≈ 0.9. Note also that the one
self-reported number that looks much better than this — pMTnet-omni's unseen-epitope AUROC 0.855 —
is not reproduced by any community benchmark, and pMTnet-omni's own authors report "further loss in
AUROC" under more stringent train/test dissimilarity filtering and conclude the model "is yet to
improve". Do not plan against 0.855.

### B.3 Does TCR-repertoire-aware scoring help rank neoantigens?

**Without the patient's TCR repertoire: no. There is nothing to use.**

Every method in Section B.2 requires a TCR as input. Without patient TCRs, the only repertoire-aware
construct available is a *population prior* on T-cell availability, and the only such priors in the
literature are the MHC amplitude / agretopicity term and the IEDB foreignness term — both of which
test null on T-cell-validated benchmarks (Section D.5). I found no validated repertoire-free
"available-T-cell-precursor" score. If you do not have TILs or PBMC TCR-seq, TCR-aware scoring is
not available to you in any evidence-backed form.

**With the patient's repertoire: suggestive positive evidence, all small-n, none definitive.**

- **TIL TCRβ-informed prioritisation** — eLife 2024, DOI [10.7554/eLife.94658](https://doi.org/10.7554/eLife.94658).
  TCRβ sequencing of TILs from **28 colorectal cancer patients**. A combined model (pHLA binding +
  pHLA-TCR binding) versus pHLA-only. At >95% specificity: **PPV 65.9% vs 22.2% and 47.1%**;
  sensitivity 39.7% vs 5.9% and 19.7%; NPV 87.7% vs 82.0% and 84.2%. At >99% specificity: PPV 62.3%
  vs 11.8% and 39.4%. Validated by long-peptide ELISpot in **4 patients**, where the combined model's
  top candidates outperformed NetMHCpan. This is the strongest primary result in this category, and
  the ELISpot validation is n = 4.
- **NEXNEO** (CGTP + epiTCR), conference abstract, **n = 18** CRC patients with ELISpot readout
  (mutant IFN-γ ≥ 2× wild-type plus 4-1BB⁺CD8⁺ expansion). Outperformed NetMHCpan in 15/18 patients;
  median PPV **40% vs 16%**; median ranking coverage 0.39 vs 0.21; captured 50% more immunogenic
  peptides in top-ranked candidates. **But: χ² p = 0.16, not statistically significant.** The authors
  report this honestly; secondary summaries describing it as "significantly enhances" are wrong.
- **VACINUS TCR** — Exp. Mol. Med. 2024, DOI [10.1038/s12276-024-01259-2](https://doi.org/10.1038/s12276-024-01259-2).
  Uses pMTnet v1.0 to score tumour-reactive TIL CDR3β against neoepitope-MHC, splitting candidates
  into tier 1 / non-tier 1. Tier-1 vaccination controlled tumour growth in vivo; non-tier-1 did not.
  The authors' own listed limitations: they do not know why non-tier-1 failed, and "only the publicly
  available pMTnet tool was utilized … A more precise model is needed."
- **DeepPROTECTNeo** — BMC Biology 2026, DOI [10.1186/s12915-026-02725-1](https://doi.org/10.1186/s12915-026-02725-1)
  (preprint: bioRxiv [10.1101/2025.01.04.631301](https://doi.org/10.1101/2025.01.04.631301)).
  End-to-end WES/WGS → variant calling → HLA typing → pMHC binding → variant-driven TCR repertoire
  mining → cross-attention TCR-epitope model. Under a strict TCR-split: mean AUROC **0.7856**,
  AUPRC 0.7932, claimed 4–5% over six SOTA predictors. Recovered **18 of 34** validated
  high-affinity neoepitopes in a patient cohort. Note: TCR-split, *not* epitope-split, so this number
  says nothing about unseen-epitope generalisation, which is the regime neoantigens live in.

**Honest assessment.** Directionally, four independent groups report that adding patient TCR
information improves neoantigen PPV. The effect sizes are large (PPV roughly 2–3×) but the cohorts
are 18–28 patients with 4-patient functional validation, one key result is non-significant, and all
of them are built on TCR-pMHC predictors that the community benchmarks show fail on unseen epitopes
— which every genuine neoantigen is. A recent review (J. Biomed. Sci. 2026,
DOI [10.1186/s12929-026-01286-3](https://doi.org/10.1186/s12929-026-01286-3)) states the position
correctly: "TCR-informed models are most appropriately used as prioritization or
hypothesis-generating tools, with functional validation remaining important."

If you have patient TIL TCR-seq, this is the *one* place in this report where I would say a
prospective experiment is justified — and I would scope it as: take your top ~30 MHCflurry-ranked
candidates, run AF3-TD/TCRdock against the top ~50 expanded TIL clonotypes, and treat the result as
a re-ordering hint with a measured hit rate, not as a filter. Budget ~5 days of A6000 time per 1,500
complexes.

---

## C. Protein language models

### C.1 Do PLM embeddings beat one-hot/BLOSUM for short 8–11mers?

**The concern you raised is substantiated by the cleanest available evidence.**

**USMPep** (BMC Bioinformatics 2021, DOI [10.1186/s12859-021-04289-z](https://doi.org/10.1186/s12859-021-04289-z))
ran the controlled experiment: language-model pretraining on *protein* data versus on
*proteasome-cleaved peptide* data versus training from scratch, all feeding the same MHC-I
downstream regressor (trained on the MHCflurry-2018 dataset, evaluated on IEDB16_I).

Findings, in their words:

- Peptide-domain pretraining "performs best, generally performing slightly better than the
  corresponding model trained from scratch."
- "**In contrast, pretraining on protein data in general even leads to a loss in performance
  compared to training from scratch.**"
- Their stated mechanism is exactly your hypothesis: "the language modeling task on peptide data
  poses additional difficulties compared to language modeling on protein data as the sequences are
  comparably short and the model thus cannot build up a lot of context."
- All differences were small and "mostly remain consistent within error bars".

**Where PLMs do report gains, and why the gains are not what they look like:**

- **Fine-tuning, not embeddings.** Front. Bioinformatics 2023
  (DOI [10.3389/fbinf.2023.1207380](https://doi.org/10.3389/fbinf.2023.1207380)) reports beating
  NetMHCpan on NetMHCpan-4.1's own test set — but the authors state directly that frozen ESM
  embeddings "may not be optimal for the specific MHC task and input format", that MHC
  pseudo-sequences "were not available in the ESM pre-training data", and that they had to fine-tune
  all parameters plus add a graph attention network. The gain is attributable to supervised
  fine-tuning on MHC data, not to the pretrained representation as a feature.
- **ESMCBA** (arXiv [2507.13077](https://arxiv.org/abs/2507.13077);
  [github.com/sermare/ESMCBA](https://github.com/sermare/ESMCBA); weights at
  `huggingface.co/smares/ESMCBA`). ESM Cambrian 300M + continued masked-LM pretraining on
  HLA-associated epitopes (two input formats: epitope alone, epitope ⊕ HLA heavy chain), then
  fine-tuned for IC50. Reported median Spearman across 25 common HLA alleles: **ESMCBA 0.62**,
  NetMHCpan 4.1 0.56, MHCflurry 0.49, HLApollo 0.44, HLAthena 0.37, MHCnuggets 0.22.
  **Three critical qualifications.** (i) The task is *quantitative IC50 regression*, and the authors
  deliberately trained only on high-quality functional affinity data, "avoiding mass spectrometry
  biases that are inherited by existing methods". MHCflurry 2.0 *presentation* — what you use — is
  built on exactly that MS eluted-ligand signal. So this is not evidence that ESMCBA would beat your
  presentation score on your ranking task; it is a different target with a different ground truth,
  and the comparison is structurally unfavourable to MHCflurry. (ii) Continued pretraining helped most
  for alleles with **moderate** data (500–2,000 peptides), improving correlations by ~0.10 — i.e. the
  gain is in the data-poor regime. (iii) Preprint, single group, not independently replicated.
- **TransHLA** (GigaScience 2025): an ESM-2 probe with a predicted-contact-map branch, for the
  allele-free question "is this peptide presentable at all". IEDB test split: 84.72% accuracy /
  91.95 AUC (class I). The ablation is genuinely informative: replacing the ESM-2 sequence embedding
  with a randomly-initialised one drops class I accuracy to **73.37%**, while removing the contact-map
  branch costs ~1 point. So the embedding matters *for this task* — but the task is presentability,
  not immunogenicity, not allele-specific, and it inherits IEDB's allele and pathogen biases.
- **Chemical-fingerprint controls** (two 2025–2026 papers on non-canonical/citrullinated peptides):
  residue-level chemical fingerprints "matched the performance of sequence-based encodings (BLOSUM62
  and one-hot) for canonical peptides", with advantages only for modified residues at anchor
  positions. Peptide-level (position-destroying) fingerprints performed poorly. Reads as: for
  canonical 9-mers, the encoding choice is close to a wash.

**The decisive evidence is from the immunogenicity task itself.** Two large T-cell-validated studies
(details in D.5) show that the bottleneck is nowhere near the encoding:

- **ICERFIRE** (NAR Cancer 2024, DOI [10.1093/narcan/zcae002](https://doi.org/10.1093/narcan/zcae002))
  uses a baseline of **amino-acid composition + %Rank** and reaches AUC 0.719 (CEDAR) / 0.693
  (PRIME). Every additional feature family tested — physico-chemical properties, antigen expression,
  self-similarity, BLOSUM/codon mutation scores, foreignness, agretopicity — moved this by ≤ 0.006
  AUC on CEDAR, and the feature that helped on CEDAR *hurt* on PRIME (0.693 → 0.658, p = 0.048).
- **IMPROVE** (Front. Immunol. 2024, DOI [10.3389/fimmu.2024.1360281](https://doi.org/10.3389/fimmu.2024.1360281)):
  full model AUC 0.630, feature-reduced "simple" model 0.643, CEDAR independent set 0.624 versus
  NetMHCpan RankEL alone at 0.583 (p = 0.009).

Neoepitope immunogenicity prediction sits at **AUC 0.62–0.73 on independent data**. Swapping BLOSUM
for ESM-2 is not going to move that. **Recommendation: do not invest in PLM embeddings for your
8–11mer scoring.** The one place PLM work is worth tracking is not embeddings-as-features but
**structure distillation** — TAPIR3 (ESM fine-tuned to predict Chai-1's pTM/ipTM/iPAE) hit
macro-AUC₀.₁ 0.538 on IMMREP25 unseen peptides at a fraction of AF3's cost, and outperformed native
Chai-1 ranking. That is a real and interesting result. Note the contrary datapoint from the same
contest: the Gowthaman team's XGBoost on AF3 metrics **plus** 1,280-d ESM-2 embeddings of
CDR3α/CDR3β/peptide was *worse* than the AF3-metrics-only variant.

### C.2 2025–2026 foundation models for immunopeptidomics / TCR-pMHC

| Model | What it is | Status | Source |
|---|---|---|---|
| **ESMCBA** | ESM Cambrian 300M, continued MLM on epitopes / epitope⊕HLA, fine-tuned for IC50 | Weights + code public | arXiv [2507.13077](https://arxiv.org/abs/2507.13077) |
| **TCRBinder** | 4 encoders: TCRα/β Roformers masked-LM-pretrained on 5.3M unpaired full-length TCRs (12 blocks, 12 heads, hidden 640, FFN 2560, ~20 h on 4× A100) + ESM2 for antigen and HLA pseudo-sequence | Published | PLoS Comput. Biol., [10.1371/journal.pcbi.1014396](https://doi.org/10.1371/journal.pcbi.1014396) |
| **DecoderTCR** | Two-stage continual pretraining of ESM-2 with component-specific masking (marginal pMHC/TCR data first, then paired interactions), entropy-guided decoding for design | Published | PMC12889700 |
| **ImmSET** | Set-encoding transformer over variable-length sequences; identifies and controls a failure mode that inflated prior sequence-model results; beats fine-tuned ESM2 on the same data | Published | PMLR v297 |
| **UniPMT** | Unified peptide-MHC-TCR; up to +15% AUPRC on pMHC-TCR, +17.62% on neoantigen-specific binding | Published | Nat. Mach. Intell. 2025, [10.1038/s42256-025-01002-0](https://doi.org/10.1038/s42256-025-01002-0) |
| **TAPIR3** | ESM distillation of Chai-1 confidence metrics | IMMREP25 submission, 0.538 macro-AUC₀.₁ | bioRxiv [10.1101/2023.09.12.557285](https://doi.org/10.1101/2023.09.12.557285) (TAPIR v1) |
| **TCR-BERT, TCRLM, IgLM** | Immune-receptor-specific PLMs for repertoire-level tasks (clonotype clustering, motif discovery) | Published | see DecoderTCR related work |
| **ESM-3 for immunopeptidomics** | **I found no peer-reviewed pMHC or immunopeptidomics benchmark using ESM-3.** The only ESM-generation-3-family immunology result I could verify is ESMCBA's use of ESM Cambrian. | **Unverified** | — |
| **OmniTCR** | Described in a search synthesis as 113M autoregressive, pretrained on 328M immune records, for recognition plus conditional TCR generation | **Unverified — no primary source located** | — |

A recurring theme in DecoderTCR's own related-work framing, worth internalising: "standard
pretraining captures marginal regularities within single sequence families rather than conditional
dependencies across multi-component interfaces." That is the structural reason PLMs have not cracked
TCR-pMHC, and it applies with extra force to 9-mers.

---

## D. Self-tolerance and foreignness modelling

This is the section where the gap between the field's narrative and the field's benchmarks is
widest, and it is directly load-bearing for your pipeline's proxy features.

### D.1 The Łuksza / Balachandran neoantigen fitness model — exact model form

**Łuksza et al., "A neoantigen fitness model predicts tumour response to checkpoint blockade
immunotherapy"**, Nature 551:517–520, 2017.
DOI [10.1038/nature24473](https://doi.org/10.1038/nature24473), PMID 29132144.
Companion: **Balachandran et al.**, Nature 551, 2017, DOI [10.1038/nature24462](https://doi.org/10.1038/nature24462).

A neoantigen's fitness cost is the product of two factors:

**(1) MHC amplitude `A`** — the ratio of wild-type to mutant dissociation constants:

```
A = Kd^WT / Kd^MT
```

Mechanistic story: if the wild-type peptide binds MHC poorly (high `Kd^WT`), WT-specific T cells
escaped negative selection, so a peripheral precursor pool exists to recognise the mutant peptide
(low `Kd^MT`). If the WT peptide *is* well-presented, tolerance mechanisms should have depleted
WT-reactive TCRs, and cross-reactivity means mutant-specific TCRs are depleted too. The authors are
explicit that using `A` outperformed using mutant or wild-type affinity alone in their model.

**(2) TCR-recognition probability `R`** — the cross-reactivity partition function:

```
R(s) = (1/Z(k)) · Σ_{e ∈ IEDB} exp[ −k ( a − |s,e| ) ]

Z(k) = 1 + Σ_{e ∈ IEDB} exp[ −k ( a − |s,e| ) ]
```

where `|s,e|` is a **BLOSUM62 local alignment score** between neoantigen `s` and epitope `e`
(blastp for candidate identification, gap open −11, gap extend −1; alignment scores via
Biopython `Bio.pairwise2`), and the epitope set `e` is restricted to **IEDB human infectious-disease
class-I-restricted targets with positive immune assays**.

- `a` = **horizontal displacement** of the binding curve (the sigmoid midpoint).
- `k` = **slope/steepness** at `a`.
- The `+1` in `Z(k)` is the "no-recognition" state, which is what makes `R` a proper probability and
  makes the whole thing a partition function over a two-state (recognised / not recognised) system —
  this is the Boltzmann/thermodynamic analogy that gives the model its name.

Neoantigen quality (2017) = `A × R`. Clone fitness is the weighted effect of the dominant neoantigens
in each subclone. Predicted survival in anti-CTLA-4 melanoma (Van Allen, Snyder cohorts) and
anti-PD-1 NSCLC (Rizvi cohort).

**Łuksza et al., "Neoantigen quality predicts immunoediting in survivors of pancreatic cancer"**,
Nature 606:389–395, 2022. DOI [10.1038/s41586-022-04735-9](https://doi.org/10.1038/s41586-022-04735-9).
Code: [github.com/LukszaLab/NeoantigenEditing](https://github.com/LukszaLab/NeoantigenEditing).

Quality is redefined multiplicatively:

```
Q(p_MT, h) = R(p_MT) × D(p_MT, p_WT, h)

D = log(A) + log(C)          [with relative weight w between the two terms]
```

`R` is unchanged from 2017. The new term is **cross-reactivity distance `C`** — "the antigenic
distance required for T cells to discriminate between `p_MT` and `p_WT`", operationalised as the
**ratio of EC50 values from sigmoidal T-cell activation curves**.

How `C` was measured (this is the empirical core of the 2022 paper):

- Model `p_WT` = HLA-A\*02:01-restricted CMV epitope **NLVPMVATV (NLV)**, with 3 NLV-specific TCRs.
- Every amino acid substituted at every position to generate `p_MT` variants; cross-reactivity
  compared across a 10,000-fold peptide concentration range.
- Substitutions partitioned into highly / moderately / poorly cross-reactive, dependent on both
  position and substituted residue.
- Replicated with a weaker HLA-A\*02:01 `p_WT` from the melanoma self-antigen **gp100** plus 3
  further TCRs, supporting the claim that conserved substitution patterns define `C`.
- Pooled **1,197 TCR-`p_MT` pairs** into a composite `C`. Two factors promote cross-reactivity:
  substitutions at peptide termini, and hydrophobicity.

**Fitted parameters (Extended Data Table 1, Supplementary Methods):**

- `a` = **22.9** (midpoint of the logistic recognition function `R`)
- `k` = **1.5**, *fixed* rather than fitted, explicitly "to limit the number of parameters"
- `w` = **0.22** (relative weight of the two terms in `D`)
- `a` was chosen "to optimize the separation of survival curves", with the cohort split at the
  median value of the score.

The composite fitness model: `F_α = −σ_I · (immune recognition of high-quality neoantigens) + σ_P ·
F_P^α`, where `F_P^α` counts missense mutations in canonical PDAC drivers (KRAS, TP53, CDKN2A,
SMAD4), and `σ_I, σ_P ≥ 0` set the amplitudes of the two components plus a tree-specific
normalisation constant `F_0(T, σ_I, σ_P)`.

**What it actually predicted.** Joint multi-sample phylogenies were reconstructed for all tumours per
patient; the model predicts **the frequencies of clones propagated to recurrent tumours**, i.e. it
identifies immunoedited clones. This is a **clone-level evolutionary** prediction validated on
longitudinal PDAC long-term-survivor samples. It is **not** a per-peptide immunogenicity classifier,
and the 2022 paper does not claim it is one. Any pipeline that lifts `Q = R × D` out of this context
and uses it to rank individual candidate peptides in a single tumour biopsy is using the model
outside its validated regime.

**Three methodological cautions on the model form itself:**

1. **`R` is IEDB-version-dependent by construction.** The authors say so explicitly: "As the peptides
   in IEDB can change over time, we use the current version of IEDB and list the positive epitopes
   used (Supplementary Table 3)." Your `R` values are not comparable to theirs, or to any other
   group's, unless you pin the IEDB snapshot. This is a reproducibility hazard for any pipeline.
2. **`a` was fit to maximise survival-curve separation on the cohort.** That is a legitimate modelling
   choice for a descriptive evolutionary model and a serious overfitting hazard for a transported
   predictive feature.
3. **`k` was fixed, not fitted**, so the sigmoid's steepness — the parameter that controls how sharply
   "similar to a known epitope" converts into "recognised" — was never identified from data.

### D.2 Critiques and replication status of the fitness model

**Direct challenge to the 2017 clinical claim.** "Tumor fitness, immune exhaustion and clinical
outcomes: impact of immune checkpoint inhibitors", bioRxiv
[10.1101/679886](https://doi.org/10.1101/679886) (v3, **preprint — not peer-reviewed**). Applied the
2017 checkpoint-based fitness measures to matched checkpoint-**treatment-naive** TCGA samples where
cytolytic activity confers a known survival benefit, and reported:

- **No significant survival predictive power beyond overall TMB.**
- **No association** between checkpoint-based fitness and tumour T-cell infiltration, cytolytic
  activity (CYT), or TIL burden.
- Probing the model's core assumption in the HBV-infected liver cancer TCGA cohort, they report
  "suggestive evidence that tumor neoepitopes actually dominate viral epitopes in putative
  immunogenicity", which cuts against the viral-similarity premise of `R`.

**Be precise about what this does and does not show.** The cohort is treatment-naive, so it is not a
strict replication attempt of the anti-CTLA-4/anti-PD-1 survival claims. It shows the fitness score
carries no independent prognostic or immune-correlate signal in an independent, untreated cohort. It
is a preprint. I found **no peer-reviewed direct replication attempt** of the 2017 ICB survival
result in an independent ICB cohort, and **no replication attempt at all** of the 2022 PDAC
immunoediting result. That is itself a notable state of affairs for a pair of Nature papers nine and
four years old.

**Erosion of the upstream premise.** The fitness model's clinical claim inherits the assumption that
neoantigen-derived quantities predict ICB response. Two critiques bear on this:

- Genome Medicine 2020, DOI [10.1186/s13073-020-00729-2](https://doi.org/10.1186/s13073-020-00729-2):
  TMB is "a dubious predictor" — no better than chance in RCC, highly sensitive to variant-calling
  methodology, poor patient-level classification.
- eLife reviewed preprint [87465](https://elifesciences.org/reviewed-preprints/87465): assembled the
  largest pan-cancer ICB dataset and found "little evidence that TMB is predictive of response",
  demonstrating that comparable associations arise in **shuffled data** when multiple testing and
  confounding disease subtypes are not controlled.

### D.3 Richman 2019 self-proteome dissimilarity

**Richman, Vonderheide, Rech**, "Neoantigen Dissimilarity to the Self-Proteome Predicts
Immunogenicity and Response to Immune Checkpoint Blockade", Cell Systems 9(4):375–382.e4, 2019.
DOI [10.1016/j.cels.2019.08.009](https://doi.org/10.1016/j.cels.2019.08.009), PMID 31606370.
Tool: **antigen.garnish**, [github.com/immune-health/antigen.garnish](https://github.com/immune-health/antigen.garnish).

Demonstrated (original paper), across five clinical datasets / 318 patients:

- Dissimilarity defined as low sequence alignment of a mutant peptide *and its sub-peptides* to the
  non-mutated reference proteome (BLAST-based).
- On the Chowell et al. immunogenic/non-immunogenic peptide set: **dissimilarity AUC = 0.85**, IEDB
  similarity score AUC = 0.70, mean Kyte-Doolittle hydropathy AUC = 0.70, Atchley factor I AUC =
  0.71, and **ensemble MHC affinity AUC = 0.54** (i.e. barely better than chance).
- High-dissimilarity neoantigens were unique, enriched for hydrophobic sequences, and correlated with
  survival after PD-1 blockade in NSCLC independent of predicted MHC affinity, whereas
  affinity-threshold-defined neoantigen burden correlated poorly.

**Critiques and replication status.** I want to be exact here, because the honest answer is more
nuanced than either "replicated" or "refuted".

- **I found no formal failed-replication paper, retraction, or matched-cohort rebuttal targeting
  Richman 2019 specifically.** Subsequent reviews describe the method and its performance without
  challenging it. If a direct replication failure exists, I did not find it.
- **However, the same underlying quantity tests null on larger T-cell-validated benchmarks.** IMPROVE
  (D.5) evaluated self-similarity (mutant vs normal peptide, kernel distance) across its broad-scale
  T-cell-recognition dataset: **p = 0.24**, not significant — and still not significant when
  stratified into "conserved binders" (mutation outside anchors, WT also presented) versus "improved
  binders" (anchor mutation), which is precisely the stratification the tolerance argument predicts
  should matter.
- **ICERFIRE found self-similarity to be dataset-dependent**: it entered the optimal PRIME model but
  not the optimal CEDAR model, and models optimised on one dataset generalised badly to the other.
- **There is a specific, documented mechanism by which dissimilarity-to-self AUCs get inflated, and it
  applies to the Chowell-style evaluation set.** ICERFIRE showed that adding viral peptides to
  neoepitope training data progressively inflated tryptophan feature importance (up to 8% at 81.8%
  viral data) while **dropping neoepitope AUC to ~60%**, and concluded that the previously-reported
  tryptophan enrichment in immunogenic peptides "might be attributed to the addition of viral data
  rather than a specific enrichment in immunogenic neo-epitopes". Any classifier distinguishing
  "immunogenic" peptides that are substantially viral in origin from "non-immunogenic" peptides that
  are substantially self-derived will score highly on a *dissimilarity-to-human-proteome* feature for
  a trivial reason. Richman's own AUC 0.85 on such a set is consistent with a real tolerance effect
  *and* consistent with this artefact, and the original paper does not distinguish them.

**Verdict on D.3:** self-proteome dissimilarity is **not refuted** but is **unproven as a
transportable feature**. It performs well on curated viral-vs-self evaluation sets and null on
T-cell-validated neoepitope-only sets. Given that your candidates are all neoepitopes, the
neoepitope-only benchmarks are the relevant ones.

### D.4 Thymic / central tolerance modelling: no usable tool exists

**Direct answer: there is no published method that maps thymic proteasome activity or mTEC
expression onto "does this neoepitope have an available T-cell repertoire". I looked specifically and
found nothing. This is a genuine gap, not a hidden literature.**

What exists is theoretical repertoire modelling, which is interesting and *undermines* rather than
supports the amplitude term:

- **"Sparse, random sampling is sufficient for central tolerance"** (bioRxiv
  [10.1101/2025.12.09.693230](https://doi.org/10.1101/2025.12.09.693230), **preprint**) is the most
  quantitatively grounded and the most relevant. Using MHC-abundance estimates, peptide-presentation
  estimates, and single-cell expression from **2,000 human thymic mTECs**, they derive:
  - Each mTEC displays **100–3,000 unique peptides** distributed over roughly **250,000 pMHC slots**.
  - Each thymocyte interacts with roughly **240 mTECs** during a ~2-day medullary dwell.
  - Therefore each T cell encounters **20,000–200,000 unique self peptides** in the thymus.
  - **Random, sparse sampling of ~5% of the unique self peptidome (~30,000 peptides) was sufficient
    to protect most of peripheral self.**

  **Why this matters for your pipeline.** The mechanistic justification for the MHC amplitude term
  `A` (and for agretopicity/DAI generally) is "if the WT peptide is presented, WT-reactive T cells
  were deleted". This result says the thymic self-peptide set each thymocyte actually sees is a small
  *random* subsample. Deletion of any specific WT-peptide-reactive clone is therefore stochastic
  rather than deterministic, and tolerance is maintained collectively rather than per-peptide. That is
  a principled reason to expect a per-peptide `A`/DAI feature to be *noisy in expectation* — which is
  exactly what the benchmarks in D.5 find (DAI p = 0.96).

- Earlier theory, useful for framing, not for features: Yates, "Theories and Quantification of Thymic
  Selection" (Front. Immunol., PMC3912788); "The effects of thymic selection on the range of T cell
  cross-reactivity" (PMC1857316 — simulates 30,000 self peptides across 3 MHC types with positive and
  negative selection thresholds); "Influence of correlated antigen presentation on T-cell negative
  selection in the thymus" (J. R. Soc. Interface 2018,
  DOI [10.1098/rsif.2018.0311](https://doi.org/10.1098/rsif.2018.0311) — finds correlated mTEC
  peptide co-expression patterns matter for negatively selecting low-degeneracy thymocytes, and
  quantifies the escape probability of autoreactive thymocytes).

None of these produce a per-neoepitope score. None use thymoproteasome (β5t / PSMB11) cleavage
specificity, which would be the obvious mechanistic handle and which I found **no** predictive work
on in a neoantigen context.

### D.5 TCR cross-reactivity-to-self as a filter, and the null results on foreignness/agretopicity

**Tools that exist:**

- **ARDitox** — J. Cancer Res. Clin. Oncol. 2025,
  DOI [10.1007/s00432-025-06330-7](https://doi.org/10.1007/s00432-025-06330-7). Predicts off-target
  epitopes (OTEs) for a **given TCR**: proteome-wide search, then filters to predicted binders
  (MHCflurry, ≤2000 nM) with predicted surface presentation (in-house MS-trained model,
  ≥ threshold probability), then assigns a safety score. Designed for TCR-T **safety** assessment,
  not neoantigen ranking, and **requires a TCR as input**. A useful incidental result: 16
  frameshift-derived neoepitopes had **10× fewer** putative OTEs than 16 TAA epitopes (336 vs 3,911),
  and no frameshift OTE had a safety score < 3 — an argument for frameshift neoepitopes over TAAs in
  TCR therapy, contingent on escaping nonsense-mediated decay.
- **NeoPrecis / NeoPrecis-Immuno** — PMC12932759. Refines Łuksza's cross-reactivity distance by
  incorporating MHC-binding motifs, trained on 11,530 MHC-I and 2,610 MHC-II data points. Their
  framing of why CRD is preferable to peptide-only immunogenicity models (PRIME, DeepNeo) is the best
  statement of the tolerance rationale in the recent literature. **But note their own internal null:
  "component contribution analysis showed that the sigmoid scaling factor had minimal impact on
  prediction performance", which they attribute to limited training-set size.** That is an independent
  group finding that the specific sigmoid functional form inherited from Łuksza contributes little.
- **Nothing predicts deletion or anergy of a neoepitope-specific T-cell population.** No tool answers
  "would this neoepitope's T cells have been deleted or rendered anergic". **Not available.**

**Now the central null results.** These are the numbers I would most want you to act on.

**IMPROVE** — Front. Immunol. 15:1360281, 2024.
DOI [10.3389/fimmu.2024.1360281](https://doi.org/10.3389/fimmu.2024.1360281), PMC11021644. Broad-scale
validation of T-cell recognition; feature-by-feature significance testing of immunogenic vs
non-immunogenic neoepitopes:

| Feature | Definition | p-value |
|---|---|---|
| **Foreignness** | Foreignness score, computed with the `antigen.garnish` function | **p = 0.24** |
| **DAI** | Differential agretopicity index (mutant/WT binding ratio) | **p = 0.96** |
| **SelfSim** | Self-similarity, mutant vs normal peptide (kernel distance) | **p = 0.24** |
| Mutation position | Position of mutation in peptide | 10-mer gap p = 0.01 |
| CelPrev | Cellular prevalence (Sequenza / PyClone) | p = 0.016 |
| Expression | Expression level (Kallisto) | p = 0.16 |
| VarAlFreq | Variant allele frequency | p = 0.72 |
| PrioScore | Priority score | p = 0.088 |
| (hydrophobicity / peptide sequence features) | — | p = 2.9 × 10⁻⁹ |

Model performance: cross-validated AUC 0.630 (full) / 0.643 (simple); independent CEDAR AUC 0.624
and AUC₀.₁ 0.0114, versus NetMHCpan RankEL alone at AUC 0.583 / AUC₀.₁ 0.0109 (p = 0.009). So the
whole feature-engineering enterprise buys ~0.04 AUC over a single presentation rank, and the three
tolerance/foreignness features contribute **nothing statistically detectable**, while peptide
sequence features (hydrophobicity, composition) dominate at p = 2.9 × 10⁻⁹.

**ICERFIRE** — NAR Cancer 6(1):zcae002, 2024.
DOI [10.1093/narcan/zcae002](https://doi.org/10.1093/narcan/zcae002).

- Foreignness computed with `antigen.garnish` v2.3.1 (Docker install per the authors' instructions).
- The optimal CEDAR model *did* include foreignness (along with inverted-IC positional weighting,
  Boman index, BLOSUM mutation score, WT %Rank, antigen expression) and beat baseline (p = 0.003).
  **But foreignness accounted for only 1.7% of feature importance, "ranking at the bottom five least
  important features along with the amino acids cysteine, histidine, asparagine and tryptophan", and
  that CEDAR-optimal model "demonstrated a significant performance drop when evaluated on the
  alternate PRIME dataset, suggesting poor power of generalisation."**
- On agretopicity: their rank-based redefinition (ICORE mutant rank / WT ICORE rank, scaled by the
  absolute difference) reached **AUC 0.589**; the classical affinity-ratio agretopicity reached
  **AUC 0.539**. Both are marginal.
- Their overall conclusion on the field's headline features is worth quoting in spirit: the large
  divergence in optimal features and weighting between CEDAR and PRIME "highlight[s] the challenges
  of generalisation and the potential limitations of dataset-specific models".
- They also document severe data pathologies you should know about if you benchmark on CEDAR/IEDB:
  HLA-A\*02:01 alone is 38.35% of peptides and the top 10 alleles are 84.47%; immunogenicity labels
  collapse heterogeneous assay types; and **several peptides annotated as neoepitopes are not**
  (e.g. SLLMWITQV is an anchor-optimised analogue of a NY-ESO-1 peptide, ELAGIGILTV is an
  anchor-optimised variant of EAAGIGILTV — both engineered for binding, not tumour mutations).

**ITSNdb / "Unraveling tumor specific neoantigen immunogenicity prediction"** — Front. Immunol. 2023,
DOI [10.3389/fimmu.2023.1094236](https://doi.org/10.3389/fimmu.2023.1094236). Built a database with
experimentally confirmed cell processing, pMHC-I binding, **and** positive/negative immunogenicity —
i.e. all negatives are confirmed MHC-I binders, which removes the usual confound.

- Best DAI variant: **ratio-DAI anchor position AUC 0.59**; non-anchor 0.56, with non-anchor DAI
  clustered near 1 as expected.
- "**In contrast to previous reports, no significant difference was found between immunogenic and
  non-immunogenic TSNs according to DAI values (P = 0.25).**"
- "Our findings indicate that relative BA over neoantigens with anchor mutations fails in detecting
  immunogenicity, in opposition to previous reports that assume as immunogenic those neoantigens
  whose WT counterpart was predicted as no binder, since no tolerance process would have been carried
  against it." That sentence is a direct rejection of the mechanistic premise behind the amplitude
  term.

**"MHC-I binding affinity derived metrics fail to predict tumor specific neoantigen immunogenicity"**
(same group). Additional findings: wild-type counterpart peptides were themselves frequently
classified as immunogenic; tumour neoantigen burden computed over **WT counterpart** predictions was
*also* associated with ICB response; and TNB failed to beat TMB as an ICB response predictor. Their
conclusion: "**the concept of 'foreignness' is not a predictor of immunogenicity, in agreement with
previous reports concluding that similarity, or not, to viral peptides was poorly informative to
predict/prioritize neoantigen immunogenicity.**" They also identify the root cause: immunogenicity
predictors were trained with positives drawn largely from non-human proteins and negatives drawn by
random sampling of the human proteome, so "they were mainly developed to evaluate foreign peptides
instead of self-derived ones such as neoantigens".

### D.6 TESLA: what the two-axis framing actually showed

**Wells et al. (TESLA Consortium)**, "Key Parameters of Tumor Epitope Immunogenicity Revealed Through
a Consortium Approach Improve Neoantigen Prediction", Cell 183(3):818–834.e13, 2020.
DOI [10.1016/j.cell.2020.09.015](https://doi.org/10.1016/j.cell.2020.09.015), PMID 33038342,
PMC7652061.

Design: 25+ groups each predicted immunogenic epitopes from shared tumour sequencing data; **608
epitopes** were then assayed for T-cell binding in patient-matched samples. Validated in an
independent cohort of **310 epitopes**.

**Presentation axis (three hard thresholds):**

- MHC binding affinity stronger than **34 nM**
- Tumour abundance greater than **33 TPM**
- pMHC binding stability greater than **1.4 h**

**Recognition axis — exactly what was shown.** Among the **29 pMHC that passed all three presentation
criteria (12 immunogenic, 17 non-immunogenic)**:

- Agretopicity (mutant/WT binding-affinity ratio) and foreignness (TCR recognition probability from
  IEDB homology, i.e. Łuksza's `R`) were **independent of the presentation-associated parameters**.
- The majority had agretopicity **> 0.1** and foreignness **< 10⁻¹⁶**. Smaller subsets had values
  orders of magnitude lower (agretopicity, "**group 1**") or higher (foreignness, "**group 2**").
- **The two groups were mutually exclusive (p = 0.005, binomial test).** This is the actual content of
  the "two-axis" framing: low agretopicity and high foreignness identify *disjoint* sets of
  recognised peptides, so they should be combined as a **logical OR**, not multiplied or summed.
- "Recognition" := low agretopicity **OR** high foreignness. Recognised peptides were enriched among
  immunogenic pMHC: **p = 0.003 (Fisher's exact), odds ratio 14.3**, robust to threshold choice.
- Recognition is defined as a property of the peptide alone and is explicitly "not a priori
  associated with the immunogenicity of that peptide in a particular patient".

**Full model:** filtered out **98% of non-immunogenic peptides at precision > 0.70**, retaining
approximately **45% of immunogenic** peptides. Pipelines that prioritised the model's features had
superior performance, and modifying pipelines to use them improved performance.

**The skeptical read you asked for.** The 98% / precision-0.70 headline is a real and useful result,
and most of its work is done by the three *presentation* thresholds — affinity, abundance, stability
— which are uncontroversial and which your pipeline substantially already implements via MHCflurry
presentation plus expression. The **recognition** axis, which is the part that motivates foreignness
and agretopicity features, rests on **n = 29 pMHC (12 vs 17)**. An odds ratio of 14.3 at that sample
size has an extremely wide confidence interval; the mutual-exclusivity p = 0.005 is computed on
subsets of a 29-element set; and the whole thing is conditional on having already passed three hard
filters, so it does not license using agretopicity or foreignness as general-purpose scoring
features across an unfiltered candidate list. Most importantly: **the larger, later,
T-cell-validated benchmarks (IMPROVE n in the thousands, ICERFIRE on CEDAR and PRIME, ITSNdb) do not
reproduce either axis as a significant predictor.** TESLA's recognition finding should be read as a
well-executed hypothesis-generating observation that has not survived scale-up.

---

## E. What I could not verify

Stated explicitly so you do not treat absence of contradiction as confirmation:

- **"ImmunoBind"** — no publication, preprint or repository found under this name. It may be
  misremembered, renamed, or very recent. The real methods in that space are NetTCR-struc, STAG,
  UniPMT, TRAP, TCRBinder, ImmSET and DecoderTCR.
- **OmniTCR** (113M autoregressive, 328M immune records) — appeared only in a search synthesis; no
  primary source located.
- **Chai-1 weight licensing for commercial use** — current status not verified. Check before use.
- **Boltz-1 / Boltz-2 license text** — the papers and commentary position them as open alternatives
  to AF3 and Chai-1 with training code available; I did not read the current repository license.
  Verify before relying on it commercially.
- **Any ESM-3 based immunopeptidomics or pMHC benchmark** — none found. ESMCBA uses ESM Cambrian.
- **Any study correlating pMHC-only confidence metrics (pLDDT/PAE/ipTM) with measured
  immunogenicity** — none found. This is the specific gap that your proposed AF3 pMHC filter would
  be betting on, and it is unfilled.
- **A peer-reviewed direct replication attempt of Łuksza 2017's ICB survival claim in an independent
  ICB cohort, or of Łuksza 2022's PDAC immunoediting result** — none found.
- **A formal failed replication of Richman 2019** — none found. The evidence against it is indirect
  (null results for the same quantity on different, arguably more appropriate, datasets).
- DOIs with the `10.64898/` prefix (IMMREP25; the AF3 characterisation preprint) are reproduced as
  reported by the sources and are preprints; resolve them before citing formally.

---

## F. Concrete recommendations for the pipeline

**1. Do not implement the AF3 pMHC secondary filter as a scoring component.** No evidence links pMHC
confidence to immunogenicity; the only head-to-head on binding puts out-of-box AF confidence *below*
NetMHCpan; a fine-tuned version reaches only parity; AF3 weights and the AlphaFold Server terms both
block the use case commercially, with the Server terms specifically prohibiting automated
peptide-binding prediction.

**2. Repurpose the interface as a modelling-quality gate, keeping the demote-only contract.** The
demote-only design is correct in shape. Make the semantics honest: demote when the pMHC is
*unmodellable* (peptide-core pLDDT below threshold, ambiguous register), because TFold's
pLDDT↔RMSD correlation (ρ = 0.42–0.55) supports exactly that inference and nothing stronger. Log the
confidence value; do not let it act as a proxy for immunogenicity. If you want this cheaply and with
a commercially usable license, use the AF2-based TFold/`alphafold_finetune` template protocol
(Apache-2.0), not AF3.

**3. If you want structure to earn its place, target TCR-pMHC and require patient TCRs.** This is the
only structural application with benchmark support: IMMREP25 macro-AUC₀.₁ 0.601 for AF3-TD/TCRdock
with pLDDT, against near-random for everything sequence-based. Scope it as a shortlist re-ranker:
≤30 top candidates × ≤50 expanded TIL clonotypes. Budget from measured numbers: ~5 days per 1,500
complexes on one A6000. Use interface-restricted confidence (RL-ipTM or I-pLDDT), apply TCRdock's
row/column double-normalisation to remove TCR- and pMHC-intrinsic confidence, and do **not** post-hoc
rescore with a TCR-specific statistical potential — TCRen reduced performance versus AF3's own
confidence in IMMREP25.

**4. Do not add IEDB foreignness or agretopicity/DAI as scoring features, and consider removing them
if present.** On the best available T-cell-validated benchmark (IMPROVE), foreignness p = 0.24, DAI
p = 0.96, self-similarity p = 0.24. On ICERFIRE, foreignness contributes 1.7% of feature importance
and sits in the bottom five. On ITSNdb, DAI shows no significant difference (P = 0.25) with best
AUC 0.59. If you keep foreignness for continuity with published pipelines, pin the IEDB snapshot,
document the `a` and `k` parameters you use, and do not weight it meaningfully.

**5. If you keep any mutant-vs-wild-type construct, use the two ranks as separate features rather
than their ratio.** ICERFIRE found that the classical affinity-ratio agretopicity (AUC 0.539) was
beaten by a rank-based, difference-scaled reformulation (AUC 0.589), and that including the WT %Rank
as its own feature was more useful than the aggregated ratio. Both are marginal, but the
decomposition is strictly better than the ratio.

**6. If you use TESLA's framing, combine agretopicity and foreignness as a logical OR, never as a
product or sum** — they identified mutually exclusive peptide sets (p = 0.005). And treat that
finding as n = 29, conditional on already passing affinity < 34 nM, abundance > 33 TPM and
stability > 1.4 h.

**7. Spend effort where the benchmarks say the signal is.** IMPROVE, ICERFIRE and TESLA all converge
on the same unglamorous list: mutant %Rank / presentation score, WT %Rank as a separate feature,
peptide hydrophobicity and amino-acid composition (IMPROVE: p = 2.9 × 10⁻⁹, the strongest single
effect in that study), antigen expression (TPM), pMHC binding stability, and clonality / cellular
prevalence (p = 0.016). Also worth implementing: **APE-Gen2.0** for the one genuine capability gap
structure fills — post-translationally modified and non-canonically-extended peptides, which no
sequence predictor covers.

**8. Do not invest in PLM embeddings for 8–11mer scoring.** USMPep's controlled comparison shows
generic protein pretraining *degrades* MHC-I downstream performance relative to training from
scratch; peptide-domain pretraining gives marginal, within-error-bar gains; chemical fingerprints
merely match BLOSUM62 for canonical peptides. The reported PLM wins come from full supervised
fine-tuning on MHC data, not from embeddings-as-features, and ESMCBA's headline gains are on
quantitative IC50 regression with MS data deliberately excluded — a different target from your
presentation score. The one PLM direction worth watching is **structure distillation** (TAPIR3,
0.538 on IMMREP25 unseen peptides at a small fraction of AF3's cost), not representation transfer.

**9. Calibrate expectations for the whole pipeline.** The independent-test-set ceiling for neoepitope
immunogenicity prediction is AUC ≈ 0.62–0.73. A single presentation rank gets you to ≈ 0.58. Nothing
in the structural, TCR or PLM literature as of September 2026 changes that ceiling for the
no-patient-TCR setting. Any design that assumes a large step change is assuming a result that has
not been published.

---

## G. Reference list

**pMHC structure**

1. Accurate modeling of peptide-MHC structures with AlphaFold (TFold). *Structure*, 2024. [10.1016/j.str.2023.11.011](https://doi.org/10.1016/j.str.2023.11.011); PMC10872456
2. PANDORA: A Fast, Anchor-Restrained Modelling Protocol for Peptide:MHC Complexes. *Front. Immunol.*, 2022. [10.3389/fimmu.2022.878762](https://doi.org/10.3389/fimmu.2022.878762)
3. Motmaen A, Dauparas J, Baek M, Abedi MH, Baker D, Bradley P. Peptide-binding specificity prediction using fine-tuned protein structure prediction networks. *PNAS*, 2023. Code: [github.com/phbradley/alphafold_finetune](https://github.com/phbradley/alphafold_finetune)
4. APE-Gen2.0. *J. Chem. Inf. Model.* 64(5):1730–1750, 2024. PMID 38415656; PMC10936522; [apegen.kavrakilab.org](https://apegen.kavrakilab.rice.edu/)
5. DockTope. *Sci. Rep.* 5:18413, 2015. [10.1038/srep18413](https://doi.org/10.1038/srep18413). **Server offline since 2024-03-26.**
6. NeoaPred. *Bioinformatics* 40(9):btae547, 2024. [10.1093/bioinformatics/btae547](https://doi.org/10.1093/bioinformatics/btae547); [github.com/Dulab2020/NeoaPred](https://github.com/Dulab2020/NeoaPred)
7. MHCIIFold-GNN ("Learned Geometry, Predicted Binding"). bioRxiv [10.1101/2025.07.10.664203](https://doi.org/10.1101/2025.07.10.664203)

**Structure prediction models and licensing**

8. AlphaFold 3 repository and Model Parameters Terms of Use. [github.com/google-deepmind/alphafold3](https://github.com/google-deepmind/alphafold3)
9. AlphaFold Server Additional Terms of Service. [gstatic.com/alphafoldserver/.../AlphaFold-Server-Additional-Terms-of-Service.pdf](https://www.gstatic.com/alphafoldserver/app/app/routes/terms_text/AlphaFold-Server-Additional-Terms-of-Service.pdf)
10. Abramson J et al. Accurate structure prediction of biomolecular interactions with AlphaFold 3. *Nature* 630:493–500, 2024
11. Boltz-1. bioRxiv 2024.11.19.624167
12. Boltz-2. PMC12262699
13. Chai-1. bioRxiv 2024.10.10.615955
14. Benchmarking AlphaFold and related deep learning approaches for modeling antibody and TCR antigen recognition. PMC13370930
15. A Unified Framework for TCR-pMHC Structural Model Assessment. 2026 preprint

**TCR-pMHC**

16. Bradley P. Structure-based prediction of T cell receptor:peptide-MHC interactions. *eLife* 12:e82813, 2023. [10.7554/eLife.82813](https://doi.org/10.7554/eLife.82813); [github.com/phbradley/TCRdock](https://github.com/phbradley/TCRdock)
17. Yin R et al. TCRmodel2. *Nucleic Acids Res.* 51(W1):W569–W576, 2023. [10.1093/nar/gkad356](https://doi.org/10.1093/nar/gkad356); [tcrmodel.ibbr.umd.edu](https://tcrmodel.ibbr.umd.edu)
18. NetTCR-2.2. *eLife* 2024. [10.7554/eLife.93934](https://doi.org/10.7554/eLife.93934); [github.com/mnielLab/NetTCR-2.2](https://github.com/mnielLab/NetTCR-2.2)
19. ERGO-II. [github.com/mapo9/ERGO-II](https://github.com/mapo9/ERGO-II)
20. TITAN. [github.com/PaccMann/TITAN](https://github.com/PaccMann/TITAN)
21. TULIP. *PNAS* 121, 2024. [10.1073/pnas.2316401121](https://doi.org/10.1073/pnas.2316401121); PMC11181096; [github.com/barthelemymp/TULIP-TCR](https://github.com/barthelemymp/TULIP-TCR)
22. STAPLER. bioRxiv [10.1101/2023.04.25.538237](https://doi.org/10.1101/2023.04.25.538237)
23. pMTnet-omni. *Nat. Commun.* 2026. [10.1038/s41467-026-73396-3](https://doi.org/10.1038/s41467-026-73396-3); preprint bioRxiv [10.1101/2023.12.01.569599](https://doi.org/10.1101/2023.12.01.569599)
24. UniPMT. *Nat. Mach. Intell.* 2025. [10.1038/s42256-025-01002-0](https://doi.org/10.1038/s42256-025-01002-0)
25. TRAP. *Chem. Sci.* 16:9881, 2025. [10.1039/D4SC08141B](https://doi.org/10.1039/D4SC08141B)
26. TCRBinder. *PLoS Comput. Biol.* [10.1371/journal.pcbi.1014396](https://doi.org/10.1371/journal.pcbi.1014396)
27. STAG. *IEEE TCBB*, Jan 2025. [10.1109/TCBBIO.2024.3504235](https://doi.org/10.1109/TCBBIO.2024.3504235)
28. Deleuran SN, Nielsen M. NetTCR-struc. *Front. Immunol.* 16:1616328, 2025. [10.3389/fimmu.2025.1616328](https://doi.org/10.3389/fimmu.2025.1616328)
29. ImmSET. *PMLR* v297; arXiv 2603.26994
30. DecoderTCR. PMC12889700

**Generalisation benchmarks (the negative results)**

31. Grazioli F et al. On TCR binding predictors failing to generalize to unseen peptides. *Front. Immunol.* 13:1014256, 2022. [10.3389/fimmu.2022.1014256](https://doi.org/10.3389/fimmu.2022.1014256)
32. Deng L et al. Performance comparison of TCR-pMHC prediction tools reveals a strong data dependency. *Front. Immunol.* 14:1128326, 2023. [10.3389/fimmu.2023.1128326](https://doi.org/10.3389/fimmu.2023.1128326); PMC10152969
33. Nielsen M et al. Lessons learned from the IMMREP23 TCR-epitope prediction challenge. *ImmunoInformatics*, 2024. [10.1016/j.immuno.2024.100045](https://doi.org/10.1016/j.immuno.2024.100045)
34. Assessment of computational methods in predicting TCR-epitope binding recognition. *Nat. Methods*, 2025. [10.1038/s41592-025-02910-0](https://doi.org/10.1038/s41592-025-02910-0)
35. IMMREP25: Unseen Peptides. bioRxiv [10.64898/2026.03.30.715276](https://doi.org/10.64898/2026.03.30.715276) (preprint)
36. Noakes et al. Characterising AlphaFold 3's ability to predict T cell antigen specificity. bioRxiv [10.64898/2026.07.08.737208](https://doi.org/10.64898/2026.07.08.737208) (preprint)

**TCR-repertoire-aware neoantigen prioritisation**

37. The T cell receptor β chain repertoire of tumor infiltrating lymphocytes improves neoantigen prediction and prioritization. *eLife* 2024. [10.7554/eLife.94658](https://doi.org/10.7554/eLife.94658)
38. VACINUS. *Exp. Mol. Med.* 2024. [10.1038/s12276-024-01259-2](https://doi.org/10.1038/s12276-024-01259-2)
39. DeepPROTECTNeo. *BMC Biol.* 2026. [10.1186/s12915-026-02725-1](https://doi.org/10.1186/s12915-026-02725-1)
40. NEXNEO (209RO). Conference abstract; n = 18; PPV 40% vs 16%, p = 0.16
41. Artificial intelligence for translational personalized neoantigen cancer vaccine development. *J. Biomed. Sci.* 2026. [10.1186/s12929-026-01286-3](https://doi.org/10.1186/s12929-026-01286-3)

**Protein language models**

42. USMPep: universal sequence models for MHC binding affinity prediction. *BMC Bioinformatics*, 2021. [10.1186/s12859-021-04289-z](https://doi.org/10.1186/s12859-021-04289-z)
43. Improved prediction of MHC-peptide binding using protein language models. *Front. Bioinform.* 2023. [10.3389/fbinf.2023.1207380](https://doi.org/10.3389/fbinf.2023.1207380)
44. Mares S. Continued domain-specific pre-training of protein language models for pMHC-I binding prediction (ESMCBA). arXiv [2507.13077](https://arxiv.org/abs/2507.13077); [github.com/sermare/ESMCBA](https://github.com/sermare/ESMCBA)
45. TransHLA. *GigaScience*, 2025
46. NetMHCpan 4.0 / 4.1. PMC5679736; [10.1093/nar/gkaa379](https://doi.org/10.1093/nar/gkaa379)
47. BigMHC. *Nat. Mach. Intell.*, 2023. [10.1038/s42256-023-00694-6](https://doi.org/10.1038/s42256-023-00694-6)

**Self-tolerance, foreignness, and the immunogenicity benchmarks**

48. Łuksza M et al. A neoantigen fitness model predicts tumour response to checkpoint blockade immunotherapy. *Nature* 551:517–520, 2017. [10.1038/nature24473](https://doi.org/10.1038/nature24473); PMID 29132144
49. Balachandran VP et al. *Nature* 551, 2017. [10.1038/nature24462](https://doi.org/10.1038/nature24462)
50. Łuksza M et al. Neoantigen quality predicts immunoediting in survivors of pancreatic cancer. *Nature* 606:389–395, 2022. [10.1038/s41586-022-04735-9](https://doi.org/10.1038/s41586-022-04735-9); code [github.com/LukszaLab/NeoantigenEditing](https://github.com/LukszaLab/NeoantigenEditing)
51. Tumor fitness, immune exhaustion and clinical outcomes: impact of immune checkpoint inhibitors. bioRxiv [10.1101/679886](https://doi.org/10.1101/679886) (**preprint**)
52. Richman LP, Vonderheide RH, Rech AJ. Neoantigen Dissimilarity to the Self-Proteome Predicts Immunogenicity and Response to Immune Checkpoint Blockade. *Cell Syst.* 9(4):375–382.e4, 2019. [10.1016/j.cels.2019.08.009](https://doi.org/10.1016/j.cels.2019.08.009); PMID 31606370; [github.com/immune-health/antigen.garnish](https://github.com/immune-health/antigen.garnish)
53. Wells DK et al. (TESLA Consortium). Key Parameters of Tumor Epitope Immunogenicity Revealed Through a Consortium Approach Improve Neoantigen Prediction. *Cell* 183(3):818–834.e13, 2020. [10.1016/j.cell.2020.09.015](https://doi.org/10.1016/j.cell.2020.09.015); PMID 33038342; PMC7652061
54. IMPROVE: a feature model to predict neoepitope immunogenicity through broad-scale validation of T-cell recognition. *Front. Immunol.* 15:1360281, 2024. [10.3389/fimmu.2024.1360281](https://doi.org/10.3389/fimmu.2024.1360281); PMC11021644
55. ICERFIRE / A large-scale study of peptide features defining immunogenicity of cancer neo-epitopes. *NAR Cancer* 6(1):zcae002, 2024. [10.1093/narcan/zcae002](https://doi.org/10.1093/narcan/zcae002)
56. Unraveling tumor specific neoantigen immunogenicity prediction: a comprehensive analysis (ITSNdb). *Front. Immunol.* 14:1094236, 2023. [10.3389/fimmu.2023.1094236](https://doi.org/10.3389/fimmu.2023.1094236)
57. MHC-I binding affinity derived metrics fail to predict tumor specific neoantigen immunogenicity
58. NeoPrecis. PMC12932759
59. ARDitox. *J. Cancer Res. Clin. Oncol.* 2025. [10.1007/s00432-025-06330-7](https://doi.org/10.1007/s00432-025-06330-7)
60. Burden of tumor mutations, neoepitopes, and other variants are weak predictors of cancer immunotherapy response and overall survival. *Genome Med.* 2020. [10.1186/s13073-020-00729-2](https://doi.org/10.1186/s13073-020-00729-2)
61. Is tumor mutational burden predictive of response to immunotherapy? *eLife* reviewed preprint [87465](https://elifesciences.org/reviewed-preprints/87465)

**Thymic selection theory**

62. Sparse, random sampling is sufficient for central tolerance. bioRxiv [10.1101/2025.12.09.693230](https://doi.org/10.1101/2025.12.09.693230) (**preprint**)
63. Yates AJ. Theories and Quantification of Thymic Selection. *Front. Immunol.*, 2014. PMC3912788
64. The effects of thymic selection on the range of T cell cross-reactivity. PMC1857316
65. Influence of correlated antigen presentation on T-cell negative selection in the thymus. *J. R. Soc. Interface*, 2018. [10.1098/rsif.2018.0311](https://doi.org/10.1098/rsif.2018.0311)
