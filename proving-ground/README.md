# Proving Ground

**We sell the exam, not the student.**

AI labs can generate endless code and speculative text. What they cannot easily build is an **eval environment that cannot be cheated** on hard scientific and statistical reasoning. 

Proving Ground is a suite of deterministic, self-grading reinforcement learning and evaluation environments for frontier AI labs. Each environment places an agent inside a realistic scientific workflow—with messy raw files, real biological data, and ambiguous candidate claims—and evaluates whether the model can detect subtle methodology bugs, reconcile conflicting data sources, and compute the rigorous scientific answer.

Every environment here is derived from a real failure mode that was caught by hand in production bioinformatics and biostatistics pipelines.

---

## The Rules Every Environment Obeys

- **Deterministic code grading, never an LLM judge:** The grader is pure Python with zero network calls and no non-deterministic language model judges. Identical inputs produce byte-identical scores every time.
- **Answer keys remain out of reach:** The agent executes in an isolated sandbox containing only input data and a task specification (`task.md`). The ground truth answer key lives in a separate `key/` directory that is never accessible inside the sandbox.
- **Immune to reward hacking:** Graders verify recomputed distributions, recovered entity sets, numerical tolerances, and precise code locations. A model cannot guess its way to a pass with superficial explanations or boilerplate summaries.
- **Mirrored fixtures defeat lazy policies:** Several environments ship with mirrored pairs (e.g., an underpowered study vs. an adequately powered one, or full leakage vs. partial leakage). A static heuristic (such as always crying "underpowered" or "contaminated") receives a zero score on the counterpart fixture.

---

## How an Environment Runs

1. The model receives a workspace of input files and `task.md`. No network access is permitted.
2. The model investigates the files, performs its computations, and writes its structured findings to `/output/submission.json`.
3. A separate grading process reads the submission, compares it against the isolated `key/`, and returns a structured result:

```json
{
  "score": 0.8,
  "passed": false,
  "checks": [
    {"name": "recomputation", "passed": false, "detail": "Observed 594 leaked pairs; expected 514 (failed to exclude held_out_from_fit.csv)"}
  ],
  "reward_hack_flags": []
}
```

---

## The Five Environments

All five environments are fully implemented, parameterized, and covered by deterministic grading suites:

| Environment | Domain | The Failure Mode It Catches | What the Grader Verifies |
| :--- | :--- | :--- | :--- |
| **`contam-sentinel`** | Data Integrity & Leakage | Benchmark evaluation examples secretly present in training sets under aliased names. | Exact set of leaked pairs identified; exclusion of labeled holdout files; recomputed clean AUROC. |
| **`power-gate`** | Statistical Inference | Mistaking a nominally significant p-value for an adequately powered study, while relying on stale design specs. | True statistical power computed; sample size for 80% power derived; data vs. design discrepancy flagged. |
| **`pvalue-police`** | Multiple Testing | Reporting raw p-values without correction across hundreds of features, or correcting after pre-filtering. | Bug localized in script; control spike-ins excluded from correction; Benjamini-Hochberg recomputed on screen analytes. |
| **`coord-drift-detector`** | Gene Annotation Versioning | Silent 1-codon coordinate shifts between reference releases altering downstream translated peptides. | Detection of coordinate drift; affected transcripts identified; correct neo-ORF peptide sequences recomputed through the 3' UTR. |
| **`allele-normalizer`** | Entity Resolution | Records silently dropped during cross-table joins due to non-standardized biological nomenclature. | Full reconciliation across serotypes and heterodimers; zero silent drops; payload preservation. |

### 1. `contam-sentinel`
- **The problem:** A predictor claims state-of-the-art accuracy on an evaluation set. However, benchmark pairs were secretly present in the training corpus under different naming conventions.
- **The task:** The model must normalize identifiers, locate the exact overlapping pairs across panels, recognize that 80 pairs listed in `held_out_from_fit.csv` were legitimately held out, and recompute the clean AUROC (**Area Under the Receiver Operating Characteristic curve**—a standard classification metric ranging from 0.5 for random chance to 1.0 for perfect ranking) over the 972 genuine test pairs.
- **Key metrics:** Canonical true leak count is **514** (out of 594 apparent overlaps); clean AUROC is approximately **0.780**.

### 2. `power-gate`
- **The problem:** A clinical study reports an association with $p = 0.014$ and claims the trial was adequately powered. In reality, the analysis relied on post-hoc power at the observed effect, and the pre-registered study design document lists outdated population figures that conflict with delivered data.
- **The task:** The model must read the actual counts from `population.tsv` and `sample.tsv` rather than blindly trusting `study_design.json`, compute true **statistical power** (the probability that a statistical test will detect a genuine effect of a given magnitude), calculate the minimum sample size required to reach 80% power, and reject the headline claim.
- **Current key:**
  - True power: **~0.418** (0.418091).
  - Minimum sample size for 80% power: **235** (against an actual sample of 60).
  - Design vs. data discrepancy: Pre-registered group size is **120**, but actual delivered count is **115**.
- **Mirrored fixture:** An `adequately_powered` mirror ($N=5000, n=1200$, observed=52, power $>0.99$) tests whether the agent actually evaluates power or merely hedges reflexively.

### 3. `pvalue-police`
- **The problem:** A drug target screen tests hundreds of analytes and claims dozens of "significant hits" using uncorrected p-values, or applies correction only after cherry-picking candidates with $p < 0.05$. Furthermore, non-biological spike-in controls were mixed into the hypothesis pool.
- **The task:** The model must inspect `bin/select_candidates.py` to localize the methodology error, cross-reference `raw/sample_sheet.csv` to separate experimental screen analytes from control spike-ins, and apply the **Benjamini-Hochberg (BH) procedure** (a statistical method that controls the False Discovery Rate across multiple simultaneous hypothesis tests) strictly to the screen hypotheses.
- **Key metrics:** The reported list had 41 uncorrected/mis-filtered features; the true Benjamini-Hochberg screen-only set is exactly **4** features.

### 4. `coord-drift-detector`
- **The problem:** A variant-annotation pipeline translates neo-open reading frames (neo-ORFs) using cached coordinates generated against an older genome release, while the current annotation file reflects a newer release where the annotated start moved downstream by one codon.
- **The task:** The agent must compare the cached offsets against the official annotation, verify the sequence context in the FASTA file, detect which transcripts experience bona fide frame shifts versus neutral end-coordinate differences, and recompute the correct downstream peptide sequences extending through the 3' UTR.
- **Key metrics:** Drift detected with magnitude 1 codon downstream; 4 affected peptide sequences recomputed matching ground truth.

### 5. `allele-normalizer`
- **The problem:** Two clinical datasets recording values for immune alleles must be merged. Naive string matching discards over half the rows due to mixed nomenclature: legacy serological aliases (e.g., `A2` vs `HLA-A*02:01`), omitted gene prefixes, differing field resolutions, and single- versus double-prefixed heterodimers (`HLA-DQA1*01:02-HLA-DQB1*05:01` vs `HLA-DQA1*01:02-DQB1*05:01`).
- **The task:** Normalize all alleles into standard two-field **HLA** (Human Leukocyte Antigen—the cell-surface proteins that present peptide fragments to T cells) nomenclature, execute a complete join without dropping valid records, apply specified unit conversions, and ignore outdated entries in poisoned alias tables.
- **Key metrics:** All 30 shared alleles successfully joined; 0 silent drops; payload scaling verified to $10^{-6}$ precision.

---

## Determinism & Verification

A core design requirement is absolute determinism. Every grader has been subjected to stress testing across hundreds of repeated runs.

```bash
# Run the determinism verification harness across all environments
python3 verify_determinism.py --runs 100
```

**Verification result:**
- **42 reference submissions × 100 runs each: all identical, good/bad verdicts correct.**
- Zero score drift, zero flakiness, zero dependency on ambient environment state.

---

## Honest Flagship Results

We evaluate our environments against frontier LLMs (`grok-4.7-xhigh`, `composer-2.5-fast`, `gpt-5.6-sol-max`, `gemini-4-argon-high`). We report results transparently:

- **`contam-sentinel` is a verified, genuine flagship reasoning failure:** When tasked with auditing train/test overlap, flagship models (such as Cursor Grok 4.6 and Composer) failed by taking the naive union across panels (594 pairs) and failing to cross-check `held_out_from_fit.csv`, which explicitly flagged 80 pairs as excluded from model training. The correct expert answer is 514.
- **The other four environments are fair, deterministic, and currently solved by flagships:** Earlier iterations of `power-gate`, `pvalue-police`, `coord-drift-detector`, and `allele-normalizer` caused flagship failures due to unstated conventions or grader ambiguities (e.g., silently recounting files without instructing the model to reconcile data, or unstated string formatting conventions). Once each environment was hardened in **v10** to make the instructions scientifically fair and self-contained, current frontier models (`grok-4.7-xhigh` and `composer-2.5-fast`) successfully solved them.
- We do **not** claim "5/5 models fail." Our product is the cheat-proof, deterministic grading harness that accurately evaluates whether a model can perform rigorous scientific and statistical auditing.

---

## Running Environments Locally

Requires Python 3.10+. No external packages required for the core test harness:

```bash
# 1. Build fixtures for an environment
python3 envs/contam_sentinel/build.py
python3 envs/power_gate/build.py

# 2. Grade a submission
python3 envs/contam_sentinel/grade.py --submission path/to/submission.json --fixture canonical
python3 envs/power_gate/grade.py --submission path/to/submission.json --fixture canonical

# 3. Verify deterministic grading across all reference submissions
python3 verify_determinism.py --runs 100
```

---

## Repository Layout

```
proving-ground/
  contract/              Shared evaluation contract and result types
    core.py              Check schemas, scoring logic, and gate assertions
    determinism.py       Determinism testing framework
  envs/                  The five evaluation environments
    contam_sentinel/     Train/test contamination audit
    power_gate/          Statistical power and study design reconciliation
    pvalue_police/       Multiple testing correction and control exclusion
    coord_drift_detector/ Annotation release coordinate drift
    allele_normalizer/   HLA entity resolution and nomenclature join
  verify_determinism.py  Automated determinism verifier (42 submissions x 100 runs)
  run_models.py          Harness for executing frontier models against environments
  HARDENING_LOG.md       Changelog documenting environment hardening from v0 to v10
```
