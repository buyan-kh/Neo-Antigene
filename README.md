# Neo Antigene

**Deterministic, self-grading science and statistics evaluation environments for AI labs, built from real bugs caught in a bioinformatics pipeline.**

---

## Overview

Frontier AI labs need evaluation environments for complex scientific and statistical reasoning that cannot be gamed. While language models can easily generate plausible-sounding text, evaluating whether a model has genuinely verified an analysis, caught subtle methodology flaws, or avoided data contamination requires rigorous, deterministic ground truth.

**Neo Antigene** provides deterministic, cheat-proof evaluation and reinforcement-learning environments where an AI agent is placed inside a realistic scientific workspace, supplied with messy raw data files (e.g., variant call files, gene annotations, plate reader outputs, and study manifests), and tasked with auditing claims or recomputing correct findings.

All grading is performed by pure Python verifiers against isolated answer keys—with **zero LLM-as-judge non-determinism** and **no network access**.

The primary product lives in [`proving-ground/`](proving-ground/):

| Environment | Focus Area | Failure Mode Caught |
| :--- | :--- | :--- |
| [`contam-sentinel`](proving-ground#1-contam-sentinel) | Data Integrity & Leakage | Hidden train/test contamination across aliased entity names and labeled holdout sets. |
| [`power-gate`](proving-ground#2-power-gate) | Statistical Inference | Conflating statistical significance with statistical power, and relying on stale design specs. |
| [`pvalue-police`](proving-ground#3-pvalue-police) | Multiple Testing | Reporting uncorrected nominal p-values across hundreds of features and failing to exclude controls. |
| [`coord-drift-detector`](proving-ground#4-coord-drift-detector) | Annotation Versioning | Silent 1-codon coordinate shifts between genome releases corrupting translated peptide sequences. |
| [`allele-normalizer`](proving-ground#5-allele-normalizer) | Entity Resolution | Dropping records during cross-table merges due to non-standardized biological nomenclature. |

See the [`proving-ground/README.md`](proving-ground/README.md) for full benchmark documentation, verification results, and honest model evaluation tables.

---

## Origins

The codebase in `src/neoantigene` was originally developed as an end-to-end computational biology pipeline for prioritizing tumor neoantigens (mutated protein fragments presented on tumor cells) from patient somatic mutations, RNA expression, and immune HLA types.

While building and validating that pipeline against published datasets and external tools, automated multi-agent audits uncovered pervasive methodological breakdowns across the field:
- Test benchmarks secretly contaminated with training data ([CEDAR / NetMHCpan](docs/HARNESS.md)).
- High-throughput discovery pipelines that omitted multiple-testing corrections or applied them after pre-filtering.
- Published studies claiming statistical power based on post-hoc estimates while conflicting with raw population data.
- Genome release updates that silently shifted translation reading frames.

Rather than treating these as one-off pipeline bugs, they became the foundation for our core product: **self-grading, reward-hack-resistant environments designed to evaluate whether frontier AI models can perform rigorous scientific auditing.** The ranking pipeline and its validation suite now serve as the real-world source material and testbed for these environments.

---

## The Underlying Engine

The pipeline code under `src/neoantigene` demonstrates the full data processing workflows that supply our evaluation environments:

```
somatic VCF + RNA expression + HLA type
        |
        v
  somatic filtering        clonal, expressed, not germline
        |
        v
  peptide enumeration      mutant + positionally matched wild-type
        |
        v
  presentation prediction  MHCflurry 2.x, pluggable
        |
        v
  ranking                  explicit logistic model over normalized features
        |
        v
  ranked shortlist -> evaluation harness -> audit
```

Implemented features include:
- Ingestion and somatic filtering with purity-corrected cancer cell fraction (CCF).
- Peptide enumeration including frameshift and neo-open reading frame (neo-ORF) tails.
- MHC presentation prediction via MHCflurry 2.x.
- A transparent, inspectable nine-feature scoring model with literature-justified weights.
- An automated evaluation harness ([`docs/HARNESS.md`](docs/HARNESS.md)) that checks benchmark statistical power before scoring.

---

## Installation

Python 3.12 managed with [uv](https://docs.astral.sh/uv/):

```bash
# Provision interpreter and dependencies
uv sync

# Real presentation predictions (optional)
uv sync --extra presentation
uv run mhcflurry-downloads fetch models_class1_presentation
```

Run CLI commands via `uv run`:

```bash
uv run neoantigene --help
```

---

## Quick Start

### 1. Run the Evaluation Environments (Proving Ground)

No external dependencies beyond Python 3.10+:

```bash
# Verify determinism across all 5 environments (42 reference submissions x 100 runs)
python3 proving-ground/verify_determinism.py --runs 100

# Run a specific environment grader
python3 proving-ground/envs/power_gate/grade.py \
  --submission proving-ground/envs/power_gate/reference/canonical/good.json \
  --fixture canonical
```

### 2. Run the Underlying Bioinformatic Pipeline

Smoke test on the bundled synthetic example:

```bash
uv run neoantigene rank data/examples/sample.yaml \
  --backend null --allow-null-backend \
  --out-dir results
```

Outputs:
- `<sample>.ranked.tsv` — Ranked peptide shortlist.
- `<sample>.features.json` — Scored peptide-HLA pairs with normalized feature values.
- `<sample>.run.json` — Reproducibility run manifest with input SHA-256 digests.

---

## Evaluation Harness

The repository includes a general-purpose evaluation harness for auditing ranking algorithms against experimental labels:

```bash
uv run neoantigene compare labels.tsv \
  -m ours=ranked.tsv -m baseline=baseline.tsv \
  --baseline baseline --k 20 -o report.md
```

Before reporting retrieval metrics (such as precision@k or concordance), the harness computes whether the assayed pool possesses sufficient statistical power to detect an effect above chance. If a dataset is provably underpowered at any outcome, the harness flags it directly. Protocol and accounting details are in [`docs/HARNESS.md`](docs/HARNESS.md).

---

## Development

```bash
# Run test suite
uv run pytest

# Formatting and linting
uv run ruff check .
uv run ruff format .
uv run mypy            # strict type-checking on src and tests
```

### Test Suite Structure

| File | What It Covers |
| :--- | :--- |
| `test_models.py` | Pydantic boundary validation and data immutability |
| `test_peptides.py`, `test_vcf.py`, `test_clonality.py` | Ingestion, somatic filtering, and peptide enumeration |
| `test_scoring.py`, `test_pipeline.py` | Scoring functions, gating, purity correction, determinism |
| `test_run.py`, `test_config.py` | Reproducibility manifests, config loading, and logging |
| `test_eval.py` | Precision@k and concordance against baselines |
| `test_cli.py` | Command line interface validation |

All test fixtures use synthetic or public benchmark data; no patient-identifiable sequencing data is stored in the repository.

---

## Repository Structure

```
Neo-Antigene/
  proving-ground/    Deterministic evaluation environments for AI labs
    contract/        Grader schemas and determinism harnesses
    envs/            The five environments (contam-sentinel, power-gate, etc.)
    verify_determinism.py Automated 100-run determinism verification suite
  src/neoantigene/   Bioinformatics engine (source material for the envs)
    variants/        Somatic variant ingestion and purity-corrected CCF
    peptides/        Proteome index and mutant peptide enumeration
    presentation/    MHC presentation prediction backends
    scoring/         Feature normalization and scoring model
    assays/          Ground-truth schema and evaluation metrics
  docs/              Audit records, benchmark accounting, and harness specs
    AUDIT.md         Source findings from multi-agent pipeline audits
    HARNESS.md       Evaluation harness methodology and statistical power rules
    BENCHMARK.md     Accounting for published validation cases
  tests/             Unit and integration tests for the core engine
```
