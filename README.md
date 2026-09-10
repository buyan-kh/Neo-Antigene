# Neo Antigene

Rank tumor vaccine targets so more of them actually trigger real T cells, and
prove it with lab assays.

Neo Antigene takes a patient's somatic variants, RNA expression and HLA type,
and returns a ranked peptide shortlist optimized for one number: the fraction
of top-ranked peptides that come back positive in a T-cell assay. That number,
measured against NetMHCpan-class baselines on the same assayed pool, is the
KPI the whole system is built around.

Not a consumer product. Customers are labs and biotechs using the ranking
tool, and pharma licensing the platform once the hit rate beats standard
tools.

## What it does today

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
  ranked shortlist -> synthesis -> T-cell assay -> labels -> refit
```

Implemented: ingestion, somatic filtering with purity-corrected CCF, peptide
enumeration, presentation prediction, ranking, assay request sheets, the KPI
metric, weight refitting, and active-learning batch selection.

Not implemented: the FASTQ/BAM front end (`neoantigene.fastq`), the AlphaFold 3
pMHC secondary filter (`neoantigene.structure`), and frameshift / neo-ORF
peptides. Each raises a clear error rather than silently degrading.

## Install

Python 3.12, managed with [uv](https://docs.astral.sh/uv/). `uv sync` provisions
the interpreter, installs from the committed `uv.lock`, and needs no
pre-existing virtualenv.

```bash
uv sync

# real presentation predictions
uv sync --extra presentation
uv run mhcflurry-downloads fetch models_class1_presentation
```

Run anything through `uv run`, which keeps the environment in sync with the
lockfile:

```bash
uv run neoantigene --help
```

## Quick start

Smoke test on the bundled example, using the development stand-in predictor:

```bash
uv run neoantigene rank data/examples/sample.yaml \
  --backend null --allow-null-backend \
  --out-dir results
```

A real run needs MHCflurry and a full Ensembl proteome:

```bash
uv run neoantigene rank my-sample.yaml --out-dir results
```

Outputs land in `results/`:

- `<sample>.ranked.tsv` — the shortlist a lab acts on
- `<sample>.features.json` — every scored peptide-HLA pair with its features,
  which is what `refit` and `evaluate` consume later
- `<sample>.run.json` — the run manifest

## Reproducibility

Every run gets a sortable run id (`20260909T184233Z-1a2b3c4d`) that is stamped
on every log line and written into both output files. The run manifest records
the package version, interpreter, platform, backend, fully resolved config, a
digest of that config, and a SHA-256 of every input file.

That means a shortlist a lab spent money on can always be traced back to the
exact inputs and settings that produced it, and a partner auditing the platform
can verify it independently. Pass `--json-logs` for structured output.

## Inputs

A sample manifest points at everything for one patient:

```yaml
sample_id: PDAC-001
cancer_type: pancreatic_adenocarcinoma
tumor_purity: 0.6

hla: [HLA-A*11:01, HLA-A*02:01, HLA-B*07:02, HLA-C*08:02]

processed:
  somatic_vcf: somatic.vep.vcf.gz
  expression_tsv: rna_tpm.tsv
  normal_expression_tsv: gtex_median_tpm.tsv   # optional
  proteome_fasta: Homo_sapiens.GRCh38.pep.all.fa.gz
  tumor_sample_name: TUMOR
```

**Somatic VCF** must be VEP-annotated — Neo Antigene reads protein
consequences rather than re-deriving them:

```bash
vep --input_file somatic.vcf --format vcf --vcf --symbol --terms SO \
    --canonical --biotype --transcript_version --cache --offline \
    --output_file somatic.vep.vcf
```

If you annotate elsewhere, use `variant_tsv:` instead and supply
`chrom, pos, ref, alt, transcript, protein_position, amino_acids, consequence`
(see `data/examples/variants.tsv`).

**Proteome** must be the Ensembl release matching your VEP cache. A mismatch
surfaces as a `ReferenceMismatch` on the first variant, which is deliberate —
silently ranking peptides built from the wrong reference is worse than
failing.

## How ranking works

The score is `sigmoid(bias + Σ wᵢxᵢ)` over seven features, each normalized to
`[0, 1]` and oriented so higher is better:

| feature | prior weight | what it captures |
| --- | --- | --- |
| `presentation` | 2.5 | MHCflurry presentation score for this peptide-allele pair |
| `clonality` | 1.5 | purity- and ploidy-corrected cancer cell fraction |
| `expression` | 1.2 | tumor TPM of the source transcript |
| `agretopicity` | 0.5 | WT/mutant affinity ratio — did the mutation create the binding event |
| `tumor_selectivity` | 0.4 | inverse of healthy-tissue expression |
| `foreignness` | 0.25 | BLOSUM62 dissimilarity to the wild-type peptide at TCR contacts |
| `hydrophobicity` | 0.15 | Kyte-Doolittle over TCR-facing residues |

Weights are scaled by strength of evidence. The top three are close to
necessary conditions for a response; the bottom four are proxies with thinner
support. That gap is deliberate — with comparable weights, the speculative
terms dominate the top of the shortlist, which is precisely the part a lab
synthesizes. `tests/test_eval.py` fails if that regresses.

Hard gates run before scoring and drop candidates outright: weak binders,
subclonal variants, germline-common variants, unexpressed transcripts, and any
peptide that exactly matches a sequence in the reference proteome.

Three deliberate choices:

- **There is no driver-gene term.** Driver status is a proxy for publication
  history, not for T-cell response. Ranking KRAS G12D highly because it is
  famous, rather than because it is clonal, expressed, presented and
  untolerized, is exactly the failure mode this tool exists to avoid.
- **The model is linear and inspectable.** `neoantigene explain` decomposes any
  candidate's score into per-feature contributions. A lab spending its own
  money on synthesis is entitled to know why a peptide was nominated.
- **Structure is a filter, not a score.** The pMHC structure step can demote a
  candidate, never promote one, so a folding artifact cannot override the
  assay-fit sequence model.

The shipped weights in `config/default.yaml` are **priors, not fitted values**.
They encode which features should matter and roughly how much. Replace them
with `neoantigene refit` output as soon as you have labels.

## The assay loop

```bash
# 1. shortlist
uv run neoantigene rank sample.yaml -o results

# 2. request sheet for synthesis and assay
uv run neoantigene request-assays results/PDAC-001.ranked.tsv -o assays/request.tsv

# 3. lab fills in `call` and `effect_size`, then measure the KPI
uv run neoantigene evaluate assays/completed.tsv \
  --ranked results/PDAC-001.ranked.tsv \
  --baseline results/netmhcpan.tsv \
  --k 20

# 4. refit weights once ~40 labelled peptides exist
uv run neoantigene refit assays/completed.tsv \
  -f results/PDAC-001.features.json \
  --out config/fitted.yaml
```

`evaluate` compares methods only over the intersection of peptides that were
actually assayed, so a method cannot win by having nominated an easier set.
Peptides that were ranked but never tested are excluded rather than counted as
negatives.

For the next round, `--select active` trades some expected hits for
information gain, which is what makes the labels worth more than the peptides:

```bash
uv run neoantigene rank sample.yaml --select active --exploitation 0.6 --top-n 24
```

## First wedge

Pancreatic adenocarcinoma, peptides, one standardized T-cell assay.

`neoantigene.wedge.pdac` catalogs recurrent PDAC variants with documented HLA
restrictions — KRAS G12D on `HLA-C*08:02` (`GADGVGKSA`, `GADGVGKSAL`) and
`HLA-A*11:01` (`VVGADGVGK`), G12V on `HLA-A*11:01` (`VVGAVGVGK`). Shared
antigens make the assay loop cheap: the same synthesized peptide can be tested
against many donors, so ground truth accumulates faster than it would with
private neoantigens.

```bash
uv run neoantigene catalog --hla "A*11:01,C*08:02"
```

The catalog is a stocking and screening aid. It does not feed the score.

## Development

```bash
uv run pytest
uv run ruff check .
uv run ruff format .
uv run mypy            # strict, on src and tests
```

CI runs all four on every push and pull request, plus a check that no
sequencing data has been committed outside `data/examples/`.

Conventions:

- **Library first, CLI second.** `neoantigene.cli` parses arguments and prints;
  it contains no ranking logic.
- **Stages are pure functions.** `ingest -> filter -> present -> rank -> export`,
  each taking and returning models. Nothing mutates its input.
- **Pydantic at every boundary.** Models are frozen with their invariants
  enforced at construction, so a peptide with a bad residue or an out-of-range
  percentile fails immediately rather than three stages later.
- **No web UI, no workflow engine, no services.** Not until a lab is actually
  uploading.

### Test layers

| file | what it covers |
| --- | --- |
| `test_models.py` | boundary validation — malformed data cannot travel |
| `test_peptides.py`, `test_vcf.py`, `test_clonality.py` | ingestion and enumeration, including published KRAS epitopes |
| `test_scoring.py`, `test_pipeline.py` | ranking, gating, purity, determinism |
| `test_run.py`, `test_config.py` | reproducibility and config loading |
| `test_eval.py` | precision@k against baselines (`-m eval`) |
| `test_cli.py` | command wiring and guard rails |

Fixture patients are generated, never real — see `tests/fixtures/synthetic.py`.
`NullBackend` is a deterministic hash-based stand-in that lets the pipeline and
tests run without TensorFlow or downloaded weights. It has no biological
validity, and the CLI refuses it unless you pass `--allow-null-backend`.

### Benchmarking

`neoantigene benchmark` re-ranks a sample with every built-in baseline
(`binding_only`, `presentation_only`, `expression_only`, `arbitrary`) and
reports hit rate for each over the same assayed pool:

```bash
uv run neoantigene benchmark sample.yaml assays/completed.tsv --k 20
```

`binding_only` is the one that matters commercially — it is the
NetMHCpan-class workflow of sorting by predicted binding and taking the top N.
If Neo Antigene cannot beat it on real assay data, there is no product.

## Layout

```
src/neoantigene/
  models.py      frozen pydantic models exchanged between stages
  config.py      pipeline configuration (YAML or TOML)
  run.py         run ids, logging, reproducibility manifest
  pipeline.py    stage composition
  io/            manifests, VEP VCF, expression, output
  variants/      somatic gating, purity-corrected CCF
  peptides/      proteome index, mutant peptide enumeration
  presentation/  MHCflurry backend, registry, dev stand-in
  scoring/       feature assembly, immunogenicity terms, ranking
  assays/        ground-truth schema, validation-rate KPI
  learning/      training set assembly, weight refit, active batch selection
  wedge/         PDAC shared-antigen catalog
  structure/     pMHC secondary filter interface (not implemented)
  fastq/         raw sequencing front end (not implemented)
```

## Data handling

`.gitignore` blocks VCF, BAM and FASTQ files and `data/samples/` by default,
and CI fails if any are committed outside `data/examples/`. Patient-derived
data must not be committed. All test fixtures are synthetic.
