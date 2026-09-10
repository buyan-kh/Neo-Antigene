"""Deterministic synthetic patients and assay outcomes.

Two things are generated here:

  - `build_cohort` writes a proteome, per-patient variant/expression tables
    and manifests to disk. This exercises the real ingestion path rather than
    hand-building model objects.
  - `simulate_assay_results` invents ground truth under an explicit
    generative model, so ranking quality can be measured before real assay
    data exists.

The generative model is an assumption, not a finding. It says a peptide
elicits a response when it is presented AND clonal AND expressed. Evaluation
tests built on it therefore check that the ranker recovers a stated
multi-factor truth — not that the biology is correct.
"""

from __future__ import annotations

import math
import random
from collections.abc import Sequence
from dataclasses import dataclass
from pathlib import Path

import yaml

from neoantigene.assays.schema import AssayCall, AssayResult, AssayType
from neoantigene.models import AMINO_ACIDS, ScoredCandidate

RESIDUES = "".join(sorted(AMINO_ACIDS))

DEFAULT_HLA = ("HLA-A*02:01", "HLA-A*11:01", "HLA-B*07:02", "HLA-C*08:02")

#: Coefficients of the synthetic truth. Presentation matters most, but a
#: subclonal or unexpressed peptide is not rescued by strong binding.
TRUTH_PRESENTATION = 3.4
TRUTH_CLONALITY = 2.6
TRUTH_EXPRESSION = 2.0
#: Tuned so roughly one candidate in eight validates. Real neoantigen hit
#: rates are low, and a benchmark run at a 50% base rate has no headroom to
#: distinguish between methods.
TRUTH_BIAS = -7.0

#: Probability a truly-positive peptide is still called negative, and vice
#: versa. Real assays are noisy; a benchmark that assumes they are not will
#: overstate every method.
FALSE_NEGATIVE_RATE = 0.10
FALSE_POSITIVE_RATE = 0.03


@dataclass(frozen=True)
class SyntheticPatient:
    sample_id: str
    manifest_path: Path
    variant_count: int


@dataclass(frozen=True)
class SyntheticCohort:
    root: Path
    proteome_path: Path
    patients: tuple[SyntheticPatient, ...]

    def patient(self, index: int = 0) -> SyntheticPatient:
        return self.patients[index]


def _protein(rng: random.Random, length: int) -> str:
    return "".join(rng.choice(RESIDUES) for _ in range(length))


def _write_proteome(path: Path, genes: Sequence[tuple[str, str]]) -> None:
    lines: list[str] = []
    for index, (symbol, sequence) in enumerate(genes, start=1):
        lines.append(
            f">ENSP9{index:010d}.1 pep synthetic:SYNTH:1:1:1:1 "
            f"gene:ENSG9{index:010d}.1 transcript:ENST9{index:010d}.1 "
            f"gene_biotype:protein_coding transcript_biotype:protein_coding "
            f"gene_symbol:{symbol}"
        )
        for offset in range(0, len(sequence), 60):
            lines.append(sequence[offset : offset + 60])
    path.write_text("\n".join(lines) + "\n")


def build_cohort(
    root: Path,
    patients: int = 2,
    genes: int = 12,
    protein_length: int = 240,
    variants_per_patient: int = 14,
    seed: int = 20260909,
) -> SyntheticCohort:
    """Write a complete synthetic cohort to `root` and return its layout."""
    rng = random.Random(seed)
    root.mkdir(parents=True, exist_ok=True)

    symbols = [f"SYN{index:03d}" for index in range(1, genes + 1)]
    sequences = [_protein(rng, protein_length) for _ in symbols]
    proteome_path = root / "proteome.fa"
    _write_proteome(proteome_path, list(zip(symbols, sequences, strict=True)))

    # Expression is shared across the cohort: a third of genes are quiet, which
    # is what lets the expression term separate candidates at all.
    tumor_tpm = {
        symbol: (rng.uniform(0.0, 0.8) if index % 3 == 0 else rng.uniform(5.0, 400.0))
        for index, symbol in enumerate(symbols)
    }
    normal_tpm = {symbol: rng.uniform(0.0, 60.0) for symbol in symbols}

    built: list[SyntheticPatient] = []
    for patient_index in range(patients):
        sample_id = f"SYN-PDAC-{patient_index + 1:03d}"
        directory = root / sample_id
        directory.mkdir(parents=True, exist_ok=True)

        _write_expression(directory / "expression.tsv", symbols, tumor_tpm, len(symbols))
        _write_expression(directory / "normal_expression.tsv", symbols, normal_tpm, len(symbols))
        count = _write_variants(
            directory / "variants.tsv", rng, symbols, sequences, variants_per_patient
        )

        manifest_path = directory / "sample.yaml"
        manifest_path.write_text(
            yaml.safe_dump(
                {
                    "sample_id": sample_id,
                    "cancer_type": "pancreatic_adenocarcinoma",
                    "tumor_purity": 0.6,
                    "hla": list(DEFAULT_HLA),
                    "processed": {
                        "variant_tsv": "variants.tsv",
                        "expression_tsv": "expression.tsv",
                        "normal_expression_tsv": "normal_expression.tsv",
                        "proteome_fasta": "../proteome.fa",
                    },
                },
                sort_keys=False,
            )
        )
        built.append(SyntheticPatient(sample_id, manifest_path, count))

    return SyntheticCohort(root, proteome_path, tuple(built))


def _write_expression(
    path: Path, symbols: Sequence[str], values: dict[str, float], count: int
) -> None:
    rows = ["transcript_id\tgene\ttpm"]
    for index, symbol in enumerate(symbols[:count], start=1):
        rows.append(f"ENST9{index:010d}\t{symbol}\t{values[symbol]:.3f}")
    path.write_text("\n".join(rows) + "\n")


def _write_variants(
    path: Path,
    rng: random.Random,
    symbols: Sequence[str],
    sequences: Sequence[str],
    count: int,
) -> int:
    header = (
        "chrom\tpos\tref\talt\tgene\ttranscript\tconsequence\tprotein_position\t"
        "amino_acids\tdna_vaf\trna_vaf\ttumor_depth\tcopy_number\tpopulation_af\tfilter"
    )
    rows = [header]
    used: set[tuple[int, int]] = set()
    written = 0
    attempts = 0

    while written < count and attempts < count * 20:
        attempts += 1
        gene_index = rng.randrange(len(symbols))
        sequence = sequences[gene_index]
        # Keep away from the termini so every requested peptide length fits.
        position = rng.randrange(20, len(sequence) - 20)
        if (gene_index, position) in used:
            continue
        used.add((gene_index, position))

        reference = sequence[position]
        alternate = rng.choice([r for r in RESIDUES if r != reference])
        # A quarter of variants are markedly subclonal, so clonality has
        # something to discriminate on.
        vaf = rng.uniform(0.05, 0.14) if written % 4 == 0 else rng.uniform(0.30, 0.48)

        rows.append(
            "\t".join(
                [
                    str(gene_index + 1),
                    str(100_000 + position * 3),
                    "C",
                    "T",
                    symbols[gene_index],
                    f"ENST9{gene_index + 1:010d}",
                    "missense_variant",
                    str(position + 1),
                    f"{reference}/{alternate}",
                    f"{vaf:.4f}",
                    f"{vaf * rng.uniform(0.7, 1.2):.4f}",
                    str(rng.randrange(60, 200)),
                    "2",
                    ".",
                    "PASS",
                ]
            )
        )
        written += 1

    path.write_text("\n".join(rows) + "\n")
    return written


def truth_probability(features: dict[str, float]) -> float:
    """Probability of a positive assay under the synthetic generative model."""
    logit = (
        TRUTH_BIAS
        + TRUTH_PRESENTATION * features.get("presentation", 0.0)
        + TRUTH_CLONALITY * features.get("clonality", 0.0)
        + TRUTH_EXPRESSION * features.get("expression", 0.0)
    )
    return 1.0 / (1.0 + math.exp(-logit))


def simulate_assay_results(
    scored: Sequence[ScoredCandidate],
    sample_id: str,
    seed: int = 7,
    assay: AssayType = AssayType.IFNG_ELISPOT,
    indeterminate_every: int = 0,
) -> list[AssayResult]:
    """Invent assay calls for already-scored candidates.

    Labels depend only on the *features*, never on the model's score, so a
    method cannot be rewarded for agreeing with itself.
    """
    rng = random.Random(seed)
    results: list[AssayResult] = []
    for index, item in enumerate(scored):
        truth = rng.random() < truth_probability(item.features)
        flipped = rng.random() < (FALSE_NEGATIVE_RATE if truth else FALSE_POSITIVE_RATE)
        positive = truth != flipped

        if indeterminate_every and index % indeterminate_every == 0:
            call = AssayCall.INDETERMINATE
        else:
            call = AssayCall.POSITIVE if positive else AssayCall.NEGATIVE

        results.append(
            AssayResult(
                sample_id=sample_id,
                peptide=item.candidate.mutant_peptide,
                allele=item.allele,
                assay=assay,
                call=call,
                effect_size=round(rng.uniform(40.0, 500.0), 1) if positive else 0.0,
                replicate_count=3,
            )
        )
    return results
