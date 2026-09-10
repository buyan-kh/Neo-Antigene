"""Serialization of ranked shortlists."""

from __future__ import annotations

import json
import math
from collections.abc import Sequence
from pathlib import Path
from typing import Any

import polars as pl

from ..models import ScoredCandidate

COLUMNS: tuple[str, ...] = (
    "rank",
    "run_id",
    "sample_id",
    "gene",
    "hgvsp",
    "transcript",
    "variant",
    "variant_class",
    "mutant_peptide",
    "wildtype_peptide",
    "length",
    "mutation_position",
    "allele",
    "score",
    "presentation_score",
    "affinity_nm",
    "affinity_percentile",
    "wt_affinity_nm",
    "agretopicity",
    "ccf",
    "dna_vaf",
    "tpm",
    "rna_vaf",
    "normal_tpm",
    "gate_failures",
)

_SCHEMA: dict[str, pl.DataType] = {
    "rank": pl.Int64(),
    "run_id": pl.Utf8(),
    "sample_id": pl.Utf8(),
    "gene": pl.Utf8(),
    "hgvsp": pl.Utf8(),
    "transcript": pl.Utf8(),
    "variant": pl.Utf8(),
    "variant_class": pl.Utf8(),
    "mutant_peptide": pl.Utf8(),
    "wildtype_peptide": pl.Utf8(),
    "length": pl.Int64(),
    "mutation_position": pl.Int64(),
    "allele": pl.Utf8(),
    "score": pl.Float64(),
    "presentation_score": pl.Float64(),
    "affinity_nm": pl.Float64(),
    "affinity_percentile": pl.Float64(),
    "wt_affinity_nm": pl.Float64(),
    "agretopicity": pl.Float64(),
    "ccf": pl.Float64(),
    "dna_vaf": pl.Float64(),
    "tpm": pl.Float64(),
    "rna_vaf": pl.Float64(),
    "normal_tpm": pl.Float64(),
    "gate_failures": pl.Utf8(),
}


def _clean(value: float | None) -> float | None:
    """Drop NaN sentinels so they serialize as empty rather than `NaN`."""
    if value is None or math.isnan(value):
        return None
    return value


def to_records(
    scored: Sequence[ScoredCandidate], sample_id: str, run_id: str = ""
) -> list[dict[str, Any]]:
    records: list[dict[str, Any]] = []
    for index, item in enumerate(scored, start=1):
        variant = item.candidate.variant
        wildtype = item.wildtype_call
        records.append(
            {
                "rank": index,
                "run_id": run_id,
                "sample_id": sample_id,
                "gene": variant.gene,
                "hgvsp": variant.hgvsp_short,
                "transcript": variant.transcript,
                "variant": variant.key,
                "variant_class": variant.variant_class.value,
                "mutant_peptide": item.candidate.mutant_peptide,
                "wildtype_peptide": item.candidate.wildtype_peptide,
                "length": item.candidate.length,
                "mutation_position": item.candidate.mutation_position,
                "allele": item.allele,
                "score": round(item.score, 6),
                "presentation_score": item.mutant_call.presentation_score,
                "affinity_nm": item.mutant_call.affinity_nm,
                "affinity_percentile": item.mutant_call.affinity_percentile,
                "wt_affinity_nm": wildtype.affinity_nm if wildtype else None,
                "agretopicity": _clean(item.features.get("agretopicity_raw")),
                "ccf": variant.ccf,
                "dna_vaf": variant.dna_vaf,
                "tpm": _clean(item.features.get("tpm_raw")),
                "rna_vaf": variant.rna_vaf,
                "normal_tpm": _clean(item.features.get("normal_tpm_raw")),
                "gate_failures": ";".join(g.value for g in item.gate_failures),
            }
        )
    return records


def write_tsv(
    scored: Sequence[ScoredCandidate], sample_id: str, path: Path, run_id: str = ""
) -> pl.DataFrame:
    records = to_records(scored, sample_id, run_id)
    frame = pl.DataFrame(records, schema=_SCHEMA).select(COLUMNS)
    path.parent.mkdir(parents=True, exist_ok=True)
    frame.write_csv(path, separator="\t")
    return frame


def write_features_json(
    scored: Sequence[ScoredCandidate], sample_id: str, path: Path, run_id: str = ""
) -> Path:
    payload = [
        {
            "run_id": run_id,
            "sample_id": sample_id,
            "candidate_id": item.candidate.id,
            "peptide": item.candidate.mutant_peptide,
            "allele": item.allele,
            "score": item.score,
            "features": item.features,
            "gate_failures": [g.value for g in item.gate_failures],
        }
        for item in scored
    ]
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2))
    return path
