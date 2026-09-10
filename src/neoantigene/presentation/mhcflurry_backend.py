"""MHCflurry 2.x presentation backend.

Install with the `presentation` extra, then fetch model weights once:

    uv sync --extra presentation
    uv run mhcflurry-downloads fetch models_class1_presentation
"""

from __future__ import annotations

import math
from collections.abc import Sequence
from typing import Any

from ..models import MAX_PEPTIDE_LENGTH, MIN_PEPTIDE_LENGTH, PresentationCall
from .base import PresentationBackend


class MHCflurryBackend(PresentationBackend):
    name = "mhcflurry"

    def __init__(self) -> None:
        try:
            from mhcflurry import Class1PresentationPredictor
        except ImportError as exc:  # pragma: no cover - depends on optional extra
            raise ImportError(
                "mhcflurry is not installed. Install with: uv sync --extra presentation "
                "and run: uv run mhcflurry-downloads fetch models_class1_presentation"
            ) from exc
        self._predictor = Class1PresentationPredictor.load()

    def supports_length(self, length: int) -> bool:
        return MIN_PEPTIDE_LENGTH <= length <= MAX_PEPTIDE_LENGTH

    def predict(
        self,
        peptides: Sequence[str],
        alleles: Sequence[str],
        n_flanks: Sequence[str] | None = None,
        c_flanks: Sequence[str] | None = None,
    ) -> list[PresentationCall]:
        if not peptides or not alleles:
            return []

        # One single-allele "sample" per allele gives the full peptide x allele
        # matrix; passing the genotype directly would collapse to a best allele.
        genotypes = {allele: [allele] for allele in alleles}
        frame = self._predictor.predict(
            peptides=list(peptides),
            alleles=genotypes,
            n_flanks=list(n_flanks) if n_flanks else None,
            c_flanks=list(c_flanks) if c_flanks else None,
            include_affinity_percentile=True,
            verbose=0,
        )
        records: list[dict[str, Any]] = frame.to_dict("records")
        return [
            PresentationCall(
                peptide=str(row["peptide"]),
                allele=str(row["sample_name"]),
                presentation_score=_finite(row.get("presentation_score")),
                affinity_nm=_finite(row.get("affinity")),
                affinity_percentile=_finite(row.get("affinity_percentile")),
                processing_score=_finite(row.get("processing_score")),
            )
            for row in records
        ]


def _finite(value: Any) -> float | None:
    if value is None:
        return None
    try:
        result = float(value)
    except (TypeError, ValueError):
        return None
    return result if math.isfinite(result) else None
