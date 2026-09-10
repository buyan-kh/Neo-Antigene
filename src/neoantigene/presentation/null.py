"""Deterministic stand-in predictor.

DEVELOPMENT ONLY. Produces reproducible pseudo-scores so the pipeline, tests
and example data run without TensorFlow or downloaded model weights. It has no
biological validity and must never be used to pick peptides for synthesis; the
CLI refuses it unless `--allow-null-backend` is passed.
"""

from __future__ import annotations

import hashlib
import math
from collections.abc import Sequence

from ..models import PresentationCall
from .base import PresentationBackend


class NullBackend(PresentationBackend):
    name = "null"

    def predict(
        self,
        peptides: Sequence[str],
        alleles: Sequence[str],
        n_flanks: Sequence[str] | None = None,
        c_flanks: Sequence[str] | None = None,
    ) -> list[PresentationCall]:
        return [_call(peptide, allele) for peptide in peptides for allele in alleles]


def _call(peptide: str, allele: str) -> PresentationCall:
    unit = _unit_hash(f"{peptide}|{allele}")
    return PresentationCall(
        peptide=peptide,
        allele=allele,
        presentation_score=round(unit, 6),
        affinity_nm=round(10 ** (1.0 + 4.0 * (1.0 - unit)), 3),
        affinity_percentile=round(100.0 * (1.0 - unit) ** 2, 4),
        processing_score=round(_unit_hash(peptide), 6),
    )


def _unit_hash(key: str) -> float:
    digest = hashlib.sha256(key.encode()).digest()
    raw = int.from_bytes(digest[:8], "big") / float(1 << 64)
    return 1.0 / (1.0 + math.exp(-6.0 * (raw - 0.5)))
