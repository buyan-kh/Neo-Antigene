"""Batch selection for the next round of synthesis and assay.

Assay capacity, not compute, is the binding constraint. A batch chosen purely
by score tells you almost nothing new; a batch chosen purely by uncertainty
wastes the round on peptides nobody would ever nominate. `select_batch` mixes
the two and then enforces diversity so the round covers distinct hypotheses.
"""

from __future__ import annotations

from collections.abc import Sequence

from pydantic import BaseModel

from ..models import ScoredCandidate


class BatchItem(BaseModel):
    candidate: ScoredCandidate
    exploit: float
    explore: float
    acquisition: float

    @property
    def reason(self) -> str:
        return "exploit" if self.exploit >= self.explore else "explore"


def uncertainty(score: float) -> float:
    """Peaks at 0.5, where the model is least sure and a label is worth most."""
    return 1.0 - 2.0 * abs(score - 0.5)


def select_batch(
    scored: Sequence[ScoredCandidate],
    batch_size: int,
    exploitation_weight: float = 0.6,
    max_per_variant: int = 2,
    max_per_allele: int = 0,
) -> list[BatchItem]:
    """Rank by acquisition value, then cap redundancy per variant and allele.

    `exploitation_weight` = 1.0 reproduces the plain shortlist; 0.0 is pure
    uncertainty sampling. The default leans toward delivering hits this round
    while still buying information for the next model fit. `max_per_allele`
    of 0 means unlimited.
    """
    if not 0.0 <= exploitation_weight <= 1.0:
        raise ValueError(f"exploitation_weight must be in [0, 1], got {exploitation_weight}")
    if batch_size <= 0:
        raise ValueError(f"batch_size must be positive, got {batch_size}")

    items = [
        BatchItem(
            candidate=entry,
            exploit=entry.score,
            explore=uncertainty(entry.score),
            acquisition=(
                exploitation_weight * entry.score
                + (1.0 - exploitation_weight) * uncertainty(entry.score)
            ),
        )
        for entry in scored
        if entry.passed
    ]
    items.sort(key=lambda i: (-i.acquisition, i.candidate.candidate.id, i.candidate.allele))

    per_variant: dict[str, int] = {}
    per_allele: dict[str, int] = {}
    selected: list[BatchItem] = []
    for item in items:
        variant_key = item.candidate.candidate.variant.key
        allele = item.candidate.allele
        if max_per_variant and per_variant.get(variant_key, 0) >= max_per_variant:
            continue
        if max_per_allele and per_allele.get(allele, 0) >= max_per_allele:
            continue
        per_variant[variant_key] = per_variant.get(variant_key, 0) + 1
        per_allele[allele] = per_allele.get(allele, 0) + 1
        selected.append(item)
        if len(selected) >= batch_size:
            break
    return selected
