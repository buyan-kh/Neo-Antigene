"""Ranking.

The score is an explicit logistic model over normalized features. It is
deliberately not a black box: the shortlist a lab acts on has to be auditable,
and the same functional form is what `learning.train` refits once assay labels
come back.

Nothing here mutates its input; `score_all` returns new `ScoredCandidate`
objects.
"""

from __future__ import annotations

import math
from collections.abc import Sequence

from ..config import OutputConfig, ScoringWeights
from ..models import ScoredCandidate, VariantClass

SCORED_FEATURES: tuple[str, ...] = (
    "clonality",
    "expression",
    "presentation",
    "agretopicity",
    "tumor_selectivity",
    "mutation_exposure",
    "wt_dissimilarity",
    "self_dissimilarity",
    "hydrophobicity",
)


def weight_vector(weights: ScoringWeights) -> dict[str, float]:
    payload = weights.model_dump()
    return {name: float(payload[name]) for name in SCORED_FEATURES}


def logit(features: dict[str, float], weights: ScoringWeights) -> float:
    total = weights.bias
    for name, weight in weight_vector(weights).items():
        total += weight * features.get(name, 0.0)
    return total


def score_features(features: dict[str, float], weights: ScoringWeights) -> float:
    return 1.0 / (1.0 + math.exp(-logit(features, weights)))


def score_all(scored: Sequence[ScoredCandidate], weights: ScoringWeights) -> list[ScoredCandidate]:
    return [item.with_score(score_features(item.features, weights)) for item in scored]


def shortlist(
    scored: Sequence[ScoredCandidate],
    config: OutputConfig,
) -> list[ScoredCandidate]:
    """Rank, then cap per-variant redundancy.

    Without the cap a single strong variant fills the list with its own
    overlapping registers, which spends synthesis budget on one hypothesis.
    Ties break on candidate id so a run is reproducible.

    Frameshifts get a higher cap, because downstream of the shift the extra
    peptides are separate hypotheses rather than registers of one. See
    `OutputConfig.max_per_frameshift_variant`.
    """
    eligible = [s for s in scored if s.passed or config.include_failed]
    ordered = sorted(eligible, key=lambda s: (-s.score, s.candidate.id, s.allele))

    per_variant: dict[str, int] = {}
    selected: list[ScoredCandidate] = []
    for item in ordered:
        key = item.candidate.variant.key
        if per_variant.get(key, 0) >= _variant_cap(item, config):
            continue
        per_variant[key] = per_variant.get(key, 0) + 1
        selected.append(item)
        if len(selected) >= config.top_n:
            break
    return selected


def _variant_cap(item: ScoredCandidate, config: OutputConfig) -> int:
    if item.candidate.variant.variant_class is VariantClass.FRAMESHIFT:
        return config.max_per_frameshift_variant
    return config.max_per_variant


def contributions(features: dict[str, float], weights: ScoringWeights) -> dict[str, float]:
    """Per-feature logit contribution, for explaining a given rank."""
    return {
        name: weight * features.get(name, 0.0) for name, weight in weight_vector(weights).items()
    }
