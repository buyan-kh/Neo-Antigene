"""Batch selection for the next round of synthesis and assay.

Assay capacity, not compute, is the binding constraint. A batch chosen purely
by score tells you almost nothing new; a batch chosen purely by uncertainty
wastes the round on peptides nobody would ever nominate. `select_batch` mixes
the two and then enforces diversity so the round covers distinct hypotheses.

There are two different uncertainties here and only one of them is worth
paying for.

`outcome_uncertainty` is how unsure the model is about *this peptide's*
result, and it peaks at a score of 0.5 by construction. It does not depend on
the weights at all, so it is unchanged by every label collected so far: a
peptide the model has confidently learned to score 0.5 looks exactly as
informative as one it has never had evidence about. Buying a label on the
first teaches nothing.

`weight_uncertainty` is how much the fitted weights would move if the data had
come out differently, measured as the spread of a candidate's score across an
ensemble of refits on resampled data. That is the quantity a label actually
reduces, so it is what an information-seeking batch should chase. It needs an
ensemble, which needs labels, so runs before the first refit necessarily fall
back to the outcome measure — documented here rather than silently, because
the fallback is a weaker criterion and not a substitute.
"""

from __future__ import annotations

import math
from collections.abc import Mapping, Sequence

from pydantic import BaseModel

from ..config import ScoringWeights
from ..models import ScoredCandidate
from ..scoring.rank import score_features

#: Largest achievable standard deviation of a probability, reached when an
#: ensemble splits evenly between 0 and 1. Dividing by it puts
#: `weight_uncertainty` on the same [0, 1] scale as the exploit term so
#: `exploitation_weight` keeps meaning the same thing either way.
_MAX_SCORE_SPREAD = 0.5


class BatchItem(BaseModel):
    candidate: ScoredCandidate
    exploit: float
    explore: float
    acquisition: float

    @property
    def reason(self) -> str:
        return "exploit" if self.exploit >= self.explore else "explore"


def outcome_uncertainty(score: float) -> float:
    """Peaks at 0.5, where this peptide's result is least predictable.

    Independent of how well-determined the weights are, so it does not measure
    what a label would teach. See the module docstring.
    """
    return 1.0 - 2.0 * abs(score - 0.5)


def weight_uncertainty(
    features: Mapping[str, float],
    ensemble: Sequence[ScoringWeights],
) -> float:
    """Spread of this candidate's score across an ensemble of refits, in [0, 1].

    High when the members disagree, which is exactly when the label would pin
    down a coefficient. Zero for a single-member ensemble, which carries no
    information about weight uncertainty.
    """
    if len(ensemble) < 2:
        return 0.0
    scores = [score_features(dict(features), weights) for weights in ensemble]
    mean = sum(scores) / len(scores)
    variance = sum((score - mean) ** 2 for score in scores) / (len(scores) - 1)
    return min(1.0, math.sqrt(variance) / _MAX_SCORE_SPREAD)


def select_batch(
    scored: Sequence[ScoredCandidate],
    batch_size: int,
    exploitation_weight: float = 0.6,
    max_per_variant: int = 2,
    max_per_allele: int = 0,
    ensemble: Sequence[ScoringWeights] | None = None,
) -> list[BatchItem]:
    """Rank by acquisition value, then cap redundancy per variant and allele.

    `exploitation_weight` = 1.0 reproduces the plain shortlist; 0.0 is pure
    uncertainty sampling. The default leans toward delivering hits this round
    while still buying information for the next model fit. `max_per_allele`
    of 0 means unlimited.

    Pass `ensemble` — from `learning.train.bootstrap_weights` — to score the
    explore term by how much the weights are still undetermined. Without it
    the term falls back to outcome uncertainty, which is all that can be
    computed before any labels exist.
    """
    if not 0.0 <= exploitation_weight <= 1.0:
        raise ValueError(f"exploitation_weight must be in [0, 1], got {exploitation_weight}")
    if batch_size <= 0:
        raise ValueError(f"batch_size must be positive, got {batch_size}")

    def explore_value(entry: ScoredCandidate) -> float:
        if ensemble is not None and len(ensemble) >= 2:
            return weight_uncertainty(entry.features, ensemble)
        return outcome_uncertainty(entry.score)

    items = [
        BatchItem(
            candidate=entry,
            exploit=entry.score,
            explore=explore_value(entry),
            acquisition=(
                exploitation_weight * entry.score
                + (1.0 - exploitation_weight) * explore_value(entry)
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
