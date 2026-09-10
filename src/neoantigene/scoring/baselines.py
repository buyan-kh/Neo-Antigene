"""Comparator rankers.

The product claim is "our shortlist validates more often than what you use
today". That claim is only measurable against something, so the baselines are
first-class code with tests, not a notebook cell.

`BINDING_ONLY` is the one that matters commercially: it is the NetMHCpan-class
workflow of sorting candidates by predicted binding and taking the top N. If
Neo Antigene cannot beat it on assay hit rate, there is no product.
"""

from __future__ import annotations

import hashlib
from collections.abc import Callable, Sequence
from typing import Final

from ..models import ScoredCandidate

Ranker = Callable[[Sequence[ScoredCandidate]], dict[str, float]]


def _rank_by(
    scored: Sequence[ScoredCandidate], key: Callable[[ScoredCandidate], float]
) -> dict[str, float]:
    return {item.key: key(item) for item in scored}


def binding_only(scored: Sequence[ScoredCandidate]) -> dict[str, float]:
    """Sort by predicted binding strength alone.

    Uses affinity percentile when available (lower is stronger, so it is
    negated), falling back to raw nM. This is the standard-of-care baseline.
    """

    def strength(item: ScoredCandidate) -> float:
        call = item.mutant_call
        if call.affinity_percentile is not None:
            return -call.affinity_percentile
        if call.affinity_nm is not None:
            return -call.affinity_nm
        return float("-inf")

    return _rank_by(scored, strength)


def presentation_only(scored: Sequence[ScoredCandidate]) -> dict[str, float]:
    """Sort by the predictor's combined presentation score alone."""

    def strength(item: ScoredCandidate) -> float:
        score = item.mutant_call.presentation_score
        return score if score is not None else 0.0

    return _rank_by(scored, strength)


def expression_only(scored: Sequence[ScoredCandidate]) -> dict[str, float]:
    """Sort by source-transcript expression alone."""
    return _rank_by(scored, lambda item: item.features.get("expression", 0.0))


def arbitrary(scored: Sequence[ScoredCandidate], seed: str = "neoantigene") -> dict[str, float]:
    """A deterministic arbitrary order: the floor any real method must clear.

    Seeded by candidate key rather than by a PRNG so the ordering is stable
    across runs and machines, which keeps the eval suite from flaking.
    """

    def pseudo(item: ScoredCandidate) -> float:
        digest = hashlib.sha256(f"{seed}|{item.key}".encode()).digest()
        return int.from_bytes(digest[:8], "big") / float(1 << 64)

    return _rank_by(scored, pseudo)


def neoantigene(scored: Sequence[ScoredCandidate]) -> dict[str, float]:
    """The full model, for symmetry with the baselines."""
    return _rank_by(scored, lambda item: item.score)


BASELINES: Final[dict[str, Ranker]] = {
    "binding_only": binding_only,
    "presentation_only": presentation_only,
    "expression_only": expression_only,
    "arbitrary": arbitrary,
}

ALL_RANKERS: Final[dict[str, Ranker]] = {"neoantigene": neoantigene, **BASELINES}
