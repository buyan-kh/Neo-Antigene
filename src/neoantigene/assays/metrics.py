"""The KPI: validation rate of top-ranked neoantigens.

`validation_rate_at_k` is the number the whole company is judged on. It is
only meaningful when the same assayed peptide set is re-ranked by each method,
so `compare_rankings` restricts every method to the intersection of what they
all nominated.
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence

from pydantic import BaseModel

from .schema import AssayResult

BASELINE_METHOD_EXCLUDED = "neoantigene"


class ValidationSummary(BaseModel):
    method: str
    k: int
    assayed: int
    positives: int
    hits_in_top_k: int
    validation_rate: float
    baseline_rate: float | None = None

    @property
    def lift(self) -> float | None:
        if not self.baseline_rate:
            return None
        return self.validation_rate / self.baseline_rate

    def describe(self) -> str:
        line = (
            f"{self.method:<20} k={self.k:<4} "
            f"validated {self.hits_in_top_k}/{self.k} = {self.validation_rate:.1%}  "
            f"(assayed pool {self.assayed}, positives {self.positives})"
        )
        lift = self.lift
        if lift is not None and self.method == BASELINE_METHOD_EXCLUDED:
            line += f"  lift x{lift:.2f}"
        return line


def labelled(results: Sequence[AssayResult]) -> list[AssayResult]:
    return [r for r in results if r.label is not None]


def validation_rate_at_k(
    results: Sequence[AssayResult],
    ranking: Mapping[str, float],
    k: int,
    method: str = BASELINE_METHOD_EXCLUDED,
) -> ValidationSummary:
    """Fraction of the top-k ranked peptides that scored positive in assay.

    `ranking` maps `AssayResult.key` to a score, higher is better. Peptides
    that were ranked but never assayed are excluded rather than counted as
    negatives: not testing a peptide is not evidence against it.
    """
    pool = [r for r in labelled(results) if r.key in ranking]
    if not pool:
        return ValidationSummary(
            method=method, k=0, assayed=0, positives=0, hits_in_top_k=0, validation_rate=0.0
        )

    # Tie-break on key so equal scores produce a stable, reproducible order.
    ordered = sorted(pool, key=lambda r: (-ranking[r.key], r.key))
    top = ordered[:k]
    hits = sum(1 for r in top if r.label == 1)
    return ValidationSummary(
        method=method,
        k=len(top),
        assayed=len(pool),
        positives=sum(1 for r in pool if r.label == 1),
        hits_in_top_k=hits,
        validation_rate=hits / len(top) if top else 0.0,
    )


def compare_rankings(
    results: Sequence[AssayResult],
    rankings: Mapping[str, Mapping[str, float]],
    k: int,
    baseline: str | None = None,
) -> list[ValidationSummary]:
    """Head-to-head validation rate for each method over the same assayed pool.

    Restricting to the intersection is the point: a method cannot win by
    having nominated a different, easier set of peptides.
    """
    if not rankings:
        return []
    shared = set.intersection(*(set(r) for r in rankings.values()))
    pool = [r for r in labelled(results) if r.key in shared]

    summaries = [
        validation_rate_at_k(pool, ranking, k, method=name) for name, ranking in rankings.items()
    ]

    reference = _reference_rate(summaries, baseline)
    if reference is None:
        return summaries
    return [s.model_copy(update={"baseline_rate": reference}) for s in summaries]


def _reference_rate(summaries: Sequence[ValidationSummary], baseline: str | None) -> float | None:
    if baseline is not None:
        match = next((s for s in summaries if s.method == baseline), None)
        return None if match is None else match.validation_rate
    fallback = next((s for s in summaries if s.method != BASELINE_METHOD_EXCLUDED), None)
    return None if fallback is None else fallback.validation_rate


class RetrievalSummary(BaseModel):
    """The four numbers a published benchmark is reported on.

    Precision and recall are both here because they answer different questions
    and diverge sharply when positives are scarce. Precision@10 is what a lab
    feels: of ten peptides synthesized, how many worked. Recall@10 is what a
    paper asks: of the epitopes that exist, how many did the method surface.
    With three positives in a pool, precision@10 cannot exceed 0.3 however
    perfect the ranking, so quoting it alone understates a good method.
    """

    method: str
    assayed: int
    positives: int
    recall_at_10: float
    recall_at_20: float
    precision_at_10: float
    hits_at_10: int
    hits_at_20: int
    best_positive_rank: int | None

    def describe(self) -> str:
        best = "none ranked" if self.best_positive_rank is None else f"#{self.best_positive_rank}"
        return (
            f"{self.method:<20} "
            f"recall@10 {self.recall_at_10:>6.1%} ({self.hits_at_10}/{self.positives})  "
            f"recall@20 {self.recall_at_20:>6.1%} ({self.hits_at_20}/{self.positives})  "
            f"precision@10 {self.precision_at_10:>6.1%}  "
            f"best {best}"
        )


def positive_ranks(
    results: Sequence[AssayResult],
    ranking: Mapping[str, float],
) -> list[int]:
    """1-based ranks of every assay-positive peptide within the assayed pool.

    Scoped to peptides that were both ranked and assayed, which is the same
    rule `validation_rate_at_k` uses.
    """
    pool = [r for r in labelled(results) if r.key in ranking]
    ordered = sorted(pool, key=lambda r: (-ranking[r.key], r.key))
    return [index for index, result in enumerate(ordered, start=1) if result.label == 1]


def retrieval_at_k(
    results: Sequence[AssayResult],
    ranking: Mapping[str, float],
    method: str = BASELINE_METHOD_EXCLUDED,
) -> RetrievalSummary:
    """Recall@10, recall@20, precision@10 and the best positive's rank."""
    pool = [r for r in labelled(results) if r.key in ranking]
    positives = sum(1 for r in pool if r.label == 1)
    ranks = positive_ranks(results, ranking)

    hits_10 = sum(1 for r in ranks if r <= 10)
    hits_20 = sum(1 for r in ranks if r <= 20)
    depth_10 = min(10, len(pool))

    return RetrievalSummary(
        method=method,
        assayed=len(pool),
        positives=positives,
        recall_at_10=hits_10 / positives if positives else 0.0,
        recall_at_20=hits_20 / positives if positives else 0.0,
        precision_at_10=hits_10 / depth_10 if depth_10 else 0.0,
        hits_at_10=hits_10,
        hits_at_20=hits_20,
        best_positive_rank=min(ranks) if ranks else None,
    )


def compare_retrieval(
    results: Sequence[AssayResult],
    rankings: Mapping[str, Mapping[str, float]],
) -> list[RetrievalSummary]:
    """`retrieval_at_k` for every method over one shared assayed pool."""
    if not rankings:
        return []
    shared = set.intersection(*(set(r) for r in rankings.values()))
    pool = [r for r in labelled(results) if r.key in shared]
    return [retrieval_at_k(pool, ranking, method=name) for name, ranking in rankings.items()]


def enrichment_by_feature(
    results: Sequence[AssayResult],
    features: Mapping[str, dict[str, float]],
    feature_name: str,
) -> dict[str, float]:
    """Mean feature value in validated vs non-validated peptides.

    Cheap diagnostic for which ranking terms are actually carrying signal
    before there is enough data to refit weights.
    """
    positive: list[float] = []
    negative: list[float] = []
    for result in labelled(results):
        values = features.get(result.key)
        if values is None or feature_name not in values:
            continue
        (positive if result.label == 1 else negative).append(values[feature_name])
    return {
        "positive_mean": sum(positive) / len(positive) if positive else float("nan"),
        "negative_mean": sum(negative) / len(negative) if negative else float("nan"),
        "n_positive": float(len(positive)),
        "n_negative": float(len(negative)),
    }
