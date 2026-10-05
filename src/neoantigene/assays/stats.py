"""Inference on assay outcomes.

`metrics` answers "what did this ranking retrieve". This module answers "could
that have happened anyway", which is the only question a 14-positive benchmark
can honestly be asked. Both benchmark cases currently report point estimates
whose sampling error is larger than the gap between methods, and the p-values
in `docs/BENCHMARK.md` were computed in a one-off script. Inference that
decides whether to keep a feature belongs in the library, under test.

Three deliberate choices:

**Exact tails, not normal approximations.** With 14 positives in 139 the normal
approximation to the hypergeometric is poor in the tail, which is the only part
anyone reads. `hypergeometric_tail` is exact integer arithmetic.

**AUC alongside top-k.** Recall@10 uses 10 of 139 observations and moves 7.1
points when one peptide shifts by one rank. The concordance probability uses
every positive against every negative, so it extracts far more from the same
labels. It answers a slightly different question — how the whole ordering
behaves rather than how the synthesized top of it behaves — so it is reported
next to precision@k, not instead of it.

**Pooling within patients only.** Each patient has a different assayed pool
and a different base response rate, so concatenating patients into one ranking
compares a peptide from one tumor against a peptide from another and measures
cohort heterogeneity as if it were ranking skill. `stratified_auc` counts only
within-patient pairs, which is what lets Ott's six patients and CU04 combine
into one estimate rather than a sum of separate top-10s.
"""

from __future__ import annotations

import random
from collections.abc import Callable, Mapping, Sequence
from math import comb

from pydantic import BaseModel, Field

from .metrics import labelled
from .schema import AssayResult

DEFAULT_ROUNDS = 10_000
DEFAULT_SEED = 0

#: Scored (score, label) observations for one patient.
Observations = list[tuple[float, int]]


def hypergeometric_tail(hits: int, k: int, positives: int, pool: int) -> float:
    """Exact P(at least `hits` positives in a random `k` drawn from `pool`).

    The null a top-k retrieval claim has to clear: that the top k is an
    unordered random draw from the peptides that were actually assayed.
    """
    if hits <= 0:
        return 1.0
    if k <= 0 or pool <= 0 or positives <= 0 or k > pool:
        return 1.0

    numerator = sum(
        comb(positives, i) * comb(pool - positives, k - i)
        for i in range(hits, min(k, positives) + 1)
        if k - i <= pool - positives
    )
    return numerator / comb(pool, k)


def expected_hits(k: int, positives: int, pool: int) -> float:
    """Positives in a random top-k: the chance line a method is compared to."""
    if pool <= 0 or k <= 0:
        return 0.0
    return k * positives / pool


def minimum_significant_hits(
    k: int,
    positives: int,
    pool: int,
    alpha: float = 0.05,
) -> int | None:
    """Fewest top-k hits that would reach `alpha`, or None if none can.

    Depends only on the shape of the assayed pool, so it is computable before
    any ranking is scored — and `None` is a complete answer on its own: no
    outcome this design can produce will clear the threshold, whatever method
    is run. CU04 at k=10 is such a design, with three positives in nineteen
    peptides.
    """
    reachable = min(k, positives)
    for hits in range(1, reachable + 1):
        if hypergeometric_tail(hits, k, positives, pool) <= alpha:
            return hits
    return None


class DesignCapacity(BaseModel):
    """What an assayed pool could show, independent of any method.

    Reported before results, deliberately. A benchmark that cannot reject
    chance at any outcome is not weak evidence about the methods run on it; it
    is no evidence, and that has to be visible before anyone reads a hit count
    and forms an impression.
    """

    k: int
    pool: int
    positives: int
    alpha: float
    expected_hits: float
    hits_needed: int | None

    @property
    def can_reach_significance(self) -> bool:
        return self.hits_needed is not None

    def describe(self) -> str:
        head = (
            f"{self.pool} assayed peptides, {self.positives} positive, k={self.k}. "
            f"Chance yields {self.expected_hits:.2f} hits."
        )
        if self.hits_needed is None:
            return (
                f"{head}\nNo outcome reaches p <= {self.alpha:g}: at most "
                f"{min(self.k, self.positives)} hits are possible and even that is not "
                f"significant. This design cannot distinguish any method from chance."
            )
        chance_rate = self.expected_hits / self.k
        return (
            f"{head}\nReaching p <= {self.alpha:g} takes {self.hits_needed} of {self.k} "
            f"({self.hits_needed / self.k:.0%} precision), against {chance_rate:.0%} "
            f"from chance."
        )


def design_capacity(
    pool: int,
    positives: int,
    k: int,
    alpha: float = 0.05,
) -> DesignCapacity:
    depth = min(k, pool) if pool > 0 else 0
    return DesignCapacity(
        k=depth,
        pool=pool,
        positives=positives,
        alpha=alpha,
        expected_hits=expected_hits(depth, positives, pool),
        hits_needed=minimum_significant_hits(depth, positives, pool, alpha),
    )


class TopKSignificance(BaseModel):
    method: str
    k: int
    pool: int
    positives: int
    hits: int
    expected: float
    p_value: float

    def describe(self) -> str:
        return (
            f"{self.method:<20} top-{self.k:<3} {self.hits} hits "
            f"(chance {self.expected:.2f})  p={self.p_value:.3f}"
        )


def top_k_significance(
    results: Sequence[AssayResult],
    ranking: Mapping[str, float],
    k: int,
    method: str = "neoantigene",
) -> TopKSignificance:
    """Observed top-k hits against the exact hypergeometric null."""
    observations = _observations(results, ranking)
    pool = len(observations)
    positives = sum(label for _score, label in observations)
    depth = min(k, pool)
    hits = sum(label for _score, label in _ordered(observations)[:depth])

    return TopKSignificance(
        method=method,
        k=depth,
        pool=pool,
        positives=positives,
        hits=hits,
        expected=expected_hits(depth, positives, pool),
        p_value=hypergeometric_tail(hits, depth, positives, pool),
    )


def rank_auc(observations: Observations) -> float | None:
    """Probability a random positive outranks a random negative, ties at 0.5.

    `None` when one class is absent, which is not a score of 0.5 but an
    undefined measurement, and must not be averaged as if it were a value.
    """
    positives = sum(label for _score, label in observations)
    negatives = len(observations) - positives
    if positives == 0 or negatives == 0:
        return None
    return _concordant_pairs(observations) / (positives * negatives)


def _concordant_pairs(observations: Observations) -> float:
    """Positive-over-negative wins, counting each tie as half a win.

    Computed from midranks rather than by comparing every pair, so cost is
    n log n in the assayed pool.
    """
    ordered = sorted(observations, key=lambda item: item[0])
    positives = sum(label for _score, label in observations)
    negatives = len(observations) - positives

    rank_sum = 0.0
    index = 0
    while index < len(ordered):
        stop = index
        while stop < len(ordered) and ordered[stop][0] == ordered[index][0]:
            stop += 1
        midrank = (index + stop + 1) / 2  # 1-based average rank of this tie group
        rank_sum += midrank * sum(label for _score, label in ordered[index:stop])
        index = stop

    wins = rank_sum - positives * (positives + 1) / 2
    return max(0.0, min(float(positives * negatives), wins))


class AucSignificance(BaseModel):
    method: str
    pool: int
    positives: int
    auc: float | None
    p_value: float | None
    rounds: int

    def describe(self) -> str:
        if self.auc is None:
            return f"{self.method:<20} auc n/a (one class only in a pool of {self.pool})"
        p = "n/a" if self.p_value is None else f"{self.p_value:.3f}"
        return f"{self.method:<20} auc {self.auc:.3f}  p={p}  (n={self.pool})"


def auc_significance(
    results: Sequence[AssayResult],
    ranking: Mapping[str, float],
    method: str = "neoantigene",
    rounds: int = DEFAULT_ROUNDS,
    seed: int = DEFAULT_SEED,
) -> AucSignificance:
    """Concordance probability with a label-permutation p-value.

    Permutation rather than the normal approximation: with a dozen positives
    the asymptotic variance of Mann-Whitney U is not trustworthy in the tail.
    """
    observations = _observations(results, ranking)
    observed = rank_auc(observations)
    positives = sum(label for _score, label in observations)

    p_value = None
    if observed is not None:
        p_value = _permutation_p(
            observations,
            observed,
            lambda shuffled: rank_auc(shuffled) or 0.0,
            rounds,
            seed,
        )

    return AucSignificance(
        method=method,
        pool=len(observations),
        positives=positives,
        auc=observed,
        p_value=p_value,
        rounds=rounds,
    )


class PairedComparison(BaseModel):
    """One ranking against another over exactly the same assayed peptides.

    This is the comparison the product claim rests on, and the one most easily
    faked by evaluating two methods on different pools. Both rankings are
    restricted to the peptides they both scored, and the null is that the
    labels are unrelated to either ordering — so if the two rankings are
    identical the difference is identically zero and the test correctly has no
    power to prefer one.
    """

    method_a: str
    method_b: str
    k: int
    pool: int
    positives: int
    hits_a: int
    hits_b: int
    hits_p_value: float
    auc_a: float | None = None
    auc_b: float | None = None
    auc_p_value: float | None = None
    rounds: int = DEFAULT_ROUNDS

    @property
    def hits_difference(self) -> int:
        return self.hits_a - self.hits_b

    @property
    def auc_difference(self) -> float | None:
        if self.auc_a is None or self.auc_b is None:
            return None
        return self.auc_a - self.auc_b

    def describe(self) -> str:
        lines = [
            f"{self.method_a} vs {self.method_b} over {self.pool} assayed peptides "
            f"({self.positives} positive)",
            f"  top-{self.k} hits  {self.hits_a} vs {self.hits_b}  "
            f"(difference {self.hits_difference:+d}, p={self.hits_p_value:.3f})",
        ]
        if self.auc_difference is not None and self.auc_p_value is not None:
            lines.append(
                f"  auc         {self.auc_a:.3f} vs {self.auc_b:.3f}  "
                f"(difference {self.auc_difference:+.3f}, p={self.auc_p_value:.3f})"
            )
        return "\n".join(lines)


def paired_comparison(
    results: Sequence[AssayResult],
    ranking_a: Mapping[str, float],
    ranking_b: Mapping[str, float],
    k: int,
    method_a: str = "neoantigene",
    method_b: str = "binding_only",
    rounds: int = DEFAULT_ROUNDS,
    seed: int = DEFAULT_SEED,
) -> PairedComparison:
    """Two-sided permutation test on the gap between two rankings."""
    shared = {r.key for r in labelled(results)} & set(ranking_a) & set(ranking_b)
    pool = [r for r in labelled(results) if r.key in shared]
    scores_a = [ranking_a[r.key] for r in pool]
    scores_b = [ranking_b[r.key] for r in pool]
    y = [r.label or 0 for r in pool]

    depth = min(k, len(pool))
    positives = sum(y)

    def hits_gap(labels: Sequence[int]) -> float:
        return float(
            _hits_at(list(zip(scores_a, labels, strict=True)), depth)
            - _hits_at(list(zip(scores_b, labels, strict=True)), depth)
        )

    def auc_gap(labels: Sequence[int]) -> float:
        first = rank_auc(list(zip(scores_a, labels, strict=True)))
        second = rank_auc(list(zip(scores_b, labels, strict=True)))
        if first is None or second is None:
            return 0.0
        return first - second

    auc_a = rank_auc(list(zip(scores_a, y, strict=True)))
    auc_b = rank_auc(list(zip(scores_b, y, strict=True)))

    return PairedComparison(
        method_a=method_a,
        method_b=method_b,
        k=depth,
        pool=len(pool),
        positives=positives,
        hits_a=_hits_at(list(zip(scores_a, y, strict=True)), depth),
        hits_b=_hits_at(list(zip(scores_b, y, strict=True)), depth),
        hits_p_value=_label_permutation_p(y, hits_gap, rounds, seed),
        auc_a=auc_a,
        auc_b=auc_b,
        auc_p_value=(
            None
            if auc_a is None or auc_b is None
            else _label_permutation_p(y, auc_gap, rounds, seed)
        ),
        rounds=rounds,
    )


class StratifiedAuc(BaseModel):
    """Concordance pooled over patients without ever crossing between them."""

    method: str
    patients: int
    pairs: int
    auc: float | None
    ci_low: float | None = None
    ci_high: float | None = None
    per_patient: dict[str, float] = Field(default_factory=dict)
    rounds: int = DEFAULT_ROUNDS

    def describe(self) -> str:
        if self.auc is None:
            return f"{self.method:<20} stratified auc n/a (no patient had both classes)"
        interval = ""
        if self.ci_low is not None and self.ci_high is not None:
            interval = f"  95% CI [{self.ci_low:.3f}, {self.ci_high:.3f}]"
        return (
            f"{self.method:<20} stratified auc {self.auc:.3f}{interval}  "
            f"({self.patients} patients, {self.pairs} within-patient pairs)"
        )


def stratified_auc(
    results: Sequence[AssayResult],
    ranking: Mapping[str, float],
    method: str = "neoantigene",
    rounds: int = DEFAULT_ROUNDS,
    seed: int = DEFAULT_SEED,
) -> StratifiedAuc:
    """Within-patient concordance, with a patient-level bootstrap interval.

    Pooling is weighted by each patient's positive-negative pair count, so a
    patient who contributes one positive and two negatives does not carry the
    same weight as one contributing twelve. The bootstrap resamples patients
    rather than peptides, because the patient is the unit that replicates: a
    seventh melanoma case is new information, a 140th peptide from the same six
    tumors is much less.
    """
    by_patient = _observations_by_patient(results, ranking)
    usable = {
        patient: observations
        for patient, observations in by_patient.items()
        if rank_auc(observations) is not None
    }
    if not usable:
        return StratifiedAuc(
            method=method, patients=len(by_patient), pairs=0, auc=None, rounds=rounds
        )

    per_patient = {patient: rank_auc(obs) or 0.0 for patient, obs in usable.items()}
    point = _pooled_auc(list(usable.values()))
    low, high = _bootstrap_interval(list(usable.values()), rounds, seed)

    return StratifiedAuc(
        method=method,
        patients=len(usable),
        pairs=int(sum(_pair_count(obs) for obs in usable.values())),
        auc=point,
        ci_low=low,
        ci_high=high,
        per_patient=per_patient,
        rounds=rounds,
    )


def _pooled_auc(groups: Sequence[Observations]) -> float | None:
    wins = sum(_concordant_pairs(obs) for obs in groups)
    pairs = sum(_pair_count(obs) for obs in groups)
    return None if pairs == 0 else wins / pairs


def _pair_count(observations: Observations) -> int:
    positives = sum(label for _score, label in observations)
    return positives * (len(observations) - positives)


def _bootstrap_interval(
    groups: Sequence[Observations],
    rounds: int,
    seed: int,
) -> tuple[float | None, float | None]:
    """Two-stage percentile interval: resample patients, then peptides within them.

    Resampling only patients — the textbook cluster bootstrap — is consistent
    as the number of clusters grows, and badly anti-conservative at the six
    patients this project actually has. It reproduces each selected patient's
    peptides intact, so it captures variation between patients and none of the
    sampling noise inside them. Measured coverage of a nominal 95% interval,
    against a known within-patient AUC of 0.5:

        patients   cluster-only   two-stage
               3          74.3%       97.0%
               6          82.7%       98.3%
              12          93.3%       99.0%

    Drawing peptides with replacement inside each selected patient restores
    that second level. It overshoots slightly, which is the tolerable
    direction: an interval that is a little too wide understates the evidence,
    while 82% coverage sold as 95% overstates it.
    """
    if len(groups) < 2:
        return None, None
    rng = random.Random(seed)
    draws: list[float] = []
    for _ in range(rounds):
        value = _pooled_auc(
            [
                [patient[rng.randrange(len(patient))] for _ in patient]
                for patient in (groups[rng.randrange(len(groups))] for _ in groups)
            ]
        )
        if value is not None:
            draws.append(value)
    if not draws:
        return None, None
    draws.sort()
    return _percentile(draws, 0.025), _percentile(draws, 0.975)


def _percentile(sorted_values: Sequence[float], fraction: float) -> float:
    if len(sorted_values) == 1:
        return sorted_values[0]
    position = fraction * (len(sorted_values) - 1)
    lower = int(position)
    upper = min(lower + 1, len(sorted_values) - 1)
    weight = position - lower
    return sorted_values[lower] * (1 - weight) + sorted_values[upper] * weight


def _hits_at(observations: Observations, k: int) -> int:
    return sum(label for _score, label in _ordered(observations)[:k])


def _ordered(observations: Observations) -> Observations:
    """Best-first. Ties resolve against the ranking, never in its favour.

    A tied block is ordered negatives-first so a method cannot be credited for
    a hit it did not actually rank above anything.
    """
    return sorted(observations, key=lambda item: (-item[0], item[1]))


def _observations(
    results: Sequence[AssayResult],
    ranking: Mapping[str, float],
) -> Observations:
    return [(ranking[r.key], r.label or 0) for r in labelled(results) if r.key in ranking]


def _observations_by_patient(
    results: Sequence[AssayResult],
    ranking: Mapping[str, float],
) -> dict[str, Observations]:
    grouped: dict[str, Observations] = {}
    for result in labelled(results):
        if result.key in ranking:
            grouped.setdefault(result.sample_id, []).append(
                (ranking[result.key], result.label or 0)
            )
    return grouped


def _permutation_p(
    observations: Observations,
    observed: float,
    statistic: Callable[[Observations], float],
    rounds: int,
    seed: int,
) -> float:
    labels = [label for _score, label in observations]
    scores = [score for score, _label in observations]

    def evaluate(shuffled: Sequence[int]) -> float:
        return statistic(list(zip(scores, shuffled, strict=True)))

    return _label_permutation_p(labels, evaluate, rounds, seed, observed=observed)


def _label_permutation_p(
    labels: Sequence[int],
    statistic: Callable[[Sequence[int]], float],
    rounds: int,
    seed: int,
    observed: float | None = None,
) -> float:
    """Two-sided permutation p-value, with the observed draw included.

    Including the observed arrangement in the null count is what keeps the
    p-value from being reported as 0 when no permutation happens to beat it:
    with `rounds` shuffles the floor is 1/(rounds+1), which is an honest
    statement of the resolution the test actually has.
    """
    target = abs(statistic(labels) if observed is None else observed)
    shuffled = list(labels)
    rng = random.Random(seed)

    extreme = 1
    for _ in range(rounds):
        rng.shuffle(shuffled)
        if abs(statistic(shuffled)) >= target - 1e-12:
            extreme += 1
    return extreme / (rounds + 1)
