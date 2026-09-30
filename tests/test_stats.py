"""Inference on assay outcomes.

These tests exist because the numbers this module produces are the ones that
decide whether a feature stays in the model. Two things are checked: that the
closed-form statistics agree with published reference values, and that the
resampling ones are calibrated — a permutation test that rejects at 20% when
it claims 5% would manufacture exactly the false confidence the rest of this
project works to avoid.
"""

import random

import pytest

from neoantigene.assays.schema import AssayCall, AssayResult, AssayType
from neoantigene.assays.stats import (
    auc_significance,
    expected_hits,
    hypergeometric_tail,
    paired_comparison,
    rank_auc,
    stratified_auc,
    top_k_significance,
)


def _results(labels_by_patient: dict[str, list[int]]) -> list[AssayResult]:
    out: list[AssayResult] = []
    for patient, labels in labels_by_patient.items():
        for index, label in enumerate(labels):
            out.append(
                AssayResult(
                    sample_id=patient,
                    peptide=f"{patient}PEP{index:04d}",
                    allele="HLA-A*02:01",
                    assay=AssayType.IFNG_ELISPOT,
                    call=AssayCall.POSITIVE if label else AssayCall.NEGATIVE,
                )
            )
    return out


def _cohort(
    patients: int,
    per_patient: int,
    rate: float,
    rng: random.Random,
    signal: float = 0.0,
) -> tuple[list[AssayResult], dict[str, float], dict[str, float]]:
    """A cohort plus two rankings. `signal` is how much ranking A really works."""
    results: list[AssayResult] = []
    good: dict[str, float] = {}
    other: dict[str, float] = {}
    for patient in range(patients):
        for index in range(per_patient):
            score = rng.random()
            label = rng.random() < rate + signal * (score - 0.5)
            result = AssayResult(
                sample_id=f"PT{patient}",
                peptide=f"P{patient}PEP{index:04d}",
                allele="HLA-A*02:01",
                assay=AssayType.IFNG_ELISPOT,
                call=AssayCall.POSITIVE if label else AssayCall.NEGATIVE,
            )
            results.append(result)
            good[result.key] = score
            other[result.key] = rng.random()
    return results, good, other


class TestHypergeometricTail:
    @pytest.mark.parametrize(
        ("hits", "k", "positives", "pool", "expected"),
        [
            # The exact values docs/BENCHMARK.md reports, which were previously
            # computed in a one-off script.
            (3, 10, 14, 139, 0.064),
            (2, 10, 14, 139, 0.265),
            (1, 10, 14, 139, 0.667),
            (3, 10, 3, 19, 0.124),
            (1, 10, 3, 19, 0.913),
        ],
    )
    def test_it_reproduces_the_published_benchmark_p_values(
        self, hits, k, positives, pool, expected
    ):
        assert hypergeometric_tail(hits, k, positives, pool) == pytest.approx(expected, abs=5e-4)

    def test_finding_nothing_is_never_evidence(self):
        assert hypergeometric_tail(0, 10, 14, 139) == 1.0

    def test_finding_every_positive_is_strong_evidence(self):
        assert hypergeometric_tail(7, 7, 7, 200) < 1e-6

    def test_it_is_monotone_in_hits(self):
        tails = [hypergeometric_tail(h, 20, 14, 139) for h in range(1, 10)]
        assert tails == sorted(tails, reverse=True)

    def test_degenerate_inputs_do_not_claim_significance(self):
        assert hypergeometric_tail(1, 0, 14, 139) == 1.0
        assert hypergeometric_tail(1, 10, 0, 139) == 1.0
        assert hypergeometric_tail(1, 200, 14, 139) == 1.0

    def test_expected_hits_is_the_chance_line(self):
        assert expected_hits(10, 14, 139) == pytest.approx(1.0072, abs=1e-4)
        assert expected_hits(10, 3, 19) == pytest.approx(1.5789, abs=1e-4)


class TestRankAuc:
    def test_a_perfect_ordering_scores_one(self):
        assert rank_auc([(0.9, 1), (0.8, 1), (0.2, 0), (0.1, 0)]) == 1.0

    def test_a_reversed_ordering_scores_zero(self):
        assert rank_auc([(0.9, 0), (0.8, 0), (0.2, 1), (0.1, 1)]) == 0.0

    def test_all_ties_score_one_half(self):
        assert rank_auc([(0.5, 1), (0.5, 0), (0.5, 1), (0.5, 0)]) == 0.5

    def test_a_partial_tie_gets_half_credit(self):
        # One positive above one negative, one pair tied.
        assert rank_auc([(0.9, 1), (0.5, 1), (0.5, 0), (0.1, 0)]) == pytest.approx(0.875)

    def test_one_class_is_undefined_rather_than_one_half(self):
        """A pool with no negatives has no measurable ordering."""
        assert rank_auc([(0.9, 1), (0.8, 1)]) is None
        assert rank_auc([(0.9, 0), (0.8, 0)]) is None
        assert rank_auc([]) is None


class TestTopKSignificance:
    def test_it_scores_the_top_of_the_ranking(self):
        results = _results({"PT1": [1, 1, 0, 0, 0, 0, 0, 0]})
        scores = [8.0, 7.0, 6.0, 5.0, 4.0, 3.0, 2.0, 1.0]
        ranking = {r.key: s for r, s in zip(results, scores, strict=True)}
        summary = top_k_significance(results, ranking, k=2)

        assert summary.hits == 2
        assert summary.pool == 8
        assert summary.positives == 2
        assert summary.p_value == pytest.approx(1 / 28)

    def test_k_is_clamped_to_the_assayed_pool(self):
        results = _results({"PT1": [1, 0, 0]})
        ranking = {r.key: 1.0 for r in results}
        assert top_k_significance(results, ranking, k=50).k == 3

    def test_ties_are_not_resolved_in_the_ranking_s_favour(self):
        """Every score equal must not be reported as a lucky top-k."""
        results = _results({"PT1": [1, 0, 0, 0]})
        ranking = {r.key: 0.5 for r in results}
        assert top_k_significance(results, ranking, k=1).hits == 0


class TestPairedComparison:
    def test_identical_rankings_cannot_prefer_either(self):
        results, ranking, _ = _cohort(3, 20, 0.2, random.Random(0))
        comparison = paired_comparison(results, ranking, ranking, k=10, rounds=200)

        assert comparison.hits_difference == 0
        assert comparison.auc_difference == 0.0
        assert comparison.hits_p_value == 1.0

    def test_a_real_advantage_is_detected(self):
        results, good, noise = _cohort(4, 40, 0.2, random.Random(1), signal=0.9)
        comparison = paired_comparison(results, good, noise, k=20, rounds=500, seed=1)

        assert comparison.auc_difference is not None
        assert comparison.auc_difference > 0.2
        assert comparison.auc_p_value is not None
        assert comparison.auc_p_value < 0.05

    def test_it_restricts_both_methods_to_the_shared_pool(self):
        """A method must not be scored on peptides the other never ranked."""
        results = _results({"PT1": [1, 1, 0, 0]})
        full = {r.key: 1.0 for r in results}
        partial = {results[0].key: 1.0, results[2].key: 0.0}
        assert paired_comparison(results, full, partial, k=10, rounds=50).pool == 2

    def test_the_p_value_floor_reflects_the_rounds_used(self):
        results, good, noise = _cohort(4, 40, 0.25, random.Random(2), signal=1.4)
        comparison = paired_comparison(results, good, noise, k=20, rounds=100, seed=2)
        assert comparison.auc_p_value is not None
        assert comparison.auc_p_value >= 1 / 101


class TestStratifiedAuc:
    def test_it_only_counts_within_patient_pairs(self):
        """Two patients, each perfectly ordered, but on disjoint score ranges.

        Pooled naively the low-scoring patient's positives sit below the
        high-scoring patient's negatives and the AUC drops. Stratified, both
        patients are perfect and so is the pooled estimate.
        """
        results = _results({"PT1": [1, 0], "PT2": [1, 0]})
        bands = {"PT1": (0.2, 0.1), "PT2": (0.9, 0.8)}
        ranking = {
            r.key: bands[r.sample_id][0] if r.label else bands[r.sample_id][1] for r in results
        }
        naive = rank_auc([(ranking[r.key], r.label or 0) for r in results])

        assert stratified_auc(results, ranking, rounds=50).auc == 1.0
        assert naive is not None and naive < 1.0

    def test_patients_with_one_class_are_dropped_not_scored_as_half(self):
        results = _results({"PT1": [1, 0, 0], "PT2": [0, 0, 0]})
        ranking = {r.key: 1.0 - i / 10 for i, r in enumerate(results)}
        summary = stratified_auc(results, ranking, rounds=50)

        assert summary.patients == 1
        assert set(summary.per_patient) == {"PT1"}

    def test_pooling_weights_patients_by_their_pair_count(self):
        """A patient contributing 1x2 pairs must not outweigh one contributing 2x4."""
        results = _results({"SMALL": [1, 0, 0], "BIG": [1, 1, 0, 0, 0, 0]})
        by_key = {r.key: r for r in results}
        # SMALL ordered perfectly, BIG ordered backwards.
        ranking = {}
        for key, result in by_key.items():
            if result.sample_id == "SMALL":
                ranking[key] = 1.0 if result.label else 0.0
            else:
                ranking[key] = 0.0 if result.label else 1.0
        summary = stratified_auc(results, ranking, rounds=50)

        assert summary.per_patient["SMALL"] == 1.0
        assert summary.per_patient["BIG"] == 0.0
        # 2 wins out of 2 + 8 pairs, not the unweighted mean of 0.5.
        assert summary.pairs == 10
        assert summary.auc == pytest.approx(0.2)

    def test_no_usable_patient_yields_no_estimate(self):
        results = _results({"PT1": [0, 0], "PT2": [0, 0]})
        summary = stratified_auc(results, {r.key: 1.0 for r in results}, rounds=50)
        assert summary.auc is None
        assert "n/a" in summary.describe()

    def test_the_interval_brackets_a_strong_effect(self):
        results, good, _ = _cohort(6, 30, 0.2, random.Random(3), signal=0.9)
        summary = stratified_auc(results, good, rounds=400, seed=3)

        assert summary.auc is not None
        assert summary.ci_low is not None and summary.ci_high is not None
        assert summary.ci_low < summary.auc < summary.ci_high
        assert summary.ci_low > 0.5, "a strong planted effect should exclude chance"


class TestCalibration:
    """The property that makes the rest of the module usable.

    Slow by nature — a calibration check needs many trials of many
    permutations — so these run under the `eval` marker with the other
    quality benchmarks rather than on every save.
    """

    @pytest.mark.eval
    def test_the_permutation_test_rejects_at_about_its_nominal_rate(self):
        rng = random.Random(21)
        rejections = 0
        trials = 150
        for _ in range(trials):
            results, _good, noise = _cohort(3, 30, 0.15, rng)
            p = auc_significance(results, noise, rounds=200, seed=rng.randrange(10**6)).p_value
            assert p is not None
            rejections += p < 0.05
        rate = rejections / trials
        assert rate < 0.12, (
            f"rejected the null in {rate:.1%} of trials where the ranking was pure noise; "
            f"a test claiming 5% must not be far above it"
        )

    @pytest.mark.eval
    def test_the_stratified_interval_covers_chance_when_there_is_no_signal(self):
        rng = random.Random(22)
        covered = 0
        trials = 100
        for _ in range(trials):
            results, _good, noise = _cohort(6, 30, 0.2, rng)
            summary = stratified_auc(results, noise, rounds=200, seed=rng.randrange(10**6))
            if summary.ci_low is not None and summary.ci_high is not None:
                covered += summary.ci_low <= 0.5 <= summary.ci_high
        coverage = covered / trials
        assert coverage >= 0.90, (
            f"a nominal 95% interval contained the true value in {coverage:.1%} of trials; "
            f"under-coverage here would be sold as confidence"
        )
