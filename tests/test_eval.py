"""Ranking quality benchmarks.

These encode the product thesis: on a cohort whose truth depends on
presentation *and* clonality *and* expression, a ranker that uses all three
should put more validated peptides in the top k than a NetMHCpan-class ranker
that sorts on binding alone.

What this does and does not prove:

  - It does prove the ranking machinery, feature assembly and metric recover
    a known multi-factor signal, and that the loop composes end to end. It
    also guards a real regression: give the weakly-evidenced sequence
    features too much weight and they dominate the top of the shortlist.
  - It does not prove the biology. Truth here is simulated (see
    `tests.fixtures.synthetic`). Only real assay data settles that, which is
    what `neoantigene evaluate` is for.

Results are pooled over patients and label-noise seeds. A single patient at
k=10 is ten observations, which is not enough to separate methods and would
make CI flake; pooling is the difference between a benchmark and an anecdote.
"""

from collections.abc import Sequence
from typing import Any

import pytest

from neoantigene.assays.metrics import compare_rankings, validation_rate_at_k
from neoantigene.io.manifest import Sample
from neoantigene.pipeline import run as run_pipeline
from neoantigene.presentation.null import NullBackend
from neoantigene.scoring.baselines import ALL_RANKERS, BASELINES

from .fixtures import simulate_assay_results

BASELINE = "binding_only"
LABEL_SEEDS = (11, 12, 13)
K = 40

pytestmark = pytest.mark.eval


@pytest.fixture(scope="module")
def replicates(cohort, benchmark_config):
    """One (results, rankings) pair per patient per label-noise seed."""
    out = []
    for patient in cohort.patients:
        sample = Sample.load(patient.manifest_path)
        scored = run_pipeline(sample, benchmark_config, backend=NullBackend()).all_scored
        sample_id = sample.sample_id
        rankings = {
            name: {f"{sample_id}|{key}": value for key, value in ranker(scored).items()}
            for name, ranker in ALL_RANKERS.items()
        }
        for seed in LABEL_SEEDS:
            out.append((simulate_assay_results(scored, sample_id, seed=seed), rankings))
    return out


def mean_rates(replicates: Sequence[Any], k: int) -> dict[str, float]:
    totals: dict[str, list[float]] = {name: [] for name in ALL_RANKERS}
    for results, rankings in replicates:
        for summary in compare_rankings(results, rankings, k, baseline=BASELINE):
            totals[summary.method].append(summary.validation_rate)
    return {name: sum(values) / len(values) for name, values in totals.items()}


def test_the_benchmark_has_power_to_distinguish_methods(replicates):
    """A base rate near 0 or near 1 would make every method look identical."""
    assert len(replicates) == 6
    for results, _ in replicates:
        assert len(results) >= 500
        base_rate = sum(1 for r in results if r.label == 1) / len(results)
        assert 0.10 < base_rate < 0.40


def test_full_model_beats_binding_only_baseline(replicates):
    rates = mean_rates(replicates, K)
    assert rates["neoantigene"] > rates[BASELINE], (
        f"neoantigene {rates['neoantigene']:.1%} did not beat {BASELINE} {rates[BASELINE]:.1%}"
    )


def test_full_model_beats_every_single_feature_baseline(replicates):
    rates = mean_rates(replicates, K)
    for name in BASELINES:
        assert rates["neoantigene"] > rates[name], (
            f"lost to {name}: {rates['neoantigene']:.1%} vs {rates[name]:.1%}"
        )


def test_lift_over_baseline_is_material(replicates):
    rates = mean_rates(replicates, K)
    lift = rates["neoantigene"] / rates[BASELINE]
    assert lift > 1.2, f"lift over {BASELINE} was only x{lift:.2f}"


def test_advantage_holds_across_cutoffs(replicates):
    """A win at one k that vanishes at another is a fluke, not a result."""
    for k in (10, 20, 40, 80, 160):
        rates = mean_rates(replicates, k)
        assert rates["neoantigene"] > rates[BASELINE], (
            f"lost at k={k}: {rates['neoantigene']:.1%} vs {rates[BASELINE]:.1%}"
        )


def test_speculative_features_do_not_dominate_the_top_of_the_list(replicates):
    """Regression guard on the shipped priors.

    Weighting agretopicity/foreignness/hydrophobicity comparably to
    presentation makes noise float to the top, which is worst exactly where it
    matters most — the handful of peptides a lab actually synthesizes.
    """
    rates = mean_rates(replicates, 10)
    assert rates["neoantigene"] > rates[BASELINE]


def test_every_method_is_scored_on_the_same_pool(replicates):
    for results, rankings in replicates:
        summaries = compare_rankings(results, rankings, K, baseline=BASELINE)
        assert len({s.assayed for s in summaries}) == 1
        assert len({s.positives for s in summaries}) == 1


def test_arbitrary_ranking_approximates_the_base_rate(replicates):
    """Sanity check on the metric itself: a meaningless order scores near chance."""
    rates = mean_rates(replicates, 160)
    base_rates = [
        sum(1 for r in results if r.label == 1) / len(results) for results, _ in replicates
    ]
    expected = sum(base_rates) / len(base_rates)
    assert abs(rates["arbitrary"] - expected) < 0.10


def test_precision_at_k_degrades_as_k_grows(replicates):
    """With a real signal, the densest hits are at the top of the list."""
    top = mean_rates(replicates, 20)["neoantigene"]
    deep = mean_rates(replicates, 400)["neoantigene"]
    assert top > deep


def test_benchmark_is_reproducible(replicates):
    """Same inputs, same seed, same numbers — CI cannot flake on this."""
    first = mean_rates(replicates, K)
    assert first == mean_rates(replicates, K)


def test_metric_handles_a_ranking_over_untested_peptides(replicates):
    results, rankings = replicates[0]
    tested = {r.key for r in results}
    padded = dict(rankings["neoantigene"]) | {"OTHER|ZZZZZZZZZ|HLA-A*02:01": 99.0}
    summary = validation_rate_at_k(results, padded, K)
    assert summary.assayed == len([r for r in results if r.label is not None])
    assert summary.k == K
    assert tested
