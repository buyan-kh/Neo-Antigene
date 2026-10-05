"""The evaluation harness.

What is under test is not arithmetic — `test_stats.py` covers that — but
whether the report can be misread as a win it did not earn. Each test below
drives the harness to one verdict and asserts the verdict is the honest one,
including the two cases where the honest answer is "this tells you nothing".

The other property that matters is that it works on rankings this package did
not produce. A harness only usable on its author's own output measures
nothing, so none of these fixtures run the pipeline.
"""

import pytest

from neoantigene.assays.harness import SourceFile, build_report
from neoantigene.assays.schema import AssayCall, AssayResult, AssayType

ALLELE = "HLA-B*18:01"


def _labels(spec: dict[str, list[int]]) -> list[AssayResult]:
    out: list[AssayResult] = []
    for patient, calls in spec.items():
        for index, label in enumerate(calls):
            out.append(
                AssayResult(
                    sample_id=patient,
                    peptide=f"{patient}PEP{index:03d}",
                    allele=ALLELE,
                    assay=AssayType.IFNG_ELISPOT,
                    call=AssayCall.POSITIVE if label else AssayCall.NEGATIVE,
                )
            )
    return out


def _perfect(results: list[AssayResult]) -> dict[str, float]:
    return {r.key: (1.0 if r.label else 0.0) for r in results}


def _inverted(results: list[AssayResult]) -> dict[str, float]:
    return {r.key: (0.0 if r.label else 1.0) for r in results}


def _noise(results: list[AssayResult], seed: int = 0) -> dict[str, float]:
    import random

    rng = random.Random(seed)
    return {r.key: rng.random() for r in results}


def _report(results, rankings, **kwargs):
    defaults = {
        "run_id": "20260101T000000Z-test",
        "package_version": "0.0.0-test",
        "rounds": 400,
        "k": 10,
    }
    return build_report(results, rankings, **{**defaults, **kwargs})


class TestVerdict:
    def test_a_pool_that_cannot_reject_chance_says_so_first(self):
        """CU04-shaped: 19 peptides, 3 positives, nothing is reachable."""
        results = _labels({"CU04": [1, 1, 1] + [0] * 16})
        report = _report(results, {"a": _noise(results), "b": _noise(results, 1)}, baseline="b")

        assert report.summary().startswith("INCONCLUSIVE")
        assert not report.supports_a_claim
        assert not report.capacity.can_reach_significance

    def test_a_perfect_ranking_on_an_underpowered_pool_is_not_sold_as_a_win(self):
        """The case that would be most tempting to quote.

        A perfect ordering still cannot produce a significant precision@10 on
        19 peptides with 3 positives, so the verdict must scope the claim to
        the ordering test and forbid the hit rate.

        The comparator is the exactly-reversed ranking rather than a random
        one, to make the ordering gap maximal and the test deterministic. With
        3 positives against 16 negatives there are only 48 pairs, and a
        perfect-versus-random comparison lands near p = 0.06 — correct
        behaviour for a pool this size, but too close to the threshold to
        assert on.
        """
        results = _labels({"CU04": [1, 1, 1] + [0] * 16})
        report = _report(
            results, {"ours": _perfect(results), "theirs": _inverted(results)}, baseline="theirs"
        )

        assert report.retrieval[0].recall_at_10 == 1.0, "fixture must be a perfect ranking"
        assert report.summary().startswith("ORDERING ONLY")
        assert "quote no hit rate" in report.summary()
        assert any("Do not quote a hit rate" in note for note in report.caveats())

    def test_a_powered_pool_with_no_method_above_chance_reports_a_result(self):
        """Absence of signal in a pool that could have shown it is a finding."""
        results = _labels({f"PT{p}": [1] * 14 + [0] * 125 for p in range(1)})
        report = _report(results, {"a": _inverted(results)}, k=20)

        assert report.capacity.can_reach_significance
        assert report.summary().startswith("NO METHOD CLEARS CHANCE")
        assert "not a missing measurement" in report.summary()

    def test_above_chance_is_distinguished_from_above_baseline(self):
        """Two near-identical strong rankings clear chance but cannot separate."""
        results = _labels({"PT1": [1] * 14 + [0] * 125})
        shared = _perfect(results)
        report = _report(results, {"ours": shared, "theirs": dict(shared)}, baseline="theirs", k=20)

        assert report.above_chance
        assert not report.separated
        assert report.summary().startswith("ABOVE CHANCE, NOT ABOVE BASELINE")
        assert "not the product claim" in report.summary()

    def test_a_real_separation_is_reported_with_a_replication_caveat(self):
        results = _labels({f"PT{p}": [1] * 5 + [0] * 45 for p in range(4)})
        report = _report(
            results, {"ours": _perfect(results), "theirs": _noise(results)}, baseline="theirs", k=20
        )

        assert report.separated == ["ours"]
        assert report.summary().startswith("SEPARATION FOUND")
        assert "independent cohort" in report.summary()


class TestScoping:
    def test_methods_are_restricted_to_the_peptides_they_all_scored(self):
        """A method must not be credited on peptides a rival never ranked."""
        results = _labels({"PT1": [1, 1, 0, 0, 0, 0]})
        full = _perfect(results)
        partial = {results[0].key: 1.0, results[2].key: 0.0}
        report = _report(results, {"full": full, "partial": partial}, baseline="partial")

        assert report.capacity.pool == 2
        assert all(item.assayed == 2 for item in report.retrieval)

    def test_an_unknown_baseline_is_ignored_rather_than_silently_substituted(self):
        results = _labels({"PT1": [1, 1] + [0] * 10})
        report = _report(results, {"a": _perfect(results)}, baseline="nope")

        assert report.baseline is None
        assert report.comparisons == []
        assert any("No baseline given" in note for note in report.caveats())

    def test_it_requires_at_least_one_ranking(self):
        with pytest.raises(ValueError, match="at least one ranking"):
            _report(_labels({"PT1": [1, 0]}), {})

    def test_indeterminate_labels_are_counted_but_never_scored(self):
        results = _labels({"PT1": [1, 0, 0, 0]})
        results.append(
            AssayResult(
                sample_id="PT1",
                peptide="PT1PEP999",
                allele=ALLELE,
                assay=AssayType.IFNG_ELISPOT,
                call=AssayCall.INDETERMINATE,
            )
        )
        report = _report(results, {"a": {r.key: 1.0 for r in results}})

        assert report.label_counts["indeterminate"] == 1
        assert report.capacity.pool == 4


class TestProvenance:
    def test_input_digests_are_carried_into_the_report(self):
        results = _labels({"PT1": [1, 0, 0, 0]})
        sources = [SourceFile(role="labels", path="/x/labels.tsv", digest="sha256:abc")]
        report = _report(results, {"a": _perfect(results)}, sources=sources)

        assert report.sources[0].digest == "sha256:abc"
        assert "sha256:abc" in report.to_markdown()

    def test_the_markdown_leads_with_the_verdict(self):
        results = _labels({"CU04": [1, 1, 1] + [0] * 16})
        report = _report(results, {"a": _noise(results)})
        body = report.to_markdown()

        heading, verdict = body.splitlines()[0], body.splitlines()[2]
        assert heading.startswith("# ")
        assert verdict.strip("*") == report.summary()

    def test_capacity_is_stated_before_any_result_table(self):
        """Order matters: a reader who sees a hit count first has anchored."""
        results = _labels({"CU04": [1, 1, 1] + [0] * 16})
        body = _report(results, {"a": _perfect(results)}).to_markdown()

        assert body.index("What this design could detect") < body.index("## Retrieval")
        assert body.index("## Retrieval") < body.index("## Reading this honestly")

    def test_small_cohorts_are_flagged_in_the_caveats(self):
        results = _labels({"PT1": [1, 1] + [0] * 20})
        notes = " ".join(_report(results, {"a": _perfect(results)}).caveats())

        assert "positives" in notes
        assert "patient" in notes
        assert "never assayed are excluded" in notes
