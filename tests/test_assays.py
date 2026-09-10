import pytest
from pydantic import ValidationError

from neoantigene.assays.metrics import (
    compare_rankings,
    enrichment_by_feature,
    validation_rate_at_k,
)
from neoantigene.assays.schema import AssayCall, AssayResult, AssayType, read_results


def result(peptide: str, call: AssayCall, allele: str = "HLA-A*02:01") -> AssayResult:
    return AssayResult(
        sample_id="S1",
        peptide=peptide,
        allele=allele,
        assay=AssayType.IFNG_ELISPOT,
        call=call,
    )


class TestSchema:
    def test_label_mapping(self):
        assert result("AAAAAAAAA", AssayCall.POSITIVE).label == 1
        assert result("AAAAAAAAA", AssayCall.NEGATIVE).label == 0
        assert result("AAAAAAAAA", AssayCall.INDETERMINATE).label is None

    def test_allele_is_normalized(self):
        assert result("AAAAAAAAA", AssayCall.POSITIVE, allele="a0201").allele == "HLA-A*02:01"

    def test_unknown_call_is_rejected(self):
        with pytest.raises(ValidationError):
            AssayResult(
                sample_id="S1",
                peptide="AAAAAAAAA",
                assay=AssayType.IFNG_ELISPOT,
                call="probably?",  # type: ignore[arg-type]
            )

    def test_key_handles_missing_allele(self):
        record = AssayResult(
            sample_id="S1",
            peptide="AAAAAAAAA",
            assay=AssayType.IFNG_ELISPOT,
            call=AssayCall.POSITIVE,
        )
        assert record.key == "S1|AAAAAAAAA|*"

    def test_round_trips_through_tsv(self, tmp_path):
        path = tmp_path / "results.tsv"
        path.write_text(
            "sample_id\tpeptide\tallele\tassay\tcall\teffect_size\treplicate_count\n"
            "S1\tAAAAAAAAA\tA0201\tifng_elispot\tpositive\t120.5\t3\n"
            "S1\tCCCCCCCCC\t.\tifng_elispot\tnegative\t\t\n"
        )
        records = read_results(path)
        assert [r.label for r in records] == [1, 0]
        assert records[0].allele == "HLA-A*02:01"
        assert records[1].allele is None
        assert records[0].effect_size == 120.5


class TestValidationRate:
    def test_counts_only_the_top_k(self):
        results = [
            result("AAAAAAAAA", AssayCall.POSITIVE),
            result("CCCCCCCCC", AssayCall.POSITIVE),
            result("DDDDDDDDD", AssayCall.NEGATIVE),
            result("EEEEEEEEE", AssayCall.NEGATIVE),
        ]
        ranking = {r.key: score for r, score in zip(results, [0.9, 0.8, 0.7, 0.6], strict=True)}
        summary = validation_rate_at_k(results, ranking, k=2)
        assert summary.hits_in_top_k == 2
        assert summary.validation_rate == 1.0
        assert summary.assayed == 4
        assert summary.positives == 2

    def test_indeterminate_calls_are_excluded_not_counted_as_negative(self):
        results = [
            result("AAAAAAAAA", AssayCall.POSITIVE),
            result("CCCCCCCCC", AssayCall.INDETERMINATE),
        ]
        summary = validation_rate_at_k(results, {r.key: 0.5 for r in results}, k=2)
        assert summary.assayed == 1
        assert summary.validation_rate == 1.0

    def test_untested_peptides_are_not_counted_as_negatives(self):
        results = [result("AAAAAAAAA", AssayCall.POSITIVE)]
        ranking = {results[0].key: 0.9, "S1|ZZZZZZZZZ|HLA-A*02:01": 0.95}
        summary = validation_rate_at_k(results, ranking, k=2)
        assert summary.k == 1
        assert summary.validation_rate == 1.0

    def test_empty_pool_is_not_a_crash(self):
        summary = validation_rate_at_k([], {}, k=10)
        assert summary.assayed == 0
        assert summary.validation_rate == 0.0

    def test_ties_are_broken_deterministically(self):
        results = [
            result("AAAAAAAAA", AssayCall.POSITIVE),
            result("CCCCCCCCC", AssayCall.NEGATIVE),
        ]
        ranking = dict.fromkeys((r.key for r in results), 0.5)
        rates = {validation_rate_at_k(results, ranking, k=1).validation_rate for _ in range(5)}
        assert len(rates) == 1


class TestComparison:
    @pytest.fixture
    def head_to_head(self):
        results = [
            result("AAAAAAAAA", AssayCall.POSITIVE),
            result("CCCCCCCCC", AssayCall.NEGATIVE),
            result("DDDDDDDDD", AssayCall.POSITIVE),
            result("EEEEEEEEE", AssayCall.NEGATIVE),
        ]
        keys = [r.key for r in results]
        return results, {
            "neoantigene": dict(zip(keys, [0.9, 0.2, 0.85, 0.1], strict=True)),
            "netmhcpan": dict(zip(keys, [0.1, 0.9, 0.2, 0.85], strict=True)),
        }

    def test_restricts_to_the_shared_assayed_pool(self, head_to_head):
        results, rankings = head_to_head
        summaries = {s.method: s for s in compare_rankings(results, rankings, k=2)}
        assert summaries["neoantigene"].validation_rate == 1.0
        assert summaries["netmhcpan"].validation_rate == 0.0
        assert summaries["neoantigene"].assayed == summaries["netmhcpan"].assayed

    def test_lift_is_reported_against_the_baseline(self, head_to_head):
        results, rankings = head_to_head
        rankings["netmhcpan"] = dict.fromkeys(rankings["netmhcpan"], 0.5)
        summaries = {
            s.method: s for s in compare_rankings(results, rankings, k=2, baseline="netmhcpan")
        }
        assert summaries["neoantigene"].baseline_rate == summaries["netmhcpan"].validation_rate

    def test_a_method_cannot_win_by_nominating_an_easier_set(self):
        """Peptides only one method ranked are dropped from both."""
        results = [
            result("AAAAAAAAA", AssayCall.POSITIVE),
            result("CCCCCCCCC", AssayCall.NEGATIVE),
        ]
        keys = [r.key for r in results]
        rankings = {
            "neoantigene": {keys[0]: 0.9},
            "netmhcpan": {keys[0]: 0.5, keys[1]: 0.9},
        }
        summaries = {s.method: s for s in compare_rankings(results, rankings, k=5)}
        assert {s.assayed for s in summaries.values()} == {1}

    def test_no_rankings_is_not_a_crash(self):
        assert compare_rankings([], {}, k=5) == []


def test_enrichment_by_feature_separates_positives_from_negatives():
    results = [
        result("AAAAAAAAA", AssayCall.POSITIVE),
        result("CCCCCCCCC", AssayCall.NEGATIVE),
        result("DDDDDDDDD", AssayCall.INDETERMINATE),
    ]
    features = {
        results[0].key: {"presentation": 0.9},
        results[1].key: {"presentation": 0.2},
        results[2].key: {"presentation": 0.99},
    }
    summary = enrichment_by_feature(results, features, "presentation")
    assert summary["positive_mean"] == pytest.approx(0.9)
    assert summary["negative_mean"] == pytest.approx(0.2)
    # The indeterminate row must not contribute to either arm.
    assert summary["n_positive"] == 1.0
    assert summary["n_negative"] == 1.0
