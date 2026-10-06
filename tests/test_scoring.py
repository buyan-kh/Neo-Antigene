import math

import pytest

from neoantigene.config import OutputConfig, ScoringWeights
from neoantigene.models import (
    GateFailure,
    PeptideCandidate,
    PresentationCall,
    ScoredCandidate,
    Variant,
    VariantClass,
)
from neoantigene.scoring.baselines import arbitrary, binding_only, expression_only
from neoantigene.scoring.features import (
    NEUTRAL,
    agretopicity_feature,
    expression_feature,
    presentation_feature,
    tumor_selectivity_feature,
    wt_dissimilarity_feature,
)
from neoantigene.scoring.immunogenicity import (
    agretopicity,
    dissimilarity_to_wildtype,
    tcr_contact_hydrophobicity,
)
from neoantigene.scoring.rank import contributions, logit, score_all, score_features, shortlist


def make_scored(
    position: int,
    score: float = 0.5,
    peptide: str = "AIVLIVLLV",
    allele: str = "HLA-A*02:01",
    features: dict[str, float] | None = None,
    gates: tuple[GateFailure, ...] = (),
    variant_class: VariantClass = VariantClass.MISSENSE,
    **call_kwargs,
) -> ScoredCandidate:
    variant = Variant(
        chrom="1",
        pos=position,
        ref="A",
        alt="T",
        gene="G",
        transcript="ENST1",
        variant_class=variant_class,
        protein_start=10,
        protein_end=10,
        aa_ref="A",
        aa_alt="T",
    )
    candidate = PeptideCandidate(variant=variant, mutant_peptide=peptide, mutation_position=3)
    return ScoredCandidate(
        candidate=candidate,
        allele=allele,
        mutant_call=PresentationCall(peptide=peptide, allele=allele, **call_kwargs),
        features=features or {},
        score=score,
        gate_failures=gates,
    )


def _value(measurement: float | None) -> float:
    """Narrow an optional measurement, failing loudly if it is absent."""
    assert measurement is not None
    return measurement


class TestImmunogenicity:
    def test_agretopicity_rewards_mutation_created_binding(self):
        assert agretopicity(50.0, 5000.0) == 100.0
        assert agretopicity(50.0, None) is None
        assert agretopicity(None, 5000.0) is None

    def test_missing_wildtype_is_undefined_not_maximal(self):
        """No positional counterpart means the quantity cannot be measured.

        It previously returned 1.0, the maximum. Together with the same default
        on `agretopicity` that handed every indel, frameshift and neo-ORF
        peptide +0.65 logit on no evidence.
        """
        assert dissimilarity_to_wildtype("SIINFEKL", None) is None
        assert dissimilarity_to_wildtype("SIINFEKL", "SIINFEKLA") is None

    def test_dissimilarity_ignores_anchor_only_changes(self):
        # Position 2 is an MHC anchor, not a TCR contact.
        assert dissimilarity_to_wildtype("SAINFEKL", "SIINFEKL") == 0.0
        assert _value(dissimilarity_to_wildtype("SIIWFEKL", "SIINFEKL")) > 0.0

    def test_radical_substitution_scores_above_conservative(self):
        conservative = _value(dissimilarity_to_wildtype("SIIDFEKL", "SIINFEKL"))
        radical = _value(dissimilarity_to_wildtype("SIIWFEKL", "SIINFEKL"))
        assert radical > conservative

    def test_single_radical_change_is_not_diluted(self):
        """Averaging over substituted positions only, so one big change counts."""
        assert _value(dissimilarity_to_wildtype("SIIWFEKL", "SIINFEKL")) > 0.5

    def test_hydrophobicity_is_bounded_and_ordered(self):
        hydrophobic = tcr_contact_hydrophobicity("AIVLIVLLV")
        hydrophilic = tcr_contact_hydrophobicity("ARDEKRDEK")
        assert 0.0 <= hydrophilic < hydrophobic <= 1.0


class TestFeatureTransforms:
    def test_expression_feature_saturates(self):
        assert expression_feature(None) == 0.5
        assert expression_feature(0.0) == 0.0
        assert expression_feature(10.0) < expression_feature(1000.0)
        assert expression_feature(10**9) == 1.0

    def test_tumor_selectivity_penalizes_normal_expression(self):
        assert tumor_selectivity_feature(0.0) == 1.0
        assert tumor_selectivity_feature(500.0) < tumor_selectivity_feature(1.0)
        assert tumor_selectivity_feature(None) == 0.5

    def test_presentation_feature_prefers_score_then_percentile(self):
        call = PresentationCall(
            peptide="AIVLIVLLV",
            allele="HLA-A*02:01",
            presentation_score=0.9,
            affinity_percentile=50.0,
        )
        assert presentation_feature(call) == pytest.approx(0.9)

        strong = PresentationCall(
            peptide="AIVLIVLLV", allele="HLA-A*02:01", affinity_percentile=0.05
        )
        weak = PresentationCall(peptide="AIVLIVLLV", allele="HLA-A*02:01", affinity_percentile=40.0)
        assert presentation_feature(strong) > presentation_feature(weak)

    def test_presentation_feature_falls_back_to_affinity(self):
        call = PresentationCall(peptide="AIVLIVLLV", allele="HLA-A*02:01", affinity_nm=20.0)
        assert presentation_feature(call) > 0.5

    def test_presentation_feature_is_zero_without_any_signal(self):
        assert presentation_feature(PresentationCall(peptide="A", allele="A0201")) == 0.0

    def test_agretopicity_feature_is_monotonic_and_bounded(self):
        values = [agretopicity_feature(r) for r in (0.01, 0.1, 1.0, 10.0, 1000.0)]
        assert values == sorted(values)
        assert all(0.0 <= v <= 1.0 for v in values)
        assert agretopicity_feature(1.0) == pytest.approx(0.5)

    def test_an_absent_comparator_scores_neutral_in_both_affected_features(self):
        """The imputation that would have made frameshifts self-validating.

        Both features defaulted to 1.0, their maximum, so a peptide with no
        wild-type counterpart collected 0.5 + 0.15 = 0.65 logit for free.
        """
        assert agretopicity_feature(None) == NEUTRAL
        assert wt_dissimilarity_feature(None) == NEUTRAL

    def test_the_absence_itself_is_reported_so_a_refit_can_price_it(self):
        weights = ScoringWeights()
        free = weights.agretopicity * 1.0 + weights.wt_dissimilarity * 1.0
        neutral = (weights.agretopicity + weights.wt_dissimilarity) * NEUTRAL
        assert free - neutral == pytest.approx(0.325)


class TestScoring:
    def test_score_is_bounded_and_increases_with_presentation(self):
        weights = ScoringWeights()
        low = score_features({"presentation": 0.1}, weights)
        high = score_features({"presentation": 1.0}, weights)
        assert 0.0 < low < high < 1.0

    def test_no_driver_gene_term_exists(self):
        """The absence of this feature is a product decision, not an oversight."""
        fields = set(ScoringWeights.model_fields)
        assert not {f for f in fields if "driver" in f or "gene" in f}

    def test_contributions_sum_with_bias_to_the_logit(self):
        weights = ScoringWeights()
        features = dict.fromkeys(("clonality", "expression", "presentation"), 0.4)
        total = weights.bias + sum(contributions(features, weights).values())
        assert total == pytest.approx(logit(features, weights))
        assert score_features(features, weights) == pytest.approx(1 / (1 + math.exp(-total)))

    def test_score_all_returns_new_objects(self):
        original = make_scored(1, score=0.0, features={"presentation": 0.9})
        rescored = score_all([original], ScoringWeights())
        assert original.score == 0.0
        assert rescored[0].score > 0.0
        assert rescored[0] is not original


class TestShortlist:
    def test_caps_candidates_per_variant(self):
        scored = [
            make_scored(variant_index, score=0.9 - index * 0.01, peptide=f"AIVLIVLL{residue}")
            for variant_index in (1, 2, 3)
            for index, residue in enumerate("VILMF")
        ]
        selected = shortlist(scored, OutputConfig(top_n=50, max_per_variant=2))
        assert len(selected) == 6
        counts: dict[str, int] = {}
        for item in selected:
            counts[item.candidate.variant.key] = counts.get(item.candidate.variant.key, 0) + 1
        assert set(counts.values()) == {2}

    def test_a_frameshift_gets_its_own_larger_cap(self):
        """Downstream of a shift, extra peptides are distinct hypotheses.

        The substitution cap exists to stop one variant filling the list with
        overlapping registers of a single hypothesis. That reasoning inverts
        for a neo-ORF tail, where the extra peptides are independent epitopes,
        so applying the substitution cap would discard the poly-epitope
        structure that makes frameshifts worth enumerating.
        """
        peptides = [f"AIVLIVLL{residue}" for residue in "VILMFWY"]
        missense = [make_scored(1, score=0.9 - i * 0.01, peptide=p) for i, p in enumerate(peptides)]
        frameshift = [
            make_scored(2, score=0.9 - i * 0.01, peptide=p, variant_class=VariantClass.FRAMESHIFT)
            for i, p in enumerate(peptides)
        ]
        config = OutputConfig(top_n=50, max_per_variant=2, max_per_frameshift_variant=5)

        selected = shortlist(missense + frameshift, config)
        by_class: dict[VariantClass, int] = {}
        for item in selected:
            key = item.candidate.variant.variant_class
            by_class[key] = by_class.get(key, 0) + 1

        assert by_class[VariantClass.MISSENSE] == 2
        assert by_class[VariantClass.FRAMESHIFT] == 5

    def test_lowering_the_frameshift_cap_recovers_the_old_behaviour(self):
        peptides = [f"AIVLIVLL{residue}" for residue in "VILMF"]
        frameshift = [
            make_scored(1, score=0.9 - i * 0.01, peptide=p, variant_class=VariantClass.FRAMESHIFT)
            for i, p in enumerate(peptides)
        ]
        config = OutputConfig(top_n=50, max_per_variant=2, max_per_frameshift_variant=2)
        assert len(shortlist(frameshift, config)) == 2

    def test_respects_top_n(self):
        scored = [make_scored(i, score=0.5) for i in range(1, 20)]
        assert len(shortlist(scored, OutputConfig(top_n=5, max_per_variant=1))) == 5

    def test_excludes_gated_candidates_by_default(self):
        passing = make_scored(1, score=0.9)
        gated = make_scored(2, score=0.99, gates=(GateFailure.SELF_PEPTIDE,))
        selected = shortlist([passing, gated], OutputConfig(top_n=10, max_per_variant=2))
        assert [s.candidate.variant.key for s in selected] == [passing.candidate.variant.key]

    def test_include_failed_reinstates_gated_candidates(self):
        passing = make_scored(1, score=0.9)
        gated = make_scored(2, score=0.99, gates=(GateFailure.SELF_PEPTIDE,))
        selected = shortlist(
            [passing, gated], OutputConfig(top_n=10, max_per_variant=2, include_failed=True)
        )
        assert len(selected) == 2

    def test_ties_break_deterministically(self):
        scored = [make_scored(i, score=0.5) for i in range(1, 10)]
        first = shortlist(scored, OutputConfig(top_n=5, max_per_variant=1))
        second = shortlist(list(reversed(scored)), OutputConfig(top_n=5, max_per_variant=1))
        assert [s.key for s in first] == [s.key for s in second]


class TestBaselines:
    def test_binding_only_prefers_lower_percentile(self):
        strong = make_scored(1, affinity_percentile=0.1)
        weak = make_scored(2, peptide="AIVLIVLLI", affinity_percentile=30.0)
        ranking = binding_only([strong, weak])
        assert ranking[strong.key] > ranking[weak.key]

    def test_binding_only_falls_back_to_affinity(self):
        strong = make_scored(1, affinity_nm=5.0)
        weak = make_scored(2, peptide="AIVLIVLLI", affinity_nm=8000.0)
        ranking = binding_only([strong, weak])
        assert ranking[strong.key] > ranking[weak.key]

    def test_expression_only_reads_the_feature(self):
        high = make_scored(1, features={"expression": 0.9})
        low = make_scored(2, peptide="AIVLIVLLI", features={"expression": 0.1})
        ranking = expression_only([high, low])
        assert ranking[high.key] > ranking[low.key]

    def test_arbitrary_is_deterministic_but_unrelated_to_score(self):
        peptides = [f"AIVLIVL{a}{b}" for a in "VILM" for b in "VILM"]
        scored = [
            make_scored(index, score=index / 20, peptide=peptide)
            for index, peptide in enumerate(peptides, start=1)
        ]
        first = arbitrary(scored)
        assert first == arbitrary(list(reversed(scored)))
        by_arbitrary = [s.key for s in sorted(scored, key=lambda s: -first[s.key])]
        by_score = [s.key for s in sorted(scored, key=lambda s: -s.score)]
        assert by_arbitrary != by_score
