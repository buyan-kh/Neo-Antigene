"""Boundary-validation tests.

The point of putting pydantic at every stage boundary is that malformed data
cannot travel. These tests pin that behaviour, because a silent coercion here
becomes a wrong shortlist three stages later.
"""

import pytest
from pydantic import ValidationError

from neoantigene.models import (
    GateFailure,
    PeptideCandidate,
    PresentationCall,
    ScoredCandidate,
    Variant,
    VariantClass,
)


def make_variant(**overrides) -> Variant:
    payload = {
        "chrom": "12",
        "pos": 25245350,
        "ref": "C",
        "alt": "T",
        "gene": "KRAS",
        "transcript": "ENST00000256078",
        "variant_class": VariantClass.MISSENSE,
        "protein_start": 12,
        "protein_end": 12,
        "aa_ref": "G",
        "aa_alt": "D",
    }
    return Variant(**{**payload, **overrides})


def make_candidate(**overrides) -> PeptideCandidate:
    payload = {
        "variant": make_variant(),
        "mutant_peptide": "GADGVGKSA",
        "wildtype_peptide": "GAGGVGKSA",
        "mutation_position": 3,
    }
    return PeptideCandidate(**{**payload, **overrides})


class TestVariant:
    def test_models_are_frozen(self):
        variant = make_variant()
        with pytest.raises(ValidationError):
            variant.ccf = 0.5  # type: ignore[misc]

    def test_model_copy_is_the_mutation_path(self):
        variant = make_variant()
        updated = variant.model_copy(update={"ccf": 0.9})
        assert variant.ccf is None
        assert updated.ccf == 0.9

    @pytest.mark.parametrize("field", ["dna_vaf", "rna_vaf", "ccf", "population_af"])
    @pytest.mark.parametrize("value", [-0.1, 1.4])
    def test_fractions_must_be_in_unit_interval(self, field, value):
        with pytest.raises(ValidationError):
            make_variant(**{field: value})

    def test_protein_span_must_be_ordered(self):
        with pytest.raises(ValidationError, match="precedes"):
            make_variant(protein_start=40, protein_end=12)

    def test_position_must_be_positive(self):
        with pytest.raises(ValidationError):
            make_variant(pos=0)

    def test_amino_acids_reject_junk(self):
        with pytest.raises(ValidationError, match="non-residue"):
            make_variant(aa_alt="Z1")

    def test_stop_and_unknown_residues_are_accepted_on_variants(self):
        assert make_variant(aa_alt="*").aa_alt == "*"
        assert make_variant(aa_alt="X").aa_alt == "X"

    def test_extra_fields_are_rejected(self):
        with pytest.raises(ValidationError):
            make_variant(is_driver_gene=True)

    def test_hgvsp_renders_each_variant_class(self):
        assert make_variant().hgvsp_short == "KRAS p.G12D"
        deletion = make_variant(
            variant_class=VariantClass.INFRAME_DEL,
            protein_start=41,
            protein_end=43,
            aa_ref="LAV",
            aa_alt="",
        )
        assert deletion.hgvsp_short == "KRAS p.LAV41_43del"
        insertion = make_variant(
            variant_class=VariantClass.INFRAME_INS,
            protein_start=41,
            protein_end=42,
            aa_ref="",
            aa_alt="AC",
        )
        assert insertion.hgvsp_short == "KRAS p.41_42insAC"


class TestPeptideCandidate:
    def test_length_is_derived_not_stored(self):
        assert make_candidate().length == 9

    def test_peptides_reject_non_residues(self):
        with pytest.raises(ValidationError, match="non-residue"):
            make_candidate(mutant_peptide="GADGVGKS*")

    def test_peptide_length_bounds_are_enforced(self):
        with pytest.raises(ValidationError):
            make_candidate(mutant_peptide="GADGVG", wildtype_peptide=None)
        with pytest.raises(ValidationError):
            make_candidate(mutant_peptide="G" * 16, wildtype_peptide=None)

    def test_wildtype_must_share_the_register(self):
        with pytest.raises(ValidationError, match="register"):
            make_candidate(wildtype_peptide="GAGGVGKSAL")

    def test_wildtype_identical_to_mutant_is_rejected(self):
        with pytest.raises(ValidationError, match="identical"):
            make_candidate(wildtype_peptide="GADGVGKSA")

    def test_mutation_position_must_fall_inside_the_peptide(self):
        with pytest.raises(ValidationError, match="outside peptide"):
            make_candidate(mutation_position=12)


class TestPresentationCall:
    def test_allele_is_normalized_on_construction(self):
        call = PresentationCall(peptide="GADGVGKSA", allele="a0201")
        assert call.allele == "HLA-A*02:01"

    def test_invalid_allele_is_rejected(self):
        with pytest.raises(ValidationError):
            PresentationCall(peptide="GADGVGKSA", allele="banana")

    def test_percentile_range_is_enforced(self):
        with pytest.raises(ValidationError):
            PresentationCall(peptide="GADGVGKSA", allele="A0201", affinity_percentile=140.0)

    def test_presentation_score_is_a_probability(self):
        with pytest.raises(ValidationError):
            PresentationCall(peptide="GADGVGKSA", allele="A0201", presentation_score=1.2)

    def test_negative_affinity_is_rejected(self):
        with pytest.raises(ValidationError):
            PresentationCall(peptide="GADGVGKSA", allele="A0201", affinity_nm=-5.0)


class TestScoredCandidate:
    def _scored(self, **overrides) -> ScoredCandidate:
        payload = {
            "candidate": make_candidate(),
            "allele": "HLA-A*02:01",
            "mutant_call": PresentationCall(peptide="GADGVGKSA", allele="HLA-A*02:01"),
        }
        return ScoredCandidate(**{**payload, **overrides})

    def test_with_score_returns_a_new_object(self):
        original = self._scored()
        updated = original.with_score(0.8)
        assert original.score == 0.0
        assert updated.score == 0.8

    def test_score_must_be_a_probability(self):
        with pytest.raises(ValidationError):
            self._scored(score=1.5)

    def test_passed_reflects_gate_failures(self):
        assert self._scored().passed
        assert not self._scored(gate_failures=(GateFailure.SELF_PEPTIDE,)).passed

    def test_gate_failures_must_be_known_reasons(self):
        with pytest.raises(ValidationError):
            self._scored(gate_failures=("vibes",))
