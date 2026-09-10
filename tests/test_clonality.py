import pytest

from neoantigene.models import Variant, VariantClass
from neoantigene.variants.clonality import cancer_cell_fraction, with_ccf


def variant(**overrides) -> Variant:
    payload = {
        "chrom": "1",
        "pos": 100,
        "ref": "A",
        "alt": "T",
        "gene": "G",
        "transcript": "ENST1",
        "variant_class": VariantClass.MISSENSE,
        "protein_start": 10,
        "protein_end": 10,
        "aa_ref": "A",
        "aa_alt": "T",
    }
    return Variant(**{**payload, **overrides})


def test_diploid_pure_tumor_vaf_half_is_clonal():
    assert cancer_cell_fraction(vaf=0.5, purity=1.0) == pytest.approx(1.0)


def test_low_purity_inflates_ccf():
    low_purity = cancer_cell_fraction(vaf=0.2, purity=0.3)
    high_purity = cancer_cell_fraction(vaf=0.2, purity=0.9)
    assert low_purity is not None and high_purity is not None
    assert low_purity > high_purity


def test_ccf_is_clamped_to_one():
    assert cancer_cell_fraction(vaf=0.9, purity=0.4) == 1.0


def test_missing_vaf_returns_none():
    assert cancer_cell_fraction(vaf=None, purity=0.5) is None


@pytest.mark.parametrize("purity", [0.0, -0.2, 1.5])
def test_invalid_purity_raises(purity):
    with pytest.raises(ValueError, match="purity"):
        cancer_cell_fraction(vaf=0.3, purity=purity)


def test_non_positive_multiplicity_raises():
    with pytest.raises(ValueError, match="multiplicity"):
        cancer_cell_fraction(vaf=0.3, purity=0.5, multiplicity=0.0)


def test_amplified_locus_raises_ccf_for_same_vaf():
    diploid = cancer_cell_fraction(vaf=0.25, purity=0.6, tumor_copy_number=2)
    amplified = cancer_cell_fraction(vaf=0.25, purity=0.6, tumor_copy_number=6)
    assert amplified is not None and diploid is not None
    assert amplified > diploid


def test_with_ccf_does_not_mutate_the_input():
    original = variant(dna_vaf=0.3)
    annotated = with_ccf(original, purity=0.6)
    assert original.ccf is None
    assert annotated.ccf is not None
    assert annotated.dna_vaf == original.dna_vaf


def test_with_ccf_handles_missing_vaf():
    assert with_ccf(variant(), purity=0.6).ccf is None
