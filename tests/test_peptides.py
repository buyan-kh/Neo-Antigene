import pytest

from neoantigene.models import Variant, VariantClass
from neoantigene.peptides.generate import (
    FrameshiftNotSupported,
    ReferenceMismatch,
    TranscriptNotFound,
    build_mutant_protein,
    generate_for_variant,
)

KRAS_TRANSCRIPT = "ENST00000256078"
SYNTH_TRANSCRIPT = "ENST90000000001"


def kras_g12d(**overrides) -> Variant:
    payload = {
        "chrom": "12",
        "pos": 25245350,
        "ref": "C",
        "alt": "T",
        "gene": "KRAS",
        "transcript": KRAS_TRANSCRIPT,
        "variant_class": VariantClass.MISSENSE,
        "protein_start": 12,
        "protein_end": 12,
        "aa_ref": "G",
        "aa_alt": "D",
    }
    return Variant(**{**payload, **overrides})


def synth_deletion() -> Variant:
    return Variant(
        chrom="1",
        pos=1000420,
        ref="CCTGGCTGTG",
        alt="C",
        gene="SYNTHA",
        transcript=SYNTH_TRANSCRIPT,
        variant_class=VariantClass.INFRAME_DEL,
        protein_start=41,
        protein_end=43,
        aa_ref="LAV",
        aa_alt="",
    )


class TestMutantProtein:
    def test_substitutes_a_single_residue(self, proteome):
        wildtype = proteome.get(KRAS_TRANSCRIPT)
        mutant, start, end = build_mutant_protein(kras_g12d(), wildtype)
        assert (start, end) == (11, 12)
        assert mutant[11] == "D"
        assert len(mutant) == len(wildtype)
        assert mutant[:11] == wildtype[:11]

    def test_reference_mismatch_is_rejected(self, proteome):
        with pytest.raises(ReferenceMismatch, match="does not match"):
            build_mutant_protein(kras_g12d(aa_ref="W"), proteome.get(KRAS_TRANSCRIPT))

    def test_position_beyond_transcript_is_rejected(self, proteome):
        variant = kras_g12d(protein_start=9999, protein_end=9999)
        with pytest.raises(ReferenceMismatch, match="outside transcript"):
            build_mutant_protein(variant, proteome.get(KRAS_TRANSCRIPT))

    def test_frameshift_is_explicitly_unsupported(self, proteome):
        variant = kras_g12d(variant_class=VariantClass.FRAMESHIFT)
        with pytest.raises(FrameshiftNotSupported, match="mutant CDS"):
            build_mutant_protein(variant, proteome.get(KRAS_TRANSCRIPT))

    def test_deletion_shortens_the_protein_at_the_junction(self, proteome):
        wildtype = proteome.get(SYNTH_TRANSCRIPT)
        mutant, start, end = build_mutant_protein(synth_deletion(), wildtype)
        assert len(mutant) == len(wildtype) - 3
        assert start == end == 40

    def test_pure_insertion_lengthens_without_deleting(self, proteome):
        wildtype = proteome.get(SYNTH_TRANSCRIPT)
        variant = Variant(
            chrom="1",
            pos=1000420,
            ref="C",
            alt="CGCTGCA",
            gene="SYNTHA",
            transcript=SYNTH_TRANSCRIPT,
            variant_class=VariantClass.INFRAME_INS,
            protein_start=41,
            protein_end=42,
            aa_ref="",
            aa_alt="AC",
        )
        mutant, start, end = build_mutant_protein(variant, wildtype)
        assert len(mutant) == len(wildtype) + 2
        assert (start, end) == (41, 43)
        assert mutant[:41] == wildtype[:41]
        assert mutant[43:] == wildtype[41:]


class TestGeneration:
    def test_missing_transcript_is_an_error(self, proteome):
        with pytest.raises(TranscriptNotFound):
            generate_for_variant(kras_g12d(transcript="ENSTNOPE", gene="NOPE"), proteome, [9])

    def test_missense_yields_every_register(self, proteome):
        candidates = generate_for_variant(kras_g12d(), proteome, lengths=[9])
        assert len(candidates) == 9
        assert all("D" in c.mutant_peptide for c in candidates)
        assert {c.mutation_position for c in candidates} == set(range(1, 10))

    def test_known_public_epitopes_are_recovered(self, proteome):
        candidates = generate_for_variant(kras_g12d(), proteome, lengths=[9, 10])
        peptides = {c.mutant_peptide for c in candidates}
        # HLA-C*08:02 and HLA-A*11:01 restricted KRAS G12D epitopes from the
        # adoptive-transfer literature.
        assert {"GADGVGKSA", "GADGVGKSAL", "VVGADGVGK"} <= peptides

    def test_wildtype_counterpart_differs_at_exactly_one_position(self, proteome):
        for candidate in generate_for_variant(kras_g12d(), proteome, lengths=[9]):
            assert candidate.wildtype_peptide is not None
            differences = [
                (a, b)
                for a, b in zip(candidate.mutant_peptide, candidate.wildtype_peptide, strict=True)
                if a != b
            ]
            assert differences == [("D", "G")]

    def test_flanks_come_from_the_mutant_protein(self, proteome):
        candidates = generate_for_variant(kras_g12d(), proteome, lengths=[9], flank_length=5)
        first = next(c for c in candidates if c.mutant_peptide == "GADGVGKSA")
        assert first.n_flank == "KLVVV"
        assert first.c_flank == "LTIQL"

    def test_peptides_are_unique_across_lengths(self, proteome):
        candidates = generate_for_variant(kras_g12d(), proteome, lengths=[8, 9, 10, 11])
        peptides = [c.mutant_peptide for c in candidates]
        assert len(peptides) == len(set(peptides))

    def test_deletion_windows_all_straddle_the_junction(self, proteome):
        variant = synth_deletion()
        wildtype = proteome.get(SYNTH_TRANSCRIPT)
        mutant, _, _ = build_mutant_protein(variant, wildtype)

        candidates = generate_for_variant(variant, proteome, lengths=[9])
        assert candidates
        for candidate in candidates:
            assert candidate.mutant_peptide in mutant
            # A window that did not span the junction would still exist in WT.
            assert candidate.mutant_peptide not in wildtype
            # No positional wild-type counterpart exists for a length change.
            assert candidate.wildtype_peptide is None

    def test_generation_does_not_mutate_the_variant(self, proteome):
        variant = kras_g12d()
        before = variant.model_dump()
        generate_for_variant(variant, proteome, lengths=[9])
        assert variant.model_dump() == before
