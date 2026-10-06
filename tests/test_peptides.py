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


#: An invented novel tail. No claim is made about it being a real sequence —
#: the tests assert structure, not biology. Deliberately free of residues that
#: would collide with the synthetic reference, so "not in wildtype" is a
#: meaningful assertion.
TAIL = "WYFMKWYFHKWYFQKWYFRKWYFNK"

#: A longer tail with no 9-mer in common with `TAIL`, so doubling the novel
#: sequence genuinely doubles the candidate count instead of being collapsed
#: by the generator's deduplication.
LONGER_TAIL = TAIL + "AGCDEAGCTEAGCSEAGCVEAGCIE"


def synth_frameshift(**overrides) -> Variant:
    payload = {
        "chrom": "1",
        "pos": 1000420,
        "ref": "CT",
        "alt": "C",
        "gene": "SYNTHA",
        "transcript": SYNTH_TRANSCRIPT,
        "variant_class": VariantClass.FRAMESHIFT,
        "protein_start": 41,
        "protein_end": 41,
        "aa_ref": "L",
        "aa_alt": "",
        "downstream_protein": TAIL,
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

    def test_frameshift_without_its_novel_tail_is_refused(self, proteome):
        """The tail is not derivable from the reference, so it must be supplied."""
        variant = kras_g12d(variant_class=VariantClass.FRAMESHIFT)
        with pytest.raises(FrameshiftNotSupported, match="novel tail"):
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


class TestFrameshift:
    """Neo-ORF enumeration.

    The payoff is poly-epitope: one frameshift can produce several
    independently immunogenic peptides, which is why windows must be allowed
    to sit entirely inside the novel tail rather than being required to
    contain the whole altered interval. Frameshifts were 0.8% of candidate
    variants in the Ott benchmark and ~17% of its validated positives.
    """

    def test_the_tail_replaces_everything_from_the_shift(self, proteome):
        wildtype = proteome.get(SYNTH_TRANSCRIPT)
        variant = synth_frameshift()
        mutant, start, end = build_mutant_protein(variant, wildtype)

        assert start == variant.protein_start - 1
        assert end == len(mutant)
        assert mutant[:start] == wildtype[:start]
        assert mutant[start:] == TAIL

    @pytest.mark.parametrize("terminator", ["*", "X"])
    def test_the_tail_is_truncated_at_the_first_terminator(self, proteome, terminator):
        """Past a stop there is no protein; past an `X` there is no known one."""
        variant = synth_frameshift(downstream_protein=f"{TAIL}{terminator}GGKKWW")
        mutant, _, _ = build_mutant_protein(variant, proteome.get(SYNTH_TRANSCRIPT))
        assert mutant.endswith(TAIL)
        assert terminator not in mutant

    def test_a_tail_identical_to_the_reference_is_rejected(self, proteome):
        """If the annotated tail matches the reference it encodes no shift."""
        wildtype = proteome.get(SYNTH_TRANSCRIPT)
        variant = synth_frameshift(downstream_protein=wildtype[40:60])
        with pytest.raises(ReferenceMismatch, match="encodes no frameshift"):
            build_mutant_protein(variant, wildtype)

    def test_a_shift_beyond_the_transcript_is_rejected(self, proteome):
        variant = synth_frameshift(protein_start=9999, protein_end=9999)
        with pytest.raises(ReferenceMismatch, match="outside transcript"):
            build_mutant_protein(variant, proteome.get(SYNTH_TRANSCRIPT))

    def test_every_peptide_carries_novel_sequence(self, proteome):
        wildtype = proteome.get(SYNTH_TRANSCRIPT)
        candidates = generate_for_variant(synth_frameshift(), proteome, lengths=[9])

        assert candidates
        for candidate in candidates:
            assert candidate.mutant_peptide not in wildtype

    def test_peptides_wholly_inside_the_tail_are_emitted(self, proteome):
        """The neo-ORF epitopes. A containment rule would drop all of these."""
        candidates = generate_for_variant(synth_frameshift(), proteome, lengths=[9])
        wholly_novel = [c for c in candidates if c.mutant_peptide in TAIL]

        assert len(wholly_novel) >= len(TAIL) - 8, (
            "every 9-mer window inside the novel tail should be a candidate"
        )

    def test_the_peptide_count_scales_with_the_tail_not_the_peptide_length(self, proteome):
        """The poly-epitope property, and what the per-variant cap must not erase.

        A substitution yields exactly one window per register, so its count is
        the peptide length and nothing else. A frameshift's count tracks how
        much novel sequence there is, so doubling the tail roughly doubles the
        candidates.
        """
        missense = generate_for_variant(kras_g12d(), proteome, lengths=[9])
        short = generate_for_variant(synth_frameshift(), proteome, lengths=[9])
        longer = generate_for_variant(
            synth_frameshift(downstream_protein=LONGER_TAIL), proteome, lengths=[9]
        )

        assert len(missense) == 9, "a substitution yields one window per register"
        assert len(short) > len(missense)
        assert len(longer) >= 2 * len(short) - 9

    def test_no_peptide_gets_a_wildtype_counterpart(self, proteome):
        """Nothing downstream of a shift has a positional counterpart."""
        for candidate in generate_for_variant(synth_frameshift(), proteome, lengths=[9, 10]):
            assert candidate.wildtype_peptide is None

    def test_a_same_length_frameshift_still_gets_no_counterpart(self, proteome):
        """Equal protein lengths must not be mistaken for a shared register."""
        wildtype = proteome.get(SYNTH_TRANSCRIPT)
        shift_at = 41
        tail = "".join("WYFMK"[i % 5] for i in range(len(wildtype) - shift_at + 1))
        variant = synth_frameshift(downstream_protein=tail)

        mutant, _, _ = build_mutant_protein(variant, wildtype)
        assert len(mutant) == len(wildtype), "fixture must produce an equal-length protein"
        for candidate in generate_for_variant(variant, proteome, lengths=[9]):
            assert candidate.wildtype_peptide is None

    def test_mutation_position_stays_inside_the_peptide(self, proteome):
        """A peptide downstream of the shift is entirely novel, so position 1."""
        candidates = generate_for_variant(synth_frameshift(), proteome, lengths=[9])
        for candidate in candidates:
            assert 1 <= candidate.mutation_position <= candidate.length
        assert any(c.mutation_position == 1 for c in candidates)

    def test_flanks_come_from_the_mutant_protein(self, proteome):
        candidates = generate_for_variant(synth_frameshift(), proteome, lengths=[9])
        inside = next(c for c in candidates if c.mutant_peptide in TAIL)
        assert inside.n_flank or inside.c_flank
        combined = f"{inside.n_flank}{inside.mutant_peptide}{inside.c_flank}"
        mutant, _, _ = build_mutant_protein(synth_frameshift(), proteome.get(SYNTH_TRANSCRIPT))
        assert combined in mutant


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
