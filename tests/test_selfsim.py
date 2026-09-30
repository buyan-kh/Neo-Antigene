"""Nearest-self search.

The property that matters is the seeding guarantee: for any single substitution
in a peptide of length >= 8, at least one 4-mer avoids the mutated position, so
the wild-type counterpart is always found. If that breaks, the
`self_dissimilarity` feature silently starts calling ordinary self-adjacent
peptides foreign, which is the flattering direction and therefore the dangerous
one.
"""

from pathlib import Path

import pytest

from neoantigene.matrices import BLOSUM62
from neoantigene.peptides import selfsim
from neoantigene.peptides.proteome import ProteomeIndex
from neoantigene.peptides.reference import MIN_REAL_PROTEOME_PROTEINS, provenance
from neoantigene.peptides.selfsim import SEED_LENGTH, SelfProteome

EXAMPLES = Path(__file__).resolve().parents[1] / "data" / "examples"

# One invented protein, used only as a self background. No claim is made about
# it being a real sequence; the tests assert relationships, not biology.
BACKGROUND = "MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQANLVSTQLLLNGDNADILIGLLQKMTEYKLVVVGAGGVGKSALTIQ"


@pytest.fixture(scope="module")
def reference() -> SelfProteome:
    return SelfProteome([BACKGROUND])


def substitute(peptide: str, position: int) -> str:
    """Change one residue, 1-based, to something that is definitely different."""
    replacement = "W" if peptide[position - 1] != "W" else "G"
    return peptide[: position - 1] + replacement + peptide[position:]


class TestNearestSimilarity:
    def test_an_exact_self_peptide_scores_one(self, reference):
        assert reference.nearest_similarity([BACKGROUND[5:14]])[BACKGROUND[5:14]] == 1.0

    def test_an_unrelated_peptide_scores_below_a_self_adjacent_one(self, reference):
        native = BACKGROUND[20:29]
        mutant = substitute(native, 5)
        scores = reference.nearest_similarity([native, mutant, "WWWWWWWWW"])
        assert scores[native] > scores[mutant] > scores["WWWWWWWWW"]

    @pytest.mark.parametrize("length", [8, 9, 10, 11])
    def test_every_single_substitution_still_finds_its_wildtype(self, reference, length):
        """The seeding guarantee, checked at every position of every length."""
        native = BACKGROUND[30 : 30 + length]
        for position in range(1, length + 1):
            mutant = substitute(native, position)
            score = reference.nearest_similarity([mutant])[mutant]
            assert score > 0.0, (
                f"no seed hit for a single substitution at P{position} of a {length}mer, "
                f"so the wild-type counterpart was missed"
            )

    def test_similarity_is_bounded(self, reference):
        peptides = [BACKGROUND[i : i + 9] for i in range(0, 40, 3)]
        peptides += [substitute(p, 4) for p in peptides]
        for value in reference.nearest_similarity(peptides).values():
            assert 0.0 <= value <= 1.0

    def test_a_conservative_substitution_stays_closer_to_self(self, reference):
        """BLOSUM62, not identity: I->V must read as more self than I->P."""
        native = BACKGROUND[25:34]
        position = native.index("I") + 1
        conservative = native[: position - 1] + "V" + native[position:]
        radical = native[: position - 1] + "P" + native[position:]
        assert BLOSUM62[("I", "V")] > BLOSUM62[("I", "P")]
        scores = reference.nearest_similarity([conservative, radical])
        assert scores[conservative] > scores[radical]

    def test_no_match_can_span_two_proteins(self):
        """A window straddling a boundary would invent a self peptide."""
        left, right = "MKTAYIAKQ", "RQISFVKSH"
        reference = SelfProteome([left, right])
        straddling = left[-4:] + right[:5]
        assert reference.nearest_similarity([straddling])[straddling] < 1.0

    def test_empty_and_degenerate_inputs(self, reference):
        assert reference.nearest_similarity([]) == {}
        assert reference.nearest_similarity([""]) == {}
        assert SelfProteome([]).nearest_similarity(["AGGVGKSAL"]) == {"AGGVGKSAL": 0.0}

    def test_seed_length_guarantees_one_clean_seed(self):
        """Why SEED_LENGTH is 4 and not 5, as an executable argument."""
        for length in (8, 9, 10, 11):
            for position in range(1, length + 1):
                starts = range(1, length - SEED_LENGTH + 2)
                clean = [s for s in starts if not s <= position <= s + SEED_LENGTH - 1]
                assert clean, f"length {length}, mutation at P{position} has no clean seed"


class TestExactMatches:
    """The self-peptide hard gate, batched.

    `exact_matches` replaced a per-candidate substring scan over the whole
    reference. Batching is only worth doing if it answers identically, so the
    equivalence test below is the one that matters; the rest pin the edges
    where a prefilter-and-verify scheme could plausibly go wrong.
    """

    @pytest.mark.parametrize("length", [8, 9, 10, 11])
    def test_a_real_self_peptide_is_found_at_every_offset(self, reference, length):
        for start in range(len(BACKGROUND) - length + 1):
            peptide = BACKGROUND[start : start + length]
            assert reference.exact_matches([peptide]) == {peptide}

    def test_a_single_substitution_is_not_a_self_match(self, reference):
        native = BACKGROUND[20:29]
        assert reference.exact_matches([substitute(native, 5)]) == set()

    def test_it_agrees_with_the_per_candidate_scan_it_replaced(self, reference):
        """The equivalence that licenses the optimization."""
        index = ProteomeIndex()
        index._by_protein = {"BACKGROUND": BACKGROUND}

        queries: list[str] = []
        for length in (8, 9, 10, 11):
            for start in range(0, len(BACKGROUND) - length + 1, 3):
                native = BACKGROUND[start : start + length]
                queries += [native, substitute(native, 1), substitute(native, length)]
        queries += ["WWWWWWWWW", "GGGGGGGG", BACKGROUND[:11]]

        expected = {q for q in queries if index.contains_peptide(q)}
        assert reference.exact_matches(queries) == expected
        assert expected, "the fixture must contain real self matches for this to mean anything"

    def test_no_match_can_span_two_proteins(self):
        left, right = "MKTAYIAKQ", "RQISFVKSH"
        proteome = SelfProteome([left, right])
        straddling = left[-4:] + right[:5]
        assert proteome.exact_matches([straddling]) == set()
        assert proteome.exact_matches([left]) == {left}

    def test_peptides_spanning_a_chunk_boundary_are_still_found(self, monkeypatch):
        """The proteome is streamed, so every boundary offset must be covered.

        Forced to a tiny chunk so the boundaries fall inside the fixture. A
        peptide whose prefix starts just before a chunk edge is the case an
        off-by-one in the stride would silently drop.
        """
        monkeypatch.setattr(selfsim, "_CHUNK_RESIDUES", 4)
        proteome = SelfProteome([BACKGROUND])
        for length in (8, 9, 10, 11):
            peptides = [
                BACKGROUND[start : start + length] for start in range(len(BACKGROUND) - length + 1)
            ]
            assert proteome.exact_matches(peptides) == set(peptides)

    def test_non_standard_residues_never_match_the_separator(self):
        """`X` encodes like the inter-protein separator, so it must be rejected."""
        proteome = SelfProteome(["MKTAYIAKQ", "RQISFVKSH"])
        assert proteome.exact_matches(["MKTAYIAKX", "XXXXXXXX"]) == set()

    def test_empty_and_degenerate_inputs(self, reference):
        assert reference.exact_matches([]) == set()
        assert reference.exact_matches([""]) == set()
        assert SelfProteome([]).exact_matches(["AGGVGKSAL"]) == set()

    def test_duplicate_queries_are_answered_once(self, reference):
        peptide = BACKGROUND[10:19]
        assert reference.exact_matches([peptide] * 5) == {peptide}


class TestStubDetection:
    def test_the_bundled_example_is_reported_as_a_stub(self):
        report = provenance(EXAMPLES / "proteome.mini.fa")
        assert report.is_stub
        assert report.proteins < MIN_REAL_PROTEOME_PROTEINS
        assert report.caveat is not None
        assert "fetch-proteome" in report.caveat

    def test_the_index_agrees_with_the_standalone_check(self):
        index = ProteomeIndex.from_fasta(EXAMPLES / "proteome.mini.fa")
        assert index.is_stub
        assert index.protein_count == provenance(EXAMPLES / "proteome.mini.fa").proteins
