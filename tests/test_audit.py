"""Reconciliation of published epitope claims.

This is an oracle, so the tests are about whether it ever gives a confident
wrong answer. Three failure modes would each manufacture a finding:

- calling a claim contradicted when it merely could not be checked
- reconciling against the wrong isoform because the first one was tried
- losing a gene to a rename and reporting it as absent

Each has a test below. The fixture proteome is invented; the assertions are
about arithmetic, not biology, except where real published KRAS epitopes are
used as a positive control.
"""

import json

import pytest

from neoantigene.audit import (
    EpitopeClaim,
    SymbolIndex,
    Verdict,
    parse_protein_change,
    reconcile,
    reconcile_claim,
    write_report,
)
from neoantigene.models import VariantClass

#: Real KRAS 4B protein sequence prefix, so the published G12D epitopes are a
#: genuine positive control rather than a construction.
KRAS = (
    "MTEYKLVVVGAGGVGKSALTIQLIQNHFVDEYDPTIEDSYRKQVVIDGETCLLDILDTAGQEEYSAMRDQYMRTGEG"
    "FLCVFAINNTKSFEDIHHYREQIKRVKDSEDVPMVLVGNKCDLPSRTVDTKQAQDLARSYGIPFIETSAKTRQGVDDAF"
)
OTHER = "MKWVTFISLLLLFSSAYSRGVFRRDTHKSEIAHRFKDLGEENFKALVLIAFAQYLQQCPFD"


@pytest.fixture
def index(tmp_path):
    """Two KRAS isoforms plus a gene reachable only through a rename.

    The short isoform is a truncation that still carries G at 12, so it is a
    valid but worse answer than the full one — which is what makes isoform
    ordering testable.
    """
    fasta = tmp_path / "pep.all.fa"
    fasta.write_text(
        f">ENSP1 pep chromosome:GRCh38:12:1:2:1 gene:ENSG1 "
        f"transcript:ENST_SHORT.1 gene_symbol:KRAS\n{KRAS[:40]}\n"
        f">ENSP2 pep chromosome:GRCh38:12:1:2:1 gene:ENSG1 "
        f"transcript:ENST_FULL.1 gene_symbol:KRAS\n{KRAS}\n"
        f">ENSP3 pep chromosome:GRCh38:3:1:2:1 gene:ENSG2 "
        f"transcript:ENST_ACP3.1 gene_symbol:ACP3\n{OTHER}\n"
    )
    return SymbolIndex.from_fasta(fasta)


def claim(peptide: str, gene: str = "KRAS", change: str = "p.G12D", **kwargs) -> EpitopeClaim:
    return EpitopeClaim(peptide=peptide, gene=gene, protein_change=change, **kwargs)


class TestProteinChange:
    @pytest.mark.parametrize("text", ["p.G12D", "G12D", "p.G12D "])
    def test_missense_is_read_with_or_without_the_prefix(self, text):
        change = parse_protein_change(text)
        assert change is not None
        assert (change.start, change.aa_ref, change.aa_alt) == (12, "G", "D")
        assert change.variant_class is VariantClass.MISSENSE

    @pytest.mark.parametrize("text", ["p.S754fs", "S754fs", "p.S754Qfs*12"])
    def test_frameshift_is_matched_before_missense(self, text):
        """`p.S754Qfs*12` would otherwise fall through to no match at all."""
        change = parse_protein_change(text)
        assert change is not None
        assert change.variant_class is VariantClass.FRAMESHIFT
        assert change.start == 754

    def test_deletions_and_insertions_are_read(self):
        single = parse_protein_change("p.L41del")
        spanning = parse_protein_change("p.L41_V43del")
        inserted = parse_protein_change("p.L41_V42insAC")

        assert single is not None and single.variant_class is VariantClass.INFRAME_DEL
        assert spanning is not None and (spanning.start, spanning.end) == (41, 43)
        assert inserted is not None and inserted.aa_alt == "AC"

    @pytest.mark.parametrize("text", ["", "garbage", "p.12", "p.G12", "p.L43_V41del"])
    def test_unreadable_changes_return_none(self, text):
        assert parse_protein_change(text) is None

    def test_an_absurd_deletion_span_is_rejected(self):
        assert parse_protein_change("p.A1_B9999del") is None


class TestReconciliation:
    def test_published_kras_epitopes_reconcile(self, index):
        """Positive control: epitopes from the adoptive-transfer literature."""
        for peptide in ("GADGVGKSA", "GADGVGKSAL", "VVGADGVGK"):
            result = reconcile_claim(claim(peptide), index)
            assert result.verdict is Verdict.DERIVABLE, f"{peptide}: {result.detail}"

    def test_it_resolves_to_an_isoform_long_enough_to_contain_the_epitope(self, index):
        """A truncation also carrying G12 must not be preferred over the full one."""
        result = reconcile_claim(claim("GADGVGKSAL"), index)
        assert result.transcript == "ENST_FULL"

    def test_a_peptide_the_variant_cannot_produce_is_contradicted(self, index):
        result = reconcile_claim(claim("WWWWWWWWW"), index)
        assert result.verdict is Verdict.PEPTIDE_NOT_GENERATED
        assert not result.verdict.is_unverifiable

    def test_a_wrong_reference_residue_is_reported_as_isoform_inconsistency(self, index):
        """The CASP1 p.P172S shape: no isoform carries that residue there."""
        result = reconcile_claim(claim("GADGVGKSA", change="p.P12D"), index)
        assert result.verdict is Verdict.NO_CONSISTENT_ISOFORM
        assert "carries 'P' at 12" in result.detail

    def test_a_position_past_the_protein_is_not_a_reconciliation(self, index):
        result = reconcile_claim(claim("GADGVGKSA", change="p.G9999D"), index)
        assert result.verdict is Verdict.NO_CONSISTENT_ISOFORM

    def test_an_absent_gene_is_distinguished_from_a_bad_coordinate(self, index):
        result = reconcile_claim(claim("GADGVGKSA", gene="NOTAGENE"), index)
        assert result.verdict is Verdict.GENE_NOT_FOUND

    def test_a_frameshift_without_its_tail_is_unverifiable_not_contradicted(self, index):
        """The distinction that stops this oracle manufacturing findings."""
        result = reconcile_claim(claim("ILMHGLVSL", change="p.S754fs"), index)

        assert result.verdict is Verdict.TAIL_REQUIRED
        assert result.verdict.is_unverifiable
        assert not result.verdict.is_verified

    def test_a_frameshift_with_its_tail_is_checked(self, index):
        tail = "WYFMKWYFHKWYFQKWYFRK"
        supplied = claim("WYFMKWYFH", change="p.G12fs", downstream_protein=tail)
        assert reconcile_claim(supplied, index).verdict is Verdict.DERIVABLE

        wrong = claim("AAAAAAAAA", change="p.G12fs", downstream_protein=tail)
        assert reconcile_claim(wrong, index).verdict is Verdict.PEPTIDE_NOT_GENERATED

    def test_an_unreadable_change_is_unverifiable(self, index):
        result = reconcile_claim(claim("GADGVGKSA", change="not a change"), index)
        assert result.verdict is Verdict.UNPARSEABLE_CHANGE
        assert result.verdict.is_unverifiable

    def test_peptides_of_unusual_length_are_still_checked(self, index):
        """A 12-mer is outside the default lengths but was still claimed."""
        result = reconcile_claim(claim("VVVGADGVGKSA"), index)
        assert result.verdict is Verdict.DERIVABLE


class TestSymbolRenames:
    def test_a_retired_symbol_is_missing_until_the_rename_is_loaded(self, index, tmp_path):
        """On the Ott benchmark this step recovered 284 variants."""
        retired = claim("KWVTFISLL", gene="ACPP", change="p.K2W")
        assert reconcile_claim(retired, index).verdict is Verdict.GENE_NOT_FOUND

        aliases = tmp_path / "aliases.json"
        aliases.write_text(json.dumps({"ACPP": "ACP3"}))
        assert index.load_aliases(aliases) == 1

        assert reconcile_claim(retired, index).transcript == "ENST_ACP3"

    def test_a_rename_cannot_invent_a_lookup_target(self, index, tmp_path):
        """Only renames pointing at a gene actually in the FASTA are kept."""
        aliases = tmp_path / "aliases.json"
        aliases.write_text(json.dumps({"ACPP": "NOT_IN_FASTA", "OLD": "ACP3"}))
        assert index.load_aliases(aliases) == 1

    def test_a_rename_never_overrides_a_present_symbol(self, index, tmp_path):
        aliases = tmp_path / "aliases.json"
        aliases.write_text(json.dumps({"KRAS": "ACP3"}))
        assert index.load_aliases(aliases) == 0
        assert reconcile_claim(claim("GADGVGKSA"), index).transcript == "ENST_FULL"


class TestReport:
    def test_unverifiable_is_counted_apart_from_contradicted(self, index):
        report = reconcile(
            [
                claim("GADGVGKSA"),
                claim("WWWWWWWWW"),
                claim("ILMHGLVSL", change="p.S754fs"),
                claim("GADGVGKSA", change="garbage"),
            ],
            index,
        )
        body = report.describe()

        assert "derivable          1" in body
        assert "contradicted       1" in body
        assert "unverifiable       2" in body
        assert "Unverifiable is not contradicted" in body

    def test_an_empty_run_says_so_rather_than_dividing_by_zero(self, index):
        assert reconcile([], index).describe() == "no claims checked"

    def test_verdicts_round_trip_to_tsv(self, index, tmp_path):
        report = reconcile([claim("GADGVGKSA"), claim("WWWWWWWWW")], index)
        path = write_report(report, tmp_path / "verdicts.tsv")
        lines = path.read_text().strip().splitlines()

        assert lines[0].split("\t")[0] == "peptide"
        assert len(lines) == 3
        assert "derivable" in lines[1] or "derivable" in lines[2]
