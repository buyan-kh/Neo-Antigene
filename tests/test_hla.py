import pytest

from neoantigene.hla import (
    HLAError,
    class_i_only,
    gene_of,
    normalize_allele,
    parse_hla_string,
    read_hla_file,
)


@pytest.mark.parametrize(
    "raw,expected",
    [
        ("A0201", "HLA-A*02:01"),
        ("A*02:01", "HLA-A*02:01"),
        ("HLA-A02:01", "HLA-A*02:01"),
        ("HLA-A*02:01:01:02", "HLA-A*02:01"),
        ("  c0802 ", "HLA-C*08:02"),
        ("B*57:01", "HLA-B*57:01"),
        ("DRB1*04:01", "HLA-DRB1*04:01"),
    ],
)
def test_normalize_allele(raw, expected):
    assert normalize_allele(raw) == expected


@pytest.mark.parametrize("raw", ["not-an-allele", "", "A*02", "banana"])
def test_normalize_allele_rejects_garbage(raw):
    with pytest.raises(HLAError):
        normalize_allele(raw)


def test_parse_hla_string_dedupes_and_preserves_order():
    assert parse_hla_string("A0201, A*02:01 B0702") == ["HLA-A*02:01", "HLA-B*07:02"]


def test_gene_of():
    assert gene_of("HLA-A*02:01") == "A"
    assert gene_of("HLA-DRB1*04:01") == "DRB1"


def test_class_i_only_drops_class_ii():
    alleles = ["HLA-A*02:01", "HLA-DRB1*04:01", "HLA-C*08:02"]
    assert class_i_only(alleles) == ["HLA-A*02:01", "HLA-C*08:02"]


def test_read_hla_file_plain_list(tmp_path):
    path = tmp_path / "hla.txt"
    path.write_text("# genotype\nA0201\nB*07:02\n\nC08:02\n")
    assert read_hla_file(path) == ["HLA-A*02:01", "HLA-B*07:02", "HLA-C*08:02"]


def test_read_hla_file_tsv_with_allele_column(tmp_path):
    path = tmp_path / "hla.tsv"
    path.write_text("sample\tallele\nS1\tA0201\nS1\tA0201\nS1\tB0702\n")
    assert read_hla_file(path) == ["HLA-A*02:01", "HLA-B*07:02"]
