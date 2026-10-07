import pytest

from neoantigene.io.vcf import VCFError, read_variant_tsv, read_vep_vcf
from neoantigene.models import VariantClass


@pytest.fixture()
def variants(examples_dir):
    return list(read_vep_vcf(examples_dir / "somatic.vep.vcf", tumor_sample="TUMOR"))


def test_only_protein_altering_consequences_are_returned(variants):
    assert len(variants) == 6
    assert all(v.variant_class is not VariantClass.OTHER for v in variants)
    assert not any(v.protein_start == 60 for v in variants)


def test_missense_fields_are_parsed(variants):
    g12d = next(v for v in variants if v.aa_alt == "D")
    assert (g12d.gene, g12d.transcript) == ("KRAS", "ENST00000256078")
    assert g12d.variant_class is VariantClass.MISSENSE
    assert (g12d.protein_start, g12d.protein_end) == (12, 12)
    assert (g12d.aa_ref, g12d.aa_alt) == ("G", "D")
    assert g12d.hgvsp_short == "KRAS p.G12D"


def test_vaf_is_read_from_the_named_tumor_sample(variants):
    g12d = next(v for v in variants if v.aa_alt == "D")
    assert g12d.dna_vaf == pytest.approx(54 / 142)
    assert g12d.tumor_depth == 142


def test_defaulting_to_the_last_sample_matches_explicit_selection(examples_dir):
    implicit = list(read_vep_vcf(examples_dir / "somatic.vep.vcf"))
    explicit = list(read_vep_vcf(examples_dir / "somatic.vep.vcf", tumor_sample="TUMOR"))
    assert [v.dna_vaf for v in implicit] == [v.dna_vaf for v in explicit]


def test_selecting_the_normal_sample_gives_zero_vaf(examples_dir):
    normal = list(read_vep_vcf(examples_dir / "somatic.vep.vcf", tumor_sample="NORMAL"))
    g12d = next(v for v in normal if v.aa_alt == "D")
    assert g12d.dna_vaf == 0.0


def test_unknown_sample_name_is_an_error(examples_dir):
    with pytest.raises(VCFError):
        list(read_vep_vcf(examples_dir / "somatic.vep.vcf", tumor_sample="MISSING"))


def test_inframe_deletion_positions_and_alleles(variants):
    deletion = next(v for v in variants if v.variant_class is VariantClass.INFRAME_DEL)
    assert (deletion.protein_start, deletion.protein_end) == (41, 43)
    assert (deletion.aa_ref, deletion.aa_alt) == ("LAV", "")
    assert deletion.hgvsp_short == "SYNTHA p.LAV41_43del"


def test_filters_and_population_frequency_are_carried(variants):
    weak = next(v for v in variants if v.protein_start == 52)
    assert weak.filters == ("weak_evidence",)
    common = next(v for v in variants if v.protein_start == 58)
    assert common.population_af == pytest.approx(0.0231)


def test_frameshift_sequence_wins_when_downstream_omits_the_altered_residue(tmp_path):
    """The whole-protein field includes the residue DownstreamProtein can drop."""
    path = tmp_path / "frameshift.vcf"
    path.write_text(
        "##fileformat=VCFv4.2\n"
        '##VEP="v114.2" ensembl-io=114 ensembl=114\n'
        '##INFO=<ID=CSQ,Number=.,Type=String,Description="Consequence annotations '
        "from Ensembl VEP. Format: Allele|Consequence|SYMBOL|Feature|"
        "Protein_position|Amino_acids|DownstreamProtein|ProteinLengthChange|"
        'FrameshiftSequence">\n'
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
        "1\t100\t.\tCT\tC\t.\tPASS\t"
        "CSQ=C|frameshift_variant|SYNTHA|ENST90000000001|41|L/-|WYFMKWYFH|4|"
        "MMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMNWYFMKWYFH\n"
    )
    variant = next(read_vep_vcf(path))
    assert variant.vep_release == 114
    assert variant.downstream_protein == "NWYFMKWYFH"
    assert variant.protein_length_change is None
    assert variant.variant_class is VariantClass.FRAMESHIFT


def test_downstream_protein_keeps_its_length_change(tmp_path):
    path = tmp_path / "downstream-only.vcf"
    path.write_text(
        "##fileformat=VCFv4.2\n"
        '##VEP="v114.2" ensembl=114\n'
        '##INFO=<ID=CSQ,Number=.,Type=String,Description="Consequence annotations '
        "from Ensembl VEP. Format: Allele|Consequence|SYMBOL|Feature|"
        'Protein_position|Amino_acids|DownstreamProtein|ProteinLengthChange">\n'
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
        "1\t100\t.\tCT\tC\t.\tPASS\t"
        "CSQ=C|frameshift_variant|SYNTHA|ENST90000000001|41|L/N|WYFMKWYFH|4\n"
    )
    variant = next(read_vep_vcf(path))
    assert variant.downstream_protein == "WYFMKWYFH"
    assert variant.protein_length_change == 4
    assert variant.aa_alt == "N"


def test_vep_version_is_the_release_when_the_cache_line_is_absent(tmp_path):
    path = tmp_path / "version-only.vcf"
    path.write_text(
        "##fileformat=VCFv4.2\n"
        '##VEP="v111"\n'
        '##INFO=<ID=CSQ,Number=.,Type=String,Description="Consequence annotations '
        "from Ensembl VEP. Format: Allele|Consequence|SYMBOL|Feature|"
        'Protein_position|Amino_acids">\n'
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
        "12\t25245350\t.\tC\tT\t.\tPASS\t"
        "CSQ=T|missense_variant|KRAS|ENST00000256078|12|G/D\n"
    )
    variant = next(read_vep_vcf(path))
    assert variant.vep_release == 111
    assert variant.protein_length_change is None


def test_a_frameshift_sequence_is_sliced_and_does_not_keep_the_length_change(tmp_path):
    path = tmp_path / "pvac.vcf"
    path.write_text(
        "##fileformat=VCFv4.2\n"
        '##VEP="v115" ensembl=115\n'
        '##INFO=<ID=CSQ,Number=.,Type=String,Description="Consequence annotations '
        "from Ensembl VEP. Format: Allele|Consequence|SYMBOL|Feature|"
        'Protein_position|Amino_acids|FrameshiftSequence|ProteinLengthChange">\n'
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
        "1\t100\t.\tCT\tC\t.\tPASS\t"
        "CSQ=C|frameshift_variant|SYNTHA|ENST90000000001|41|L/-|"
        "MMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMWYFMKWYFH|4\n"
    )
    variant = next(read_vep_vcf(path))
    assert variant.downstream_protein == "WYFMKWYFH"
    assert variant.protein_length_change is None
    assert variant.vep_release == 115


def test_vcf_without_csq_header_is_rejected(tmp_path):
    path = tmp_path / "plain.vcf"
    path.write_text(
        "##fileformat=VCFv4.2\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
        "1\t100\t.\tA\tT\t.\tPASS\t.\n"
    )
    with pytest.raises(VCFError, match="CSQ"):
        list(read_vep_vcf(path))


def test_vcf_and_tsv_paths_agree(examples_dir):
    from_vcf = list(read_vep_vcf(examples_dir / "somatic.vep.vcf", tumor_sample="TUMOR"))
    from_tsv = list(read_variant_tsv(examples_dir / "variants.tsv"))

    def identity(variant):
        return (variant.gene, variant.protein_start, variant.aa_ref, variant.aa_alt)

    tsv_identities = {identity(v) for v in from_tsv}
    assert tsv_identities <= {identity(v) for v in from_vcf}
