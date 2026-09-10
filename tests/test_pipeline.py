import pytest

from neoantigene.config import PipelineConfig
from neoantigene.io.manifest import Sample
from neoantigene.io.writers import to_records, write_tsv
from neoantigene.models import GateFailure
from neoantigene.pipeline import resolve_config
from neoantigene.pipeline import run as run_pipeline
from neoantigene.presentation.null import NullBackend


@pytest.fixture
def result(example_sample, dev_config):
    return run_pipeline(example_sample, dev_config, backend=NullBackend())


class TestEndToEnd:
    def test_produces_a_ranked_shortlist(self, result):
        assert result.report.variants_read == 5
        # The VAF 0.08 variant is subclonal at 60% purity and must be dropped.
        assert result.report.variants_kept == 4
        assert "subclonal" in result.report.filter_report.summary()
        assert result.report.peptides_generated > 0
        assert result.ranked
        scores = [item.score for item in result.ranked]
        assert scores == sorted(scores, reverse=True)

    def test_is_deterministic(self, example_sample, dev_config):
        first = run_pipeline(example_sample, dev_config, backend=NullBackend())
        second = run_pipeline(example_sample, dev_config, backend=NullBackend())
        assert [(i.key, i.score) for i in first.ranked] == [(i.key, i.score) for i in second.ranked]

    def test_every_hla_allele_is_evaluated(self, result, example_sample):
        assert {item.allele for item in result.all_scored} == set(example_sample.hla)

    def test_shortlist_is_a_subset_of_everything_scored(self, result):
        assert {i.key for i in result.ranked} <= {i.key for i in result.all_scored}

    def test_self_peptide_gate_keeps_wildtype_sequences_out(self, result):
        assert all(item.features["self_match"] == 0.0 for item in result.ranked)

    def test_gate_counts_are_reported(self, result):
        assert result.report.pairs_gated
        assert all(isinstance(k, GateFailure) for k in result.report.pairs_gated)


class TestPurity:
    def test_run_does_not_mutate_the_config_it_was_given(self, example_sample, dev_config):
        before = dev_config.model_dump()
        run_pipeline(example_sample, dev_config, backend=NullBackend())
        assert dev_config.model_dump() == before

    def test_manifest_purity_overrides_config(self, example_sample, dev_config):
        dev_config.variant_filters.tumor_purity = 0.1
        resolved = resolve_config(example_sample, dev_config)
        assert resolved.variant_filters.tumor_purity == example_sample.tumor_purity
        assert dev_config.variant_filters.tumor_purity == 0.1

    def test_config_without_sample_purity_is_returned_unchanged(self, dev_config, examples_dir):
        sample = Sample.load(examples_dir / "sample.yaml").model_copy(update={"tumor_purity": None})
        assert resolve_config(sample, dev_config) is dev_config


class TestVcfInput:
    def test_vcf_path_runs_end_to_end(self, examples_dir, dev_config, tmp_path):
        import yaml

        manifest = tmp_path / "sample.yaml"
        manifest.write_text(
            yaml.safe_dump(
                {
                    "sample_id": "VCF-001",
                    "tumor_purity": 0.6,
                    "hla": ["HLA-A*11:01", "HLA-C*08:02"],
                    "processed": {
                        "somatic_vcf": str(examples_dir / "somatic.vep.vcf"),
                        "expression_tsv": str(examples_dir / "expression.tsv"),
                        "proteome_fasta": str(examples_dir / "proteome.mini.fa"),
                        "tumor_sample_name": "TUMOR",
                    },
                }
            )
        )
        result = run_pipeline(Sample.load(manifest), dev_config, backend=NullBackend())
        assert result.report.variants_read == 6
        summary = result.report.filter_report.summary()
        assert "germline_common" in summary
        assert "not_pass" in summary
        assert result.ranked


class TestValidation:
    def test_missing_hla_is_rejected(self, example_sample, dev_config):
        sample = example_sample.model_copy(update={"hla": []})
        with pytest.raises(ValueError, match="no HLA alleles"):
            run_pipeline(sample, dev_config, backend=NullBackend())

    def test_missing_proteome_is_rejected(self, example_sample, dev_config):
        processed = example_sample.processed.model_copy(update={"proteome_fasta": None})
        sample = example_sample.model_copy(update={"processed": processed})
        with pytest.raises(ValueError, match="proteome_fasta"):
            run_pipeline(sample, dev_config, backend=NullBackend())


class TestOutput:
    def test_round_trips_to_tsv(self, result, tmp_path):
        path = tmp_path / "ranked.tsv"
        frame = write_tsv(result.ranked, result.sample.sample_id, path, run_id="r1")
        assert path.exists()
        assert frame.height == len(result.ranked)
        assert frame["rank"].to_list() == list(range(1, len(result.ranked) + 1))
        assert set(frame["run_id"].to_list()) == {"r1"}

    def test_records_carry_the_sample_and_peptide(self, result):
        records = to_records(result.ranked, result.sample.sample_id)
        assert all(r["sample_id"] == result.sample.sample_id for r in records)
        assert all(r["mutant_peptide"] for r in records)

    def test_nan_sentinels_are_written_as_empty(self, result, tmp_path):
        path = tmp_path / "ranked.tsv"
        write_tsv(result.ranked, result.sample.sample_id, path)
        assert "NaN" not in path.read_text()


class TestSyntheticCohort:
    def test_runs_on_every_generated_patient(self, cohort, benchmark_config):
        for patient in cohort.patients:
            sample = Sample.load(patient.manifest_path)
            result = run_pipeline(sample, benchmark_config, backend=NullBackend())
            assert result.report.variants_read == patient.variant_count
            assert result.all_scored

    def test_disabled_gates_leave_candidates_ungated(self, benchmark_run):
        assert all(item.passed for item in benchmark_run.all_scored)


def test_empty_candidate_set_returns_an_empty_result(examples_dir, tmp_path):
    import yaml

    empty = tmp_path / "variants.tsv"
    empty.write_text(
        "chrom\tpos\tref\talt\tgene\ttranscript\tconsequence\tprotein_position\tamino_acids\n"
    )
    manifest = tmp_path / "sample.yaml"
    manifest.write_text(
        yaml.safe_dump(
            {
                "sample_id": "EMPTY",
                "hla": ["HLA-A*02:01"],
                "processed": {
                    "variant_tsv": "variants.tsv",
                    "proteome_fasta": str(examples_dir / "proteome.mini.fa"),
                },
            }
        )
    )
    config = PipelineConfig()
    config.presentation.backend = "null"
    result = run_pipeline(Sample.load(manifest), config, backend=NullBackend())
    assert result.ranked == []
    assert result.all_scored == []
    assert result.manifest.finished_at is not None
