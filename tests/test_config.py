import pytest
import yaml
from pydantic import ValidationError

from neoantigene.config import PipelineConfig, ScoringWeights, discover, read_config_file


class TestValidation:
    def test_defaults_are_valid(self):
        weights = PipelineConfig().weights
        # Well-evidenced features must outweigh the speculative ones.
        assert weights.presentation > weights.agretopicity
        assert weights.clonality > weights.wt_dissimilarity
        assert weights.expression > weights.hydrophobicity

    def test_unknown_keys_are_rejected(self):
        """A typo in a config file must fail loudly, not be silently ignored."""
        with pytest.raises(ValidationError):
            PipelineConfig.model_validate({"weights": {"presentaton": 1.0}})

    def test_out_of_range_values_are_rejected(self):
        with pytest.raises(ValidationError):
            PipelineConfig.model_validate({"variant_filters": {"tumor_purity": 0.0}})
        with pytest.raises(ValidationError):
            PipelineConfig.model_validate({"variant_filters": {"min_ccf": 1.4}})
        with pytest.raises(ValidationError):
            PipelineConfig.model_validate({"output": {"top_n": 0}})

    def test_assignment_is_validated(self):
        config = PipelineConfig()
        with pytest.raises(ValidationError):
            config.output.top_n = -5

    def test_peptide_lengths_are_normalized_on_demand(self):
        config = PipelineConfig.model_validate({"peptides": {"lengths": [10, 9, 9, 8]}})
        assert config.peptides.sorted_lengths == [8, 9, 10]


class TestLoading:
    def test_yaml_round_trip(self, tmp_path):
        original = PipelineConfig()
        original.weights = ScoringWeights(presentation=3.3)
        path = original.dump(tmp_path / "config.yaml")
        assert PipelineConfig.load(path).weights.presentation == 3.3

    def test_toml_is_supported(self, tmp_path):
        path = tmp_path / "neoantigene.toml"
        path.write_text("[weights]\npresentation = 4.0\n\n[output]\ntop_n = 7\n")
        config = PipelineConfig.load(path)
        assert config.weights.presentation == 4.0
        assert config.output.top_n == 7

    def test_empty_yaml_yields_defaults(self, tmp_path):
        path = tmp_path / "config.yaml"
        path.write_text("")
        assert PipelineConfig.load(path) == PipelineConfig()

    def test_non_mapping_yaml_is_rejected(self, tmp_path):
        path = tmp_path / "config.yaml"
        path.write_text("- a\n- b\n")
        with pytest.raises(ValueError, match="mapping"):
            read_config_file(path)

    def test_no_config_anywhere_yields_defaults(self, tmp_path):
        assert PipelineConfig.load(search_from=tmp_path) == PipelineConfig()

    def test_discovery_prefers_the_first_match(self, tmp_path):
        (tmp_path / "config").mkdir()
        (tmp_path / "config" / "default.yaml").write_text("{}")
        assert discover(tmp_path) == tmp_path / "config" / "default.yaml"

        (tmp_path / "neoantigene.yaml").write_text("{}")
        assert discover(tmp_path) == tmp_path / "neoantigene.yaml"

    def test_shipped_default_config_matches_code_defaults(self):
        """`config/default.yaml` documents the priors; drift would mislead."""
        from pathlib import Path

        path = Path(__file__).resolve().parents[1] / "config" / "default.yaml"
        loaded = PipelineConfig.model_validate(yaml.safe_load(path.read_text()))
        assert loaded == PipelineConfig()
