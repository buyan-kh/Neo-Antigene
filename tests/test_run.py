"""Reproducibility tests.

A shortlist that cannot be traced back to the inputs and config that produced
it is not something a lab can act on or a partner can audit.
"""

import json
import logging
import re

from neoantigene.config import PipelineConfig
from neoantigene.pipeline import run as run_pipeline
from neoantigene.presentation.null import NullBackend
from neoantigene.run import (
    RunManifest,
    config_digest,
    configure_logging,
    file_digest,
    new_run_id,
)


class TestRunId:
    def test_is_sortable_and_unique(self):
        ids = [new_run_id() for _ in range(5)]
        assert len(set(ids)) == 5
        assert all(re.fullmatch(r"\d{8}T\d{6}Z-[0-9a-f]{8}", i) for i in ids)


class TestDigests:
    def test_file_digest_is_stable_and_content_sensitive(self, tmp_path):
        path = tmp_path / "a.txt"
        path.write_text("hello")
        first = file_digest(path)
        assert first == file_digest(path)
        assert first.startswith("sha256:")

        path.write_text("hello!")
        assert file_digest(path) != first

    def test_config_digest_ignores_key_order(self):
        left = PipelineConfig.model_validate({"output": {"top_n": 5}, "weights": {"bias": -1.0}})
        right = PipelineConfig.model_validate({"weights": {"bias": -1.0}, "output": {"top_n": 5}})
        assert config_digest(left) == config_digest(right)

    def test_config_digest_changes_with_weights(self):
        base = PipelineConfig()
        changed = base.model_copy(
            update={"weights": base.weights.model_copy(update={"presentation": 9.9})}
        )
        assert config_digest(base) != config_digest(changed)


class TestManifest:
    def test_records_inputs_config_and_environment(self, example_sample, dev_config):
        manifest = RunManifest.start(
            sample_id=example_sample.sample_id,
            config=dev_config,
            backend="null",
            inputs=example_sample.input_paths(),
        )
        assert manifest.sample_id == example_sample.sample_id
        assert manifest.backend == "null"
        assert manifest.package_version
        assert manifest.python_version.startswith("3.")
        assert manifest.inputs
        assert all(i.digest.startswith("sha256:") and i.bytes > 0 for i in manifest.inputs)
        assert manifest.config == dev_config

    def test_finish_records_counts_without_mutating(self, example_sample, dev_config):
        manifest = RunManifest.start(
            example_sample.sample_id, dev_config, "null", example_sample.input_paths()
        )
        finished = manifest.finish({"shortlisted": 8})
        assert manifest.finished_at is None
        assert finished.finished_at is not None
        assert finished.counts == {"shortlisted": 8}

    def test_written_manifest_reloads(self, example_sample, dev_config, tmp_path):
        manifest = RunManifest.start(
            example_sample.sample_id, dev_config, "null", example_sample.input_paths()
        ).finish({"shortlisted": 3})
        path = manifest.write(tmp_path / "run.json")
        reloaded = RunManifest.model_validate_json(path.read_text())
        assert reloaded == manifest

    def test_pipeline_attaches_a_manifest_with_real_counts(self, example_sample, dev_config):
        result = run_pipeline(example_sample, dev_config, backend=NullBackend(), run_id="fixed-id")
        manifest = result.manifest
        assert manifest.run_id == "fixed-id"
        assert manifest.finished_at is not None
        assert manifest.counts["shortlisted"] == len(result.ranked)
        assert manifest.counts["variants_read"] == result.report.variants_read
        assert manifest.config == result.config

    def test_manifest_config_reflects_sample_purity_override(self, example_sample, dev_config):
        result = run_pipeline(example_sample, dev_config, backend=NullBackend())
        assert result.manifest.config.variant_filters.tumor_purity == example_sample.tumor_purity

    def test_same_inputs_give_the_same_config_digest(self, example_sample, dev_config):
        first = run_pipeline(example_sample, dev_config, backend=NullBackend())
        second = run_pipeline(example_sample, dev_config, backend=NullBackend())
        assert first.manifest.config_digest == second.manifest.config_digest
        assert first.manifest.run_id != second.manifest.run_id


class TestLogging:
    def test_records_carry_the_run_id(self, capsys):
        configure_logging("run-abc")
        logging.getLogger("neoantigene.test").info("hello")
        line = capsys.readouterr().err.strip().splitlines()[-1]
        assert "[run-abc]" in line
        assert "hello" in line

    def test_library_logs_do_not_leak_to_the_root_logger(self, capsys):
        """A library that reconfigures the host application's logging is rude."""
        configure_logging("run-xyz")
        assert logging.getLogger("neoantigene").propagate is False

    def test_json_logs_are_parseable(self, capsys):
        configure_logging("run-json", json_logs=True)
        logging.getLogger("neoantigene.test").info("structured")
        payload = json.loads(capsys.readouterr().err.strip().splitlines()[-1])
        assert payload["run_id"] == "run-json"
        assert payload["message"] == "structured"
        assert payload["level"] == "INFO"
