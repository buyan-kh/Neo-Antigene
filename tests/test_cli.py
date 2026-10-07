"""CLI smoke tests.

The CLI is a thin shell, so these assert wiring and guard rails rather than
ranking behaviour: that commands exit 0, write the files they claim to, and
refuse the dangerous option.
"""

import csv
import json

import pytest
from typer.testing import CliRunner

from neoantigene.assays.schema import AssayCall
from neoantigene.cli import app

from .fixtures import simulate_assay_results

runner = CliRunner()


@pytest.fixture
def ranked(tmp_path, examples_dir):
    out = tmp_path / "results"
    result = runner.invoke(
        app,
        [
            "rank",
            str(examples_dir / "sample.yaml"),
            "--backend",
            "null",
            "--allow-null-backend",
            "--out-dir",
            str(out),
        ],
    )
    assert result.exit_code == 0, result.output
    return out, result


class TestRank:
    def test_writes_all_three_artifacts(self, ranked):
        out, _ = ranked
        assert (out / "EXAMPLE-PDAC-001.ranked.tsv").exists()
        assert (out / "EXAMPLE-PDAC-001.features.json").exists()
        assert (out / "EXAMPLE-PDAC-001.run.json").exists()

    def test_reports_a_run_id_that_matches_the_manifest(self, ranked):
        out, result = ranked
        manifest = json.loads((out / "EXAMPLE-PDAC-001.run.json").read_text())
        assert manifest["run_id"] in result.output
        assert manifest["counts"]["shortlisted"] > 0

    def test_shortlist_rows_carry_the_run_id(self, ranked):
        out, _ = ranked
        with open(out / "EXAMPLE-PDAC-001.ranked.tsv") as handle:
            rows = list(csv.DictReader(handle, delimiter="\t"))
        assert rows
        assert len({r["run_id"] for r in rows}) == 1

    def test_null_backend_requires_acknowledgement(self, tmp_path, examples_dir):
        result = runner.invoke(
            app,
            [
                "rank",
                str(examples_dir / "sample.yaml"),
                "--backend",
                "null",
                "--out-dir",
                str(tmp_path),
            ],
        )
        assert result.exit_code != 0
        assert "development only" in result.output

    def test_missing_manifest_is_rejected(self, tmp_path):
        result = runner.invoke(app, ["rank", str(tmp_path / "nope.yaml")])
        assert result.exit_code != 0

    def test_active_selection_respects_batch_size(self, tmp_path, examples_dir):
        out = tmp_path / "active"
        result = runner.invoke(
            app,
            [
                "rank",
                str(examples_dir / "sample.yaml"),
                "--backend",
                "null",
                "--allow-null-backend",
                "--select",
                "active",
                "--top-n",
                "4",
                "--out-dir",
                str(out),
            ],
        )
        assert result.exit_code == 0, result.output
        with open(out / "EXAMPLE-PDAC-001.ranked.tsv") as handle:
            assert len(list(csv.DictReader(handle, delimiter="\t"))) == 4


class TestAssayLoop:
    def test_request_sheet_has_blank_call_column(self, ranked, tmp_path):
        out, _ = ranked
        request = tmp_path / "request.tsv"
        result = runner.invoke(
            app,
            ["request-assays", str(out / "EXAMPLE-PDAC-001.ranked.tsv"), "-o", str(request)],
        )
        assert result.exit_code == 0, result.output
        with open(request) as handle:
            rows = list(csv.DictReader(handle, delimiter="\t"))
        assert rows
        assert all(row["call"] == "" for row in rows)
        assert all(row["predicted_rank"] for row in rows)

    def test_evaluate_reports_a_validation_rate(self, ranked, tmp_path):
        out, _ = ranked
        completed = tmp_path / "completed.tsv"
        _write_completed_from_ranked(out / "EXAMPLE-PDAC-001.ranked.tsv", completed)

        result = runner.invoke(
            app,
            [
                "evaluate",
                str(completed),
                "--ranked",
                str(out / "EXAMPLE-PDAC-001.ranked.tsv"),
                "--k",
                "4",
            ],
        )
        assert result.exit_code == 0, result.output
        assert "neoantigene" in result.output
        assert "validated" in result.output

    def test_evaluate_stats_adds_inference_without_removing_the_point_estimate(
        self, ranked, tmp_path
    ):
        out, _ = ranked
        completed = tmp_path / "completed.tsv"
        _write_completed_from_ranked(out / "EXAMPLE-PDAC-001.ranked.tsv", completed)

        result = runner.invoke(
            app,
            [
                "evaluate",
                str(completed),
                "--ranked",
                str(out / "EXAMPLE-PDAC-001.ranked.tsv"),
                "--k",
                "4",
                "--stats",
            ],
        )
        assert result.exit_code == 0, result.output
        assert "validated" in result.output
        assert "chance" in result.output
        assert "p=" in result.output
        assert "auc" in result.output

    def test_compare_audits_arbitrary_rankings_and_writes_a_report(self, ranked, tmp_path):
        """The harness must work on files it did not produce."""
        out, _ = ranked
        completed = tmp_path / "completed.tsv"
        _write_completed_from_ranked(out / "EXAMPLE-PDAC-001.ranked.tsv", completed)

        ranking = out / "EXAMPLE-PDAC-001.ranked.tsv"
        report = tmp_path / "report.md"
        result = runner.invoke(
            app,
            [
                "compare",
                str(completed),
                "-m",
                f"ours={ranking}",
                "-m",
                f"theirs={ranking}",
                "--baseline",
                "theirs",
                "--k",
                "4",
                "--rounds",
                "200",
                "-o",
                str(report),
            ],
        )
        assert result.exit_code in (0, 2), result.output
        body = report.read_text()
        assert body.startswith("# Neoantigen ranking evaluation")
        assert "What this design could detect" in body
        assert "sha256:" in body
        assert "Reading this honestly" in body

    def test_compare_rejects_a_malformed_method_argument(self, ranked, tmp_path):
        out, _ = ranked
        completed = tmp_path / "completed.tsv"
        _write_completed_from_ranked(out / "EXAMPLE-PDAC-001.ranked.tsv", completed)

        result = runner.invoke(
            app, ["compare", str(completed), "-m", str(out / "EXAMPLE-PDAC-001.ranked.tsv")]
        )
        assert result.exit_code != 0
        assert "NAME=path" in result.output

    def test_compare_rejects_a_baseline_that_is_not_a_given_method(self, ranked, tmp_path):
        out, _ = ranked
        completed = tmp_path / "completed.tsv"
        _write_completed_from_ranked(out / "EXAMPLE-PDAC-001.ranked.tsv", completed)

        result = runner.invoke(
            app,
            [
                "compare",
                str(completed),
                "-m",
                f"ours={out / 'EXAMPLE-PDAC-001.ranked.tsv'}",
                "--baseline",
                "netmhcpan",
            ],
        )
        assert result.exit_code != 0
        assert "not one of the given methods" in result.output

    def test_refit_refuses_insufficient_data(self, ranked, tmp_path):
        out, _ = ranked
        completed = tmp_path / "completed.tsv"
        _write_completed_from_ranked(out / "EXAMPLE-PDAC-001.ranked.tsv", completed)

        result = runner.invoke(
            app,
            [
                "refit",
                str(completed),
                "-f",
                str(out / "EXAMPLE-PDAC-001.features.json"),
            ],
        )
        assert result.exit_code == 1
        assert "cannot refit" in result.output

    def test_refit_recovers_the_generative_model(self, tmp_path, cohort, benchmark_run):
        """End-to-end learning loop on synthetic truth."""
        from neoantigene.io.writers import write_features_json

        sample_id = benchmark_run.sample.sample_id
        features_path = write_features_json(
            benchmark_run.all_scored, sample_id, tmp_path / "features.json"
        )

        results = simulate_assay_results(benchmark_run.all_scored, sample_id, seed=3)
        completed = tmp_path / "completed.tsv"
        _write_results(results, completed)

        fitted = tmp_path / "fitted.yaml"
        result = runner.invoke(
            app,
            ["refit", str(completed), "-f", str(features_path), "--out", str(fitted)],
        )
        assert result.exit_code == 0, result.output
        assert fitted.exists()

        import yaml

        weights = yaml.safe_load(fitted.read_text())["weights"]
        # Truth was built from presentation, clonality and expression, so those
        # three must dominate the refit.
        driven = {"presentation", "clonality", "expression"}
        strongest = sorted((k for k in weights if k != "bias"), key=lambda k: -abs(weights[k]))[:3]
        assert driven == set(strongest), f"refit emphasised {strongest}"


class TestInspection:
    def test_explain_decomposes_a_score(self, ranked):
        out, _ = ranked
        with open(out / "EXAMPLE-PDAC-001.ranked.tsv") as handle:
            peptide = next(csv.DictReader(handle, delimiter="\t"))["mutant_peptide"]

        result = runner.invoke(
            app, ["explain", str(out / "EXAMPLE-PDAC-001.features.json"), peptide]
        )
        assert result.exit_code == 0, result.output
        assert "presentation" in result.output
        assert "bias" in result.output

    def test_explain_unknown_peptide_exits_nonzero(self, ranked):
        out, _ = ranked
        result = runner.invoke(
            app, ["explain", str(out / "EXAMPLE-PDAC-001.features.json"), "WWWWWWWWW"]
        )
        assert result.exit_code == 1

    def test_catalog_marks_hla_matches(self):
        result = runner.invoke(app, ["catalog", "--hla", "A*11:01,C*08:02"])
        assert result.exit_code == 0, result.output
        assert "* KRAS G12D" in result.output
        assert "  KRAS G12R" in result.output

    def test_init_config_writes_loadable_yaml(self, tmp_path):
        from neoantigene.config import PipelineConfig

        target = tmp_path / "conf.yaml"
        result = runner.invoke(app, ["init-config", "-o", str(target)])
        assert result.exit_code == 0, result.output
        assert PipelineConfig.load(target) == PipelineConfig()


def test_benchmark_compares_against_every_baseline(tmp_path, cohort, benchmark_run):
    sample_id = benchmark_run.sample.sample_id
    completed = tmp_path / "completed.tsv"
    _write_results(simulate_assay_results(benchmark_run.all_scored, sample_id, seed=5), completed)

    config = tmp_path / "bench.yaml"
    config.write_text(
        "variant_filters:\n  min_ccf: 0.0\n  min_dna_vaf: 0.0\n  min_tumor_depth: 0\n"
        "expression:\n  min_tpm: 0.0\n  min_rna_vaf: null\n"
        # Quoted: bare `null` is YAML's None, which is not a backend name.
        'presentation:\n  backend: "null"\n  max_affinity_percentile: 100.0\n'
        "peptides:\n  lengths: [9]\n  drop_self_matching: false\n"
        "output:\n  top_n: 1000\n  max_per_variant: 1000\n"
    )

    result = runner.invoke(
        app,
        [
            "benchmark",
            str(cohort.patient(0).manifest_path),
            str(completed),
            "--config",
            str(config),
            "--backend",
            "null",
            "--allow-null-backend",
            "--k",
            "30",
        ],
    )
    assert result.exit_code == 0, result.output
    for method in ("neoantigene", "binding_only", "arbitrary"):
        assert method in result.output


def _write_completed_from_ranked(ranked_tsv, out) -> None:
    with open(ranked_tsv) as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    lines = ["sample_id\tpeptide\tallele\tassay\tcall"]
    for index, row in enumerate(rows):
        call = "positive" if index % 2 == 0 else "negative"
        lines.append(
            f"{row['sample_id']}\t{row['mutant_peptide']}\t{row['allele']}\tifng_elispot\t{call}"
        )
    out.write_text("\n".join(lines) + "\n")


def _write_results(results, out) -> None:
    lines = ["sample_id\tpeptide\tallele\tassay\tcall\teffect_size"]
    for record in results:
        effect = "" if record.call is AssayCall.INDETERMINATE else record.effect_size
        lines.append(
            f"{record.sample_id}\t{record.peptide}\t{record.allele}\t"
            f"{record.assay.value}\t{record.call.value}\t{effect}"
        )
    out.write_text("\n".join(lines) + "\n")
