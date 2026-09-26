"""API behaviour, with the emphasis on the things that could mislead a viewer.

The demo's whole value is that a partner can trust what is on screen, so the
tests that matter most here are the ones asserting the API refuses to invent
anything: no benchmark without labels, no dev-backend scores without a banner,
no silent stub proteome.
"""

from __future__ import annotations

import io

from fastapi.testclient import TestClient

from neoantigene.io.writers import COLUMNS
from neoantigene.presentation.registry import DEVELOPMENT_BACKENDS
from neoantigene.scoring.rank import SCORED_FEATURES


class TestDiscovery:
    def test_health(self, client: TestClient) -> None:
        assert client.get("/api/health").json()["status"] == "ok"

    def test_backends_reports_development_ones_separately(self, client: TestClient) -> None:
        payload = client.get("/api/backends").json()
        assert set(payload["development"]) == set(DEVELOPMENT_BACKENDS)
        assert "mhcflurry" in payload["available"]

    def test_reference_reports_the_bundled_stub_honestly(self, client: TestClient) -> None:
        payload = client.get("/api/reference").json()
        assert payload["bundled_stub_proteins"] < 10
        assert "fetch-proteome" in payload["fetch_command"]

    def test_weights_match_the_scored_feature_set(self, client: TestClient) -> None:
        payload = client.get("/api/weights").json()
        assert set(payload) == {"bias", *SCORED_FEATURES}

    def test_example_exposes_the_bundled_manifest(self, client: TestClient) -> None:
        payload = client.get("/api/example").json()
        assert payload["sample_id"]
        assert payload["hla"]
        assert payload["variant_tsv"]


class TestGuardRails:
    def test_development_backend_is_refused_without_acknowledgement(
        self, client: TestClient
    ) -> None:
        """Mirrors the CLI's --allow-null-backend guard."""
        response = client.post("/api/runs", json={"source": "example", "backend": "null"})
        assert response.status_code == 400
        assert "meaningless" in response.json()["detail"]

    def test_unknown_backend_is_rejected(self, client: TestClient) -> None:
        response = client.post("/api/runs", json={"source": "example", "backend": "nope"})
        assert response.status_code == 400

    def test_acknowledged_development_run_carries_a_banner(
        self, client: TestClient, finished_run: str
    ) -> None:
        warnings = client.get(f"/api/runs/{finished_run}").json()["warnings"]
        assert warnings["development_backend"] is True
        assert "meaningless" in warnings["development_detail"]

    def test_a_stub_proteome_is_disclosed_on_the_run(
        self, client: TestClient, finished_run: str
    ) -> None:
        warnings = client.get(f"/api/runs/{finished_run}").json()["warnings"]
        assert warnings["proteome_is_stub"] is True
        assert "not a complete human proteome" in warnings["proteome_detail"]

    def test_unknown_run_is_404(self, client: TestClient) -> None:
        assert client.get("/api/runs/nope").status_code == 404
        assert client.get("/api/runs/nope/results").status_code == 404

    def test_uploads_must_name_a_workspace(self, client: TestClient) -> None:
        response = client.post(
            "/api/runs",
            json={"source": "upload", "sample": {"sample_id": "X", "hla": ["HLA-A*02:01"]}},
        )
        assert response.status_code == 422


class TestProgress:
    def test_events_come_from_the_pipeline_itself(
        self, client: TestClient, finished_run: str
    ) -> None:
        state = client.get(f"/api/runs/{finished_run}").json()
        messages = " ".join(event["message"] for event in state["events"])
        assert "protein-altering variants" in messages
        assert "somatic filtering" in messages
        assert "peptide-allele pairs" in messages
        assert "shortlisted" in messages

    def test_report_counts_are_populated(self, client: TestClient, finished_run: str) -> None:
        report = client.get(f"/api/runs/{finished_run}").json()["report"]
        assert report["variants_read"] > 0
        assert report["peptides_generated"] > 0
        assert report["pairs_scored"] > 0


class TestResults:
    def test_columns_are_exactly_the_writer_columns(
        self, client: TestClient, finished_run: str
    ) -> None:
        """The UI must not invent or reorder columns."""
        payload = client.get(f"/api/runs/{finished_run}/results").json()
        assert payload["columns"] == list(COLUMNS)

    def test_rows_carry_every_column(self, client: TestClient, finished_run: str) -> None:
        payload = client.get(f"/api/runs/{finished_run}/results").json()
        assert payload["rows"]
        for row in payload["rows"]:
            assert set(COLUMNS) <= set(row)

    def test_mutant_and_wildtype_peptides_differ_where_present(
        self, client: TestClient, finished_run: str
    ) -> None:
        for row in client.get(f"/api/runs/{finished_run}/results").json()["rows"]:
            if row["wildtype_peptide"]:
                assert row["wildtype_peptide"] != row["mutant_peptide"]
                assert len(row["wildtype_peptide"]) == len(row["mutant_peptide"])

    def test_gated_rows_name_their_gate(self, client: TestClient, finished_run: str) -> None:
        payload = client.get(f"/api/runs/{finished_run}/results").json()
        gated = [r for r in payload["rows"] if r["gate_failures"]]
        for row in gated:
            assert row["gate_failures"] in payload["gate_failures"] or ";" in row["gate_failures"]

    def test_excluding_gated_rows_shortens_the_table(
        self, client: TestClient, finished_run: str
    ) -> None:
        with_gated = client.get(f"/api/runs/{finished_run}/results").json()["total"]
        without = client.get(
            f"/api/runs/{finished_run}/results", params={"include_gated": False}
        ).json()["total"]
        assert without <= with_gated

    def test_downloads_are_served(self, client: TestClient, finished_run: str) -> None:
        tsv = client.get(f"/api/runs/{finished_run}/download/ranked.tsv")
        assert tsv.status_code == 200
        assert tsv.text.splitlines()[0].split("\t") == list(COLUMNS)
        assert client.get(f"/api/runs/{finished_run}/download/features.json").status_code == 200
        assert client.get(f"/api/runs/{finished_run}/download/nope").status_code == 404


class TestExplain:
    def test_contributions_reconstruct_the_logit(
        self, client: TestClient, finished_run: str
    ) -> None:
        """The side-sheet chart must add up, or it is decoration."""
        row = client.get(f"/api/runs/{finished_run}/results").json()["rows"][0]
        payload = client.get(
            f"/api/runs/{finished_run}/explain",
            params={"peptide": row["mutant_peptide"], "allele": row["allele"]},
        ).json()

        total = payload["bias"] + sum(c["contribution"] for c in payload["contributions"])
        assert abs(total - payload["logit"]) < 1e-9

    def test_every_scored_feature_appears(self, client: TestClient, finished_run: str) -> None:
        row = client.get(f"/api/runs/{finished_run}/results").json()["rows"][0]
        payload = client.get(
            f"/api/runs/{finished_run}/explain",
            params={"peptide": row["mutant_peptide"], "allele": row["allele"]},
        ).json()
        assert {c["feature"] for c in payload["contributions"]} == set(SCORED_FEATURES)

    def test_unknown_peptide_is_404(self, client: TestClient, finished_run: str) -> None:
        response = client.get(
            f"/api/runs/{finished_run}/explain",
            params={"peptide": "WWWWWWWWW", "allele": "HLA-A*02:01"},
        )
        assert response.status_code == 404


class TestBenchmark:
    """The benchmark endpoint must refuse to produce a number it cannot justify."""

    def test_unparseable_labels_are_rejected(self, client: TestClient, finished_run: str) -> None:
        response = client.post(
            f"/api/runs/{finished_run}/benchmark",
            files={
                "labels": (
                    "labels.tsv",
                    io.BytesIO(b"not\ta\tlabel\tfile\n"),
                    "text/tab-separated-values",
                )
            },
        )
        assert response.status_code == 422

    def test_labels_that_match_nothing_are_rejected(
        self, client: TestClient, finished_run: str
    ) -> None:
        """No silent empty benchmark: a zero would read as a real result."""
        sheet = (
            "sample_id\tpeptide\tallele\tassay\tcall\n"
            "OTHER\tWWWWWWWWW\tHLA-A*02:01\tifng_elispot\tpositive\n"
        )
        response = client.post(
            f"/api/runs/{finished_run}/benchmark",
            files={
                "labels": ("labels.tsv", io.BytesIO(sheet.encode()), "text/tab-separated-values")
            },
        )
        assert response.status_code == 422
        assert "match a" in response.json()["detail"]

    def test_matching_labels_produce_the_four_numbers(
        self, client: TestClient, finished_run: str
    ) -> None:
        """Labels here are synthetic and only exercise the arithmetic path."""
        rows = client.get(f"/api/runs/{finished_run}/results").json()["rows"][:6]
        lines = ["sample_id\tpeptide\tallele\tassay\tcall"]
        for index, row in enumerate(rows):
            call = "positive" if index % 3 == 0 else "negative"
            lines.append(
                f"{row['sample_id']}\t{row['mutant_peptide']}\t{row['allele']}\tifng_elispot\t{call}"
            )
        sheet = "\n".join(lines) + "\n"

        payload = client.post(
            f"/api/runs/{finished_run}/benchmark",
            files={
                "labels": ("labels.tsv", io.BytesIO(sheet.encode()), "text/tab-separated-values")
            },
        ).json()

        methods = {entry["method"] for entry in payload["retrieval"]}
        assert "neoantigene" in methods
        assert "binding_only" in methods
        for entry in payload["retrieval"]:
            assert 0.0 <= entry["recall_at_10"] <= 1.0
            assert 0.0 <= entry["precision_at_10"] <= 1.0
            assert entry["positives"] == payload["source"]["positives"]
        assert payload["source"]["filename"] == "labels.tsv"
        assert payload["pool_size"] > 0

    def test_every_method_is_scored_on_the_same_pool(
        self, client: TestClient, finished_run: str
    ) -> None:
        rows = client.get(f"/api/runs/{finished_run}/results").json()["rows"][:8]
        lines = ["sample_id\tpeptide\tallele\tassay\tcall"]
        for index, row in enumerate(rows):
            call = "positive" if index % 2 == 0 else "negative"
            lines.append(
                f"{row['sample_id']}\t{row['mutant_peptide']}\t{row['allele']}\tifng_elispot\t{call}"
            )
        payload = client.post(
            f"/api/runs/{finished_run}/benchmark",
            files={
                "labels": (
                    "labels.tsv",
                    io.BytesIO(("\n".join(lines) + "\n").encode()),
                    "text/tab-separated-values",
                )
            },
        ).json()
        assert len({entry["assayed"] for entry in payload["retrieval"]}) == 1
        assert len({entry["positives"] for entry in payload["retrieval"]}) == 1
