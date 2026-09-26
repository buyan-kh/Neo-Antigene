"""Test setup for the API app.

Deliberately not wired into the repo's root `pytest` run: `testpaths` in
pyproject.toml stays `tests`, so CI keeps testing the package alone and this
optional app cannot break it. Run these explicitly:

    uv run pytest apps/api/tests
"""

from __future__ import annotations

import sys
from collections.abc import Iterator
from pathlib import Path

import pytest

APP_ROOT = Path(__file__).resolve().parents[1]
if str(APP_ROOT) not in sys.path:
    sys.path.insert(0, str(APP_ROOT))

from fastapi.testclient import TestClient  # noqa: E402
from neoantigene_api.main import app  # noqa: E402


@pytest.fixture
def client() -> Iterator[TestClient]:
    with TestClient(app) as test_client:
        yield test_client


@pytest.fixture
def finished_run(client: TestClient) -> str:
    """One completed run on the bundled example, using the dev backend.

    The dev backend keeps the suite fast and independent of downloaded
    MHCflurry weights. Its scores are meaningless, so nothing here asserts on
    score *values* — only on plumbing, shapes and guard rails.
    """
    response = client.post(
        "/api/runs",
        json={
            "source": "example",
            "backend": "null",
            "acknowledge_development_backend": True,
        },
    )
    assert response.status_code == 202, response.text
    run_id = response.json()["run_id"]

    for _ in range(600):
        state = client.get(f"/api/runs/{run_id}").json()
        if state["status"] in {"succeeded", "failed"}:
            break
    assert state["status"] == "succeeded", state.get("error")
    return run_id
