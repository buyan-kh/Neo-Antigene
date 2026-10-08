"""The Vercel service entrypoint must be the same FastAPI app the demo runs."""

from __future__ import annotations

from fastapi import FastAPI


def test_vercel_entrypoint_is_the_api() -> None:
    import app as vercel_app

    assert isinstance(vercel_app.app, FastAPI)
    paths = {getattr(route, "path", None) for route in vercel_app.app.routes}
    assert "/api/health" in paths
    assert "/api/runs" in paths
