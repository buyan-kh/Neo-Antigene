"""ASGI entrypoint for the Vercel ``api`` service.

``scripts/dev.sh`` still launches ``neoantigene_api.main:app`` with uvicorn's
``--app-dir``. This module exists because the Vercel service root is the
repository (where ``pyproject.toml`` and ``src/`` live), and the FastAPI
package itself lives in ``apps/api``.
"""

from __future__ import annotations

import sys
from pathlib import Path
from typing import Any


def _load() -> Any:
    api_dir = str(Path(__file__).resolve().parent / "apps" / "api")
    if api_dir not in sys.path:
        sys.path.insert(0, api_dir)
    from neoantigene_api.main import app as fastapi_app

    return fastapi_app


app = _load()
