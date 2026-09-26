"""FastAPI wrapper around the `neoantigene` package.

Lives outside `src/neoantigene` on purpose. The package is the product and has
no web dependencies; this app depends on the package and the package never
depends on it.
"""

from __future__ import annotations

__all__ = ["app"]


def __getattr__(name: str) -> object:
    # Imported lazily so `neoantigene_api` stays importable for tooling that
    # only wants the module docstring, without pulling in FastAPI.
    if name == "app":
        from .main import app

        return app
    raise AttributeError(name)
