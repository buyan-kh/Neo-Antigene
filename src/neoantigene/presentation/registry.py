from __future__ import annotations

from collections.abc import Callable

from .base import PresentationBackend
from .null import NullBackend

BackendFactory = Callable[[], PresentationBackend]

_BACKENDS: dict[str, BackendFactory] = {}

#: Backends that produce meaningless scores and must be opted into explicitly.
DEVELOPMENT_BACKENDS = frozenset({"null"})


def register(name: str, factory: BackendFactory) -> None:
    _BACKENDS[name] = factory


def available() -> list[str]:
    return sorted(_BACKENDS)


def get_backend(name: str) -> PresentationBackend:
    if name not in _BACKENDS:
        raise KeyError(f"unknown presentation backend {name!r}; available: {available()}")
    return _BACKENDS[name]()


def _mhcflurry() -> PresentationBackend:
    from .mhcflurry_backend import MHCflurryBackend

    return MHCflurryBackend()


register("mhcflurry", _mhcflurry)
register("null", NullBackend)
