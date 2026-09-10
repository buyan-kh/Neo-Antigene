"""Run identity and reproducibility.

Every ranking run gets an id and writes a manifest recording exactly what
produced the shortlist: package version, resolved config, backend, and a
digest of every input file. A lab that spent money synthesizing peptides is
entitled to reconstruct the run that nominated them, and a pharma partner
auditing the platform will ask for precisely this.
"""

from __future__ import annotations

import hashlib
import json
import logging
import platform
import sys
from collections.abc import Sequence
from datetime import UTC, datetime
from pathlib import Path
from typing import Any
from uuid import uuid4

from pydantic import BaseModel, Field

from . import __version__
from .config import PipelineConfig

logger = logging.getLogger(__name__)

_DIGEST_CHUNK = 1 << 20


def file_digest(path: Path) -> str:
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        while chunk := handle.read(_DIGEST_CHUNK):
            digest.update(chunk)
    return f"sha256:{digest.hexdigest()}"


def config_digest(config: PipelineConfig) -> str:
    payload = json.dumps(config.model_dump(mode="json"), sort_keys=True).encode()
    return f"sha256:{hashlib.sha256(payload).hexdigest()}"


def new_run_id() -> str:
    """Sortable, human-legible run id: `20260909T184233Z-1a2b3c4d`."""
    stamp = datetime.now(UTC).strftime("%Y%m%dT%H%M%SZ")
    return f"{stamp}-{uuid4().hex[:8]}"


class InputDigest(BaseModel):
    path: str
    digest: str
    bytes: int


class RunManifest(BaseModel):
    """Everything needed to explain or reproduce one ranking run."""

    run_id: str
    started_at: datetime
    finished_at: datetime | None = None
    sample_id: str
    package_version: str
    python_version: str
    platform: str
    backend: str
    config_digest: str
    config: PipelineConfig
    inputs: list[InputDigest] = Field(default_factory=list)
    counts: dict[str, int] = Field(default_factory=dict)

    @classmethod
    def start(
        cls,
        sample_id: str,
        config: PipelineConfig,
        backend: str,
        inputs: Sequence[Path],
        run_id: str | None = None,
    ) -> RunManifest:
        return cls(
            run_id=run_id or new_run_id(),
            started_at=datetime.now(UTC),
            sample_id=sample_id,
            package_version=__version__,
            python_version=sys.version.split()[0],
            platform=platform.platform(),
            backend=backend,
            config_digest=config_digest(config),
            config=config,
            inputs=[
                InputDigest(
                    path=str(path),
                    digest=file_digest(path),
                    bytes=path.stat().st_size,
                )
                for path in inputs
            ],
        )

    def finish(self, counts: dict[str, int]) -> RunManifest:
        return self.model_copy(update={"finished_at": datetime.now(UTC), "counts": dict(counts)})

    def write(self, path: Path) -> Path:
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(self.model_dump_json(indent=2))
        return path


class RunIdFilter(logging.Filter):
    """Stamps every record with the active run id."""

    def __init__(self, run_id: str) -> None:
        super().__init__()
        self.run_id = run_id

    def filter(self, record: logging.LogRecord) -> bool:
        record.run_id = self.run_id
        return True


def configure_logging(run_id: str, verbose: bool = False, json_logs: bool = False) -> None:
    handler: logging.Handler = logging.StreamHandler(sys.stderr)
    handler.addFilter(RunIdFilter(run_id))
    handler.setFormatter(_JsonFormatter() if json_logs else _TextFormatter())

    root = logging.getLogger("neoantigene")
    root.handlers.clear()
    root.addHandler(handler)
    root.setLevel(logging.DEBUG if verbose else logging.INFO)
    root.propagate = False


class _TextFormatter(logging.Formatter):
    def __init__(self) -> None:
        super().__init__("%(asctime)s %(levelname)-7s [%(run_id)s] %(name)s: %(message)s")


class _JsonFormatter(logging.Formatter):
    def format(self, record: logging.LogRecord) -> str:
        payload: dict[str, Any] = {
            "ts": datetime.fromtimestamp(record.created, UTC).isoformat(),
            "level": record.levelname,
            "run_id": getattr(record, "run_id", None),
            "logger": record.name,
            "message": record.getMessage(),
        }
        if record.exc_info:
            payload["exception"] = self.formatException(record.exc_info)
        return json.dumps(payload)
