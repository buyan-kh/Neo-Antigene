"""In-memory run registry.

Runs take real time, so a POST starts one on a worker thread and returns an id
the client polls. State lives in a dict: this is a demo service, and a database
would add operational surface without making anything on screen more true.

Restarting the API loses run history. That is the documented trade.
"""

from __future__ import annotations

import logging
import shutil
import tempfile
import threading
import traceback
from collections.abc import Iterator
from concurrent.futures import ThreadPoolExecutor
from datetime import UTC, datetime
from pathlib import Path
from typing import Final

from neoantigene.config import PipelineConfig
from neoantigene.io.manifest import Sample
from neoantigene.io.writers import write_features_json, write_tsv
from neoantigene.models import ScoredCandidate
from neoantigene.pipeline import RunResult
from neoantigene.pipeline import run as run_pipeline
from neoantigene.presentation.registry import DEVELOPMENT_BACKENDS, get_backend
from neoantigene.run import new_run_id

from . import progress
from .schemas import RunState, RunStatus, RunWarnings

logger = logging.getLogger(__name__)

#: One run at a time. MHCflurry holds a TensorFlow graph and the machine doing
#: the screen share is also running a browser.
MAX_CONCURRENT_RUNS: Final[int] = 1


class RunNotFound(KeyError):
    pass


class Run:
    """One pipeline execution and everything the API needs to answer about it."""

    def __init__(self, run_id: str, out_dir: Path) -> None:
        self.state = RunState(
            run_id=run_id,
            status=RunStatus.QUEUED,
            created_at=datetime.now(UTC),
        )
        self.out_dir = out_dir
        self.result: RunResult | None = None
        self.ranked_tsv: Path | None = None
        self.features_json: Path | None = None
        self.lock = threading.Lock()

    def scored(self) -> list[ScoredCandidate]:
        return list(self.result.all_scored) if self.result else []


class RunStore:
    """Registry of runs, their workspaces and their outputs."""

    def __init__(self, root: Path | None = None) -> None:
        self._root = root or Path(tempfile.mkdtemp(prefix="neoantigene-api-"))
        self._runs: dict[str, Run] = {}
        self._lock = threading.Lock()
        self._pool = ThreadPoolExecutor(
            max_workers=MAX_CONCURRENT_RUNS, thread_name_prefix="neoantigene-run"
        )
        progress.install(self._record_event)

    @property
    def root(self) -> Path:
        return self._root

    def workspace(self, workspace_id: str) -> Path:
        path = self._root / "workspaces" / workspace_id
        path.mkdir(parents=True, exist_ok=True)
        return path

    def new_workspace(self) -> tuple[str, Path]:
        workspace_id = new_run_id()
        return workspace_id, self.workspace(workspace_id)

    def get(self, run_id: str) -> Run:
        with self._lock:
            run = self._runs.get(run_id)
        if run is None:
            raise RunNotFound(run_id)
        return run

    def list_runs(self) -> list[RunState]:
        with self._lock:
            runs = list(self._runs.values())
        return sorted((r.state for r in runs), key=lambda s: s.created_at, reverse=True)

    def submit(
        self,
        sample: Sample,
        config: PipelineConfig,
        warnings: RunWarnings,
    ) -> str:
        run_id = new_run_id()
        run = Run(run_id, self._root / "runs" / run_id)
        run.state.sample_id = sample.sample_id
        run.state.backend = config.presentation.backend
        run.state.warnings = warnings
        run.state.config = config
        with self._lock:
            self._runs[run_id] = run
        self._pool.submit(self._execute, run, sample, config)
        return run_id

    def _record_event(self, run_id: str, event: progress.ProgressEvent) -> None:
        with self._lock:
            run = self._runs.get(run_id)
        if run is None:
            return
        with run.lock:
            run.state.events.append(event)

    def _execute(self, run: Run, sample: Sample, config: PipelineConfig) -> None:
        with run.lock:
            run.state.status = RunStatus.RUNNING
            run.state.started_at = datetime.now(UTC)

        try:
            with progress.attributed_to(run.state.run_id):
                result = run_pipeline(
                    sample,
                    config,
                    backend=get_backend(config.presentation.backend),
                    run_id=run.state.run_id,
                )
                run.out_dir.mkdir(parents=True, exist_ok=True)
                ranked_tsv = run.out_dir / f"{sample.sample_id}.ranked.tsv"
                features_json = run.out_dir / f"{sample.sample_id}.features.json"
                write_tsv(result.ranked, sample.sample_id, ranked_tsv, run_id=run.state.run_id)
                write_features_json(
                    result.all_scored,
                    sample.sample_id,
                    features_json,
                    run_id=run.state.run_id,
                )
                result.manifest.write(run.out_dir / f"{sample.sample_id}.run.json")
        except Exception as exc:
            logger.exception("run %s failed", run.state.run_id)
            with run.lock:
                run.state.status = RunStatus.FAILED
                run.state.finished_at = datetime.now(UTC)
                # The real message, not a generic apology: a failed run on a
                # screen share is a debugging session, not an error page.
                run.state.error = str(exc) or repr(exc)
                run.state.error_type = type(exc).__name__
                run.state.events.append(
                    progress.ProgressEvent(
                        at=datetime.now(UTC),
                        level="ERROR",
                        stage="pipeline",
                        message="".join(traceback.format_exception_only(type(exc), exc)).strip(),
                    )
                )
            return

        with run.lock:
            run.result = result
            run.ranked_tsv = ranked_tsv
            run.features_json = features_json
            run.state.report = result.report
            run.state.config = result.config
            run.state.shortlisted = len(result.ranked)
            run.state.warnings = _refresh_warnings(run.state.warnings, result)
            run.state.status = RunStatus.SUCCEEDED
            run.state.finished_at = datetime.now(UTC)

    def shutdown(self) -> None:
        self._pool.shutdown(wait=False, cancel_futures=True)
        shutil.rmtree(self._root, ignore_errors=True)

    def __iter__(self) -> Iterator[Run]:
        with self._lock:
            return iter(list(self._runs.values()))


def _refresh_warnings(warnings: RunWarnings, result: RunResult) -> RunWarnings:
    """Replace pre-run guesses with what the run actually observed."""
    report = result.report
    updated = warnings.model_copy()
    updated.proteome_is_stub = report.proteome_is_stub
    updated.proteome_proteins = report.proteome_proteins
    if report.proteome_is_stub:
        updated.proteome_detail = (
            f"The reference proteome holds {report.proteome_proteins} protein sequences, "
            "not a complete human proteome. The self-peptide gate cannot fire and the "
            "self-dissimilarity feature was skipped rather than guessed."
        )
    else:
        updated.proteome_detail = None
    updated.expression_available = report.expression_available
    if not report.expression_available:
        updated.expression_detail = (
            "No RNA expression table was supplied, so the expression and tumor "
            "selectivity features are neutral for every candidate and contribute "
            "nothing to the ranking."
        )
    else:
        updated.expression_detail = None
    return updated


def development_warning(backend: str) -> tuple[bool, str | None]:
    """The same guard the CLI applies, phrased for a banner."""
    if backend not in DEVELOPMENT_BACKENDS:
        return False, None
    return True, (
        f"The {backend!r} backend is a deterministic hash-based stand-in with no "
        "biological validity. Every score on this page is meaningless and must not be "
        "read as a prediction."
    )
