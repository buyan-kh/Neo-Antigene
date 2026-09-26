"""Turning the pipeline's own log stream into progress events.

`neoantigene.pipeline.run` is a synchronous function with no progress callback,
and adding one would mean changing the library to suit a demo. It does, however,
log every stage it completes. So the API attaches a handler to the `neoantigene`
logger and republishes those records as progress events.

The consequence worth knowing: progress reflects what the pipeline actually
said it did, not a script of what a frontend hopes it will do. If a stage stops
logging, the step disappears from the UI rather than showing a fabricated tick.
"""

from __future__ import annotations

import logging
import threading
from collections.abc import Callable
from datetime import UTC, datetime

from pydantic import BaseModel

LIBRARY_LOGGER = "neoantigene"

_active = threading.local()


class ProgressEvent(BaseModel):
    """One line the pipeline emitted while running."""

    at: datetime
    level: str
    stage: str
    message: str


class _Router(logging.Handler):
    """Routes each record to whichever run's thread produced it."""

    def __init__(self, sink: Callable[[str, ProgressEvent], None]) -> None:
        super().__init__(level=logging.INFO)
        self._sink = sink

    def emit(self, record: logging.LogRecord) -> None:
        run_id = getattr(_active, "run_id", None)
        if run_id is None:
            return
        try:
            message = record.getMessage()
        except Exception:  # pragma: no cover - defensive, never break a run to log
            return
        self._sink(
            run_id,
            ProgressEvent(
                at=datetime.fromtimestamp(record.created, UTC),
                level=record.levelname,
                stage=record.name.removeprefix(f"{LIBRARY_LOGGER}."),
                message=message,
            ),
        )


def install(sink: Callable[[str, ProgressEvent], None]) -> _Router:
    """Attach the router to the library logger, replacing any prior one.

    Deliberately does not call `neoantigene.run.configure_logging`, which
    clears handlers and assumes one run per process.
    """
    logger = logging.getLogger(LIBRARY_LOGGER)
    for existing in [h for h in logger.handlers if isinstance(h, _Router)]:
        logger.removeHandler(existing)
    handler = _Router(sink)
    logger.addHandler(handler)
    logger.setLevel(logging.INFO)
    logger.propagate = True
    return handler


class attributed_to:
    """Marks the calling thread's log records as belonging to one run."""

    def __init__(self, run_id: str) -> None:
        self._run_id = run_id
        self._previous: str | None = None

    def __enter__(self) -> None:
        self._previous = getattr(_active, "run_id", None)
        _active.run_id = self._run_id

    def __exit__(self, *_: object) -> None:
        _active.run_id = self._previous
