from __future__ import annotations

import gzip
from collections.abc import Iterator
from pathlib import Path
from typing import IO


def _open(path: Path) -> IO[str]:
    if str(path).endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path)


def read_fasta(path: Path) -> Iterator[tuple[str, str]]:
    """Yield (header, sequence) pairs. Header excludes the leading '>'."""
    header: str | None = None
    chunks: list[str] = []
    with _open(path) as handle:
        for raw_line in handle:
            line = raw_line.rstrip("\n")
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(chunks)
                header = line[1:]
                chunks = []
            else:
                chunks.append(line.strip())
    if header is not None:
        yield header, "".join(chunks)
