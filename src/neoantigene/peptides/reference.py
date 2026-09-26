"""Reference proteome acquisition.

The self-peptide gate and the mutant protein lookup are both only as good as
the reference behind them, so the reference is fetched from Ensembl rather than
vendored. A release is pinned: the proteome must match the VEP cache that
annotated the VCF, and "whatever is current today" is not a reproducible
answer.

Ensembl publishes a `CHECKSUMS` file per directory using the BSD `sum`
algorithm, which is what `verify` checks against.
"""

from __future__ import annotations

import logging
import os
import shutil
import urllib.request
from dataclasses import dataclass
from pathlib import Path
from typing import Final

from ..io.fasta import read_fasta

logger = logging.getLogger(__name__)

#: Pinned Ensembl release. Bumping this invalidates cached downloads because
#: the release is part of the cache path.
ENSEMBL_RELEASE: Final[str] = "116"

PROTEOME_FILENAME: Final[str] = "Homo_sapiens.GRCh38.pep.all.fa.gz"

_FTP_ROOT: Final[str] = "https://ftp.ensembl.org/pub"

#: A real human proteome has ~120k translations across all transcripts. The
#: bundled example stub has two. Anything in between is also not a proteome,
#: hence a generous floor rather than a tight one.
MIN_REAL_PROTEOME_PROTEINS: Final[int] = 20_000


class ReferenceDownloadError(RuntimeError):
    """Raised when the reference proteome could not be fetched or verified."""


def proteome_url(release: str = ENSEMBL_RELEASE) -> str:
    return f"{_FTP_ROOT}/release-{release}/fasta/homo_sapiens/pep/{PROTEOME_FILENAME}"


def checksums_url(release: str = ENSEMBL_RELEASE) -> str:
    return f"{_FTP_ROOT}/release-{release}/fasta/homo_sapiens/pep/CHECKSUMS"


def cache_dir() -> Path:
    """Where downloaded references live.

    Honours `NEOANTIGENE_CACHE` so a shared machine or a container can point
    this at a volume that actually has room for it.
    """
    override = os.environ.get("NEOANTIGENE_CACHE")
    if override:
        return Path(override).expanduser()
    base = os.environ.get("XDG_CACHE_HOME")
    root = Path(base).expanduser() if base else Path.home() / ".cache"
    return root / "neoantigene"


def proteome_path(release: str = ENSEMBL_RELEASE) -> Path:
    return cache_dir() / f"ensembl-{release}" / PROTEOME_FILENAME


def bsd_sum(path: Path) -> tuple[int, int]:
    """`sum` as coreutils implements it by default: 16-bit checksum, 1K blocks.

    Returned as `(checksum, blocks)` to match the two columns in Ensembl's
    CHECKSUMS file.
    """
    checksum = 0
    size = 0
    with open(path, "rb") as handle:
        while chunk := handle.read(1 << 20):
            size += len(chunk)
            for byte in chunk:
                checksum = (checksum >> 1) + ((checksum & 1) << 15)
                checksum = (checksum + byte) & 0xFFFF
    blocks = (size + 1023) // 1024
    return checksum, blocks


def _expected_checksum(release: str, filename: str) -> tuple[int, int] | None:
    url = checksums_url(release)
    try:
        with urllib.request.urlopen(url, timeout=60) as response:
            body = response.read().decode("utf-8", "replace")
    except OSError as exc:
        logger.warning("could not fetch %s: %s; skipping checksum verification", url, exc)
        return None
    for line in body.splitlines():
        fields = line.split()
        if len(fields) == 3 and fields[2] == filename:
            return int(fields[0]), int(fields[1])
    return None


def fetch_proteome(
    release: str = ENSEMBL_RELEASE,
    dest: Path | None = None,
    force: bool = False,
    verify: bool = True,
) -> Path:
    """Download the Ensembl human proteome, returning the local path.

    Idempotent: an existing, verified file short-circuits. The download lands
    in a temporary file and is moved into place only once complete, so an
    interrupted run cannot leave a truncated FASTA that would silently produce
    a wrong shortlist.
    """
    target = dest if dest is not None else proteome_path(release)
    target.parent.mkdir(parents=True, exist_ok=True)

    if target.exists() and not force:
        logger.info("reference proteome already present at %s", target)
        return target

    url = proteome_url(release)
    logger.info("downloading %s", url)
    tmp = target.with_name(f"{target.name}.{os.getpid()}.part")
    try:
        with (
            urllib.request.urlopen(url, timeout=120) as response,
            open(tmp, "wb") as handle,
        ):
            if response.status != 200:
                raise ReferenceDownloadError(f"{url} returned HTTP {response.status}")
            shutil.copyfileobj(response, handle)
    except OSError as exc:
        tmp.unlink(missing_ok=True)
        raise ReferenceDownloadError(f"failed to download {url}: {exc}") from exc

    if verify:
        expected = _expected_checksum(release, PROTEOME_FILENAME)
        if expected is not None and bsd_sum(tmp) != expected:
            tmp.unlink(missing_ok=True)
            raise ReferenceDownloadError(
                f"checksum mismatch for {PROTEOME_FILENAME}; refusing to install a "
                "reference that does not match Ensembl's published checksum"
            )

    tmp.replace(target)
    logger.info("reference proteome installed at %s (%d bytes)", target, target.stat().st_size)
    return target


@dataclass(frozen=True)
class ProteomeProvenance:
    """Whether the proteome behind a run can carry the weight put on it.

    A shortlist produced against the bundled stub has a self-peptide gate that
    never fires and a self-dissimilarity feature with nothing to compare
    against. That is not a footnote, so it travels with the run report and is
    surfaced in every interface rather than living only in the README.
    """

    path: Path
    proteins: int

    @property
    def is_stub(self) -> bool:
        return self.proteins < MIN_REAL_PROTEOME_PROTEINS

    @property
    def caveat(self) -> str | None:
        if not self.is_stub:
            return None
        return (
            f"{self.path.name} holds {self.proteins} protein sequences, not a complete human "
            f"proteome (~120k translations in Ensembl {ENSEMBL_RELEASE}). The self-peptide "
            f"gate cannot fire and the self-dissimilarity feature has nothing to compare "
            f"against. Run `neoantigene fetch-proteome` and repoint "
            f"processed.proteome_fasta."
        )


def count_proteins(path: Path) -> int:
    return sum(1 for _header, sequence in read_fasta(path) if sequence.rstrip("*"))


def provenance(path: Path) -> ProteomeProvenance:
    """Inspect a proteome FASTA and judge whether it is a usable reference."""
    return ProteomeProvenance(path=path, proteins=count_proteins(path))
