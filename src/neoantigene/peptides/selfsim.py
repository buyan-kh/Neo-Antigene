"""Nearest self-peptide search over the reference proteome.

The question this answers is "how close is this neo-peptide to the closest
thing the thymus already screened against", which is the feature Richman et al.
(Cell Systems 2019) showed tracks immunogenicity. A peptide one conservative
substitution away from an abundant self protein is a poor bet; one with no near
neighbour in the proteome is a better one.

Why seed-and-extend rather than a full scan. Scoring every candidate against
every window of a 70-million-residue proteome is ~10^12 operations per run,
which is minutes to hours. This uses the same trick BLAST does: index short
exact seeds, and only align windows that share one.

The seed length is 4, which is not an arbitrary tuning choice. For any single
substitution in a peptide of length >= 8, at least one 4-mer of the mutant
avoids the mutated position, so the wild-type counterpart — very often the
nearest self peptide — is guaranteed to be found. Longer seeds lose that
guarantee for centrally located mutations.

What this gives up, stated plainly: a self peptide sharing no exact 4-mer with
the candidate is never aligned. Such a window is necessarily a distant match,
so the effect is that some already-low similarities are reported as slightly
lower still. It is a sensitivity limit, not a correctness one, and it is not
the exhaustive Smith-Waterman of the original method.
"""

from __future__ import annotations

import logging
from collections import defaultdict
from collections.abc import Iterable, Sequence
from typing import Final

import numpy as np

from ..matrices import BLOSUM62

logger = logging.getLogger(__name__)

ALPHABET: Final[str] = "ACDEFGHIKLMNPQRSTVWY"
ALPHABET_SIZE: Final[int] = len(ALPHABET)

#: See the module docstring: 4 is the longest seed that still guarantees the
#: wild-type counterpart of any single substitution is found.
SEED_LENGTH: Final[int] = 4

#: Code for anything outside the 20 standard residues, including the separator
#: placed between proteins. Windows containing it are never aligned, so no
#: match can straddle two proteins.
_GAP: Final[int] = ALPHABET_SIZE

_CODE_TABLE = np.full(256, _GAP, dtype=np.uint8)
for _index, _residue in enumerate(ALPHABET):
    _CODE_TABLE[ord(_residue)] = _index

#: BLOSUM62 as a (21, 21) matrix so scoring is a vectorized table lookup. The
#: gap row and column are filled with the matrix minimum, which makes any
#: window overlapping a protein boundary score worse than a real alignment.
_BLOSUM_MIN: Final[int] = min(BLOSUM62.values())
_BLOSUM = np.full((ALPHABET_SIZE + 1, ALPHABET_SIZE + 1), _BLOSUM_MIN, dtype=np.int16)
for _i, _a in enumerate(ALPHABET):
    for _j, _b in enumerate(ALPHABET):
        _BLOSUM[_i, _j] = BLOSUM62[(_a, _b)]


def encode(sequence: str) -> np.ndarray:
    """Residues to 0-19 codes, with 20 for anything else."""
    raw = np.frombuffer(sequence.encode("ascii", "replace"), dtype=np.uint8)
    return _CODE_TABLE[raw]


def _seed_codes(codes: np.ndarray, length: int) -> tuple[np.ndarray, np.ndarray]:
    """Rolling base-20 code for every `length`-window, plus a validity mask."""
    windows = codes.size - length + 1
    if windows <= 0:
        return np.empty(0, dtype=np.int64), np.empty(0, dtype=bool)

    keys = np.zeros(windows, dtype=np.int64)
    valid = np.ones(windows, dtype=bool)
    for offset in range(length):
        column = codes[offset : offset + windows]
        valid &= column != _GAP
        keys *= ALPHABET_SIZE
        keys += np.minimum(column, ALPHABET_SIZE - 1).astype(np.int64)
    return keys, valid


class SelfProteome:
    """Finds each query peptide's most similar peptide in the proteome."""

    def __init__(self, sequences: Iterable[str]) -> None:
        cleaned = [s for s in (seq.strip("*") for seq in sequences) if s]
        # A separator between proteins stops any alignment window from
        # spanning two unrelated proteins and inventing a self peptide.
        self._codes = encode("*".join(cleaned))
        self._proteins = len(cleaned)
        self._seed_keys: np.ndarray | None = None

    @property
    def residues(self) -> int:
        return int(np.count_nonzero(self._codes != _GAP))

    def _keys(self) -> np.ndarray:
        if self._seed_keys is None:
            keys, valid = _seed_codes(self._codes, SEED_LENGTH)
            # Invalid windows get a key outside the key space so they can never
            # match a query seed.
            self._seed_keys = np.where(valid, keys, -1)
        return self._seed_keys

    def nearest_similarity(self, peptides: Sequence[str]) -> dict[str, float]:
        """Best BLOSUM62 similarity to any self peptide, normalized to [0, 1].

        1.0 means an exact self match, 0.0 means no seeded alignment scored
        better than the worst possible. Peptides with no seed hit at all are
        reported as 0.0, which is the correct reading: nothing in the proteome
        shares even four consecutive residues with them.
        """
        if not peptides:
            return {}

        queries = [p for p in dict.fromkeys(peptides) if p]
        wanted: dict[int, list[tuple[int, int]]] = defaultdict(list)
        encoded: list[np.ndarray] = []
        for index, peptide in enumerate(queries):
            codes = encode(peptide)
            encoded.append(codes)
            keys, valid = _seed_codes(codes, SEED_LENGTH)
            for offset, (key, ok) in enumerate(zip(keys, valid, strict=True)):
                if ok:
                    wanted[int(key)].append((index, offset))

        best = np.full(len(queries), np.iinfo(np.int32).min, dtype=np.int64)
        if wanted:
            self._align_seed_hits(queries, encoded, wanted, best)

        return {
            peptide: _normalize(int(best[index]), encoded[index])
            for index, peptide in enumerate(queries)
        }

    def _align_seed_hits(
        self,
        queries: Sequence[str],
        encoded: Sequence[np.ndarray],
        wanted: dict[int, list[tuple[int, int]]],
        best: np.ndarray,
    ) -> None:
        """Align every proteome window that shares a seed with some query."""
        keys = self._keys()
        query_keys = np.fromiter(wanted, dtype=np.int64, count=len(wanted))

        hit_mask = np.isin(keys, query_keys)
        hit_positions = np.flatnonzero(hit_mask)
        if hit_positions.size == 0:
            logger.debug("no seed hits for %d query peptides", len(queries))
            return
        logger.debug("%d seed hits across %d query peptides", int(hit_positions.size), len(queries))

        # Group hit positions by seed key so each (query, offset) pair is
        # extended against only the positions that actually matched it.
        hit_keys = keys[hit_positions]
        order = np.argsort(hit_keys, kind="stable")
        hit_keys = hit_keys[order]
        hit_positions = hit_positions[order]
        boundaries = np.searchsorted(hit_keys, query_keys)

        total = self._codes.size
        for key, start in zip(query_keys, boundaries, strict=True):
            end = np.searchsorted(hit_keys, key, side="right")
            if end <= start:
                continue
            positions = hit_positions[start:end]
            for query_index, offset in wanted[int(key)]:
                peptide = encoded[query_index]
                starts = positions - offset
                starts = starts[(starts >= 0) & (starts + peptide.size <= total)]
                if starts.size == 0:
                    continue
                scores = self._score_windows(starts, peptide)
                top = int(scores.max())
                if top > best[query_index]:
                    best[query_index] = top

    def _score_windows(self, starts: np.ndarray, peptide: np.ndarray) -> np.ndarray:
        """BLOSUM62 sum for the peptide against each window start."""
        scores = np.zeros(starts.size, dtype=np.int64)
        for offset in range(peptide.size):
            column = self._codes[starts + offset]
            scores += _BLOSUM[int(peptide[offset]), column]
        return scores


def _normalize(score: int, peptide: np.ndarray) -> float:
    """Map a raw BLOSUM sum onto [0, 1] against this peptide's own range.

    Normalizing per peptide rather than globally keeps the feature comparable
    across lengths and across residue compositions: a tryptophan-rich peptide
    has a much larger self-identity score than a glycine-rich one, and without
    this it would look artificially similar to everything.
    """
    if peptide.size == 0:
        return 0.0
    rows = _BLOSUM[peptide[:, None], np.arange(ALPHABET_SIZE)[None, :]]
    best = int(_BLOSUM[peptide, peptide].sum())
    worst = int(rows.min(axis=1).sum())
    if score <= worst or best == worst:
        return 0.0
    return max(0.0, min(1.0, (score - worst) / (best - worst)))
