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

#: Prefix length used to prefilter exact self matches. Unrelated to SEED_LENGTH,
#: and chosen for the opposite reason: seeding wants short keys so near matches
#: are still found, whereas exact matching wants the longest key every candidate
#: shares so that almost nothing survives the prefilter. 8 is the minimum
#: supported peptide length, so every candidate has a prefix this long.
EXACT_PREFIX_LENGTH: Final[int] = 8

#: Residues per chunk when streaming the proteome. Keeps the rolling-key
#: intermediates in cache-friendly blocks instead of allocating one int64 array
#: per key over the whole 157 Mb reference.
_CHUNK_RESIDUES: Final[int] = 1 << 23

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


def _is_standard(peptide: str) -> bool:
    """True when every residue is one of the 20 standard amino acids.

    Anything else encodes to the same code as the inter-protein separator, so
    admitting it would let a peptide "match" self across a protein boundary.
    """
    return all(residue in ALPHABET for residue in peptide)


def _prefix_key(codes: np.ndarray, length: int) -> int:
    """Rolling key for the first `length` residues, or -1 if there is no valid one."""
    keys, valid = _seed_codes(codes[:length], length)
    if keys.size == 0 or not bool(valid[0]):
        return -1
    return int(keys[0])


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

    def exact_matches(self, peptides: Sequence[str]) -> set[str]:
        """Which of `peptides` occur verbatim in the reference proteome.

        This is the self-peptide hard gate. A candidate identical to a sequence
        the thymus already screened against is not a neoantigen, whatever its
        predicted affinity.

        Answering it one peptide at a time costs a substring scan over ~157 Mb
        per candidate, which dominates a real run. Batching inverts the loop:
        the proteome is streamed once and every candidate is tested against
        each position at the same time, so the reference is read once per run
        rather than once per candidate.

        The prefilter is what makes that cheap. An exact match requires the
        candidate's first `EXACT_PREFIX_LENGTH` residues to match, and there
        are 20^8 such prefixes, so for a query set of tens of thousands only a
        few hundred proteome positions survive and need full comparison.
        """
        clean = [p for p in dict.fromkeys(peptides) if p and _is_standard(p)]
        if not clean:
            return set()

        by_prefix_length: dict[int, list[str]] = defaultdict(list)
        for peptide in clean:
            by_prefix_length[min(EXACT_PREFIX_LENGTH, len(peptide))].append(peptide)

        found: set[str] = set()
        for prefix_length, group in by_prefix_length.items():
            found |= self._exact_with_prefix(group, prefix_length)
        return found

    def _exact_with_prefix(self, queries: Sequence[str], prefix_length: int) -> set[str]:
        encoded = {peptide: encode(peptide) for peptide in queries}
        by_prefix: dict[int, list[str]] = defaultdict(list)
        for peptide in queries:
            key = _prefix_key(encoded[peptide], prefix_length)
            if key >= 0:
                by_prefix[key].append(peptide)
        if not by_prefix:
            return set()

        prefix_keys = np.array(sorted(by_prefix), dtype=np.int64)
        total = self._codes.size
        found: set[str] = set()
        pending = {peptide for group in by_prefix.values() for peptide in group}

        start = 0
        while start + prefix_length <= total and pending:
            stop = min(start + _CHUNK_RESIDUES + prefix_length - 1, total)
            keys, valid = _seed_codes(self._codes[start:stop], prefix_length)
            keys = np.where(valid, keys, -1)

            # prefix_keys is sorted, so searchsorted finds each chunk key's
            # candidate slot without sorting the chunk itself.
            slots = np.searchsorted(prefix_keys, keys)
            np.clip(slots, 0, prefix_keys.size - 1, out=slots)
            offsets = np.flatnonzero(prefix_keys[slots] == keys)

            for offset in offsets[offsets < _CHUNK_RESIDUES].tolist():
                for peptide in by_prefix[int(keys[offset])]:
                    if peptide in pending and self._matches_at(start + offset, encoded[peptide]):
                        found.add(peptide)
                        pending.discard(peptide)
            start += _CHUNK_RESIDUES

        return found

    def _matches_at(self, position: int, codes: np.ndarray) -> bool:
        end = position + codes.size
        if end > self._codes.size:
            return False
        return bool(np.array_equal(self._codes[position:end], codes))

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
