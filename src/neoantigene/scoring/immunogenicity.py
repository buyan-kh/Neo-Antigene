"""Sequence-intrinsic immunogenicity terms.

These are proxies, not models. They exist so the ranking has something to say
about "will a T cell see this" beyond "will it be presented", and every one of
them is expected to be replaced by coefficients fit on real assay outcomes.
"""

from __future__ import annotations

BLOSUM_MATCH_SCORE = 4.0
BLOSUM_RANGE = 8.0

KYTE_DOOLITTLE = {
    "A": 1.8,
    "R": -4.5,
    "N": -3.5,
    "D": -3.5,
    "C": 2.5,
    "Q": -3.5,
    "E": -3.5,
    "G": -0.4,
    "H": -3.2,
    "I": 4.5,
    "L": 3.8,
    "K": -3.9,
    "M": 1.9,
    "F": 2.8,
    "P": -1.6,
    "S": -0.8,
    "T": -0.7,
    "W": -0.9,
    "Y": -1.3,
    "V": 4.2,
}

_BLOSUM62_ORDER = "ARNDCQEGHILKMFPSTWYV"
_BLOSUM62_ROWS = [
    [4, -1, -2, -2, 0, -1, -1, 0, -2, -1, -1, -1, -1, -2, -1, 1, 0, -3, -2, 0],
    [-1, 5, 0, -2, -3, 1, 0, -2, 0, -3, -2, 2, -1, -3, -2, -1, -1, -3, -2, -3],
    [-2, 0, 6, 1, -3, 0, 0, 0, 1, -3, -3, 0, -2, -3, -2, 1, 0, -4, -2, -3],
    [-2, -2, 1, 6, -3, 0, 2, -1, -1, -3, -4, -1, -3, -3, -1, 0, -1, -4, -3, -3],
    [0, -3, -3, -3, 9, -3, -4, -3, -3, -1, -1, -3, -1, -2, -3, -1, -1, -2, -2, -1],
    [-1, 1, 0, 0, -3, 5, 2, -2, 0, -3, -2, 1, 0, -3, -1, 0, -1, -2, -1, -2],
    [-1, 0, 0, 2, -4, 2, 5, -2, 0, -3, -3, 1, -2, -3, -1, 0, -1, -3, -2, -2],
    [0, -2, 0, -1, -3, -2, -2, 6, -2, -4, -4, -2, -3, -3, -2, 0, -2, -2, -3, -3],
    [-2, 0, 1, -1, -3, 0, 0, -2, 8, -3, -3, -1, -2, -1, -2, -1, -2, -2, 2, -3],
    [-1, -3, -3, -3, -1, -3, -3, -4, -3, 4, 2, -3, 1, 0, -3, -2, -1, -3, -1, 3],
    [-1, -2, -3, -4, -1, -2, -3, -4, -3, 2, 4, -2, 2, 0, -3, -2, -1, -2, -1, 1],
    [-1, 2, 0, -1, -3, 1, 1, -2, -1, -3, -2, 5, -1, -3, -1, 0, -1, -3, -2, -2],
    [-1, -1, -2, -3, -1, 0, -2, -3, -2, 1, 2, -1, 5, 0, -2, -1, -1, -1, -1, 1],
    [-2, -3, -3, -3, -2, -3, -3, -3, -1, 0, 0, -3, 0, 6, -4, -2, -2, 1, 3, -1],
    [-1, -2, -2, -1, -3, -1, -1, -2, -2, -3, -3, -1, -2, -4, 7, -1, -1, -4, -3, -2],
    [1, -1, 1, 0, -1, 0, 0, 0, -1, -2, -2, 0, -1, -2, -1, 4, 1, -3, -2, -2],
    [0, -1, 0, -1, -1, -1, -1, -2, -2, -1, -1, -1, -1, -2, -1, 1, 5, -2, -2, 0],
    [-3, -3, -4, -4, -2, -2, -3, -2, -2, -3, -2, -3, -1, 1, -4, -3, -2, 11, 2, -3],
    [-2, -2, -2, -3, -2, -1, -2, -3, 2, -1, -1, -2, -1, 3, -3, -2, -2, 2, 7, -1],
    [0, -3, -3, -3, -1, -2, -2, -3, -3, 3, 1, -2, 1, -1, -2, -2, 0, -3, -1, 4],
]
BLOSUM62 = {
    (a, b): _BLOSUM62_ROWS[i][j]
    for i, a in enumerate(_BLOSUM62_ORDER)
    for j, b in enumerate(_BLOSUM62_ORDER)
}


def anchor_positions(length: int) -> set[int]:
    """1-based canonical MHC class I anchor positions."""
    if length < 8:
        return set()
    return {2, length}


def tcr_contact_positions(length: int) -> list[int]:
    """1-based positions that face the T-cell receptor rather than the groove."""
    anchors = anchor_positions(length)
    return [p for p in range(1, length + 1) if p not in anchors and p != 1]


def tcr_contact_hydrophobicity(peptide: str) -> float:
    """Mean Kyte-Doolittle score over TCR-facing residues, scaled to [0, 1].

    Hydrophobic TCR-contact residues are one of the few sequence features that
    has held up across immunogenicity datasets.
    """
    positions = tcr_contact_positions(len(peptide))
    values = [KYTE_DOOLITTLE.get(peptide[p - 1], 0.0) for p in positions]
    if not values:
        return 0.5
    mean = sum(values) / len(values)
    return _clip((mean + 4.5) / 9.0)


def agretopicity(
    mutant_affinity_nm: float | None, wildtype_affinity_nm: float | None
) -> float | None:
    """Differential agretopicity index: WT affinity / mutant affinity.

    A high ratio means the mutation created the binding event, so the WT
    counterpart was never presented and the T-cell repertoire was never
    tolerized against it.
    """
    if not mutant_affinity_nm or mutant_affinity_nm <= 0:
        return None
    if not wildtype_affinity_nm or wildtype_affinity_nm <= 0:
        return None
    return wildtype_affinity_nm / mutant_affinity_nm


def dissimilarity_to_wildtype(mutant: str, wildtype: str | None) -> float:
    """Mean BLOSUM62 dissimilarity over substituted TCR-facing positions, in [0, 1].

    A conservative substitution buried against the MHC groove is invisible to
    a T cell; a radical one at P5 is the whole point. Averaging over only the
    substituted positions keeps a single radical change from being diluted by
    the unchanged residues around it.
    """
    if wildtype is None or len(mutant) != len(wildtype):
        # No positional WT counterpart exists (indel / neo-ORF): the surface is
        # novel by construction.
        return 1.0
    positions = tcr_contact_positions(len(mutant))
    substituted = [p for p in positions if mutant[p - 1] != wildtype[p - 1]]
    if not substituted:
        return 0.0
    total = 0.0
    for position in substituted:
        a, b = mutant[position - 1], wildtype[position - 1]
        score = BLOSUM62.get((a, b))
        total += 1.0 if score is None else _clip((BLOSUM_MATCH_SCORE - score) / BLOSUM_RANGE)
    return _clip(total / len(substituted))


def _clip(value: float, low: float = 0.0, high: float = 1.0) -> float:
    return max(low, min(high, value))
