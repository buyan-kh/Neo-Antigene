"""Sequence-intrinsic immunogenicity terms.

These are proxies, not models. They exist so the ranking has something to say
about "will a T cell see this" beyond "will it be presented", and every one of
them is expected to be replaced by coefficients fit on real assay outcomes.
"""

from __future__ import annotations

from ..matrices import BLOSUM62, KYTE_DOOLITTLE

BLOSUM_MATCH_SCORE = 4.0
BLOSUM_RANGE = 8.0


def anchor_positions(length: int) -> set[int]:
    """1-based canonical MHC class I anchor positions."""
    if length < 8:
        return set()
    return {2, length}


def tcr_contact_positions(length: int) -> list[int]:
    """1-based positions that face the T-cell receptor rather than the groove."""
    anchors = anchor_positions(length)
    return [p for p in range(1, length + 1) if p not in anchors and p != 1]


#: P1 sits partly in the A pocket and partly under the TCR, so a mutation there
#: is neither clearly buried nor clearly visible. It gets the neutral value
#: rather than being forced into one of the two clean cases.
_AMBIGUOUS_POSITION = 1
_ANCHOR_EXPOSURE = 0.0
_AMBIGUOUS_EXPOSURE = 0.5
_TCR_EXPOSURE = 1.0


def mutation_exposure(position: int, length: int) -> float:
    """Whether the mutated residue points at the TCR or into the MHC groove.

    Returned on [0, 1], higher meaning more TCR-facing.

    Read the weight on this term before trusting it. The mechanistic story is
    that a substitution at an anchor (P2 and the C-terminus for class I)
    changes *whether* the peptide is presented while the surface a T cell sees
    stays the self surface, whereas a substitution at a solvent-exposed
    position changes that surface directly.

    The literature does not support the directional claim. Capietto et al.
    (J Exp Med 2020) establish that mutation position matters, but what they
    show is that position determines *which* affinity metric predicts
    immunogenicity — absolute affinity for non-anchor mutations, affinity
    relative to wild-type for anchor mutations — and they state explicitly that
    "both anchor and nonanchor mutated peptides contain cases that show CD8
    responses". TESLA (Wells et al., Cell 2020) went further and found mutation
    position was not useful for filtering at all.

    So this is kept near zero, computed and reported for a future refit rather
    than trusted to order a shortlist today.
    """
    if length <= 0 or not 1 <= position <= length:
        return _AMBIGUOUS_EXPOSURE
    if position in anchor_positions(length):
        return _ANCHOR_EXPOSURE
    if position == _AMBIGUOUS_POSITION:
        return _AMBIGUOUS_EXPOSURE
    return _TCR_EXPOSURE


def tcr_contact_hydrophobicity(peptide: str) -> float:
    """Mean Kyte-Doolittle score over TCR-facing residues, scaled to [0, 1].

    Shipped with a weight of zero, deliberately. The two best sources disagree
    on the sign of this effect, and they disagree in the setting that matters
    here. Chowell et al. (PNAS 2015) found "a strong bias toward hydrophobic
    amino acids at T-cell receptor contact residues within immunogenic
    epitopes", validated in vivo — but on viral and self epitopes. TESLA
    (Wells et al., Cell 2020), working on human tumor neoepitopes, found
    immunogenic pMHC were significantly *less* hydrophobic (p = 0.04) and that
    hydrophobicity was not useful for filtering.

    Guessing a sign here would be inventing a result. The value is still
    computed and written to `features.json` so `neoantigene refit` can settle
    it from real labels, which is the only thing that will.
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
