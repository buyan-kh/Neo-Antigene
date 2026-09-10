"""Mutant peptide enumeration.

For each variant we build the mutant protein, then emit every window of the
requested lengths that overlaps the altered residues, paired with the
positionally matched wild-type peptide where one exists. The WT counterpart is
what makes agretopicity computable, so it is carried through rather than
recomputed later.

Supported today: missense, inframe insertion, inframe deletion.
Not yet supported: frameshift and stop-loss neo-ORFs. Those need the mutant
CDS (or VEP's Downstream plugin output), which the current input contract does
not carry — see `FrameshiftNotSupported`.
"""

from __future__ import annotations

from collections.abc import Iterator, Sequence

from ..models import AMINO_ACIDS, PeptideCandidate, Variant, VariantClass
from .proteome import ProteomeIndex


class PeptideGenerationError(ValueError):
    pass


class FrameshiftNotSupported(PeptideGenerationError):
    pass


class ReferenceMismatch(PeptideGenerationError):
    pass


class TranscriptNotFound(PeptideGenerationError):
    pass


def build_mutant_protein(variant: Variant, wildtype: str) -> tuple[str, int, int]:
    """Return (mutant_protein, altered_start, altered_end) in mutant coordinates.

    The altered interval is half-open. For pure deletions it is empty and marks
    the junction, which `_windows` widens to the flanking residues.
    """
    if variant.variant_class is VariantClass.FRAMESHIFT:
        raise FrameshiftNotSupported(
            f"{variant.hgvsp_short}: frameshift peptides require the mutant CDS"
        )

    if variant.aa_ref:
        start = variant.protein_start - 1
        end = variant.protein_end
    else:
        # Pure insertion: VEP reports the flanking residues, e.g. `-/AB` at 10-11.
        start = variant.protein_start
        end = start

    if start < 0 or end > len(wildtype):
        raise ReferenceMismatch(
            f"{variant.hgvsp_short}: protein position {variant.protein_start}-"
            f"{variant.protein_end} outside transcript of length {len(wildtype)}"
        )

    observed = wildtype[start:end]
    if variant.aa_ref and observed != variant.aa_ref:
        raise ReferenceMismatch(
            f"{variant.hgvsp_short}: reference residue {observed!r} at "
            f"{variant.protein_start} does not match annotated {variant.aa_ref!r}; "
            "proteome build and annotation are likely mismatched"
        )

    mutant = wildtype[:start] + variant.aa_alt + wildtype[end:]
    return mutant, start, start + len(variant.aa_alt)


def _windows(
    protein_length: int,
    altered_start: int,
    altered_end: int,
    length: int,
) -> Iterator[int]:
    """Yield window start offsets of `length` that overlap the altered interval.

    For a zero-length interval (a pure deletion) the window must straddle the
    junction, so it has to contain the residues on both sides of it.
    """
    if length > protein_length:
        return
    if altered_end > altered_start:
        lowest = altered_end - length
        highest = altered_start
    else:
        lowest = altered_start + 1 - length
        highest = altered_start - 1
    lowest = max(0, lowest)
    highest = min(highest, protein_length - length)
    yield from range(lowest, highest + 1)


def _is_clean(peptide: str) -> bool:
    return bool(peptide) and all(residue in AMINO_ACIDS for residue in peptide)


def _mutation_position(peptide_start: int, altered_start: int, altered_end: int) -> int:
    """1-based offset of the first altered residue within the peptide."""
    anchor = altered_start if altered_end > altered_start else max(altered_start - 1, 0)
    return anchor - peptide_start + 1


def generate_for_variant(
    variant: Variant,
    proteome: ProteomeIndex,
    lengths: Sequence[int],
    flank_length: int = 10,
) -> list[PeptideCandidate]:
    wildtype = proteome.get(variant.transcript, variant.gene)
    if wildtype is None:
        raise TranscriptNotFound(
            f"{variant.hgvsp_short}: transcript {variant.transcript} not in proteome"
        )

    mutant, altered_start, altered_end = build_mutant_protein(variant, wildtype)
    same_register = len(mutant) == len(wildtype)

    candidates: list[PeptideCandidate] = []
    seen: set[str] = set()
    for length in lengths:
        for start in _windows(len(mutant), altered_start, altered_end, length):
            mutant_peptide = mutant[start : start + length]
            if not _is_clean(mutant_peptide) or mutant_peptide in seen:
                continue
            seen.add(mutant_peptide)

            wildtype_peptide: str | None = None
            if same_register:
                counterpart = wildtype[start : start + length]
                if _is_clean(counterpart) and counterpart != mutant_peptide:
                    wildtype_peptide = counterpart

            candidates.append(
                PeptideCandidate(
                    variant=variant,
                    mutant_peptide=mutant_peptide,
                    wildtype_peptide=wildtype_peptide,
                    mutation_position=_mutation_position(start, altered_start, altered_end),
                    n_flank=_trim_n_flank(mutant[max(0, start - flank_length) : start]),
                    c_flank=_trim_c_flank(mutant[start + length : start + length + flank_length]),
                )
            )
    return candidates


def _trim_n_flank(flank: str) -> str:
    """Keep the residues adjacent to the peptide, stopping at any non-residue.

    Filtering rather than truncating would splice non-adjacent residues
    together and hand the processing predictor a sequence that does not exist.
    """
    for index in range(len(flank) - 1, -1, -1):
        if flank[index] not in AMINO_ACIDS:
            return flank[index + 1 :]
    return flank


def _trim_c_flank(flank: str) -> str:
    for index, residue in enumerate(flank):
        if residue not in AMINO_ACIDS:
            return flank[:index]
    return flank
