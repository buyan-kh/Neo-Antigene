"""Mutant peptide enumeration.

For each variant we build the mutant protein, then emit every window of the
requested lengths that overlaps the altered residues, paired with the
positionally matched wild-type peptide where one exists. The WT counterpart is
what makes agretopicity computable, so it is carried through rather than
recomputed later.

Supported today: missense, inframe insertion, inframe deletion, and frameshift
neo-ORFs where the annotation carries the novel tail. Stop-loss is still
unsupported; VEP maps it to `OTHER` and it never reaches here.

Frameshifts differ from substitutions in three ways that matter here, and all
three are handled explicitly rather than by letting the substitution logic
stretch to cover them:

**The novel sequence is not in the reference.** A frameshift reads through in
a different frame, so its C-terminal tail exists in no proteome. It has to
arrive on the variant as `downstream_protein`; absent that this module raises
rather than inventing one.

**Every residue after the shift is altered.** A substitution's altered interval
is one or two residues, and a window is required to contain all of them. A
frameshift's altered interval runs to the stop codon, so that rule would admit
nothing. Windows are instead required only to *overlap* the altered interval,
which is what lets peptides lying entirely inside the novel tail be emitted —
and those are the neo-ORF epitopes. One frameshift can therefore yield many
peptides, which is correct: Roudko et al. (2020) showed single recurrent MSI
frameshifts producing several independently immunogenic epitopes.

**There is no positionally matched wild-type.** Downstream of the shift there
is nothing to compare against, so `wildtype_peptide` is None and the features
that need it go neutral. Never 1.0 — see `scoring.features.agretopicity_feature`
for why that default would have made this entire variant class look good for
free.

**Where the tail starts depends on the annotator.** A tail is joined at
`protein_start` unless the record carries an Ensembl release of 114 or later
and `ProteinLengthChange`. That pair says where the tail starts in the mutant
protein. It does not contain a residue the plugin omitted. An insertion
between codons omits a reference residue, which is copied from the proteome.
A deletion that omits the new residue is built only when `Amino_acids` states
that residue; otherwise this module raises rather than writing the reference
residue in its place. Releases 112 and 113 report a `ProteinLengthChange`
that cannot place the tail, so it is ignored there.
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


#: VEP writes the stop codon as `*` and an unresolved residue as `X`. A novel
#: tail ends at the first of either: past a stop there is no protein, and past
#: an `X` the sequence is not one we can claim to know.
_TAIL_TERMINATORS = "*X"

#: Ensembl release at which `Downstream.pm` started reporting
#: `ProteinLengthChange` as the full mutant peptide length minus the reference
#: length. That definition, with the tail, is the join. The codon-offset
#: change that can drop a residue landed in release 112, but 112 and 113 kept
#: the older length-change definition, which assumes the tail already starts
#: at `protein_start` and so cannot correct the drop.
DOWNSTREAM_LENGTH_CHANGE_RELEASE = 114


def build_mutant_protein(variant: Variant, wildtype: str) -> tuple[str, int, int]:
    """Return (mutant_protein, altered_start, altered_end) in mutant coordinates.

    The altered interval is half-open. For pure deletions it is empty and marks
    the junction, which `_windows` widens to the flanking residues. For
    frameshifts it runs from the shift to the end of the mutant protein.
    """
    if variant.variant_class is VariantClass.FRAMESHIFT:
        return _build_frameshift_protein(variant, wildtype)

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


def _build_frameshift_protein(variant: Variant, wildtype: str) -> tuple[str, int, int]:
    tail = _truncate_tail(variant.downstream_protein or "")
    if not tail:
        raise FrameshiftNotSupported(
            f"{variant.hgvsp_short}: frameshift peptides need the novel tail, which the "
            f"reference proteome cannot supply. Annotate with VEP's Downstream plugin and "
            f"carry DownstreamProtein, or supply the mutant protein from pVACtools' "
            f"Frameshift.pm sliced at protein_start"
        )

    mutant, start = _assemble_frameshift(variant, wildtype, tail)
    return mutant, start, len(mutant)


def _assemble_frameshift(variant: Variant, wildtype: str, tail: str) -> tuple[str, int]:
    """Return (mutant protein, 0-based index of the first altered residue).

    Without a usable length change the tail is taken to begin at
    `protein_start`. That matches pVACtools sliced at that position, and
    Downstream.pm through Ensembl 111.

    From Ensembl 114, `ProteinLengthChange` is the mutant length minus the
    reference length, so it says which mutant residue the tail starts at.
    When that is later than `protein_start`, the missing residue is not in
    the tail. Copying it from the reference is right for an insertion
    between codons, where the skipped residue was never changed, and wrong
    for a deletion, where the skipped residue is a new amino acid. A deletion
    is built only from the residue `Amino_acids` actually states.
    """
    annotated = variant.protein_start - 1
    if annotated < 0 or annotated > len(wildtype):
        raise ReferenceMismatch(
            f"{variant.hgvsp_short}: frameshift at protein position "
            f"{variant.protein_start} is outside transcript of length {len(wildtype)}"
        )

    tail_start = _tail_start(variant, wildtype, tail)
    if tail_start is None or tail_start == variant.protein_start:
        _reject_reference_tail(variant, wildtype, annotated, tail)
        return wildtype[:annotated] + tail, annotated

    if tail_start < variant.protein_start or tail_start > len(wildtype) + 1:
        raise ReferenceMismatch(
            f"{variant.hgvsp_short}: ProteinLengthChange {variant.protein_length_change} "
            f"places the novel tail at residue {tail_start}, which does not follow "
            f"protein position {variant.protein_start}"
        )

    gap = tail_start - variant.protein_start
    if not variant.aa_ref:
        # Nothing in the reference was replaced. The residues the tail
        # started after are still the reference residues.
        _reject_reference_tail(variant, wildtype, tail_start - 1, tail)
        return wildtype[: tail_start - 1] + tail, tail_start - 1

    bridge = _known_variant_residues(variant.aa_alt)
    observed = wildtype[annotated : annotated + len(variant.aa_ref)]
    if observed != variant.aa_ref:
        raise ReferenceMismatch(
            f"{variant.hgvsp_short}: reference residue {observed!r} at "
            f"{variant.protein_start} does not match annotated {variant.aa_ref!r}; "
            "proteome build and annotation are likely mismatched"
        )
    if len(bridge) != gap:
        raise FrameshiftNotSupported(
            f"{variant.hgvsp_short}: DownstreamProtein starts at residue {tail_start}, "
            f"so the {gap} residue(s) at position {variant.protein_start} are absent "
            f"from the tail. They are not the reference, and the amino-acid annotation "
            f"{variant.aa_alt!r} does not state them. Supply pVACtools FrameshiftSequence, "
            f"which includes the altered residue, or a DownstreamProtein that does."
        )
    return wildtype[:annotated] + bridge + tail, annotated


def _tail_start(variant: Variant, wildtype: str, tail: str) -> int | None:
    """1-based mutant index where `tail` begins, when the release makes it knowable."""
    release = variant.vep_release
    change = variant.protein_length_change
    if release is None or change is None or release < DOWNSTREAM_LENGTH_CHANGE_RELEASE:
        return None
    # length(prefix) + length(tail) = length(mutant), and ProteinLengthChange
    # is length(mutant) - length(reference).
    tail_start = len(wildtype) + change - len(tail) + 1
    if tail_start < 1 or tail_start > len(wildtype) + 1:
        raise ReferenceMismatch(
            f"{variant.hgvsp_short}: ProteinLengthChange {change} places the novel "
            f"tail at residue {tail_start}, outside a transcript of length {len(wildtype)}"
        )
    return tail_start


def _known_variant_residues(amino_acids: str) -> str:
    """Standard residues before the first unknown (`X`) or stop (`*`)."""
    known: list[str] = []
    for residue in amino_acids:
        if residue not in AMINO_ACIDS:
            break
        known.append(residue)
    return "".join(known)


def _reject_reference_tail(variant: Variant, wildtype: str, start: int, tail: str) -> None:
    if wildtype[start : start + len(tail)] == tail:
        raise ReferenceMismatch(
            f"{variant.hgvsp_short}: the annotated downstream protein is identical to the "
            f"reference from position {start + 1}, so it encodes no frameshift. "
            f"The annotation and the proteome build are likely mismatched"
        )


def _truncate_tail(tail: str) -> str:
    for index, residue in enumerate(tail):
        if residue in _TAIL_TERMINATORS:
            return tail[:index]
    return tail


def _windows(
    protein_length: int,
    altered_start: int,
    altered_end: int,
    length: int,
    overlap_only: bool = False,
) -> Iterator[int]:
    """Yield window start offsets of `length` that cover the altered interval.

    By default a window must contain the whole altered interval, which is the
    right rule for a substitution or a short indel: a peptide carrying only
    part of an insertion is a different hypothesis from one carrying all of it.

    `overlap_only` relaxes that to intersecting the interval, which is required
    once the interval is a frameshift's novel tail. Demanding containment of a
    hundred-residue tail inside a 9-mer admits nothing, and the peptides lying
    wholly inside the tail are exactly the neo-ORF epitopes.

    For a zero-length interval (a pure deletion) the window must straddle the
    junction, so it has to contain the residues on both sides of it.
    """
    if length > protein_length:
        return
    if altered_end > altered_start:
        if overlap_only:
            lowest = altered_start - length + 1
            highest = altered_end - 1
        else:
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
    """1-based offset of the first altered residue within the peptide.

    Clamped to 1 because a peptide can lie entirely downstream of a frameshift,
    in which case the first altered residue precedes the peptide and every
    residue in it is novel.
    """
    anchor = altered_start if altered_end > altered_start else max(altered_start - 1, 0)
    return max(1, anchor - peptide_start + 1)


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
    frameshift = variant.variant_class is VariantClass.FRAMESHIFT
    # A frameshift changes the reading frame, so even when the mutant protein
    # happens to come out the same length there is no positional counterpart
    # downstream of the shift. Comparing lengths alone would occasionally pair a
    # neo-ORF peptide with an unrelated wild-type one.
    same_register = not frameshift and len(mutant) == len(wildtype)

    candidates: list[PeptideCandidate] = []
    seen: set[str] = set()
    for length in lengths:
        for start in _windows(
            len(mutant), altered_start, altered_end, length, overlap_only=frameshift
        ):
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
