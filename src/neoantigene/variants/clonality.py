"""Cancer cell fraction estimation.

A subclonal neoantigen is present in only part of the tumor, so even a
perfectly immunogenic peptide leaves an escape population behind. CCF is
therefore a first-class ranking term rather than a QC annotation.
"""

from __future__ import annotations

from ..models import Variant

NORMAL_COPY_NUMBER = 2.0
DEFAULT_TUMOR_COPY_NUMBER = 2.0


def cancer_cell_fraction(
    vaf: float | None,
    purity: float,
    tumor_copy_number: float | None = None,
    multiplicity: float = 1.0,
) -> float | None:
    """Purity- and ploidy-corrected CCF.

    CCF = VAF * (purity * CN_tumor + 2 * (1 - purity)) / (purity * multiplicity)

    `multiplicity` is the number of tumor chromosome copies carrying the
    variant; without an allele-specific CN caller it defaults to 1, which is
    the conservative choice (it cannot inflate CCF).
    """
    if vaf is None or vaf <= 0:
        return None
    if not 0 < purity <= 1:
        raise ValueError(f"tumor purity must be in (0, 1], got {purity}")
    if multiplicity <= 0:
        raise ValueError(f"multiplicity must be positive, got {multiplicity}")
    copy_number = (
        tumor_copy_number
        if tumor_copy_number is not None and tumor_copy_number > 0
        else DEFAULT_TUMOR_COPY_NUMBER
    )
    observable = purity * copy_number + NORMAL_COPY_NUMBER * (1 - purity)
    ccf = vaf * observable / (purity * multiplicity)
    return max(0.0, min(1.0, ccf))


def with_ccf(variant: Variant, purity: float) -> Variant:
    """Return a copy of `variant` carrying its estimated CCF."""
    return variant.model_copy(
        update={
            "ccf": cancer_cell_fraction(
                vaf=variant.dna_vaf,
                purity=purity,
                tumor_copy_number=variant.copy_number,
            )
        }
    )
