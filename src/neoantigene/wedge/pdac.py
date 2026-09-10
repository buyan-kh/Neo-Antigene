"""PDAC wedge: recurrent pancreatic adenocarcinoma variants.

Shared (public) neoantigens are the entry point for a self-funded assay loop:
the same peptide can be tested against many donors, so ground truth
accumulates faster than it would with private neoantigens.

This catalog is a screening aid, not a ranking shortcut. Membership here does
not add score — the ranker has no "known driver" term by design. Its uses are
(a) checking whether a sample carries a variant the lab already has peptide
for, and (b) prioritising which shared antigens to stock.

Frequencies are among KRAS-mutant PDAC (KRAS itself is mutated in ~90%),
from AACR GENIE-based and clinicogenomic cohort summaries; they vary by
cohort and assay.
"""

from __future__ import annotations

from collections.abc import Sequence

from pydantic import BaseModel, ConfigDict, Field

from ..hla import Allele
from ..models import Variant


class SharedAntigen(BaseModel):
    model_config = ConfigDict(frozen=True)

    gene: str
    protein_change: str
    frequency: float | None = Field(default=None, ge=0.0, le=1.0)
    restricting_alleles: dict[Allele, tuple[str, ...]] = Field(default_factory=dict)
    note: str = ""

    @property
    def key(self) -> str:
        return f"{self.gene} {self.protein_change}"


KRAS_HOTSPOTS: tuple[SharedAntigen, ...] = (
    SharedAntigen(
        gene="KRAS",
        protein_change="G12D",
        frequency=0.41,
        restricting_alleles={
            "HLA-C*08:02": ("GADGVGKSA", "GADGVGKSAL"),
            "HLA-A*11:01": ("VVGADGVGK", "VVVGADGVGK"),
        },
        note=(
            "C*08:02 epitopes are mutant-specific: D12 forms an anchoring salt "
            "bridge that the wild-type glycine cannot, so the WT counterpart is "
            "not presented and the repertoire is untolerized."
        ),
    ),
    SharedAntigen(
        gene="KRAS",
        protein_change="G12V",
        frequency=0.32,
        restricting_alleles={"HLA-A*11:01": ("VVGAVGVGK", "VVVGAVGVGK")},
        note=(
            "A*11:01 G12V TCRs have been reported to cross-recognize a RAB7B "
            "self peptide; screen candidate TCRs against it before any "
            "engineering work."
        ),
    ),
    SharedAntigen(gene="KRAS", protein_change="G12R", frequency=0.16),
    SharedAntigen(gene="KRAS", protein_change="Q61H", frequency=0.05),
    SharedAntigen(gene="KRAS", protein_change="G12C", frequency=0.02),
)

OTHER_RECURRENT: tuple[SharedAntigen, ...] = (
    SharedAntigen(
        gene="TP53",
        protein_change="R175H",
        restricting_alleles={"HLA-A*02:01": ("HMTEVVRHC",)},
        note="Recurrent across tumor types; PDAC TP53 mutations are mostly private.",
    ),
    SharedAntigen(gene="SMAD4", protein_change="R361H"),
    SharedAntigen(gene="CDKN2A", protein_change="R58*"),
)

CATALOG: tuple[SharedAntigen, ...] = KRAS_HOTSPOTS + OTHER_RECURRENT

_BY_KEY: dict[str, SharedAntigen] = {antigen.key: antigen for antigen in CATALOG}


def lookup(gene: str, protein_change: str) -> SharedAntigen | None:
    return _BY_KEY.get(f"{gene} {protein_change}")


def protein_change_of(variant: Variant) -> str:
    return f"{variant.aa_ref}{variant.protein_start}{variant.aa_alt}"


def match_variants(variants: Sequence[Variant]) -> list[SharedAntigen]:
    """Return catalog entries covered by a sample's variants, in catalog order."""
    matched = {
        antigen.key
        for variant in variants
        if (antigen := lookup(variant.gene, protein_change_of(variant))) is not None
    }
    return [a for a in CATALOG if a.key in matched]


def eligible_by_hla(alleles: Sequence[str]) -> list[SharedAntigen]:
    """Catalog entries with a documented epitope for one of these alleles.

    Absence here means "no published restriction", not "will not work" — the
    ranker still evaluates these variants against the full genotype.
    """
    allele_set = set(alleles)
    return [a for a in CATALOG if allele_set & set(a.restricting_alleles)]
