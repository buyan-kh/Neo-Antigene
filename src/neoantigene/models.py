"""Domain models exchanged between pipeline stages.

Every object that crosses a stage boundary is a frozen pydantic model with its
invariants enforced at construction. Two reasons this is worth the ceremony:

  - A peptide with a non-residue character, an affinity percentile of 400, or
    a CCF of 1.7 is a bug that would otherwise surface as a silently wrong
    shortlist. Failing at the boundary makes it a stack trace instead.
  - Frozen models mean stages cannot mutate their inputs, so a run is a
    composition of pure functions and is reproducible by construction.

Stages return new objects rather than editing in place; use `model_copy` for
derived values (see `variants.clonality.with_ccf`).
"""

from __future__ import annotations

from enum import StrEnum
from typing import Annotated, Self

from pydantic import BaseModel, ConfigDict, Field, field_validator, model_validator

from .hla import Allele

AMINO_ACIDS = frozenset("ACDEFGHIKLMNPQRSTVWY")

#: VEP legitimately emits `*` for stop and `X` for an unknown residue in
#: `Amino_acids`. They are accepted on variants and rejected on peptides, so a
#: variant carrying one is kept for reporting but can never yield a candidate.
ANNOTATION_RESIDUES = AMINO_ACIDS | frozenset("*X")

#: Peptide lengths a class I predictor will accept.
MIN_PEPTIDE_LENGTH = 8
MAX_PEPTIDE_LENGTH = 15

Fraction = Annotated[float, Field(ge=0.0, le=1.0)]
Percentile = Annotated[float, Field(ge=0.0, le=100.0)]
NonNegativeFloat = Annotated[float, Field(ge=0.0)]


def _validate_residues(sequence: str, field: str, alphabet: frozenset[str] = AMINO_ACIDS) -> str:
    invalid = sorted(set(sequence) - alphabet)
    if invalid:
        raise ValueError(f"{field} contains non-residue characters {invalid}: {sequence!r}")
    return sequence


class VariantClass(StrEnum):
    MISSENSE = "missense"
    INFRAME_INS = "inframe_insertion"
    INFRAME_DEL = "inframe_deletion"
    FRAMESHIFT = "frameshift"
    OTHER = "other"


class Frozen(BaseModel):
    model_config = ConfigDict(frozen=True, extra="forbid", validate_assignment=True)


class Variant(Frozen):
    """A protein-altering somatic variant on one transcript."""

    chrom: str
    pos: int = Field(gt=0)
    ref: str
    alt: str
    gene: str
    transcript: str
    variant_class: VariantClass
    protein_start: int = Field(gt=0)
    protein_end: int = Field(gt=0)
    aa_ref: str = ""
    aa_alt: str = ""
    dna_vaf: Fraction | None = None
    rna_vaf: Fraction | None = None
    tumor_depth: int | None = Field(default=None, ge=0)
    copy_number: NonNegativeFloat | None = None
    population_af: Fraction | None = None
    filters: tuple[str, ...] = ()
    ccf: Fraction | None = None

    @field_validator("aa_ref", "aa_alt")
    @classmethod
    def _residues_only(cls, value: str) -> str:
        return _validate_residues(value, "amino acid", ANNOTATION_RESIDUES)

    @model_validator(mode="after")
    def _coherent_protein_span(self) -> Self:
        if self.protein_end < self.protein_start:
            raise ValueError(
                f"protein_end {self.protein_end} precedes protein_start {self.protein_start}"
            )
        return self

    @property
    def key(self) -> str:
        return f"{self.chrom}:{self.pos}:{self.ref}>{self.alt}"

    @property
    def span(self) -> str:
        if self.protein_start == self.protein_end:
            return str(self.protein_start)
        return f"{self.protein_start}_{self.protein_end}"

    @property
    def hgvsp_short(self) -> str:
        if self.aa_ref and self.aa_alt:
            return f"{self.gene} p.{self.aa_ref}{self.span}{self.aa_alt}"
        if self.aa_ref:
            return f"{self.gene} p.{self.aa_ref}{self.span}del"
        return f"{self.gene} p.{self.span}ins{self.aa_alt}"


class PeptideCandidate(Frozen):
    """A single mutant peptide in a single register, before HLA assignment."""

    variant: Variant
    mutant_peptide: str = Field(min_length=MIN_PEPTIDE_LENGTH, max_length=MAX_PEPTIDE_LENGTH)
    wildtype_peptide: str | None = None
    mutation_position: int = Field(gt=0)
    n_flank: str = ""
    c_flank: str = ""

    @field_validator("mutant_peptide", "wildtype_peptide", "n_flank", "c_flank")
    @classmethod
    def _residues_only(cls, value: str | None) -> str | None:
        return None if value is None else _validate_residues(value, "peptide")

    @model_validator(mode="after")
    def _consistent_register(self) -> Self:
        if self.mutation_position > len(self.mutant_peptide):
            raise ValueError(
                f"mutation_position {self.mutation_position} is outside peptide "
                f"of length {len(self.mutant_peptide)}"
            )
        if self.wildtype_peptide is not None:
            if len(self.wildtype_peptide) != len(self.mutant_peptide):
                raise ValueError(
                    "wildtype_peptide must share the register of mutant_peptide; "
                    f"got {len(self.wildtype_peptide)} vs {len(self.mutant_peptide)}"
                )
            if self.wildtype_peptide == self.mutant_peptide:
                raise ValueError("wildtype_peptide is identical to mutant_peptide")
        return self

    @property
    def length(self) -> int:
        return len(self.mutant_peptide)

    @property
    def id(self) -> str:
        return f"{self.variant.key}|{self.variant.transcript}|{self.mutant_peptide}"


class PresentationCall(Frozen):
    """One predictor output for a peptide-allele pair."""

    peptide: str
    allele: Allele
    presentation_score: Fraction | None = None
    affinity_nm: NonNegativeFloat | None = None
    affinity_percentile: Percentile | None = None
    processing_score: Fraction | None = None

    @property
    def key(self) -> tuple[str, str]:
        return self.peptide, self.allele


class GateFailure(StrEnum):
    """Reasons a candidate is excluded outright rather than scored down."""

    WEAK_BINDER = "weak_binder"
    LOW_PRESENTATION = "low_presentation"
    SELF_PEPTIDE = "self_peptide"
    NORMAL_TISSUE_EXPRESSED = "normal_tissue_expressed"


class ScoredCandidate(Frozen):
    """A peptide-allele pair with all ranking features resolved."""

    candidate: PeptideCandidate
    allele: Allele
    mutant_call: PresentationCall
    wildtype_call: PresentationCall | None = None
    features: dict[str, float] = Field(default_factory=dict)
    score: Fraction = 0.0
    gate_failures: tuple[GateFailure, ...] = ()

    @property
    def passed(self) -> bool:
        return not self.gate_failures

    @property
    def key(self) -> str:
        return f"{self.candidate.mutant_peptide}|{self.allele}"

    def with_score(self, score: float) -> ScoredCandidate:
        return self.model_copy(update={"score": score})
