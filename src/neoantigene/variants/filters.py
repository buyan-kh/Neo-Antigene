"""Variant-level gating.

Gates are hard drops applied before peptide enumeration. Anything that is a
matter of degree belongs in scoring, not here.

`apply_variant_filters` does not mutate its input: it returns annotated copies
of the variants that survived, plus a report of what was dropped and why.
"""

from __future__ import annotations

from collections.abc import Sequence
from enum import StrEnum

from pydantic import BaseModel, Field

from ..config import ExpressionConfig, VariantFilterConfig
from ..io.expression import ExpressionTable
from ..models import Variant
from .clonality import with_ccf


class DropReason(StrEnum):
    NOT_PASS = "not_pass"
    LOW_DEPTH = "low_depth"
    LOW_VAF = "low_vaf"
    GERMLINE_COMMON = "germline_common"
    SUBCLONAL = "subclonal"
    NO_EXPRESSION_DATA = "no_expression_data"
    NOT_EXPRESSED = "not_expressed"
    ALLELE_NOT_EXPRESSED = "allele_not_expressed"


class FilterReport(BaseModel):
    kept: int = 0
    dropped: dict[DropReason, int] = Field(default_factory=dict)

    @property
    def total_dropped(self) -> int:
        return sum(self.dropped.values())

    def summary(self) -> str:
        if not self.dropped:
            return f"kept {self.kept}, dropped 0"
        detail = ", ".join(f"{k.value}={v}" for k, v in sorted(self.dropped.items()))
        return f"kept {self.kept}, dropped {self.total_dropped} ({detail})"


def _drop_reason(
    variant: Variant,
    config: VariantFilterConfig,
    expression_config: ExpressionConfig,
    expression: ExpressionTable | None,
) -> DropReason | None:
    if config.require_pass and variant.filters and variant.filters != ("PASS",):
        return DropReason.NOT_PASS
    if variant.tumor_depth is not None and variant.tumor_depth < config.min_tumor_depth:
        return DropReason.LOW_DEPTH
    if variant.dna_vaf is not None and variant.dna_vaf < config.min_dna_vaf:
        return DropReason.LOW_VAF
    if variant.population_af is not None and variant.population_af > config.max_population_af:
        return DropReason.GERMLINE_COMMON
    if variant.ccf is not None and variant.ccf < config.min_ccf:
        return DropReason.SUBCLONAL

    if expression is not None:
        tpm = expression.lookup(variant.transcript, variant.gene)
        if tpm is None:
            if expression_config.require_rna_evidence:
                return DropReason.NO_EXPRESSION_DATA
        elif tpm < expression_config.min_tpm:
            return DropReason.NOT_EXPRESSED

    if (
        expression_config.min_rna_vaf is not None
        and variant.rna_vaf is not None
        and variant.rna_vaf < expression_config.min_rna_vaf
    ):
        return DropReason.ALLELE_NOT_EXPRESSED

    return None


def apply_variant_filters(
    variants: Sequence[Variant],
    config: VariantFilterConfig,
    expression_config: ExpressionConfig,
    expression: ExpressionTable | None = None,
) -> tuple[list[Variant], FilterReport]:
    report = FilterReport()
    kept: list[Variant] = []

    for variant in variants:
        annotated = with_ccf(variant, config.tumor_purity)
        reason = _drop_reason(annotated, config, expression_config, expression)
        if reason is not None:
            report.dropped[reason] = report.dropped.get(reason, 0) + 1
            continue
        report.kept += 1
        kept.append(annotated)

    return kept, report
