"""Feature assembly.

Every scored feature is normalized to [0, 1] and oriented so that higher is
better, which is what makes the weight vector directly interpretable and lets
`learning.train` refit it as plain logistic regression on assay labels.

Raw (unnormalized) values are kept alongside under `*_raw` keys for reporting
and for refitting on a different transform later. They are not scored.
"""

from __future__ import annotations

import math
from dataclasses import dataclass

from ..config import PipelineConfig
from ..io.expression import ExpressionTable, NormalExpressionReference
from ..models import GateFailure, PeptideCandidate, PresentationCall
from ..peptides.proteome import ProteomeIndex
from .immunogenicity import (
    agretopicity,
    dissimilarity_to_wildtype,
    mutation_exposure,
    tcr_contact_hydrophobicity,
)

NEUTRAL = 0.5

EXPRESSION_LOG_CEILING = 2.5
NORMAL_EXPRESSION_LOG_CEILING = 2.5
AGRETOPICITY_LOG_CEILING = 2.0


@dataclass(frozen=True)
class FeatureContext:
    """Everything feature assembly needs that is not on the candidate itself."""

    config: PipelineConfig
    expression: ExpressionTable | None = None
    normal_expression: NormalExpressionReference | None = None
    proteome: ProteomeIndex | None = None

    #: Nearest-self BLOSUM62 similarity per mutant peptide, precomputed in one
    #: proteome pass by `peptides.selfsim`. `None` means the search could not be
    #: run (no reference, or a stub too small to mean anything), in which case
    #: the feature is neutral rather than optimistic.
    self_similarity: dict[str, float] | None = None

    #: Mutant peptides found verbatim in the reference proteome, precomputed in
    #: one batched pass by `SelfProteome.exact_matches`. `None` falls back to
    #: scanning the proteome per candidate, which is correct but costs a full
    #: substring search each time — fine for a single `explain`, ruinous for a
    #: whole run.
    self_matches: frozenset[str] | None = None


def expression_feature(tpm: float | None) -> float:
    if tpm is None:
        return NEUTRAL
    return _clip(math.log10(max(tpm, 0.0) + 1.0) / EXPRESSION_LOG_CEILING)


def tumor_selectivity_feature(normal_tpm: float | None) -> float:
    """Higher when the source transcript is quiet in healthy tissue.

    This is the on-target/off-tumor term. It is a penalty on the gene, not on
    the peptide; peptide-level self-identity is handled as a hard gate.
    """
    if normal_tpm is None:
        return NEUTRAL
    burden = _clip(math.log10(max(normal_tpm, 0.0) + 1.0) / NORMAL_EXPRESSION_LOG_CEILING)
    return 1.0 - burden


def presentation_feature(call: PresentationCall) -> float:
    if call.presentation_score is not None:
        return _clip(call.presentation_score)
    if call.affinity_percentile is not None:
        return _clip(1.0 - math.log10(max(call.affinity_percentile, 0.01) + 1.0) / 2.0)
    if call.affinity_nm is not None:
        return _clip(1.0 - (math.log10(max(call.affinity_nm, 1.0)) - 1.0) / 3.0)
    return 0.0


def agretopicity_feature(raw: float | None) -> float:
    """Normalized WT/mutant affinity ratio, higher meaning mutation-created binding.

    `None` means there is no positionally matched wild-type peptide to compare
    against, which is true of every indel, frameshift and neo-ORF. That is
    scored NEUTRAL rather than 1.0.

    Returning 1.0 was a substantive and unsupported claim. Combined with
    `wt_dissimilarity`, which had the same default, it granted +0.65 logit —
    just under a tenth of the model's total absolute weight — to every peptide
    lacking a comparator, on no evidence. The feature it rewards does not
    survive contact with well-powered data: Litchfield's 1,008-patient
    meta-analysis put the differential agretopic index at OR 1.03 (p = 0.79),
    and IMPROVE measured a rank-based formulation at p = 0.96 over 467
    positives.

    The failure mode that default creates is specific and bad. Frameshift and
    neo-ORF peptides are genuinely enriched for real epitopes, so once
    enumeration supports them every one would arrive pre-boosted and the
    feature would appear to be working while reporting an imputation.
    """
    if raw is None:
        return NEUTRAL
    if raw <= 0:
        return 0.0
    return _clip(math.log10(raw) / AGRETOPICITY_LOG_CEILING + 0.5)


def wt_dissimilarity_feature(raw: float | None) -> float:
    """BLOSUM62 distance from the wild-type peptide, NEUTRAL when there is none.

    Same reasoning as `agretopicity_feature`, and the same reason `_gates` does
    not treat an absent comparator as a disqualification: unavailable is not
    informative in either direction.
    """
    return NEUTRAL if raw is None else _clip(raw)


def self_dissimilarity_feature(similarity: float | None) -> float:
    """Distance from the closest human self peptide, higher being more foreign.

    `None` means the nearest-self search did not run, which is a different
    thing from "nothing similar was found". Returning NEUTRAL keeps an
    unavailable feature from silently acting as a strong positive signal, at
    the cost of the feature contributing nothing.
    """
    if similarity is None:
        return NEUTRAL
    return _clip(1.0 - similarity)


def _is_self_peptide(peptide: str, context: FeatureContext) -> bool:
    if context.self_matches is not None:
        return peptide in context.self_matches
    if context.proteome is not None:
        return context.proteome.contains_peptide(peptide)
    return False


def _gates(
    candidate: PeptideCandidate,
    mutant_call: PresentationCall,
    normal_tpm: float | None,
    context: FeatureContext,
) -> tuple[list[GateFailure], bool]:
    config = context.config
    gates: list[GateFailure] = []

    percentile = mutant_call.affinity_percentile
    if percentile is not None and percentile > config.presentation.max_affinity_percentile:
        gates.append(GateFailure.WEAK_BINDER)
    if (
        mutant_call.presentation_score is not None
        and mutant_call.presentation_score < config.presentation.min_presentation_score
    ):
        gates.append(GateFailure.LOW_PRESENTATION)

    # The proteome scan is the expensive gate, so only peptides that are
    # otherwise viable pay for it.
    self_match = False
    if config.peptides.drop_self_matching and not gates:
        self_match = _is_self_peptide(candidate.mutant_peptide, context)
        if self_match:
            gates.append(GateFailure.SELF_PEPTIDE)

    ceiling = config.expression.normal_tpm_ceiling
    if ceiling is not None and normal_tpm is not None and normal_tpm > ceiling:
        gates.append(GateFailure.NORMAL_TISSUE_EXPRESSED)

    return gates, self_match


def assemble(
    candidate: PeptideCandidate,
    mutant_call: PresentationCall,
    wildtype_call: PresentationCall | None,
    context: FeatureContext,
) -> tuple[dict[str, float], tuple[GateFailure, ...]]:
    variant = candidate.variant

    tpm = (
        context.expression.lookup(variant.transcript, variant.gene)
        if context.expression is not None
        else None
    )
    normal_tpm = (
        context.normal_expression.lookup(variant.transcript, variant.gene)
        if context.normal_expression is not None
        else None
    )

    gates, self_match = _gates(candidate, mutant_call, normal_tpm, context)

    raw_agretopicity = agretopicity(
        mutant_call.affinity_nm,
        wildtype_call.affinity_nm if wildtype_call else None,
    )

    raw_self_similarity = (
        context.self_similarity.get(candidate.mutant_peptide)
        if context.self_similarity is not None
        else None
    )

    raw_wt_dissimilarity = dissimilarity_to_wildtype(
        candidate.mutant_peptide, candidate.wildtype_peptide
    )
    #: Two features go NEUTRAL together whenever a peptide has no positionally
    #: matched wild-type, so the fact of that absence is reported as its own
    #: value. It is not scored — there is no prior worth asserting about it —
    #: but it is the only way a refit can learn what missingness is worth, and
    #: the only way an audit can tell a measured 0.5 from an unavailable one.
    wildtype_missing = candidate.wildtype_peptide is None or raw_agretopicity is None

    features: dict[str, float] = {
        "clonality": variant.ccf if variant.ccf is not None else NEUTRAL,
        "expression": expression_feature(tpm),
        "presentation": presentation_feature(mutant_call),
        "agretopicity": agretopicity_feature(raw_agretopicity),
        "tumor_selectivity": tumor_selectivity_feature(normal_tpm),
        "mutation_exposure": mutation_exposure(candidate.mutation_position, candidate.length),
        "hydrophobicity": tcr_contact_hydrophobicity(candidate.mutant_peptide),
        "wt_dissimilarity": wt_dissimilarity_feature(raw_wt_dissimilarity),
        "self_dissimilarity": self_dissimilarity_feature(raw_self_similarity),
        "tpm_raw": tpm if tpm is not None else math.nan,
        "normal_tpm_raw": normal_tpm if normal_tpm is not None else math.nan,
        "agretopicity_raw": raw_agretopicity if raw_agretopicity is not None else math.nan,
        "wt_dissimilarity_raw": (
            raw_wt_dissimilarity if raw_wt_dissimilarity is not None else math.nan
        ),
        "self_similarity_raw": (
            raw_self_similarity if raw_self_similarity is not None else math.nan
        ),
        "self_match": 1.0 if self_match else 0.0,
        "wildtype_missing": 1.0 if wildtype_missing else 0.0,
    }

    return features, tuple(gates)


def _clip(value: float, low: float = 0.0, high: float = 1.0) -> float:
    return max(low, min(high, value))
