"""End-to-end ranking run.

Stage order mirrors the locked pipeline: ingest -> somatic filtering ->
peptide enumeration -> processing/presentation -> ranking. Structure-based
filtering and assay feedback attach after this, never inside it.

Each stage is a module-level function that takes models and returns models.
`run` is the composition, and is the only place that touches the filesystem
or the clock.
"""

from __future__ import annotations

import logging
from collections.abc import Mapping, Sequence
from pathlib import Path

from pydantic import BaseModel, ConfigDict, Field

from .config import PipelineConfig
from .io.expression import ExpressionTable, NormalExpressionReference
from .io.manifest import Sample
from .io.vcf import read_variant_tsv, read_vep_vcf
from .models import GateFailure, PeptideCandidate, PresentationCall, ScoredCandidate, Variant
from .peptides.generate import PeptideGenerationError, generate_for_variant
from .peptides.proteome import ProteomeIndex
from .peptides.selfsim import SelfProteome
from .presentation.base import PresentationBackend
from .presentation.registry import get_backend
from .run import RunManifest
from .scoring.features import FeatureContext, assemble
from .scoring.rank import score_all, shortlist
from .variants.filters import FilterReport, apply_variant_filters

logger = logging.getLogger(__name__)


class RunReport(BaseModel):
    variants_read: int = 0
    variants_kept: int = 0
    filter_report: FilterReport = Field(default_factory=FilterReport)
    peptides_generated: int = 0
    peptide_errors: dict[str, int] = Field(default_factory=dict)
    pairs_scored: int = 0
    pairs_gated: dict[GateFailure, int] = Field(default_factory=dict)
    shortlisted: int = 0

    #: Size of the loaded reference proteome, and whether it is too small to be
    #: a real one. A stub reference makes the self-peptide gate and the
    #: self-dissimilarity feature meaningless, which consumers must disclose.
    proteome_proteins: int = 0
    proteome_is_stub: bool = False

    #: False when the sample carried no RNA, in which case the `expression` and
    #: `tumor_selectivity` features fall back to a neutral 0.5 and contribute
    #: nothing but their bias.
    expression_available: bool = False
    normal_expression_available: bool = False

    def counts(self) -> dict[str, int]:
        """Flat integer summary, for the run manifest."""
        return {
            "variants_read": self.variants_read,
            "variants_kept": self.variants_kept,
            "peptides_generated": self.peptides_generated,
            "pairs_scored": self.pairs_scored,
            "pairs_gated": sum(self.pairs_gated.values()),
            "shortlisted": self.shortlisted,
        }


class RunResult(BaseModel):
    model_config = ConfigDict(arbitrary_types_allowed=True)

    sample: Sample
    config: PipelineConfig
    manifest: RunManifest
    ranked: list[ScoredCandidate]
    all_scored: list[ScoredCandidate]
    report: RunReport


def read_variants(sample: Sample) -> list[Variant]:
    processed = sample.processed
    if processed.somatic_vcf is not None:
        return list(read_vep_vcf(processed.somatic_vcf, processed.tumor_sample_name))
    assert processed.variant_tsv is not None  # guaranteed by ProcessedInputs
    return list(read_variant_tsv(processed.variant_tsv))


def resolve_config(sample: Sample, config: PipelineConfig) -> PipelineConfig:
    """Fold sample-level facts into the config so the run has one source of truth."""
    if sample.tumor_purity is None:
        return config
    return config.model_copy(
        update={
            "variant_filters": config.variant_filters.model_copy(
                update={"tumor_purity": sample.tumor_purity}
            )
        }
    )


def generate_peptides(
    variants: Sequence[Variant],
    proteome: ProteomeIndex,
    config: PipelineConfig,
) -> tuple[list[PeptideCandidate], dict[str, int]]:
    candidates: list[PeptideCandidate] = []
    errors: dict[str, int] = {}
    for variant in variants:
        try:
            candidates.extend(
                generate_for_variant(
                    variant,
                    proteome,
                    config.peptides.lengths,
                    config.peptides.flank_length,
                )
            )
        except PeptideGenerationError as exc:
            reason = type(exc).__name__
            errors[reason] = errors.get(reason, 0) + 1
            logger.debug("skipping %s: %s", variant.hgvsp_short, exc)
    return candidates, errors


def scan_nearest_self(
    candidates: Sequence[PeptideCandidate],
    proteome: ProteomeIndex,
    calls: Mapping[tuple[str, str], PresentationCall] | None = None,
    config: PipelineConfig | None = None,
    reference: SelfProteome | None = None,
) -> dict[str, float] | None:
    """Nearest-self similarity for plausibly presented peptides, or None if unusable.

    Skipped against a stub reference: a two-protein FASTA would report every
    peptide as maximally foreign, which is a fabricated signal rather than a
    weak one.

    Restricted to peptides that at least one allele is predicted to present.
    That is not only a saving. The search seeds on exact 4-mers, and there are
    only 20^4 of them, so a query set of tens of thousands of peptides covers
    most of the seed space and every proteome position becomes a hit — the
    index stops discriminating and the scan degrades to a full comparison.
    Scanning the gate-passing subset keeps the seeding selective, and peptides
    no allele presents are dropped downstream regardless.
    """
    if proteome.is_stub:
        logger.warning("skipping nearest-self search: reference proteome is a stub")
        return None

    peptides = _presentable_peptides(candidates, calls, config)
    if not peptides:
        logger.info("no peptides cleared presentation, skipping nearest-self search")
        return None

    reference = reference or SelfProteome(proteome.sequences())
    similarities = reference.nearest_similarity(peptides)
    logger.info(
        "scanned %d of %d peptides against %d self residues",
        len(peptides),
        len({c.mutant_peptide for c in candidates}),
        reference.residues,
    )
    return similarities


def scan_self_matches(
    candidates: Sequence[PeptideCandidate],
    proteome: ProteomeIndex,
    calls: Mapping[tuple[str, str], PresentationCall] | None = None,
    config: PipelineConfig | None = None,
    reference: SelfProteome | None = None,
) -> frozenset[str] | None:
    """Mutant peptides that occur verbatim in the reference proteome.

    Precomputed in one batched pass so the self-peptide gate reads a set
    membership instead of scanning the whole reference per candidate. `None`
    means the gate is disabled, in which case feature assembly does not
    consult it at all.

    Unlike the nearest-self search this still runs against a stub reference.
    An exact match found in two proteins is a real match; the stub only limits
    how many matches can be found, which the run report already discloses.
    """
    if config is not None and not config.peptides.drop_self_matching:
        return None

    peptides = _presentable_peptides(candidates, calls, config)
    if not peptides:
        return frozenset()

    reference = reference or SelfProteome(proteome.sequences())
    matches = reference.exact_matches(peptides)
    logger.info(
        "self-peptide gate: %d of %d presented peptides match the reference exactly",
        len(matches),
        len(peptides),
    )
    return frozenset(matches)


def _presentable_peptides(
    candidates: Sequence[PeptideCandidate],
    calls: Mapping[tuple[str, str], PresentationCall] | None,
    config: PipelineConfig | None,
) -> list[str]:
    """Mutant peptides at least one allele is predicted to present."""
    everything = sorted({c.mutant_peptide for c in candidates})
    if calls is None or config is None:
        return everything

    ceiling = config.presentation.max_affinity_percentile
    floor = config.presentation.min_presentation_score
    keep: set[str] = set()
    for (peptide, _allele), call in calls.items():
        percentile = call.affinity_percentile
        if percentile is not None and percentile > ceiling:
            continue
        score = call.presentation_score
        if score is not None and score < floor:
            continue
        keep.add(peptide)
    return sorted(keep & set(everything))


def predict_presentation(
    candidates: Sequence[PeptideCandidate],
    alleles: Sequence[str],
    backend: PresentationBackend,
    config: PipelineConfig,
) -> dict[tuple[str, str], PresentationCall]:
    """Score every distinct peptide once, mutant and wild-type together."""
    flanks: dict[str, tuple[str, str]] = {}
    for candidate in candidates:
        flanks.setdefault(candidate.mutant_peptide, (candidate.n_flank, candidate.c_flank))
    if config.presentation.score_wildtype:
        for candidate in candidates:
            wildtype = candidate.wildtype_peptide
            if wildtype is not None:
                flanks.setdefault(wildtype, (candidate.n_flank, candidate.c_flank))

    peptides = [p for p in flanks if backend.supports_length(len(p))]
    calls = backend.predict(
        peptides,
        list(alleles),
        n_flanks=[flanks[p][0] for p in peptides],
        c_flanks=[flanks[p][1] for p in peptides],
    )
    return {call.key: call for call in calls}


def build_scored(
    candidates: Sequence[PeptideCandidate],
    alleles: Sequence[str],
    calls: dict[tuple[str, str], PresentationCall],
    context: FeatureContext,
) -> list[ScoredCandidate]:
    scored: list[ScoredCandidate] = []
    for candidate in candidates:
        for allele in alleles:
            mutant_call = calls.get((candidate.mutant_peptide, allele))
            if mutant_call is None:
                continue
            wildtype_call = (
                calls.get((candidate.wildtype_peptide, allele))
                if candidate.wildtype_peptide is not None
                else None
            )
            features, gates = assemble(candidate, mutant_call, wildtype_call, context)
            scored.append(
                ScoredCandidate(
                    candidate=candidate,
                    allele=allele,
                    mutant_call=mutant_call,
                    wildtype_call=wildtype_call,
                    features=features,
                    gate_failures=gates,
                )
            )
    return scored


def run(
    sample: Sample,
    config: PipelineConfig,
    backend: PresentationBackend | None = None,
    run_id: str | None = None,
) -> RunResult:
    alleles = sample.require_hla()
    proteome_path = sample.require_proteome()
    resolved = resolve_config(sample, config)
    predictor = backend or get_backend(resolved.presentation.backend)

    manifest = RunManifest.start(
        sample_id=sample.sample_id,
        config=resolved,
        backend=predictor.name,
        inputs=sample.input_paths(),
        run_id=run_id,
    )
    logger.info(
        "run %s: sample %s, %d alleles, backend %s",
        manifest.run_id,
        sample.sample_id,
        len(alleles),
        predictor.name,
    )

    report = RunReport()
    variants = read_variants(sample)
    report.variants_read = len(variants)
    logger.info("read %d protein-altering variants", len(variants))

    expression = _load(ExpressionTable, sample.processed.expression_tsv)
    normal_expression = _load(NormalExpressionReference, sample.processed.normal_expression_tsv)
    report.expression_available = expression is not None
    report.normal_expression_available = normal_expression is not None
    if expression is None:
        logger.warning(
            "no expression table for %s: the expression feature is neutral for every "
            "candidate and cannot discriminate",
            sample.sample_id,
        )

    kept, filter_report = apply_variant_filters(
        variants, resolved.variant_filters, resolved.expression, expression
    )
    report.variants_kept = len(kept)
    report.filter_report = filter_report
    logger.info("somatic filtering: %s", filter_report.summary())

    proteome = ProteomeIndex.from_fasta(proteome_path)
    report.proteome_proteins = proteome.protein_count
    report.proteome_is_stub = proteome.is_stub
    logger.info("loaded proteome with %d transcripts", len(proteome))
    if proteome.is_stub:
        logger.warning(
            "%s holds only %d proteins, which is not a reference proteome: the "
            "self-peptide gate and self-dissimilarity feature are meaningless. Run "
            "`neoantigene fetch-proteome` and point the manifest at the result.",
            proteome_path,
            proteome.protein_count,
        )

    candidates, errors = generate_peptides(kept, proteome, resolved)
    report.peptides_generated = len(candidates)
    report.peptide_errors = errors
    logger.info("generated %d unique mutant peptides", len(candidates))
    if errors:
        logger.info("peptide generation skips: %s", errors)

    if not candidates:
        return RunResult(
            sample=sample,
            config=resolved,
            manifest=manifest.finish(report.counts()),
            ranked=[],
            all_scored=[],
            report=report,
        )

    calls = predict_presentation(candidates, alleles, predictor, resolved)
    logger.info("predicted %d peptide-allele pairs", len(calls))

    # One encoded copy of the reference, shared by both proteome scans: the
    # nearest-self search and the exact self-match gate.
    reference = SelfProteome(proteome.sequences())
    context = FeatureContext(
        config=resolved,
        expression=expression,
        normal_expression=normal_expression,
        proteome=proteome,
        self_similarity=scan_nearest_self(candidates, proteome, calls, resolved, reference),
        self_matches=scan_self_matches(candidates, proteome, calls, resolved, reference),
    )
    all_scored = score_all(build_scored(candidates, alleles, calls, context), resolved.weights)
    report.pairs_scored = len(all_scored)
    for item in all_scored:
        for gate in item.gate_failures:
            report.pairs_gated[gate] = report.pairs_gated.get(gate, 0) + 1

    ranked = shortlist(all_scored, resolved.output)
    report.shortlisted = len(ranked)
    logger.info("shortlisted %d candidates", len(ranked))

    return RunResult(
        sample=sample,
        config=resolved,
        manifest=manifest.finish(report.counts()),
        ranked=ranked,
        all_scored=all_scored,
        report=report,
    )


def _load[T: ExpressionTable](table: type[T], path: Path | None) -> T | None:
    return None if path is None else table.from_tsv(path)
