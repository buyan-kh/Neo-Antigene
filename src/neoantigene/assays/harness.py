"""A blinded evaluation harness for neoantigen rankings.

The published literature cannot currently tell whether its own methods work.
A model trained on the HLA allele alone, with no peptide, outscored all
fifteen entrants in one published benchmark. PRIME 1.0 had roughly 70% of its
training peptides inside a CEDAR-derived evaluation set, as measured by a
third party, and lost performance when that overlap was removed. IEDB's
curation manual permits a prediction to supersede a coarser experimental
restriction.

Contamination is not uniform, which is the reason to measure rather than
assume: IMPROVE's in-house set has exactly zero overlap with the CEDAR
benchmark on the same join key. An audit that presumed leakage everywhere
would be making the error it exists to detect.

This module is the measurement layer that absence implies. It takes a label
set and any number of named rankings — from this package or from any other
tool — and produces one auditable report over the peptides they all scored.

Three properties are the point:

**It works on other people's rankings.** A harness only usable on its author's
own method measures nothing. Input is a peptide-to-score mapping; where the
scores came from is irrelevant.

**It states what the design can detect before it reports what happened.**
`DesignCapacity` is computed from the pool alone and printed first. When a
pool cannot reject chance at any outcome, that appears above the hit counts
rather than in a caveat after them, because a reader who has already seen
"3 of 10" has formed an impression that no later sentence undoes.

**It cannot be read as a win it did not earn.** Every method is restricted to
the shared pool, ties resolve against the ranking, intervals come from a
bootstrap whose coverage was measured rather than assumed, and the verdict is
stated in `summary()` rather than left to the reader to assemble.
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence

from pydantic import BaseModel, Field

from .metrics import RetrievalSummary, compare_retrieval, labelled
from .schema import AssayCall, AssayResult
from .stats import (
    DEFAULT_ROUNDS,
    DEFAULT_SEED,
    AucSignificance,
    DesignCapacity,
    PairedComparison,
    StratifiedAuc,
    TopKSignificance,
    auc_significance,
    design_capacity,
    paired_comparison,
    stratified_auc,
    top_k_significance,
)


class SourceFile(BaseModel):
    """One input and its digest, so a report can be tied to exact bytes."""

    role: str
    path: str
    digest: str


class HarnessReport(BaseModel):
    run_id: str
    package_version: str
    k: int
    alpha: float
    rounds: int

    sources: list[SourceFile] = Field(default_factory=list)
    label_counts: dict[str, int] = Field(default_factory=dict)
    patients: list[str] = Field(default_factory=list)

    capacity: DesignCapacity
    retrieval: list[RetrievalSummary] = Field(default_factory=list)
    significance: list[TopKSignificance] = Field(default_factory=list)
    pooled_auc: list[AucSignificance] = Field(default_factory=list)
    stratified: list[StratifiedAuc] = Field(default_factory=list)
    comparisons: list[PairedComparison] = Field(default_factory=list)
    baseline: str | None = None

    @property
    def separated(self) -> list[str]:
        """Methods that beat the baseline on the paired ordering test."""
        return [
            c.method_a
            for c in self.comparisons
            if c.auc_p_value is not None and c.auc_p_value <= self.alpha
        ]

    @property
    def above_chance(self) -> list[str]:
        return [s.method for s in self.significance if s.p_value <= self.alpha]

    @property
    def supports_a_claim(self) -> bool:
        """Whether any analysis in this report could sustain a claim at all."""
        return bool(self.above_chance or self.separated)

    def summary(self) -> str:
        """The one-sentence reading, stated by the harness rather than inferred.

        Written to be quotable as-is, because the common failure is a reader
        lifting a number out of a table and dropping the condition attached to
        it.

        The top-k and ordering analyses are reported separately on purpose.
        They have very different power — a top-k cut uses k of n observations
        while concordance uses every positive-negative pair — so a pool can be
        structurally unable to support a precision@k claim while the ordering
        test still resolves a difference. Collapsing those into one verdict
        would be wrong in one direction or the other.
        """
        underpowered = not self.capacity.can_reach_significance
        if underpowered and not self.separated:
            return (
                f"INCONCLUSIVE: at k={self.k} no outcome this pool can produce reaches "
                f"p <= {self.alpha:g}, and the ordering test separated nothing either. "
                f"This is not evidence about any method evaluated here."
            )
        if underpowered:
            return (
                f"ORDERING ONLY: precision@{self.k} cannot be significant on this pool at "
                f"any outcome, so quote no hit rate from it. {', '.join(self.separated)} "
                f"did separate from {self.baseline} on the paired ordering test, which "
                f"uses every pair rather than the top {self.k}. Confirm on more patients."
            )
        if not self.above_chance:
            return (
                f"NO METHOD CLEARS CHANCE: none reached p <= {self.alpha:g} at k={self.k}. "
                f"The pool could have shown it — {self.capacity.hits_needed} hits were "
                f"needed — so this is a result, not a missing measurement."
            )
        if self.baseline and not self.separated:
            return (
                f"ABOVE CHANCE, NOT ABOVE BASELINE: {', '.join(self.above_chance)} cleared "
                f"chance, but no method separated from {self.baseline} on a paired test. "
                f"Beating chance is not the product claim."
            )
        return (
            f"SEPARATION FOUND: {', '.join(self.separated)} separated from {self.baseline} "
            f"at p <= {self.alpha:g}. Confirm on an independent cohort before quoting."
        )

    def to_markdown(self) -> str:
        lines = [
            "# Neoantigen ranking evaluation",
            "",
            f"**{self.summary()}**",
            "",
            f"- run id: `{self.run_id}`",
            f"- neoantigene: `{self.package_version}`",
            f"- k: {self.k}, alpha: {self.alpha:g}, permutations: {self.rounds:,}",
            "",
            "## Inputs",
            "",
            "| role | path | digest |",
            "| --- | --- | --- |",
        ]
        lines += [f"| {s.role} | `{s.path}` | `{s.digest}` |" for s in self.sources]

        counts = ", ".join(f"{name} {count}" for name, count in sorted(self.label_counts.items()))
        lines += [
            "",
            f"Labels: {counts}. Patients: {len(self.patients)}"
            f"{' (' + ', '.join(self.patients) + ')' if self.patients else ''}.",
            "",
            "## What this design could detect",
            "",
            "Computed from the assayed pool alone, before any ranking was scored.",
            "",
            "```",
            self.capacity.describe(),
            "```",
            "",
            "## Retrieval",
            "",
            "| method | recall@10 | recall@20 | precision@10 | best positive rank |",
            "| --- | ---: | ---: | ---: | ---: |",
        ]
        for item in self.retrieval:
            best = "none" if item.best_positive_rank is None else f"#{item.best_positive_rank}"
            lines.append(
                f"| {item.method} | {item.recall_at_10:.1%} ({item.hits_at_10}/{item.positives}) "
                f"| {item.recall_at_20:.1%} ({item.hits_at_20}/{item.positives}) "
                f"| {item.precision_at_10:.1%} | {best} |"
            )

        lines += [
            "",
            "## Against chance",
            "",
            f"Exact hypergeometric tail at k={self.k}.",
            "",
            "| method | hits | chance | p |",
            "| --- | ---: | ---: | ---: |",
        ]
        lines += [
            f"| {s.method} | {s.hits} | {s.expected:.2f} | {s.p_value:.3f} |"
            for s in self.significance
        ]

        lines += [
            "",
            "## Ordering quality",
            "",
            "Concordance over every positive-negative pair, which uses far more of the",
            "labels than a top-k cut. `stratified` counts only within-patient pairs, so",
            "no peptide is ever compared against one from a different tumor.",
            "",
            "| method | pooled auc | p | stratified auc | 95% CI | pairs |",
            "| --- | ---: | ---: | ---: | :---: | ---: |",
        ]
        by_method = {s.method: s for s in self.stratified}
        for pooled in self.pooled_auc:
            strat = by_method.get(pooled.method)
            auc = "n/a" if pooled.auc is None else f"{pooled.auc:.3f}"
            p = "n/a" if pooled.p_value is None else f"{pooled.p_value:.3f}"
            s_auc = "n/a" if strat is None or strat.auc is None else f"{strat.auc:.3f}"
            interval = "n/a"
            if strat is not None and strat.ci_low is not None and strat.ci_high is not None:
                interval = f"[{strat.ci_low:.3f}, {strat.ci_high:.3f}]"
            pairs = "0" if strat is None else f"{strat.pairs}"
            lines.append(f"| {pooled.method} | {auc} | {p} | {s_auc} | {interval} | {pairs} |")

        if self.comparisons:
            lines += ["", f"## Paired against `{self.baseline}`", ""]
            lines += [
                "Permutation test on the gap between two rankings over the peptides they",
                "both scored. This is the comparison a product claim rests on: clearing",
                "chance is not the same as beating the tool a lab already uses.",
                "",
            ]
            for comparison in self.comparisons:
                lines += ["```", comparison.describe(), "```", ""]

        lines += ["## Reading this honestly", ""] + [f"- {note}" for note in self.caveats()]
        return "\n".join(lines) + "\n"

    def caveats(self) -> list[str]:
        notes: list[str] = []
        if not self.capacity.can_reach_significance:
            notes.append(
                f"**No precision@{self.k} figure from this report can be significant**, "
                f"whatever it reads. Do not quote a hit rate from it. The ordering "
                f"columns are the only ones with power on a pool this size."
            )
        if self.capacity.positives < 100:
            notes.append(
                f"{self.capacity.positives} positives. Detecting a doubling of "
                f"precision@k over a baseline takes on the order of 100, so a gap of one "
                f"or two peptides here is within sampling error."
            )
        if len(self.patients) < 3:
            notes.append(
                f"{len(self.patients)} patient(s). Per-patient response rates vary "
                f"severalfold, and with fewer than three the stratified interval has no "
                f"between-patient variance to estimate."
            )
        if self.baseline is None:
            notes.append(
                "No baseline given, so nothing here speaks to whether these rankings "
                "beat sorting by predicted presentation."
            )
        notes.append(
            "Peptides ranked but never assayed are excluded, not counted as negatives. "
            "Not testing a peptide is not evidence against it."
        )
        notes.append(
            "If any ranking used a predictor fine-tuned on the same labels, this report "
            "measures that overlap and not the method. Record predictor versions."
        )
        return notes


def build_report(
    results: Sequence[AssayResult],
    rankings: Mapping[str, Mapping[str, float]],
    run_id: str,
    package_version: str,
    k: int = 20,
    alpha: float = 0.05,
    baseline: str | None = None,
    sources: Sequence[SourceFile] = (),
    rounds: int = DEFAULT_ROUNDS,
    seed: int = DEFAULT_SEED,
) -> HarnessReport:
    """Evaluate every ranking over the peptides all of them scored.

    Takes loaded data rather than paths so it is testable without a filesystem
    and callable on rankings this package did not produce.
    """
    if not rankings:
        raise ValueError("at least one ranking is required")

    shared = set.intersection(*(set(r) for r in rankings.values()))
    pool = [r for r in labelled(results) if r.key in shared]
    positives = sum(1 for r in pool if r.label == 1)
    depth = min(k, len(pool))

    reference = baseline if baseline in rankings else None
    comparisons = [
        paired_comparison(
            pool,
            ranking,
            rankings[reference],
            k=depth,
            method_a=name,
            method_b=reference,
            rounds=rounds,
            seed=seed,
        )
        for name, ranking in rankings.items()
        if reference is not None and name != reference
    ]

    return HarnessReport(
        run_id=run_id,
        package_version=package_version,
        k=depth,
        alpha=alpha,
        rounds=rounds,
        sources=list(sources),
        label_counts=_label_counts(results),
        patients=sorted({r.sample_id for r in pool}),
        capacity=design_capacity(len(pool), positives, depth, alpha),
        retrieval=compare_retrieval(pool, rankings),
        significance=[
            top_k_significance(pool, ranking, depth, method=name)
            for name, ranking in rankings.items()
        ],
        pooled_auc=[
            auc_significance(pool, ranking, method=name, rounds=rounds, seed=seed)
            for name, ranking in rankings.items()
        ],
        stratified=[
            stratified_auc(pool, ranking, method=name, rounds=rounds, seed=seed)
            for name, ranking in rankings.items()
        ],
        comparisons=comparisons,
        baseline=reference,
    )


def _label_counts(results: Sequence[AssayResult]) -> dict[str, int]:
    counts = dict.fromkeys((call.value for call in AssayCall), 0)
    for result in results:
        counts[result.call.value] += 1
    return counts
