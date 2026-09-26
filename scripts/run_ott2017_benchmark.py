#!/usr/bin/env python
"""Run the Ott et al. 2017 benchmark and report the numbers, good or bad.

    uv run --extra presentation --extra plots python scripts/run_ott2017_benchmark.py

Structure enforces the blinding. `rank_everything()` takes only the sample
manifests and the shipped config; it has no parameter through which a label
could reach it. `score_rankings()` is the first function that opens
`validated.tsv`, and by then every ranking is frozen. The run prints a line
saying so, with the label file's digest, so the ordering is visible in the log
rather than merely asserted in prose.

Nothing here is tuned. The config is `config/default.yaml` as shipped, and the
weights in it were frozen from the literature before this script first ran.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import sys
from dataclasses import dataclass
from math import comb
from pathlib import Path
from typing import Any, Final

from neoantigene.assays.metrics import (
    RetrievalSummary,
    compare_rankings,
    compare_retrieval,
    labelled,
    retrieval_at_k,
)
from neoantigene.assays.schema import AssayResult, read_results
from neoantigene.config import PipelineConfig
from neoantigene.io.manifest import Sample
from neoantigene.models import ScoredCandidate
from neoantigene.pipeline import run as run_pipeline
from neoantigene.presentation.registry import get_backend
from neoantigene.run import configure_logging, new_run_id
from neoantigene.scoring.baselines import BASELINES

REPO_ROOT: Final[Path] = Path(__file__).resolve().parents[1]
CASE_ROOT: Final[Path] = REPO_ROOT / "data" / "benchmarks" / "ott2017"
DEFAULT_CONFIG: Final[Path] = REPO_ROOT / "config" / "default.yaml"

#: `binding_only` is the MHCflurry-only comparator: sort by predicted binding
#: and take the top N, which is the standard workflow this tool has to beat.
HEADLINE_BASELINE: Final[str] = "binding_only"


@dataclass
class PatientRun:
    sample_id: str
    scored: list[ScoredCandidate]
    variants_read: int
    variants_kept: int
    peptides: int
    pairs: int
    proteome_is_stub: bool


def rank_everything(
    case_root: Path, config_path: Path, backend_name: str
) -> tuple[list[PatientRun], PipelineConfig]:
    """Rank every patient. Sees manifests and config only -- never labels."""
    config = PipelineConfig.load(config_path)
    config.presentation.backend = backend_name
    # Keep gated candidates in the scored set: an assayed peptide that fails a
    # gate must still receive a rank, or the comparison quietly excludes the
    # cases where the gates were wrong.
    config.output.include_failed = True

    runs: list[PatientRun] = []
    directories = sorted((case_root / "inputs").glob("patient-*"))
    if not directories:
        raise SystemExit(
            f"no patient inputs under {case_root / 'inputs'}\n"
            "Run `uv run python scripts/build_ott2017_benchmark.py` first."
        )

    for directory in directories:
        sample = Sample.load(directory / "sample.yaml")
        print(f"  ranking {sample.sample_id} ({len(sample.hla)} alleles) ...", flush=True)
        result = run_pipeline(
            sample,
            config,
            backend=get_backend(config.presentation.backend),
            run_id=new_run_id(),
        )
        report = result.report
        runs.append(
            PatientRun(
                sample_id=sample.sample_id,
                scored=list(result.all_scored),
                variants_read=report.variants_read,
                variants_kept=report.variants_kept,
                peptides=report.peptides_generated,
                pairs=report.pairs_scored,
                proteome_is_stub=report.proteome_is_stub,
            )
        )
        print(
            f"    {report.variants_read} variants -> {report.peptides_generated} peptides "
            f"-> {report.pairs_scored} peptide-HLA pairs",
            flush=True,
        )
    return runs, config


def build_rankings(runs: list[PatientRun]) -> dict[str, dict[str, float]]:
    """Pooled rankings, keyed as `assays.metrics` expects.

    Scores are on one scale across patients, so pooling is legitimate and
    buys the statistical power a single 27-peptide pool does not have.
    """
    rankings: dict[str, dict[str, float]] = {"neoantigene": {}}
    for name in BASELINES:
        rankings[name] = {}

    for run in runs:
        prefix = f"{run.sample_id}|"
        for item in run.scored:
            rankings["neoantigene"][f"{prefix}{item.key}"] = item.score
        for name, ranker in BASELINES.items():
            for key, value in ranker(run.scored).items():
                rankings[name][f"{prefix}{key}"] = value
    return rankings


def score_rankings(
    label_path: Path,
    rankings: dict[str, dict[str, float]],
    k: int,
) -> tuple[list[AssayResult], list[RetrievalSummary], list[Any]]:
    """The first and only place labels are read."""
    digest = hashlib.sha256(label_path.read_bytes()).hexdigest()
    print(f"\n  opening labels now that every ranking is frozen: {label_path.name}")
    print(f"  label file sha256 = {digest}")

    results = read_results(label_path)
    return (
        results,
        compare_retrieval(results, rankings),
        compare_rankings(results, rankings, k, baseline=HEADLINE_BASELINE),
    )


def chance_level(pool: int, positives: int, k: int, hits: int) -> dict[str, float]:
    """Exact hypergeometric tail: how often would blind luck do this well?

    With a pool this small, a one-peptide swing moves recall@10 by eight
    points, so raw metric differences between methods are not interpretable on
    their own. Drawing k of `pool` without replacement is exactly the null of
    "ranking carries no information", and `comb` gives the tail in closed form
    with no sampling error and no scipy dependency.
    """
    if pool <= 0 or positives <= 0 or k <= 0:
        return {"expected_hits": 0.0, "p_value": 1.0}
    draw = min(k, pool)
    reachable = min(positives, draw)
    tail = sum(
        comb(positives, h) * comb(pool - positives, draw - h) / comb(pool, draw)
        for h in range(min(hits, reachable), reachable + 1)
    )
    return {
        "expected_hits": positives * draw / pool,
        "p_value": min(1.0, tail),
    }


def coverage(results: list[AssayResult], ranking: dict[str, float]) -> dict[str, int]:
    """How much of the assayed set the ranker actually produced.

    An assayed peptide the pipeline never generated is not a miss to be hidden;
    it is a scope limitation to be counted.
    """
    assayed = labelled(results)
    matched = [r for r in assayed if r.key in ranking]
    return {
        "labels_in_file": len(results),
        "labels_usable": len(assayed),
        "labels_scored": len(matched),
        "labels_unscored": len(assayed) - len(matched),
        "positives_usable": sum(1 for r in assayed if r.label == 1),
        "positives_scored": sum(1 for r in matched if r.label == 1),
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--case", type=Path, default=CASE_ROOT)
    parser.add_argument("--config", type=Path, default=DEFAULT_CONFIG)
    parser.add_argument("--backend", default="mhcflurry")
    parser.add_argument("--k", type=int, default=20)
    parser.add_argument("--out", type=Path, default=REPO_ROOT / "results" / "ott2017")
    parser.add_argument("--no-plots", action="store_true")
    args = parser.parse_args()

    configure_logging(new_run_id())

    print("Ott et al. 2017 benchmark (doi:10.1038/nature22991)")
    print(f"config  : {args.config}  [shipped defaults, not tuned]")
    print(f"backend : {args.backend}\n")

    print("Stage 1 -- ranking, with no access to labels")
    runs, config = rank_everything(args.case, args.config, args.backend)
    rankings = build_rankings(runs)

    if any(run.proteome_is_stub for run in runs):
        print("\n  WARNING: a stub proteome was used; self-dissimilarity is meaningless here")

    print("\nStage 2 -- scoring the rankings")
    label_path = args.case / "validated.tsv"
    results, retrieval, validation = score_rankings(label_path, rankings, args.k)
    stats = coverage(results, rankings["neoantigene"])

    print("\n  assay pool")
    for key, value in stats.items():
        print(f"    {key:<20} {value}")

    print(f"\n  retrieval over the {stats['labels_scored']} assayed, scored peptide-HLA pairs")
    print(f"  {'method':<20} {'recall@10':>10} {'recall@20':>10} {'prec@10':>9} {'best rank':>10}")
    for summary in retrieval:
        best = "-" if summary.best_positive_rank is None else str(summary.best_positive_rank)
        print(
            f"  {summary.method:<20} {summary.recall_at_10:>9.1%} {summary.recall_at_20:>9.1%} "
            f"{summary.precision_at_10:>8.1%} {best:>10}"
        )

    print(f"\n  validation rate at k={args.k} (precision among the top k)")
    for summary in validation:
        print(f"    {summary.describe()}")

    pool, positives = stats["labels_scored"], stats["positives_scored"]
    print("\n  against blind chance (exact hypergeometric, top 10)")
    print(f"    a random ordering would place {positives * 10 / max(pool, 1):.2f} positives in 10")
    chance: dict[str, dict[str, float]] = {}
    for summary in retrieval:
        test = chance_level(pool, positives, 10, summary.hits_at_10)
        chance[summary.method] = test
        verdict = "not distinguishable from chance" if test["p_value"] > 0.05 else "above chance"
        print(
            f"    {summary.method:<20} {summary.hits_at_10} hits  "
            f"p={test['p_value']:.3f}  {verdict}"
        )

    ours = next((s for s in retrieval if s.method == "neoantigene"), None)
    base = next((s for s in retrieval if s.method == HEADLINE_BASELINE), None)
    verdict = _verdict(ours, base)
    print(f"\n  verdict: {verdict}")

    per_patient = _per_patient(runs, results)
    print("\n  per patient (recall@10 / positives scored)")
    for row in per_patient:
        print(
            f"    {row['sample_id']:<16} neoantigene {row['neoantigene_recall_at_10']:>6.1%}   "
            f"{HEADLINE_BASELINE} {row['baseline_recall_at_10']:>6.1%}   "
            f"positives {row['positives_scored']}/{row['positives']}"
        )

    args.out.mkdir(parents=True, exist_ok=True)
    payload = {
        "case": "ott2017",
        "doi": "10.1038/nature22991",
        "backend": args.backend,
        "config_digest": hashlib.sha256(args.config.read_bytes()).hexdigest(),
        "label_digest": hashlib.sha256(label_path.read_bytes()).hexdigest(),
        "weights": config.weights.model_dump(),
        "coverage": stats,
        "retrieval": [s.model_dump() for s in retrieval],
        "validation": [s.model_dump() for s in validation],
        "per_patient": per_patient,
        "verdict": verdict,
        "chance": chance,
        "patients": [
            {
                "sample_id": run.sample_id,
                "variants_read": run.variants_read,
                "variants_kept": run.variants_kept,
                "peptides": run.peptides,
                "pairs": run.pairs,
            }
            for run in runs
        ],
    }
    results_json = args.out / "benchmark.json"
    results_json.write_text(json.dumps(payload, indent=2))
    print(f"\n  wrote {results_json}")

    if not args.no_plots:
        _plots(runs, results, rankings, args.out)
    return 0


def _verdict(ours: RetrievalSummary | None, base: RetrievalSummary | None) -> str:
    """State the result plainly, including when we lose."""
    if ours is None or base is None:
        return "could not compare: a method produced no ranking over the assayed pool"
    if ours.recall_at_20 > base.recall_at_20:
        return (
            f"Neo Antigene beat {HEADLINE_BASELINE} on recall@20 "
            f"({ours.recall_at_20:.1%} vs {base.recall_at_20:.1%})"
        )
    if ours.recall_at_20 < base.recall_at_20:
        return (
            f"Neo Antigene LOST to {HEADLINE_BASELINE} on recall@20 "
            f"({ours.recall_at_20:.1%} vs {base.recall_at_20:.1%})"
        )
    return (
        f"Neo Antigene tied {HEADLINE_BASELINE} on recall@20 at {ours.recall_at_20:.1%}; "
        f"best rank {ours.best_positive_rank} vs {base.best_positive_rank}"
    )


def _per_patient(runs: list[PatientRun], results: list[AssayResult]) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for run in runs:
        prefix = f"{run.sample_id}|"
        subset = [r for r in results if r.sample_id == run.sample_id]
        ours = {f"{prefix}{item.key}": item.score for item in run.scored}
        base = {
            f"{prefix}{key}": value
            for key, value in BASELINES[HEADLINE_BASELINE](run.scored).items()
        }
        mine = retrieval_at_k(subset, ours, method="neoantigene")
        theirs = retrieval_at_k(subset, base, method=HEADLINE_BASELINE)
        rows.append(
            {
                "sample_id": run.sample_id,
                "assayed": mine.assayed,
                "positives": sum(1 for r in labelled(subset) if r.label == 1),
                "positives_scored": mine.positives,
                "neoantigene_recall_at_10": mine.recall_at_10,
                "baseline_recall_at_10": theirs.recall_at_10,
                "neoantigene_best_rank": mine.best_positive_rank,
                "baseline_best_rank": theirs.best_positive_rank,
            }
        )
    return rows


def _plots(
    runs: list[PatientRun],
    results: list[AssayResult],
    rankings: dict[str, dict[str, float]],
    out: Path,
) -> None:
    try:
        from neoantigene.plots import PlottingUnavailable, recall_curve, score_distribution
        from neoantigene.plots import validated_ranks as plot_validated_ranks
    except ImportError:
        print("  plots skipped: neoantigene.plots unavailable")
        return

    every = [item for run in runs for item in run.scored]
    ordered = {
        name: rankings[name] for name in ("neoantigene", HEADLINE_BASELINE) if name in rankings
    }
    ordered.update({k: v for k, v in rankings.items() if k not in ordered})

    try:
        print(f"  {score_distribution(every, out / 'score_distribution.png')}")
        print(f"  {plot_validated_ranks(results, ordered, out / 'validated_ranks.png')}")
        print(f"  {recall_curve(results, ordered, out / 'recall_curve.png')}")
    except PlottingUnavailable as exc:
        print(f"  plots skipped: {exc}")


if __name__ == "__main__":
    sys.exit(main())
