#!/usr/bin/env python
"""Recompute benchmark metrics and redraw plots from a saved run.

    uv run --extra plots python scripts/replot_benchmark.py

Reads `results/ott2017/rankings.json` written by `run_ott2017_benchmark.py`
and re-derives every retrieval number from it. Two uses: redrawing figures
without a 24-minute MHCflurry rerun, and letting someone else verify that the
published numbers follow from the published scores and the published label
file. It has no access to the pipeline, so it cannot change a ranking -- only
re-measure one.

Score distribution is not redrawn here; it needs per-candidate features rather
than the ranking alone.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path
from typing import Final

from neoantigene.assays.metrics import compare_retrieval
from neoantigene.assays.schema import read_results

REPO_ROOT: Final[Path] = Path(__file__).resolve().parents[1]


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--results", type=Path, default=REPO_ROOT / "results" / "ott2017")
    parser.add_argument(
        "--labels",
        type=Path,
        default=REPO_ROOT / "data" / "benchmarks" / "ott2017" / "validated.tsv",
    )
    args = parser.parse_args()

    cache = args.results / "rankings.json"
    if not cache.exists():
        raise SystemExit(
            f"no saved rankings at {cache}\nRun scripts/run_ott2017_benchmark.py first."
        )

    rankings: dict[str, dict[str, float]] = json.loads(cache.read_text())
    results = read_results(args.labels)

    # Refuse to draw from a cache that shares no keys with the label file.
    # Matplotlib will happily render an axis with nothing on it, and an empty
    # figure that looks like a real one is worse than no figure at all.
    assayed = {r.key for r in results}
    overlap = {
        method: len(assayed.intersection(ranking)) for method, ranking in rankings.items()
    }
    if not any(overlap.values()):
        raise SystemExit(
            f"{cache} has no keys in common with {args.labels.name}.\n"
            f"  methods found: {', '.join(rankings) or 'none'}\n"
            "Keys must be `sample_id|peptide|HLA-A*02:01`. Re-run "
            "scripts/run_ott2017_benchmark.py to regenerate the cache."
        )
    print("  assayed pairs per method: " + ", ".join(f"{m}={n}" for m, n in overlap.items()))

    print(f"  {'method':<20} {'recall@10':>10} {'recall@20':>10} {'prec@10':>9} {'best rank':>10}")
    for summary in compare_retrieval(results, rankings):
        best = "-" if summary.best_positive_rank is None else str(summary.best_positive_rank)
        print(
            f"  {summary.method:<20} {summary.recall_at_10:>9.1%} {summary.recall_at_20:>9.1%} "
            f"{summary.precision_at_10:>8.1%} {best:>10}"
        )

    from neoantigene.plots import PlottingUnavailable, recall_curve, validated_ranks

    try:
        print(f"  {validated_ranks(results, rankings, args.results / 'validated_ranks.png')}")
        print(f"  {recall_curve(results, rankings, args.results / 'recall_curve.png')}")
    except PlottingUnavailable as exc:
        print(f"  plots skipped: {exc}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
