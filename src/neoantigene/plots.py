"""Benchmark figures.

Two plots, because two are what the numbers actually need: the shape of the
score distribution, and where the experimentally validated peptides landed in
each ranking. Everything else is decoration.

matplotlib is an optional dependency (`uv sync --extra plots`), imported inside
the functions so the core package stays importable without it.
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from pathlib import Path
from typing import TYPE_CHECKING, Any, Final

from .assays.metrics import labelled
from .assays.schema import AssayResult
from .models import ScoredCandidate

if TYPE_CHECKING:  # pragma: no cover
    from matplotlib.figure import Figure

#: Muted, colour-blind-safe, and deliberately not a rainbow.
INK: Final[str] = "#27272a"
MUTED: Final[str] = "#a1a1aa"
ACCENT: Final[str] = "#2563eb"
POSITIVE: Final[str] = "#15803d"
NEGATIVE: Final[str] = "#b91c1c"


class PlottingUnavailable(RuntimeError):
    """Raised when matplotlib is not installed."""


def _pyplot() -> Any:
    try:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError as exc:  # pragma: no cover - depends on the environment
        raise PlottingUnavailable(
            "plots need matplotlib; install it with `uv sync --extra plots`"
        ) from exc
    return plt


def _style(axes: Any) -> None:
    """Thin axes, no gridlines, no box. The data carries the figure."""
    for side in ("top", "right"):
        axes.spines[side].set_visible(False)
    for side in ("left", "bottom"):
        axes.spines[side].set_linewidth(0.6)
        axes.spines[side].set_color(MUTED)
    axes.tick_params(colors=INK, labelsize=8, width=0.6, length=3)
    axes.set_axisbelow(True)


def score_distribution(
    scored: Sequence[ScoredCandidate],
    path: Path,
    title: str = "Score distribution",
) -> Path:
    """Histogram of every scored peptide-HLA pair, split by whether it passed gating.

    The useful thing to see here is that the distribution is not bimodal at the
    ends: a ranker whose scores pile up at 0 and 1 is not ranking, it is
    thresholding.
    """
    plt = _pyplot()
    passed = [s.score for s in scored if s.passed]
    gated = [s.score for s in scored if not s.passed]

    figure, axes = plt.subplots(figsize=(6.4, 3.2), dpi=160)
    bins = 40
    if gated:
        axes.hist(gated, bins=bins, range=(0, 1), color=MUTED, label=f"gated ({len(gated)})")
    if passed:
        axes.hist(passed, bins=bins, range=(0, 1), color=ACCENT, label=f"passed ({len(passed)})")

    axes.set_xlabel("score")
    axes.set_ylabel("peptide-HLA pairs")
    axes.set_title(title, fontsize=10, color=INK, loc="left")
    if passed or gated:
        axes.legend(frameon=False, fontsize=8)
    _style(axes)
    return _save(figure, path)


def validated_ranks(
    results: Sequence[AssayResult],
    rankings: Mapping[str, Mapping[str, float]],
    path: Path,
    title: str = "Where the validated peptides landed",
    max_rank: int | None = None,
) -> Path:
    """Rank position of every assay-positive peptide, one row per method.

    This is the plot that can embarrass the ranker, which is the point: if the
    positives are scattered uniformly the method is not working, and that is
    visible here in a way a single recall number hides.
    """
    plt = _pyplot()
    assayed = labelled(results)
    positives = {r.key for r in assayed if r.label == 1}
    pool = {r.key for r in assayed}

    methods = list(rankings)
    figure, axes = plt.subplots(figsize=(6.4, 0.7 * len(methods) + 1.6), dpi=160)

    # Ranks run over the assayed peptides a method actually scored, which is
    # smaller than the label file whenever the pipeline did not generate one.
    # Labelling the axis with the full label count would overstate the pool.
    scored = {key for ranking in rankings.values() for key in ranking if key in pool}
    ceiling = max_rank or max(len(scored), 1)
    for row, method in enumerate(methods):
        ranks = _ranks_of_positives(rankings[method], pool, positives)
        axes.scatter(
            ranks,
            [row] * len(ranks),
            s=42,
            marker="|",
            linewidths=1.4,
            color=POSITIVE if ranks else MUTED,
        )
        best = min(ranks) if ranks else None
        label = f"best rank {best}" if best is not None else "no positives ranked"
        axes.annotate(
            label,
            xy=(ceiling, row),
            xytext=(4, 0),
            textcoords="offset points",
            va="center",
            fontsize=7,
            color=MUTED,
        )

    axes.set_yticks(range(len(methods)))
    axes.set_yticklabels(methods, fontsize=8)
    axes.set_xlim(0.5, ceiling + 0.5)
    axes.set_ylim(-0.6, len(methods) - 0.4)
    axes.invert_yaxis()
    axes.set_xlabel(f"rank among the {len(scored)} assayed peptides scored (lower is better)")
    axes.set_title(title, fontsize=10, color=INK, loc="left")
    _style(axes)
    return _save(figure, path)


def recall_curve(
    results: Sequence[AssayResult],
    rankings: Mapping[str, Mapping[str, float]],
    path: Path,
    title: str = "Validated peptides recovered by depth",
) -> Path:
    """Cumulative count of validated positives against list depth.

    A win that only exists at one k is a fluke. This shows the whole curve so
    the reader can see whether the advantage holds.
    """
    plt = _pyplot()
    assayed = labelled(results)
    positives = {r.key for r in assayed if r.label == 1}
    pool = {r.key for r in assayed}

    figure, axes = plt.subplots(figsize=(6.4, 3.4), dpi=160)
    # Depth and the chance line must both run over the peptides that were
    # actually ranked. Using the full label count would stretch the axis past
    # the end of every list and flatten the chance line below true chance.
    scored = {key for ranking in rankings.values() for key in ranking if key in pool}
    ranked_positives = len(positives & scored)
    depths = range(1, max(len(scored), 1) + 1)
    palette = [ACCENT, INK, MUTED, POSITIVE, NEGATIVE]

    for index, (method, ranking) in enumerate(rankings.items()):
        ranks = sorted(_ranks_of_positives(ranking, pool, positives))
        recovered = [sum(1 for r in ranks if r <= depth) for depth in depths]
        axes.step(
            list(depths),
            recovered,
            where="post",
            linewidth=1.4 if index == 0 else 1.0,
            color=palette[index % len(palette)],
            label=method,
        )

    axes.plot(
        list(depths),
        [ranked_positives * depth / max(len(scored), 1) for depth in depths],
        linestyle=(0, (3, 3)),
        linewidth=0.9,
        color=MUTED,
        label="chance",
    )
    axes.set_xlabel("depth into the ranked list")
    axes.set_ylabel("validated peptides recovered")
    axes.set_title(title, fontsize=10, color=INK, loc="left")
    axes.legend(frameon=False, fontsize=8, loc="lower right")
    _style(axes)
    return _save(figure, path)


def _ranks_of_positives(
    ranking: Mapping[str, float],
    pool: set[str],
    positives: set[str],
) -> list[int]:
    """1-based ranks of the validated positives, over the assayed pool only.

    Restricting to the assayed pool is the same rule the metrics use: a method
    cannot be rewarded or punished for peptides nobody tested.
    """
    scored = [(key, ranking[key]) for key in pool if key in ranking]
    ordered = sorted(scored, key=lambda item: (-item[1], item[0]))
    return [index for index, (key, _) in enumerate(ordered, start=1) if key in positives]


def _save(figure: Figure, path: Path) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    figure.tight_layout()
    figure.savefig(path, bbox_inches="tight", transparent=False)
    figure.clf()
    return path
