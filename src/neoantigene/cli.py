"""Neo Antigene command line.

Thin shell over the library: every command parses arguments, calls one or two
library functions, and prints. Nothing here contains ranking logic.
"""

from __future__ import annotations

import csv
import json
from collections.abc import Mapping, Sequence
from enum import StrEnum
from pathlib import Path
from typing import Annotated, Any

import typer

from . import __version__
from .assays.harness import SourceFile, build_report
from .assays.metrics import compare_rankings, validation_rate_at_k
from .assays.schema import AssayResult, AssayType, read_results, write_request
from .assays.stats import paired_comparison, stratified_auc, top_k_significance
from .config import PipelineConfig
from .hla import parse_hla_string
from .io.manifest import Sample
from .io.writers import write_features_json, write_tsv
from .learning.active import select_batch
from .learning.dataset import build_training_set, load_feature_records
from .learning.train import InsufficientData, compare_to_prior, fit_weights
from .models import ScoredCandidate
from .peptides.reference import (
    ENSEMBL_RELEASE,
    ReferenceDownloadError,
    fetch_proteome,
    provenance,
)
from .pipeline import run as run_pipeline
from .presentation.registry import DEVELOPMENT_BACKENDS, available, get_backend
from .run import configure_logging, file_digest, new_run_id
from .scoring.baselines import BASELINES
from .scoring.rank import contributions
from .wedge import pdac

app = typer.Typer(add_completion=False, no_args_is_help=True, help=__doc__)


class SelectMode(StrEnum):
    SCORE = "score"
    ACTIVE = "active"


ManifestArg = Annotated[Path, typer.Argument(help="Sample manifest YAML.", exists=True)]
ConfigOpt = Annotated[Path | None, typer.Option("--config", "-c", exists=True)]
OutDirOpt = Annotated[Path, typer.Option("--out-dir", "-o")]


@app.command()
def rank(
    manifest: ManifestArg,
    out_dir: OutDirOpt = Path("results"),
    config_path: ConfigOpt = None,
    backend: Annotated[str | None, typer.Option("--backend", help=f"One of {available()}")] = None,
    top_n: Annotated[int | None, typer.Option("--top-n", min=1)] = None,
    select: Annotated[SelectMode, typer.Option("--select")] = SelectMode.SCORE,
    exploitation: Annotated[float, typer.Option("--exploitation", min=0.0, max=1.0)] = 0.6,
    allow_null_backend: Annotated[bool, typer.Option("--allow-null-backend")] = False,
    json_logs: Annotated[bool, typer.Option("--json-logs")] = False,
    verbose: Annotated[bool, typer.Option("--verbose", "-v")] = False,
) -> None:
    """Rank neoantigen candidates for one sample."""
    run_id = new_run_id()
    configure_logging(run_id, verbose=verbose, json_logs=json_logs)

    sample = Sample.load(manifest)
    config = PipelineConfig.load(config_path)
    if backend is not None:
        config.presentation.backend = backend
    if top_n is not None:
        config.output.top_n = top_n

    if config.presentation.backend in DEVELOPMENT_BACKENDS and not allow_null_backend:
        raise typer.BadParameter(
            f"the {config.presentation.backend!r} backend produces meaningless scores "
            "and is for development only; pass --allow-null-backend to acknowledge this"
        )

    result = run_pipeline(
        sample,
        config,
        backend=get_backend(config.presentation.backend),
        run_id=run_id,
    )

    ranked = _select(result.all_scored, result.ranked, select, exploitation, result.config)

    out_dir.mkdir(parents=True, exist_ok=True)
    tsv_path = out_dir / f"{sample.sample_id}.ranked.tsv"
    write_tsv(ranked, sample.sample_id, tsv_path, run_id=run_id)
    write_features_json(
        result.all_scored,
        sample.sample_id,
        out_dir / f"{sample.sample_id}.features.json",
        run_id=run_id,
    )
    manifest_path = result.manifest.write(out_dir / f"{sample.sample_id}.run.json")

    report = result.report
    typer.echo(f"run id             : {run_id}")
    typer.echo(f"variants read      : {report.variants_read}")
    typer.echo(f"somatic filtering  : {report.filter_report.summary()}")
    typer.echo(f"peptides generated : {report.peptides_generated}")
    if report.peptide_errors:
        typer.echo(f"peptide skips      : {report.peptide_errors}")
    typer.echo(f"peptide-HLA pairs  : {report.pairs_scored}")
    if report.pairs_gated:
        gated = {gate.value: count for gate, count in report.pairs_gated.items()}
        typer.echo(f"gated pairs        : {gated}")
    typer.echo(f"shortlisted        : {len(ranked)} -> {tsv_path}")
    typer.echo(f"run manifest       : {manifest_path}")

    matches = pdac.match_variants([s.candidate.variant for s in ranked])
    if matches:
        typer.echo("shared antigens    : " + ", ".join(m.key for m in matches))


def _select(
    all_scored: list[ScoredCandidate],
    ranked: list[ScoredCandidate],
    mode: SelectMode,
    exploitation: float,
    config: PipelineConfig,
) -> list[ScoredCandidate]:
    if mode is SelectMode.SCORE:
        return ranked
    batch = select_batch(
        all_scored,
        batch_size=config.output.top_n,
        exploitation_weight=exploitation,
        max_per_variant=config.output.max_per_variant,
    )
    return [item.candidate for item in batch]


@app.command("request-assays")
def request_assays(
    ranked_tsv: Annotated[Path, typer.Argument(help="Ranked shortlist TSV.", exists=True)],
    out: Annotated[Path, typer.Option("--out", "-o")],
    assay: Annotated[AssayType, typer.Option("--assay")] = AssayType.IFNG_ELISPOT,
) -> None:
    """Emit a synthesis/assay request sheet for the lab to fill in."""
    records = _read_tsv(ranked_tsv)
    write_request(records, out, assay=assay)
    typer.echo(f"wrote {len(records)} peptides to {out}")


@app.command()
def evaluate(
    results_tsv: Annotated[Path, typer.Argument(help="Completed assay results.", exists=True)],
    ranked_tsv: Annotated[Path, typer.Option("--ranked", exists=True)],
    baseline_tsv: Annotated[
        Path | None,
        typer.Option("--baseline", exists=True, help="Comparator ranking over the same peptides."),
    ] = None,
    k: Annotated[int, typer.Option("--k", min=1)] = 20,
    stats: Annotated[
        bool,
        typer.Option(
            "--stats",
            help="Add exact hypergeometric tails, rank AUC and a paired permutation test.",
        ),
    ] = False,
) -> None:
    """Validation rate of the top-k ranked peptides, against a baseline."""
    results = read_results(results_tsv)
    rankings: dict[str, dict[str, float]] = {"neoantigene": _ranking_from_tsv(ranked_tsv)}
    if baseline_tsv is not None:
        rankings[baseline_tsv.stem] = _ranking_from_tsv(baseline_tsv)
        summaries = compare_rankings(results, rankings, k)
    else:
        summaries = [validation_rate_at_k(results, rankings["neoantigene"], k)]

    for summary in summaries:
        typer.echo(summary.describe())

    if stats:
        baseline_name = baseline_tsv.stem if baseline_tsv is not None else None
        _echo_significance(results, rankings, k, baseline=baseline_name)


@app.command()
def compare(
    results_tsv: Annotated[Path, typer.Argument(help="Ground-truth assay results.", exists=True)],
    method: Annotated[
        list[str],
        typer.Option(
            "--method",
            "-m",
            help="NAME=path.tsv ranking to evaluate. Repeatable. Any tool's output.",
        ),
    ],
    baseline: Annotated[
        str | None,
        typer.Option("--baseline", help="Which --method name to run paired tests against."),
    ] = None,
    k: Annotated[int, typer.Option("--k", min=1)] = 20,
    alpha: Annotated[float, typer.Option("--alpha", min=1e-6, max=0.5)] = 0.05,
    rounds: Annotated[int, typer.Option("--rounds", min=100)] = 10_000,
    out: Annotated[
        Path | None, typer.Option("-o", "--out", help="Write a markdown report.")
    ] = None,
) -> None:
    """Audit any set of rankings against one label set.

    Accepts rankings from any tool, not just this one, and reports what the
    assayed pool could have detected before reporting what it did. A pool that
    cannot reject chance at any outcome says so at the top of the report.

        neoantigene compare labels.tsv \\
          -m ours=ranked.tsv -m netmhcpan=theirs.tsv --baseline netmhcpan
    """
    run_id = new_run_id()
    configure_logging(run_id)

    rankings: dict[str, dict[str, float]] = {}
    sources = [SourceFile(role="labels", path=str(results_tsv), digest=file_digest(results_tsv))]
    for entry in method:
        name, _, raw_path = entry.partition("=")
        if not name or not raw_path:
            raise typer.BadParameter(f"expected NAME=path.tsv, got {entry!r}")
        path = Path(raw_path)
        if not path.exists():
            raise typer.BadParameter(f"{name}: {path} does not exist")
        rankings[name] = _ranking_from_tsv(path)
        sources.append(SourceFile(role=f"ranking:{name}", path=str(path), digest=file_digest(path)))

    if baseline is not None and baseline not in rankings:
        raise typer.BadParameter(
            f"--baseline {baseline!r} is not one of the given methods: {sorted(rankings)}"
        )

    report = build_report(
        read_results(results_tsv),
        rankings,
        run_id=run_id,
        package_version=__version__,
        k=k,
        alpha=alpha,
        baseline=baseline,
        sources=sources,
        rounds=rounds,
    )

    if out is not None:
        out.parent.mkdir(parents=True, exist_ok=True)
        out.write_text(report.to_markdown())
        typer.echo(f"report -> {out}")
    else:
        typer.echo(report.to_markdown())

    # Exit 2 when nothing in the report could sustain a claim, so a pipeline or
    # a CI job cannot mistake an inconclusive benchmark for a passing one.
    if not report.supports_a_claim:
        raise typer.Exit(code=2)


@app.command()
def benchmark(
    manifest: ManifestArg,
    results_tsv: Annotated[Path, typer.Argument(help="Completed assay results.", exists=True)],
    config_path: ConfigOpt = None,
    backend: Annotated[str | None, typer.Option("--backend")] = None,
    k: Annotated[int, typer.Option("--k", min=1)] = 20,
    allow_null_backend: Annotated[bool, typer.Option("--allow-null-backend")] = False,
    stats: Annotated[
        bool,
        typer.Option(
            "--stats",
            help="Add exact hypergeometric tails, rank AUC and a paired permutation test.",
        ),
    ] = False,
) -> None:
    """Re-rank one sample with every built-in baseline and compare hit rates.

    Every method is scored over the intersection of peptides they all
    nominated, so none can win by picking an easier set.
    """
    run_id = new_run_id()
    configure_logging(run_id)

    sample = Sample.load(manifest)
    config = PipelineConfig.load(config_path)
    if backend is not None:
        config.presentation.backend = backend
    if config.presentation.backend in DEVELOPMENT_BACKENDS and not allow_null_backend:
        raise typer.BadParameter("pass --allow-null-backend to benchmark with a dev backend")

    result = run_pipeline(sample, config, run_id=run_id)
    results = read_results(results_tsv)

    prefix = f"{sample.sample_id}|"
    rankings: dict[str, dict[str, float]] = {
        "neoantigene": {f"{prefix}{item.key}": item.score for item in result.all_scored}
    }
    for name, ranker in BASELINES.items():
        rankings[name] = {
            f"{prefix}{key}": value for key, value in ranker(result.all_scored).items()
        }

    for summary in compare_rankings(results, rankings, k, baseline="binding_only"):
        typer.echo(summary.describe())

    if stats:
        _echo_significance(results, rankings, k, baseline="binding_only")


@app.command()
def refit(
    results_tsv: Annotated[Path, typer.Argument(help="Completed assay results.", exists=True)],
    features: Annotated[
        list[Path], typer.Option("--features", "-f", help="features.json files.", exists=True)
    ],
    config_path: ConfigOpt = None,
    out: Annotated[Path | None, typer.Option("--out", help="Write refitted config here.")] = None,
) -> None:
    """Refit ranking weights on returned assay labels."""
    config = PipelineConfig.load(config_path)
    dataset = build_training_set(read_results(results_tsv), load_feature_records(features))

    try:
        report = fit_weights(dataset)
    except InsufficientData as exc:
        typer.echo(f"cannot refit: {exc}")
        raise typer.Exit(code=1) from exc

    typer.echo(report.describe())
    typer.echo("")
    for line in compare_to_prior(config.weights, report.weights):
        typer.echo(line)

    if out is not None:
        config.weights = report.weights
        config.dump(out)
        typer.echo(f"\nwrote refitted config to {out}")


@app.command()
def explain(
    features_json: Annotated[Path, typer.Argument(help="features.json.", exists=True)],
    peptide: Annotated[str, typer.Argument(help="Mutant peptide sequence.")],
    config_path: ConfigOpt = None,
) -> None:
    """Break a candidate's score into per-feature contributions."""
    config = PipelineConfig.load(config_path)
    entries = [e for e in json.loads(features_json.read_text()) if e["peptide"] == peptide.upper()]
    if not entries:
        typer.echo(f"no candidate with peptide {peptide}")
        raise typer.Exit(code=1)

    for entry in sorted(entries, key=lambda e: -e["score"]):
        typer.echo(f"\n{peptide} / {entry['allele']}  score={entry['score']:.4f}")
        if entry["gate_failures"]:
            typer.echo(f"  gated: {', '.join(entry['gate_failures'])}")
        typer.echo(f"  {'bias':<18} {config.weights.bias:+.3f}")
        ordered = sorted(
            contributions(entry["features"], config.weights).items(),
            key=lambda kv: -abs(kv[1]),
        )
        for name, value in ordered:
            typer.echo(f"  {name:<18} {value:+.3f}   (x={entry['features'][name]:.3f})")


@app.command()
def catalog(
    hla: Annotated[str | None, typer.Option("--hla", help="Comma-separated genotype.")] = None,
) -> None:
    """Show the PDAC shared-antigen catalog."""
    eligible: set[str] = set()
    if hla is not None:
        eligible = {a.key for a in pdac.eligible_by_hla(parse_hla_string(hla))}

    for antigen in pdac.CATALOG:
        frequency = f"{antigen.frequency:.0%}" if antigen.frequency else "-"
        marker = "*" if antigen.key in eligible else " "
        restrictions = ", ".join(antigen.restricting_alleles) or "no published restriction"
        typer.echo(f"{marker} {antigen.key:<16} {frequency:>5}  {restrictions}")
    if hla is not None:
        typer.echo("\n* has a published epitope for one of the supplied alleles")


@app.command("fetch-proteome")
def fetch_proteome_command(
    release: Annotated[str, typer.Option("--release", help="Ensembl release.")] = ENSEMBL_RELEASE,
    out: Annotated[Path | None, typer.Option("--out", help="Destination file.")] = None,
    force: Annotated[bool, typer.Option("--force", help="Re-download if present.")] = False,
) -> None:
    """Download the Ensembl human proteome the self-peptide gate needs."""
    try:
        path = fetch_proteome(release=release, dest=out, force=force)
    except ReferenceDownloadError as exc:
        typer.echo(f"download failed: {exc}")
        raise typer.Exit(code=1) from exc

    report = provenance(path)
    typer.echo(f"proteome : {path}")
    typer.echo(f"proteins : {report.proteins}")
    if report.caveat:
        typer.echo(f"warning  : {report.caveat}")
    else:
        typer.echo("point processed.proteome_fasta at that path")


@app.command("init-config")
def init_config(
    out: Annotated[Path, typer.Option("--out", "-o")] = Path("config/default.yaml"),
) -> None:
    """Write a config file populated with the default priors."""
    typer.echo(f"wrote {PipelineConfig().dump(out)}")


def _read_tsv(path: Path) -> list[dict[str, Any]]:
    with open(path) as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def _echo_significance(
    results: Sequence[AssayResult],
    rankings: Mapping[str, Mapping[str, float]],
    k: int,
    baseline: str | None = None,
) -> None:
    """Report whether the retrieval numbers above could have happened anyway.

    Printed under `--stats` rather than always, because it costs permutations
    and because the point estimates are what a lab reads. The paired test is
    the one that matters: it asks whether the full model beats its own
    strongest input feature on the same peptides.
    """
    typer.echo("")
    for name, ranking in rankings.items():
        typer.echo(top_k_significance(results, ranking, k, method=name).describe())
    typer.echo("")
    for name, ranking in rankings.items():
        typer.echo(stratified_auc(results, ranking, method=name).describe())

    reference = baseline if baseline in rankings else None
    if reference is None or reference == "neoantigene":
        return
    typer.echo("")
    typer.echo(
        paired_comparison(
            results,
            rankings["neoantigene"],
            rankings[reference],
            k=k,
            method_b=reference,
        ).describe()
    )


def _ranking_from_tsv(path: Path) -> dict[str, float]:
    ranking: dict[str, float] = {}
    for row in _read_tsv(path):
        peptide = row.get("mutant_peptide") or row.get("peptide")
        score = row.get("score")
        if not peptide or not score:
            continue
        allele = row.get("allele") or "*"
        ranking[f"{row['sample_id']}|{peptide}|{allele}"] = float(score)
    return ranking


if __name__ == "__main__":
    app()
