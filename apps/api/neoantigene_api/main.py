"""HTTP surface over the `neoantigene` package.

A thin wrapper by design: every endpoint parses a request, calls one or two
library functions, and serializes the result. No ranking, scoring, peptide
generation or metric arithmetic happens here, and none of it is reimplemented
in the frontend either.
"""

from __future__ import annotations

import json
from collections.abc import AsyncIterator
from contextlib import asynccontextmanager
from pathlib import Path
from typing import Annotated, Any, Final

from fastapi import Depends, FastAPI, File, HTTPException, Query, Request, UploadFile, status
from fastapi.middleware.cors import CORSMiddleware
from fastapi.responses import FileResponse
from pydantic import ValidationError

from neoantigene import __version__
from neoantigene.assays.metrics import (
    compare_rankings,
    compare_retrieval,
    labelled,
    positive_ranks,
)
from neoantigene.assays.schema import AssayCall, AssayResult, read_results
from neoantigene.config import PipelineConfig
from neoantigene.io.manifest import Sample
from neoantigene.io.writers import COLUMNS, to_records
from neoantigene.models import ScoredCandidate
from neoantigene.peptides import reference
from neoantigene.presentation.registry import DEVELOPMENT_BACKENDS, available
from neoantigene.scoring.baselines import BASELINES
from neoantigene.scoring.rank import contributions, logit, weight_vector

from .runs import Run, RunNotFound, RunStore, development_warning
from .schemas import (
    BackendInfo,
    BenchmarkReport,
    Contribution,
    Explanation,
    FieldError,
    LabelSource,
    ReferenceInfo,
    ResultsPage,
    RunCreated,
    RunRequest,
    RunState,
    RunWarnings,
    UploadedFile,
    Workspace,
)

REPO_ROOT: Final[Path] = Path(__file__).resolve().parents[3]
EXAMPLE_MANIFEST: Final[Path] = REPO_ROOT / "data" / "examples" / "sample.yaml"
EXAMPLE_PROTEOME: Final[Path] = REPO_ROOT / "data" / "examples" / "proteome.mini.fa"
DEFAULT_CONFIG: Final[Path] = REPO_ROOT / "config" / "default.yaml"

#: Upload ceiling. Generous for a variant TSV, far below a real proteome.
MAX_UPLOAD_BYTES: Final[int] = 64 * 1024 * 1024

#: Starlette deprecated its own names for these two; the numbers are stable.
UNPROCESSABLE: Final[int] = 422
TOO_LARGE: Final[int] = 413


@asynccontextmanager
async def lifespan(application: FastAPI) -> AsyncIterator[None]:
    """One run store per application lifespan, torn down with it."""
    application.state.store = RunStore()
    try:
        yield
    finally:
        application.state.store.shutdown()
        application.state.store = None


app = FastAPI(
    title="Neo Antigene API",
    version=__version__,
    summary="Ranks mutant peptides for personalized cancer vaccines.",
    lifespan=lifespan,
)

app.add_middleware(
    CORSMiddleware,
    allow_origins=["http://localhost:3000", "http://127.0.0.1:3000"],
    allow_methods=["*"],
    allow_headers=["*"],
)


def store(request: Request) -> RunStore:
    existing = getattr(request.app.state, "store", None)
    if existing is None:
        # Only reachable if the app is served without its lifespan running.
        request.app.state.store = existing = RunStore()
    return existing


StoreDep = Annotated[RunStore, Depends(store)]


def _run_or_404(runs: RunStore, run_id: str) -> Run:
    try:
        return runs.get(run_id)
    except RunNotFound as exc:
        raise HTTPException(status.HTTP_404_NOT_FOUND, f"unknown run {run_id}") from exc


@app.get("/api/health")
def health() -> dict[str, str]:
    return {"status": "ok", "version": __version__}


@app.get("/api/backends", response_model=BackendInfo)
def backends() -> BackendInfo:
    """Which predictors are installed, and whether MHCflurry has its weights."""
    ready, detail = _mhcflurry_state()
    names = available()
    default = "mhcflurry" if "mhcflurry" in names and ready else names[0]
    return BackendInfo(
        available=names,
        development=sorted(DEVELOPMENT_BACKENDS),
        default=default,
        mhcflurry_ready=ready,
        mhcflurry_detail=detail,
    )


def _mhcflurry_state() -> tuple[bool, str]:
    """Distinguish "not installed" from "installed but no weights downloaded"."""
    try:
        from mhcflurry import Class1PresentationPredictor
    except ImportError:
        return False, "mhcflurry is not installed; `uv sync --extra presentation`"
    try:
        Class1PresentationPredictor.load()
    except Exception as exc:
        return False, (
            f"mhcflurry is installed but its models did not load ({exc}); run "
            "`uv run mhcflurry-downloads fetch models_class1_presentation`"
        )
    return True, "models_class1_presentation loaded"


@app.get("/api/reference", response_model=ReferenceInfo)
def reference_info() -> ReferenceInfo:
    """Reference proteome state, so the UI can warn before a run rather than after."""
    cached = reference.proteome_path()
    info = ReferenceInfo(
        bundled_stub_path=str(EXAMPLE_PROTEOME.relative_to(REPO_ROOT)),
        bundled_stub_proteins=reference.count_proteins(EXAMPLE_PROTEOME),
        cached_release=reference.ENSEMBL_RELEASE,
        fetch_command="uv run neoantigene fetch-proteome",
    )
    if cached.exists():
        report = reference.provenance(cached)
        info.cached_path = str(cached)
        info.cached_proteins = report.proteins
        info.cached_is_stub = report.is_stub
    return info


@app.get("/api/config", response_model=PipelineConfig)
def default_config() -> PipelineConfig:
    path = DEFAULT_CONFIG if DEFAULT_CONFIG.exists() else None
    return PipelineConfig.load(path)


@app.get("/api/weights")
def weights() -> dict[str, float]:
    """The scored feature weights, which are priors rather than fitted values."""
    config = default_config()
    return {"bias": config.weights.bias, **weight_vector(config.weights)}


@app.get("/api/example")
def example_sample() -> dict[str, Any]:
    """The bundled example manifest, for the one-click path."""
    if not EXAMPLE_MANIFEST.exists():
        raise HTTPException(status.HTTP_404_NOT_FOUND, "bundled example is missing")
    sample = Sample.load(EXAMPLE_MANIFEST)
    return {
        "path": str(EXAMPLE_MANIFEST.relative_to(REPO_ROOT)),
        "sample_id": sample.sample_id,
        "cancer_type": sample.cancer_type,
        "tumor_purity": sample.tumor_purity,
        "hla": list(sample.hla),
        "variant_tsv": _name(sample.processed.variant_tsv),
        "somatic_vcf": _name(sample.processed.somatic_vcf),
        "expression_tsv": _name(sample.processed.expression_tsv),
        "normal_expression_tsv": _name(sample.processed.normal_expression_tsv),
        "proteome_fasta": _name(sample.processed.proteome_fasta),
    }


def _name(path: Path | None) -> str | None:
    return None if path is None else path.name


@app.post("/api/uploads", response_model=Workspace)
async def upload(
    runs: StoreDep,
    files: Annotated[list[UploadFile], File()],
    workspace: Annotated[str | None, Query()] = None,
) -> Workspace:
    """Stash uploaded inputs in a workspace directory keyed by id."""
    workspace_id, directory = (
        (workspace, runs.workspace(workspace)) if workspace else runs.new_workspace()
    )
    stored: list[UploadedFile] = []
    for item in files:
        name = Path(item.filename or "unnamed").name
        payload = await item.read()
        if len(payload) > MAX_UPLOAD_BYTES:
            raise HTTPException(
                TOO_LARGE,
                f"{name} is larger than {MAX_UPLOAD_BYTES // (1024 * 1024)} MB",
            )
        (directory / name).write_bytes(payload)
        stored.append(UploadedFile(name=name, bytes=len(payload)))

    existing = [
        UploadedFile(name=p.name, bytes=p.stat().st_size)
        for p in sorted(directory.iterdir())
        if p.is_file()
    ]
    return Workspace(id=workspace_id, files=existing or stored)


@app.post("/api/runs", response_model=RunCreated, status_code=status.HTTP_202_ACCEPTED)
def start_run(runs: StoreDep, request: RunRequest) -> RunCreated:
    """Begin a ranking run and return its id. Poll `/api/runs/{id}` for progress."""
    config = default_config()
    if request.backend:
        config.presentation.backend = request.backend
    if request.top_n:
        config.output.top_n = request.top_n

    if config.presentation.backend not in available():
        raise HTTPException(
            status.HTTP_400_BAD_REQUEST,
            f"unknown backend {config.presentation.backend!r}; available: {available()}",
        )

    is_development, development_detail = development_warning(config.presentation.backend)
    if is_development and not request.acknowledge_development_backend:
        # Same refusal as the CLI's --allow-null-backend guard.
        raise HTTPException(
            status.HTTP_400_BAD_REQUEST,
            f"the {config.presentation.backend!r} backend produces meaningless scores and is "
            "for development only; set acknowledge_development_backend to proceed",
        )

    sample = _resolve_sample(runs, request)
    warnings = RunWarnings(
        development_backend=is_development,
        development_detail=development_detail,
    )
    return RunCreated(run_id=runs.submit(sample, config, warnings))


def _resolve_sample(runs: RunStore, request: RunRequest) -> Sample:
    if request.source == "example":
        if not EXAMPLE_MANIFEST.exists():
            raise HTTPException(status.HTTP_404_NOT_FOUND, "bundled example is missing")
        return Sample.load(EXAMPLE_MANIFEST)

    if request.sample is None:
        raise HTTPException(UNPROCESSABLE, "sample is required for uploads")
    if not request.workspace:
        raise HTTPException(UNPROCESSABLE, "workspace is required")

    directory = runs.workspace(request.workspace)
    payload = _manifest_payload(request.sample, directory)
    try:
        # Reuse the library's own manifest validation rather than duplicating it.
        return Sample.model_validate(payload)
    except ValidationError as exc:
        raise HTTPException(
            UNPROCESSABLE,
            detail=[e.model_dump() for e in _field_errors(exc)],
        ) from exc


def _manifest_payload(sample: Any, directory: Path) -> dict[str, Any]:
    """Build a manifest dict, resolving named uploads to workspace paths."""

    def resolve(name: str | None) -> str | None:
        if not name:
            return None
        candidate = directory / Path(name).name
        if not candidate.exists():
            raise HTTPException(
                UNPROCESSABLE,
                detail=[
                    FieldError(
                        field=name, message=f"{Path(name).name} was not uploaded to this workspace"
                    ).model_dump()
                ],
            )
        return str(candidate)

    proteome = sample.proteome
    if proteome == "stub":
        proteome_path: str | None = str(EXAMPLE_PROTEOME)
    elif proteome == "ensembl":
        cached = reference.proteome_path()
        if not cached.exists():
            raise HTTPException(
                UNPROCESSABLE,
                detail=[
                    FieldError(
                        field="proteome",
                        message=(
                            "the Ensembl proteome is not downloaded yet; run "
                            "`uv run neoantigene fetch-proteome` or choose the bundled stub"
                        ),
                    ).model_dump()
                ],
            )
        proteome_path = str(cached)
    else:
        proteome_path = resolve(proteome)

    return {
        "sample_id": sample.sample_id,
        "cancer_type": sample.cancer_type,
        "tumor_purity": sample.tumor_purity,
        "hla": sample.hla,
        "processed": {
            "somatic_vcf": resolve(sample.somatic_vcf),
            "variant_tsv": resolve(sample.variant_tsv),
            "expression_tsv": resolve(sample.expression_tsv),
            "normal_expression_tsv": resolve(sample.normal_expression_tsv),
            "proteome_fasta": proteome_path,
            "tumor_sample_name": sample.tumor_sample_name,
        },
    }


def _field_errors(exc: ValidationError) -> list[FieldError]:
    """Map pydantic locations onto the form fields the frontend rendered."""
    aliases = {
        "processed.somatic_vcf": "somatic_vcf",
        "processed.variant_tsv": "variant_tsv",
        "processed.expression_tsv": "expression_tsv",
        "processed.normal_expression_tsv": "normal_expression_tsv",
        "processed.proteome_fasta": "proteome",
        "processed.tumor_sample_name": "tumor_sample_name",
        "processed": "variant_tsv",
    }
    errors: list[FieldError] = []
    for error in exc.errors():
        location = ".".join(str(part) for part in error["loc"] if part != "__root__")
        field = aliases.get(location, location or "sample_id")
        errors.append(FieldError(field=field, message=error["msg"]))
    return errors


@app.get("/api/runs", response_model=list[RunState])
def list_runs(runs: StoreDep) -> list[RunState]:
    return runs.list_runs()


@app.get("/api/runs/{run_id}", response_model=RunState)
def run_state(runs: StoreDep, run_id: str) -> RunState:
    return _run_or_404(runs, run_id).state


@app.get("/api/runs/{run_id}/results", response_model=ResultsPage)
def results(
    runs: StoreDep,
    run_id: str,
    include_gated: Annotated[bool, Query()] = True,
) -> ResultsPage:
    """The ranked table, built with the same writer the TSV uses."""
    run = _run_or_404(runs, run_id)
    if run.result is None:
        raise HTTPException(status.HTTP_409_CONFLICT, f"run {run_id} has no results yet")

    sample_id = run.result.sample.sample_id
    shortlist = list(run.result.ranked)
    if include_gated:
        shortlist = _with_gated(shortlist, run.result.all_scored, run.result.config.output.top_n)

    rows = to_records(shortlist, sample_id, run_id=run_id)
    return ResultsPage(
        run_id=run_id,
        sample_id=sample_id,
        columns=list(COLUMNS),
        rows=rows,
        total=len(rows),
        genes=sorted({str(r["gene"]) for r in rows if r.get("gene")}),
        alleles=sorted({str(r["allele"]) for r in rows if r.get("allele")}),
        gate_failures=sorted(
            {gate for r in rows for gate in str(r.get("gate_failures") or "").split(";") if gate}
        ),
    )


def _with_gated(
    ranked: list[ScoredCandidate],
    all_scored: list[ScoredCandidate],
    top_n: int,
) -> list[ScoredCandidate]:
    """Append the best gated candidates so the UI can show what was excluded.

    Kept after the shortlist and never renumbered into it: a gated peptide is
    not a ranked recommendation, it is an explanation of one.
    """
    shortlisted = {item.key for item in ranked}
    gated = sorted(
        (item for item in all_scored if not item.passed and item.key not in shortlisted),
        key=lambda item: (-item.score, item.candidate.id, item.allele),
    )
    return ranked + gated[:top_n]


@app.get("/api/runs/{run_id}/explain", response_model=Explanation)
def explain(
    runs: StoreDep,
    run_id: str,
    peptide: Annotated[str, Query(min_length=1)],
    allele: Annotated[str, Query(min_length=1)],
) -> Explanation:
    """Per-feature score decomposition, from `scoring.rank.contributions`."""
    run = _run_or_404(runs, run_id)
    if run.result is None:
        raise HTTPException(status.HTTP_409_CONFLICT, f"run {run_id} has no results yet")

    key = f"{peptide.upper()}|{allele}"
    match = next((item for item in run.result.all_scored if item.key == key), None)
    if match is None:
        raise HTTPException(status.HTTP_404_NOT_FOUND, f"no candidate {key} in run {run_id}")

    config = run.result.config
    weights = weight_vector(config.weights)
    breakdown = contributions(match.features, config.weights)
    return Explanation(
        peptide=match.candidate.mutant_peptide,
        wildtype_peptide=match.candidate.wildtype_peptide,
        allele=match.allele,
        gene=match.candidate.variant.gene,
        hgvsp=match.candidate.variant.hgvsp_short,
        score=match.score,
        bias=config.weights.bias,
        logit=logit(match.features, config.weights),
        mutation_position=match.candidate.mutation_position,
        gate_failures=[gate.value for gate in match.gate_failures],
        contributions=[
            Contribution(
                feature=name,
                weight=weights[name],
                value=match.features.get(name, 0.0),
                contribution=value,
            )
            for name, value in sorted(breakdown.items(), key=lambda kv: -abs(kv[1]))
        ],
    )


@app.get("/api/runs/{run_id}/download/{artifact}")
def download(runs: StoreDep, run_id: str, artifact: str) -> FileResponse:
    """Serve the files the library's writers already produced."""
    run = _run_or_404(runs, run_id)
    paths = {"ranked.tsv": run.ranked_tsv, "features.json": run.features_json}
    if artifact not in paths:
        raise HTTPException(status.HTTP_404_NOT_FOUND, f"unknown artifact {artifact}")
    path = paths[artifact]
    if path is None or not path.exists():
        raise HTTPException(status.HTTP_409_CONFLICT, f"run {run_id} has not written {artifact}")
    return FileResponse(
        path,
        filename=path.name,
        media_type="text/tab-separated-values" if artifact.endswith(".tsv") else "application/json",
    )


@app.post("/api/runs/{run_id}/benchmark", response_model=BenchmarkReport)
async def benchmark(
    runs: StoreDep,
    run_id: str,
    labels: Annotated[UploadFile, File()],
    k: Annotated[int, Query(ge=1)] = 20,
) -> BenchmarkReport:
    """Score this run against the baselines using an uploaded assay-label TSV.

    Labels arrive only here, after the run has finished. Nothing in the ranking
    path can see them, which is what makes the comparison worth reading.
    """
    run = _run_or_404(runs, run_id)
    if run.result is None:
        raise HTTPException(status.HTTP_409_CONFLICT, f"run {run_id} has no results yet")

    filename = Path(labels.filename or "labels.tsv").name
    destination = runs.workspace(f"labels-{run_id}") / filename
    destination.write_bytes(await labels.read())

    try:
        results_list = read_results(destination)
    except (ValidationError, ValueError, KeyError) as exc:
        raise HTTPException(
            UNPROCESSABLE,
            f"could not read {filename} as an assay result TSV: {exc}",
        ) from exc

    if not results_list:
        raise HTTPException(UNPROCESSABLE, f"{filename} has no rows")

    sample_id = run.result.sample.sample_id
    rankings = _rankings(run, sample_id)
    pool = [r for r in labelled(results_list) if r.key in rankings["neoantigene"]]

    if not pool:
        raise HTTPException(
            UNPROCESSABLE,
            (
                f"none of the {len(results_list)} labelled peptides in {filename} match a "
                f"scored peptide-HLA pair in this run. Labels must use sample_id "
                f"{sample_id!r} and the same allele notation."
            ),
        )

    return BenchmarkReport(
        run_id=run_id,
        source=_label_source(filename, results_list, len(pool)),
        retrieval=compare_retrieval(results_list, rankings),
        validation=compare_rankings(results_list, rankings, k, baseline="binding_only"),
        k=k,
        positive_ranks={
            name: positive_ranks(results_list, ranking) for name, ranking in rankings.items()
        },
        pool_size=len(pool),
    )


def _rankings(run: Run, sample_id: str) -> dict[str, dict[str, float]]:
    """Neo Antigene plus every built-in comparator, keyed as the metrics expect."""
    scored = run.scored()
    prefix = f"{sample_id}|"
    rankings: dict[str, dict[str, float]] = {
        "neoantigene": {f"{prefix}{item.key}": item.score for item in scored}
    }
    for name, ranker in BASELINES.items():
        rankings[name] = {f"{prefix}{key}": value for key, value in ranker(scored).items()}
    return rankings


def _label_source(filename: str, results_list: list[AssayResult], matched: int) -> LabelSource:
    assays: dict[str, int] = {}
    for item in results_list:
        assays[item.assay.value] = assays.get(item.assay.value, 0) + 1
    return LabelSource(
        filename=filename,
        assayed=len(results_list),
        positives=sum(1 for r in results_list if r.call is AssayCall.POSITIVE),
        negatives=sum(1 for r in results_list if r.call is AssayCall.NEGATIVE),
        indeterminate=sum(1 for r in results_list if r.call is AssayCall.INDETERMINATE),
        assays=assays,
        matched_to_run=matched,
    )


@app.get("/api/runs/{run_id}/features")
def features(runs: StoreDep, run_id: str) -> list[dict[str, Any]]:
    """Raw feature records, the same payload `features.json` contains."""
    run = _run_or_404(runs, run_id)
    if run.features_json is None or not run.features_json.exists():
        raise HTTPException(status.HTTP_409_CONFLICT, f"run {run_id} has no features yet")
    payload: list[dict[str, Any]] = json.loads(run.features_json.read_text())
    return payload
