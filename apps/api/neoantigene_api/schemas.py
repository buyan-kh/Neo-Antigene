"""HTTP request and response shapes.

Every response model either wraps a library model directly or is built from one,
so the OpenAPI schema the frontend generates its types from stays tied to the
Python definitions rather than drifting from them.
"""

from __future__ import annotations

from datetime import datetime
from enum import StrEnum
from typing import Any, Literal

from pydantic import BaseModel, Field

from neoantigene.assays.metrics import RetrievalSummary, ValidationSummary
from neoantigene.config import PipelineConfig
from neoantigene.pipeline import RunReport

from .progress import ProgressEvent


class RunStatus(StrEnum):
    QUEUED = "queued"
    RUNNING = "running"
    SUCCEEDED = "succeeded"
    FAILED = "failed"


class BackendInfo(BaseModel):
    """Which presentation predictors this deployment can actually use."""

    available: list[str]
    development: list[str]
    default: str
    mhcflurry_ready: bool
    mhcflurry_detail: str


class ReferenceInfo(BaseModel):
    """State of the reference proteome, which gates two features and a hard gate."""

    bundled_stub_path: str
    bundled_stub_proteins: int
    cached_release: str
    cached_path: str | None = None
    cached_proteins: int | None = None
    cached_is_stub: bool | None = None
    fetch_command: str


class FieldError(BaseModel):
    """A manifest validation failure, addressed to one form field."""

    field: str
    message: str


class SampleInput(BaseModel):
    """Mirrors `data/examples/sample.yaml`.

    File fields name a file previously uploaded to a workspace, not a server
    path, so a request cannot reach outside its own upload directory.
    """

    sample_id: str = Field(min_length=1)
    cancer_type: str | None = None
    tumor_purity: float | None = None
    hla: list[str] = Field(default_factory=list)

    somatic_vcf: str | None = None
    variant_tsv: str | None = None
    expression_tsv: str | None = None
    normal_expression_tsv: str | None = None
    tumor_sample_name: str | None = None

    #: `stub` uses the bundled two-protein example, `ensembl` the fetched
    #: release, or name an uploaded FASTA.
    proteome: str = "ensembl"


class RunRequest(BaseModel):
    source: Literal["example", "upload"] = "example"
    workspace: str | None = None
    sample: SampleInput | None = None
    backend: str | None = None
    top_n: int | None = Field(default=None, ge=1)
    #: Mirrors the CLI's `--allow-null-backend` guard.
    acknowledge_development_backend: bool = False


class UploadedFile(BaseModel):
    name: str
    bytes: int


class Workspace(BaseModel):
    id: str
    files: list[UploadedFile]


class RunCreated(BaseModel):
    run_id: str


class RunWarnings(BaseModel):
    """Caveats that must be visible in any interface, not just the README."""

    proteome_is_stub: bool = False
    proteome_proteins: int = 0
    proteome_detail: str | None = None
    development_backend: bool = False
    development_detail: str | None = None
    expression_available: bool = True
    expression_detail: str | None = None


class RunState(BaseModel):
    run_id: str
    status: RunStatus
    sample_id: str | None = None
    backend: str | None = None
    created_at: datetime
    started_at: datetime | None = None
    finished_at: datetime | None = None
    events: list[ProgressEvent] = Field(default_factory=list)
    report: RunReport | None = None
    warnings: RunWarnings = Field(default_factory=RunWarnings)
    config: PipelineConfig | None = None
    error: str | None = None
    error_type: str | None = None
    shortlisted: int = 0


class ResultsPage(BaseModel):
    """The ranked table, in exactly the column order the TSV uses."""

    run_id: str
    sample_id: str
    columns: list[str]
    rows: list[dict[str, Any]]
    total: int
    genes: list[str]
    alleles: list[str]
    gate_failures: list[str]


class Contribution(BaseModel):
    feature: str
    weight: float
    value: float
    contribution: float


class Explanation(BaseModel):
    """The per-feature decomposition behind one candidate's score."""

    peptide: str
    wildtype_peptide: str | None
    allele: str
    gene: str
    hgvsp: str
    score: float
    bias: float
    logit: float
    mutation_position: int
    gate_failures: list[str]
    contributions: list[Contribution]


class BenchmarkRequest(BaseModel):
    run_id: str
    k: int = Field(default=20, ge=1)


class LabelSource(BaseModel):
    """Provenance of the loaded labels, shown next to every benchmark number."""

    filename: str
    assayed: int
    positives: int
    negatives: int
    indeterminate: int
    assays: dict[str, int]
    matched_to_run: int


class BenchmarkReport(BaseModel):
    run_id: str
    source: LabelSource
    retrieval: list[RetrievalSummary]
    validation: list[ValidationSummary]
    k: int
    positive_ranks: dict[str, list[int]]
    pool_size: int
