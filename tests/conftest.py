from pathlib import Path

import pytest

from neoantigene.config import (
    ExpressionConfig,
    OutputConfig,
    PeptideConfig,
    PipelineConfig,
    PresentationConfig,
    VariantFilterConfig,
)
from neoantigene.io.manifest import Sample
from neoantigene.peptides.proteome import ProteomeIndex
from neoantigene.pipeline import RunResult
from neoantigene.pipeline import run as run_pipeline
from neoantigene.presentation.null import NullBackend

from .fixtures import SyntheticCohort, build_cohort

EXAMPLES = Path(__file__).resolve().parents[1] / "data" / "examples"


@pytest.fixture(scope="session")
def examples_dir() -> Path:
    return EXAMPLES


@pytest.fixture(scope="session")
def proteome() -> ProteomeIndex:
    return ProteomeIndex.from_fasta(EXAMPLES / "proteome.mini.fa")


@pytest.fixture
def dev_config() -> PipelineConfig:
    """Config for pipeline tests: deterministic backend, default gates left on."""
    config = PipelineConfig()
    config.presentation = PresentationConfig(backend="null")
    return config


@pytest.fixture
def example_sample(examples_dir: Path) -> Sample:
    return Sample.load(examples_dir / "sample.yaml")


@pytest.fixture(scope="session")
def cohort(tmp_path_factory: pytest.TempPathFactory) -> SyntheticCohort:
    """A synthetic multi-patient cohort written to a temporary directory."""
    return build_cohort(tmp_path_factory.mktemp("cohort"), variants_per_patient=30)


@pytest.fixture(scope="session")
def benchmark_config() -> PipelineConfig:
    """Config for ranking-quality benchmarks.

    Hard gates are disabled on purpose. The question a benchmark answers is
    "which ranking function orders candidates better", and a gate that removes
    a candidate before either method sees it only obscures that.
    """
    return PipelineConfig(
        variant_filters=VariantFilterConfig(min_ccf=0.0, min_dna_vaf=0.0, min_tumor_depth=0),
        # The expression gate has to come off too, or low-expression candidates
        # are removed before either ranker can be judged on them.
        expression=ExpressionConfig(min_tpm=0.0, min_rna_vaf=None),
        presentation=PresentationConfig(backend="null", max_affinity_percentile=100.0),
        peptides=PeptideConfig(lengths=[9], drop_self_matching=False),
        output=OutputConfig(top_n=2000, max_per_variant=1000),
    )


@pytest.fixture(scope="session")
def benchmark_run(cohort: SyntheticCohort, benchmark_config: PipelineConfig) -> RunResult:
    sample = Sample.load(cohort.patient(0).manifest_path)
    return run_pipeline(sample, benchmark_config, backend=NullBackend(), run_id="test-run")
