"""Pipeline configuration.

Scoring weights live in config rather than code because they are refit from
assay outcomes (see `neoantigene.learning.train`). The shipped defaults are
priors, not fitted values.

Config is discovered relative to the working directory, never from a path
baked into the package, so an installed wheel and a source checkout behave the
same. Both YAML and TOML are accepted.
"""

from __future__ import annotations

import tomllib
from pathlib import Path
from typing import Any

import yaml
from pydantic import BaseModel, ConfigDict, Field

#: Searched in order, relative to the working directory, when no explicit
#: config path is given.
DEFAULT_CONFIG_FILENAMES = (
    Path("neoantigene.yaml"),
    Path("neoantigene.toml"),
    Path("config/default.yaml"),
)


class ConfigModel(BaseModel):
    model_config = ConfigDict(extra="forbid", validate_assignment=True)


class VariantFilterConfig(ConfigModel):
    require_pass: bool = True
    min_tumor_depth: int = Field(default=20, ge=0)
    min_dna_vaf: float = Field(default=0.05, ge=0.0, le=1.0)
    max_population_af: float = Field(default=0.0001, ge=0.0, le=1.0)
    min_ccf: float = Field(default=0.6, ge=0.0, le=1.0)
    tumor_purity: float = Field(default=0.5, gt=0.0, le=1.0)


class ExpressionConfig(ConfigModel):
    min_tpm: float = Field(default=1.0, ge=0.0)
    min_rna_vaf: float | None = Field(default=0.02, ge=0.0, le=1.0)
    require_rna_evidence: bool = False
    normal_tpm_ceiling: float | None = Field(default=None, ge=0.0)


class PeptideConfig(ConfigModel):
    lengths: list[int] = Field(default_factory=lambda: [8, 9, 10, 11])
    flank_length: int = Field(default=10, ge=0, le=30)
    drop_self_matching: bool = True

    @property
    def sorted_lengths(self) -> list[int]:
        return sorted(set(self.lengths))


class PresentationConfig(ConfigModel):
    backend: str = "mhcflurry"
    max_affinity_percentile: float = Field(default=2.0, gt=0.0, le=100.0)
    min_presentation_score: float = Field(default=0.0, ge=0.0, le=1.0)
    score_wildtype: bool = True


class ScoringWeights(ConfigModel):
    """Coefficients of a logistic model over normalized features.

    Weights are scaled by how well evidenced each feature is, not by how
    interesting it sounds. Presentation, clonality and expression are close to
    necessary conditions for a T-cell response and carry most of the weight.
    The sequence-intrinsic terms are proxies with weak literature support, so
    they are kept small on purpose: given them comparable weight and they
    dominate the very top of the shortlist, which is the part a lab actually
    synthesizes. `tests/test_eval.py` guards against that regression.

    Deliberately absent: any "known driver gene" or literature-prominence term.
    Driver status is a proxy for publication history, not for T-cell response,
    and biases the shortlist toward antigens that are already known to fail.

    These are priors. Replace them with `neoantigene refit` output.
    """

    # Chosen so a middling candidate lands near 0.5 rather than saturating at
    # 0.99, which keeps scores readable and keeps the active-learning
    # uncertainty term from collapsing.
    bias: float = -3.0

    presentation: float = 2.5
    clonality: float = 1.5
    expression: float = 1.2

    agretopicity: float = 0.5
    tumor_selectivity: float = 0.4
    self_dissimilarity: float = 0.25
    wt_dissimilarity: float = 0.15
    mutation_exposure: float = 0.1

    #: Zero on purpose: Chowell (PNAS 2015) and TESLA (Cell 2020) disagree on
    #: the sign of this effect in tumor neoepitopes. See
    #: `scoring.immunogenicity.tcr_contact_hydrophobicity` and docs/CITATIONS.md.
    hydrophobicity: float = 0.0


class OutputConfig(ConfigModel):
    top_n: int = Field(default=50, gt=0)
    max_per_variant: int = Field(default=2, gt=0)

    #: A frameshift is allowed more of the shortlist than a substitution,
    #: because it has more independent epitopes to offer rather than more
    #: registers of the same one.
    #:
    #: `max_per_variant` exists to stop one strong variant filling the list
    #: with overlapping windows of a single hypothesis. Downstream of a
    #: frameshift that reasoning inverts: Roudko et al. (2020) showed single
    #: recurrent MSI frameshifts producing several independently immunogenic
    #: epitopes from the novel tail, so capping at two discards distinct
    #: hypotheses rather than redundant ones.
    #:
    #: 6 is a deliberate choice, not a fitted one. The evidence says "several",
    #: not a number. Lower it to `max_per_variant` to recover the previous
    #: behaviour.
    max_per_frameshift_variant: int = Field(default=6, gt=0)

    include_failed: bool = False


class PipelineConfig(ConfigModel):
    variant_filters: VariantFilterConfig = Field(default_factory=VariantFilterConfig)
    expression: ExpressionConfig = Field(default_factory=ExpressionConfig)
    peptides: PeptideConfig = Field(default_factory=PeptideConfig)
    presentation: PresentationConfig = Field(default_factory=PresentationConfig)
    weights: ScoringWeights = Field(default_factory=ScoringWeights)
    output: OutputConfig = Field(default_factory=OutputConfig)

    @classmethod
    def load(cls, path: Path | None = None, search_from: Path | None = None) -> PipelineConfig:
        resolved = path if path is not None else discover(search_from)
        if resolved is None:
            return cls()
        return cls.model_validate(read_config_file(resolved))

    def dump(self, path: Path) -> Path:
        path.parent.mkdir(parents=True, exist_ok=True)
        payload = self.model_dump(mode="json")
        with open(path, "w") as handle:
            yaml.safe_dump(payload, handle, sort_keys=False)
        return path


def discover(search_from: Path | None = None) -> Path | None:
    base = search_from or Path.cwd()
    for candidate in DEFAULT_CONFIG_FILENAMES:
        path = base / candidate
        if path.is_file():
            return path
    return None


def read_config_file(path: Path) -> dict[str, Any]:
    if path.suffix == ".toml":
        with open(path, "rb") as binary:
            return tomllib.load(binary)
    with open(path) as handle:
        loaded = yaml.safe_load(handle)
    if loaded is None:
        return {}
    if not isinstance(loaded, dict):
        raise ValueError(
            f"{path}: expected a mapping at the top level, got {type(loaded).__name__}"
        )
    return loaded
