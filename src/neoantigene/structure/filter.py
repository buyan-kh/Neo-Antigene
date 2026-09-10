"""pMHC structure filter interface (not yet implemented).

Intended contract once wired to AlphaFold 3 or an equivalent pMHC predictor:

  1. Take the top-N shortlist only, after sequence-based ranking.
  2. Fold each peptide-MHC pair; read out interface confidence and whether the
     mutated residue points toward the TCR rather than into the groove.
  3. Return a verdict that can demote or annotate a candidate. It must not be
     able to promote one, because doing so would let a structural artifact
     override the assay-fit sequence model.
"""

from __future__ import annotations

from collections.abc import Sequence

from pydantic import BaseModel, ConfigDict

from ..models import ScoredCandidate


class StructureVerdict(BaseModel):
    model_config = ConfigDict(frozen=True)

    candidate_id: str
    allele: str
    interface_confidence: float | None = None
    mutation_solvent_exposed: bool | None = None
    demote: bool = False
    note: str = ""


class StructureFilter:
    """Placeholder for the AlphaFold 3 pMHC secondary filter."""

    def __init__(self, model: str = "alphafold3") -> None:
        self.model = model

    def evaluate(self, shortlist: Sequence[ScoredCandidate]) -> list[StructureVerdict]:
        raise NotImplementedError(
            "pMHC structure filtering is not implemented. It is a post-ranking "
            "step and is not required for a ranking run."
        )
