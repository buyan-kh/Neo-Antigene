from __future__ import annotations

from abc import ABC, abstractmethod
from collections.abc import Sequence

from ..models import MAX_PEPTIDE_LENGTH, MIN_PEPTIDE_LENGTH, PresentationCall


class PresentationBackend(ABC):
    """Predicts peptide-MHC presentation for every peptide x allele pair.

    Backends are swappable so that the ranking layer never depends on a
    specific predictor, and so that a retrained in-house model can be
    benchmarked against the incumbent on the same assay-labelled set.
    """

    name: str = "base"

    @abstractmethod
    def predict(
        self,
        peptides: Sequence[str],
        alleles: Sequence[str],
        n_flanks: Sequence[str] | None = None,
        c_flanks: Sequence[str] | None = None,
    ) -> list[PresentationCall]: ...

    def supports_length(self, length: int) -> bool:
        return MIN_PEPTIDE_LENGTH <= length <= MAX_PEPTIDE_LENGTH
