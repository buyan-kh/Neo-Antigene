"""Structure-based secondary filtering.

Deliberately a re-ranker over an already-shortlisted set, never a primary
score. Structure prediction is slow, and pMHC models mostly recapitulate what
the sequence-based presentation predictor already said; its value is catching
implausible geometry in the handful of peptides heading to synthesis.

Not implemented — see `structure.filter`.
"""

from .filter import StructureFilter, StructureVerdict

__all__ = ["StructureFilter", "StructureVerdict"]
