"""Oracle-backed checks on published claims.

Every check here answers with arithmetic rather than with a model, which is
what makes it safe to point a large number of agents at the literature: a
claim is accepted because a deterministic tool confirmed it, not because an
agent found it plausible. See `audit.reconcile` for the reasoning.
"""

from .reconcile import (
    EpitopeClaim,
    ProteinChange,
    Reconciliation,
    ReconciliationReport,
    SymbolIndex,
    Verdict,
    parse_protein_change,
    read_claims,
    reconcile,
    reconcile_claim,
    write_report,
)

__all__ = [
    "EpitopeClaim",
    "ProteinChange",
    "Reconciliation",
    "ReconciliationReport",
    "SymbolIndex",
    "Verdict",
    "parse_protein_change",
    "read_claims",
    "reconcile",
    "reconcile_claim",
    "write_report",
]
