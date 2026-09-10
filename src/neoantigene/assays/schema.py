"""Ground-truth assay records.

This schema is the product. Everything upstream is a hypothesis generator;
these rows are the only thing that says whether a ranking was right.

One row per (sample, peptide, allele, assay) with an explicit call. Ambiguous
results stay ambiguous — coercing `indeterminate` to `negative` silently
inflates every downstream metric.
"""

from __future__ import annotations

import csv
from collections.abc import Iterator, Sequence
from datetime import date
from enum import StrEnum
from pathlib import Path
from typing import Any

from pydantic import BaseModel, ConfigDict, Field

from ..hla import Allele

REQUEST_COLUMNS: tuple[str, ...] = (
    "sample_id",
    "peptide",
    "allele",
    "assay",
    "call",
    "effect_size",
    "replicate_count",
    "run_date",
    "operator",
    "notes",
    "predicted_rank",
    "predicted_score",
)


class AssayType(StrEnum):
    IFNG_ELISPOT = "ifng_elispot"
    ICS = "intracellular_cytokine_staining"
    MULTIMER = "pmhc_multimer"
    ACTIVATION_INDUCED_MARKER = "aim"


class AssayCall(StrEnum):
    POSITIVE = "positive"
    NEGATIVE = "negative"
    INDETERMINATE = "indeterminate"


class AssayResult(BaseModel):
    model_config = ConfigDict(extra="forbid", frozen=True)

    sample_id: str = Field(min_length=1)
    peptide: str = Field(min_length=1)
    allele: Allele | None = None
    assay: AssayType
    call: AssayCall
    effect_size: float | None = Field(
        default=None,
        description="Assay-native magnitude, e.g. background-subtracted SFU/1e6 cells.",
    )
    replicate_count: int | None = Field(default=None, ge=1)
    run_date: date | None = None
    operator: str | None = None
    notes: str | None = None
    predicted_rank: int | None = Field(default=None, ge=1)
    predicted_score: float | None = None

    @property
    def key(self) -> str:
        return f"{self.sample_id}|{self.peptide}|{self.allele or '*'}"

    @property
    def label(self) -> int | None:
        """1 positive, 0 negative, None if the assay could not call it."""
        match self.call:
            case AssayCall.POSITIVE:
                return 1
            case AssayCall.NEGATIVE:
                return 0
            case AssayCall.INDETERMINATE:
                return None


def read_results(path: Path) -> list[AssayResult]:
    return list(iter_results(path))


def iter_results(path: Path) -> Iterator[AssayResult]:
    with open(path) as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            payload: dict[str, Any] = {
                key: (value if value not in ("", ".") else None)
                for key, value in row.items()
                if key
            }
            yield AssayResult.model_validate(payload)


def write_request(
    ranked_records: Sequence[dict[str, Any]],
    path: Path,
    assay: AssayType = AssayType.IFNG_ELISPOT,
) -> Path:
    """Emit a synthesis/assay request sheet from a ranked shortlist.

    Pre-filled with rank and score so the returned sheet joins back without a
    separate key file — the lab edits `call` and `effect_size`.
    """
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=REQUEST_COLUMNS, delimiter="\t")
        writer.writeheader()
        for record in ranked_records:
            writer.writerow(
                {
                    "sample_id": record["sample_id"],
                    "peptide": record["mutant_peptide"],
                    "allele": record["allele"],
                    "assay": assay.value,
                    "call": "",
                    "effect_size": "",
                    "replicate_count": "",
                    "run_date": "",
                    "operator": "",
                    "notes": "",
                    "predicted_rank": record["rank"],
                    "predicted_score": record["score"],
                }
            )
    return path
