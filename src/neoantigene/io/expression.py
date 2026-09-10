"""Tumor RNA expression and normal-tissue reference expression."""

from __future__ import annotations

import csv
from pathlib import Path
from typing import Self

from .ids import strip_version


class ExpressionTable:
    """Transcript- and gene-level TPM lookup.

    Input is a TSV with a `tpm` column plus at least one of
    `transcript_id` / `gene` (or `gene_name`, `gene_id`).
    """

    def __init__(
        self,
        by_transcript: dict[str, float] | None = None,
        by_gene: dict[str, float] | None = None,
    ) -> None:
        self.by_transcript = by_transcript or {}
        self.by_gene = by_gene or {}

    @classmethod
    def from_tsv(cls, path: Path) -> Self:
        by_transcript: dict[str, float] = {}
        by_gene: dict[str, float] = {}
        with open(path) as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            fields = {name.lower(): name for name in (reader.fieldnames or [])}
            tpm_column = fields.get("tpm")
            if tpm_column is None:
                raise ValueError(f"{path}: expression table needs a 'tpm' column")
            transcript_column = fields.get("transcript_id") or fields.get("transcript")
            gene_column = fields.get("gene") or fields.get("gene_name") or fields.get("gene_id")
            if transcript_column is None and gene_column is None:
                raise ValueError(f"{path}: expression table needs a transcript or gene column")

            for row in reader:
                raw_tpm = (row.get(tpm_column) or "").strip()
                if not raw_tpm:
                    continue
                tpm = float(raw_tpm)
                if tpm < 0:
                    raise ValueError(f"{path}: negative TPM {tpm} in row {row}")
                if transcript_column:
                    key = strip_version(row[transcript_column])
                    if key:
                        by_transcript[key] = tpm
                if gene_column:
                    key = strip_version(row[gene_column])
                    if key:
                        by_gene[key] = max(by_gene.get(key, 0.0), tpm)
        return cls(by_transcript, by_gene)

    def lookup(self, transcript: str | None, gene: str | None) -> float | None:
        if transcript:
            value = self.by_transcript.get(strip_version(transcript))
            if value is not None:
                return value
        if gene:
            return self.by_gene.get(strip_version(gene))
        return None

    def __len__(self) -> int:
        return len(self.by_transcript) + len(self.by_gene)


class NormalExpressionReference(ExpressionTable):
    """Median expression across healthy tissues (e.g. a GTEx summary table).

    Used as an on-target/off-tumor risk term, not as a hard filter: a gene that
    is highly expressed in healthy tissue is both a tolerance risk and a
    toxicity risk, but the mutant peptide may still be tumor-restricted.
    """
