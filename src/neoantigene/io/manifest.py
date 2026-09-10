"""Sample manifest.

The manifest describes one patient/sample and is the single entry point for
both supported front ends:

  - `processed`: somatic VCF + expression table + HLA calls (implemented today)
  - `raw`: tumor/normal FASTQ or BAM (consumed by `neoantigene.fastq`, not yet
    implemented) which produces the `processed` fields

Ranking always runs off the `processed` block, so adding the raw front end
later does not change anything downstream.
"""

from __future__ import annotations

from pathlib import Path
from typing import Self

import yaml
from pydantic import BaseModel, ConfigDict, Field, model_validator

from ..hla import Allele


class ManifestModel(BaseModel):
    model_config = ConfigDict(extra="forbid")


class RawInputs(ManifestModel):
    """Unaligned/aligned sequencing inputs. Reserved for the FASTQ front end."""

    tumor_dna: list[Path] = Field(default_factory=list)
    normal_dna: list[Path] = Field(default_factory=list)
    tumor_rna: list[Path] = Field(default_factory=list)
    tumor_bam: Path | None = None
    normal_bam: Path | None = None


class ProcessedInputs(ManifestModel):
    somatic_vcf: Path | None = None
    variant_tsv: Path | None = None
    expression_tsv: Path | None = None
    normal_expression_tsv: Path | None = None
    proteome_fasta: Path | None = None
    tumor_sample_name: str | None = None

    @model_validator(mode="after")
    def _one_variant_source(self) -> Self:
        if bool(self.somatic_vcf) == bool(self.variant_tsv):
            raise ValueError("set exactly one of somatic_vcf or variant_tsv")
        return self

    @property
    def variant_source(self) -> Path:
        source = self.somatic_vcf or self.variant_tsv
        assert source is not None  # guaranteed by _one_variant_source
        return source


class Sample(ManifestModel):
    sample_id: str = Field(min_length=1)
    cancer_type: str | None = None
    tumor_purity: float | None = Field(default=None, gt=0.0, le=1.0)
    hla: list[Allele] = Field(default_factory=list)
    processed: ProcessedInputs
    raw: RawInputs = Field(default_factory=RawInputs)

    @classmethod
    def load(cls, path: Path) -> Sample:
        """Read a manifest, resolving relative paths against its own location."""
        with open(path) as handle:
            payload = yaml.safe_load(handle) or {}
        sample = cls.model_validate(payload)
        return sample._with_resolved_paths(Path(path).resolve().parent)

    def _with_resolved_paths(self, base: Path) -> Sample:
        return self.model_copy(
            update={
                "processed": _resolve_block(self.processed, base),
                "raw": _resolve_block(self.raw, base),
            }
        )

    def require_proteome(self) -> Path:
        if self.processed.proteome_fasta is None:
            raise ValueError(f"{self.sample_id}: processed.proteome_fasta is required")
        return self.processed.proteome_fasta

    def require_hla(self) -> list[str]:
        if not self.hla:
            raise ValueError(f"{self.sample_id}: no HLA alleles in manifest")
        return self.hla

    def input_paths(self) -> list[Path]:
        """Every existing input file, for run-manifest digests."""
        candidates: list[Path] = []
        for block in (self.processed, self.raw):
            for _, value in block:
                if isinstance(value, Path):
                    candidates.append(value)
                elif isinstance(value, list):
                    candidates.extend(p for p in value if isinstance(p, Path))
        return [p for p in candidates if p.exists()]


def _resolve_block[T: ManifestModel](block: T, base: Path) -> T:
    updates: dict[str, Path | list[Path]] = {}
    for name, value in block:
        if isinstance(value, Path):
            updates[name] = _resolve(value, base)
        elif isinstance(value, list) and value and all(isinstance(p, Path) for p in value):
            updates[name] = [_resolve(p, base) for p in value]
    return block.model_copy(update=updates)


def _resolve(path: Path, base: Path) -> Path:
    return path if path.is_absolute() else (base / path).resolve()
