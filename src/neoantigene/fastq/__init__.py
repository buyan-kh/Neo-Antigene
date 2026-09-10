"""Raw-sequencing front end (not yet implemented).

Neo Antigene currently starts from a VEP-annotated somatic VCF, an expression
table and an HLA genotype. This package is the reserved home for producing
those three artifacts from tumor/normal FASTQ or BAM:

  align -> somatic calling (tumor/normal) -> RNA quantification -> HLA typing

It is a separate package on purpose. The ranking engine is the defensible
part of the product; read alignment and somatic calling are commodity steps
that most customer labs already run, and several will insist on their own.
Keeping the boundary at `Sample.processed` means adopting this front end never
changes ranking behaviour.
"""

from __future__ import annotations

from pathlib import Path

from ..io.manifest import ProcessedInputs, Sample


def run_front_end(sample: Sample, work_dir: Path) -> ProcessedInputs:
    """Produce `ProcessedInputs` from `sample.raw`."""
    raise NotImplementedError(
        "the FASTQ/BAM front end is not implemented; populate manifest "
        "`processed:` fields with an externally produced VCF, expression "
        "table and HLA genotype"
    )
