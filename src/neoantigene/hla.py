"""HLA allele parsing and normalization.

Lives at package top level rather than under `io` because allele identity is a
domain concept that the models validate against, and `io` imports the models.
"""

from __future__ import annotations

import csv
import re
from pathlib import Path
from typing import Annotated

from pydantic import AfterValidator

_ALLELE_RE = re.compile(
    r"^(?:HLA-)?([ABC]|D[PQR][AB]1?)\*?(\d{2,3}):?(\d{2,3})(?::\d{2,3})*[A-Z]?$",
    re.IGNORECASE,
)

SUPPORTED_CLASS_I = ("A", "B", "C")


class HLAError(ValueError):
    pass


def normalize_allele(raw: str) -> str:
    """Return a canonical `HLA-A*02:01` string.

    Accepts `A0201`, `A*02:01`, `HLA-A02:01`, and 6/8-digit forms, which are
    truncated to two fields since presentation predictors are trained at that
    resolution.
    """
    token = raw.strip()
    if not token:
        raise HLAError("empty allele")
    match = _ALLELE_RE.match(token)
    if not match:
        raise HLAError(f"unrecognized HLA allele: {raw!r}")
    gene, field1, field2 = match.groups()
    return f"HLA-{gene.upper()}*{field1}:{field2}"


#: Any string field that must hold a canonical HLA allele. Normalizes on the
#: way in, so downstream code can compare alleles with `==`.
Allele = Annotated[str, AfterValidator(normalize_allele)]


def gene_of(allele: str) -> str:
    return allele.split("-", 1)[1].split("*", 1)[0]


def class_i_only(alleles: list[str]) -> list[str]:
    return [a for a in alleles if gene_of(a) in SUPPORTED_CLASS_I]


def parse_hla_string(value: str) -> list[str]:
    """Parse a comma/space separated genotype string."""
    tokens = [t for t in re.split(r"[,\s;]+", value) if t]
    return dedupe([normalize_allele(t) for t in tokens])


def read_hla_file(path: Path) -> list[str]:
    """Read a genotype from a one-allele-per-line list or a TSV with an `allele` column."""
    with open(path) as handle:
        first = handle.readline()
        handle.seek(0)
        if "\t" in first and "allele" in first.lower():
            reader = csv.DictReader(handle, delimiter="\t")
            column = next(c for c in reader.fieldnames or [] if c.lower() == "allele")
            return dedupe([normalize_allele(row[column]) for row in reader if row[column]])
        alleles: list[str] = []
        for line in handle:
            stripped = line.strip()
            if not stripped or stripped.startswith("#"):
                continue
            alleles.extend(parse_hla_string(stripped))
        return dedupe(alleles)


def dedupe(alleles: list[str]) -> list[str]:
    seen: set[str] = set()
    out: list[str] = []
    for allele in alleles:
        if allele not in seen:
            seen.add(allele)
            out.append(allele)
    return out
