from __future__ import annotations

import re
from typing import Optional

_ALLELE_RE = re.compile(
    r"^(?:HLA[-_]?)?([A-Z]{1,4}[0-9]?)\*?([0-9]{2}):?([0-9]{2,3})(?::[0-9]{2,3}){0,2}$"
)


def normalize_peptide(raw: object) -> Optional[str]:
    if not isinstance(raw, str):
        return None
    cleaned = "".join(raw.split()).upper()
    if not cleaned or not cleaned.isalpha():
        return None
    return cleaned


def normalize_allele(raw: object) -> Optional[str]:
    if not isinstance(raw, str):
        return None
    cleaned = "".join(raw.split()).upper().replace("_", "-")
    match = _ALLELE_RE.match(cleaned)
    if not match:
        return None
    gene, field_one, field_two = match.groups()
    return "HLA-{}*{}:{}".format(gene, field_one, field_two)


def normalize_pair(peptide: object, allele: object) -> Optional[tuple]:
    pep = normalize_peptide(peptide)
    allele_norm = normalize_allele(allele)
    if pep is None or allele_norm is None:
        return None
    return (pep, allele_norm)
