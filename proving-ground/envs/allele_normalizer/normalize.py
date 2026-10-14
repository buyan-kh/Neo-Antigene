from __future__ import annotations

import csv
import re
from pathlib import Path
from typing import Dict, Optional

_FORMAT_RE = re.compile(
    r"^(?:HLA[-_]?)?([A-Z]{1,4}[0-9]?)\*?([0-9]{2}):?([0-9]{2,3})(?::[0-9]{2,3}){0,2}([NLSCAQGP]?)$"
)


def load_aliases(path: Path) -> Dict[str, str]:
    mapping: Dict[str, str] = {}
    with open(path, "r", encoding="utf-8") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            canonical = row["canonical_symbol"].strip()
            previous = row.get("prev_symbol", "").strip()
            if previous:
                mapping[previous.upper()] = canonical
            for synonym in (row.get("synonyms") or "").split("|"):
                synonym = synonym.strip()
                if synonym:
                    mapping[synonym.upper()] = canonical
    return mapping


def canonicalize(raw: object, aliases: Dict[str, str]) -> Optional[str]:
    if not isinstance(raw, str):
        return None
    cleaned = "".join(raw.split()).upper().replace("_", "-")
    if not cleaned:
        return None

    if cleaned.count("*") >= 2 and "-" in cleaned:
        parts = [canonicalize(part, aliases) for part in cleaned.split("-")]
        if all(parts):
            return "-".join(parts)
        return None

    match = _FORMAT_RE.match(cleaned)
    if match:
        gene, field_one, field_two, suffix = match.groups()
        if gene == "CW":
            gene = "C"
        return "HLA-{}*{}:{}{}".format(gene, field_one, field_two, suffix)

    if cleaned in aliases:
        return aliases[cleaned]

    return None
