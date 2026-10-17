from __future__ import annotations

import csv
import re
from pathlib import Path
from typing import Dict, Optional

_FORMAT_RE = re.compile(
    r"^(?:HLA[-_]?)?([A-Z]{1,4}[0-9]?)\*?([0-9]{2}):?([0-9]{2,3})(?::[0-9]{2,3}){0,2}([NLSCAQGP]?)$"
)
_TOKEN_RE = re.compile(
    r"(?:HLA-)?[A-Z]{1,4}[0-9]?\*[0-9]{2}:?[0-9]{2,3}(?::[0-9]{2,3}){0,2}[NLSCAQGP]?"
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

    tokens = _TOKEN_RE.findall(cleaned)
    if len(tokens) >= 2:
        parts = [canonicalize(token, aliases) for token in tokens]
        if all(parts):
            return "-".join(sorted(parts))
        return None

    if "*" not in cleaned and "-" in cleaned:
        pieces = cleaned.split("-")
        if not (len(pieces) == 2 and pieces[0] == "HLA"):
            parts = [canonicalize(part, aliases) for part in pieces]
            if len(parts) >= 2 and all(parts):
                return "-".join(sorted(parts))
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
