from __future__ import annotations


def strip_version(identifier: str | None) -> str:
    """Drop an Ensembl-style version suffix: ENST00000256078.10 -> ENST00000256078."""
    if not identifier:
        return ""
    return identifier.strip().split(".", 1)[0]
