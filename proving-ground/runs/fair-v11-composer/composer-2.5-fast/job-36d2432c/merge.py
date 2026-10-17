#!/usr/bin/env python3
from __future__ import annotations

import csv
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[4]
sys.path.insert(0, str(ROOT))

from envs.allele_normalizer.normalize import canonicalize, load_aliases  # noqa: E402

WORKDIR = Path(__file__).resolve().parent


def read_source(path: Path, payload_col: str) -> tuple[dict[str, float], list[tuple[str, float]]]:
    rows: list[tuple[str, float]] = []
    with open(path, newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle):
            raw = row["allele"]
            rows.append((raw, float(row[payload_col])))
    by_canonical: dict[str, float] = {}
    for raw, value in rows:
        by_canonical[raw] = value  # placeholder, filled after canonicalize
    return {}, rows


def main() -> None:
    aliases = load_aliases(WORKDIR / "gene_aliases.tsv")

    payload_a: dict[str, float] = {}
    raw_a_rows: list[tuple[str, float]] = []
    with open(WORKDIR / "source_a.csv", newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle):
            raw = row["allele"]
            val = float(row["payload_a"])
            raw_a_rows.append((raw, val))
            canon = canonicalize(raw, aliases)
            if canon:
                payload_a[canon] = val

    payload_b: dict[str, float] = {}
    raw_b_rows: list[tuple[str, float]] = []
    with open(WORKDIR / "source_b.csv", newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle):
            raw = row["allele"]
            val = float(row["payload_b"]) / 100.0
            raw_b_rows.append((raw, val))
            canon = canonicalize(raw, aliases)
            if canon:
                payload_b[canon] = val

    mapping: dict[str, str] = {}
    for raw, _ in raw_a_rows + raw_b_rows:
        if raw not in mapping:
            canon = canonicalize(raw, aliases)
            if canon:
                mapping[raw] = canon

    joined = sorted(set(payload_a) & set(payload_b))
    raw_a_strings = {raw for raw, _ in raw_a_rows}
    raw_b_strings = {raw for raw, _ in raw_b_rows}
    naive_joined = raw_a_strings & raw_b_strings
    naive_dropped = len(joined) - len(naive_joined)

    merged = [
        {
            "canonical_allele": c,
            "a_payload": payload_a[c],
            "b_payload": payload_b[c],
        }
        for c in joined
    ]

    out = {
        "merged": merged,
        "joined_count": len(joined),
        "naive_dropped_count": naive_dropped,
        "mapping": dict(sorted(mapping.items())),
    }

    out_path = WORKDIR / "output" / "submission.json"
    out_path.write_text(json.dumps(out, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(f"Wrote {out_path}: joined={len(joined)} naive_dropped={naive_dropped}")


if __name__ == "__main__":
    main()
