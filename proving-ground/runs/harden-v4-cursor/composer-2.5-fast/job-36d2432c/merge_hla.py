#!/usr/bin/env python3
import csv
import json
import re
from pathlib import Path

ROOT = Path(__file__).resolve().parent

GENES = ("DRB1", "DQB1", "DQA1", "DPB1", "DPA1", "DRB3", "DRB4", "DRB5")


def load_aliases(path: Path) -> dict[str, str]:
    alias_to_canonical: dict[str, str] = {}
    with path.open() as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            canonical = row["canonical_symbol"].strip()
            for key in ("canonical_symbol", "prev_symbol", "synonyms"):
                if key == "synonyms":
                    parts = row["synonyms"].split("|") if row["synonyms"] else []
                else:
                    parts = [row[key]]
                for p in parts:
                    p = p.strip()
                    if p:
                        alias_to_canonical[p] = canonical
                        alias_to_canonical[p.lower()] = canonical
    return alias_to_canonical


def split_allele_digits(digits: str) -> tuple[str, str, str]:
    """Return field1, field2, remainder suffix (e.g. empty or 'G' or 'N')."""
    digits = digits.replace(":", "")
    m = re.match(r"^(\d{2})(\d{2})(.*)$", digits)
    if not m:
        raise ValueError(f"cannot split digits: {digits}")
    f1, f2, rest = m.group(1), m.group(2), m.group(3)
    return f1, f2, rest


def format_canonical(gene: str, f1: str, f2: str, suffix: str = "") -> str:
    gene = gene.upper()
    if not gene.startswith("HLA-"):
        gene = f"HLA-{gene}"
    return f"{gene}*{f1}:{f2}{suffix}"


def parse_to_canonical(raw: str, alias_to_canonical: dict[str, str]) -> str:
    s = raw.strip()
    if not s:
        raise ValueError("empty allele")

    if s in alias_to_canonical:
        return alias_to_canonical[s]

    lower = s.lower()
    if lower in alias_to_canonical:
        return alias_to_canonical[lower]

    # Normalize case for parsing; preserve suffix letters on two-field
    work = s.strip()

    # HLA-GENE*fields or hla-gene*fields
    m = re.match(
        r"^(?:HLA-)?([A-Za-z0-9]+)\*(\d+(?::\d+)*)([A-Za-z]?)$",
        work,
        re.IGNORECASE,
    )
    if m:
        gene, fields, trailing_suffix = m.group(1), m.group(2), m.group(3) or ""
        parts = fields.split(":")
        f1, f2 = parts[0], parts[1] if len(parts) > 1 else None
        if f2 is None:
            if len(f1) >= 4:
                f1, f2 = f1[:2], f1[2:4]
                extra = f1[4:] if len(f1) > 4 else ""
            else:
                raise ValueError(f"short field: {raw}")
        else:
            extra = trailing_suffix
            if len(parts) > 2:
                pass  # drop higher resolution fields
        # suffix on last part e.g. 01:04N in HLA-A*01:04N
        m2 = re.match(r"^(\d{2})([A-Za-z]?)$", f2)
        if m2:
            f2, sfx = m2.group(1), m2.group(2)
            extra = extra or sfx
        return format_canonical(gene, f1.zfill(2) if len(f1) == 1 else f1[:2], f2, extra)

    # HLA-A*01:04N style with N on second field
    m = re.match(
        r"^(?:HLA-)?([A-Za-z0-9]+)\*(\d{2}):(\d{2})([A-Za-z]?)$",
        work,
        re.IGNORECASE,
    )
    if m:
        gene, f1, f2, sfx = m.groups()
        return format_canonical(gene, f1, f2, sfx or "")

    # HLA-DRB115:01 or HLA-DQB106:02 — gene digits merged into locus name
    m = re.match(r"^HLA-(DRB1|DQB1|DQA1|DPB1)(\d{2}):(\d{2})$", work, re.IGNORECASE)
    if m:
        gene, f1, f2 = m.group(1), m.group(2), m.group(3)
        return format_canonical(gene, f1, f2)

    # HLA-A26:01 or HLA-C03:04 or HLA-B14:02 (no asterisk)
    m = re.match(r"^HLA-([ABC])(\d{2}):(\d{2})$", work, re.IGNORECASE)
    if m:
        return format_canonical(m.group(1), m.group(2), m.group(3))

    # HLA-B07:02 style
    m = re.match(r"^HLA-([ABC])(\d{2}):(\d{2})$", work, re.IGNORECASE)
    if m:
        return format_canonical(m.group(1), m.group(2), m.group(3))

    # HLA-A02:01
    m = re.match(r"^HLA-([ABC])(\d{2}):(\d{2})$", work, re.IGNORECASE)
    if m:
        return format_canonical(m.group(1), m.group(2), m.group(3))

    # Short: B*15:01, C*12:03, A*33:03
    m = re.match(r"^([ABC])\*(\d{2}):(\d{2})$", work, re.IGNORECASE)
    if m:
        return format_canonical(m.group(1), m.group(2), m.group(3))

    # DRB1*11:01 or DRB1*1101
    m = re.match(r"^(DRB1|DQB1|DQA1)\*(\d{2}):(\d{2})$", work, re.IGNORECASE)
    if m:
        return format_canonical(m.group(1), m.group(2), m.group(3))

    m = re.match(r"^(DRB1|DQB1|DQA1)(\d{4,})$", work, re.IGNORECASE)
    if m:
        gene, digits = m.group(1), m.group(2)
        f1, f2, _ = split_allele_digits(digits)
        return format_canonical(gene, f1, f2)

    m = re.match(r"^(DRB1|DQB1|DQA1)\*(\d{4,})$", work, re.IGNORECASE)
    if m:
        gene, digits = m.group(1), m.group(2)
        f1, f2, _ = split_allele_digits(digits)
        return format_canonical(gene, f1, f2)

    # C*0202, C*0702, B*1302
    m = re.match(r"^([ABC])\*(\d{4})$", work, re.IGNORECASE)
    if m:
        f1, f2, _ = split_allele_digits(m.group(2))
        return format_canonical(m.group(1), f1, f2)

    # A*3303, A*0104, A*2402
    m = re.match(r"^([ABC])\*(\d{4,})$", work, re.IGNORECASE)
    if m:
        f1, f2, _ = split_allele_digits(m.group(2))
        return format_canonical(m.group(1), f1, f2)

    # Compact: C0501, B0801, A0206, A2902, DRB10301, B3801, DQB1*0501 handled above
    m = re.match(r"^([ABC])(\d{4})$", work, re.IGNORECASE)
    if m:
        f1, f2, _ = split_allele_digits(m.group(2))
        return format_canonical(m.group(1), f1, f2)

    m = re.match(r"^(DRB1)(\d{4,})$", work, re.IGNORECASE)
    if m:
        f1, f2, _ = split_allele_digits(m.group(2))
        return format_canonical(m.group(1), f1, f2)

    raise ValueError(f"unparsed allele: {raw!r}")


def read_source(path: Path, payload_col: str) -> list[tuple[str, float]]:
    rows = []
    with path.open() as f:
        reader = csv.DictReader(f)
        for row in reader:
            rows.append((row["allele"], float(row[payload_col])))
    return rows


def main():
    alias_to_canonical = load_aliases(ROOT / "gene_aliases.tsv")
    mapping: dict[str, str] = {}

    def canonicalize(raw: str) -> str:
        if raw not in mapping:
            mapping[raw] = parse_to_canonical(raw, alias_to_canonical)
        return mapping[raw]

    a_rows = read_source(ROOT / "source_a.csv", "payload_a")
    b_rows = read_source(ROOT / "source_b.csv", "payload_b")

    a_by_canonical: dict[str, tuple[str, float]] = {}
    for raw, payload in a_rows:
        c = canonicalize(raw)
        if c in a_by_canonical:
            prev_raw, prev_p = a_by_canonical[c]
            if prev_p != payload:
                raise ValueError(f"A conflict {c}: {prev_raw} vs {raw}")
        else:
            a_by_canonical[c] = (raw, payload)

    b_by_canonical: dict[str, tuple[str, float]] = {}
    for raw, payload in b_rows:
        c = canonicalize(raw)
        if c in b_by_canonical:
            prev_raw, prev_p = b_by_canonical[c]
            if prev_p != payload:
                raise ValueError(f"B conflict {c}: {prev_raw} vs {raw}")
        else:
            b_by_canonical[c] = (raw, payload)

    merged = []
    naive_dropped = 0
    for c in sorted(a_by_canonical.keys() & b_by_canonical.keys()):
        raw_a, pa = a_by_canonical[c]
        raw_b, pb = b_by_canonical[c]
        if raw_a != raw_b:
            naive_dropped += 1
        merged.append(
            {"canonical_allele": c, "a_payload": pa, "b_payload": pb}
        )

    out = {
        "merged": merged,
        "joined_count": len(merged),
        "naive_dropped_count": naive_dropped,
        "mapping": mapping,
    }

    out_path = ROOT / "output" / "submission.json"
    out_path.write_text(json.dumps(out, indent=2) + "\n")
    print(f"joined_count={len(merged)} naive_dropped={naive_dropped}")
    print(f"mapping entries={len(mapping)}")


if __name__ == "__main__":
    main()
