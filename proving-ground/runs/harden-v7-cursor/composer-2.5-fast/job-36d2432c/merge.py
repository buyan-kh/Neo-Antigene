#!/usr/bin/env python3
import csv
import json
import re
from pathlib import Path
from typing import Dict, List, Optional, Tuple

ROOT = Path(__file__).parent


def load_aliases(path: Path) -> dict[str, str]:
    """Case-insensitive lookup -> canonical_symbol from gene_aliases."""
    lookup: dict[str, str] = {}
    with path.open() as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            canon = row["canonical_symbol"].strip()
            lookup[canon.lower()] = canon
            prev = row["prev_symbol"].strip()
            if prev:
                lookup[prev.lower()] = canon
            syns = row["synonyms"].strip()
            if syns:
                for syn in syns.split("|"):
                    syn = syn.strip()
                    if syn:
                        lookup[syn.lower()] = canon
    return lookup


def split_gene_and_rest(s: str) -> Optional[Tuple[str, str]]:
    """Return (gene_with_hla_prefix, rest_after_gene) or None."""
    s = s.strip()
    # Remove optional HLA- prefix for parsing
    m = re.match(r"(?i)^(HLA-)?(.+)$", s)
    if not m:
        return None
    body = m.group(2)

    # Multi-locus haplotype
    if "-" in body and re.search(r"\*.*-", body):
        return None

    genes = ["DRB1", "DQB1", "DQA1", "DRB3", "DRB4", "DRB5", "DPA1", "DPB1", "DQA2", "DQB2"]
    for g in genes:
        if re.match(rf"(?i)^{g}", body):
            rest = body[len(g) :]
            return (f"HLA-{g.upper()}", rest)

    # Single letter loci A, B, C (possibly with old Cw prefix)
    m2 = re.match(r"(?i)^C(w)?(\*|.*)?", body)
    if m2 and body[0].upper() == "C":
        # Cw*02:02 or C0501 or C*0602
        if body.upper().startswith("CW"):
            rest = body[2:]
            return ("HLA-C", rest)
        if body[0].upper() == "C":
            rest = body[1:]
            return ("HLA-C", rest)

    m3 = re.match(r"(?i)^([AB])(\*|.*)?", body)
    if m3:
        letter = m3.group(1).upper()
        rest = body[1:]
        return (f"HLA-{letter}", rest)

    return None


def insert_star(gene: str, rest: str) -> str:
    rest = rest.lstrip()
    if rest.startswith("*"):
        return gene + rest
    return gene + "*" + rest


def parse_fields_after_star(allele_part: str) -> Tuple[List[str], str]:
    """
    allele_part is like *02:01:01:01G or *1101 or *0702
    Returns (numeric fields as strings with leading zeros normalized, suffix on 2nd field)
    """
    if not allele_part.startswith("*"):
        allele_part = "*" + allele_part
    core = allele_part[1:]

    # Split on colons
    if ":" in core:
        parts = core.split(":")
        fields = []
        suffix = ""
        for i, p in enumerate(parts):
            m = re.match(r"^(\d+)([A-Za-z].*)?$", p)
            if not m:
                fields.append(p)
                continue
            num, suf = m.group(1), m.group(2) or ""
            fields.append(num.zfill(2) if len(num) <= 2 else num)
            if i == 1 and suf:
                suffix = suf
            elif i == len(parts) - 1 and len(parts) > 2 and suf and not suffix:
                # group suffix on last field (e.g. 01G on 4th field)
                suffix = suf
        return fields, suffix

    # Compact digits only: 1101, 0702, 1302, 0201 for DQB1, etc.
    m = re.match(r"^(\d+)([A-Za-z].*)?$", core)
    if not m:
        return [core], ""
    digits, suffix = m.group(1), m.group(2) or ""
    if len(digits) == 4:
        return [digits[:2].zfill(2), digits[2:].zfill(2)], suffix
    if len(digits) == 6:
        return [digits[:2].zfill(2), digits[2:4].zfill(2), digits[4:6].zfill(2)], suffix
    if len(digits) == 2:
        return [digits.zfill(2)], suffix
    # 3 digit like 206 for A0206 - actually A0206 is 4 digits 0206
    return [digits], suffix


def to_two_field(gene: str, fields: list[str], suffix: str) -> str:
    if len(fields) < 2:
        if len(fields) == 1:
            return f"{gene}*{fields[0]}"
        return f"{gene}*"
    f1 = fields[0].zfill(2) if fields[0].isdigit() else fields[0]
    f2 = fields[1].zfill(2) if fields[1].isdigit() else fields[1]
    # suffix from second field or group letter
    if suffix and not f2[-1].isalpha():
        f2 = f2 + suffix
    elif len(fields) >= 2:
        m = re.match(r"^(\d+)([A-Za-z]+)$", fields[1])
        if m:
            f2 = m.group(1).zfill(2) + m.group(2)
    return f"{gene}*{f1}:{f2}"


def normalize_hla(raw: str) -> Optional[str]:
    s = raw.strip()
    if not s:
        return None

    # Fix patterns like HLA-A02:01 -> HLA-A*02:01
    s = re.sub(
        r"(?i)^(HLA-[ABC])\s*(\d{2}:\d{2})",
        r"\1*\2",
        s,
    )
    s = re.sub(
        r"(?i)^(HLA-DRB1)(\d{2}:\d{2})",
        r"\1*\2",
        s,
    )
    s = re.sub(
        r"(?i)^(HLA-B)(\d{2}:\d{2})",
        r"\1*\2",
        s,
    )
    s = re.sub(
        r"(?i)^(HLA-C)(\d{2}:\d{2})",
        r"\1*\2",
        s,
    )
    s = re.sub(
        r"(?i)^(HLA-A)(\d{2}:\d{2})",
        r"\1*\2",
        s,
    )

    # HLA-A29:02 missing star
    s = re.sub(r"(?i)^(HLA-A)(\d{2}:\d{2})", r"\1*\2", s)

    parsed = split_gene_and_rest(s)
    if not parsed:
        return None
    gene, rest = parsed

    # DRB115:01 style already partially fixed; DRB1 without HLA
    if rest and not rest.startswith("*") and re.match(r"^\d", rest):
        rest = "*" + rest

    # B1501, A0206, C0501, DQB10201, DRB10401, A3303, B3801, DRB1*1101
    if gene.startswith("HLA-"):
        rest = rest.lstrip()
        if not rest.startswith("*"):
            if re.match(r"^\d{4,6}$", rest):
                rest = "*" + rest
            elif re.match(r"^\d", rest):
                rest = "*" + rest

    allele_part = rest
    if not allele_part.startswith("*"):
        allele_part = "*" + allele_part

    fields, suffix = parse_fields_after_star(allele_part)
    return to_two_field(gene, fields, suffix)


def canonicalize(raw: str, aliases: dict[str, str]) -> str:
    key = raw.strip().lower()
    if key in aliases:
        return aliases[key]

    norm = normalize_hla(raw)
    if norm:
        return norm

    return raw.strip()


def read_source(path: Path, payload_col: str) -> list[tuple[str, float]]:
    rows = []
    with path.open() as f:
        reader = csv.DictReader(f)
        for row in reader:
            allele = row["allele"]
            payload = float(row[payload_col])
            rows.append((allele, payload))
    return rows


def main():
    aliases = load_aliases(ROOT / "gene_aliases.tsv")
    rows_a = read_source(ROOT / "source_a.csv", "payload_a")
    rows_b = read_source(ROOT / "source_b.csv", "payload_b")

    mapping: dict[str, str] = {}
    by_canon_a: dict[str, tuple[str, float]] = {}
    by_canon_b: dict[str, tuple[str, float]] = {}

    for raw, payload in rows_a:
        canon = canonicalize(raw, aliases)
        mapping[raw] = canon
        if canon in by_canon_a and by_canon_a[canon][1] != payload:
            raise ValueError(f"A conflict for {canon}: {by_canon_a[canon]} vs {raw}, {payload}")
        by_canon_a[canon] = (raw, payload)

    for raw, payload in rows_b:
        canon = canonicalize(raw, aliases)
        mapping[raw] = canon
        if canon in by_canon_b and by_canon_b[canon][1] != payload:
            raise ValueError(f"B conflict for {canon}: {by_canon_b[canon]} vs {raw}, {payload}")
        by_canon_b[canon] = (raw, payload)

    common = sorted(set(by_canon_a) & set(by_canon_b))
    merged = []
    naive_dropped = 0
    for canon in common:
        raw_a, pa = by_canon_a[canon]
        raw_b, pb = by_canon_b[canon]
        if raw_a != raw_b:
            naive_dropped += 1
        merged.append(
            {
                "canonical_allele": canon,
                "a_payload": pa,
                "b_payload": pb,
            }
        )

    out = {
        "merged": merged,
        "joined_count": len(merged),
        "naive_dropped_count": naive_dropped,
        "mapping": mapping,
    }

    out_path = ROOT / "output" / "submission.json"
    with out_path.open("w") as f:
        json.dump(out, f, indent=2, sort_keys=False)

    print(f"joined_count={len(merged)} naive_dropped={naive_dropped}")
    for m in merged:
        ra, rb = by_canon_a[m["canonical_allele"]][0], by_canon_b[m["canonical_allele"]][0]
        flag = "!" if ra != rb else " "
        print(f"{flag} {m['canonical_allele']}: A={ra!r} B={rb!r}")

    unmapped_a = set(by_canon_a) - set(by_canon_b)
    unmapped_b = set(by_canon_b) - set(by_canon_a)
    print("A only:", unmapped_a)
    print("B only:", unmapped_b)


if __name__ == "__main__":
    main()
