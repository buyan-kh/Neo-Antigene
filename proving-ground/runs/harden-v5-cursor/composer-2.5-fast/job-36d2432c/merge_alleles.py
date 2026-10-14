#!/usr/bin/env python3
import csv
import json
import re
from pathlib import Path

ROOT = Path(__file__).resolve().parent


def load_alias_map() -> dict[str, str]:
    alias_to_canonical: dict[str, str] = {}
    with open(ROOT / "gene_aliases.tsv", newline="") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            canonical = row["canonical_symbol"].strip()
            for key in [row["canonical_symbol"], row["prev_symbol"]]:
                if key:
                    alias_to_canonical[key.strip()] = canonical
                    alias_to_canonical[key.strip().lower()] = canonical
            syns = row.get("synonyms") or ""
            for syn in syns.split("|"):
                syn = syn.strip()
                if syn:
                    alias_to_canonical[syn] = canonical
                    alias_to_canonical[syn.lower()] = canonical
    return alias_to_canonical


def two_field_canonical(gene: str, f1: str, f2: str, suffix: str = "") -> str:
    if not gene.upper().startswith("HLA-"):
        gene = "HLA-" + gene.upper()
    else:
        gene = "HLA-" + gene[4:].upper() if gene.upper().startswith("HLA-") else gene
    # Normalize gene part after HLA-
    g = gene[4:]
    if g.upper() == "CW":
        g = "C"
    gene = "HLA-" + g.upper()
    return f"{gene}*{f1.zfill(2)}:{f2.zfill(2)}{suffix}"


def parse_allele(raw: str, alias_map: dict[str, str]) -> str:
    original = raw
    s = raw.strip()

    if s in alias_map:
        return alias_map[s]
    if s.lower() in alias_map:
        return alias_map[s.lower()]

    # HLA-A*02:01G style — capture suffix on second field
    m = re.match(
        r"^(?i)(hla-)?([a-z0-9]+)[\*:](\d{1,2}):(\d{1,2})([a-z]?(?:n|g|l|q|s)?)(?::\d+)*$",
        s.replace(" ", ""),
    )
    if m:
        prefix, gene, f1, f2, suf = m.groups()
        suf = suf or ""
        if gene.upper() == "CW":
            gene = "C"
        if not prefix:
            gene = gene.upper()
        return two_field_canonical(
            f"HLA-{gene.upper()}" if prefix else gene, f1, f2, suf.upper() if len(suf) == 1 and suf.isalpha() else suf
        )

    # HLA-Cw*02:02
    m = re.match(r"^(?i)hla-cw\*(\d{1,2}):(\d{1,2})", s.replace(" ", ""))
    if m:
        return two_field_canonical("C", m.group(1), m.group(2))

    # HLA-A02:01 or HLA-B14:02 or HLA-C12:03 (missing * between gene and digits)
    m = re.match(r"^(?i)hla-([a-z0-9]+)(\d{2}):(\d{2})", s.replace(" ", ""))
    if m:
        gene, f1, f2 = m.groups()
        if gene.upper().startswith("DRB1"):
            return two_field_canonical("DRB1", f1, f2)
        if gene.upper().startswith("DQB1"):
            return two_field_canonical("DQB1", f1, f2)
        return two_field_canonical(gene, f1, f2)

    # HLA-DRB115:01 or HLA-DRB104:01
    m = re.match(r"^(?i)hla-drb1(\d{2}):(\d{2})", s.replace(" ", ""))
    if m:
        return two_field_canonical("DRB1", m.group(1), m.group(2))

    # A*29:02, B*18:01, DRB1*15:01, DQB1*03:01
    m = re.match(
        r"^(?i)([a-z]\*?|drb1\*?|dqb1\*?)(\d{1,2}):(\d{1,2})([a-z]?)(?::\d+)*$",
        s.replace(" ", ""),
    )
    if m:
        gene, f1, f2, suf = m.groups()
        gene = gene.replace("*", "").upper()
        if gene == "A":
            gene = "A"
        return two_field_canonical(gene, f1, f2, (suf or "").upper())

    # A*0104, C*0602, B*1302, DRB1*1101 (compressed 4 digits after *)
    m = re.match(r"^(?i)([a-z]\*|drb1\*|dqb1\*|c\*|b\*)(\d{4})([a-z]?)$", s.replace(" ", ""))
    if m:
        gene_part, digits, suf = m.groups()
        gene = gene_part.replace("*", "").upper()
        f1, f2 = digits[:2], digits[2:]
        return two_field_canonical(gene, f1, f2, (suf or "").upper())

    # B*0702 style without colon
    m = re.match(r"^(?i)([a-z]\*|drb1\*|dqb1\*|c\*|b\*)(\d{4})$", s.replace(" ", ""))
    if m:
        gene_part, digits = m.groups()
        gene = gene_part.replace("*", "").upper()
        f1, f2 = digits[:2], digits[2:]
        return two_field_canonical(gene, f1, f2)

    # Old 4-letter digit: B1501, C0501, A0206, B3801, A3303, DQB10201
    m = re.match(r"^(?i)(dqb1|drb1|[a-z])(\d{4})$", s.replace(" ", ""))
    if m:
        gene, digits = m.groups()
        gene = gene.upper()
        f1, f2 = digits[:2], digits[2:]
        return two_field_canonical(gene, f1, f2)

    # HLA-A29:02 missing *
    m = re.match(r"^(?i)hla-([a-z])(\d{2}):(\d{2})", s.replace(" ", ""))
    if m:
        return two_field_canonical(m.group(1), m.group(2), m.group(3))

    raise ValueError(f"Cannot parse allele: {original!r}")


def load_source(path: Path) -> list[tuple[str, float]]:
    rows = []
    with open(path, newline="") as f:
        reader = csv.DictReader(f)
        for row in reader:
            allele = row["allele"]
            payload_key = [k for k in row if k != "allele"][0]
            rows.append((allele, float(row[payload_key])))
    return rows


def main():
    alias_map = load_alias_map()
    source_a = load_source(ROOT / "source_a.csv")
    source_b = load_source(ROOT / "source_b.csv")

    mapping: dict[str, str] = {}
    for allele, _ in source_a + source_b:
        if allele not in mapping:
            mapping[allele] = parse_allele(allele, alias_map)

    a_by_canonical: dict[str, tuple[str, float]] = {}
    for allele, payload in source_a:
        canon = mapping[allele]
        if canon in a_by_canonical:
            prev_raw, prev_p = a_by_canonical[canon]
            if prev_p != payload or prev_raw != allele:
                raise ValueError(f"Conflict in A for {canon}: {prev_raw}/{prev_p} vs {allele}/{payload}")
        a_by_canonical[canon] = (allele, payload)

    b_by_canonical: dict[str, tuple[str, float]] = {}
    for allele, payload in source_b:
        canon = mapping[allele]
        if canon in b_by_canonical:
            prev_raw, prev_p = b_by_canonical[canon]
            if prev_p != payload or prev_raw != allele:
                raise ValueError(f"Conflict in B for {canon}: {prev_raw}/{prev_p} vs {allele}/{payload}")
        b_by_canonical[canon] = (allele, payload)

    common = sorted(set(a_by_canonical) & set(b_by_canonical))
    merged = []
    naive_dropped = 0
    for canon in common:
        raw_a, pa = a_by_canonical[canon]
        raw_b, pb = b_by_canonical[canon]
        if raw_a != raw_b:
            naive_dropped += 1
        merged.append(
            {"canonical_allele": canon, "a_payload": pa, "b_payload": pb}
        )

    out = {
        "merged": merged,
        "joined_count": len(merged),
        "naive_dropped_count": naive_dropped,
        "mapping": dict(sorted(mapping.items(), key=lambda x: x[0])),
    }

    out_path = ROOT / "output" / "submission.json"
    with open(out_path, "w") as f:
        json.dump(out, f, indent=2)
        f.write("\n")

    print(f"joined_count={len(merged)} naive_dropped={naive_dropped}")
    only_a = set(a_by_canonical) - set(b_by_canonical)
    only_b = set(b_by_canonical) - set(a_by_canonical)
    print("only_a:", sorted(only_a))
    print("only_b:", sorted(only_b))


if __name__ == "__main__":
    main()
