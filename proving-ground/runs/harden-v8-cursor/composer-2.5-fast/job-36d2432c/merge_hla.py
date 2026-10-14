#!/usr/bin/env python3
import csv
import json
import re
from collections import defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parent


def load_aliases(path):
    alias_to_canonical = {}
    with open(path, newline="") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            canon = row["canonical_symbol"].strip()
            for key in [row["prev_symbol"]] + row["synonyms"].split("|"):
                key = key.strip()
                if key:
                    alias_to_canonical[key.upper()] = canon
                    alias_to_canonical[key] = canon
    return alias_to_canonical


def two_field_canonical(gene: str, field1: str, field2: str, suffix: str = "") -> str:
    return f"HLA-{gene}*{field1}:{field2}{suffix}"


def parse_compact_allele(s: str):
    """Parse forms like A3201, B1801, DQB10301, DRB11101, A6801, A0206."""
    s = s.strip().upper()
    m = re.match(r"^(DQB1|DRB1|DRB|DQA1|A|B|C)(\d{3,4})$", s)
    if m:
        prefix, digits = m.group(1), m.group(2)
        gene_map = {
            "A": "A",
            "B": "B",
            "C": "C",
            "DQB1": "DQB1",
            "DQA1": "DQA1",
            "DRB1": "DRB1",
            "DRB": "DRB1",
        }
        gene = gene_map[prefix]
        if len(digits) == 4:
            f1, f2 = digits[:2], digits[2:]
        elif len(digits) == 3:
            f1, f2 = digits[0], digits[1:]
        else:
            return None
        return two_field_canonical(gene, f1, f2)

    m = re.match(r"^([A-Z]+)(\d{2,4})$", s, re.I)
    if not m:
        return None
    prefix, digits = m.group(1).upper(), m.group(2)
    gene_map = {
        "A": "A",
        "B": "B",
        "C": "C",
    }
    if prefix not in gene_map:
        return None
    gene = gene_map[prefix]
    if len(digits) == 4:
        f1, f2 = digits[:2], digits[2:]
    elif len(digits) == 3:
        f1, f2 = digits[0], digits[1:]
    else:
        return None
    return two_field_canonical(gene, f1, f2)


def normalize_allele(raw: str, alias_to_canonical: dict):
    original = raw
    s = raw.strip()
    if not s:
        return None

    if s in alias_to_canonical:
        return alias_to_canonical[s]
    su = s.upper()
    if su in alias_to_canonical:
        return alias_to_canonical[su]

    if "-" in s and re.search(r"DQA.*DQB", s, re.I):
        parts = re.findall(
            r"([A-Z0-9]+)\*(\d{2}):(\d{2})", s.upper()
        )
        if len(parts) == 2:
            return "-".join(
                two_field_canonical(g, f1, f2) for g, f1, f2 in parts
            )
        return None

    work = s.strip()
    work_upper = work.upper()

    # Already HLA-prefixed variants
    if work_upper.startswith("HLA-"):
        rest = work[4:] if work[4:5].isupper() else work[4:]
        # Fix missing *: HLA-A02:01, HLA-B14:02, HLA-DRB115:01
        m = re.match(
            r"^([A-Z0-9]+?)(\d{2}):(\d{2}(?:[A-Z])?)(.*)$",
            rest,
            re.I,
        )
        if m and "*" not in rest.split(":")[0]:
            gene_part, f1, f2, tail = m.groups()
            gene_part = gene_part.upper()
            # DRB115 -> DRB1 + 15
            if re.match(r"^DRB1?\d{2}$", gene_part, re.I):
                gm = re.match(r"^(DRB1?)(\d{2})$", gene_part, re.I)
                if gm and gm.group(1).upper() in ("DRB", "DRB1"):
                    gene_part = "DRB1"
                    f1 = gm.group(2)
            elif gene_part in ("A", "B", "C"):
                pass
            elif gene_part.startswith("DQB"):
                gene_part = "DQB1"
            suffix = ""
            if re.match(r"^[A-Z]$", f2[-1:]) and f2[-1].isalpha():
                suffix = f2[-1]
                f2 = f2[:-1]
            canon = two_field_canonical(gene_part, f1, f2, suffix)
            return canon

        # HLA-Cw*03:04 -> C
        rest_norm = re.sub(r"^C(W)?\*", "C*", rest, flags=re.I)
        m = re.match(
            r"^([A-Z0-9]+)\*(\d{2}):(\d{2}(?:[A-Z])?)",
            rest_norm,
            re.I,
        )
        if m:
            gene, f1, f2 = m.group(1).upper(), m.group(2), m.group(3)
            if gene == "CW":
                gene = "C"
            suffix = ""
            if f2[-1].isalpha():
                suffix = f2[-1]
                f2 = f2[:-1]
            return two_field_canonical(gene, f1, f2, suffix)

        m = re.match(r"^([A-Z0-9]+)(\d{2}):(\d{2})", rest, re.I)
        if m and "*" not in rest[:20]:
            gene, f1, f2 = m.group(1).upper(), m.group(2), m.group(3)
            if gene.startswith("DQB"):
                gene = "DQB1"
            return two_field_canonical(gene, f1, f2)

    # Gene*allele with optional HLA prefix missing
    work_norm = work_upper
    if not work_norm.startswith("HLA-"):
        # A*3201, B*0702, C*0802, DRB1*15:01
        m = re.match(
            r"^(?:HLA-)?([A-Z0-9]+)\*(\d{2,4})(?::(\d{2}(?:[A-Z])?))?$",
            work_norm,
        )
        if m:
            gene, d1, d2 = m.group(1), m.group(2), m.group(3)
            if gene == "CW":
                gene = "C"
            if d2 is None:
                if len(d1) == 4:
                    f1, f2 = d1[:2], d1[2:]
                elif len(d1) == 3:
                    f1, f2 = d1[0], d1[1:]
                else:
                    f1, f2 = d1, "01"
            else:
                f1, f2 = d1, d2
            suffix = ""
            if f2[-1].isalpha():
                suffix = f2[-1]
                f2 = f2[:-1]
            return two_field_canonical(gene, f1, f2, suffix)

        compact = parse_compact_allele(work_norm.replace("*", ""))
        if compact:
            return compact

        # DRB1*1101 without colon
        m = re.match(r"^(DRB1?)\*(\d{4})$", work_norm)
        if m:
            d = m.group(2)
            return two_field_canonical("DRB1", d[:2], d[2:])

    # HLA-A*full form
    m = re.match(
        r"^HLA-([A-Z0-9]+)\*(\d{2}):(\d{2}([A-Z])?)",
        work_upper.replace(" ", ""),
    )
    if m:
        gene, f1, f2, suffix = m.group(1), m.group(2), m.group(3), m.group(4) or ""
        if gene == "CW":
            gene = "C"
        return two_field_canonical(gene, f1, f2, suffix)

    return None


def load_source(path, payload_col):
    rows = []
    with open(path, newline="") as f:
        reader = csv.DictReader(f)
        for row in reader:
            raw = row["allele"]
            payload = float(row[payload_col])
            rows.append((raw, payload))
    return rows


def match_rows(rows_a, rows_b, alias_to_canonical):
    mapping = {}
    norm_a = []
    norm_b = []

    for raw, payload in rows_a:
        canon = normalize_allele(raw, alias_to_canonical)
        mapping[raw] = canon if canon else raw
        norm_a.append((raw, payload, canon))

    for raw, payload in rows_b:
        canon = normalize_allele(raw, alias_to_canonical)
        mapping[raw] = canon if canon else raw
        norm_b.append((raw, payload, canon))

    by_canon_a = defaultdict(list)
    by_canon_b = defaultdict(list)
    for item in norm_a:
        if item[2]:
            by_canon_a[item[2]].append(item)
    for item in norm_b:
        if item[2]:
            by_canon_b[item[2]].append(item)

    merged = []
    used_a = set()
    used_b = set()

    # Pass 1: exact string matches (case-sensitive as written in files)
    for i, (raw_a, pa, ca) in enumerate(norm_a):
        if not ca or i in used_a:
            continue
        for j, (raw_b, pb, cb) in enumerate(norm_b):
            if cb != ca or j in used_b:
                continue
            if raw_a == raw_b:
                merged.append(
                    {
                        "canonical_allele": ca,
                        "a_payload": pa,
                        "b_payload": pb,
                        "raw_a": raw_a,
                        "raw_b": raw_b,
                    }
                )
                used_a.add(i)
                used_b.add(j)
                break

    # Pass 2: remaining normalized matches (greedy)
    for i, (raw_a, pa, ca) in enumerate(norm_a):
        if not ca or i in used_a:
            continue
        for j, (raw_b, pb, cb) in enumerate(norm_b):
            if cb != ca or j in used_b:
                continue
            merged.append(
                {
                    "canonical_allele": ca,
                    "a_payload": pa,
                    "b_payload": pb,
                    "raw_a": raw_a,
                    "raw_b": raw_b,
                }
            )
            used_a.add(i)
            used_b.add(j)
            break

    naive_dropped = sum(1 for m in merged if m["raw_a"] != m["raw_b"])

    out_merged = [
        {
            "canonical_allele": m["canonical_allele"],
            "a_payload": m["a_payload"],
            "b_payload": m["b_payload"],
        }
        for m in sorted(merged, key=lambda x: x["canonical_allele"])
    ]

    return out_merged, len(merged), naive_dropped, mapping


def main():
    alias_to_canonical = load_aliases(ROOT / "gene_aliases.tsv")
    rows_a = load_source(ROOT / "source_a.csv", "payload_a")
    rows_b = load_source(ROOT / "source_b.csv", "payload_b")

    merged, joined_count, naive_dropped, mapping = match_rows(
        rows_a, rows_b, alias_to_canonical
    )

    # Fix mapping: use canonical string for all successfully normalized
    for raw in list(mapping.keys()):
        c = normalize_allele(raw, alias_to_canonical)
        if c:
            mapping[raw] = c

    out = {
        "merged": merged,
        "joined_count": joined_count,
        "naive_dropped_count": naive_dropped,
        "mapping": mapping,
    }

    out_path = ROOT / "output" / "submission.json"
    with open(out_path, "w") as f:
        json.dump(out, f, indent=2)
        f.write("\n")

    print(f"joined_count={joined_count} naive_dropped={naive_dropped}")
    for m in merged:
        print(m)


if __name__ == "__main__":
    main()
