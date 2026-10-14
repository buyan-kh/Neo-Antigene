#!/usr/bin/env python3
import csv
import json
import re
from pathlib import Path

ROOT = Path(__file__).parent

GENE_MAP = {
    "A": "A",
    "B": "B",
    "C": "C",
    "CW": "C",
    "DRB1": "DRB1",
    "DRB": "DRB1",
    "DQB1": "DQB1",
    "DQA1": "DQA1",
}


def load_aliases(path):
    alias_to_canonical = {}
    with open(path, newline="") as f:
        for row in csv.DictReader(f, delimiter="\t"):
            canon = row["canonical_symbol"].strip()
            for key in [row["prev_symbol"]] + row["synonyms"].split("|"):
                key = key.strip()
                if key:
                    alias_to_canonical[key.upper()] = canon
    return alias_to_canonical


def split_fields(digits):
    """Split a digit run into two-field first:second."""
    if len(digits) >= 4:
        return digits[:2], digits[2:4]
    if len(digits) == 3:
        return digits[0], digits[1:]
    if len(digits) == 2:
        return digits, "01"
    return digits, "01"


def format_allele(gene, f1, f2, suffix=""):
    gene = GENE_MAP.get(gene.upper(), gene.upper())
    if gene == "CW":
        gene = "C"
    return f"HLA-{gene}*{f1}:{f2}{suffix}"


def parse_haplotype(s):
    parts = s.split("-")
    if len(parts) != 2:
        return None
    left = canonicalize_raw(parts[0].strip(), {}, allow_alias=False)
    right = canonicalize_raw(parts[1].strip(), {}, allow_alias=False)
    if left and right:
        return f"{left}-{right}"
    return None


def canonicalize_raw(raw, aliases, allow_alias=True):
    s = raw.strip()
    upper = s.upper()

    if allow_alias and upper in aliases:
        return aliases[upper]

    if "-" in s and re.search(r"DQA|DQB", upper):
        hap = parse_haplotype(s)
        if hap:
            return hap

    # HLA-DRB115:01 (missing star between gene and allele)
    m = re.match(r"^HLA-DRB1(\d{2}):(\d{2})$", upper)
    if m:
        return format_allele("DRB1", m.group(1), m.group(2))

    # Standard HLA-GENE*ff:ss:... with optional trailing expression/group letter
    m = re.match(
        r"^HLA-(CW|[ABCDQ][A-Z0-9]*)\*(\d{2}):(\d{2})(?::\d{2})*(?::\d{2})*([A-ZNGQLS])?$",
        upper,
    )
    if m:
        gene, f1, f2, suffix = m.group(1), m.group(2), m.group(3), m.group(4) or ""
        if suffix == "G" and len(re.findall(r":\d{2}", upper)) >= 2:
            suffix = ""
        return format_allele(gene, f1, f2, suffix)

    # HLA-A02:01 / HLA-B14:02 (no star)
    m = re.match(r"^HLA-([ABC])(\d{2}):(\d{2})([A-ZNGQLS])?$", upper)
    if m:
        return format_allele(m.group(1), m.group(2), m.group(3), m.group(4) or "")

    m = re.match(r"^HLA-([BC])(\d{2}):(\d{2})$", upper)
    if m:
        return format_allele(m.group(1), m.group(2), m.group(3))

    # Gene*allele compact: DRB1*15:01, DRB1*1101, A*29:02, C*0802
    m = re.match(
        r"^(DRB1|DQB1|DQA1|DPB1|CW|A|B|C)\*(\d{2}):(\d{2})([A-ZNGQLS])?$",
        upper,
    )
    if m:
        return format_allele(m.group(1), m.group(2), m.group(3), m.group(4) or "")

    m = re.match(r"^(DRB1|DQB1|DQA1|CW|A|B|C)\*(\d{3,4})([A-ZNGQLS])?$", upper)
    if m:
        f1, f2 = split_fields(m.group(2))
        return format_allele(m.group(1), f1, f2, m.group(3) or "")

    # Compact without star: B1801, A6801, DQB10301, A0206
    m = re.match(r"^(DQB1|DRB1|DQA1|A|B|C)(\d{3,4})([A-ZNGQLS])?$", upper)
    if m:
        f1, f2 = split_fields(m.group(2))
        return format_allele(m.group(1), f1, f2, m.group(3) or "")

    m = re.match(r"^(A|B|C)(\d{4})$", upper)
    if m:
        f1, f2 = split_fields(m.group(2))
        return format_allele(m.group(1), f1, f2)

    return None


def canonicalize(raw, aliases):
    s = raw.strip()
    if s.upper() in aliases:
        return aliases[s.upper()]
    result = canonicalize_raw(s, aliases)
    if result:
        return result
    raise ValueError(f"Cannot canonicalize: {raw!r}")


def load_csv(path, payload_col):
    rows = {}
    raw_names = []
    with open(path, newline="") as f:
        reader = csv.DictReader(f)
        for row in reader:
            allele = row["allele"]
            raw_names.append(allele)
            payload = float(row[payload_col])
            rows.setdefault(allele, payload)
    return rows, raw_names


def main():
    aliases = load_aliases(ROOT / "gene_aliases.tsv")
    a_rows, a_raw_list = load_csv(ROOT / "source_a.csv", "payload_a")
    b_rows, b_raw_list = load_csv(ROOT / "source_b.csv", "payload_b")

    all_raw = list(dict.fromkeys(a_raw_list + b_raw_list))

    mapping = {raw: canonicalize(raw, aliases) for raw in all_raw}

    common = sorted(
        {mapping[r] for r in a_rows} & {mapping[r] for r in b_rows}
    )

    def pick_payload_pair(c):
        a_raws = [r for r in a_rows if mapping[r] == c]
        b_raws = [r for r in b_rows if mapping[r] == c]
        exact = [(ra, rb) for ra in a_raws for rb in b_raws if ra == rb]
        if exact:
            exact.sort(
                key=lambda pair: (
                    0 if pair[0].upper().startswith("HLA") else 1,
                    -len(pair[0]),
                    pair[0],
                )
            )
            ra, rb = exact[0]
            return a_rows[ra], b_rows[rb], a_raws, b_raws
        if len(a_raws) == 1 and len(b_raws) == 1:
            return a_rows[a_raws[0]], b_rows[b_raws[0]], a_raws, b_raws
        raise ValueError(f"Ambiguous payloads for {c}: A={a_raws} B={b_raws}")

    naive_dropped = 0
    merged = []
    for c in common:
        a_payload, b_payload, a_raws, b_raws = pick_payload_pair(c)
        has_identical = any(ra == rb for ra in a_raws for rb in b_raws)
        if not has_identical:
            naive_dropped += 1
        merged.append(
            {
                "canonical_allele": c,
                "a_payload": a_payload,
                "b_payload": b_payload,
            }
        )

    out = {
        "merged": merged,
        "joined_count": len(common),
        "naive_dropped_count": naive_dropped,
        "mapping": mapping,
    }

    out_path = ROOT / "output" / "submission.json"
    with open(out_path, "w") as f:
        json.dump(out, f, indent=2)

    print(f"joined: {len(common)}")
    print(f"naive_dropped: {naive_dropped}")


if __name__ == "__main__":
    main()
