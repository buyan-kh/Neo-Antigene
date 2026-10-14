#!/usr/bin/env python3
import csv
import json
import re
from pathlib import Path

ROOT = Path(__file__).resolve().parent
GENES = [
    "DRB1", "DRB3", "DRB4", "DRB5",
    "DQB1", "DQA1", "DPB1", "DPA1",
    "DMA", "DMB", "DOA", "DOB",
    "MICA", "MICB", "TAP1", "TAP2",
    "CW", "A", "B", "C", "E", "F", "G",
]


def load_aliases(path):
    aliases = {}
    with open(path, newline="") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            canonical = row["canonical_symbol"].strip()
            names = [canonical, row["prev_symbol"].strip()]
            syn = row.get("synonyms") or ""
            names.extend(s.strip() for s in syn.split("|") if s.strip())
            for name in names:
                aliases[name.upper()] = canonical
    return aliases


def split_haplotype(s):
    m = re.match(
        r"^((?:HLA-)?[A-Z0-9]+\*[0-9A-Z:]+)-((?:HLA-)?[A-Z0-9]+\*[0-9A-Z:]+)$",
        s,
        flags=re.IGNORECASE,
    )
    if m:
        return m.group(1), m.group(2)
    return None


def parse_molecular(s):
    text = s.strip().upper()
    if text.startswith("HLA-"):
        text = text[4:]
    gene = None
    rest = None
    for g in GENES:
        if text.startswith(g):
            gene = "C" if g == "CW" else g
            rest = text[len(g):]
            break
    if gene is None or rest is None:
        return None
    if rest.startswith("*"):
        rest = rest[1:]
    if not rest:
        return None
    suffix = ""
    m = re.search(r"([A-Z]+)$", rest)
    if m:
        suffix = m.group(1)
        rest = rest[: m.start()]
    if not rest:
        return None
    if ":" in rest:
        fields = rest.split(":")
        if not all(re.fullmatch(r"\d+", fld) for fld in fields):
            return None
    else:
        if not rest.isdigit() or len(rest) < 4 or len(rest) % 2 != 0:
            return None
        fields = [rest[i : i + 2] for i in range(0, len(rest), 2)]
    if len(fields) < 2:
        return None
    return f"HLA-{gene}*{fields[0]}:{fields[1]}{suffix}"


def canonicalize(raw, aliases):
    s = raw.strip()
    hap = split_haplotype(s)
    if hap:
        a = canonicalize(hap[0], aliases)
        b = canonicalize(hap[1], aliases)
        if a and b:
            return f"{a}-{b}"
        return None
    parsed = parse_molecular(s)
    if parsed:
        return parsed
    return aliases.get(s.upper())


def load_source(path, payload_key):
    rows = []
    with open(path, newline="") as f:
        reader = csv.DictReader(f)
        for row in reader:
            rows.append((row["allele"], float(row[payload_key])))
    return rows


def main():
    aliases = load_aliases(ROOT / "gene_aliases.tsv")
    src_a = load_source(ROOT / "source_a.csv", "payload_a")
    src_b = load_source(ROOT / "source_b.csv", "payload_b")

    mapping = {}
    a_by_canon = {}
    b_by_canon = {}
    a_raw_by_canon = {}
    b_raw_by_canon = {}
    conflicts = []

    for raw, payload in src_a:
        canon = canonicalize(raw, aliases)
        if not canon:
            raise SystemExit(f"unmapped A: {raw!r}")
        mapping[raw] = canon
        if canon in a_by_canon and a_by_canon[canon] != payload:
            conflicts.append(("A", canon, a_raw_by_canon[canon], raw, a_by_canon[canon], payload))
        a_by_canon[canon] = payload
        a_raw_by_canon.setdefault(canon, raw)

    for raw, payload in src_b:
        canon = canonicalize(raw, aliases)
        if not canon:
            raise SystemExit(f"unmapped B: {raw!r}")
        mapping[raw] = canon
        if canon in b_by_canon and b_by_canon[canon] != payload:
            conflicts.append(("B", canon, b_raw_by_canon[canon], raw, b_by_canon[canon], payload))
        b_by_canon[canon] = payload
        b_raw_by_canon.setdefault(canon, raw)

    shared = sorted(set(a_by_canon) & set(b_by_canon))
    merged = [
        {
            "canonical_allele": name,
            "a_payload": a_by_canon[name],
            "b_payload": b_by_canon[name],
        }
        for name in shared
    ]

    naive_dropped = 0
    for name in shared:
        if a_raw_by_canon[name] != b_raw_by_canon[name]:
            naive_dropped += 1

    print("CONFLICTS", conflicts)
    print("A only", sorted(set(a_by_canon) - set(b_by_canon)))
    print("B only", sorted(set(b_by_canon) - set(a_by_canon)))
    print("SHARED", shared)
    print("joined", len(shared), "naive_dropped", naive_dropped)
    print("MAPPING")
    for k, v in mapping.items():
        print(repr(k), "->", v)

    out = {
        "merged": merged,
        "joined_count": len(shared),
        "naive_dropped_count": naive_dropped,
        "mapping": mapping,
    }
    out_path = ROOT / "output" / "submission.json"
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with open(out_path, "w") as f:
        json.dump(out, f, indent=2)
        f.write("\n")
    print("wrote", out_path)


if __name__ == "__main__":
    main()
