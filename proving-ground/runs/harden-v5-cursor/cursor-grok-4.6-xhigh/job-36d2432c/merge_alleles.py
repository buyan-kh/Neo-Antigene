#!/usr/bin/env python3
import csv
import json
import re
from pathlib import Path

ROOT = Path(__file__).resolve().parent

GENES = [
    "DRB1", "DRB3", "DRB4", "DRB5", "DQA1", "DQB1", "DPA1", "DPB1",
    "DRA", "DMA", "DMB", "DOA", "DOB", "CW", "A", "B", "C", "E", "F", "G",
]
GENE_RE = re.compile(r"^(" + "|".join(GENES) + r")")
GENE_MAP = {"CW": "C"}


def load_aliases(path):
    alias_to_canonical = {}
    with open(path, newline="") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            canon = row["canonical_symbol"].strip()
            prev = (row.get("prev_symbol") or "").strip()
            syns = (row.get("synonyms") or "").strip()
            names = []
            if prev:
                names.append(prev)
            if syns:
                names.extend(s.strip() for s in syns.split("|") if s.strip())
            for name in names:
                alias_to_canonical[name] = canon
                alias_to_canonical[name.upper()] = canon
                alias_to_canonical[name.lower()] = canon
    return alias_to_canonical


def parse_fields(rest):
    """Parse allele numeric fields + optional trailing letter suffix."""
    rest = rest.strip()
    suffix = ""
    m = re.search(r"([A-Za-z]+)$", rest)
    if m:
        suffix = m.group(1).upper()
        rest = rest[: m.start()]
    rest = rest.strip()
    if not rest:
        return None, suffix
    if ":" in rest:
        parts = [p for p in rest.split(":") if p != ""]
        fields = []
        for p in parts:
            if not p.isdigit():
                return None, suffix
            fields.append(int(p))
        return fields, suffix
    if not rest.isdigit():
        return None, suffix
    if len(rest) % 2 != 0:
        # odd length: pad left so 301 -> 03:01 style serological-ish, or 3 -> 03
        rest = rest.zfill(len(rest) + 1)
    fields = [int(rest[i : i + 2]) for i in range(0, len(rest), 2)]
    return fields, suffix


def normalize(raw, aliases):
    original = raw
    s = raw.strip()
    if not s:
        raise ValueError(f"empty allele: {original!r}")

    # serological / historical aliases (exact, after strip, case-insensitive)
    if s in aliases:
        return aliases[s]
    if s.upper() in aliases:
        return aliases[s.upper()]

    t = s.upper()
    if t.startswith("HLA-"):
        t = t[4:]
    t = t.strip()

    if "*" in t:
        gene, rest = t.split("*", 1)
        gene = gene.strip()
        rest = rest.strip()
    else:
        m = GENE_RE.match(t)
        if not m:
            raise ValueError(f"cannot find gene in {original!r} -> {t!r}")
        gene = m.group(1)
        rest = t[m.end() :].strip()

    gene = GENE_MAP.get(gene, gene)
    fields, suffix = parse_fields(rest)
    if not fields or len(fields) < 1:
        raise ValueError(f"cannot parse fields from {original!r} rest={rest!r}")
    if len(fields) == 1:
        # single field: treat as allele group only — not two-field unless alias handled it
        raise ValueError(f"only one field in {original!r}: {fields}")
    f1, f2 = fields[0], fields[1]
    return f"HLA-{gene}*{f1:02d}:{f2:02d}{suffix}"


def load_source(path, payload_key):
    rows = []
    with open(path, newline="") as f:
        reader = csv.DictReader(f)
        for row in reader:
            allele = row["allele"]
            payload = float(row[payload_key])
            rows.append((allele, payload))
    return rows


def main():
    aliases = load_aliases(ROOT / "gene_aliases.tsv")
    print("ALIASES:")
    for k, v in sorted(aliases.items(), key=lambda kv: (kv[1], kv[0])):
        print(f"  {k!r} -> {v}")

    a_rows = load_source(ROOT / "source_a.csv", "payload_a")
    b_rows = load_source(ROOT / "source_b.csv", "payload_b")

    print("\nSOURCE A:")
    a_map = {}  # canonical -> (raw, payload)
    mapping = {}
    for raw, payload in a_rows:
        canon = normalize(raw, aliases)
        mapping[raw] = canon
        print(f"  {raw!r:30} -> {canon:20}  payload={payload}")
        if canon in a_map:
            old_raw, old_payload = a_map[canon]
            if old_payload != payload:
                print(f"    CONFLICT A {canon}: {old_raw}={old_payload} vs {raw}={payload}")
        a_map[canon] = (raw, payload)

    print("\nSOURCE B:")
    b_map = {}
    for raw, payload in b_rows:
        canon = normalize(raw, aliases)
        mapping[raw] = canon
        print(f"  {raw!r:30} -> {canon:20}  payload={payload}")
        if canon in b_map:
            old_raw, old_payload = b_map[canon]
            if old_payload != payload:
                print(f"    CONFLICT B {canon}: {old_raw}={old_payload} vs {raw}={payload}")
        b_map[canon] = (raw, payload)

    shared = sorted(set(a_map) & set(b_map))
    a_only = sorted(set(a_map) - set(b_map))
    b_only = sorted(set(b_map) - set(a_map))
    print("\nSHARED:")
    naive_dropped = 0
    merged = []
    for canon in shared:
        raw_a, pa = a_map[canon]
        raw_b, pb = b_map[canon]
        same = raw_a == raw_b
        if not same:
            naive_dropped += 1
        print(f"  {canon:20} A={raw_a!r:30} B={raw_b!r:30} same={same} pa={pa} pb={pb}")
        merged.append({
            "canonical_allele": canon,
            "a_payload": pa,
            "b_payload": pb,
        })

    print("\nA ONLY:", a_only)
    print("B ONLY:", b_only)
    print(f"\njoined_count={len(shared)}")
    print(f"naive_dropped_count={naive_dropped}")
    print(f"mapping size={len(mapping)}")

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
    print(f"wrote {out_path}")


if __name__ == "__main__":
    main()
