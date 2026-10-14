#!/usr/bin/env python3
import csv
import json
import re
from collections import OrderedDict
from pathlib import Path

ROOT = Path(__file__).resolve().parent

GENES = [
    "DQB1",
    "DQA1",
    "DPB1",
    "DPA1",
    "DRB1",
    "DRB3",
    "DRB4",
    "DRB5",
    "DRA",
    "DMA",
    "DMB",
    "DOA",
    "DOB",
    "MICA",
    "MICB",
    "TAP1",
    "TAP2",
    "E",
    "F",
    "G",
    "A",
    "B",
    "C",
]


def load_aliases(path):
    aliases = {}
    with open(path, newline="") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            canonical = row["canonical_symbol"].strip()
            keys = [canonical]
            prev = (row.get("prev_symbol") or "").strip()
            if prev:
                keys.append(prev)
            syn = (row.get("synonyms") or "").strip()
            if syn:
                keys.extend(s.strip() for s in syn.split("|") if s.strip())
            for k in keys:
                aliases[k.lower()] = canonical
    return aliases


def parse_allele(name, aliases):
    original = name
    s = name.strip()
    if not s:
        raise ValueError(f"empty allele name: {original!r}")

    alias_hit = aliases.get(s.lower())
    if alias_hit:
        return alias_hit, "alias"

    t = s.upper()
    if t.startswith("HLA-"):
        t = t[4:]

    gene = None
    rest = None
    if t.startswith("CW"):
        gene = "C"
        rest = t[2:]
    else:
        for g in GENES:
            if t.startswith(g):
                gene = g
                rest = t[len(g) :]
                break
    if gene is None:
        raise ValueError(f"cannot parse gene from {original!r}")

    if rest.startswith("*"):
        rest = rest[1:]
    rest = rest.strip()
    if not rest:
        raise ValueError(f"missing fields in {original!r}")

    suffix = ""
    m = re.fullmatch(r"([0-9:]+?)([A-Z])", rest)
    if m:
        rest = m.group(1)
        suffix = m.group(2)

    if ":" in rest:
        fields = rest.split(":")
        if any(not fld.isdigit() for fld in fields):
            raise ValueError(f"bad fields in {original!r}: {rest}")
    else:
        if not rest.isdigit() or len(rest) < 4 or len(rest) % 2 != 0:
            raise ValueError(f"bad compact fields in {original!r}: {rest}")
        fields = [rest[i : i + 2] for i in range(0, len(rest), 2)]

    if len(fields) < 2:
        raise ValueError(f"need two fields in {original!r}")

    f1 = fields[0].zfill(2)
    f2 = fields[1].zfill(2)
    canonical = f"HLA-{gene}*{f1}:{f2}{suffix}"
    return canonical, "parsed"


def load_source(path, payload_key):
    rows = []
    with open(path, newline="") as f:
        reader = csv.reader(f)
        header = next(reader)
        assert header[0] == "allele"
        for row in reader:
            if not row:
                continue
            rows.append({"allele": row[0], "payload": row[1]})
    return rows


def main():
    aliases = load_aliases(ROOT / "gene_aliases.tsv")
    source_a = load_source(ROOT / "source_a.csv", "payload_a")
    source_b = load_source(ROOT / "source_b.csv", "payload_b")

    mapping = OrderedDict()
    a_by_canon = OrderedDict()
    b_by_canon = OrderedDict()
    a_raw_by_canon = OrderedDict()
    b_raw_by_canon = OrderedDict()

    print("=== SOURCE A ===")
    for row in source_a:
        raw = row["allele"]
        canon, how = parse_allele(raw, aliases)
        mapping[raw] = canon
        print(f"A  {raw!r:30} -> {canon:22} ({how}) payload={row['payload']}")
        if canon in a_by_canon:
            if a_by_canon[canon] != row["payload"]:
                raise SystemExit(f"conflicting A payloads for {canon}")
            if a_raw_by_canon[canon] != raw:
                print(f"  NOTE: extra A raw name {raw!r} for {canon}")
        else:
            a_by_canon[canon] = row["payload"]
            a_raw_by_canon[canon] = raw

    print("\n=== SOURCE B ===")
    for row in source_b:
        raw = row["allele"]
        canon, how = parse_allele(raw, aliases)
        mapping[raw] = canon
        print(f"B  {raw!r:30} -> {canon:22} ({how}) payload={row['payload']}")
        if canon in b_by_canon:
            if b_by_canon[canon] != row["payload"]:
                raise SystemExit(f"conflicting B payloads for {canon}")
            if b_raw_by_canon[canon] != raw:
                print(f"  NOTE: extra B raw name {raw!r} for {canon}")
        else:
            b_by_canon[canon] = row["payload"]
            b_raw_by_canon[canon] = raw

    joined = [c for c in a_by_canon if c in b_by_canon]
    a_only = [c for c in a_by_canon if c not in b_by_canon]
    b_only = [c for c in b_by_canon if c not in a_by_canon]

    print("\n=== JOINED ===")
    naive_same = []
    naive_diff = []
    merged = []
    for c in joined:
        ra = a_raw_by_canon[c]
        rb = b_raw_by_canon[c]
        same = ra == rb
        (naive_same if same else naive_diff).append(c)
        print(
            f"{c:22} A={a_by_canon[c]:>8} B={b_by_canon[c]:>8} "
            f"rawA={ra!r} rawB={rb!r} identical={same}"
        )
        merged.append(
            {
                "canonical_allele": c,
                "a_payload": json.loads(a_by_canon[c]),
                "b_payload": json.loads(b_by_canon[c]),
            }
        )

    print("\n=== A ONLY ===")
    for c in a_only:
        print(f"{c:22} raw={a_raw_by_canon[c]!r} payload={a_by_canon[c]}")
    print("\n=== B ONLY ===")
    for c in b_only:
        print(f"{c:22} raw={b_raw_by_canon[c]!r} payload={b_by_canon[c]}")

    print("\n=== COUNTS ===")
    print("unique A canonical", len(a_by_canon))
    print("unique B canonical", len(b_by_canon))
    print("joined", len(joined))
    print("identical raw strings among joined", len(naive_same))
    print("different raw strings among joined", len(naive_diff))
    print("mapping keys", len(mapping))
    print("a only", len(a_only))
    print("b only", len(b_only))

    # Also compute naive exact-string inner join on raw names
    a_raw_set = {row["allele"] for row in source_a}
    b_raw_set = {row["allele"] for row in source_b}
    exact = a_raw_set & b_raw_set
    print("exact raw string overlap", len(exact), sorted(exact))

    out = {
        "merged": merged,
        "joined_count": len(joined),
        "naive_dropped_count": len(naive_diff),
        "mapping": mapping,
    }
    out_path = ROOT / "output" / "submission.json"
    with open(out_path, "w") as f:
        json.dump(out, f, indent=2)
        f.write("\n")
    print("wrote", out_path)


if __name__ == "__main__":
    main()
