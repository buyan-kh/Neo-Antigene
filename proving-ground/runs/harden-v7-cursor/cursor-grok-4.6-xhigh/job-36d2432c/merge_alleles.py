#!/usr/bin/env python3
import csv
import json
import re
from pathlib import Path

ROOT = Path(__file__).resolve().parent

GENES = [
    "DQA1", "DQB1", "DPA1", "DPB1",
    "DRB1", "DRB3", "DRB4", "DRB5",
    "MICA", "MICB", "TAP1", "TAP2",
    "DRA", "DMA", "DMB", "DOA", "DOB",
    "CW",
    "A", "B", "C", "E", "F", "G", "H", "J", "K", "L",
]


def load_aliases(path):
    alias_to_canonical = {}
    with open(path, newline="") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            canonical = row["canonical_symbol"].strip()
            keys = [canonical, row["prev_symbol"].strip()]
            syn = row.get("synonyms") or ""
            keys.extend(s.strip() for s in syn.split("|") if s.strip())
            for key in keys:
                alias_to_canonical[norm_key(key)] = canonical
    return alias_to_canonical


def norm_key(name):
    s = name.strip().upper()
    s = s.replace("HLA-", "")
    s = s.replace("*", "")
    s = s.replace(":", "")
    s = s.replace(" ", "")
    return s


def extract_suffix(token):
    m = re.search(r"([A-Z]+)$", token)
    if m:
        return token[: m.start()], m.group(1)
    return token, ""


def split_fields(rest):
    rest = rest.strip()
    if not rest:
        return [], ""
    rest, suffix = extract_suffix(rest)
    rest = rest.strip(":")
    if ":" in rest:
        fields = [p for p in rest.split(":") if p]
    else:
        digits = re.sub(r"\D", "", rest)
        if len(digits) % 2 == 1 and len(digits) > 1:
            # odd length: first field may be 1 digit (rare); pad pairs from the right
            fields = [digits[:1].zfill(2)] + [digits[i : i + 2] for i in range(1, len(digits), 2)]
        else:
            if len(digits) == 1 or len(digits) == 2:
                fields = [digits.zfill(2)]
            else:
                fields = [digits[i : i + 2] for i in range(0, len(digits), 2)]
    fields = [f.zfill(2) if f.isdigit() else f for f in fields]
    return fields, suffix


def split_gene(body):
    body = body.upper()
    if "*" in body:
        gene, rest = body.split("*", 1)
        return gene.replace(" ", ""), rest
    for gene in GENES:
        if body.startswith(gene) and (len(body) == len(gene) or not body[len(gene)].isalpha()):
            return gene, body[len(gene) :]
    m = re.match(r"([A-Z]+)", body)
    if m:
        return m.group(1), body[m.end() :]
    return body, ""


def canonical_gene(gene):
    gene = gene.upper().replace("HLA-", "")
    if gene == "CW":
        return "C"
    return gene


def is_haplotype(name):
    s = name.strip()
    s_wo = re.sub(r"^HLA-", "", s, flags=re.I)
    # hyphen not part of HLA- prefix, with allele-like material on both sides
    return bool(re.search(r"[A-Z0-9*]\s*-\s*(?:HLA-)?[A-Z]", s_wo, re.I))


def parse_one(name, aliases):
    raw_stripped = name.strip()
    key = norm_key(raw_stripped)

    body = re.sub(r"^HLA-", "", raw_stripped, flags=re.I).strip()
    gene, rest = split_gene(body)
    gene = canonical_gene(gene)
    fields, suffix = split_fields(rest)

    if len(fields) >= 2:
        return f"HLA-{gene}*{fields[0]}:{fields[1]}{suffix}"

    # serological / short alias
    if key in aliases:
        return aliases[key]

    sero_keys = [norm_key(gene + (fields[0] if fields else rest))]
    if gene == "C":
        sero_keys.append(norm_key("CW" + (fields[0] if fields else rest)))
        sero_keys.append(norm_key("Cw" + (fields[0].lstrip("0") if fields else rest)))
    for sk in sero_keys:
        if sk in aliases:
            return aliases[sk]

    if fields:
        return f"HLA-{gene}*{fields[0]}{suffix}"
    return f"HLA-{gene}"


def to_canonical(name, aliases):
    stripped = name.strip()
    if is_haplotype(stripped):
        body = re.sub(r"^HLA-", "", stripped, flags=re.I)
        parts = re.split(r"-(?=(?:HLA-)?[A-Z])", body, maxsplit=1)
        return "-".join(parse_one(p, aliases) for p in parts)
    return parse_one(stripped, aliases)


def load_source(path, allele_col, payload_col):
    rows = []
    with open(path, newline="") as f:
        reader = csv.DictReader(f)
        for row in reader:
            rows.append((row[allele_col], float(row[payload_col])))
    return rows


def main():
    aliases = load_aliases(ROOT / "gene_aliases.tsv")
    src_a = load_source(ROOT / "source_a.csv", "allele", "payload_a")
    src_b = load_source(ROOT / "source_b.csv", "allele", "payload_b")

    mapping = {}
    a_by_canon = {}
    a_raw_by_canon = {}
    for raw, payload in src_a:
        canon = to_canonical(raw, aliases)
        mapping[raw] = canon
        a_by_canon[canon] = payload
        a_raw_by_canon[canon] = raw

    b_by_canon = {}
    b_raw_by_canon = {}
    for raw, payload in src_b:
        canon = to_canonical(raw, aliases)
        mapping[raw] = canon
        b_by_canon[canon] = payload
        b_raw_by_canon[canon] = raw

    joined = sorted(set(a_by_canon) & set(b_by_canon))
    naive_dropped = 0
    merged = []
    print("=== JOINED ===")
    for canon in joined:
        raw_a = a_raw_by_canon[canon]
        raw_b = b_raw_by_canon[canon]
        same = raw_a == raw_b
        if not same:
            naive_dropped += 1
        rec = {
            "canonical_allele": canon,
            "a_payload": a_by_canon[canon],
            "b_payload": b_by_canon[canon],
        }
        merged.append(rec)
        flag = "SAME" if same else "DIFF"
        print(f"{flag:4} {canon:22} A={raw_a!r:30} B={raw_b!r:30} pa={a_by_canon[canon]} pb={b_by_canon[canon]}")

    print("\n=== A ONLY ===")
    for canon in sorted(set(a_by_canon) - set(b_by_canon)):
        print(f"  {canon:22} raw={a_raw_by_canon[canon]!r}")

    print("\n=== B ONLY ===")
    for canon in sorted(set(b_by_canon) - set(a_by_canon)):
        print(f"  {canon:22} raw={b_raw_by_canon[canon]!r}")

    print("\n=== MAPPING ===")
    for raw, canon in mapping.items():
        print(f"  {raw!r:35} -> {canon}")

    print(f"\njoined_count={len(joined)}")
    print(f"naive_dropped_count={naive_dropped}")
    print(f"mapping_size={len(mapping)}")
    print(f"a unique canon={len(a_by_canon)} b unique canon={len(b_by_canon)}")

    out = {
        "merged": merged,
        "joined_count": len(joined),
        "naive_dropped_count": naive_dropped,
        "mapping": mapping,
    }
    out_path = ROOT / "output" / "submission.json"
    out_path.parent.mkdir(exist_ok=True)
    with open(out_path, "w") as f:
        json.dump(out, f, indent=2)
        f.write("\n")
    print(f"wrote {out_path}")


if __name__ == "__main__":
    main()
