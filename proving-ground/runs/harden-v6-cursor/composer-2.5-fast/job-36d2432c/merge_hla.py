#!/usr/bin/env python3
import csv
import json
import re
from pathlib import Path

ROOT = Path(__file__).resolve().parent


def load_gene_aliases(path):
    """Map any alias (prev_symbol, synonyms) -> canonical_symbol from TSV."""
    alias_to_canonical = {}
    with open(path, newline="") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            canon = row["canonical_symbol"].strip()
            alias_to_canonical[canon.upper()] = canon
            prev = row["prev_symbol"].strip()
            if prev:
                alias_to_canonical[prev.upper()] = canon
            syns = row["synonyms"].strip()
            if syns:
                for s in syns.split("|"):
                    s = s.strip()
                    if s:
                        alias_to_canonical[s.upper()] = canon
    return alias_to_canonical


GENE_ALIASES = None


def parse_allele_fields(after_star: str):
    """Split allele designation into fields; preserve trailing suffix on last field."""
    after_star = after_star.strip()
    parts = after_star.split(":")
    fields = []
    for i, p in enumerate(parts):
        p = p.strip()
        if not p:
            continue
        if i == 0 and len(p) >= 2:
            fields.append(p[:2])
            rest = p[2:]
            if rest:
                fields.append(rest)
        else:
            fields.append(p)
    return fields


def suffix_from_field(field: str):
    m = re.match(r"^(\d+)([A-Za-z]+)$", field)
    if m:
        return m.group(1), m.group(2).upper()
    return field, ""


def to_two_field(gene: str, fields: list) -> str:
    """Build HLA-<gene>*dd:dd with optional suffix on second field."""
    if len(fields) < 2:
        raise ValueError(f"need at least 2 fields, got {fields}")
    f0 = fields[0].zfill(2) if fields[0].isdigit() else fields[0]
    f1, suf = suffix_from_field(fields[1])
    if f1.isdigit():
        f1 = f1.zfill(2)
    second = f1 + suf
    return f"HLA-{gene}*{f0}:{second}"


def normalize_raw(raw: str) -> str:
    global GENE_ALIASES
    original = raw
    s = raw.strip()

    if s in GENE_ALIASES:
        key = s.upper()
        if key in GENE_ALIASES:
            return GENE_ALIASES[key]

    key = s.upper()
    if key in GENE_ALIASES:
        return GENE_ALIASES[key]

    work = s
    if work.lower().startswith("hla-"):
        work = work[4:]

    work = re.sub(r"^HLA-", "", work, flags=re.I)

    combo = re.match(
        r"^(DQA1\*[\d:]+)-DQB1\*([\d:]+)$", work, re.I
    )
    if combo:
        dqb = combo.group(2)
        fields = parse_allele_fields(dqb.replace("*", ""))
        if "*" in combo.group(1):
            pass
        dqb_part = combo.group(2)
        if ":" not in dqb_part and len(re.sub(r"\D", "", dqb_part)) >= 4:
            digits = re.sub(r"[^0-9]", "", dqb_part)
            dqb_part = f"{digits[:2]}:{digits[2:4]}"
        return f"HLA-DQA1*01:02-DQB1*{dqb_part}" if "01:02" in combo.group(1) else normalize_raw(
            f"DQA1*01:02-DQB1*{dqb_part}"
        )

    combo2 = re.match(r"^DQA1\*01:02-DQB1\*([\d:]+)$", work, re.I)
    if combo2:
        dqb = combo2.group(1)
        fields = parse_allele_fields(dqb)
        return to_two_field("DQA1*01:02-DQB1", fields) if False else _combo_canonical(dqb)

    if re.match(r"^DQA1\*01:02-DQB1\*", work, re.I):
        m = re.match(r"^DQA1\*01:02-DQB1\*(.+)$", work, re.I)
        dqb_rest = m.group(1)
        fields = parse_allele_fields(dqb_rest)
        dqb2 = to_two_field("DQB1", fields).replace("HLA-DQB1*", "")
        return f"HLA-DQA1*01:02-DQB1*{dqb2}"

    work_upper = work.upper()

    m = re.match(r"^(CW?\*?|C\*?|A\*?|B\*?|DRB1\*?|DQB1\*?)(.*)$", work, re.I)
    if not m:
        m2 = re.match(r"^(DRB1|DQB1|A|B|C|CW?)(\d{2,4})$", work, re.I)
        if m2:
            gene_raw = m2.group(1).upper()
            digits = m2.group(2)
            if gene_raw.startswith("C") and gene_raw != "C":
                gene = "C"
            elif gene_raw == "DRB1":
                gene = "DRB1"
            elif gene_raw == "DQB1":
                gene = "DQB1"
            else:
                gene = gene_raw.replace("CW", "C")
            if len(digits) == 4:
                after = f"{digits[:2]}:{digits[2:]}"
            elif len(digits) == 2:
                after = f"{digits}:01"
            else:
                after = digits
            fields = parse_allele_fields(after)
            return to_two_field(gene, fields)

        m3 = re.match(r"^(A|B|C)(\d{4})$", work, re.I)
        if m3:
            gene = m3.group(1).upper()
            d = m3.group(2)
            fields = parse_allele_fields(f"{d[:2]}:{d[2:]}")
            return to_two_field(gene, fields)

        raise ValueError(f"cannot parse: {original!r}")

    gene_part = m.group(1).upper().replace("CW", "C").replace("*", "")
    rest = m.group(2).strip()

    if gene_part == "C" and rest.startswith("W"):
        rest = rest[1:]
        if rest.startswith("*"):
            rest = rest[1:]

    gene_map = {"A": "A", "B": "B", "C": "C", "DRB1": "DRB1", "DQB1": "DQB1"}
    gene = gene_map.get(gene_part.replace("W", ""), gene_part)

    if not rest:
        raise ValueError(f"no allele rest: {original!r}")

    rest = rest.lstrip("*")

    m_hla_style = re.match(r"^(\d{2}):(\d{2,}.*)$", rest)
    if m_hla_style and "*" not in work and gene_part in ("A", "B", "C", "DRB1", "DQB1"):
        rest = f"{m_hla_style.group(1)}:{m_hla_style.group(2)}"

    m_embed = re.match(r"^(\d{2})(\d{2,})$", rest.replace(":", ""))
    if ":" not in rest and m_embed and len(rest.replace(":", "")) >= 4:
        digits = rest.replace(":", "").replace(" ", "")
        if len(digits) >= 4:
            rest = f"{digits[:2]}:{digits[2:]}"

    m_drb = re.match(r"^(\d{2})(\d{2}.*)$", rest)
    if gene == "DRB1" and ":" not in rest and m_drb:
        rest = f"{m_drb.group(1)}:{m_drb.group(2)}"

    fields = parse_allele_fields(rest)
    return to_two_field(gene, fields)


def _combo_canonical(dqb_rest):
    fields = parse_allele_fields(dqb_rest)
    dqb2 = to_two_field("DQB1", fields).replace("HLA-DQB1*", "")
    return f"HLA-DQA1*01:02-DQB1*{dqb2}"


def normalize_with_hla_prefix(raw: str) -> str:
    s = raw.strip()
    if not s.lower().startswith("hla-"):
        inner = s
    else:
        inner = s[4:]

    m = re.match(r"^([A-Z0-9]+)(\d{2}):(\d{2,}.*)$", inner, re.I)
    if m and m.group(1).upper() in ("A", "B", "C") and "*" not in inner:
        gene = m.group(1).upper()
        return to_two_field(gene, parse_allele_fields(f"{m.group(2)}:{m.group(3)}"))

    m2 = re.match(r"^(DRB1)(\d{2}):(\d{2,}.*)$", inner, re.I)
    if m2 and "*" not in inner:
        return to_two_field("DRB1", parse_allele_fields(f"{m2.group(2)}:{m2.group(3)}"))

    m3 = re.match(r"^(C)(\d{2}):(\d{2,}.*)$", inner, re.I)
    if m3 and "*" not in inner and not inner.upper().startswith("CW"):
        return to_two_field("C", parse_allele_fields(f"{m3.group(2)}:{m3.group(3)}"))

    return normalize_raw(raw)


def canonicalize(raw: str) -> str:
    s = raw.strip()
    try:
        return normalize_with_hla_prefix(s)
    except ValueError:
        return normalize_raw(s)


def read_csv(path):
    rows = []
    with open(path, newline="") as f:
        reader = csv.DictReader(f)
        for row in reader:
            allele = row["allele"]
            payload_key = [k for k in row if k != "allele"][0]
            rows.append((allele, float(row[payload_key])))
    return rows


def main():
    global GENE_ALIASES
    GENE_ALIASES = load_gene_aliases(ROOT / "gene_aliases.tsv")

    a_rows = read_csv(ROOT / "source_a.csv")
    b_rows = read_csv(ROOT / "source_b.csv")

    mapping = {}
    all_raw = set()
    for allele, _ in a_rows + b_rows:
        all_raw.add(allele)

    errors = []
    for raw in sorted(all_raw):
        try:
            mapping[raw] = canonicalize(raw)
        except Exception as e:
            errors.append((raw, str(e)))

    if errors:
        print("ERRORS:", errors)
        for raw in errors:
            print(f"  {raw}")

    a_by_canon = {}
    for allele, payload in a_rows:
        c = mapping[allele]
        if c in a_by_canon and a_by_canon[c][1] != payload:
            print(f"WARN A duplicate canon {c}: {a_by_canon[c]} vs {payload}")
        a_by_canon[c] = (allele, payload)

    b_by_canon = {}
    for allele, payload in b_rows:
        c = mapping[allele]
        if c in b_by_canon and b_by_canon[c][1] != payload:
            print(f"WARN B duplicate canon {c}: {b_by_canon[c]} vs {payload}")
        b_by_canon[c] = (allele, payload)

    common = set(a_by_canon) & set(b_by_canon)
    merged = []
    naive_dropped = 0
    for c in sorted(common):
        raw_a, pa = a_by_canon[c]
        raw_b, pb = b_by_canon[c]
        if raw_a != raw_b:
            naive_dropped += 1
        merged.append(
            {"canonical_allele": c, "a_payload": pa, "b_payload": pb}
        )

    only_a = set(a_by_canon) - common
    only_b = set(b_by_canon) - common
    print(f"joined: {len(common)}")
    print(f"naive_dropped: {naive_dropped}")
    print(f"only_a ({len(only_a)}):", sorted(only_a))
    print(f"only_b ({len(only_b)}):", sorted(only_b))

    out = {
        "merged": merged,
        "joined_count": len(common),
        "naive_dropped_count": naive_dropped,
        "mapping": mapping,
    }
    out_path = ROOT / "output" / "submission.json"
    with open(out_path, "w") as f:
        json.dump(out, f, indent=2)
        f.write("\n")
    print(f"Wrote {out_path}")


if __name__ == "__main__":
    main()
