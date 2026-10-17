from __future__ import annotations

import csv
import json
import random
import sys
from pathlib import Path

ENV_DIR = Path(__file__).resolve().parent
ROOT = ENV_DIR.parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from envs.allele_normalizer.normalize import canonicalize, load_aliases  # noqa: E402

SEED = 20261014

GENES = {
    "A": ["01:01", "01:04", "02:01", "02:06", "03:01", "11:01", "23:01", "24:02", "26:01", "29:02", "31:01", "32:01", "33:03", "68:01"],
    "B": ["07:02", "08:01", "13:02", "14:02", "15:01", "18:01", "27:05", "35:01", "38:01", "40:01", "44:02", "51:01", "57:01"],
    "C": ["01:02", "02:02", "03:04", "04:01", "05:01", "06:02", "07:01", "07:02", "08:02", "12:03", "14:02", "15:02", "16:01"],
    "DRB1": ["01:01", "03:01", "04:01", "07:01", "11:01", "15:01"],
    "DQB1": ["02:01", "03:01", "03:02", "05:01", "06:02"],
}

SEROLOGICAL = {
    "HLA-A*01:01": "A1",
    "HLA-A*02:01": "A2",
    "HLA-A*03:01": "A3",
    "HLA-A*11:01": "A11",
    "HLA-A*24:02": "A24",
    "HLA-B*07:02": "B7",
    "HLA-B*08:01": "B8",
    "HLA-B*27:05": "B27",
    "HLA-B*44:02": "B44",
    "HLA-B*57:01": "B57",
    "HLA-C*07:01": "Cw7",
    "HLA-C*08:02": "Cw8",
    "HLA-DRB1*04:01": "DR4",
    "HLA-DRB1*15:01": "DR15",
    "HLA-DQB1*03:02": "DQ8",
}

FORCED_SHARED = ["HLA-A*01:04", "HLA-A*02:01", "HLA-A*02:06", "HLA-B*07:02"]
FORCED_ONLY_A = {"HLA-C*07:01": "Cw7"}
FORCED_ONLY_B = {
    "HLA-C*07:02": "C*0702",
    "HLA-A*01:04N": "HLA-A*01:04N",
    "HLA-A*02:01G": "HLA-A*02:01:01:01G",
}

HETERODIMERS = [
    (
        "HLA-A*02:01-HLA-B*07:02",
        "A2-B7",
        "HLA-B*07:02-HLA-A*02:01",
    ),
    (
        "HLA-A*02:06-HLA-B*07:02",
        "A*0206-B*0702",
        "HLA-B*07:02-HLA-A*02:06",
    ),
    (
        "HLA-DQA1*01:02-HLA-DQB1*05:01",
        "HLA-DQA1*01:02-HLA-DQB1*05:01",
        "DQB1*05:01-DQA1*01:02",
    ),
]

README_TEXT = """# Merge two allele tables

Two groups typed the same panel of HLA alleles and recorded a value for each.
The two tables were produced years apart by different teams. An alias table
from the older team's records is included. It is not guaranteed to be current.

`payload_b` is stored as a percentage of the scale used for `payload_a`.
Divide it by 100 when merging.

Files:

```
source_a.csv     allele, payload_a
source_b.csv     allele, payload_b
gene_aliases.tsv canonical_symbol, prev_symbol, synonyms (pipe separated)
```
"""


def all_canonicals():
    names = []
    for gene, fields in GENES.items():
        for field in fields:
            names.append("HLA-{}*{}".format(gene, field))
    return sorted(names)


def spell(canonical: str, style: int) -> str:
    gene, fields = canonical[4:].split("*")
    field_one, field_two = fields.split(":")[:2]
    if style == 0:
        return canonical
    if style == 1:
        return "{}*{}{}".format(gene, field_one, field_two)
    if style == 2:
        return "HLA-{}{}:{}".format(gene, field_one, field_two)
    if style == 3:
        return "{}{}{}".format(gene, field_one, field_two)
    if style == 4 and canonical in SEROLOGICAL:
        return SEROLOGICAL[canonical]
    if style == 5:
        return canonical.lower()
    if style == 6:
        return canonical + " "
    if style == 7:
        return "{}:01".format(canonical)
    if style == 8:
        return "{}:01:01".format(canonical)
    if style == 9 and gene == "C":
        return "HLA-Cw*{}:{}".format(field_one, field_two)
    return "{}*{}:{}".format(gene, field_one, field_two)


STYLE_COUNT = 10
DUPLICATE_ROWS_A = 2


def build_fixture(name: str) -> dict:
    rng = random.Random(SEED)
    canonicals = [
        c
        for c in all_canonicals()
        if c not in FORCED_SHARED and c not in FORCED_ONLY_A and c not in FORCED_ONLY_B
    ]
    rng.shuffle(canonicals)

    shared = FORCED_SHARED + canonicals[: 30 - len(FORCED_SHARED)]
    rest = canonicals[30 - len(FORCED_SHARED) :]
    only_a = list(FORCED_ONLY_A) + rest[:4]
    only_b = list(FORCED_ONLY_B) + rest[4:7]

    fixture_dir = ENV_DIR / "fixtures" / name
    inputs_dir = fixture_dir / "inputs"
    key_dir = fixture_dir / "key"
    inputs_dir.mkdir(parents=True, exist_ok=True)
    key_dir.mkdir(parents=True, exist_ok=True)

    alias_rows = []
    for canonical in sorted(SEROLOGICAL):
        alias_rows.append(
            {
                "canonical_symbol": canonical,
                "prev_symbol": SEROLOGICAL[canonical],
                "synonyms": "{}|{}".format(
                    SEROLOGICAL[canonical].lower(), SEROLOGICAL[canonical]
                ),
            }
        )
    aliases_path = inputs_dir / "gene_aliases.tsv"
    with open(aliases_path, "w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=["canonical_symbol", "prev_symbol", "synonyms"],
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows(alias_rows)
        writer.writerow(
            {
                "canonical_symbol": "HLA-A*02:01",
                "prev_symbol": "A0206",
                "synonyms": "A*0206|A206|HLA-A*02:06",
            }
        )

    aliases = load_aliases(aliases_path)

    rows_a = []
    rows_b = []
    raw_to_canonical = {}

    for index, canonical in enumerate(sorted(shared)):
        style_a = index % STYLE_COUNT
        if index % 15 < 7:
            style_b = style_a
        else:
            style_b = (index + 3) % STYLE_COUNT
            if spell(canonical, style_b) == spell(canonical, style_a):
                style_b = (style_b + 1) % STYLE_COUNT
        raw_a = spell(canonical, style_a)
        raw_b = spell(canonical, style_b)
        rows_a.append((raw_a, round(rng.uniform(0.0, 100.0), 3)))
        rows_b.append((raw_b, round(rng.uniform(0.0, 100.0), 3)))
        raw_to_canonical[raw_a] = canonical
        raw_to_canonical[raw_b] = canonical

    for index, canonical in enumerate(sorted(only_a)):
        raw = FORCED_ONLY_A.get(canonical) or spell(canonical, index % 4)
        rows_a.append((raw, round(rng.uniform(0.0, 100.0), 3)))
        raw_to_canonical[raw] = canonical

    for index, canonical in enumerate(sorted(only_b)):
        raw = FORCED_ONLY_B.get(canonical) or spell(canonical, (index + 1) % 4)
        rows_b.append((raw, round(rng.uniform(0.0, 100.0), 3)))
        raw_to_canonical[raw] = canonical

    for canonical, raw_a, raw_b in HETERODIMERS:
        rows_a.append((raw_a, round(rng.uniform(0.0, 100.0), 3)))
        rows_b.append((raw_b, round(rng.uniform(0.0, 100.0), 3)))
        raw_to_canonical[raw_a] = canonical
        raw_to_canonical[raw_b] = canonical

    distinct_a = len(rows_a)
    distinct_b = len(rows_b)
    for row in rng.sample(rows_a, DUPLICATE_ROWS_A):
        rows_a.append(row)
    rng.shuffle(rows_a)
    rng.shuffle(rows_b)

    with open(inputs_dir / "source_a.csv", "w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, lineterminator="\n")
        writer.writerow(["allele", "payload_a"])
        writer.writerows(rows_a)
    with open(inputs_dir / "source_b.csv", "w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, lineterminator="\n")
        writer.writerow(["allele", "payload_b"])
        writer.writerows(rows_b)
    (inputs_dir / "README.md").write_text(README_TEXT, encoding="utf-8")

    for raw, canonical in raw_to_canonical.items():
        assert canonicalize(raw, aliases) == canonical, (raw, canonical)

    payload_a = {raw_to_canonical[raw]: value for raw, value in rows_a}
    payload_b = {raw_to_canonical[raw]: value for raw, value in rows_b}
    joined = sorted(set(payload_a) & set(payload_b))

    raw_a_strings = {raw for raw, _ in rows_a}
    raw_b_strings = {raw for raw, _ in rows_b}
    naive_joined = raw_a_strings & raw_b_strings
    naive_dropped = len(joined) - len(naive_joined)

    key = {
        "fixture": name,
        "mapping": {raw: canonical for raw, canonical in sorted(raw_to_canonical.items())},
        "expected_joined": joined,
        "joined_count": len(joined),
        "naive_dropped_count": naive_dropped,
        "naive_joined_count": len(naive_joined),
        "payload_a": {c: payload_a[c] for c in joined},
        "payload_b": {c: round(payload_b[c] / 100.0, 6) for c in joined},
        "payload_b_raw": {c: payload_b[c] for c in joined},
        "only_a_payloads": {c: payload_a[c] for c in sorted(only_a)},
        "only_b_payloads": {c: payload_b[c] for c in sorted(only_b)},
        "source_a_rows": distinct_a,
        "source_b_rows": distinct_b,
        "source_a_raw_rows": len(rows_a),
        "source_b_raw_rows": len(rows_b),
        "only_a": sorted(only_a),
        "only_b": sorted(only_b),
        "seed": SEED,
    }
    (key_dir / "meta.json").write_text(
        json.dumps(key, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )

    write_reference_submissions(name, key)
    return {
        "fixture": name,
        "joined_count": len(joined),
        "naive_joined_count": len(naive_joined),
        "naive_dropped_count": naive_dropped,
        "source_a_rows": len(rows_a),
        "source_b_rows": len(rows_b),
    }


def write_reference_submissions(name: str, key: dict) -> None:
    reference_dir = ENV_DIR / "reference" / name
    reference_dir.mkdir(parents=True, exist_ok=True)

    merged = [
        {
            "canonical_allele": canonical,
            "a_payload": key["payload_a"][canonical],
            "b_payload": key["payload_b"][canonical],
        }
        for canonical in key["expected_joined"]
    ]

    good = {
        "merged": merged,
        "joined_count": key["joined_count"],
        "naive_dropped_count": key["naive_dropped_count"],
        "mapping": key["mapping"],
    }
    (reference_dir / "good.json").write_text(
        json.dumps(good, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )

    naive_merged = [
        row
        for row in merged
        if row["canonical_allele"] in key["expected_joined"][: key["naive_joined_count"]]
    ]
    naive = {
        "merged": naive_merged,
        "joined_count": len(naive_merged),
        "naive_dropped_count": 0,
        "mapping": {},
    }
    (reference_dir / "bad_naive_join.json").write_text(
        json.dumps(naive, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )

    swapped = [
        {
            "canonical_allele": row["canonical_allele"],
            "a_payload": row["b_payload"],
            "b_payload": row["a_payload"],
        }
        for row in merged
    ]
    cross_wired = {
        "merged": swapped,
        "joined_count": key["joined_count"],
        "naive_dropped_count": key["naive_dropped_count"],
        "mapping": key["mapping"],
    }
    (reference_dir / "bad_cross_wired_payloads.json").write_text(
        json.dumps(cross_wired, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )

    overmerged = {
        "merged": merged
        + [
            {
                "canonical_allele": "HLA-C*07:01",
                "a_payload": key["only_a_payloads"]["HLA-C*07:01"],
                "b_payload": key["only_b_payloads"]["HLA-C*07:02"],
            }
        ],
        "joined_count": key["joined_count"] + 1,
        "naive_dropped_count": key["naive_dropped_count"] + 1,
        "mapping": dict(key["mapping"], **{"C*0702": "HLA-C*07:01"}),
    }
    (reference_dir / "bad_overmerged_serotype.json").write_text(
        json.dumps(overmerged, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )

    collapsed = [row for row in merged if row["canonical_allele"] != "HLA-A*02:06"]
    alias_override = {
        "merged": collapsed,
        "joined_count": len(collapsed),
        "naive_dropped_count": key["naive_dropped_count"],
        "mapping": dict(key["mapping"], **{"A*0206": "HLA-A*02:01", "A0206": "HLA-A*02:01"}),
    }
    (reference_dir / "bad_alias_overrides_parseable_name.json").write_text(
        json.dumps(alias_override, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )

    unscaled = [
        dict(row, b_payload=key["payload_b_raw"][row["canonical_allele"]])
        for row in merged
    ]
    (reference_dir / "bad_unscaled_payload_b.json").write_text(
        json.dumps(
            dict(good, merged=unscaled),
            indent=2,
            sort_keys=True,
        )
        + "\n",
        encoding="utf-8",
    )

    inflated = {
        "merged": merged,
        "joined_count": key["joined_count"] + 7,
        "naive_dropped_count": key["naive_dropped_count"] + 7,
        "mapping": key["mapping"],
    }
    (reference_dir / "bad_inflated_counts.json").write_text(
        json.dumps(inflated, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )


def main() -> None:
    print(json.dumps(build_fixture("canonical"), indent=2))


if __name__ == "__main__":
    main()
