#!/usr/bin/env python3
import csv
import json
import re
from pathlib import Path

ROOT = Path(__file__).resolve().parent
TRAIN = ROOT / "predictor_training_set"


def normalize_allele(raw: str) -> str:
    """Canonical form: HLA-A*24:02 style."""
    s = raw.strip()
    if not s:
        return s
    s = s.replace(" ", "")
    # HLA-A*02:01 or HLA-A*02:01:01:01 (extended typing)
    m = re.match(r"^HLA-([ABC])\*(\d{2}):(\d{2})", s, re.I)
    if m:
        return f"HLA-{m.group(1).upper()}*{m.group(2)}:{m.group(3)}"
    # HLA-A11:01 or HLA-A02:01 (missing star before digits)
    m = re.match(r"^HLA-([ABC])(\d{2}):(\d{2})", s, re.I)
    if m:
        return f"HLA-{m.group(1).upper()}*{m.group(2)}:{m.group(3)}"
    # hla-a*01:01 without HLA- prefix
    m = re.match(r"^([ABC])\*(\d{2}):(\d{2})", s, re.I)
    if m:
        return f"HLA-{m.group(1).upper()}*{m.group(2)}:{m.group(3)}"
    # A*2402 or A2402 from predictions / allele_ids
    m = re.match(r"^([ABC])\*?(\d{2})(\d{2})$", s, re.I)
    if m:
        return f"HLA-{m.group(1).upper()}*{m.group(2)}:{m.group(3)}"
    return s


def extract_peptide_panel_b(sequence: str) -> str:
    """Uppercase run in sequence is the measured peptide."""
    upper = "".join(c for c in sequence if c.isupper())
    return upper


def load_allele_map() -> dict[str, str]:
    mapping = {}
    with open(TRAIN / "allele_ids.tsv", newline="") as f:
        for row in csv.DictReader(f, delimiter="\t"):
            aid = row["allele_id"].strip()
            mapping[aid] = normalize_allele(row["allele"].strip())
    return mapping


def build_training_pairs() -> tuple[set[tuple[str, str]], dict[tuple[str, str], list[str]]]:
    """Returns set of (peptide, allele) and provenance for each pair."""
    pairs: set[tuple[str, str]] = set()
    provenance: dict[tuple[str, str], list[str]] = {}

    def add(peptide: str, allele: str, source: str):
        pep = re.sub(r"\s+", "", peptide).upper()
        alle = normalize_allele(allele)
        key = (pep, alle)
        pairs.add(key)
        provenance.setdefault(key, []).append(source)

    with open(TRAIN / "panel_a.csv", newline="") as f:
        for row in csv.DictReader(f):
            add(row["peptide"], row["mhc_allele"], "panel_a.csv")

    with open(TRAIN / "panel_b.tsv", newline="") as f:
        for row in csv.DictReader(f, delimiter="\t"):
            pep = extract_peptide_panel_b(row["sequence"])
            add(pep, row["hla"], "panel_b.tsv")

    allele_map = load_allele_map()
    with open(TRAIN / "panel_c.csv", newline="") as f:
        for row in csv.DictReader(f):
            aid = row["allele_id"].strip()
            alle = allele_map.get(aid, aid)
            add(row["peptide"], alle, "panel_c.csv")

    return pairs, provenance


def load_benchmark() -> list[dict]:
    rows = []
    with open(ROOT / "benchmark.csv", newline="") as f:
        for row in csv.DictReader(f):
            rows.append(
                {
                    "peptide": row["peptide"].strip().upper(),
                    "allele": normalize_allele(row["mhc_allele"]),
                    "label": 1 if row["true_label"].strip() == "binder" else 0,
                    "true_label": row["true_label"].strip(),
                }
            )
    return rows


def load_predictions() -> dict[tuple[str, str], float]:
    pred = {}
    with open(ROOT / "predictions.csv", newline="") as f:
        for row in csv.DictReader(f):
            pep = row["peptide"].strip().upper()
            alle = normalize_allele(row["mhc_allele"])
            pred[(pep, alle)] = float(row["pred_score"])
    return pred


def auroc(y_true, y_score):
    if len(set(y_true)) < 2:
        return None
    # Mann-Whitney U / rank-based AUROC
    pairs = sorted(zip(y_score, y_true), key=lambda x: x[0])
    n_pos = sum(y_true)
    n_neg = len(y_true) - n_pos
    rank = 1
    sum_ranks_pos = 0.0
    i = 0
    while i < len(pairs):
        j = i
        while j < len(pairs) and pairs[j][0] == pairs[i][0]:
            j += 1
        avg_rank = (rank + rank + (j - i) - 1) / 2.0
        for k in range(i, j):
            if pairs[k][1] == 1:
                sum_ranks_pos += avg_rank
        rank += j - i
        i = j
    u = sum_ranks_pos - n_pos * (n_pos + 1) / 2.0
    return u / (n_pos * n_neg)


def main():
    train_pairs, provenance = build_training_pairs()
    benchmark = load_benchmark()
    predictions = load_predictions()

    # Duplicate benchmark rows
    seen_bench: dict[tuple[str, str], int] = {}
    for r in benchmark:
        k = (r["peptide"], r["allele"])
        seen_bench[k] = seen_bench.get(k, 0) + 1
    dup_bench = {k: c for k, c in seen_bench.items() if c > 1}

    leaked = []
    leak_sources: dict[tuple[str, str], list[str]] = {}
    for r in benchmark:
        k = (r["peptide"], r["allele"])
        if k in train_pairs:
            leaked.append(k)
            leak_sources[k] = list(dict.fromkeys(provenance.get(k, ["unknown"])))

    leaked_set = set(leaked)
    leaked_sorted = sorted(leaked_set)

    clean = [r for r in benchmark if (r["peptide"], r["allele"]) not in train_pairs]

    y_true = []
    y_score = []
    for r in clean:
        k = (r["peptide"], r["allele"])
        if k not in predictions:
            raise SystemExit(f"missing prediction for {k}")
        y_true.append(r["label"])
        y_score.append(predictions[k])

    auc = auroc(y_true, y_score)
    uncomputable = auc is None

    # Full benchmark AUROC for reference
    all_y = [r["label"] for r in benchmark]
    all_s = [predictions[(r["peptide"], r["allele"])] for r in benchmark]
    full_auc = auroc(all_y, all_s)

    print(f"Training pairs: {len(train_pairs)}")
    print(f"Benchmark rows: {len(benchmark)}")
    print(f"Unique benchmark pairs: {len(seen_bench)}")
    print(f"Duplicate benchmark pairs: {len(dup_bench)}")
    print(f"Leaked pairs (in training): {len(leaked_set)}")
    print(f"Clean pairs: {len(clean)}")
    print(f"Full benchmark AUROC: {full_auc:.6f}" if full_auc else "Full AUROC uncomputable")
    print(f"Clean AUROC: {auc:.6f}" if auc is not None else "Clean AUROC uncomputable")

    # Check overlap with held_out only
    held_only = set()
    with open(TRAIN / "held_out_from_fit.csv", newline="") as f:
        for row in csv.DictReader(f):
            held_only.add(
                (row["peptide"].strip().upper(), normalize_allele(row["mhc_allele"]))
            )
    leak_held = leaked_set & held_only
    print(f"Leaked that are exactly held_out_from_fit list: {len(leak_held)}")

    # Sample provenance
    for k in leaked_sorted[:5]:
        print(f"  {k} -> {leak_sources[k]}")

    mechanism = "none"
    if dup_bench:
        mechanism = "duplicate_benchmark_rows"
    elif leaked_set:
        mechanism = "train_test_overlap"

    out_dir = ROOT / "output"
    out_dir.mkdir(exist_ok=True)

    with open(out_dir / "leaked_pairs.csv", "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["peptide", "allele"])
        for p, a in leaked_sorted:
            w.writerow([p, a])

    out = {
        "contaminated": bool(leaked_set or dup_bench),
        "leaked_pair_count": len(leaked_set),
        "mechanism": mechanism,
        "trusted_metric": not bool(leaked_set or dup_bench),
        "corrected_metric": {
            "name": "auroc",
            "value": auc if auc is not None else None,
            "n_clean_pairs": len(clean),
            "uncomputable": uncomputable,
        },
    }

    with open(out_dir / "submission.json", "w") as f:
        json.dump(out, f, indent=2)
        f.write("\n")

    # Debug file
    with open(out_dir / "debug_leakage.txt", "w") as f:
        f.write(f"mechanism={mechanism}\n")
        for k in leaked_sorted:
            f.write(f"{k[0]},{k[1]} -> {leak_sources[k]}\n")


if __name__ == "__main__":
    main()
