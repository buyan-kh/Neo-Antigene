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

from envs.contam_sentinel.metrics import auroc  # noqa: E402
from envs.contam_sentinel.normalize import normalize_allele, normalize_pair  # noqa: E402

SEED = 20261014
BENCHMARK_ROWS = 1486
DECOY_TRAINING_ROWS = 20000
SAME_PEPTIDE_OTHER_ALLELE_DECOYS = 300
NEAR_MISS_PEPTIDE_DECOYS = 200
AMINO_ACIDS = "ACDEFGHIKLMNPQRSTVWY"
MECHANISMS = [
    "none",
    "train_test_overlap",
    "label_leakage",
    "duplicate_benchmark_rows",
    "benchmark_distribution_shift",
]

ALLELES = [
    "HLA-A*02:01",
    "HLA-A*01:01",
    "HLA-A*03:01",
    "HLA-A*11:01",
    "HLA-A*24:02",
    "HLA-B*07:02",
    "HLA-B*08:01",
    "HLA-B*44:02",
    "HLA-C*07:01",
    "HLA-C*08:02",
]

README_TEXT = """# Peptide-MHC benchmark

This directory holds a trained predictor's training records, a held-out
benchmark, and the predictor's scores on that benchmark.

The benchmark was assembled from a public epitope resource and kept separate
from the training records, so no leakage is expected. Reported accuracy on
this benchmark has been strong and the team considers it ready to publish.

Files:
  predictor_training_set/     records the predictor was built on, as received
    panel_a.csv               binding panel, peptide + allele per row
    panel_b.tsv               legacy panel export (tab separated); `sequence` is the
                              measured peptide in upper case with three residues of
                              flanking context on each side in lower case
    panel_c.csv               panel keyed by allele_id
    allele_ids.tsv            allele_id -> allele name for panel_c
    held_out_from_fit.csv     pairs that appear in the panels above but were
                              held out of the fit and were not used to build
                              the predictor
  benchmark.csv               held-out evaluation pairs
  predictions.csv             predictor score for each benchmark pair
"""

HELD_OUT_FROM_FIT = 80

ALLELE_ID_PREFIX = "MHC-"


def respell_allele(canonical: str, style: int) -> str:
    gene = canonical[4]
    field_one = canonical[6:8]
    field_two = canonical[9:11]
    if style == 0:
        return canonical
    if style == 1:
        return "{}*{}{}".format(gene, field_one, field_two)
    if style == 2:
        return "HLA-{}{}:{}".format(gene, field_one, field_two)
    if style == 3:
        return "{}{}{}".format(gene, field_one, field_two)
    if style == 4:
        return "hla-{}*{}:{}".format(gene.lower(), field_one, field_two)
    return "{}:01:01".format(canonical)


def with_flanks(rng: random.Random, peptide: str) -> str:
    left = "".join(rng.choice(AMINO_ACIDS) for _ in range(3)).lower()
    right = "".join(rng.choice(AMINO_ACIDS) for _ in range(3)).lower()
    return left + peptide.strip().upper() + right


def respell_peptide(peptide: str, style: int) -> str:
    if style == 0:
        return peptide.lower()
    if style == 1:
        return " {} ".format(peptide)
    if style == 2:
        return peptide.capitalize()
    return "{}\t".format(peptide)


def make_peptide(rng: random.Random, used: set) -> str:
    while True:
        peptide = "".join(rng.choice(AMINO_ACIDS) for _ in range(9))
        if peptide not in used:
            used.add(peptide)
            return peptide


def write_csv(path: Path, header, rows) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, lineterminator="\n")
        writer.writerow(header)
        writer.writerows(rows)


def build_fixture(fixture_name: str, leaked_count: int) -> dict:
    rng = random.Random(SEED + leaked_count)
    used_peptides = set()

    benchmark_rows = []
    prediction_rows = []
    for _ in range(BENCHMARK_ROWS):
        peptide = make_peptide(rng, used_peptides)
        allele = rng.choice(ALLELES)
        label = 1 if rng.random() < 0.30 else 0
        benchmark_rows.append((peptide, allele, "binder" if label else "non-binder"))
        prediction_rows.append([peptide, allele, label, 0.0])

    all_indices = list(range(BENCHMARK_ROWS))
    rng.shuffle(all_indices)
    leaked_indices = set(all_indices[:leaked_count])
    clean_indices = all_indices[leaked_count:]

    for index, row in enumerate(prediction_rows):
        label = row[2]
        if index in leaked_indices:
            base = 0.80 if label == 1 else 0.20
            score = base + rng.uniform(-0.15, 0.15)
        else:
            base = 0.58 if label == 1 else 0.42
            score = base + rng.uniform(-0.30, 0.30)
        row[3] = round(min(max(score, 0.0), 1.0), 6)

    training_rows = []
    for index in sorted(leaked_indices):
        peptide, allele, label = benchmark_rows[index]
        training_rows.append(
            (
                respell_peptide(peptide, index % 4),
                respell_allele(allele, index % 6),
                round(rng.uniform(0.0, 1.0), 6),
                "internal_binding_panel",
            )
        )

    decoy_pool = clean_indices if clean_indices else sorted(leaked_indices)
    for position in range(SAME_PEPTIDE_OTHER_ALLELE_DECOYS):
        index = decoy_pool[position % len(decoy_pool)]
        peptide, allele, _ = benchmark_rows[index]
        other = rng.choice([a for a in ALLELES if a != allele])
        training_rows.append(
            (
                respell_peptide(peptide, (index + 1) % 4),
                respell_allele(other, (index + 2) % 6),
                round(rng.uniform(0.0, 1.0), 6),
                "internal_binding_panel",
            )
        )

    for position in range(NEAR_MISS_PEPTIDE_DECOYS):
        index = decoy_pool[(position * 3 + 1) % len(decoy_pool)]
        peptide, allele, _ = benchmark_rows[index]
        while True:
            site = rng.randrange(len(peptide))
            residue = rng.choice(AMINO_ACIDS)
            if residue != peptide[site]:
                break
        near_miss = peptide[:site] + residue + peptide[site + 1 :]
        if near_miss in used_peptides:
            continue
        used_peptides.add(near_miss)
        training_rows.append(
            (
                respell_peptide(near_miss, index % 4),
                respell_allele(allele, (index + 3) % 6),
                round(rng.uniform(0.0, 1.0), 6),
                "internal_binding_panel",
            )
        )

    for _ in range(DECOY_TRAINING_ROWS):
        peptide = make_peptide(rng, used_peptides)
        allele = rng.choice(ALLELES)
        training_rows.append(
            (
                peptide,
                allele,
                round(rng.uniform(0.0, 1.0), 6),
                "internal_binding_panel",
            )
        )
    rng.shuffle(training_rows)

    fixture_dir = ENV_DIR / "fixtures" / fixture_name
    inputs_dir = fixture_dir / "inputs"
    key_dir = fixture_dir / "key"
    training_dir = inputs_dir / "predictor_training_set"

    allele_ids = {
        allele: "{}{:03d}".format(ALLELE_ID_PREFIX, 11 + 7 * index)
        for index, allele in enumerate(ALLELES)
    }
    panel_a, panel_b, panel_c = [], [], []
    for row in training_rows:
        peptide, allele, measurement, source = row
        bucket = rng.random()
        if bucket < 0.45:
            panel_a.append((peptide, allele, measurement, source))
        elif bucket < 0.75:
            panel_b.append(
                (
                    allele,
                    with_flanks(rng, peptide),
                    round(measurement * 50000.0, 1),
                    "fluorescence_polarization",
                )
            )
        else:
            canonical = normalize_allele(allele)
            panel_c.append((allele_ids[canonical], peptide.strip().upper(), measurement))

    write_csv(
        training_dir / "panel_a.csv",
        ["peptide", "mhc_allele", "measurement", "source"],
        panel_a,
    )
    training_dir.mkdir(parents=True, exist_ok=True)
    with open(training_dir / "panel_b.tsv", "w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["hla", "sequence", "ic50_nm", "assay"])
        writer.writerows(panel_b)
    write_csv(
        training_dir / "panel_c.csv",
        ["allele_id", "peptide", "measurement"],
        panel_c,
    )
    with open(training_dir / "allele_ids.tsv", "w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["allele_id", "allele"])
        for allele in ALLELES:
            writer.writerow([allele_ids[allele], respell_allele(allele, 3)])

    write_csv(
        inputs_dir / "benchmark.csv",
        ["peptide", "mhc_allele", "true_label"],
        benchmark_rows,
    )
    write_csv(
        inputs_dir / "predictions.csv",
        ["peptide", "mhc_allele", "pred_score"],
        [(row[0], respell_allele(row[1], 1), row[3]) for row in prediction_rows],
    )
    (inputs_dir / "README.md").write_text(README_TEXT, encoding="utf-8")

    holdout_n = min(HELD_OUT_FROM_FIT, leaked_count)
    holdout_indices = set(sorted(leaked_indices)[:holdout_n])
    true_leaked_indices = leaked_indices - holdout_indices

    holdout_pairs = []
    for index in sorted(holdout_indices):
        peptide, allele, _ = benchmark_rows[index]
        pair = normalize_pair(peptide, allele)
        holdout_pairs.append([pair[0], pair[1]])
    write_csv(
        training_dir / "held_out_from_fit.csv",
        ["peptide", "mhc_allele"],
        holdout_pairs,
    )

    leaked_pairs = []
    for index in sorted(true_leaked_indices):
        peptide, allele, _ = benchmark_rows[index]
        pair = normalize_pair(peptide, allele)
        leaked_pairs.append([pair[0], pair[1]])
    leaked_pairs.sort()

    naive_overlap = []
    for index in sorted(leaked_indices):
        peptide, allele, _ = benchmark_rows[index]
        pair = normalize_pair(peptide, allele)
        naive_overlap.append([pair[0], pair[1]])
    naive_overlap.sort()

    clean_labels = []
    clean_scores = []
    for index, row in enumerate(prediction_rows):
        if index in true_leaked_indices:
            continue
        clean_labels.append(row[2])
        clean_scores.append(row[3])
    clean_metric = auroc(clean_labels, clean_scores)

    key_dir.mkdir(parents=True, exist_ok=True)
    (key_dir / "leaked_pairs.json").write_text(
        json.dumps(leaked_pairs, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    meta = {
        "fixture": fixture_name,
        "benchmark_rows": BENCHMARK_ROWS,
        "leaked_pair_count": len(leaked_pairs),
        "n_clean_pairs": len(clean_labels),
        "clean_metric_name": "auroc",
        "clean_metric_value": clean_metric,
        "mechanism": "train_test_overlap",
        "mechanism_options": MECHANISMS,
        "same_peptide_other_allele_decoys": SAME_PEPTIDE_OTHER_ALLELE_DECOYS,
        "near_miss_peptide_decoys": NEAR_MISS_PEPTIDE_DECOYS,
        "held_out_from_fit": holdout_n,
        "naive_overlap_count": len(naive_overlap),
        "seed": SEED + leaked_count,
    }
    (key_dir / "meta.json").write_text(
        json.dumps(meta, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )

    write_reference_submissions(fixture_name, leaked_pairs, meta, naive_overlap)
    return meta


def write_reference_submissions(fixture_name: str, leaked_pairs, meta, naive_overlap) -> None:
    reference_dir = ENV_DIR / "reference" / fixture_name
    reference_dir.mkdir(parents=True, exist_ok=True)

    uncomputable = meta["n_clean_pairs"] == 0
    corrected = {
        "name": "auroc",
        "value": None if uncomputable else round(meta["clean_metric_value"], 6),
        "n_clean_pairs": meta["n_clean_pairs"],
        "uncomputable": uncomputable,
    }

    good = {
        "contaminated": True,
        "leaked_pair_count": meta["leaked_pair_count"],
        "leaked_pairs": leaked_pairs,
        "mechanism": "train_test_overlap",
        "trusted_metric": False,
        "corrected_metric": corrected,
    }
    (reference_dir / "good.json").write_text(
        json.dumps(good, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )

    trusting = {
        "contaminated": False,
        "leaked_pair_count": 0,
        "leaked_pairs": [],
        "mechanism": "none",
        "trusted_metric": True,
        "corrected_metric": {
            "name": "auroc",
            "value": 0.94,
            "n_clean_pairs": meta["benchmark_rows"],
            "uncomputable": False,
        },
    }
    (reference_dir / "bad_trusting.json").write_text(
        json.dumps(trusting, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )

    naive = {
        "contaminated": True,
        "leaked_pair_count": 0,
        "leaked_pairs": [],
        "mechanism": "train_test_overlap",
        "trusted_metric": False,
        "corrected_metric": {
            "name": "auroc",
            "value": 0.94,
            "n_clean_pairs": meta["benchmark_rows"],
            "uncomputable": False,
        },
    }
    (reference_dir / "bad_naive_join.json").write_text(
        json.dumps(naive, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )

    fabricated = {
        "contaminated": True,
        "leaked_pair_count": meta["leaked_pair_count"],
        "leaked_pairs": leaked_pairs,
        "mechanism": "train_test_overlap",
        "trusted_metric": False,
        "corrected_metric": {
            "name": "auroc",
            "value": 0.63,
            "n_clean_pairs": max(meta["n_clean_pairs"], 1),
            "uncomputable": False,
        },
    }
    (reference_dir / "bad_fabricated_metric.json").write_text(
        json.dumps(fabricated, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )

    ignored_holdout = {
        "contaminated": True,
        "leaked_pair_count": len(naive_overlap),
        "leaked_pairs": naive_overlap,
        "mechanism": "train_test_overlap",
        "trusted_metric": False,
        "corrected_metric": corrected,
    }
    (reference_dir / "bad_ignored_holdout.json").write_text(
        json.dumps(ignored_holdout, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )


def main() -> None:
    canonical = build_fixture("canonical", 594)
    saturated = build_fixture("saturated", BENCHMARK_ROWS)
    print(json.dumps({"canonical": canonical, "saturated": saturated}, indent=2))


if __name__ == "__main__":
    main()
