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

from envs.pvalue_police.stats_core import (  # noqa: E402
    KNOWN_SIGMA,
    benjamini_hochberg,
    bonferroni,
    p_values,
    uncorrected,
)

SEED = 20261023
FEATURE_COUNT = 400
REPLICATES = 8
TRUE_EFFECT_COUNT = 12
MISSING_REPLICATE_FEATURES = 6
SPIKE_IN_COUNT = 80
SPIKE_IN_EFFECTS = 8
EFFECT_SIZE = 1.4
ALPHA = 0.05
Q_VALUE = 0.05

BUG_CLASS = "correction_applied_after_prefiltering"
BUG_CLASSES = [
    "none",
    "no_multiple_testing_correction",
    "correction_applied_after_prefiltering",
    "wrong_test_statistic",
    "one_sided_p_value_misreported",
    "groups_mislabelled",
    "boundary_comparison_error",
]

STATS_SOURCE = '''"""Per-feature screening statistics."""

from math import erfc, sqrt

KNOWN_SIGMA = 1.0
ALPHA = 0.05
FDR_Q = 0.05


def two_sided_p(values_a, values_b):
    """Two-sided z-test p-value for a difference in means."""
    n_a, n_b = len(values_a), len(values_b)
    mean_a = sum(values_a) / n_a
    mean_b = sum(values_b) / n_b
    standard_error = KNOWN_SIGMA * sqrt(1.0 / n_a + 1.0 / n_b)
    z = (mean_a - mean_b) / standard_error
    return erfc(abs(z) / sqrt(2.0))


def adjust_fdr(p_by_feature):
    """Benjamini-Hochberg adjusted p-values."""
    ordered = sorted(p_by_feature.items(), key=lambda item: (item[1], item[0]))
    m = len(ordered)
    adjusted = {}
    running = 1.0
    for rank in range(m, 0, -1):
        feature_id, p_value = ordered[rank - 1]
        running = min(running, p_value * m / rank)
        adjusted[feature_id] = running
    return adjusted


def screen(features):
    """Score every feature, FDR-adjust, and return the ones we call significant."""
    results = {}
    for feature_id, groups in features.items():
        results[feature_id] = two_sided_p(groups["A"], groups["B"])

    adjusted = adjust_fdr(results)
    significant = [
        feature_id for feature_id, q_value in adjusted.items() if q_value <= FDR_Q
    ]

    return results, sorted(significant)
'''

RUN_SOURCE = '''"""Driver: load the screen data and print the significant features."""

import csv
import os

from stats import screen

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, os.pardir, "data", "features.csv")


def load():
    features = {}
    with open(DATA, newline="") as handle:
        for row in csv.DictReader(handle):
            feature_id = row["feature_id"]
            features.setdefault(feature_id, {"A": [], "B": []})
            features[feature_id][row["group"]].append(float(row["value"]))
    return features


def main():
    features = load()
    _, significant = screen(features)
    print("{} features screened".format(len(features)))
    print("{} significant".format(len(significant)))
    for feature_id in significant:
        print(feature_id)


if __name__ == "__main__":
    main()
'''

PREPARE_SOURCE = '''"""Select the analysis table from the raw measurements.

The raw panel is large. Features that do not show a raw difference are
dropped here so the FDR step is not dominated by noise. The analysis
driver reads only the table this script writes.
"""

import csv
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, os.pardir, "analysis"))
from stats import ALPHA, two_sided_p

RAW = os.path.join(HERE, os.pardir, "raw", "run_2024-11-02_plate_reader.csv")
OUT = os.path.join(HERE, os.pardir, "data", "features.csv")


def load(path):
    features = {}
    rows = []
    with open(path, newline="") as handle:
        for row in csv.DictReader(handle):
            rows.append(row)
            feature_id = row["feature_id"]
            features.setdefault(feature_id, {"A": [], "B": []})
            features[feature_id][row["group"]].append(float(row["value"]))
    return features, rows


def main():
    features, rows = load(RAW)
    kept = {
        feature_id
        for feature_id, groups in features.items()
        if two_sided_p(groups["A"], groups["B"]) < ALPHA
    }
    with open(OUT, "w", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=["feature_id", "group", "replicate", "value"]
        )
        writer.writeheader()
        for row in rows:
            if row["feature_id"] in kept:
                writer.writerow(row)
    print("{} features measured".format(len(features)))
    print("{} forwarded to the screen".format(len(kept)))


if __name__ == "__main__":
    main()
'''

README_TEXT = """# Feature screen

A screening analysis over a panel of candidate features. Each feature was
measured in two groups with {replicates} replicates per group (a few wells
failed QC and were dropped, so some features have one fewer). The measurement
noise is well characterised, so a z-test against a known standard deviation of
{sigma} is used. Significance calls are false-discovery-rate controlled at 5%
using Benjamini-Hochberg.

Layout:

```
bin/select_candidates.py  builds the analysis table from the raw export
analysis/stats.py         scoring and the significance call
analysis/run.py           driver that prints the significant features
raw/run_2024-11-02_plate_reader.csv   the instrument export
raw/sample_sheet.csv                  feature role (screen vs spike-in)
data/features.csv                     the table the screen is run on
results/reported_significant.csv      the current output
```

Rebuild the analysis table with `python bin/select_candidates.py`, then run
`python analysis/run.py`.
""".format(replicates=REPLICATES, sigma=KNOWN_SIGMA)


def locate_bug_span(source: str) -> dict:
    lines = source.splitlines()
    start = None
    end = None
    for index, line in enumerate(lines, start=1):
        if "kept = {" in line:
            start = index
        if "print(\"{} forwarded to the screen\".format(len(kept)))" in line:
            end = index
    return {"file": "bin/select_candidates.py", "line_start": start, "line_end": end}


def build_fixture(name: str) -> dict:
    rng = random.Random(SEED)
    fixture_dir = ENV_DIR / "fixtures" / name
    inputs_dir = fixture_dir / "inputs"
    key_dir = fixture_dir / "key"
    (inputs_dir / "analysis").mkdir(parents=True, exist_ok=True)
    (inputs_dir / "bin").mkdir(parents=True, exist_ok=True)
    (inputs_dir / "data").mkdir(parents=True, exist_ok=True)
    (inputs_dir / "raw").mkdir(parents=True, exist_ok=True)
    (inputs_dir / "results").mkdir(parents=True, exist_ok=True)
    key_dir.mkdir(parents=True, exist_ok=True)

    feature_ids = ["feat_{:04d}".format(i) for i in range(FEATURE_COUNT)]
    true_effects = set(rng.sample(feature_ids, TRUE_EFFECT_COUNT))
    short_features = set(rng.sample(feature_ids, MISSING_REPLICATE_FEATURES))
    spike_ids = ["feat_{:04d}".format(FEATURE_COUNT + i) for i in range(SPIKE_IN_COUNT)]
    spike_effects = set(rng.sample(spike_ids, SPIKE_IN_EFFECTS))

    rows = []
    for feature_id in feature_ids + spike_ids:
        is_spike = feature_id in set(spike_ids)
        if is_spike:
            shift = EFFECT_SIZE if feature_id in spike_effects else 0.0
        else:
            shift = EFFECT_SIZE if feature_id in true_effects else 0.0
        role = "spike_in" if is_spike else "screen"
        replicates_a = REPLICATES - 1 if feature_id in short_features else REPLICATES
        for replicate in range(replicates_a):
            rows.append(
                (feature_id, "A", replicate, round(rng.gauss(shift, KNOWN_SIGMA), 6), role)
            )
        for replicate in range(REPLICATES):
            rows.append(
                (feature_id, "B", replicate, round(rng.gauss(0.0, KNOWN_SIGMA), 6), role)
            )

    raw_path = inputs_dir / "raw" / "run_2024-11-02_plate_reader.csv"
    features_path = inputs_dir / "data" / "features.csv"
    with open(raw_path, "w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, lineterminator="\n")
        writer.writerow(["feature_id", "group", "replicate", "value"])
        writer.writerows([row[:4] for row in rows])
    with open(inputs_dir / "raw" / "sample_sheet.csv", "w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, lineterminator="\n")
        writer.writerow(["feature_id", "role", "plate"])
        seen = []
        for feature_id, _, _, _, role in rows:
            if feature_id in seen:
                continue
            seen.append(feature_id)
            writer.writerow([feature_id, role, "P1" if role == "screen" else "CTRL"])

    (inputs_dir / "analysis" / "stats.py").write_text(STATS_SOURCE, encoding="utf-8")
    (inputs_dir / "analysis" / "run.py").write_text(RUN_SOURCE, encoding="utf-8")
    (inputs_dir / "bin" / "select_candidates.py").write_text(PREPARE_SOURCE, encoding="utf-8")
    (inputs_dir / "README.md").write_text(README_TEXT, encoding="utf-8")
    (inputs_dir / "requirements.txt").write_text("", encoding="utf-8")

    observed_p = p_values(raw_path)
    screen_p = {f: p for f, p in observed_p.items() if f in set(feature_ids)}
    prefiltered = {f: p for f, p in observed_p.items() if p < ALPHA}
    reported = sorted(benjamini_hochberg(prefiltered, Q_VALUE))
    bh_set = sorted(benjamini_hochberg(screen_p, Q_VALUE))
    bonferroni_set = sorted(bonferroni(screen_p, ALPHA))
    bh_all = sorted(benjamini_hochberg(observed_p, Q_VALUE))
    assert bh_set != reported and set(bh_set) != set(bh_all)

    kept_ids = set(prefiltered)
    with open(features_path, "w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, lineterminator="\n")
        writer.writerow(["feature_id", "group", "replicate", "value"])
        for row in rows:
            if row[0] in kept_ids:
                writer.writerow(row[:4])

    with open(
        inputs_dir / "results" / "reported_significant.csv",
        "w",
        newline="",
        encoding="utf-8",
    ) as handle:
        writer = csv.writer(handle, lineterminator="\n")
        writer.writerow(["feature_id"])
        writer.writerows([[feature_id] for feature_id in reported])

    key = {
        "fixture": name,
        "alpha": ALPHA,
        "q_value": Q_VALUE,
        "bug_class": BUG_CLASS,
        "bug_class_options": BUG_CLASSES,
        "bug_span": locate_bug_span(PREPARE_SOURCE),
        "reported_significant": reported,
        "benjamini_hochberg": bh_set,
        "bonferroni": bonferroni_set,
        "true_effect_features": sorted(true_effects),
        "seed": SEED,
    }
    (key_dir / "meta.json").write_text(
        json.dumps(key, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )

    write_reference_submissions(name, key, bh_including_spikes=bh_all)
    return {
        "fixture": name,
        "reported": len(reported),
        "benjamini_hochberg": len(bh_set),
        "bonferroni": len(bonferroni_set),
        "true_effects": len(true_effects),
        "bug_span": key["bug_span"],
    }


def write_reference_submissions(name: str, key: dict, bh_including_spikes: list) -> None:
    reference_dir = ENV_DIR / "reference" / name
    reference_dir.mkdir(parents=True, exist_ok=True)
    span = key["bug_span"]

    good = {
        "bug_found": True,
        "bug_class": BUG_CLASS,
        "location": span,
        "corrected_method": "benjamini_hochberg",
        "corrected_significant": key["benjamini_hochberg"],
    }
    (reference_dir / "good.json").write_text(
        json.dumps(good, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )

    good_bonferroni = {
        "bug_found": True,
        "bug_class": BUG_CLASS,
        "location": span,
        "corrected_method": "bonferroni",
        "corrected_significant": key["bonferroni"],
    }
    (reference_dir / "good_bonferroni.json").write_text(
        json.dumps(good_bonferroni, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )

    echo = {
        "bug_found": False,
        "bug_class": "none",
        "location": {"file": "analysis/stats.py", "line_start": 1, "line_end": 1},
        "corrected_method": "benjamini_hochberg",
        "corrected_significant": key["reported_significant"],
    }
    (reference_dir / "bad_no_bug.json").write_text(
        json.dumps(echo, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )

    echoed_output = {
        "bug_found": True,
        "bug_class": BUG_CLASS,
        "location": span,
        "corrected_method": "benjamini_hochberg",
        "corrected_significant": key["reported_significant"],
    }
    (reference_dir / "bad_echoed_output.json").write_text(
        json.dumps(echoed_output, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )

    wrong_place = {
        "bug_found": True,
        "bug_class": BUG_CLASS,
        "location": {"file": "analysis/stats.py", "line_start": 35, "line_end": 46},
        "corrected_method": "benjamini_hochberg",
        "corrected_significant": key["benjamini_hochberg"],
    }
    (reference_dir / "bad_wrong_location.json").write_text(
        json.dumps(wrong_place, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )

    trusted_adjustment = {
        "bug_found": True,
        "bug_class": "no_multiple_testing_correction",
        "location": span,
        "corrected_method": "benjamini_hochberg",
        "corrected_significant": key["reported_significant"],
    }
    (reference_dir / "bad_trusted_team_adjustment.json").write_text(
        json.dumps(trusted_adjustment, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )

    included_spikes = {
        "bug_found": True,
        "bug_class": BUG_CLASS,
        "location": span,
        "corrected_method": "benjamini_hochberg",
        "corrected_significant": bh_including_spikes,
    }
    (reference_dir / "bad_bh_including_spike_ins.json").write_text(
        json.dumps(included_spikes, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )

    included_spikes = {
        "bug_found": True,
        "bug_class": BUG_CLASS,
        "location": span,
        "corrected_method": "benjamini_hochberg",
        "corrected_significant": bh_including_spikes,
    }
    (reference_dir / "bad_included_spike_ins.json").write_text(
        json.dumps(included_spikes, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )


def main() -> None:
    print(json.dumps(build_fixture("canonical"), indent=2))


if __name__ == "__main__":
    main()
