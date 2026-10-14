"""Select the analysis table from the raw measurements.

The raw panel is large. Features that do not show a raw difference are
dropped here so the FDR step is not dominated by noise. The analysis
driver reads only the table this script writes.
"""

import csv
import os
import sys

from stats import ALPHA, two_sided_p

HERE = os.path.dirname(os.path.abspath(__file__))
RAW = os.path.join(HERE, os.pardir, "data", "all_measured.csv")
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
