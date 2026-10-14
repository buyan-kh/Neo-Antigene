import csv
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "analysis"))
from stats import two_sided_p, adjust_fdr, FDR_Q

ROOT = os.path.join(os.path.dirname(__file__), "..")
RAW = os.path.join(ROOT, "data", "all_measured.csv")


def load(path):
    features = {}
    with open(path, newline="") as handle:
        for row in csv.DictReader(handle):
            feature_id = row["feature_id"]
            features.setdefault(feature_id, {"A": [], "B": []})
            features[feature_id][row["group"]].append(float(row["value"]))
    return features


def main():
    features = load(RAW)
    p_by_feature = {
        fid: two_sided_p(g["A"], g["B"]) for fid, g in features.items()
    }
    adjusted = adjust_fdr(p_by_feature)
    significant = sorted(
        fid for fid, q in adjusted.items() if q <= FDR_Q
    )
    print("m=", len(p_by_feature))
    print("n_sig=", len(significant))
    for fid in significant:
        print(fid)


if __name__ == "__main__":
    main()
