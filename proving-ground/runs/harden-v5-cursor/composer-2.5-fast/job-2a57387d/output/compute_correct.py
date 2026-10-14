import csv
import os
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, "analysis"))
from stats import two_sided_p, adjust_fdr, FDR_Q

RAW = os.path.join(ROOT, "raw", "run_2024-11-02_plate_reader.csv")


def load_raw(path):
    features = {}
    with open(path, newline="") as handle:
        for row in csv.DictReader(handle):
            feature_id = row["feature_id"]
            features.setdefault(feature_id, {"A": [], "B": []})
            features[feature_id][row["group"]].append(float(row["value"]))
    return features


def bonferroni_significant(p_by_feature, alpha=0.05):
    m = len(p_by_feature)
    return sorted(
        fid for fid, p in p_by_feature.items() if p * m <= alpha
    )


def bh_significant(p_by_feature, q=0.05):
    adj = adjust_fdr(p_by_feature)
    return sorted(fid for fid, qv in adj.items() if qv <= q)


def main():
    features = load_raw(RAW)
    pvals = {
        fid: two_sided_p(g["A"], g["B"]) for fid, g in features.items()
    }
    bh = bh_significant(pvals, FDR_Q)
    bonf = bonferroni_significant(pvals)
    print("n_features", len(features))
    print("BH", len(bh), bh)
    print("Bonferroni", len(bonf), bonf)

    # features with raw p < 0.05 but not in BH after full correction
    raw_sig = {fid for fid, p in pvals.items() if p < 0.05}
    kept_only = raw_sig - set(bh)
    print("raw p<0.05 not in full BH:", sorted(kept_only))

    # BH on prefiltered only (replicate bug)
    pre = {fid: p for fid, p in pvals.items() if p < 0.05}
    bh_pref = bh_significant(pre, FDR_Q)
    print("BH on prefiltered only", len(bh_pref), bh_pref)


if __name__ == "__main__":
    main()
