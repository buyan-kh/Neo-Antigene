import csv
import os
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, "analysis"))
from stats import two_sided_p, adjust_fdr, ALPHA, FDR_Q

RAW = os.path.join(ROOT, "raw", "run_2024-11-02_plate_reader.csv")


def load_features(path):
    features = {}
    with open(path, newline="") as handle:
        for row in csv.DictReader(handle):
            feature_id = row["feature_id"]
            features.setdefault(feature_id, {"A": [], "B": []})
            features[feature_id][row["group"]].append(float(row["value"]))
    return features


def main():
    features = load_features(RAW)
    p_by = {fid: two_sided_p(g["A"], g["B"]) for fid, g in features.items()}
    m = len(p_by)
    raw_sig = sum(1 for p in p_by.values() if p < ALPHA)
    raw_sig_le = sum(1 for p in p_by.values() if p <= ALPHA)
    print(f"total features: {m}")
    print(f"raw p < {ALPHA}: {raw_sig}")
    print(f"raw p <= {ALPHA}: {raw_sig_le}")

    adj = adjust_fdr(p_by)
    bh_sig = sorted(fid for fid, q in adj.items() if q <= FDR_Q)
    print(f"BH q <= {FDR_Q}: {len(bh_sig)}")
    for fid in bh_sig:
        print(fid)

    # Bonferroni
    bonf = sorted(fid for fid, p in p_by.items() if p * m <= 0.05)
    print(f"Bonferroni: {len(bonf)}")


if __name__ == "__main__":
    main()
