import csv
import os
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, "analysis"))
from stats import two_sided_p, adjust_fdr, FDR_Q

RAW = os.path.join(ROOT, "raw", "run_2024-11-02_plate_reader.csv")


def load_raw():
    features = {}
    roles = {}
    with open(RAW, newline="") as handle:
        for row in csv.DictReader(handle):
            fid = row["feature_id"]
            roles[fid] = row["role"]
            features.setdefault(fid, {"A": [], "B": []})
            features[fid][row["group"]].append(float(row["value"]))
    return features, roles


def bh_screen(feature_ids, features):
    p_by = {fid: two_sided_p(features[fid]["A"], features[fid]["B"]) for fid in feature_ids}
    adj = adjust_fdr(p_by)
    sig = sorted(fid for fid, q in adj.items() if q <= FDR_Q)
    return sig, p_by, adj


def bonferroni_screen(feature_ids, features, alpha=0.05):
    p_by = {fid: two_sided_p(features[fid]["A"], features[fid]["B"]) for fid in feature_ids}
    m = len(feature_ids)
    sig = sorted(fid for fid, p in p_by.items() if p * m <= alpha)
    return sig, p_by


def main():
    features, roles = load_raw()
    screen_ids = sorted(fid for fid in features if roles[fid] == "screen")
    spike_ids = sorted(fid for fid in features if roles[fid] == "spike_in")
    print("screen features:", len(screen_ids))
    print("spike features:", len(spike_ids))

    bh_sig, p_by, adj = bh_screen(screen_ids, features)
    bonf_sig, _ = bonferroni_screen(screen_ids, features)
    print("BH significant:", len(bh_sig))
    print("Bonferroni significant:", len(bonf_sig))

    # prefilter set (what select_candidates keeps)
    from stats import ALPHA

    kept = [fid for fid in screen_ids if p_by[fid] < ALPHA]
    print("prefilter kept screen:", len(kept))
    bh_prefilter, _, _ = bh_screen(kept, features)
    print("BH on prefiltered only:", len(bh_prefilter))

    # compare lists
    if bh_sig != bonf_sig:
        only_bh = set(bh_sig) - set(bonf_sig)
        only_bonf = set(bonf_sig) - set(bh_sig)
        print("BH only:", len(only_bh), sorted(only_bh)[:10])
        print("Bonf only:", len(only_bonf), sorted(only_bonf)[:10])

    print("\nBH list:")
    for x in bh_sig:
        print(x)


if __name__ == "__main__":
    main()
