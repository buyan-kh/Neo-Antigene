import csv
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "analysis"))
from stats import two_sided_p, adjust_fdr, FDR_Q, ALPHA

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


features = load(RAW)
p_all = {fid: two_sided_p(g["A"], g["B"]) for fid, g in features.items()}
prefilter = {fid for fid, p in p_all.items() if p < ALPHA}
adj_all = adjust_fdr(p_all)
adj_sub = adjust_fdr({fid: p_all[fid] for fid in prefilter})

reported = set()
with open(os.path.join(ROOT, "results", "reported_significant.csv")) as f:
    next(f)
    for line in f:
        reported.add(line.strip())

sig_all = {fid for fid, q in adj_all.items() if q <= FDR_Q}
sig_sub = {fid for fid, q in adj_sub.items() if q <= FDR_Q}

print("prefilter size", len(prefilter))
print("sig all", len(sig_all), sorted(sig_all))
print("sig sub", len(sig_sub), sorted(sig_sub))
print("reported matches sub", reported == sig_sub)
print("false pos vs all", sorted(reported - sig_all))
print("missed vs all", sorted(sig_all - reported))
