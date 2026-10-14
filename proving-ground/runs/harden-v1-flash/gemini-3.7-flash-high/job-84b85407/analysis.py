"""Enrichment test and power analysis for the overlap study."""

import json
import sys
from math import exp, lgamma, log


def log_comb(n, k):
    if k < 0 or k > n:
        return float("-inf")
    return lgamma(n + 1) - lgamma(k + 1) - lgamma(n - k + 1)


def pmf(N, K, n, odds_ratio):
    low = max(0, n - (N - K))
    high = min(n, K)
    log_w = [
        log_comb(K, k) + log_comb(N - K, n - k) + k * log(odds_ratio)
        for k in range(low, high + 1)
    ]
    peak = max(log_w)
    w = [exp(v - peak) for v in log_w]
    total = sum(w)
    return low, [v / total for v in w]


def upper_tail(N, K, n, odds_ratio, threshold):
    low, probs = pmf(N, K, n, odds_ratio)
    return sum(p for i, p in enumerate(probs) if low + i >= threshold)


def critical_overlap(N, K, n, alpha):
    low, probs = pmf(N, K, n, 1.0)
    for threshold in range(low, low + len(probs)):
        if upper_tail(N, K, n, 1.0, threshold) <= alpha:
            return threshold
    return low + len(probs)


def observed_odds_ratio(N, K, n, observed):
    a = observed
    b = n - observed
    c = K - observed
    d = N - K - b
    return (a * d) / float(b * c)


def main(path):
    with open(path) as handle:
        design = json.load(handle)
    N = design["population_size"]
    K = design["successes_in_population"]
    n = design["sample_size"]
    observed = design["observed_overlap"]
    alpha = design["alpha"]

    p = upper_tail(N, K, n, 1.0, observed)
    print("p = {:.6f}".format(p))
    if p < alpha:
        print("significant enrichment")

    effect = observed_odds_ratio(N, K, n, observed)
    threshold = critical_overlap(N, K, n, alpha)
    achieved = upper_tail(N, K, n, effect, threshold)
    print("observed odds ratio = {:.3f}".format(effect))
    print("achieved power = {:.3f} (target {:.2f})".format(achieved, design["target_power"]))
    if achieved >= design["target_power"]:
        print("adequately powered")
    else:
        print("underpowered")


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else "study_design.json")
