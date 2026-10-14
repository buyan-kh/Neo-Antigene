import csv
import json
import sys
from fractions import Fraction
from math import comb


def load_population(path):
    with open(path, newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    flags = {row["item_id"]: int(row["in_group"]) for row in rows}
    assert len(flags) == len(rows)
    return flags


def load_sample(path):
    with open(path, newline="") as handle:
        ids = [row["item_id"] for row in csv.DictReader(handle, delimiter="\t")]
    assert len(set(ids)) == len(ids)
    return ids


def weights(N, K, n, odds_ratio):
    low = max(0, n - (N - K))
    high = min(n, K)
    num, den = odds_ratio.numerator, odds_ratio.denominator
    w = [
        comb(K, k) * comb(N - K, n - k) * num**k * den ** (n - k)
        for k in range(low, high + 1)
    ]
    return low, w, sum(w)


def upper_tail(low, w, total, threshold):
    start = min(max(0, threshold - low), len(w))
    return Fraction(sum(w[start:]), total)


def critical_overlap(N, K, n, alpha):
    low, w, total = weights(N, K, n, Fraction(1))
    tail = 0
    threshold = low + len(w)
    for i in range(len(w) - 1, -1, -1):
        tail += w[i]
        if Fraction(tail, total) > alpha:
            break
        threshold = low + i
    return threshold


def power(N, K, n, alpha, odds_ratio):
    threshold = critical_overlap(N, K, n, alpha)
    low, w, total = weights(N, K, n, odds_ratio)
    return threshold, upper_tail(low, w, total, threshold)


def min_sample_size(N, K, alpha, odds_ratio, target):
    for n in range(1, N + 1):
        _, achieved = power(N, K, n, alpha, odds_ratio)
        if achieved >= target:
            return n, achieved
    return None, None


def analyse(N, K, n, observed, alpha, odds_ratio, target):
    low0, w0, total0 = weights(N, K, n, Fraction(1))
    p_value = upper_tail(low0, w0, total0, observed)
    threshold, achieved = power(N, K, n, alpha, odds_ratio)
    n_req, power_at_n_req = min_sample_size(N, K, alpha, odds_ratio, target)
    return {
        "N": N,
        "K": K,
        "n": n,
        "observed": observed,
        "cutoff": threshold,
        "p_value": float(p_value),
        "power": float(achieved),
        "min_n": n_req,
        "power_at_min_n": float(power_at_n_req),
        "power_at_min_n_minus_1": float(power(N, K, n_req - 1, alpha, odds_ratio)[1]),
    }


def main():
    with open("study_design.json") as handle:
        design = json.load(handle)
    flags = load_population("population.tsv")
    sample = load_sample("sample.tsv")
    assert all(item in flags for item in sample)

    alpha = Fraction(str(design["alpha"]))
    target = Fraction(str(design["target_power"]))
    assert design["claimed_effect"]["type"] == "odds_ratio"
    odds_ratio = Fraction(str(design["claimed_effect"]["value"]))

    N = len(flags)
    K = sum(flags.values())
    n = len(sample)
    observed = sum(flags[item] for item in sample)

    actual = analyse(N, K, n, observed, alpha, odds_ratio, target)
    print("from data:  ", json.dumps(actual))

    if "--compare-design" in sys.argv:
        stale = analyse(
            design["population_size"],
            design["successes_in_population"],
            design["sample_size"],
            design["observed_overlap"],
            alpha,
            odds_ratio,
            target,
        )
        print("from design:", json.dumps(stale))

    if "--write" in sys.argv:
        underpowered = actual["power"] < float(target)
        submission = {
            "power": actual["power"],
            "underpowered": underpowered,
            "min_sample_size_for_target_power": actual["min_n"],
            "p_value": actual["p_value"],
            "endorse_claim": not underpowered,
            "reason_code": "underpowered" if underpowered else "adequately_powered",
        }
        with open("output/submission.json", "w") as handle:
            json.dump(submission, handle, indent=2)
            handle.write("\n")
        print(json.dumps(submission, indent=2))


if __name__ == "__main__":
    main()
