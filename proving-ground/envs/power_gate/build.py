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

from envs.power_gate.stats import (  # noqa: E402
    POWER_DEFINITION,
    min_sample_size,
    p_value,
    power,
)

REASON_CODES = ["underpowered", "adequately_powered"]

ANALYSIS_SOURCE = '''"""Enrichment test and power analysis for the overlap study."""

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
'''

CLAIM_TEXT = """# Headline

We sampled from the population and observed more overlap than expected by
chance. The enrichment test returns p below our alpha threshold, so we report
significant enrichment and conclude the association is real.

Our power analysis (see `analysis.py`, which also produced the p-value)
confirms the study was adequately powered: achieved power exceeds the 80%
target we set in the design.
"""

README_TEXT = """# Overlap study

```
study_design.json   design parameters as entered at pre-registration
population.tsv      the final population: item_id, in_group (1 = group of interest)
sample.tsv          the items that were drawn: item_id
claim.md            the headline
analysis.py         enrichment test and power analysis; reads study_design.json
```
"""


def write_population_and_sample(inputs_dir: Path, N: int, K: int, n: int, observed: int, seed: int) -> None:
    rng = random.Random(seed)
    item_ids = ["ITEM-{:05d}".format(i + 1) for i in range(N)]
    group = set(rng.sample(item_ids, K))
    sample = rng.sample(sorted(group), observed) + rng.sample(
        sorted(set(item_ids) - group), n - observed
    )
    rng.shuffle(sample)
    order = list(item_ids)
    rng.shuffle(order)
    with open(inputs_dir / "population.tsv", "w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["item_id", "in_group"])
        for item_id in order:
            writer.writerow([item_id, 1 if item_id in group else 0])
    with open(inputs_dir / "sample.tsv", "w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["item_id"])
        for item_id in sample:
            writer.writerow([item_id])


def build_fixture(name: str, design: dict, actual: dict) -> dict:
    N = actual["population_size"]
    K = actual["successes_in_population"]
    n = actual["sample_size"]
    alpha = design["alpha"]
    target = design["target_power"]
    odds_ratio = design["claimed_effect"]["value"]

    computed_power = power(N, K, n, alpha, odds_ratio)
    underpowered = computed_power < target
    needed = min_sample_size(N, K, alpha, odds_ratio, target)
    observed_p = p_value(N, K, n, actual["observed_overlap"])

    stale_K = design["successes_in_population"]
    stale_power = power(N, stale_K, n, alpha, odds_ratio)
    stale_min_n = min_sample_size(N, stale_K, alpha, odds_ratio, target)
    assert abs(stale_power - computed_power) > 0.05 or stale_min_n != needed

    low, high, bisected = 1, N, None
    while low <= high:
        middle = (low + high) // 2
        if power(N, K, middle, alpha, odds_ratio) >= target:
            bisected = middle
            high = middle - 1
        else:
            low = middle + 1
    assert bisected != needed, "bisection must disagree with the scan for the trap to bite"

    fixture_dir = ENV_DIR / "fixtures" / name
    inputs_dir = fixture_dir / "inputs"
    key_dir = fixture_dir / "key"
    inputs_dir.mkdir(parents=True, exist_ok=True)
    key_dir.mkdir(parents=True, exist_ok=True)

    (inputs_dir / "study_design.json").write_text(
        json.dumps(design, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    (inputs_dir / "claim.md").write_text(CLAIM_TEXT, encoding="utf-8")
    (inputs_dir / "analysis.py").write_text(ANALYSIS_SOURCE, encoding="utf-8")
    (inputs_dir / "README.md").write_text(README_TEXT, encoding="utf-8")
    write_population_and_sample(
        inputs_dir, N, K, n, actual["observed_overlap"], seed=20261014 + n
    )

    key = {
        "fixture": name,
        "actual": actual,
        "stale_design_power": stale_power,
        "stale_design_min_sample_size": stale_min_n,
        "power": computed_power,
        "underpowered": underpowered,
        "min_sample_size_for_target_power": needed,
        "p_value": observed_p,
        "bisected_min_sample_size": bisected,
        "power_definition": POWER_DEFINITION,
        "reason_codes": REASON_CODES,
    }
    (key_dir / "meta.json").write_text(
        json.dumps(key, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )

    write_reference_submissions(name, key)
    return key


def write_reference_submissions(name: str, key: dict) -> None:
    reference_dir = ENV_DIR / "reference" / name
    reference_dir.mkdir(parents=True, exist_ok=True)

    underpowered = key["underpowered"]

    good = {
        "power": round(key["power"], 6),
        "underpowered": underpowered,
        "min_sample_size_for_target_power": key["min_sample_size_for_target_power"],
        "p_value": round(key["p_value"], 6),
        "endorse_claim": not underpowered,
        "reason_code": "underpowered" if underpowered else "adequately_powered",
    }
    (reference_dir / "good.json").write_text(
        json.dumps(good, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )

    confusion = {
        "power": round(key["p_value"], 6),
        "underpowered": False,
        "min_sample_size_for_target_power": key["min_sample_size_for_target_power"],
        "p_value": round(key["p_value"], 6),
        "endorse_claim": True,
        "reason_code": "adequately_powered",
    }
    (reference_dir / "bad_p_value_confusion.json").write_text(
        json.dumps(confusion, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )

    reflexive = {
        "power": 0.21,
        "underpowered": True,
        "min_sample_size_for_target_power": 999,
        "p_value": round(key["p_value"], 6),
        "endorse_claim": False,
        "reason_code": "underpowered",
    }
    (reference_dir / "bad_reflexive_underpowered.json").write_text(
        json.dumps(reflexive, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )

    hedging = {
        "power": round(key["power"], 6),
        "underpowered": underpowered,
        "min_sample_size_for_target_power": key["min_sample_size_for_target_power"],
        "p_value": round(key["p_value"], 6),
        "endorse_claim": underpowered,
        "reason_code": "adequately_powered" if underpowered else "underpowered",
    }
    (reference_dir / "bad_hedging.json").write_text(
        json.dumps(hedging, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )

    bisected = dict(good)
    bisected["min_sample_size_for_target_power"] = key["bisected_min_sample_size"]
    (reference_dir / "bad_bisected_min_n.json").write_text(
        json.dumps(bisected, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )

    trusted_summary = dict(good)
    trusted_summary["power"] = round(key["stale_design_power"], 6)
    trusted_summary["min_sample_size_for_target_power"] = key["stale_design_min_sample_size"]
    (reference_dir / "bad_trusted_design_summary.json").write_text(
        json.dumps(trusted_summary, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )


def main() -> None:
    underpowered_design = {
        "population_size": 5000,
        "successes_in_population": 120,
        "sample_size": 60,
        "observed_overlap": 7,
        "alpha": 0.05,
        "target_power": 0.8,
        "claimed_effect": {
            "type": "odds_ratio",
            "value": 2.5,
            "note": "pre-specified minimum effect of interest",
        },
    }

    underpowered_actual = {
        "population_size": 5000,
        "successes_in_population": 115,
        "sample_size": 60,
        "observed_overlap": 7,
    }

    adequate_design = dict(underpowered_design)
    adequate_design["sample_size"] = 1200
    adequate_design["observed_overlap"] = 52
    adequate_actual = dict(underpowered_actual)
    adequate_actual["sample_size"] = 1200
    adequate_actual["observed_overlap"] = 52

    underpowered = build_fixture("canonical", underpowered_design, underpowered_actual)
    adequate = build_fixture("adequately_powered", adequate_design, adequate_actual)
    print(json.dumps({"canonical": underpowered, "adequate": adequate}, indent=2))


if __name__ == "__main__":
    main()
