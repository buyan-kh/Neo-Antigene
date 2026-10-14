from __future__ import annotations

import csv
from math import erfc, sqrt
from pathlib import Path
from typing import Dict, List, Set, Tuple

KNOWN_SIGMA = 1.0


def read_features(path: Path) -> Dict[str, Dict[str, List[float]]]:
    grouped: Dict[str, Dict[str, List[float]]] = {}
    with open(path, "r", encoding="utf-8") as handle:
        for row in csv.DictReader(handle):
            feature = row["feature_id"]
            group = row["group"]
            grouped.setdefault(feature, {"A": [], "B": []})
            grouped[feature][group].append(float(row["value"]))
    return grouped


def two_sided_p(values_a: List[float], values_b: List[float]) -> float:
    n_a, n_b = len(values_a), len(values_b)
    mean_a = sum(values_a) / n_a
    mean_b = sum(values_b) / n_b
    standard_error = KNOWN_SIGMA * sqrt(1.0 / n_a + 1.0 / n_b)
    z = (mean_a - mean_b) / standard_error
    return erfc(abs(z) / sqrt(2.0))


def p_values(path: Path) -> Dict[str, float]:
    return {
        feature: two_sided_p(groups["A"], groups["B"])
        for feature, groups in read_features(path).items()
    }


def _ordered(p_by_feature: Dict[str, float]) -> List[Tuple[str, float]]:
    return sorted(p_by_feature.items(), key=lambda item: (item[1], item[0]))


def benjamini_hochberg(p_by_feature: Dict[str, float], q: float) -> Set[str]:
    ordered = _ordered(p_by_feature)
    total = len(ordered)
    largest_index = 0
    for rank, (_, p_value) in enumerate(ordered, start=1):
        if p_value <= q * rank / total:
            largest_index = rank
    return {feature for feature, _ in ordered[:largest_index]}


def bonferroni(p_by_feature: Dict[str, float], alpha: float) -> Set[str]:
    total = len(p_by_feature)
    return {
        feature
        for feature, p_value in p_by_feature.items()
        if p_value <= alpha / total
    }


def uncorrected(p_by_feature: Dict[str, float], alpha: float) -> Set[str]:
    return {
        feature for feature, p_value in p_by_feature.items() if p_value < alpha
    }
