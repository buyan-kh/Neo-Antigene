from __future__ import annotations

from functools import lru_cache
from math import exp, lgamma, log
from typing import List, Optional, Tuple

POWER_DEFINITION = "fishers_noncentral_hypergeometric"
MIN_N_DEFINITION = "smallest_n_with_power_at_or_above_target"


def _log_comb(n: int, k: int) -> float:
    if k < 0 or k > n or n < 0:
        return float("-inf")
    return lgamma(n + 1) - lgamma(k + 1) - lgamma(n - k + 1)


def pmf(N: int, K: int, n: int, odds_ratio: float) -> Tuple[int, List[float]]:
    low = max(0, n - (N - K))
    high = min(n, K)
    if high < low:
        return low, []

    log_psi = log(odds_ratio)
    log_weights = [
        _log_comb(K, k) + _log_comb(N - K, n - k) + k * log_psi
        for k in range(low, high + 1)
    ]
    peak = max(log_weights)
    weights = [exp(value - peak) for value in log_weights]
    total = sum(weights)
    return low, [value / total for value in weights]


def _upper_tails(probabilities: List[float]) -> List[float]:
    tails = [0.0] * (len(probabilities) + 1)
    running = 0.0
    for index in range(len(probabilities) - 1, -1, -1):
        running += probabilities[index]
        tails[index] = running
    return tails


def tail_at_least(N: int, K: int, n: int, odds_ratio: float, threshold: int) -> float:
    low, probabilities = pmf(N, K, n, odds_ratio)
    index = threshold - low
    if index <= 0:
        return 1.0 if probabilities else 0.0
    if index >= len(probabilities):
        return 0.0
    return _upper_tails(probabilities)[index]


def critical_value(N: int, K: int, n: int, alpha: float) -> int:
    low, probabilities = pmf(N, K, n, 1.0)
    if not probabilities:
        return low
    tails = _upper_tails(probabilities)
    for index in range(len(probabilities)):
        if tails[index] <= alpha:
            return low + index
    return low + len(probabilities)


@lru_cache(maxsize=None)
def power(N: int, K: int, n: int, alpha: float, odds_ratio: float) -> float:
    threshold = critical_value(N, K, n, alpha)
    return tail_at_least(N, K, n, odds_ratio, threshold)


@lru_cache(maxsize=None)
def min_sample_size(
    N: int, K: int, alpha: float, odds_ratio: float, target_power: float
) -> Optional[int]:
    for n in range(1, N + 1):
        if power(N, K, n, alpha, odds_ratio) >= target_power:
            return n
    return None


def p_value(N: int, K: int, n: int, observed: int) -> float:
    return tail_at_least(N, K, n, 1.0, observed)
