"""Per-feature screening statistics."""

from math import erfc, sqrt

KNOWN_SIGMA = 1.0
ALPHA = 0.05
FDR_Q = 0.05


def two_sided_p(values_a, values_b):
    """Two-sided z-test p-value for a difference in means."""
    n_a, n_b = len(values_a), len(values_b)
    mean_a = sum(values_a) / n_a
    mean_b = sum(values_b) / n_b
    standard_error = KNOWN_SIGMA * sqrt(1.0 / n_a + 1.0 / n_b)
    z = (mean_a - mean_b) / standard_error
    return erfc(abs(z) / sqrt(2.0))


def adjust_fdr(p_by_feature):
    """Benjamini-Hochberg adjusted p-values."""
    ordered = sorted(p_by_feature.items(), key=lambda item: (item[1], item[0]))
    m = len(ordered)
    adjusted = {}
    running = 1.0
    for rank in range(m, 0, -1):
        feature_id, p_value = ordered[rank - 1]
        running = min(running, p_value * m / rank)
        adjusted[feature_id] = running
    return adjusted


def screen(features):
    """Score every feature, FDR-adjust, and return the ones we call significant."""
    results = {}
    for feature_id, groups in features.items():
        results[feature_id] = two_sided_p(groups["A"], groups["B"])

    adjusted = adjust_fdr(results)
    significant = [
        feature_id for feature_id, q_value in adjusted.items() if q_value <= FDR_Q
    ]

    return results, sorted(significant)
