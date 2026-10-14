"""Per-feature screening statistics."""

from math import erfc, sqrt

KNOWN_SIGMA = 1.0
ALPHA = 0.05


def two_sided_p(values_a, values_b):
    """Two-sided z-test p-value for a difference in means."""
    n_a, n_b = len(values_a), len(values_b)
    mean_a = sum(values_a) / n_a
    mean_b = sum(values_b) / n_b
    standard_error = KNOWN_SIGMA * sqrt(1.0 / n_a + 1.0 / n_b)
    z = (mean_a - mean_b) / standard_error
    return erfc(abs(z) / sqrt(2.0))


def screen(features):
    """Score every feature and return the ones we call significant."""
    results = {}
    for feature_id, groups in features.items():
        results[feature_id] = two_sided_p(groups["A"], groups["B"])

    significant = []
    for feature_id, p_value in results.items():
        if p_value < ALPHA:
            significant.append(feature_id)

    return results, sorted(significant)
