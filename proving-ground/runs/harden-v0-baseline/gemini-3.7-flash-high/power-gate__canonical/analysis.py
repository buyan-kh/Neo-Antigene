"""Enrichment test for the overlap study."""

from math import comb


def hypergeometric_tail(N, K, n, observed):
    total = comb(N, n)
    return sum(
        comb(K, k) * comb(N - K, n - k) / total
        for k in range(observed, min(n, K) + 1)
    )


def main(design):
    tail = hypergeometric_tail(
        design["population_size"],
        design["successes_in_population"],
        design["sample_size"],
        design["observed_overlap"],
    )
    print("p = {:.6f}".format(tail))
    if tail < design["alpha"]:
        print("significant enrichment")
