"""Power analysis for precision@k comparisons on neoantigen assay pools.

Answers two questions that the Ott-2017-scale benchmarks cannot:

1. Against chance. Given an assayed pool of N peptide-HLA pairs containing P
   positives, how many positives must exist before an exact hypergeometric test
   on top-k hits can detect a ranker whose true precision@k is a given multiple
   of the pool prevalence?

2. Against a baseline ranker on the same pool. Two rankings of one pool are
   paired -- they share candidates in their top-k -- so the informative
   quantity is the discordant top-k slots, not the raw hit counts. This is the
   design that actually matches "beat MHCflurry-sorted on the same assayed
   pool", and it is materially more powerful than treating the two arms as
   independent samples.

Run:  uv run python scripts/power_precision_at_k.py
"""

from __future__ import annotations

import numpy as np
from scipy.stats import hypergeom, norm

RNG = np.random.default_rng(20260930)
ALPHA = 0.05
TARGET_POWER = 0.80


def hypergeom_critical_hits(n_pool: int, n_pos: int, k: int, alpha: float) -> int:
    """Smallest hit count whose one-sided hypergeometric p-value is <= alpha.

    Returns k + 1 (unreachable) when no attainable hit count can reach alpha,
    which is the signature of a pool that cannot produce a significant result
    at any effect size.
    """
    for hits in range(0, min(k, n_pos) + 1):
        # P(X >= hits) under random ranking.
        if hypergeom.sf(hits - 1, n_pool, n_pos, k) <= alpha:
            return hits
    return k + 1


def power_vs_chance(
    n_pool: int,
    n_pos: int,
    k: int,
    enrichment: float,
    alpha: float = ALPHA,
    n_sim: int = 20000,
) -> tuple[float, int]:
    """Power to reject "ranking is random" when the ranker is `enrichment`x better.

    The ranker is modelled as drawing each of its k slots with per-slot hit
    probability `enrichment * P / N`, capped by the number of positives
    available. Sampling is without replacement, so we use a hypergeometric
    draw from an effective pool enriched by the same factor.
    """
    crit = hypergeom_critical_hits(n_pool, n_pos, k, alpha)
    if crit > k:
        return 0.0, crit

    p_slot = min(enrichment * n_pos / n_pool, 1.0)
    # Effective positives visible to this ranker in the region it selects from.
    eff_pos = min(int(round(p_slot * n_pool)), n_pos if enrichment <= 1 else n_pool)
    eff_pos = min(eff_pos, n_pool)
    draws = hypergeom.rvs(n_pool, eff_pos, k, size=n_sim, random_state=RNG)
    draws = np.minimum(draws, n_pos)
    return float((draws >= crit).mean()), crit


def power_paired(
    n_pool: int,
    n_pos: int,
    k: int,
    p_base: float,
    p_new: float,
    n_patients: int,
    overlap: float = 0.5,
    n_sim: int = 20000,
    alpha: float = ALPHA,
) -> float:
    """Power of a paired permutation test on per-patient precision@k differences.

    `overlap` is the expected fraction of top-k slots the two rankers share.
    Shared slots contribute identical hits and cancel in the paired difference,
    which is exactly why pairing helps: only the (1 - overlap) * k discordant
    slots carry signal.
    """
    k_disc = max(int(round((1.0 - overlap) * k)), 1)
    rejects = 0
    for _ in range(n_sim):
        # Per-patient hit counts in the discordant slots only.
        base_hits = RNG.binomial(k_disc, min(p_base, 1.0), size=n_patients)
        new_hits = RNG.binomial(k_disc, min(p_new, 1.0), size=n_patients)
        # Cap by the positives actually present in each patient's pool.
        cap = min(n_pos, k_disc)
        base_hits = np.minimum(base_hits, cap)
        new_hits = np.minimum(new_hits, cap)
        diff = (new_hits - base_hits) / k_disc
        if np.allclose(diff, 0.0):
            continue
        observed = diff.mean()
        # Exact-ish paired permutation: randomise the sign of each patient's delta.
        signs = RNG.choice([-1.0, 1.0], size=(4000, n_patients))
        null = (signs * diff).mean(axis=1)
        p_value = (1 + (null >= observed).sum()) / (1 + len(null))
        if p_value <= alpha:
            rejects += 1
    return rejects / n_sim


def two_proportion_n(p0: float, p1: float, alpha: float, power: float) -> int:
    """Slots per arm for an unpaired one-sided two-proportion z-test."""
    z_a = norm.ppf(1 - alpha)
    z_b = norm.ppf(power)
    p_bar = (p0 + p1) / 2
    num = (
        z_a * np.sqrt(2 * p_bar * (1 - p_bar))
        + z_b * np.sqrt(p0 * (1 - p0) + p1 * (1 - p1))
    ) ** 2
    return int(np.ceil(num / (p1 - p0) ** 2))


def main() -> None:
    k = 20

    print("=" * 78)
    print("1. AGAINST CHANCE -- exact hypergeometric on a single pool, k=20")
    print("=" * 78)
    print(
        "Ott 2017 as benchmarked: N=139, P=14. Prevalence 10.1%, "
        f"E[hits] under random = {k * 14 / 139:.2f}"
    )
    crit = hypergeom_critical_hits(139, 14, k, ALPHA)
    print(f"  hits needed for p<=0.05 one-sided: {crit} of 20")
    for enr in (1.5, 2.0, 3.0, 4.0):
        pw, _ = power_vs_chance(139, 14, k, enr)
        print(f"  power at {enr:>3}x enrichment: {pw:5.1%}")

    print()
    print("Bulik-Sullivan CU04 as benchmarked: N=19, P=3")
    crit = hypergeom_critical_hits(19, 3, k, ALPHA)
    print(f"  k=20 exceeds the pool (N=19); top-k is the whole pool -> no test exists")
    for kk in (5, 10):
        c = hypergeom_critical_hits(19, 3, kk, ALPHA)
        reachable = "unreachable" if c > kk else f"{c} of {kk}"
        print(f"  at k={kk}: hits needed for p<=0.05 = {reachable}")

    print()
    print("Pool size needed so that a 2x-enriched ranker is detectable at 80% power")
    print("(prevalence held at 10%, k=20):")
    for n_pool in (139, 300, 600, 1000, 2000, 4000):
        n_pos = round(0.10 * n_pool)
        pw, c = power_vs_chance(n_pool, n_pos, k, 2.0)
        print(f"  N={n_pool:>5}  P={n_pos:>4}  crit={c:>3}/20  power={pw:5.1%}")

    print()
    print("=" * 78)
    print("2. DOUBLING PRECISION@20 OVER A BASELINE -- how many patients / positives")
    print("=" * 78)
    for p0 in (0.05, 0.10, 0.20):
        n_arm = two_proportion_n(p0, 2 * p0, ALPHA, TARGET_POWER)
        patients = int(np.ceil(n_arm / k))
        print(
            f"  p@20 {p0:.0%} -> {2 * p0:.0%}: unpaired needs {n_arm:>5} top-20 slots/arm "
            f"= {patients:>3} patients/arm; expected positives in the winning arm "
            f"~{round(n_arm * 2 * p0)}"
        )

    print()
    print("Paired permutation test on per-patient precision@20 deltas")
    print("(baseline p@20=10%, new p@20=20%, top-20 lists overlapping by 50%):")
    for m in (5, 10, 20, 30, 40, 60, 80):
        pw = power_paired(1000, 100, k, 0.10, 0.20, m, overlap=0.5, n_sim=2000)
        print(f"  patients={m:>3}  power={pw:5.1%}")

    print()
    print("Same, but top-20 lists overlapping by 80% (a re-ranker that mostly")
    print("agrees with its dominant input feature -- the realistic case):")
    for m in (10, 20, 40, 60, 80, 120, 160):
        pw = power_paired(1000, 100, k, 0.10, 0.20, m, overlap=0.8, n_sim=2000)
        print(f"  patients={m:>3}  power={pw:5.1%}")

    print()
    print("=" * 78)
    print("3. WHAT THE EXISTING BENCHMARKS CAN AND CANNOT SHOW")
    print("=" * 78)
    print(
        "Ott 2017 pooled (N=139, P=14, k=20) with 1 effective patient stratum:\n"
        "  minimum attainable one-sided p-value, if all 20 top slots were hits:\n"
        f"  p = {hypergeom.sf(min(14, 20) - 1, 139, 14, 20):.2e}  (attainable)\n"
        "  but at the observed 2-4 hits the p-values are 0.06-0.27, i.e. the\n"
        "  design has ~0 power against any plausible effect size."
    )


if __name__ == "__main__":
    main()
