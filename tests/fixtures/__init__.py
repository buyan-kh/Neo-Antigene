"""Synthetic test data.

Nothing in here is derived from a real patient, and nothing in the repository
ever should be. Cohorts are generated from a seed so that failures are
reproducible and the suite does not flake.
"""

from .synthetic import (
    SyntheticCohort,
    SyntheticPatient,
    build_cohort,
    simulate_assay_results,
)

__all__ = [
    "SyntheticCohort",
    "SyntheticPatient",
    "build_cohort",
    "simulate_assay_results",
]
