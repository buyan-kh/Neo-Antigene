"""Weight refitting.

The property under test is that a cross-validated AUC means what it says.
Rows within a patient are not independent — `clonality`, `expression` and
`tumor_selectivity` are one value per patient whenever purity or RNA is
missing, which is the situation in both benchmark cases — so a row-wise split
lets the model recognize a held-out patient and replay its base response rate.
That reads as skill and is not skill.

Two tests carry the weight: one where the labels have no peptide-level signal
at all and the estimate must stay near chance, and one where they do and it
must still be found. A fix that only passed the first would be a fix that
reports 0.5 forever.
"""

import numpy as np
import pytest
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import StratifiedKFold, cross_val_score

from neoantigene.assays.schema import AssayCall, AssayResult, AssayType
from neoantigene.learning.dataset import TrainingSet, build_training_set
from neoantigene.learning.train import (
    MIN_GROUPS_FOR_CV,
    InsufficientData,
    fit_weights,
)
from neoantigene.scoring.rank import SCORED_FEATURES

PRESENTATION = SCORED_FEATURES.index("presentation")
PATIENT_CONSTANT = [
    SCORED_FEATURES.index(name) for name in ("clonality", "expression", "tumor_selectivity")
]


def _dataset(
    patients: int,
    per_patient: int,
    rates: list[float] | None = None,
    signal: float = 0.0,
    seed: int = 0,
) -> TrainingSet:
    """A cohort whose response rate varies by patient, with tunable real signal.

    `signal` is how much `presentation` actually drives the label. At 0 the
    label depends only on which patient the row came from, so the only way to
    score above chance on a held-out patient is to have leaked.
    """
    rng = np.random.default_rng(seed)
    rates = rates or [0.02, 0.05, 0.10, 0.15, 0.25, 0.35][:patients]

    dataset = TrainingSet()
    for patient in range(patients):
        # Constant across every row of this patient, as when purity and RNA are
        # unavailable and the features fall back to one sample-level value.
        constants = rng.uniform(0.2, 0.8, size=len(PATIENT_CONSTANT))
        for _ in range(per_patient):
            features = rng.uniform(0.0, 1.0, size=len(SCORED_FEATURES))
            for slot, value in zip(PATIENT_CONSTANT, constants, strict=True):
                features[slot] = value
            probability = rates[patient] + signal * (features[PRESENTATION] - 0.5)
            dataset.x.append([float(v) for v in features])
            dataset.y.append(int(rng.random() < probability))
            dataset.groups.append(f"PT{patient}")
            dataset.keys.append(f"PT{patient}|PEP{len(dataset.y)}|HLA-A*02:01")
    return dataset


def _rowwise_auc(dataset: TrainingSet, folds: int = 5) -> float:
    """The estimate a row-wise split would have reported, for comparison."""
    x = np.asarray(dataset.x, dtype=float)
    y = np.asarray(dataset.y, dtype=int)
    model = LogisticRegression(C=0.3, class_weight="balanced", max_iter=2000)
    splitter = StratifiedKFold(n_splits=folds, shuffle=True, random_state=0)
    return float(cross_val_score(model, x, y, cv=splitter, scoring="roc_auc").mean())


class TestGroupedCrossValidation:
    def test_patient_identity_alone_does_not_become_signal(self):
        """Labels driven only by patient must not produce a skilful-looking AUC."""
        dataset = _dataset(patients=6, per_patient=40, signal=0.0, seed=0)
        report = fit_weights(dataset)

        assert report.cv_auc is not None
        assert report.cv_auc < 0.6, (
            f"grouped cv_auc {report.cv_auc:.3f} on labels with no peptide-level signal: "
            f"patient information is still leaking across the split"
        )
        assert report.cv_auc < _rowwise_auc(dataset), (
            "grouping by patient did not reduce the estimate, so either the split is not "
            "grouping or this fixture no longer reproduces the leak"
        )

    def test_real_peptide_level_signal_is_still_recovered(self):
        """The grouped split must not have cost all of its power."""
        dataset = _dataset(patients=6, per_patient=40, signal=0.9, seed=1)
        report = fit_weights(dataset)

        assert report.cv_auc is not None
        assert report.cv_auc > 0.65, (
            f"grouped cv_auc {report.cv_auc:.3f} failed to find a strong planted "
            f"presentation effect"
        )

    def test_the_scheme_and_patient_count_are_reported(self):
        report = fit_weights(_dataset(patients=6, per_patient=40, seed=2))
        assert report.groups == 6
        assert "grouped" in report.cv_scheme
        assert "patient" in report.cv_scheme
        assert "patients=6" in report.describe()

    def test_too_few_patients_yields_no_estimate_rather_than_a_leaky_one(self):
        """One patient cannot be held out, so there is no honest number to give."""
        dataset = _dataset(patients=2, per_patient=40, rates=[0.2, 0.3], seed=3)
        report = fit_weights(dataset)

        assert report.groups < MIN_GROUPS_FOR_CV
        assert report.cv_auc is None
        assert report.cv_scheme == "none"
        assert "generalization is unmeasured" in report.describe()

    def test_weights_are_still_fitted_when_cv_is_impossible(self):
        report = fit_weights(_dataset(patients=2, per_patient=40, rates=[0.2, 0.3], seed=4))
        assert report.n == 80
        assert set(report.coefficients) == set(SCORED_FEATURES)

    def test_insufficient_data_is_still_refused(self):
        with pytest.raises(InsufficientData):
            fit_weights(_dataset(patients=3, per_patient=4, rates=[0.3, 0.3, 0.3], seed=5))


class TestGroupProvenance:
    def test_build_training_set_carries_the_patient_of_origin(self):
        results = [
            AssayResult(
                sample_id=sample,
                peptide=peptide,
                allele="HLA-A*02:01",
                assay=AssayType.IFNG_ELISPOT,
                call=AssayCall.POSITIVE if peptide.endswith("K") else AssayCall.NEGATIVE,
            )
            for sample, peptide in (
                ("PT1", "VVGADGVGK"),
                ("PT1", "GADGVGKSA"),
                ("PT2", "VVGAVGVGK"),
            )
        ]
        features = {r.key: dict.fromkeys(SCORED_FEATURES, 0.5) for r in results}
        dataset = build_training_set(results, features)

        assert dataset.groups == ["PT1", "PT1", "PT2"]
        assert dataset.group_count == 2
        assert len(dataset.groups) == len(dataset)

    def test_indeterminate_rows_do_not_desynchronize_groups(self):
        """Skipped rows must be skipped in every parallel list at once."""
        results = [
            AssayResult(
                sample_id=sample,
                peptide=peptide,
                allele="HLA-A*02:01",
                assay=AssayType.IFNG_ELISPOT,
                call=call,
            )
            for sample, peptide, call in (
                ("PT1", "VVGADGVGK", AssayCall.POSITIVE),
                ("PT2", "GADGVGKSA", AssayCall.INDETERMINATE),
                ("PT3", "VVGAVGVGK", AssayCall.NEGATIVE),
            )
        ]
        features = {r.key: dict.fromkeys(SCORED_FEATURES, 0.5) for r in results}
        dataset = build_training_set(results, features)

        assert dataset.groups == ["PT1", "PT3"]
        assert len(dataset.x) == len(dataset.y) == len(dataset.groups) == 2
