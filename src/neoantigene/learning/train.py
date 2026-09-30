"""Refit ranking weights on assay outcomes.

Same functional form as the shipped scorer, so a refit is a drop-in config
change rather than a model swap. Regularization is strong by default because
early datasets are tiny and heavily biased toward previously top-ranked
peptides.

Cross-validation splits by patient, never by row. Rows within one patient are
not independent: `clonality`, `expression` and `tumor_selectivity` collapse to
the same value for every candidate whenever purity or RNA is unavailable, which
is exactly the situation in both benchmark cases, and overlapping registers of
one variant differ by a residue. A row-wise split therefore lets the model
recognize the held-out patient from its own training rows and predict that
patient's base response rate, which reads as skill. On six simulated patients
whose labels carried no peptide-level signal at all, a row-wise split reported
cv_auc 0.643 where the grouped split correctly reported 0.494.
"""

from __future__ import annotations

import logging
import random
from typing import Any

from pydantic import BaseModel, Field

from ..config import ScoringWeights
from .dataset import TrainingSet

logger = logging.getLogger(__name__)

MIN_EXAMPLES = 40
MIN_PER_CLASS = 8

#: Below this many patients a grouped split cannot be formed, so no honest
#: cross-validated estimate is available from this data at all.
MIN_GROUPS_FOR_CV = 3


class InsufficientData(RuntimeError):
    pass


class FitReport(BaseModel):
    n: int
    positives: int
    weights: ScoringWeights
    cv_auc: float | None = None
    coefficients: dict[str, float] = Field(default_factory=dict)

    #: How `cv_auc` was estimated, and over how many patients. Recorded because
    #: a cross-validated AUC is uninterpretable without knowing what was held
    #: out; `none` means it could not be estimated honestly and `cv_auc` is
    #: None rather than optimistic.
    cv_scheme: str = "none"
    groups: int = 0

    def describe(self) -> str:
        auc = f"{self.cv_auc:.3f}" if self.cv_auc is not None else "n/a"
        lines = [
            f"n={self.n} positives={self.positives} patients={self.groups} "
            f"cv_auc={auc} ({self.cv_scheme})",
            "coefficients:",
        ]
        for name, value in sorted(self.coefficients.items(), key=lambda kv: -abs(kv[1])):
            lines.append(f"  {name:<18} {value:+.3f}")
        if self.cv_auc is None:
            lines.append(
                f"  no cross-validated estimate: needs >= {MIN_GROUPS_FOR_CV} patients, "
                f"have {self.groups}. The weights are fitted; their generalization is "
                f"unmeasured."
            )
        return "\n".join(lines)


def _check_sufficient(dataset: TrainingSet) -> None:
    if len(dataset) < MIN_EXAMPLES:
        raise InsufficientData(
            f"need at least {MIN_EXAMPLES} labelled peptides to refit, have {len(dataset)}"
        )
    if dataset.positives < MIN_PER_CLASS or dataset.negatives < MIN_PER_CLASS:
        raise InsufficientData(
            f"need at least {MIN_PER_CLASS} of each class, have "
            f"{dataset.positives} positive / {dataset.negatives} negative"
        )


def fit_weights(
    dataset: TrainingSet,
    regularization_c: float = 0.3,
    cv_folds: int = 5,
) -> FitReport:
    _check_sufficient(dataset)

    try:
        import numpy as np
        from sklearn.linear_model import LogisticRegression
    except ImportError as exc:  # pragma: no cover - depends on optional extra
        raise ImportError(
            "refitting requires the learning extra: uv sync --extra learning"
        ) from exc

    x = np.asarray(dataset.x, dtype=float)
    y = np.asarray(dataset.y, dtype=int)
    groups = np.asarray(dataset.groups) if len(dataset.groups) == len(dataset) else None

    model = LogisticRegression(C=regularization_c, class_weight="balanced", max_iter=2000)
    cv_auc, cv_scheme = _cross_validate(model, x, y, groups, cv_folds)

    model.fit(x, y)
    coefficients = {
        name: float(value)
        for name, value in zip(dataset.feature_names, model.coef_[0], strict=True)
    }

    return FitReport(
        n=len(dataset),
        positives=dataset.positives,
        weights=ScoringWeights(bias=float(model.intercept_[0]), **coefficients),
        cv_auc=cv_auc,
        coefficients=coefficients,
        cv_scheme=cv_scheme,
        groups=dataset.group_count,
    )


def _cross_validate(
    model: Any,
    x: Any,
    y: Any,
    groups: Any,
    cv_folds: int,
) -> tuple[float | None, str]:
    """Held-out AUC with every patient's rows kept on one side of the split.

    Returns `(None, "none")` rather than a row-wise estimate when there are too
    few patients to group by. A leaky number here is worse than no number: it
    would be quoted as evidence the refit generalizes, at the exact moment the
    first labels arrive and nobody has an independent check.
    """
    import numpy as np
    from sklearn.base import clone
    from sklearn.metrics import roc_auc_score
    from sklearn.model_selection import StratifiedGroupKFold

    distinct = 0 if groups is None else len(set(groups.tolist()))
    if distinct < MIN_GROUPS_FOR_CV:
        logger.warning(
            "not cross-validating: %d patient(s) in the training set, need at least %d to "
            "hold one out. Splitting by row instead would inflate AUC, because features "
            "that are constant within a patient let the model identify the held-out rows.",
            distinct,
            MIN_GROUPS_FOR_CV,
        )
        return None, "none"

    folds = min(cv_folds, distinct)
    splitter = StratifiedGroupKFold(n_splits=folds, shuffle=True, random_state=0)

    scores: list[float] = []
    for train_index, test_index in splitter.split(x, y, groups):
        if len(set(y[train_index].tolist())) < 2 or len(set(y[test_index].tolist())) < 2:
            # A fold with one class in either half has no defined AUC.
            continue
        fitted = clone(model).fit(x[train_index], y[train_index])
        probabilities = fitted.predict_proba(x[test_index])[:, 1]
        scores.append(float(roc_auc_score(y[test_index], probabilities)))

    if not scores:
        logger.warning(
            "not cross-validating: no patient-grouped fold had both classes on both sides"
        )
        return None, "none"

    return float(np.mean(scores)), f"grouped {len(scores)}-fold by patient"


def bootstrap_weights(
    dataset: TrainingSet,
    rounds: int = 200,
    regularization_c: float = 0.3,
    seed: int = 0,
) -> list[ScoringWeights]:
    """An ensemble of refits, each on a resample of the labelled data.

    The spread across this ensemble is what `learning.active` needs to know
    which label would actually teach it something. A single fit cannot say
    that: its score is a point estimate, and the distance of that point from
    0.5 measures how uncertain the *outcome* is, not how uncertain the
    *weights* are. A peptide can sit at 0.5 because the model is confidently
    ambivalent about it, which is the least informative case, not the most.

    Resampling is at the patient level whenever more than one patient is
    available, for the same reason cross-validation groups by patient: the
    dominant uncertainty in a cohort of six is which six, not which peptides
    within them. With one patient there is no cluster structure to resample,
    so rows are drawn instead and the resulting spread understates the true
    uncertainty — it describes only this tumor.
    """
    _check_sufficient(dataset)

    try:
        import numpy as np
        from sklearn.linear_model import LogisticRegression
    except ImportError as exc:  # pragma: no cover - depends on optional extra
        raise ImportError(
            "refitting requires the learning extra: uv sync --extra learning"
        ) from exc

    x = np.asarray(dataset.x, dtype=float)
    y = np.asarray(dataset.y, dtype=int)
    indices = _resample_indices(dataset, rounds, seed)

    ensemble: list[ScoringWeights] = []
    for rows in indices:
        labels = y[rows]
        if len(set(labels.tolist())) < 2:
            # A resample with one class cannot be fitted; dropping it biases the
            # ensemble slightly toward balanced draws, which is preferable to
            # inventing a member.
            continue
        model = LogisticRegression(C=regularization_c, class_weight="balanced", max_iter=2000)
        model.fit(x[rows], labels)
        ensemble.append(
            ScoringWeights(
                bias=float(model.intercept_[0]),
                **{
                    name: float(value)
                    for name, value in zip(dataset.feature_names, model.coef_[0], strict=True)
                },
            )
        )

    if len(ensemble) < 2:
        raise InsufficientData(
            f"could not fit a usable ensemble: only {len(ensemble)} of {rounds} resamples "
            f"contained both classes"
        )
    return ensemble


def _resample_indices(dataset: TrainingSet, rounds: int, seed: int) -> list[list[int]]:
    """Row indices for each bootstrap round, clustered by patient when possible."""
    rng = random.Random(seed)
    by_patient: dict[str, list[int]] = {}
    for row, patient in enumerate(dataset.groups):
        by_patient.setdefault(patient, []).append(row)

    if len(by_patient) < 2:
        rows = list(range(len(dataset)))
        return [[rng.choice(rows) for _ in rows] for _ in range(rounds)]

    patients = list(by_patient)
    draws: list[list[int]] = []
    for _ in range(rounds):
        selected: list[int] = []
        for _ in patients:
            selected.extend(by_patient[rng.choice(patients)])
        draws.append(selected)
    return draws


def compare_to_prior(prior: ScoringWeights, fitted: ScoringWeights) -> list[str]:
    """Human-readable diff, for reviewing a refit before adopting it."""
    prior_values = prior.model_dump()
    fitted_values = fitted.model_dump()
    return [
        f"{name:<18} {prior_values[name]:+.3f} -> {fitted_values[name]:+.3f}  "
        f"({fitted_values[name] - prior_values[name]:+.3f})"
        for name in sorted(prior_values)
    ]
