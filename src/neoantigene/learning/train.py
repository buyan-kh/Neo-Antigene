"""Refit ranking weights on assay outcomes.

Same functional form as the shipped scorer, so a refit is a drop-in config
change rather than a model swap. Regularization is strong by default because
early datasets are tiny and heavily biased toward previously top-ranked
peptides.
"""

from __future__ import annotations

from pydantic import BaseModel, Field

from ..config import ScoringWeights
from .dataset import TrainingSet

MIN_EXAMPLES = 40
MIN_PER_CLASS = 8


class InsufficientData(RuntimeError):
    pass


class FitReport(BaseModel):
    n: int
    positives: int
    weights: ScoringWeights
    cv_auc: float | None = None
    coefficients: dict[str, float] = Field(default_factory=dict)

    def describe(self) -> str:
        auc = f"{self.cv_auc:.3f}" if self.cv_auc is not None else "n/a"
        lines = [f"n={self.n} positives={self.positives} cv_auc={auc}", "coefficients:"]
        for name, value in sorted(self.coefficients.items(), key=lambda kv: -abs(kv[1])):
            lines.append(f"  {name:<18} {value:+.3f}")
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
        from sklearn.model_selection import StratifiedKFold, cross_val_score
    except ImportError as exc:  # pragma: no cover - depends on optional extra
        raise ImportError(
            "refitting requires the learning extra: uv sync --extra learning"
        ) from exc

    x = np.asarray(dataset.x, dtype=float)
    y = np.asarray(dataset.y, dtype=int)

    model = LogisticRegression(C=regularization_c, class_weight="balanced", max_iter=2000)

    cv_auc: float | None = None
    folds = min(cv_folds, dataset.positives, dataset.negatives)
    if folds >= 2:
        splitter = StratifiedKFold(n_splits=folds, shuffle=True, random_state=0)
        cv_auc = float(cross_val_score(model, x, y, cv=splitter, scoring="roc_auc").mean())

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
    )


def compare_to_prior(prior: ScoringWeights, fitted: ScoringWeights) -> list[str]:
    """Human-readable diff, for reviewing a refit before adopting it."""
    prior_values = prior.model_dump()
    fitted_values = fitted.model_dump()
    return [
        f"{name:<18} {prior_values[name]:+.3f} -> {fitted_values[name]:+.3f}  "
        f"({fitted_values[name] - prior_values[name]:+.3f})"
        for name in sorted(prior_values)
    ]
