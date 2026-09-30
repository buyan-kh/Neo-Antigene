"""Join predicted features to assay outcomes.

The training set only ever contains peptides that were actually tested, which
means it is systematically drawn from the top of previous rankings. That
sampling bias is the reason active learning (rather than plain retraining) is
in the loop — see `neoantigene.learning.active`.
"""

from __future__ import annotations

import json
from collections.abc import Sequence
from pathlib import Path

from pydantic import BaseModel, Field

from ..assays.schema import AssayResult
from ..scoring.rank import SCORED_FEATURES


class TrainingSet(BaseModel):
    keys: list[str] = Field(default_factory=list)
    feature_names: tuple[str, ...] = SCORED_FEATURES
    x: list[list[float]] = Field(default_factory=list)
    y: list[int] = Field(default_factory=list)

    #: Patient of origin per row, parallel to `x` and `y`. Cross-validation
    #: must split on this rather than on rows: several features are constant
    #: within a patient whenever purity or RNA is unavailable, so a random
    #: split lets the model identify the held-out patient from its training
    #: rows and score its known response rate. See `learning.train`.
    groups: list[str] = Field(default_factory=list)

    def __len__(self) -> int:
        return len(self.y)

    @property
    def positives(self) -> int:
        return sum(self.y)

    @property
    def negatives(self) -> int:
        return len(self.y) - self.positives

    @property
    def group_count(self) -> int:
        return len(set(self.groups))


def load_feature_records(paths: Sequence[Path]) -> dict[str, dict[str, float]]:
    """Index `features.json` outputs by `sample_id|peptide|allele`.

    Also registers an allele-agnostic `sample|peptide|*` key so results
    recorded without a restricting allele still join.
    """
    index: dict[str, dict[str, float]] = {}
    for path in paths:
        entries = json.loads(Path(path).read_text())
        for entry in entries:
            features = entry["features"]
            index[f"{entry['sample_id']}|{entry['peptide']}|{entry['allele']}"] = features
            index.setdefault(f"{entry['sample_id']}|{entry['peptide']}|*", features)
    return index


def build_training_set(
    results: Sequence[AssayResult],
    features: dict[str, dict[str, float]],
    feature_names: tuple[str, ...] = SCORED_FEATURES,
) -> TrainingSet:
    dataset = TrainingSet(feature_names=feature_names)
    for result in results:
        label = result.label
        if label is None:
            continue
        values = features.get(result.key)
        if values is None:
            continue
        dataset.keys.append(result.key)
        dataset.x.append([float(values.get(name, 0.0)) for name in feature_names])
        dataset.y.append(label)
        dataset.groups.append(result.sample_id)
    return dataset
