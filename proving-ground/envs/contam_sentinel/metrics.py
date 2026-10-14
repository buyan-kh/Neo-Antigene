from __future__ import annotations

from typing import List, Optional, Sequence, Tuple


def auroc(labels: Sequence[int], scores: Sequence[float]) -> Optional[float]:
    if len(labels) != len(scores) or not labels:
        return None
    positives = sum(1 for label in labels if label == 1)
    negatives = len(labels) - positives
    if positives == 0 or negatives == 0:
        return None

    ordered: List[Tuple[float, int]] = sorted(
        zip(scores, labels), key=lambda item: item[0]
    )
    ranks = [0.0] * len(ordered)
    index = 0
    while index < len(ordered):
        stop = index
        while stop + 1 < len(ordered) and ordered[stop + 1][0] == ordered[index][0]:
            stop += 1
        average_rank = (index + stop + 2) / 2.0
        for position in range(index, stop + 1):
            ranks[position] = average_rank
        index = stop + 1

    rank_sum = sum(
        rank for rank, (_, label) in zip(ranks, ordered) if label == 1
    )
    return (rank_sum - positives * (positives + 1) / 2.0) / float(
        positives * negatives
    )
