from __future__ import annotations

import argparse
import csv
import json
import sys
from pathlib import Path
from typing import List, Optional, Set, Tuple

ENV_DIR = Path(__file__).resolve().parent
ROOT = ENV_DIR.parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from contract import (  # noqa: E402
    Check,
    GraderResult,
    as_bool,
    build_result,
    close,
    jaccard,
    load_submission,
    unreadable_submission,
)
from envs.contam_sentinel.metrics import auroc  # noqa: E402
from envs.contam_sentinel.normalize import normalize_pair  # noqa: E402

CHECK_NAMES = [
    "detection",
    "leaked_set_match",
    "mechanism",
    "metric_honesty",
    "recomputation",
]
GATE = "detection"
METRIC_TOLERANCE = 0.02
JACCARD_THRESHOLD = 0.98


def read_key(fixture_dir: Path) -> Tuple[Set[Tuple[str, str]], dict]:
    key_dir = fixture_dir / "key"
    with open(key_dir / "leaked_pairs.json", "r", encoding="utf-8") as handle:
        raw_pairs = json.load(handle)
    with open(key_dir / "meta.json", "r", encoding="utf-8") as handle:
        meta = json.load(handle)
    pairs = {(item[0], item[1]) for item in raw_pairs}
    return pairs, meta


def recompute_clean_metric(
    fixture_dir: Path, leaked: Set[Tuple[str, str]]
) -> Tuple[Optional[float], int]:
    inputs_dir = fixture_dir / "inputs"

    labels_by_pair = {}
    with open(inputs_dir / "benchmark.csv", "r", encoding="utf-8") as handle:
        for row in csv.DictReader(handle):
            pair = normalize_pair(row["peptide"], row["mhc_allele"])
            if pair is None:
                continue
            labels_by_pair[pair] = 1 if row["true_label"].strip() == "binder" else 0

    clean_labels: List[int] = []
    clean_scores: List[float] = []
    with open(inputs_dir / "predictions.csv", "r", encoding="utf-8") as handle:
        for row in csv.DictReader(handle):
            pair = normalize_pair(row["peptide"], row["mhc_allele"])
            if pair is None or pair in leaked or pair not in labels_by_pair:
                continue
            clean_labels.append(labels_by_pair[pair])
            clean_scores.append(float(row["pred_score"]))

    return auroc(clean_labels, clean_scores), len(clean_labels)


def collect_submitted_pairs(
    submission: dict, submission_path: Path
) -> Tuple[Set[Tuple[str, str]], int, List[str]]:
    flags: List[str] = []
    raw_pairs = submission.get("leaked_pairs")

    if raw_pairs is None:
        side_file = submission_path.parent / "leaked_pairs.csv"
        if side_file.exists():
            raw_pairs = []
            with open(side_file, "r", encoding="utf-8") as handle:
                reader = csv.reader(handle)
                for index, row in enumerate(reader):
                    if index == 0 and row and not row[0].strip().isalpha():
                        continue
                    if len(row) >= 2:
                        raw_pairs.append([row[0], row[1]])

    if not isinstance(raw_pairs, list):
        return set(), 0, ["leaked_pairs missing or not a list"]

    normalized: Set[Tuple[str, str]] = set()
    malformed = 0
    for item in raw_pairs:
        if isinstance(item, (list, tuple)) and len(item) >= 2:
            pair = normalize_pair(item[0], item[1])
        elif isinstance(item, dict):
            pair = normalize_pair(item.get("peptide"), item.get("mhc_allele"))
        else:
            pair = None
        if pair is None:
            malformed += 1
            continue
        normalized.add(pair)

    if malformed:
        flags.append("{} malformed leaked_pairs entries ignored".format(malformed))
    return normalized, len(raw_pairs), flags


def grade(submission_path: Path, fixture_dir: Path) -> GraderResult:
    submission = load_submission(submission_path)
    if submission is None:
        return unreadable_submission(CHECK_NAMES, GATE)

    key_pairs, meta = read_key(fixture_dir)
    flags: List[str] = []

    contaminated = as_bool(submission.get("contaminated"))
    detection_ok = contaminated is True
    checks = [
        Check(
            "detection",
            detection_ok,
            "contaminated={}".format(submission.get("contaminated")),
        )
    ]

    submitted_pairs, raw_count, pair_flags = collect_submitted_pairs(
        submission, submission_path
    )
    flags.extend(pair_flags)

    overlap = jaccard(submitted_pairs, key_pairs)
    declared_count = submission.get("leaked_pair_count")
    count_agrees = declared_count == len(submitted_pairs)
    if not count_agrees and isinstance(declared_count, int):
        flags.append(
            "leaked_pair_count={} disagrees with {} submitted pairs".format(
                declared_count, len(submitted_pairs)
            )
        )
    checks.append(
        Check(
            "leaked_set_match",
            overlap >= JACCARD_THRESHOLD and count_agrees,
            "jaccard={:.4f} submitted={} key={} count_agrees={}".format(
                overlap, len(submitted_pairs), len(key_pairs), count_agrees
            ),
        )
    )

    checks.append(
        Check(
            "mechanism",
            submission.get("mechanism") == meta["mechanism"],
            "mechanism={}".format(submission.get("mechanism")),
        )
    )

    trusted = as_bool(submission.get("trusted_metric"))
    checks.append(
        Check(
            "metric_honesty",
            trusted is False,
            "trusted_metric={}".format(submission.get("trusted_metric")),
        )
    )

    recomputed_value, clean_count = recompute_clean_metric(fixture_dir, key_pairs)
    corrected = submission.get("corrected_metric")
    corrected = corrected if isinstance(corrected, dict) else {}
    claimed_value = corrected.get("value")
    claimed_clean = corrected.get("n_clean_pairs")
    claimed_uncomputable = corrected.get("uncomputable") is True or claimed_value is None

    if recomputed_value is None:
        recompute_ok = claimed_uncomputable and claimed_clean == clean_count
        detail = "no clean pairs remain; expected an uncomputable report"
        if not claimed_uncomputable and claimed_value is not None:
            flags.append("fabricated a metric with {} clean pairs".format(clean_count))
    else:
        recompute_ok = (
            not claimed_uncomputable
            and close(claimed_value, recomputed_value, METRIC_TOLERANCE)
            and claimed_clean == clean_count
        )
        detail = "recomputed={:.6f} claimed={} n_clean expected={} got={}".format(
            recomputed_value, claimed_value, clean_count, claimed_clean
        )
    checks.append(Check("recomputation", recompute_ok, detail))

    return build_result(
        checks,
        gates=[GATE],
        required_to_pass=CHECK_NAMES,
        flags=flags,
    )


def main() -> None:
    parser = argparse.ArgumentParser(description="Grade a contam-sentinel submission.")
    parser.add_argument("--submission", required=True, type=Path)
    parser.add_argument("--fixture", default="canonical")
    args = parser.parse_args()
    fixture_dir = ENV_DIR / "fixtures" / args.fixture
    print(grade(args.submission, fixture_dir).to_json())


if __name__ == "__main__":
    main()
