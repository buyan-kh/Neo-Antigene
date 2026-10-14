from __future__ import annotations

import argparse
import csv
import json
import sys
from pathlib import Path
from typing import List, Set

ENV_DIR = Path(__file__).resolve().parent
ROOT = ENV_DIR.parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from contract import (  # noqa: E402
    Check,
    GraderResult,
    build_result,
    load_submission,
    unreadable_submission,
)
from envs.pvalue_police.stats_core import (  # noqa: E402
    benjamini_hochberg,
    bonferroni,
    p_values,
)

CHECK_NAMES = [
    "detection",
    "localization",
    "method",
    "corrected_set",
    "non_triviality",
]
GATE = "detection"
ACCEPTED_METHODS = {"benjamini_hochberg", "bonferroni"}


def read_reported(fixture_dir: Path) -> Set[str]:
    path = fixture_dir / "inputs" / "results" / "reported_significant.csv"
    with open(path, "r", encoding="utf-8") as handle:
        return {row["feature_id"] for row in csv.DictReader(handle)}


def grade(submission_path: Path, fixture_dir: Path) -> GraderResult:
    submission = load_submission(submission_path)
    if submission is None:
        return unreadable_submission(CHECK_NAMES, GATE)

    with open(fixture_dir / "key" / "meta.json", "r", encoding="utf-8") as handle:
        key = json.load(handle)

    screen_ids = set()
    with open(fixture_dir / "inputs" / "raw" / "sample_sheet.csv", "r", encoding="utf-8") as handle:
        for row in csv.DictReader(handle):
            if row["role"].strip() == "screen":
                screen_ids.add(row["feature_id"])
    observed_p = p_values(fixture_dir / "inputs" / "raw" / "run_2024-11-02_plate_reader.csv")
    observed_p = {feature: value for feature, value in observed_p.items() if feature in screen_ids}
    recomputed = {
        "benjamini_hochberg": benjamini_hochberg(observed_p, key["q_value"]),
        "bonferroni": bonferroni(observed_p, key["alpha"]),
    }

    flags: List[str] = []
    checks = []

    detection_ok = submission.get("bug_class") == key["bug_class"]
    checks.append(
        Check("detection", detection_ok, "bug_class={}".format(submission.get("bug_class")))
    )

    span = key["bug_span"]
    location = submission.get("location")
    location = location if isinstance(location, dict) else {}
    same_file = location.get("file") == span["file"]
    try:
        start = int(location.get("line_start"))
        end = int(location.get("line_end"))
        overlaps = start <= span["line_end"] and end >= span["line_start"]
    except (TypeError, ValueError):
        overlaps = False
    checks.append(
        Check(
            "localization",
            bool(same_file and overlaps),
            "expected {} lines {}-{} got {}".format(
                span["file"], span["line_start"], span["line_end"], location
            ),
        )
    )

    method = submission.get("corrected_method")
    method_ok = method in ACCEPTED_METHODS
    checks.append(Check("method", method_ok, "corrected_method={}".format(method)))

    submitted_set = submission.get("corrected_significant")
    submitted_set = (
        {str(item) for item in submitted_set}
        if isinstance(submitted_set, list)
        else set()
    )

    if method_ok:
        expected = recomputed[method]
        if set(key[method]) != expected:
            flags.append("committed key disagreed with recomputation for " + method)
        set_ok = submitted_set == expected
        detail = "expected {} features got {}".format(len(expected), len(submitted_set))
    else:
        set_ok = False
        detail = "no accepted correction method declared"
    checks.append(Check("corrected_set", set_ok, detail))

    reported = read_reported(fixture_dir)
    non_trivial = submitted_set != reported and bool(submitted_set)
    if submitted_set == reported:
        flags.append("echoed the tool's uncorrected output as the corrected set")
    checks.append(
        Check(
            "non_triviality",
            non_trivial,
            "reported={} submitted={}".format(len(reported), len(submitted_set)),
        )
    )

    return build_result(
        checks,
        gates=[GATE],
        required_to_pass=["detection", "localization", "corrected_set"],
        flags=flags,
    )


def main() -> None:
    parser = argparse.ArgumentParser(description="Grade a pvalue-police submission.")
    parser.add_argument("--submission", required=True, type=Path)
    parser.add_argument("--fixture", default="canonical")
    args = parser.parse_args()
    print(grade(args.submission, ENV_DIR / "fixtures" / args.fixture).to_json())


if __name__ == "__main__":
    main()
