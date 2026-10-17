from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path
from typing import Dict, List

ENV_DIR = Path(__file__).resolve().parent
ROOT = ENV_DIR.parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from envs.allele_normalizer.normalize import canonicalize  # noqa: E402
from contract import (  # noqa: E402
    Check,
    GraderResult,
    as_int,
    build_result,
    close,
    jaccard,
    load_submission,
    unreadable_submission,
)

CHECK_NAMES = [
    "join_correctness",
    "payload_integrity",
    "counts",
    "mapping_validity",
    "no_silent_drops",
]
GATE = "join_correctness"
JACCARD_THRESHOLD = 0.99
PAYLOAD_TOLERANCE = 1e-6


def allele_id(name) -> str | None:
    if not isinstance(name, str):
        return None
    return canonicalize(name.strip(), {})


def grade(submission_path: Path, fixture_dir: Path) -> GraderResult:
    submission = load_submission(submission_path)
    if submission is None:
        return unreadable_submission(CHECK_NAMES, GATE)

    with open(fixture_dir / "key" / "meta.json", "r", encoding="utf-8") as handle:
        key = json.load(handle)

    flags: List[str] = []
    checks = []

    merged = submission.get("merged")
    merged = merged if isinstance(merged, list) else []
    merged_by_allele: Dict[str, dict] = {}
    for row in merged:
        if isinstance(row, dict) and isinstance(row.get("canonical_allele"), str):
            identified = allele_id(row["canonical_allele"])
            if identified:
                merged_by_allele[identified] = row

    expected = {allele_id(item) for item in key["expected_joined"]}
    overlap = jaccard(set(merged_by_allele), expected)
    checks.append(
        Check(
            "join_correctness",
            overlap >= JACCARD_THRESHOLD,
            "jaccard={:.4f} submitted={} expected={}".format(
                overlap, len(merged_by_allele), len(expected)
            ),
        )
    )

    payload_a = {allele_id(k): v for k, v in key["payload_a"].items()}
    payload_b = {allele_id(k): v for k, v in key["payload_b"].items()}
    mismatched = 0
    for allele in expected:
        row = merged_by_allele.get(allele)
        if row is None:
            mismatched += 1
            continue
        if not close(row.get("a_payload"), payload_a[allele], PAYLOAD_TOLERANCE):
            mismatched += 1
        elif not close(row.get("b_payload"), payload_b[allele], PAYLOAD_TOLERANCE):
            mismatched += 1
    if mismatched:
        flags.append("{} joined rows have wrong or swapped payloads".format(mismatched))
    checks.append(
        Check(
            "payload_integrity",
            mismatched == 0,
            "{} of {} rows mispaired".format(mismatched, len(expected)),
        )
    )

    joined_ok = as_int(submission.get("joined_count")) == key["joined_count"]
    dropped_ok = (
        as_int(submission.get("naive_dropped_count")) == key["naive_dropped_count"]
    )
    checks.append(
        Check(
            "counts",
            joined_ok and dropped_ok,
            "joined expected={} got={} | naive_dropped expected={} got={}".format(
                key["joined_count"],
                submission.get("joined_count"),
                key["naive_dropped_count"],
                submission.get("naive_dropped_count"),
            ),
        )
    )

    mapping = submission.get("mapping")
    mapping = mapping if isinstance(mapping, dict) else {}
    mapping = {str(k).strip(): v for k, v in mapping.items()}
    key_mapping = {k.strip(): v for k, v in key["mapping"].items()}
    missing = [raw for raw in key_mapping if raw not in mapping]
    wrong = [
        raw
        for raw, canonical in key_mapping.items()
        if raw in mapping and allele_id(mapping[raw]) != allele_id(canonical)
    ]
    checks.append(
        Check(
            "mapping_validity",
            not missing and not wrong,
            "{} unmapped, {} incorrect of {}".format(
                len(missing), len(wrong), len(key_mapping)
            ),
        )
    )

    declared_joined = as_int(submission.get("joined_count")) or 0
    conservation_a = declared_joined + len(key["only_a"]) == key["source_a_rows"]
    conservation_b = declared_joined + len(key["only_b"]) == key["source_b_rows"]
    checks.append(
        Check(
            "no_silent_drops",
            conservation_a and conservation_b,
            "a: {}+{}=={} b: {}+{}=={}".format(
                declared_joined,
                len(key["only_a"]),
                key["source_a_rows"],
                declared_joined,
                len(key["only_b"]),
                key["source_b_rows"],
            ),
        )
    )

    if declared_joined > key["joined_count"]:
        flags.append("claimed more joins than exist")

    return build_result(
        checks,
        gates=[GATE],
        required_to_pass=["join_correctness", "counts", "payload_integrity", "mapping_validity"],
        flags=flags,
    )


def main() -> None:
    parser = argparse.ArgumentParser(description="Grade an allele-normalizer submission.")
    parser.add_argument("--submission", required=True, type=Path)
    parser.add_argument("--fixture", default="canonical")
    args = parser.parse_args()
    print(grade(args.submission, ENV_DIR / "fixtures" / args.fixture).to_json())


if __name__ == "__main__":
    main()
