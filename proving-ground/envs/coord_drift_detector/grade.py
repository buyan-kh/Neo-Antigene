from __future__ import annotations

import argparse
import csv
import json
import sys
from pathlib import Path
from typing import Dict, List

ENV_DIR = Path(__file__).resolve().parent
ROOT = ENV_DIR.parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from contract import (  # noqa: E402
    Check,
    GraderResult,
    as_bool,
    as_int,
    build_result,
    load_submission,
    unreadable_submission,
)
from envs.coord_drift_detector.translate import (  # noqa: E402
    downstream_peptide,
    read_fasta,
)

CHECK_NAMES = [
    "detection",
    "affected_set",
    "corrected_peptides",
    "offsets",
    "non_triviality",
]
GATE = "detection"


def parse_vcf(path: Path) -> List[dict]:
    variants = []
    with open(path, "r", encoding="utf-8") as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            info = dict(
                item.split("=", 1) for item in fields[7].split(";") if "=" in item
            )
            variants.append(
                {
                    "variant_id": fields[2],
                    "ref": fields[3],
                    "alt": fields[4],
                    "transcript": info["TRANSCRIPT"],
                    "cds_pos": int(info["CDSPOS"]) - 1,
                }
            )
    return variants


def parse_cds_spans(path: Path) -> Dict[str, tuple]:
    origin = 0
    spans = {}
    with open(path, "r", encoding="utf-8") as handle:
        for line in handle:
            if line.startswith("##transcript-origin"):
                origin = int(line.split()[-1])
                continue
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9 or fields[2] != "CDS":
                continue
            info = dict(
                item.split("=", 1) for item in fields[8].split(";") if "=" in item
            )
            parent = info.get("Parent", "")
            transcript_id = parent.split(":")[-1] if parent else fields[0]
            phase = int(fields[7]) if fields[7].isdigit() else 0
            spans[transcript_id] = (int(fields[3]) - origin, int(fields[4]) - origin, phase)
    return spans


def recompute(fixture_dir: Path) -> Dict[str, dict]:
    inputs = fixture_dir / "inputs"
    transcripts = read_fasta(inputs / "transcripts.fasta")
    spans = parse_cds_spans(inputs / "annotation.gff3")
    result = {}
    for variant in parse_vcf(inputs / "variants.vcf"):
        start, end, phase = spans[variant["transcript"]]
        annotated = transcripts[variant["transcript"]][start - 1 : end]
        cds = annotated[phase:]
        cds = cds[: len(cds) - (len(cds) % 3)]
        peptide, offset = downstream_peptide(
            cds, variant["cds_pos"] - phase, variant["ref"], variant["alt"]
        )
        result[variant["variant_id"]] = {
            "downstream_peptide": peptide,
            "start_offset": offset,
        }
    return result


def read_reported(fixture_dir: Path) -> Dict[str, dict]:
    path = fixture_dir / "inputs" / "reported_neoorfs.csv"
    reported = {}
    with open(path, "r", encoding="utf-8") as handle:
        for row in csv.DictReader(handle):
            reported[row["variant_id"]] = {
                "downstream_peptide": row["downstream_peptide"],
                "start_offset": int(row["start_offset"]),
            }
    return reported


def grade(submission_path: Path, fixture_dir: Path) -> GraderResult:
    submission = load_submission(submission_path)
    if submission is None:
        return unreadable_submission(CHECK_NAMES, GATE)

    with open(fixture_dir / "key" / "meta.json", "r", encoding="utf-8") as handle:
        key = json.load(handle)

    truth = recompute(fixture_dir)
    for variant_id, expected in key["corrected"].items():
        if truth[variant_id] != expected:
            return build_result(
                [Check(name, False, "key disagreed with recomputation") for name in CHECK_NAMES],
                gates=[GATE],
                required_to_pass=CHECK_NAMES,
                flags=["committed key is not self-consistent"],
            )

    flags: List[str] = []
    checks = []

    detected = as_bool(submission.get("drift_detected")) is True
    mechanism_ok = submission.get("mechanism") == key["mechanism"]
    checks.append(
        Check(
            "detection",
            detected and mechanism_ok,
            "drift_detected={} mechanism={}".format(
                submission.get("drift_detected"), submission.get("mechanism")
            ),
        )
    )

    affected = submission.get("affected_variants")
    affected_set = {str(item) for item in affected} if isinstance(affected, list) else set()
    expected_affected = set(key["affected_variants"])
    checks.append(
        Check(
            "affected_set",
            affected_set == expected_affected,
            "expected {} got {}".format(len(expected_affected), len(affected_set)),
        )
    )

    corrected = submission.get("corrected")
    corrected_rows = {}
    if isinstance(corrected, list):
        for row in corrected:
            if isinstance(row, dict) and isinstance(row.get("variant_id"), str):
                corrected_rows[row["variant_id"]] = row

    peptide_mismatch = 0
    offset_mismatch = 0
    for variant_id in expected_affected:
        row = corrected_rows.get(variant_id)
        expected = truth[variant_id]
        if row is None:
            peptide_mismatch += 1
            offset_mismatch += 1
            continue
        if row.get("downstream_peptide") != expected["downstream_peptide"]:
            peptide_mismatch += 1
        if as_int(row.get("start_offset")) != expected["start_offset"]:
            offset_mismatch += 1

    checks.append(
        Check(
            "corrected_peptides",
            peptide_mismatch == 0 and bool(expected_affected),
            "{} of {} peptides wrong".format(peptide_mismatch, len(expected_affected)),
        )
    )
    checks.append(
        Check(
            "offsets",
            offset_mismatch == 0 and bool(expected_affected),
            "{} of {} offsets wrong".format(offset_mismatch, len(expected_affected)),
        )
    )

    reported = read_reported(fixture_dir)
    echoed = 0
    for variant_id in expected_affected:
        row = corrected_rows.get(variant_id)
        if row and row.get("downstream_peptide") == reported[variant_id]["downstream_peptide"]:
            echoed += 1
    if echoed:
        flags.append("{} corrected peptides simply echo the buggy output".format(echoed))
    checks.append(
        Check(
            "non_triviality",
            echoed == 0 and bool(corrected_rows),
            "{} of {} echoed the existing wrong output".format(
                echoed, len(expected_affected)
            ),
        )
    )

    return build_result(
        checks,
        gates=[GATE],
        required_to_pass=["detection", "affected_set", "corrected_peptides"],
        flags=flags,
    )


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Grade a coord-drift-detector submission."
    )
    parser.add_argument("--submission", required=True, type=Path)
    parser.add_argument("--fixture", default="canonical")
    args = parser.parse_args()
    print(grade(args.submission, ENV_DIR / "fixtures" / args.fixture).to_json())


if __name__ == "__main__":
    main()
