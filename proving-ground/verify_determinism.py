from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from contract.determinism import DEFAULT_RUNS, sweep_reference_submissions  # noqa: E402
from envs.allele_normalizer import grade as allele_grade  # noqa: E402
from envs.contam_sentinel import grade as contam_grade  # noqa: E402
from envs.coord_drift_detector import grade as coord_grade  # noqa: E402
from envs.power_gate import grade as power_grade  # noqa: E402
from envs.pvalue_police import grade as pvalue_grade  # noqa: E402

ENVIRONMENTS = [
    ("contam-sentinel", contam_grade),
    ("power-gate", power_grade),
    ("pvalue-police", pvalue_grade),
    ("allele-normalizer", allele_grade),
    ("coord-drift-detector", coord_grade),
]


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Run every grader repeatedly and assert identical output."
    )
    parser.add_argument("--runs", type=int, default=DEFAULT_RUNS)
    args = parser.parse_args()

    failures = 0
    total = 0
    for name, module in ENVIRONMENTS:
        env_dir = Path(module.__file__).resolve().parent
        fixtures_root = env_dir / "fixtures"
        reference_root = env_dir / "reference"
        for fixture_dir in sorted(p for p in fixtures_root.iterdir() if p.is_dir()):
            reference_dir = reference_root / fixture_dir.name
            if not reference_dir.is_dir():
                continue
            results = sweep_reference_submissions(
                module.grade, reference_dir, fixture_dir, args.runs
            )
            for submission_name, ok, detail in results:
                total += 1
                expected_pass = submission_name.startswith("good")
                verdict = None
                if ok:
                    verdict = json.loads(detail)["passed"]
                    if verdict is not expected_pass:
                        ok = False
                        detail = "expected passed={} got passed={}".format(
                            expected_pass, verdict
                        )
                status = "ok " if ok else "FAIL"
                print(
                    "{} {}/{}/{} passed={}".format(
                        status, name, fixture_dir.name, submission_name, verdict
                    )
                )
                if not ok:
                    failures += 1
                    print("     {}".format(detail))

    print(
        "\n{} submissions x {} runs each: {}".format(
            total, args.runs, "all identical, good/bad verdicts correct" if failures == 0 else
            "{} failures".format(failures)
        )
    )
    return 1 if failures else 0


if __name__ == "__main__":
    raise SystemExit(main())
