from __future__ import annotations

from pathlib import Path
from typing import Callable, List, Tuple

DEFAULT_RUNS = 100


def check_repeatable(
    grade: Callable[[Path, Path], "object"],
    submission_path: Path,
    fixture_dir: Path,
    runs: int = DEFAULT_RUNS,
) -> Tuple[bool, str]:
    first = grade(submission_path, fixture_dir).to_json()
    for index in range(1, runs):
        again = grade(submission_path, fixture_dir).to_json()
        if again != first:
            return False, "run {} differed from run 0".format(index)
    return True, first


def sweep_reference_submissions(
    grade: Callable[[Path, Path], "object"],
    reference_dir: Path,
    fixture_dir: Path,
    runs: int = DEFAULT_RUNS,
) -> List[Tuple[str, bool, str]]:
    results = []
    for submission_path in sorted(reference_dir.glob("*.json")):
        ok, detail = check_repeatable(grade, submission_path, fixture_dir, runs)
        results.append((submission_path.name, ok, detail))
    return results
