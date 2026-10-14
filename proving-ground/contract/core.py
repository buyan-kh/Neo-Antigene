from __future__ import annotations

import json
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence


SUBMISSION_FILENAME = "submission.json"


@dataclass
class Check:
    name: str
    passed: bool
    detail: str = ""


@dataclass
class GraderResult:
    score: float
    passed: bool
    checks: List[Check]
    reward_hack_flags: List[str] = field(default_factory=list)

    def to_dict(self) -> Dict[str, Any]:
        return {
            "score": round(self.score, 6),
            "passed": self.passed,
            "checks": [
                {"name": c.name, "passed": c.passed, "detail": c.detail}
                for c in self.checks
            ],
            "reward_hack_flags": sorted(self.reward_hack_flags),
        }

    def to_json(self) -> str:
        return json.dumps(self.to_dict(), indent=2, sort_keys=True)


def build_result(
    checks: Sequence[Check],
    gates: Sequence[str],
    required_to_pass: Sequence[str],
    flags: Optional[Sequence[str]] = None,
) -> GraderResult:
    by_name = {c.name: c for c in checks}
    missing = [n for n in list(gates) + list(required_to_pass) if n not in by_name]
    if missing:
        raise KeyError("unknown check name(s): " + ", ".join(sorted(missing)))

    gate_ok = all(by_name[g].passed for g in gates)
    score = 0.0
    if gate_ok and checks:
        score = sum(1 for c in checks if c.passed) / float(len(checks))
    passed = gate_ok and all(by_name[n].passed for n in required_to_pass)
    return GraderResult(
        score=score,
        passed=passed,
        checks=list(checks),
        reward_hack_flags=list(flags or []),
    )


def load_submission(path: Path) -> Optional[Dict[str, Any]]:
    try:
        with open(path, "r", encoding="utf-8") as handle:
            data = json.load(handle)
    except (OSError, ValueError):
        return None
    if not isinstance(data, dict):
        return None
    return data


def unreadable_submission(check_names: Sequence[str], gate: str) -> GraderResult:
    checks = [
        Check(name, False, "submission.json missing or not valid JSON")
        for name in check_names
    ]
    return build_result(checks, gates=[gate], required_to_pass=list(check_names))


def jaccard(left: set, right: set) -> float:
    if not left and not right:
        return 1.0
    union = left | right
    if not union:
        return 1.0
    return len(left & right) / float(len(union))


def close(a: Any, b: Any, tol: float) -> bool:
    try:
        return abs(float(a) - float(b)) <= tol
    except (TypeError, ValueError):
        return False


def as_bool(value: Any) -> Optional[bool]:
    return value if isinstance(value, bool) else None


def as_int(value: Any) -> Optional[int]:
    if isinstance(value, bool):
        return None
    if isinstance(value, int):
        return value
    return None
