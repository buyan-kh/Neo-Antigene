from __future__ import annotations

import argparse
import hashlib
import importlib
import json
import shutil
import subprocess
import sys
import time
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, List, Optional

ROOT = Path(__file__).resolve().parent
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

ENV_MODULES = {
    "contam-sentinel": "envs.contam_sentinel",
    "power-gate": "envs.power_gate",
    "pvalue-police": "envs.pvalue_police",
    "allele-normalizer": "envs.allele_normalizer",
    "coord-drift-detector": "envs.coord_drift_detector",
}

OUTPUT_NOTE = (
    "\n\n---\n"
    "Operational note: write your answer to ./output/submission.json "
    "(relative to this working directory), not to /output/. "
    "The directory ./output/ already exists. "
    "You may create and run scripts in this directory."
)


def env_dir(env_name: str) -> Path:
    module = importlib.import_module(ENV_MODULES[env_name])
    return Path(module.__file__).resolve().parent


def grade_submission(env_name: str, submission: Path, fixture_dir: Path):
    module = importlib.import_module(ENV_MODULES[env_name] + ".grade")
    return module.grade(submission, fixture_dir)


def stage_sandbox(env_name: str, fixture: str, destination: Path) -> None:
    source = env_dir(env_name)
    inputs = source / "fixtures" / fixture / "inputs"
    if destination.exists():
        shutil.rmtree(destination)
    (destination / "output").mkdir(parents=True, exist_ok=True)
    for item in sorted(inputs.iterdir()):
        if item.is_file():
            shutil.copy2(item, destination / item.name)
        else:
            shutil.copytree(item, destination / item.name)
    shutil.copy2(source / "task.md", destination / "TASK.md")


def assert_no_key(sandbox: Path) -> None:
    offenders = [
        str(path)
        for path in sandbox.rglob("*")
        if path.is_file()
        and (path.name == "meta.json" or "leaked_pairs" in path.name or path.parent.name == "key")
    ]
    if offenders:
        raise RuntimeError("key material leaked into sandbox: " + ", ".join(offenders))


def run_one(
    model: str,
    env_name: str,
    fixture: str,
    run_root: Path,
    timeout: int,
) -> Dict:
    safe_model = model.replace("/", "_")
    digest = hashlib.sha1("{}:{}".format(env_name, fixture).encode("utf-8")).hexdigest()[:8]
    sandbox = run_root / safe_model / "job-{}".format(digest)
    stage_sandbox(env_name, fixture, sandbox)
    assert_no_key(sandbox)

    prompt = (sandbox / "TASK.md").read_text(encoding="utf-8") + OUTPUT_NOTE

    command = [
        "cursor-agent",
        "--print",
        "--force",
        "--model",
        model,
        "--workspace",
        str(sandbox),
        prompt,
    ]

    started = time.time()
    timed_out = False
    try:
        completed = subprocess.run(
            command,
            cwd=str(sandbox),
            capture_output=True,
            text=True,
            timeout=timeout,
        )
        stdout, stderr, returncode = completed.stdout, completed.stderr, completed.returncode
    except subprocess.TimeoutExpired as exc:
        stdout = exc.stdout.decode("utf-8", "replace") if isinstance(exc.stdout, bytes) else (exc.stdout or "")
        stderr = exc.stderr.decode("utf-8", "replace") if isinstance(exc.stderr, bytes) else (exc.stderr or "")
        returncode = -1
        timed_out = True
    elapsed = time.time() - started

    (sandbox / "agent_stdout.txt").write_text(stdout or "", encoding="utf-8")
    (sandbox / "agent_stderr.txt").write_text(stderr or "", encoding="utf-8")

    submission = sandbox / "output" / "submission.json"
    record: Dict = {
        "model": model,
        "env": env_name,
        "fixture": fixture,
        "elapsed_seconds": round(elapsed, 1),
        "timed_out": timed_out,
        "agent_returncode": returncode,
        "submission_written": submission.exists(),
        "sandbox": str(sandbox),
    }

    if submission.exists():
        fixture_dir = env_dir(env_name) / "fixtures" / fixture
        try:
            result = grade_submission(env_name, submission, fixture_dir)
            payload = result.to_dict()
            (sandbox / "grader.json").write_text(
                json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8"
            )
            record["score"] = payload["score"]
            record["passed"] = payload["passed"]
            record["failed_checks"] = [
                check["name"] for check in payload["checks"] if not check["passed"]
            ]
            record["reward_hack_flags"] = payload["reward_hack_flags"]
        except Exception as exc:  # noqa: BLE001
            record["score"] = 0.0
            record["passed"] = False
            record["grader_error"] = "{}: {}".format(type(exc).__name__, exc)
            record["failed_checks"] = ["grader_error"]
            record["reward_hack_flags"] = []
    else:
        record["score"] = 0.0
        record["passed"] = False
        record["failed_checks"] = ["no_submission"]
        record["reward_hack_flags"] = []

    return record


def append_record(results_path: Path, record: Dict) -> None:
    with open(results_path, "a", encoding="utf-8") as handle:
        handle.write(json.dumps(record, sort_keys=True) + "\n")


def main() -> int:
    parser = argparse.ArgumentParser(description="Run models against environments.")
    parser.add_argument("--models", required=True)
    parser.add_argument(
        "--envs",
        required=True,
        help="comma separated env:fixture pairs, e.g. contam-sentinel:canonical",
    )
    parser.add_argument("--run-dir", default=None)
    parser.add_argument("--timeout", type=int, default=2400)
    args = parser.parse_args()

    models = [m.strip() for m in args.models.split(",") if m.strip()]
    targets = []
    for item in args.envs.split(","):
        item = item.strip()
        if not item:
            continue
        env_name, _, fixture = item.partition(":")
        targets.append((env_name, fixture or "canonical"))

    stamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")
    run_root = (Path(args.run_dir) if args.run_dir else ROOT / "runs" / stamp).resolve()
    run_root.mkdir(parents=True, exist_ok=True)
    results_path = run_root / "results.jsonl"

    print("run dir: {}".format(run_root), flush=True)
    for model in models:
        for env_name, fixture in targets:
            print(
                "\n>>> {} | {}:{}".format(model, env_name, fixture),
                flush=True,
            )
            record = run_one(model, env_name, fixture, run_root, args.timeout)
            append_record(results_path, record)
            print(
                "    score={} passed={} elapsed={}s failed={}".format(
                    record.get("score"),
                    record.get("passed"),
                    record.get("elapsed_seconds"),
                    ",".join(record.get("failed_checks") or []) or "-",
                ),
                flush=True,
            )
    print("\nresults: {}".format(results_path))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
