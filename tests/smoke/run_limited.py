"""Check exit codes for smoke calculations with deliberately small cycle budgets."""

import json
import subprocess
import sys
from pathlib import Path

from mlmm.cli.completion import completion_code


def check_completion(code, payload):
    if payload.get("execution_status") != "completed" or "status" in payload:
        raise ValueError("Limited run did not finish execution with the public result schema")
    scientific = payload.get("scientific_status")
    if scientific not in {"success", "partial", "failed"}:
        raise ValueError("Limited run is missing its scientific outcome")
    expected = completion_code(payload)
    if code != expected:
        raise ValueError(f"CLI exit {code} disagrees with scientific_status={scientific!r}")
    return scientific


def main():
    module, command, *arguments = sys.argv[1:]
    out_dir = None
    for index, argument in enumerate(arguments):
        if argument in {"-o", "--out-dir"}:
            out_dir = Path(arguments[index + 1])
        elif argument.startswith("--out-dir="):
            out_dir = Path(argument.split("=", 1)[1])
    if out_dir is None:
        raise SystemExit("Limited smoke run requires an explicit --out-dir")
    aggregate = command in {"all", "path-search"}
    if not aggregate and "--out-json" not in arguments:
        arguments.append("--out-json")
    code = subprocess.run([sys.executable, "-m", module, command, *arguments]).returncode
    filename = "summary.json" if aggregate else "result.json"
    payload = json.loads((out_dir / filename).read_text())
    scientific = check_completion(code, payload)
    print(f"[smoke limited] cli_exit={code} execution_status=completed scientific_status={scientific}")


if __name__ == "__main__":
    main()

