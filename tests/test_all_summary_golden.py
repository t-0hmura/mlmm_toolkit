"""Golden root ``summary.json`` for every ``mlmm all`` exit that writes it.

Each case runs the real ``all`` callback with stubbed children (see
``_all_stub_harness.py``) and compares the exit code, the child commands that
ran and the final root ``summary.json`` (every key and value) with
``tests/data/all_summary_golden/<case>.json``. Values that change from run to
run (run id, temporary directory, package version, host environment) are
replaced by placeholders first; floats are compared to 1e-9 relative.

After an intended output change, regenerate with
``MLMM_UPDATE_SUMMARY_GOLDEN=1 pytest tests/test_all_summary_golden.py`` and
review the JSON diff.
"""

from __future__ import annotations

import json
import math
import os
from pathlib import Path

import pytest

from _all_stub_harness import TS_NOT_CONVERGED, Behaviour, run_all

GOLDEN_DIR = Path(__file__).resolve().parent / "data" / "all_summary_golden"
PLACEHOLDERS = {
    "run_id": "<run_id>",
    "environment": "<environment>",
    "mlmm_toolkit_version": "<version>",
}

_BRIDGE = {
    "index": 1, "tag": "seg_001_bridge", "kind": "bridge",
    "barrier_kcal": 0.5, "delta_kcal": 0.1, "bond_changes": "", "converged": True,
}

# case: (input frames, extra CLI options, child behaviour)
CASES = {
    # MEP route, single-pass path-opt.
    "mep_path_opt": ((0.0, 2.0), (), None),
    "mep_path_opt_tsopt": ((0.0, 2.0), ("--tsopt",), None),
    "mep_path_opt_tsopt_thermo_dft": ((0.0, 2.0), ("--tsopt", "--thermo", "--dft"), None),
    "mep_stopped_before_irc": ((0.0, 2.0), ("--tsopt",), Behaviour(tsopt=dict(TS_NOT_CONVERGED))),
    "mep_endpoint_failure": ((0.0, 2.0), ("--tsopt",), Behaviour(opt={"P": "error"})),
    # MEP route, recursive path-search (--refine-path).
    "refine_path": ((0.0, 2.0), ("--refine-path",), None),
    "refine_path_tsopt": ((0.0, 2.0), ("--refine-path", "--tsopt"), None),
    "early_no_segments": ((0.0, 2.0), ("--refine-path", "--tsopt"),
                          Behaviour(path_search_segments=[])),
    "early_no_reactive_segment": ((0.0, 2.0), ("--refine-path", "--tsopt"),
                                  Behaviour(path_search_segments=[dict(_BRIDGE)])),
    # TS-only route (one input with --tsopt).
    "tsonly": ((1.0,), ("--tsopt",), None),
    "tsonly_thermo_dft": ((1.0,), ("--tsopt", "--thermo", "--dft"), None),
    "tsonly_stopped_before_irc": ((1.0,), ("--tsopt",), Behaviour(tsopt=dict(TS_NOT_CONVERGED))),
    "tsonly_endpoint_failure": ((1.0,), ("--tsopt",), Behaviour(opt={"P": "error"})),
}


def _normalize(value, roots):
    if isinstance(value, dict):
        return {
            key: PLACEHOLDERS[key] if key in PLACEHOLDERS else _normalize(item, roots)
            for key, item in value.items()
        }
    if isinstance(value, list):
        return [_normalize(item, roots) for item in value]
    if isinstance(value, str):
        for root in roots:
            value = value.replace(root, "<tmp>")
    return value


def _diff(expected, actual, path="$"):
    if isinstance(expected, dict) and isinstance(actual, dict):
        lines = [f"{path}.{key}: missing" for key in sorted(set(expected) - set(actual))]
        lines += [f"{path}.{key}: unexpected" for key in sorted(set(actual) - set(expected))]
        for key in sorted(set(expected) & set(actual)):
            lines += _diff(expected[key], actual[key], f"{path}.{key}")
        return lines
    if isinstance(expected, list) and isinstance(actual, list):
        if len(expected) != len(actual):
            return [f"{path}: {len(expected)} items expected, {len(actual)} found"]
        lines = []
        for index, (left, right) in enumerate(zip(expected, actual)):
            lines += _diff(left, right, f"{path}[{index}]")
        return lines
    if isinstance(expected, float) and isinstance(actual, float):
        if math.isclose(expected, actual, rel_tol=1e-9, abs_tol=1e-12):
            return []
        return [f"{path}: {expected!r} != {actual!r}"]
    if type(expected) is not type(actual) or expected != actual:
        return [f"{path}: {expected!r} != {actual!r}"]
    return []


def snapshot(run, tmp_path):
    roots = sorted({str(tmp_path), str(Path(tmp_path).resolve())}, key=len, reverse=True)
    return {
        "exit_code": run.result.exit_code,
        "children": run.child_names(),
        "summary": _normalize(run.summary, roots),
    }


@pytest.mark.parametrize("case", list(CASES))
def test_all_root_summary_matches_golden(tmp_path, monkeypatch, case):
    frames, extra, behaviour = CASES[case]
    run = run_all(tmp_path, monkeypatch, frames=frames, extra=extra, behaviour=behaviour)
    assert (run.out / "summary.json").is_file(), run.result.output[-3000:]
    actual = snapshot(run, tmp_path)
    golden = GOLDEN_DIR / f"{case}.json"
    if os.environ.get("MLMM_UPDATE_SUMMARY_GOLDEN") == "1":
        GOLDEN_DIR.mkdir(parents=True, exist_ok=True)
        golden.write_text(
            json.dumps(actual, indent=1, ensure_ascii=False, sort_keys=True) + "\n",
            encoding="utf-8",
        )
    expected = json.loads(golden.read_text(encoding="utf-8"))
    differences = _diff(expected, actual)
    assert not differences, "\n".join(differences[:60])
