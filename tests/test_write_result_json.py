"""Unit tests for mlmm.core.utils.write_result_json."""

from __future__ import annotations

import json
import tempfile
from pathlib import Path

from mlmm.core.utils import (
    write_result_json,
    RESULT_JSON_SCHEMA_VERSION,
    RESULT_JSON_STATUS_VALUES,
)


def test_writes_result_and_summary_mirror() -> None:
    with tempfile.TemporaryDirectory() as d:
        path = write_result_json(
            Path(d),
            {"status": "success", "energy_hartree": -1.5},
            command="opt",
            elapsed_seconds=1.234,
        )
        assert path is not None
        assert (Path(d) / "result.json").exists()
        assert (Path(d) / "summary.json").exists()
        r = json.loads((Path(d) / "result.json").read_text())
        s = json.loads((Path(d) / "summary.json").read_text())
        assert r == s
        assert r["status"] == "success"
        assert r["schema_version"] == RESULT_JSON_SCHEMA_VERSION
        assert r["command"] == "opt"
        assert r["elapsed_seconds"] == 1.234


def test_summary_json_input_skips_mirror() -> None:
    with tempfile.TemporaryDirectory() as d:
        write_result_json(
            Path(d),
            {"status": "success"},
            command="all",
            filename="summary.json",
        )
        assert (Path(d) / "summary.json").exists()
        assert not (Path(d) / "result.json").exists()


def test_disable_summary_mirror() -> None:
    with tempfile.TemporaryDirectory() as d:
        write_result_json(
            Path(d),
            {"status": "success"},
            command="opt",
            also_write_summary_json=False,
        )
        assert (Path(d) / "result.json").exists()
        assert not (Path(d) / "summary.json").exists()


def test_status_enum_documented() -> None:
    # Catch typo regressions in the documented status values.
    assert RESULT_JSON_STATUS_VALUES == ("success", "partial", "error", "unknown")
