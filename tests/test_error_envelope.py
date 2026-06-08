"""Unit tests for the structured error envelope produced by _write_error_json."""

from __future__ import annotations

import json
import tempfile
from pathlib import Path

from mlmm.cli.decorators import _write_error_json


class _CustomOptError(RuntimeError):
    pass


def test_error_envelope_includes_class_chain() -> None:
    with tempfile.TemporaryDirectory() as d:
        _write_error_json(
            Path(d),
            "opt",
            _CustomOptError("optimization diverged"),
            "opt-stage",
            time_start=None,
        )
        r = json.loads((Path(d) / "result.json").read_text())
        assert r["status"] == "error"
        assert r["error"] == "optimization diverged"
        assert r["error_type"] == "_CustomOptError"
        chain = r["error_class_chain"]
        assert "RuntimeError" in chain
        assert "Exception" in chain
        assert "BaseException" in chain
        assert chain[0] == "_CustomOptError"
        assert r["error_module"] == _CustomOptError.__module__
        assert r["error_label"] == "opt-stage"


def test_error_envelope_mirrors_to_summary_json() -> None:
    with tempfile.TemporaryDirectory() as d:
        _write_error_json(
            Path(d),
            "tsopt",
            ValueError("malformed input"),
            "tsopt-stage",
        )
        r = json.loads((Path(d) / "result.json").read_text())
        s = json.loads((Path(d) / "summary.json").read_text())
        assert r == s
