"""Regression tests for summary_log payload handling and rendering."""

from __future__ import annotations

import sys
from pathlib import Path

import pytest

pytestmark = pytest.mark.skipif(
    sys.version_info < (3, 11),
    reason="mlmm requires Python >= 3.11",
)


def test_normalize_summary_payload_sets_defaults():
    from mlmm.summary_log import normalize_summary_payload

    payload = normalize_summary_payload({"pipeline_mode": "path-search"})
    assert payload["pipeline_mode"] == "path-search"
    assert payload["root_out_dir"] == "-"
    assert payload["path_module_dir"] == "-"
    assert payload["segments"] == []
    assert payload["energy_diagrams"] == []


def test_write_summary_log_accepts_empty_payload(tmp_path: Path):
    from mlmm.summary_log import write_summary_log

    dest = tmp_path / "summary.log"
    write_summary_log(dest, {})

    text = dest.read_text(encoding="utf-8")
    assert "mlmm summary.log" in text
    assert "missing payload keys replaced with defaults" in text


def test_write_summary_log_renders_segment_section(tmp_path: Path):
    from mlmm.summary_log import write_summary_log

    out_root = tmp_path / "run"
    out_root.mkdir(parents=True, exist_ok=True)

    payload = {
        "root_out_dir": str(out_root),
        "path_module_dir": "path_search",
        "pipeline_mode": "path-search",
        "segments": [
            {
                "index": 1,
                "tag": "seg_01",
                "kind": "seg",
                "barrier_kcal": 12.3,
                "delta_kcal": -1.2,
                "bond_changes": "Broken: C1-O1",
            }
        ],
        "energy_diagrams": [],
    }

    dest = out_root / "summary.log"
    write_summary_log(dest, payload)

    text = dest.read_text(encoding="utf-8")
    assert "[2] Segment-level MEP summary" in text
    assert "Segment 01 [seg]" in text
    assert "Broken: C1-O1" in text
