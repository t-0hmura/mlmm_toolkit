"""Regression tests for the executable documentation smoke sanitizer."""

from __future__ import annotations

import importlib.util
from pathlib import Path


def _load_smoke_module():
    path = Path(__file__).parents[1] / ".github" / "scripts" / "smoke_docs_commands.py"
    spec = importlib.util.spec_from_file_location("mlmm_smoke_docs_commands", path)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_dynamic_smoke_preserves_all_input_mode(tmp_path: Path) -> None:
    smoke = _load_smoke_module()
    fixture = smoke._prepare_fixture_files(tmp_path)

    endpoint = smoke._sanitize_all_args(
        ["all", "-i", "R.pdb", "P.pdb", "-c", "LIG"], fixture
    )
    staged_scan = smoke._sanitize_all_args(
        [
            "all",
            "-i",
            "R.pdb",
            "-c",
            "LIG",
            "--scan-lists",
            "[(1, 2, 1.5)]",
        ],
        fixture,
    )

    endpoint_i = endpoint.index("-i")
    assert endpoint[endpoint_i + 1 : endpoint_i + 3] == [
        str(fixture["r_pdb"]),
        str(fixture["p_pdb"]),
    ]
    scan_i = staged_scan.index("-i")
    assert staged_scan[scan_i + 1] == str(fixture["r_pdb"])
    assert staged_scan[scan_i + 2] == "-c"
    assert staged_scan[staged_scan.index("--scan-lists") + 1] == "[(1,2,1.5)]"
