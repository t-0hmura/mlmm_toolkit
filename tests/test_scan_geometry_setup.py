"""The 1D scan reaches geometry setup with the YAML geom.coord_type (default cart)."""

from __future__ import annotations

from pathlib import Path

import pytest
from click.testing import CliRunner

from mlmm.cli import cli as root_cli


class _StopAtGeometry(RuntimeError):
    pass


def _write_layered_structure(path: Path) -> None:
    # Atoms 1-4 are ML (B = 0) and atoms 5-6 are movable MM (B = 10).
    lines = [
        f"HETATM{i:5d}  C{i:<2d} LIG A   1    {float(i - 1):8.3f}   0.000   0.000"
        f"  1.00{0.0 if i <= 4 else 10.0:6.2f}           C\n"
        for i in range(1, 7)
    ]
    path.write_text("".join(lines) + "END\n", encoding="utf-8")


def _invoke_scan(tmp_path: Path, monkeypatch, *extra: str):
    from mlmm.workflows import scan

    seen: dict[str, str] = {}

    def _fake_geom_loader(_path, *, coord_type, **_kwargs):
        seen["coord_type"] = coord_type
        raise _StopAtGeometry

    monkeypatch.setattr(scan, "geom_loader", _fake_geom_loader)
    structure = tmp_path / "system.pdb"
    _write_layered_structure(structure)
    parm = tmp_path / "system.parm7"
    parm.write_text("not reached\n", encoding="utf-8")
    result = CliRunner().invoke(
        root_cli,
        [
            "scan", "-i", str(structure), "--parm", str(parm), "-q", "0", "-m", "1",
            "--scan-lists", "[(1,2,1.2)]", "--out-dir", str(tmp_path / "out"), *extra,
        ],
    )
    return result, seen


@pytest.mark.parametrize(
    ("yaml_text", "expected"),
    [(None, "cart"), ("geom:\n  coord_type: dlc\n", "dlc"), ("geom:\n  coord_type: cart\n", "cart")],
)
def test_scan_geometry_uses_yaml_coord_type(
    tmp_path: Path, monkeypatch, yaml_text, expected
) -> None:
    extra: list[str] = []
    if yaml_text is not None:
        config = tmp_path / "config.yaml"
        config.write_text(yaml_text, encoding="utf-8")
        extra = ["--config", str(config)]

    result, seen = _invoke_scan(tmp_path, monkeypatch, *extra)

    assert seen.get("coord_type") == expected, result.output
