"""ML membership precedence is explicit model > indices > B-factor fallback."""

from __future__ import annotations

from pathlib import Path

from mlmm.core.utils import (
    read_bfactors_from_pdb,
    resolve_ml_layer_assignment,
)


def _line(serial: int, name: str, bfactor: float) -> str:
    return (
        f"ATOM  {serial:5d} {name:^4s} LIG A   1    "
        f"{float(serial):8.3f}{0.0:8.3f}{0.0:8.3f}"
        f"  1.00{bfactor:6.2f}           C\n"
    )


def test_explicit_model_wins_over_valid_bfactor_membership(tmp_path: Path) -> None:
    source = tmp_path / "full.pdb"
    model = tmp_path / "explicit.pdb"
    lines = [_line(1, "C1", 0.0), _line(2, "C2", 10.0), _line(3, "C3", 20.0)]
    source.write_text("".join(lines) + "END\n", encoding="utf-8")
    model.write_text(lines[2] + "END\n", encoding="utf-8")
    cfg: dict = {}

    resolved, layer_info = resolve_ml_layer_assignment(
        source_path=source,
        out_dir_path=tmp_path / "out",
        model_pdb=model,
        model_indices=None,
        detect_layer=True,
        hess_cutoff=None,
        movable_cutoff=None,
        calc_cfg=cfg,
    )

    assert resolved == model
    assert layer_info is not None
    assert cfg["model_pdb"] == str(model)
    assert cfg["use_bfactor_layers"] is True


def test_indices_win_when_model_is_absent(tmp_path: Path) -> None:
    source = tmp_path / "full.pdb"
    lines = [_line(1, "C1", 0.0), _line(2, "C2", 10.0), _line(3, "C3", 20.0)]
    source.write_text("".join(lines) + "END\n", encoding="utf-8")
    cfg: dict = {}

    resolved, _ = resolve_ml_layer_assignment(
        source_path=source,
        out_dir_path=tmp_path / "out",
        model_pdb=None,
        model_indices=[1],
        detect_layer=True,
        hess_cutoff=None,
        movable_cutoff=None,
        calc_cfg=cfg,
    )

    assert len(read_bfactors_from_pdb(resolved)) == 1
    assert "C2" in resolved.read_text(encoding="utf-8")


def test_bfactor_is_used_only_without_explicit_membership(tmp_path: Path) -> None:
    source = tmp_path / "full.pdb"
    lines = [_line(1, "C1", 0.0), _line(2, "C2", 10.0), _line(3, "C3", 20.0)]
    source.write_text("".join(lines) + "END\n", encoding="utf-8")
    cfg: dict = {}

    resolved, layer_info = resolve_ml_layer_assignment(
        source_path=source,
        out_dir_path=tmp_path / "out",
        model_pdb=None,
        model_indices=None,
        detect_layer=True,
        hess_cutoff=None,
        movable_cutoff=None,
        calc_cfg=cfg,
    )

    assert layer_info is not None
    assert layer_info["ml_indices"] == [0]
    assert len(read_bfactors_from_pdb(resolved)) == 1
