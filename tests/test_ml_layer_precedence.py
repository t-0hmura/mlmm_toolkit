"""ML membership precedence is explicit model > indices > B-factor fallback."""

from __future__ import annotations

from pathlib import Path

import click
import numpy as np
import pytest

from mlmm.core.utils import (
    apply_layer_freeze_constraints,
    read_bfactors_from_pdb,
    resolve_ml_layer_assignment,
)
from mlmm.backends.mlmm_calc import MLMMCore, mlmm


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
        protected_inputs=(),
    )

    assert resolved == model
    assert layer_info is not None
    assert cfg["model_pdb"] == str(model)
    assert cfg["use_bfactor_layers"] is True
    assert layer_info["ml_indices"] == [2]
    assert layer_info["frozen_indices"] == []
    assert layer_info["unassigned_indices"] == [0]

    geom_cfg: dict = {}
    assert apply_layer_freeze_constraints(geom_cfg, cfg, layer_info) == []


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
        protected_inputs=(),
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
        protected_inputs=(),
    )

    assert layer_info is not None
    assert layer_info["ml_indices"] == [0]
    assert len(read_bfactors_from_pdb(resolved)) == 1


def test_bfactor_model_generation_preserves_protected_target(
    tmp_path: Path,
) -> None:
    source = tmp_path / "full.pdb"
    source.write_text(
        _line(1, "C1", 0.0) + _line(2, "C2", 10.0) + "END\n",
        encoding="utf-8",
    )
    out_dir = tmp_path / "out"
    out_dir.mkdir()
    protected = out_dir / "model_from_bfactor.pdb"
    original = b"protected input"
    protected.write_bytes(original)

    with pytest.raises(click.ClickException, match="collides"):
        resolve_ml_layer_assignment(
            source_path=source,
            out_dir_path=out_dir,
            model_pdb=None,
            model_indices=None,
            detect_layer=True,
            hess_cutoff=None,
            movable_cutoff=None,
            calc_cfg={},
            protected_inputs=(protected,),
        )

    assert protected.read_bytes() == original


def test_bfactor_model_generation_replaces_unprotected_prior_output(
    tmp_path: Path,
) -> None:
    source = tmp_path / "full.pdb"
    source.write_text(
        _line(1, "C1", 0.0) + _line(2, "C2", 10.0) + "END\n",
        encoding="utf-8",
    )
    out_dir = tmp_path / "out"
    out_dir.mkdir()
    prior = out_dir / "model_from_bfactor.pdb"
    prior.write_text("stale generated output\n", encoding="utf-8")

    resolved, layer_info = resolve_ml_layer_assignment(
        source_path=source,
        out_dir_path=out_dir,
        model_pdb=None,
        model_indices=None,
        detect_layer=True,
        hess_cutoff=None,
        movable_cutoff=None,
        calc_cfg={},
        protected_inputs=(),
    )

    assert resolved == prior
    assert layer_info is not None
    assert "C1" in prior.read_text(encoding="utf-8")


def test_distance_movable_cutoff_bounds_default_hessian_targets() -> None:
    core = object.__new__(MLMMCore)
    core.selection_indices = [0]
    core._explicit_hess_mm_atoms = None
    core._explicit_movable_mm_atoms = None
    core._explicit_frozen_mm_atoms = None
    core.use_bfactor_layers = False
    core.hess_cutoff = None
    core.movable_cutoff = 1.5

    core._compute_layer_indices(
        np.array(
            [
                [0.0, 0.0, 0.0],
                [1.0, 0.0, 0.0],
                [3.0, 0.0, 0.0],
            ]
        )
    )

    assert core.hess_mm_indices == [1]
    assert core.movable_mm_indices == []
    assert core.frozen_layer_indices == [2]
    assert core.movable_indices == [0, 1]


def test_calculator_attach_preserves_distance_layer_freezes() -> None:
    class FakeCore:
        frozen_layer_indices = [3]
        freeze_atoms = [3]

        def _update_active_dof_mappings(self) -> None:
            self.updated = True

    calculator = object.__new__(mlmm)
    calculator.core = FakeCore()

    calculator.freeze_atoms = [1]

    assert calculator.freeze_atoms == [1, 3]
    assert calculator.core.updated is True
