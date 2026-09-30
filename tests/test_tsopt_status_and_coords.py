"""tsopt/opt keep an explicit coord_type, and a failed final TS energy is a command error."""

from __future__ import annotations

import json
from types import SimpleNamespace

import numpy as np
import pytest
from click.testing import CliRunner

from pysisyphus.Geometry import Geometry
from mlmm.workflows import freq
from mlmm.workflows import opt as opt_workflow
from mlmm.workflows import tsopt

_WATER = ("O", "H", "H")
_WATER_BOHR = np.array([0.0, 0.0, 0.0, 1.8, 0.0, 0.0, -0.4, 1.7, 0.0])


class _LayerlessCalculator:
    def __init__(self, **_kwargs):
        self.core = SimpleNamespace(
            ml_indices=[], hess_mm_indices=[], movable_mm_indices=[],
            frozen_layer_indices=[],
        )


def _prepared_inputs(tmp_path):
    source, parm, xyz = (tmp_path / name for name in ("input.pdb", "input.parm7", "input.xyz"))
    source.write_text("Input parsing is replaced at the prepared-structure boundary.\n")
    parm.write_text("No Amber parser or MM engine is constructed.\n")
    ang = _WATER_BOHR.reshape(-1, 3) * 0.529177
    xyz.write_text("3\nwater\n" + "".join(
        f"{sym} {x:.6f} {y:.6f} {z:.6f}\n" for sym, (x, y, z) in zip(_WATER, ang)
    ))
    prepared = SimpleNamespace(original_path=source, source_path=source,
                               geom_path=xyz, cleanup=lambda: None)
    return source, parm, prepared


@pytest.fixture
def dimer_cli(monkeypatch, tmp_path):
    """Run tsopt --opt-mode grad up to a fake HessianDimer that records its kwargs."""
    captured = {}

    class FakeDimer:
        def __init__(self, fn, **kwargs):
            captured.update(kwargs)
            self.geom = Geometry(_WATER, _WATER_BOHR.copy())
            self.calc_kwargs = dict(kwargs["calc_kwargs"])
            self.is_converged, self.is_stalled, self.stop_reason = True, False, ""
            self._cycles_spent, self.flatten_skip_reason = 1, None
            self.n_imaginary_modes = self.n_negative_modes = 1
            self.imaginary_frequencies_cm = [-500.0]
            self.hessian_status, self.hessian_error = "completed", None
            self.rigid_projection_info = {}

        def run(self):
            pass

    source, parm, prepared = _prepared_inputs(tmp_path)
    monkeypatch.chdir(tmp_path)
    monkeypatch.setattr(tsopt, "prepare_input_structure", lambda *_a, **_k: prepared)
    monkeypatch.setattr(tsopt, "resolve_charge_spin_or_raise", lambda *_a, **_k: (0, 1))
    monkeypatch.setattr(tsopt, "resolve_ml_layer_assignment", lambda **_k: (source, None))
    for module in (tsopt, freq):
        monkeypatch.setattr(module, "mlmm", _LayerlessCalculator)
    monkeypatch.setattr(tsopt, "HessianDimer", FakeDimer)

    def invoke(*extra):
        out_dir = tmp_path / "ts"
        result = CliRunner().invoke(tsopt.cli, [
            "-i", str(source), "--parm", str(parm), "-q", "0", "-m", "1",
            "-o", str(out_dir), "--opt-mode", "grad", "--active-dof-mode", "all",
            "--no-flatten", "--no-dump", "--no-convert-files", "--out-json", *extra,
        ])
        report_path = out_dir / "result.json"
        report = json.loads(report_path.read_text()) if report_path.exists() else None
        return result, report, captured

    return invoke


@pytest.mark.parametrize("coord_type", ["redund", "dlc"])
def test_tsopt_dimer_keeps_explicit_internal_coord_type(dimer_cli, monkeypatch, coord_type):
    monkeypatch.setattr(tsopt, "_calc_energy", lambda *_a, **_k: -1.0)
    result, report, captured = dimer_cli("--coord-type", coord_type)
    assert result.exit_code == 0, result.output
    assert captured["geom_kwargs"]["coord_type"] == coord_type
    assert report["optimization_status"] == "converged"


def test_tsopt_dimer_default_coord_type_is_cart(dimer_cli, monkeypatch):
    monkeypatch.setattr(tsopt, "_calc_energy", lambda *_a, **_k: -1.0)
    result, report, captured = dimer_cli()
    assert result.exit_code == 0, result.output
    assert captured["geom_kwargs"]["coord_type"] == "cart"
    assert report["energy_hartree"] == -1.0


def test_tsopt_dimer_final_energy_failure_is_a_command_error(dimer_cli, monkeypatch):
    def failed_energy(*_args, **_kwargs):
        raise RuntimeError("injected final energy failure")

    monkeypatch.setattr(tsopt, "_calc_energy", failed_energy)
    result, report, _ = dimer_cli()
    # A converged optimizer without a final energy must not be reported as a usable TS.
    assert result.exit_code == 1, result.output
    assert report["execution_status"] == "failed"
    assert report["error_type"] == "RuntimeError"


class _HarmonicWater(_LayerlessCalculator):
    """Harmonic O-H/H-H springs around a slightly displaced start geometry."""

    _BONDS = ((0, 1, 1.9), (0, 2, 1.9), (1, 2, 2.9))

    def _energy_gradient(self, coords):
        xyz = np.asarray(coords, dtype=float).reshape(-1, 3)
        energy, grad = 0.0, np.zeros_like(xyz)
        for i, j, r0 in self._BONDS:
            vec = xyz[i] - xyz[j]
            r = np.linalg.norm(vec)
            energy += 0.5 * (r - r0) ** 2
            grad[i] += (r - r0) * vec / r
            grad[j] -= (r - r0) * vec / r
        return energy, grad.ravel()

    def get_energy(self, atoms, coords):
        return {"energy": self._energy_gradient(coords)[0]}

    def get_forces(self, atoms, coords):
        energy, grad = self._energy_gradient(coords)
        return {"energy": energy, "forces": -grad}


def test_opt_lbfgs_runs_in_explicit_dlc(tmp_path, monkeypatch):
    loaded = []

    def load_geometry(_path, **kwargs):
        loaded.append(Geometry(_WATER, _WATER_BOHR.copy(), **kwargs))
        return loaded[-1]

    def fake_layers(*, source_path, calc_cfg, **_kwargs):
        calc_cfg["model_pdb"] = str(source_path)
        return source_path, None

    source, parm, prepared = _prepared_inputs(tmp_path)
    monkeypatch.setattr(opt_workflow, "prepare_input_structure", lambda *_a: prepared)
    monkeypatch.setattr(opt_workflow, "load_pdb_atom_metadata", lambda *_a: [])
    monkeypatch.setattr(opt_workflow, "resolve_charge_spin_or_raise", lambda *_a, **_k: (0, 1))
    monkeypatch.setattr(opt_workflow, "resolve_ml_layer_assignment", fake_layers)
    monkeypatch.setattr(opt_workflow, "geom_loader", load_geometry)
    monkeypatch.setattr(opt_workflow, "mlmm", _HarmonicWater)
    out_dir = tmp_path / "opt"
    result = CliRunner().invoke(opt_workflow.cli, [
        "-i", str(source), "--parm", str(parm), "-q", "0", "-m", "1",
        "--opt-mode", "grad", "--coord-type", "dlc", "--max-cycles", "50",
        "--no-convert-files", "--out-json", "--out-dir", str(out_dir),
    ])
    assert result.exit_code == 0, result.output
    assert [geom.coord_type for geom in loaded] == ["dlc"]
    assert json.loads((out_dir / "result.json").read_text())["optimization_status"] == "converged"
