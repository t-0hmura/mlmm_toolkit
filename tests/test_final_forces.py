"""result.json final forces come from the final geometry, the final force and
step fields agree with run.log, and convergence_thresholds are the thresholds
of the optimizer that ran last."""

from __future__ import annotations

import json
import re
from types import SimpleNamespace

import numpy as np
import pytest
import torch
from click.testing import CliRunner

from pysisyphus.Geometry import Geometry
from pysisyphus.calculators.Calculator import Calculator
from pysisyphus.optimizers.LBFGS import LBFGS
from pysisyphus.optimizers.Optimizer import CONV_THRESHS
from mlmm.cli import cli as root_cli
from mlmm.core import calc_eval
from mlmm.core.utils import optimizer_final_forces
from mlmm.io import hessian_cache
from mlmm.workflows import freq, tsopt
from mlmm.workflows import opt as opt_workflow


class _QuadraticCalculator(Calculator):
    def __init__(self, out_dir):
        super().__init__(out_dir=out_dir, check_mem=False)

    @staticmethod
    def _results(coords):
        coords = np.asarray(coords, dtype=float)
        return {"energy": float(coords @ coords), "forces": -2.0 * coords}

    def get_energy(self, atoms, coords, **prepare_kwargs):
        return self._results(coords)

    def get_forces(self, atoms, coords, **prepare_kwargs):
        return self._results(coords)


def _run_lbfgs(tmp_path, atoms, coords, freeze_atoms=None, **kwargs):
    geom = Geometry(
        atoms, np.asarray(coords, dtype=float), coord_type="cart",
        freeze_atoms=freeze_atoms or [],
    )
    geom.set_calculator(_QuadraticCalculator(tmp_path))
    opt = LBFGS(geom, line_search=False, out_dir=tmp_path, **kwargs)
    opt.run()
    return geom, opt


def test_max_cycles_stop_reports_final_geometry_forces(tmp_path, capsys) -> None:
    geom, opt = _run_lbfgs(
        tmp_path, ["H", "H"], [1.0, 0.0, 0.0, 0.0, 0.5, 0.0],
        max_cycles=1, max_step=0.1,
    )
    assert not opt.is_converged
    # A step followed the only evaluation, so the history is one geometry behind.
    assert not np.array_equal(geom.cart_coords, opt.cart_coords[-1])

    final_forces = -2.0 * geom.cart_coords
    expected = (
        float(np.abs(final_forces).max()),
        float(np.sqrt(np.mean(final_forces**2))),
    )
    assert optimizer_final_forces(opt) == expected
    assert expected[0] != opt.max_forces[-1]
    # The JSON values are the ones the final summary prints.
    out = capsys.readouterr().out
    assert f"max(forces,cartesian): {expected[0]:.6f} hartree/bohr" in out
    assert f"rms(forces,cartesian): {expected[1]:.6f} hartree/bohr" in out


def test_max_cycles_stop_keeps_active_dof_basis_with_frozen_atoms(tmp_path) -> None:
    geom, opt = _run_lbfgs(
        tmp_path, ["H", "H"], [1.0, 0.0, 0.0, 0.0, 0.5, 0.0],
        freeze_atoms=[1], max_cycles=1, max_step=0.1,
    )
    active_forces = (-2.0 * geom.cart_coords)[:3]
    assert optimizer_final_forces(opt) == (
        float(np.abs(active_forces).max()),
        float(np.sqrt(np.mean(active_forces**2))),
    )


def test_converged_run_keeps_history_values(tmp_path) -> None:
    geom, opt = _run_lbfgs(
        tmp_path, ["H"], [1.0e-4, 0.0, 0.0], thresh="gau_loose", max_cycles=5,
    )
    assert opt.is_converged
    assert opt.max_forces[-1] > 0.0
    assert optimizer_final_forces(opt) == (opt.max_forces[-1], opt.rms_forces[-1])


def test_final_forces_fall_back_to_history_without_final_summary() -> None:
    stopped = SimpleNamespace(
        final_max_force=None, final_rms_force=None,
        max_forces=[0.3, 0.2], rms_forces=[0.1, 0.05],
    )
    assert optimizer_final_forces(stopped) == (0.2, 0.05)
    assert optimizer_final_forces(SimpleNamespace(max_forces=[], rms_forces=[])) is None
    assert optimizer_final_forces(None) is None


# E = 1/2 sum_i k_i (x_i - x0_i)^2 (Hartree, Bohr); one negative k makes a saddle.
_X0 = np.array([
    [0.0, 0.0, 0.0], [2.0, 0.0, 0.0], [0.0, 2.0, 0.0],
    [0.0, 0.0, 2.0], [2.0, 2.0, 2.0], [3.0, 1.0, 2.0],
]).ravel()
_START = _X0.copy()
_START[[9, 10, 15, 16]] += [0.4, 0.5, -0.6, 0.7]
_SADDLE = np.arange(1.0, 19.0) * 0.05
_SADDLE[9] *= -1.0


class _QuadraticSaddle:
    """ML/MM stand-in: atoms 0-2 frozen MM, 3 and 5 ML, 4 movable MM."""

    curvatures = _SADDLE

    def __init__(self, **kwargs):
        self.freeze_atoms = list(kwargs.get("freeze_atoms", []))
        self.core = SimpleNamespace(
            ml_indices=[3, 5], hess_mm_indices=[4], movable_mm_indices=[4],
            frozen_layer_indices=[0, 1, 2], hess_active_atoms=list(range(6)),
            mlmm_links=[], _ml_backend=None,
        )

    def get_energy(self, atoms, coords):
        d = np.asarray(coords, dtype=float).reshape(-1) - _X0
        return {"energy": float(0.5 * np.dot(self.curvatures * d, d))}

    def get_forces(self, atoms, coords):
        d = np.asarray(coords, dtype=float).reshape(-1) - _X0
        return {**self.get_energy(atoms, coords), "forces": -self.curvatures * d}

    def get_hessian(self, atoms, coords):
        return {
            **self.get_forces(atoms, coords),
            "hessian": torch.diag(torch.as_tensor(self.curvatures)),
        }


class _QuadraticMinimum(_QuadraticSaddle):
    curvatures = np.abs(_SADDLE)


def _run_cli(tmp_path, monkeypatch, module, calculator, command, *args):
    monkeypatch.chdir(tmp_path)
    monkeypatch.setattr(hessian_cache, "_cache", {})
    for patched in (module, freq, calc_eval):
        monkeypatch.setattr(patched, "mlmm", calculator)
    monkeypatch.setattr(module, "mlmm_mm_only", lambda _core, **kw: calculator(**kw))
    source, parm, xyz, config = (
        tmp_path / name for name in ("input.pdb", "input.parm7", "input.xyz", "cpu.yaml")
    )
    source.write_text("Input parsing is replaced at the prepared-structure boundary.\n")
    parm.write_text("No Amber parser or MM engine is constructed.\n")
    xyz.write_text("6\nquadratic\n" + "".join(
        f"H {x:.8f} {y:.8f} {z:.8f}\n"
        for x, y, z in _START.reshape(-1, 3) * 0.529177210903
    ))
    config.write_text(
        "calc:\n  ml_device: cpu\n"
        + ("hessian_dimer:\n  device: cpu\n" if command == "tsopt" else "")
    )
    prepared = SimpleNamespace(original_path=source, source_path=source,
                               geom_path=xyz, cleanup=lambda: None)

    def layers(*, source_path, calc_cfg, **_kwargs):
        calc_cfg["model_pdb"] = str(source_path)
        return source_path, None

    monkeypatch.setattr(module, "prepare_input_structure", lambda *_a, **_k: prepared)
    monkeypatch.setattr(module, "resolve_charge_spin_or_raise", lambda *_a, **_k: (0, 1))
    monkeypatch.setattr(module, "resolve_ml_layer_assignment", layers)
    monkeypatch.setattr(
        module, "geom_loader", lambda _path, **kw: Geometry(["H"] * 6, _START.copy(), **kw)
    )
    if module is opt_workflow:
        monkeypatch.setattr(module, "load_pdb_atom_metadata", lambda *_a: [])
    out_dir = tmp_path / "out"
    result = CliRunner().invoke(root_cli, [
        command, "-i", str(source), "--parm", str(parm), "-q", "0", "-m", "1",
        "--config", str(config), "--out-dir", str(out_dir), "--max-cycles", "3",
        "--print-every", "1", "--no-flatten", "--no-dump", "--no-convert-files",
        "--out-json", *args,
    ])
    report = json.loads((out_dir / "result.json").read_text(encoding="utf-8"))
    assert result.exit_code == 1, result.output
    assert report["execution_status"] == "completed"
    assert report["optimization_status"] == "not_converged"
    return report, (out_dir / "run.log").read_text(encoding="utf-8")


def _final_fields(report):
    return tuple(
        f"{report[key]:.6f}"
        for key in ("final_max_force", "final_rms_force", "final_max_step", "final_rms_step")
    )


def _thresholds(preset):
    """The numeric thresholds of a named preset, as result.json writes them."""
    keys = ("max_force_thresh", "rms_force_thresh", "max_step_thresh", "rms_step_thresh")
    return {key: float(value) for key, value in zip(keys, CONV_THRESHS[preset])}


def _last_final_summary(log):
    """(max, rms) forces printed in the last Final summary."""
    return (
        re.findall(r"max\(forces,cartesian\): (\S+)", log)[-1],
        re.findall(r"rms\(forces,cartesian\): (\S+)", log)[-1],
    )


def _last_cycle_row(log):
    """(max force, rms force, max step, rms step) of the last cycle row.

    The microiteration table has one more column (micro steps).
    """
    rows = [
        line.split() for line in log.splitlines()
        if len(line.split()) in (7, 8) and line.split()[0].isdigit()
    ]
    return tuple(value.rstrip("*") for value in rows[-1][2:6])


def test_tsopt_dimer_final_fields_match_run_log(tmp_path, monkeypatch):
    mode = np.zeros((6, 3))
    mode[3, 0] = 1.0
    monkeypatch.setattr(
        tsopt, "_calc_full_hessian_torch",
        lambda *_args, **_kwargs: torch.diag(torch.as_tensor(_SADDLE)),
    )
    monkeypatch.setattr(tsopt, "_mode_direction_by_root", lambda *_a, **_k: mode.copy())
    report, log = _run_cli(
        tmp_path, monkeypatch, tsopt, _QuadraticSaddle, "tsopt",
        "--opt-mode", "dimer", "--active-dof-mode", "partial",
    )
    assert "Number of cycles exceeded!" in log
    # Forces of the final geometry, as the last LBFGS segment's Final summary
    # prints them; the step that produced it, as the last cycle row prints it.
    assert _final_fields(report)[:2] == _last_final_summary(log)
    assert _final_fields(report)[2:] == _last_cycle_row(log)[2:]
    # The cycle budget ends in the loose loop, so its LBFGS thresholds were used.
    assert "Reached --max-cycles budget after loose loop" in log
    assert report["convergence_thresholds"] == _thresholds("gau_loose")


@pytest.mark.parametrize("opt_mode", ["hess", "rsirfo", "trim"])
def test_tsopt_hessian_final_fields_match_run_log(tmp_path, monkeypatch, opt_mode):
    report, log = _run_cli(
        tmp_path, monkeypatch, tsopt, _QuadraticSaddle, "tsopt",
        "--opt-mode", opt_mode, "--no-microiter", "--active-dof-mode", "all",
        "--freeze-atoms", "1,2,3",
    )
    assert "Number of cycles exceeded!" in log
    assert _final_fields(report)[:2] == _last_final_summary(log)
    assert _final_fields(report)[2:] == _last_cycle_row(log)[2:]
    assert report["convergence_thresholds"] == _thresholds(report["thresh"])


@pytest.mark.parametrize("opt_mode", ["hess", "rsirfo", "trim"])
def test_tsopt_microiteration_final_fields_match_last_macro_row(
    tmp_path, monkeypatch, opt_mode,
):
    report, log = _run_cli(
        tmp_path, monkeypatch, tsopt, _QuadraticSaddle, "tsopt",
        "--opt-mode", opt_mode, "--microiter", "--active-dof-mode", "all",
        "--freeze-atoms", "1,2,3",
    )
    assert report["microiteration"]["used"] is True
    assert "[microiter] Reached max macro iterations (3)." in log
    # No Final summary on this path: the values are the last macro (ONIOM) cycle's.
    assert _final_fields(report) == _last_cycle_row(log)
    assert report["convergence_thresholds"] == _thresholds(report["thresh"])


@pytest.mark.parametrize("microiter", ["--microiter", "--no-microiter"])
def test_opt_hess_final_fields_match_run_log(tmp_path, monkeypatch, microiter):
    report, log = _run_cli(
        tmp_path, monkeypatch, opt_workflow, _QuadraticMinimum, "opt",
        "--opt-mode", "hess", microiter,
    )
    if microiter == "--microiter":
        assert report["microiteration"]["used"] is True
        assert _final_fields(report) == _last_cycle_row(log)
    else:
        assert "microiteration" not in report
        assert _final_fields(report)[:2] == _last_final_summary(log)
        assert _final_fields(report)[2:] == _last_cycle_row(log)[2:]
    # The macro (outer ONIOM) optimizer's thresholds on the microiteration path.
    assert report["convergence_thresholds"] == _thresholds(report["thresh"])
