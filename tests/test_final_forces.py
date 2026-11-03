"""result.json final forces come from the final geometry."""

from __future__ import annotations

from types import SimpleNamespace

import numpy as np

from pysisyphus.Geometry import Geometry
from pysisyphus.calculators.Calculator import Calculator
from pysisyphus.optimizers.LBFGS import LBFGS
from mlmm.core.utils import optimizer_final_forces


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
