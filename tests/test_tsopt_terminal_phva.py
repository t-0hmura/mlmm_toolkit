"""Terminal PHVA runs after convergence or a plateau, not max-cycle expiry."""

from pathlib import Path

import numpy as np
import torch

from mlmm.workflows import tsopt


class _Geometry:
    atomic_numbers = np.array([1])
    atoms = ["H"]
    cart_coords = np.zeros(3)
    coords3d = np.zeros((1, 3))
    within_partial_hessian = None


def _runner(tmp_path, monkeypatch, *, stalled):
    runner = object.__new__(tsopt.HessianDimer)
    runner.geom = _Geometry()
    runner.dump = False
    runner.optim_all_path = tmp_path / "optimization_all_trj.xyz"
    runner.mode_path = tmp_path / ".dimer_mode.dat"
    runner.out_dir = tmp_path
    runner.vib_dir = tmp_path / "vib"
    runner.freeze_atoms = []
    runner.masses_au_t = torch.ones(1)
    runner.device = torch.device("cpu")
    runner.root = 0
    runner.thresh_loose = "baker"
    runner.thresh = "baker"
    runner.flatten_max_iter = 1
    runner.flatten_loop_bofill = False
    runner.neg_freq_thresh_cm = 5.0
    runner.tr_projection = "none"
    runner.rigid_projection_info = {}
    runner.max_total_cycles = 1
    runner._cycles_spent = 0
    runner.is_converged = False
    runner.is_stalled = False
    runner.stop_reason = ""
    runner.n_imaginary_modes = None
    runner.imaginary_frequencies_cm = []
    runner.saddle_order_verified = False
    runner.skip_final_freq = stalled
    runner.ml_only_hessian_dimer = False
    runner.calc_kwargs_partial = {}
    runner.calc_kwargs_full = {}
    runner.calc_kwargs_ml_only = {}
    runner.analysis_active_atoms = [0]
    runner.source_path = None

    hessian_calls = []
    mode_exports = []

    def fake_hessian(*args, **kwargs):
        hessian_calls.append(True)
        return torch.eye(3)

    def fake_loop(_threshold, **kwargs):
        runner._cycles_spent = runner.max_total_cycles
        if stalled:
            runner.is_stalled = True
            runner.stop_reason = "energy plateau"
        return 1, False, False

    runner._compact_hessian_to_computed_coverage = lambda value: value
    runner._resolve_hessian_active_subspace = (
        lambda _hessian, _n_atoms: ([0], np.ones(3, dtype=bool))
    )
    runner._dimer_loop = fake_loop
    monkeypatch.setattr(tsopt, "_calc_full_hessian_torch", fake_hessian)
    monkeypatch.setattr(
        tsopt,
        "_mode_direction_by_root",
        lambda *args, **kwargs: np.ones((1, 3)),
    )
    monkeypatch.setattr(
        tsopt,
        "_reconcile_hessian_analysis_basis",
        lambda hessian, *args, **kwargs: (hessian, [0], [0], "full"),
    )
    monkeypatch.setattr(
        tsopt,
        "_modes_from_Hact_embedded",
        lambda *args, **kwargs: (
            np.array([-100.0, 20.0, 30.0]),
            torch.eye(3),
        ),
    )
    monkeypatch.setattr(
        tsopt,
        "_write_all_imag_modes",
        lambda *args, **kwargs: mode_exports.append(True) or 1,
    )
    monkeypatch.setattr(tsopt, "_clear_cuda_cache", lambda: None)
    monkeypatch.setattr(
        tsopt,
        "write",
        lambda path, _atoms: Path(path).write_text("final\n", encoding="utf-8"),
    )
    return runner, hessian_calls, mode_exports


def test_dimer_max_cycles_skips_terminal_phva(monkeypatch, tmp_path):
    runner, hessian_calls, mode_exports = _runner(
        tmp_path, monkeypatch, stalled=False
    )

    runner.run()

    assert len(hessian_calls) == 1
    assert mode_exports == []
    assert runner.n_imaginary_modes is None
    assert not runner.vib_dir.exists()


def test_dimer_plateau_runs_terminal_phva_even_if_skip_requested(
    monkeypatch, tmp_path,
):
    runner, hessian_calls, mode_exports = _runner(
        tmp_path, monkeypatch, stalled=True
    )

    runner.run()

    assert len(hessian_calls) == 2
    assert mode_exports == [True]
    assert runner.n_imaginary_modes == 1
    assert runner.is_stalled is True
