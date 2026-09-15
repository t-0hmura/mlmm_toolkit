"""The first exact minimum Hessian can reveal a smaller working space."""

import numpy as np
import pytest
import torch

from pysisyphus.Geometry import Geometry
from pysisyphus.calculators.Calculator import Calculator
from pysisyphus.optimizers.RFOptimizer import RFOptimizer
from pysisyphus.normal_modes import _frequencies_cm_and_modes, resolved_imaginary_mask


class PartialHarmonic(Calculator):
    """A stationary six-DOF harmonic PES reporting an atom-zero Hessian block."""

    def __init__(self, reference, out_dir):
        super().__init__(out_dir=out_dir, check_mem=False)
        self.reference = np.asarray(reference).copy()
        self.hessian_calls = 0

    def get_forces(self, atoms, coords, **kwargs):
        delta = np.asarray(coords) - self.reference
        # Neither a Hessian nor partial-space metadata is available until
        # the optimizer explicitly asks for the first physical Hessian.
        return {"energy": float(0.5 * delta @ delta), "forces": -delta}

    get_energy = get_forces

    def get_hessian(self, atoms, coords, **kwargs):
        self.hessian_calls += 1
        result = self.get_forces(atoms, coords)
        result.update(
            hessian=np.eye(3),
            within_partial_hessian={
                "active_atoms": np.array([0], dtype=int),
                "active_dofs": np.array([0, 1, 2], dtype=int),
                "active_n_dof": 3,
                "full_n_dof": 6,
            },
        )
        return result


def test_positive_exact_refresh_rebuilds_proposal_in_new_partial_space(tmp_path):
    coords = np.array([-0.5, 0.0, 0.0, 0.5, 0.0, 0.0])
    geom = Geometry(["H", "H"], coords.copy(), coord_type="cart", freeze_atoms=[])
    calculator = PartialHarmonic(coords, tmp_path)
    geom.set_calculator(calculator)
    opt = RFOptimizer(
        geom,
        hessian_init="unit",
        hessian_recalc=500,
        trust_radius=0.1,
        trust_min=1e-4,
        trust_max=0.1,
        thresh="baker",
        max_cycles=1,
        gdiis=False,
        line_search=False,
        out_dir=tmp_path,
        dump=False,
    )

    opt.run()

    assert calculator.hessian_calls == 1
    assert opt.using_active_dofs
    assert opt.cur_H.shape == (3, 3)
    assert opt.is_converged
    assert opt._minimum_curvature_valid
    assert opt._minimum_matches_current_geometry()
    assert len(opt.steps) == len(opt.predicted_energy_changes) == 1
    np.testing.assert_array_equal(opt.steps[0], np.zeros(6))
    np.testing.assert_array_equal(geom.cart_coords, coords)
    assert opt.predicted_energy_changes[0] == 0.0
    assert geom._hessian is None


@pytest.mark.parametrize("backend", ["numpy", "torch", "cuda"])
def test_minimum_phva_respects_ordered_partial_hessian_metadata(backend):
    if backend == "cuda" and not torch.cuda.is_available():
        pytest.skip("CUDA is unavailable")
    coords = np.array(
        [[1., 0., 0.], [0., 0., 0.], [0., 1., 0.], [0., 0., 1.]]
    )
    geom = Geometry(["H"] * 4, coords.ravel(), coord_type="cart", freeze_atoms=[1, 3])
    geom.tr_projection = "constrained"
    geom.within_partial_hessian = {
        "active_atoms": np.array([2, 0], dtype=int),
        "active_dofs": np.array([6, 7, 8, 0, 1, 2], dtype=int),
        "active_n_dof": 6,
        "full_n_dof": 12,
    }
    # Geometry explicitly supports the declared atom order, not only sorted
    # freeze-complement order. The surviving rigid motion is rotation about z.
    active_atoms, active_dofs, _, _ = geom._validated_partial_hessian_layout()
    np.testing.assert_array_equal(active_atoms, [2, 0])
    np.testing.assert_array_equal(active_dofs, [6, 7, 8, 0, 1, 2])
    true_rigid = np.array([-1., 0., 0., 0., 1., 0.]) / np.sqrt(2.)
    wrong_rigid = np.array([0., 1., 0., -1., 0., 0.]) / np.sqrt(2.)
    ordered_hessian = (
        np.eye(6) - np.outer(true_rigid, true_rigid)
        - 2. * np.outer(wrong_rigid, wrong_rigid)
    )

    opt = RFOptimizer.__new__(RFOptimizer)
    opt.geometry = geom
    opt.H = ordered_hessian.copy()
    if backend != "numpy":
        opt.H = torch.as_tensor(
            opt.H, dtype=torch.float64, device="cuda" if backend == "cuda" else "cpu"
        )
    opt.small_eigval_thresh = 1e-8
    opt.log = lambda *_: None
    opt._hessian_system(np.zeros(12))
    actual_frequencies, actual_modes = opt._mw_frequencies_and_modes()

    # Independent oracle: explicitly permute [atom2, atom0] to canonical
    # [atom0, atom2] before invoking the public-format frequency helper.
    canonical_dofs = [3, 4, 5, 0, 1, 2]
    canonical_hessian = ordered_hessian[np.ix_(canonical_dofs, canonical_dofs)]
    expected_frequencies, expected_modes = _frequencies_cm_and_modes(
        torch.tensor(canonical_hessian, dtype=torch.float64),
        [1, 1, 1, 1],
        coords.copy(),
        torch.device("cpu"),
        freeze_idx=[1, 3],
        tr_projection="constrained",
    )
    assert np.count_nonzero(resolved_imaginary_mask(expected_frequencies)) == 1
    assert np.count_nonzero(resolved_imaginary_mask(actual_frequencies)) == 1
    np.testing.assert_allclose(actual_frequencies, expected_frequencies, rtol=1e-10, atol=1e-6)
    physical_hessian = (
        opt.cur_H.detach().cpu().numpy() if isinstance(opt.cur_H, torch.Tensor)
        else opt.cur_H
    )
    np.testing.assert_array_equal(physical_hessian, ordered_hessian)
    actual_negative = actual_modes[0].detach().cpu().numpy()
    expected_negative = expected_modes[0].detach().cpu().numpy()
    assert abs(actual_negative @ expected_negative) == pytest.approx(1.)
    np.testing.assert_array_equal(geom.within_partial_hessian["active_atoms"], [2, 0])
