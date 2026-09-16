"""Terminal Hessian metadata must not reinterpret an existing TS proposal."""
import numpy as np
import pytest
import torch

from pysisyphus._array import as_numpy

from pysisyphus.Geometry import Geometry
from pysisyphus.calculators.Calculator import Calculator
from pysisyphus.tsoptimizers.RSPRFOptimizer import RSPRFOptimizer
from pysisyphus.tsoptimizers.RSIRFOptimizer import RSIRFOptimizer
from pysisyphus.tsoptimizers.TRIM import TRIM


OPTIMIZERS = [RSPRFOptimizer, RSIRFOptimizer, TRIM]


class FirstPartialHessian(Calculator):
    """First physical H discovers atom zero; the saved initial model is full."""

    def __init__(self, center, out_dir):
        super().__init__(out_dir=out_dir, check_mem=False)
        self.center = np.asarray(center).copy()
        self.diagonal = np.array([-4., -2., 10., 1., 1., 1.])
        self.hessian_points = []
        self.force_points = []

    def get_forces(self, atoms, coords, **kwargs):
        self.force_points.append(np.asarray(coords).copy())
        delta = np.asarray(coords) - self.center
        gradient = self.diagonal * delta
        return {"energy": float(.5 * delta @ gradient), "forces": -gradient}

    get_energy = get_forces

    def get_hessian(self, atoms, coords, **kwargs):
        self.hessian_points.append(np.asarray(coords).copy())
        result = self.get_forces(atoms, coords)
        result.update(
            hessian=np.diag(self.diagonal[:3]),
            within_partial_hessian={
                "active_atoms": np.array([0]), "active_dofs": np.array([0, 1, 2]),
                "active_n_dof": 3, "full_n_dof": 6,
            },
        )
        return result


@pytest.mark.parametrize("cls", OPTIMIZERS)
def test_terminal_full_model_to_first_partial_hessian_keeps_stationary_candidate(tmp_path, cls):
    center = np.array([-.5, 0., 0., .5, 0., 0.])
    geometry = Geometry(["H", "H"], center.copy(), coord_type="cart", freeze_atoms=[])
    calculator = FirstPartialHessian(center, tmp_path)
    geometry.set_calculator(calculator)
    model = np.eye(6)
    model[0, 0] = -1.
    model_file = tmp_path / "saved_full_model.dat"
    np.savetxt(model_file, model)
    optimizer = cls(
        geometry, hessian_init=str(model_file), hessian_recalc=500,
        reference_mode=np.array([1., 0., 0., 0., 0., 0.]),
        max_cycles=3, thresh="baker", trust_radius=.1, trust_max=.1,
        saddle_recovery_max_cycles=0, flatten_enabled=False,
        dump=False, out_dir=tmp_path,
    )
    assert geometry.within_partial_hessian is None
    optimizer.run()

    assert optimizer.is_converged
    assert optimizer.saddle_recovery_steps == 0
    assert not optimizer._saddle_recovery_active
    assert len(calculator.hessian_points) == 1
    assert optimizer.exact_saddle_checks == 1
    assert optimizer.cur_H.shape == (3, 3)
    np.testing.assert_array_equal(optimizer.active_dof_indices, [0, 1, 2])
    np.testing.assert_array_equal(optimizer.cur_H, np.diag([-4., -2., 10.]))
    assert optimizer._last_exact_n_imaginary == 2
    assert len(optimizer.steps) == len(optimizer.predicted_energy_changes) == 1
    np.testing.assert_array_equal(optimizer.steps[0], np.zeros(6))
    assert optimizer.predicted_energy_changes[0] == 0.
    np.testing.assert_array_equal(geometry.cart_coords, center)
    for point in calculator.force_points + calculator.hessian_points:
        np.testing.assert_array_equal(point, center)


@pytest.mark.parametrize("cls", OPTIMIZERS)
@pytest.mark.parametrize("initial_backend,terminal_backend", [
    ("numpy", "numpy"), ("numpy", "torch"), ("numpy", "cuda"),
    ("torch", "numpy"), ("cuda", "numpy"),
])
@pytest.mark.parametrize("new_order", [np.arange(6), np.array([3, 4, 5, 0, 1, 2])], ids=["same", "permuted"])
def test_terminal_same_size_map_retains_physical_proposal_and_prediction(
    tmp_path, monkeypatch, cls, new_order, initial_backend, terminal_backend,
):
    """Nonzero, distinct components discriminate an otherwise silent reorder.

    The real proposal solvers run. Only the terminal refresh is replaced by its
    H/map effect, so this checks bookkeeping without another physical H call.
    The actual physical refresh/no-motion route is covered above.
    """
    if "cuda" in (initial_backend, terminal_backend) and not torch.cuda.is_available():
        pytest.skip("CUDA is unavailable")

    def represented(array, backend):
        if backend == "numpy":
            return array.copy()
        return torch.as_tensor(array, dtype=torch.float64,
                               device="cuda" if backend == "cuda" else "cpu")

    geometry = Geometry(["H", "H"], [-.5, 0., 0., .5, 0., 0.], coord_type="cart")
    optimizer = cls(
        geometry, hessian_init="unit", roots=[0], verify_saddle=False,
        min_line_search=False, max_line_search=False, trust_radius=.1,
        trust_max=.1, max_micro_cycles=50, saddle_recovery_max_cycles=0,
        dump=False, out_dir=tmp_path,
    )
    optimizer._using_active_dofs = True
    optimizer._active_dof_indices = np.arange(6)
    values = np.array([-1., 2., 3., 4., 5., 6.])
    gradient = np.array([1., 3., 5., 7., 11., 13.]) * 1e-6
    model = np.diag(values)
    exact_full = np.diag([-4., -2., 10., 8., 6., 4.])
    optimizer.H = optimizer.cur_H = represented(model, initial_backend)
    optimizer.forces = [-gradient.copy()]
    monkeypatch.setattr(optimizer, "housekeeping", lambda: (
        0., represented(gradient, initial_backend), represented(model, initial_backend),
        represented(values, initial_backend), represented(np.eye(6), initial_backend), False,
    ))
    proposals = []

    def refresh(step):
        proposals.append(as_numpy(optimizer.full_from_active(step)).copy())
        # Mutate the array in place to ensure the pre-validation map was copied.
        optimizer._active_dof_indices[:] = new_order
        optimizer.H = optimizer.cur_H = represented(
            exact_full[np.ix_(new_order, new_order)], terminal_backend
        )

    monkeypatch.setattr(optimizer, "validate_terminal_saddle_for_step", refresh)
    returned = optimizer.optimize()
    assert len(proposals) == 1
    assert np.count_nonzero(proposals[0]) == 6
    np.testing.assert_allclose(returned, proposals[0], rtol=0., atol=1e-15)
    numerator = gradient @ returned + .5 * returned @ exact_full @ returned
    expected = numerator if cls is TRIM else numerator / (1. + returned @ returned)
    assert len(optimizer.predicted_energy_changes) == 1
    assert optimizer.predicted_energy_changes[0] == pytest.approx(expected, rel=1e-12, abs=1e-25)
    assert optimizer.saddle_recovery_steps == 0
    assert geometry.calculator is None


@pytest.mark.parametrize("device", ["cpu", "cuda"])
@pytest.mark.parametrize("changed_basis", [False, True], ids=["same-map", "new-partial-map"])
def test_rsirf_saved_numpy_model_to_physical_tensor_hessian(
    tmp_path, monkeypatch, device, changed_basis,
):
    """Exercise actual terminal H acquisition with one H and no motion."""
    if device == "cuda" and not torch.cuda.is_available():
        pytest.skip("CUDA is unavailable")

    class TensorPartialHessian(FirstPartialHessian):
        def get_hessian(self, atoms, coords, **kwargs):
            result = super().get_hessian(atoms, coords, **kwargs)
            result["hessian"] = torch.as_tensor(
                result["hessian"], dtype=torch.float64, device=device
            )
            return result

    center = np.array([-.5, 0., 0., .5, 0., 0.])
    geometry = Geometry(["H", "H"], center.copy(), coord_type="cart",
                        freeze_atoms=[] if changed_basis else [1])
    calculator = TensorPartialHessian(center, tmp_path)
    geometry.set_calculator(calculator)
    model = np.eye(6)
    model[0, 0] = -1.
    model_file = tmp_path / "saved_numpy_model.dat"
    np.savetxt(model_file, model)
    optimizer = RSIRFOptimizer(
        geometry, hessian_init=str(model_file), hessian_recalc=500,
        reference_mode=np.array([1., 0., 0., 0., 0., 0.]),
        max_cycles=3, thresh="baker", trust_radius=.1, trust_max=.1,
        saddle_recovery_max_cycles=0, flatten_enabled=False,
        dump=False, out_dir=tmp_path,
    )
    previous = optimizer.validate_terminal_saddle_for_step
    before_maps = []

    def observe(step):
        assert isinstance(optimizer.cur_H, np.ndarray)
        indices = optimizer.active_dof_indices
        before_maps.append(None if indices is None else indices.copy())
        previous(step)

    monkeypatch.setattr(optimizer, "validate_terminal_saddle_for_step", observe)
    optimizer.run()
    assert optimizer.is_converged
    assert len(before_maps) == 1
    if changed_basis:
        assert before_maps[0] is None
    else:
        np.testing.assert_array_equal(before_maps[0], [0, 1, 2])
    assert isinstance(optimizer.cur_H, torch.Tensor)
    assert optimizer.cur_H.device.type == device
    torch.testing.assert_close(optimizer.cur_H,
        torch.diag(torch.tensor([-4., -2., 10.], dtype=torch.float64, device=device)))
    assert len(calculator.hessian_points) == optimizer.exact_saddle_checks == 1
    assert optimizer.saddle_recovery_steps == 0
    assert not optimizer._saddle_recovery_active
    assert len(optimizer.steps) == len(optimizer.predicted_energy_changes) == 1
    np.testing.assert_array_equal(optimizer.steps[0], np.zeros(6))
    assert optimizer.predicted_energy_changes[0] == 0.
    np.testing.assert_array_equal(geometry.cart_coords, center)
    for point in calculator.force_points + calculator.hessian_points:
        np.testing.assert_array_equal(point, center)
