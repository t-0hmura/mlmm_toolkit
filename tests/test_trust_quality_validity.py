"""Scalar energy-quality validity without changing finite trust policy."""
import warnings

import numpy as np
import pytest

from pysisyphus.Geometry import Geometry
from pysisyphus.calculators.Calculator import Calculator
from pysisyphus.optimizers.RFOptimizer import RFOptimizer
from pysisyphus.optimizers.exceptions import OptimizationError


SCALARS = [
    pytest.param(float, id="python-float"),
    pytest.param(np.float64, id="numpy-scalar"),
    pytest.param(lambda value: np.array(value, dtype=np.float64), id="numpy-zero-dimensional"),
]


def make_optimizer(tmp_path, *, band=False, reject_uphill=False):
    geometry = Geometry(["H"], np.array([1.0, 0.0, 0.0]), coord_type="cart")
    optimizer = RFOptimizer(
        geometry, hessian_init="unit", max_cycles=1, dump=False,
        trust_radius=0.01, trust_min=1e-8, trust_max=1.0, trust_band=band,
        reject_uphill=reject_uphill, line_search=False, gdiis=False,
        out_dir=tmp_path,
    )
    optimizer.forces = [np.zeros(3), np.zeros(3)]
    optimizer.steps = [np.array([0.01, 0.0, 0.0])]
    optimizer.coords = [np.zeros(3), optimizer.steps[0].copy()]
    optimizer.cart_coords = [point.copy() for point in optimizer.coords]
    return optimizer


def remember_histories(optimizer):
    names = ("energies", "predicted_energy_changes", "forces", "steps", "coords", "cart_coords")
    return {
        name: (getattr(optimizer, name), list(getattr(optimizer, name)),
               [np.array(value, copy=True) for value in getattr(optimizer, name)])
        for name in names
    }


def assert_histories_unchanged(optimizer, before):
    for name, (original_list, original_items, copies) in before.items():
        current = getattr(optimizer, name)
        assert current is original_list
        assert len(current) == len(original_items)
        for value, original, saved in zip(current, original_items, copies):
            assert value is original
            np.testing.assert_array_equal(value, saved)


@pytest.mark.parametrize("scalar", SCALARS)
@pytest.mark.parametrize("band", [False, True])
@pytest.mark.parametrize("prediction", [0.0, -0.0], ids=["positive-zero", "negative-zero"])
@pytest.mark.parametrize("actual", [0.0, 1e-5, -1e-5])
def test_exact_zero_prediction_holds_without_division(
    tmp_path, monkeypatch, scalar, band, prediction, actual
):
    optimizer = make_optimizer(tmp_path, band=band)
    optimizer.energies = [scalar(0.0), scalar(actual)]
    optimizer.predicted_energy_changes = [scalar(prediction)]
    before = remember_histories(optimizer)
    messages = []
    monkeypatch.setattr(optimizer, "log", messages.append)

    def forbidden_set(*args, **kwargs):
        pytest.fail("Unavailable quality must not reach coefficient-based policy.")

    monkeypatch.setattr(optimizer, "set_new_trust_radius", forbidden_set)
    with warnings.catch_warnings():
        warnings.simplefilter("error", RuntimeWarning)
        with np.errstate(all="raise"):
            assert not optimizer.update_trust_radius()
    assert optimizer.trust_radius == 0.01
    assert any("zero predicted change" in message for message in messages)
    assert_histories_unchanged(optimizer, before)


@pytest.mark.parametrize("scalar", SCALARS)
@pytest.mark.parametrize("band", [False, True])
@pytest.mark.parametrize("operand", ["previous", "current", "predicted"])
@pytest.mark.parametrize("invalid", [np.nan, np.inf, -np.inf])
def test_nonfinite_energy_quality_fails_before_radius_mutation(
    tmp_path, monkeypatch, scalar, band, operand, invalid
):
    optimizer = make_optimizer(tmp_path, band=band)
    values = dict(previous=0.0, current=1e-5, predicted=1e-5)
    values[operand] = invalid
    optimizer.energies = [scalar(values["previous"]), scalar(values["current"])]
    optimizer.predicted_energy_changes = [scalar(values["predicted"])]
    before = remember_histories(optimizer)
    calls = []
    monkeypatch.setattr(optimizer, "set_new_trust_radius", lambda *a, **k: calls.append((a, k)))
    with warnings.catch_warnings():
        warnings.simplefilter("error", RuntimeWarning)
        with np.errstate(all="raise"), pytest.raises(OptimizationError, match="Non-finite"):
            optimizer.update_trust_radius()
    assert not calls
    assert optimizer.trust_radius == 0.01
    assert_histories_unchanged(optimizer, before)


@pytest.mark.parametrize("scalar", SCALARS)
@pytest.mark.parametrize("band", [False, True])
@pytest.mark.parametrize(
    "previous,current,prediction,message",
    [
        (-1e308, 1e308, 1.0, "Non-finite actual energy change"),
        (0.0, 1e308, 1e-308, "ratio overflow from finite changes"),
    ],
    ids=["finite-energies-change-overflow", "finite-changes-ratio-overflow"],
)
def test_derived_overflow_is_an_explicit_numerical_error(
    tmp_path, scalar, band, previous, current, prediction, message
):
    optimizer = make_optimizer(tmp_path, band=band)
    optimizer.energies = [scalar(previous), scalar(current)]
    optimizer.predicted_energy_changes = [scalar(prediction)]
    before = remember_histories(optimizer)
    with warnings.catch_warnings():
        warnings.simplefilter("error", RuntimeWarning)
        with np.errstate(all="raise"), pytest.raises(OptimizationError, match=message):
            optimizer.update_trust_radius()
    assert optimizer.trust_radius == 0.01
    assert_histories_unchanged(optimizer, before)


@pytest.mark.parametrize("scalar", SCALARS)
@pytest.mark.parametrize("band", [False, True])
@pytest.mark.parametrize("invalid", [np.nan, np.inf, -np.inf])
def test_direct_nonfinite_coefficient_is_not_a_policy_input(tmp_path, scalar, band, invalid):
    optimizer = make_optimizer(tmp_path, band=band)
    optimizer.energies = [0.0, 0.0]
    optimizer.predicted_energy_changes = [0.0]
    before = remember_histories(optimizer)
    with pytest.raises(ValueError, match="coefficient must be finite"):
        optimizer.set_new_trust_radius(scalar(invalid), 0.01)
    assert optimizer.trust_radius == 0.01
    assert_histories_unchanged(optimizer, before)


@pytest.mark.parametrize("scalar", SCALARS)
@pytest.mark.parametrize("band", [False, True])
@pytest.mark.parametrize(
    "actual,prediction,default_radius,band_radius",
    [
        (-8.03394e-22, -8.03394e-22, 0.02, 0.0115),
        (1e-300, 1e-300, 0.02, 0.0115),
        (-1e-300, -1e-300, 0.02, 0.0115),
        (-9e-5, -1e-5, 0.02, 0.0065),
        (9e-5, 1e-5, 0.02, 0.0065),
        (9.99978e-15 * 10001.22, 9.99978e-15, 0.02, 0.0065),
        (1e308, 1.0, 0.02, 0.0065),
        (0.0, 1e-5, 0.0025, 0.0065),
        (1e-5, -1e-5, 0.0025, 0.0065),
    ],
    ids=["signed-cancellation-is-nonzero", "tiny-positive", "tiny-negative",
         "reliable-large-lowering", "reliable-large-raising", "finite-sensitive-ratio",
         "largest-scale-finite-ratio", "zero-actual-finite-prediction", "opposite-signs"],
)
def test_finite_nonzero_quality_retains_existing_policy(
    tmp_path, scalar, band, actual, prediction, default_radius, band_radius
):
    optimizer = make_optimizer(tmp_path, band=band)
    optimizer.energies = [scalar(0.0), scalar(actual)]
    optimizer.predicted_energy_changes = [scalar(prediction)]
    before = remember_histories(optimizer)
    unexpected = optimizer.update_trust_radius()
    assert bool(unexpected) == (actual > 0 and prediction < 0)
    assert optimizer.trust_radius == pytest.approx(
        band_radius if band else default_radius, rel=1e-14, abs=0.0
    )
    assert_histories_unchanged(optimizer, before)


def test_existing_maximum_energy_increase_check_remains(tmp_path):
    optimizer = make_optimizer(tmp_path)
    optimizer.energies = [0.0, 0.1]
    optimizer.predicted_energy_changes = [-0.1]
    optimizer.max_energy_incr = 0.01
    before = remember_histories(optimizer)
    with pytest.raises(OptimizationError, match="Actual energy change too high"):
        optimizer.update_trust_radius()
    assert optimizer.trust_radius == 0.01
    assert_histories_unchanged(optimizer, before)


class QuadraticCalculator(Calculator):
    @staticmethod
    def results_at(coords):
        coords = np.asarray(coords, dtype=float)
        return {"energy": float(coords @ coords), "forces": -2.0 * coords}

    def get_energy(self, atoms, coords, **kwargs):
        return self.results_at(coords)

    def get_forces(self, atoms, coords, **kwargs):
        return self.results_at(coords)

    def get_hessian(self, atoms, coords, **kwargs):
        result = self.results_at(coords)
        result["hessian"] = 2.0 * np.eye(len(coords))
        return result


@pytest.mark.parametrize("scalar", SCALARS)
@pytest.mark.parametrize("band", [False, True])
def test_zero_prediction_does_not_bypass_optional_uphill_rollback(tmp_path, scalar, band):
    optimizer = make_optimizer(tmp_path, band=band, reject_uphill=True)
    geometry = optimizer.geometry
    geometry.set_calculator(QuadraticCalculator(out_dir=tmp_path, check_mem=False))
    accepted = geometry.coords.copy()
    energy = geometry.energy
    forces = geometry.forces.copy()
    step = np.array([1.0, 0.0, 0.0])
    geometry.coords = accepted + step
    optimizer.coords = [accepted.copy(), geometry.coords.copy()]
    optimizer.cart_coords = [point.copy() for point in optimizer.coords]
    optimizer.energies = [energy]
    optimizer.forces = [forces]
    optimizer.steps = [step]
    optimizer.predicted_energy_changes = [scalar(-0.0)]
    optimizer.image_inds = [[0], [0]]
    optimizer.image_nums = [1, 1]
    optimizer.cur_cycle = 1
    optimizer.H = np.eye(3)
    _, _, hessian, _, _, resetted = optimizer.housekeeping()
    np.testing.assert_array_equal(geometry.coords, accepted)
    np.testing.assert_array_equal(hessian, np.eye(3))
    assert resetted
    assert optimizer.trust_radius == 0.0025
    assert optimizer.rejected_uphill_steps == 1
    assert optimizer.energies == [energy]
    np.testing.assert_array_equal(optimizer.forces, [forces])
    assert len(optimizer.coords) == len(optimizer.cart_coords) == 1
    assert optimizer.steps == optimizer.predicted_energy_changes == []
