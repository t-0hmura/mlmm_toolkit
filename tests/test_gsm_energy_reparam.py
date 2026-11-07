"""Energy-weighted reparametrization keeps the CI and avoids spline excursions."""

from types import SimpleNamespace

import numpy as np
import pytest

from pysisyphus.cos.GrowingString import GrowingString


def string_fixture(points, param="energy", fully_grown=True, climbing=(2,)):
    density = np.linspace(0.0, 1.0, len(points))
    weighted = np.array([0.0, 0.001, 0.002, 0.8, 0.9, 1.0])
    obj = SimpleNamespace(
        coords=points.flatten().copy(), coords_length=3, param=param,
        fully_grown=fully_grown,
        images=[SimpleNamespace(coords=p.copy()) for p in points],
        get_climbing_indices=lambda: climbing, log=lambda *_: None,
        get_cur_param_density=lambda kind=None: weighted if kind == "energy" else density,
    )
    obj.spline = lambda: GrowingString.spline(obj)
    return obj


def test_energy_reparam_does_not_throw_a_compressed_straight_path_outside_its_endpoints():
    points = np.column_stack((np.arange(6.0), np.zeros((6, 2))))
    obj = string_fixture(points)
    desired = np.linspace(0.0, 1.0, 6)
    original_desired = desired.copy()
    GrowingString.reparam_cart(obj, desired)
    actual = obj.coords.reshape(-1, 3)
    assert np.isfinite(actual).all()
    assert actual[:, 0].min() >= 0.0
    assert actual[:, 0].max() <= 5.0
    np.testing.assert_array_equal(actual[[0, 2, 5]], points[[0, 2, 5]])
    np.testing.assert_array_equal(actual[:, 1:], 0.0)
    np.testing.assert_array_equal(desired, original_desired)


@pytest.mark.parametrize("param,fully_grown", [("equi", True), ("energy", False)])
def test_equi_and_growth_keep_the_existing_cubic_interpolation(param, fully_grown):
    x = np.linspace(0.0, 1.0, 6)
    points = np.column_stack((x, x**2, np.zeros(6)))
    obj = string_fixture(points, param, fully_grown)
    desired = np.array([0.0, 0.08, 0.22, 0.5, 0.9, 1.0])
    expected = np.column_stack((desired, desired**2, np.zeros(6)))
    expected[2] = points[2]
    GrowingString.reparam_cart(obj, desired)
    np.testing.assert_allclose(obj.coords.reshape(-1, 3), expected, atol=1e-13)
