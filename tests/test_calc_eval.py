import math

import pytest

from mlmm.core.calc_eval import calc_energy
from mlmm.workflows.freq import _calc_full_hessian_torch


class _Geometry:
    atoms = ["H"]
    cart_coords = [0.0, 0.0, 0.0]


class _Calculator:
    def __init__(self, result) -> None:
        self.result = result

    def get_energy(self, atoms, coords):
        assert atoms == _Geometry.atoms
        assert coords == _Geometry.cart_coords
        return self.result


@pytest.mark.parametrize(
    "result",
    [
        {},
        {"energy": math.nan},
        {"energy": math.inf},
        {"energy": -math.inf},
    ],
)
def test_calc_energy_rejects_missing_or_nonfinite_energy(result) -> None:
    with pytest.raises((KeyError, ValueError)):
        calc_energy(_Geometry(), {}, calc=_Calculator(result))


def test_calc_energy_returns_finite_energy() -> None:
    assert calc_energy(
        _Geometry(),
        {},
        calc=_Calculator({"energy": -1.25}),
    ) == pytest.approx(-1.25)


@pytest.mark.parametrize(
    "result",
    [
        {"hessian": [[1.0]]},
        {"hessian": [[1.0]], "energy": math.nan},
        {"hessian": [[1.0]], "energy": math.inf},
    ],
)
def test_frequency_hessian_rejects_missing_or_nonfinite_energy(result) -> None:
    class Calculator:
        def get_hessian(self, atoms, coords):
            return result

    with pytest.raises((KeyError, ValueError)):
        _calc_full_hessian_torch(
            _Geometry(),
            {},
            device="cpu",
            calculator=Calculator(),
        )
