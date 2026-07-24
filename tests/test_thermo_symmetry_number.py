from __future__ import annotations

import math

import click
import numpy as np
import pytest

from mlmm.core.defaults import THERMO_KW
from mlmm.workflows.freq import (
    cli,
    _validated_symmetry_number,
    _validated_thermo_condition,
)
from thermoanalysis.QCData import QCData
from thermoanalysis.constants import KBAU
from thermoanalysis.thermo import thermochemistry


def _water_thermochemistry(symmetry_number: int):
    qc = QCData(
        {
            "coords3d": np.array(
                [[0.0, 0.0, 0.0], [0.757, 0.586, 0.0], [-0.757, 0.586, 0.0]]
            ),
            "wavenumbers": np.array(
                [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1595.0, 3657.0, 3756.0]
            ),
            "scf_energy": -76.4,
            "masses": np.array([15.999, 1.008, 1.008]),
            "mult": 1,
        },
        point_group="c1",
        mult=1,
    )
    qc.symmetry_number = symmetry_number
    return thermochemistry(qc, 298.15, pressure=101325.0)


def test_symmetry_number_cli_is_positive_and_defaults_to_one() -> None:
    param = next(p for p in cli.params if p.name == "symmetry_number")
    assert isinstance(param.type, click.IntRange)
    assert param.type.min == 1
    assert param.default == THERMO_KW["symmetry_number"] == 1


@pytest.mark.parametrize("value", [False, 0, -1, 1.5, "2", None])
def test_yaml_symmetry_number_rejects_non_positive_integers(value) -> None:
    with pytest.raises(click.UsageError, match="integer greater than or equal to 1"):
        _validated_symmetry_number(value)


@pytest.mark.parametrize(
    "value", [None, "bad", False, True, 0, -1, math.nan, math.inf, -math.inf]
)
def test_thermochemistry_conditions_reject_non_positive_or_nonfinite(value) -> None:
    with pytest.raises(click.UsageError, match="finite number greater than zero"):
        _validated_thermo_condition(value, name="temperature")


@pytest.mark.parametrize("value, expected", [(1, 1.0), ("2.5", 2.5), (298.15, 298.15)])
def test_thermochemistry_conditions_accept_positive_finite_values(
    value, expected
) -> None:
    assert _validated_thermo_condition(value, name="pressure_atm") == expected


def test_rotational_symmetry_number_changes_entropy_and_gibbs_exactly_once() -> None:
    sigma1 = _water_thermochemistry(1)
    sigma2 = _water_thermochemistry(2)
    expected = KBAU * 298.15 * math.log(2.0)
    np.testing.assert_allclose(
        float(sigma2.G - sigma1.G), expected, rtol=1.0e-10, atol=2.0e-14
    )
    np.testing.assert_allclose(
        float(sigma1.S_tot - sigma2.S_tot),
        KBAU * math.log(2.0),
        rtol=1.0e-12,
        atol=1.0e-15,
    )
