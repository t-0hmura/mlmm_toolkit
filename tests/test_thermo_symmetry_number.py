from __future__ import annotations

import math

import click
import numpy as np
import pytest

from mlmm.core.defaults import THERMO_KW
from mlmm.workflows.freq import (
    cli,
    _symmetry_number_source,
    _validated_symmetry_number,
    _validated_thermo_condition,
)
from thermoanalysis.QCData import (
    QCData,
    detect_point_group_and_symmetry_number,
    symmetry_number_from_point_group,
)
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


def test_symmetry_number_is_automatic_without_a_public_cli_option() -> None:
    assert all(p.name != "symmetry_number" for p in cli.params)
    assert THERMO_KW["symmetry_number"] is None


@pytest.mark.parametrize("value", [False, 0, -1, 1.5, "2"])
def test_yaml_symmetry_number_rejects_non_positive_integers(value) -> None:
    with pytest.raises(click.UsageError, match="integer greater than or equal to 1"):
        _validated_symmetry_number(value)


def test_yaml_none_selects_automatic_symmetry_detection() -> None:
    assert _validated_symmetry_number(None) is None


@pytest.mark.parametrize(
    ("config_has_value", "override_has_value", "resolved_value", "expected"),
    [
        (False, False, None, "auto"),
        (True, False, 2, "config"),
        (True, True, 3, "override"),
        (True, True, None, "auto"),
    ],
)
def test_symmetry_number_source_tracks_configuration_precedence(
    config_has_value, override_has_value, resolved_value, expected
) -> None:
    assert _symmetry_number_source(
        config_has_value=config_has_value,
        override_has_value=override_has_value,
        resolved_value=resolved_value,
    ) == expected


@pytest.mark.parametrize(
    ("atomic_numbers", "coords3d", "expected"),
    [
        (
            [8, 1, 1],
            [[0.0, 0.0, 0.0], [0.757, 0.586, 0.0], [-0.757, 0.586, 0.0]],
            ("C2v", 2, "auto"),
        ),
        (
            [6, 1, 1, 1, 1],
            [
                [0.0, 0.0, 0.0],
                [1.0, 1.0, 1.0],
                [-1.0, -1.0, 1.0],
                [-1.0, 1.0, -1.0],
                [1.0, -1.0, -1.0],
            ],
            ("Td", 12, "auto"),
        ),
        ([1, 1], [[0.0, 0.0, -0.37], [0.0, 0.0, 0.37]], ("Dinfh", 2, "auto")),
        ([1], [[0.0, 0.0, 0.0]], ("Kh", 1, "auto")),
    ],
)
def test_automatic_symmetry_detection(atomic_numbers, coords3d, expected) -> None:
    point_group, symmetry_number, source = detect_point_group_and_symmetry_number(
        atomic_numbers,
        coords3d,
    )
    assert (point_group, symmetry_number, source) == expected


def test_automatic_symmetry_detection_has_conservative_fallback() -> None:
    assert detect_point_group_and_symmetry_number([0], [[0.0, 0.0, 0.0]]) == (
        "C1",
        1,
        "auto-fallback",
    )


def test_detector_failure_is_recorded_as_fallback(monkeypatch) -> None:
    import pymsym

    class BrokenContext:
        def __init__(self, *args, **kwargs):
            raise RuntimeError("detector failed")

    monkeypatch.setattr(pymsym, "Context", BrokenContext)
    assert detect_point_group_and_symmetry_number(
        [1, 1], [[0.0, 0.0, -0.37], [0.0, 0.0, 0.37]]
    ) == ("C1", 1, "auto-fallback")


@pytest.mark.parametrize(
    ("point_group", "expected"),
    [("T", 12), ("O", 24), ("I", 60), ("S6", 3), ("D3h", 6), ("C4v", 4)],
)
def test_point_group_uses_external_rotational_symmetry_number(
    point_group: str, expected: int
) -> None:
    assert symmetry_number_from_point_group(point_group) == expected


@pytest.mark.parametrize("point_group", ["C0", "C0v", "D0", "D0h", "S0"])
def test_zero_order_point_group_is_rejected(point_group: str) -> None:
    with pytest.raises(ValueError, match="order must be positive"):
        symmetry_number_from_point_group(point_group)


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
