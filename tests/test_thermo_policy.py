"""P05: the two retained thermochemistry policies are explicit and serialized.

Both entry points keep their historical numerical behaviour:

* standalone ``mlmm freq``  -> QRRHO, rotor cutoff 100 cm^-1,
  frequency scale 1, ZPE scale 1, NO imaginary inversion, NO positive floor;
* bundled ``Geometry.get_thermoanalysis`` -> the same, but inverts imaginaries
  from -15 cm^-1 and floors positive frequencies below 25 cm^-1.

The policies are named immutable objects (``WORKFLOW_THERMO_POLICY`` /
``GEOMETRY_THERMO_POLICY``); every value is passed to ``thermochemistry``
explicitly and serialized so a payload always states the effective policy.
"""

from __future__ import annotations

import warnings

import numpy as np
import pytest

from thermoanalysis.QCData import QCData
from thermoanalysis.thermo import thermochemistry
from thermoanalysis.config import (
    ThermoPolicy,
    WORKFLOW_THERMO_POLICY,
    GEOMETRY_THERMO_POLICY,
    ROTOR_CUT_DEFAULT,
)


T_K = 298.15
P_PA = 101325.0

# Water-like nonlinear fixture reused across the policy checks.
COORDS = [[0.0, 0.0, 0.0], [0.757, 0.586, 0.0], [-0.757, 0.586, 0.0]]
MASSES = [15.999, 1.008, 1.008]
SCF = -76.4

SERIALIZED_KEYS = {
    "kind",
    "rotor_cutoff_cm",
    "frequency_scale",
    "zpe_scale",
    "invert_imag_from_cm",
    "positive_frequency_floor_cm",
}


def _qc(wavenumbers):
    return QCData(
        {
            "coords3d": np.asarray(COORDS, dtype=float),
            "wavenumbers": np.asarray(wavenumbers, dtype=float),
            "scf_energy": SCF,
            "masses": np.asarray(MASSES, dtype=float),
            "mult": 1,
        },
        point_group="c1",
        mult=1,
    )


def _run(wavenumbers, **kw):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        return thermochemistry(_qc(wavenumbers), T_K, pressure=P_PA, **kw)


# ---- the named policies encode exactly the retained defaults --------------


def test_policy_field_values_are_the_retained_defaults():
    assert WORKFLOW_THERMO_POLICY == ThermoPolicy(
        kind="qrrho",
        rotor_cutoff_cm=ROTOR_CUT_DEFAULT,
        frequency_scale=1.0,
        zpe_scale=1.0,
        invert_imag_from_cm=0.0,
        positive_frequency_floor_cm=0.0,
    )
    assert GEOMETRY_THERMO_POLICY == ThermoPolicy(
        kind="qrrho",
        rotor_cutoff_cm=ROTOR_CUT_DEFAULT,
        frequency_scale=1.0,
        zpe_scale=1.0,
        invert_imag_from_cm=-15.0,
        positive_frequency_floor_cm=25.0,
    )


def test_policy_is_immutable():
    with pytest.raises(Exception):
        WORKFLOW_THERMO_POLICY.invert_imag_from_cm = -15.0  # type: ignore[misc]


def test_serialized_dict_uses_stable_field_keys():
    for policy in (WORKFLOW_THERMO_POLICY, GEOMETRY_THERMO_POLICY):
        assert set(policy.as_dict().keys()) == SERIALIZED_KEYS
    assert WORKFLOW_THERMO_POLICY.as_dict()["invert_imag_from_cm"] == 0.0
    assert WORKFLOW_THERMO_POLICY.as_dict()["positive_frequency_floor_cm"] == 0.0
    assert GEOMETRY_THERMO_POLICY.as_dict()["invert_imag_from_cm"] == -15.0
    assert GEOMETRY_THERMO_POLICY.as_dict()["positive_frequency_floor_cm"] == 25.0


# ---- identical frequencies reproduce each entry point's pre-change numbers --


def _assert_thermo_equal(a, b):
    for field in ("ZPE", "U_tot", "H", "G", "dG", "S_tot", "c_tot"):
        np.testing.assert_allclose(
            float(getattr(a, field)), float(getattr(b, field)), rtol=0.0, atol=0.0,
            err_msg=field,
        )


def test_workflow_policy_reproduces_bare_library_call():
    """WORKFLOW policy == the historical ``thermochemistry(qc, T, pressure=p)``."""
    wn = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1595.0, 3657.0, 3756.0]
    pre_change = _run(wn)  # bare call: all library defaults
    via_policy = _run(wn, **WORKFLOW_THERMO_POLICY.thermochemistry_kwargs())
    _assert_thermo_equal(via_policy, pre_change)


def test_geometry_policy_reproduces_invert15_cutoff25_call():
    """GEOMETRY policy == the historical invert_imags=-15 / cutoff=25 call."""
    wn = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1595.0, 3657.0, 3756.0]
    pre_change = _run(wn, invert_imags=-15.0, cutoff=25.0)
    via_policy = _run(wn, **GEOMETRY_THERMO_POLICY.thermochemistry_kwargs())
    _assert_thermo_equal(via_policy, pre_change)


# ---- behavioural falsifier: -10 retained by workflow, inverted by geometry --


def test_minus10_retained_by_workflow_inverted_by_geometry():
    """A -10 cm^-1 mode: workflow keeps it imaginary (dropped); geometry inverts it."""
    wn = [-10.0, 200.0]  # vibrational-only (size != 3N -> no rigid-mode drop)

    workflow = _run(wn, **WORKFLOW_THERMO_POLICY.thermochemistry_kwargs())
    # No inversion, no floor: the -10 stays imaginary and is discarded, only 200
    # survives as a real vibration.
    np.testing.assert_array_equal(
        np.asarray(workflow.wavenumbers, dtype=float), np.array([200.0])
    )

    geometry = _run(wn, **GEOMETRY_THERMO_POLICY.thermochemistry_kwargs())
    # -10 is inverted to +10, then floored to 25; 200 is untouched.
    np.testing.assert_array_equal(
        np.asarray(geometry.wavenumbers, dtype=float), np.array([25.0, 200.0])
    )


def test_positive_10_floored_only_by_geometry_policy():
    """A +10 cm^-1 mode is left at 10 by workflow but floored to 25 by geometry."""
    wn = [10.0, 200.0]

    workflow = _run(wn, **WORKFLOW_THERMO_POLICY.thermochemistry_kwargs())
    np.testing.assert_array_equal(
        np.asarray(workflow.wavenumbers, dtype=float), np.array([10.0, 200.0])
    )

    geometry = _run(wn, **GEOMETRY_THERMO_POLICY.thermochemistry_kwargs())
    np.testing.assert_array_equal(
        np.asarray(geometry.wavenumbers, dtype=float), np.array([25.0, 200.0])
    )


# ---- Geometry entry point records the effective policy on its result -------


def test_geometry_get_thermoanalysis_records_policy():
    """Geometry.get_thermoanalysis attaches the serialized geometry policy."""
    from pysisyphus.Geometry import Geometry

    atoms = ["O", "H", "H"]
    coords = np.array(COORDS, dtype=float).reshape(-1)
    geom = Geometry(atoms, coords)

    # Full 3N Hessian -> real modes; the numeric values do not affect the policy
    # record (get_thermoanalysis attaches the effective policy unconditionally).
    rng = np.random.default_rng(7)
    raw = rng.standard_normal((9, 9))
    geom.cart_hessian = raw.T @ raw + 0.5 * np.eye(9)
    geom.energy = SCF

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        geom.get_thermoanalysis()

    assert getattr(geom, "thermo_policy", None) == GEOMETRY_THERMO_POLICY.as_dict()
