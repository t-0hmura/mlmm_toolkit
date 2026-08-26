"""Regression tests for unbiased scan energies and ML/MM provenance."""

from __future__ import annotations

import numpy as np

from mlmm.core.utils import calculator_provenance, unbiased_energy_hartree


class _InternalCoordinateGeometry:
    atoms = ("H", "H")
    coords = np.array([99.0])
    coords3d = np.array([[0.0, 0.0, 0.0], [1.4, 0.0, 0.0]])


class _BareCalculator:
    def __init__(self) -> None:
        self.received = None

    def get_energy(self, atoms, coords):
        self.received = np.asarray(coords).copy()
        return {"energy": 1.2345}


def test_unbiased_energy_uses_cartesian_coordinates_for_internal_geometry() -> None:
    calculator = _BareCalculator()

    energy = unbiased_energy_hartree(_InternalCoordinateGeometry(), calculator)

    assert energy == 1.2345
    np.testing.assert_array_equal(
        calculator.received, _InternalCoordinateGeometry.coords3d
    )


def test_calculator_provenance_resolves_backend_specific_model() -> None:
    provenance = calculator_provenance(
        {
            "backend": "orb",
            "orb_model": "orb-test",
            "orb_precision": "float64",
            "mm_backend": "openmm",
            "link_atom_method": "fixed",
            "use_cmap": True,
        }
    )

    assert provenance == {
        "mlip_backend": "orb",
        "mlip_model": "orb-test",
        "mlip_model_label": "ORB-test",
        "mlip_task": None,
        "mlip_precision": "fp64",
        "mm_backend": "openmm",
        "link_atom_method": "fixed",
        "use_cmap": True,
    }


def test_calculator_provenance_labels_custom_factory() -> None:
    provenance = calculator_provenance(
        {
            "backend": "custom",
            "calc_file": "/tmp/my_calc.py",
            "calc_factory": "build",
        }
    )
    assert provenance["mlip_backend"] == "custom"
    assert provenance["mlip_model"] == "my_calc.py:build"
    assert provenance["mlip_precision"] is None
