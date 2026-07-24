"""M55: the opt-in IRC endpoint guard must decide on an EXACT Hessian at the
exact current coordinates, not on the integrator's quasi-Newton ``mw_hessian``.

The guard is exercised at the method level with a stub geometry — the same
"direct current-source probe" shape the audit used — because a full EulerPC run
needs a real MLIP calculator.  ``IRC.__new__`` bypasses ``__init__`` so the
predicate/helper methods can be driven in isolation.
"""

from __future__ import annotations

import numpy as np

from pysisyphus.Geometry import Geometry
from pysisyphus.irc.IRC import IRC


class _StubGeom:
    """Geometry whose exact Hessian is fixed and whose reads are counted."""

    is_analytical_2d = True  # bypasses the rigid TR projection in _project_active

    def __init__(self, coords, exact_hessian):
        self._coords = np.asarray(coords, dtype=float)
        self._exact_hessian = exact_hessian
        self.hessian_calls = 0
        self.seen_coords = []
        self._hessian = "quasi-newton-placeholder"

    @property
    def cart_coords(self):
        return self._coords

    def get_energy_and_cart_hessian_at(self, coords):
        self.hessian_calls += 1
        self.seen_coords.append(np.asarray(coords, dtype=float).copy())
        return {"hessian": self._exact_hessian}


def _guard(coords, exact_hessian, quasi_hessian):
    irc = IRC.__new__(IRC)
    irc.geometry = _StubGeom(coords, exact_hessian)
    irc._act_dofs = np.arange(6)
    irc._act_atoms = np.arange(2)
    irc.mm_inv2 = np.eye(6)
    irc.mw_hessian = quasi_hessian
    return irc


def test_positive_quasi_newton_with_negative_exact_blocks_convergence() -> None:
    irc = _guard(np.zeros(6), -np.eye(6), np.eye(6))

    # The legacy predicate would pass on the positive quasi-Newton matrix.
    assert irc._mw_hessian_is_pos_def() is True
    # The exact guard consults the negative exact Hessian and rejects.
    assert irc._exact_endpoint_is_pos_def() is False
    assert irc.geometry.hessian_calls == 1


def test_negative_quasi_newton_with_positive_exact_allows_convergence() -> None:
    irc = _guard(np.zeros(6), 2.0 * np.eye(6), -np.eye(6))

    assert irc._mw_hessian_is_pos_def() is False
    assert irc._exact_endpoint_is_pos_def() is True
    assert irc.geometry.hessian_calls == 1


def test_exact_hessian_requested_at_convergence_coordinates() -> None:
    coords = np.array([0.1, 0.2, 0.3, 0.4, 0.5, 0.6])
    irc = _guard(coords, 2.0 * np.eye(6), np.eye(6))

    irc._exact_endpoint_is_pos_def()
    np.testing.assert_array_equal(irc.geometry.seen_coords[0], coords)


def test_inconsistent_exact_shape_fails_closed() -> None:
    # n_full = n_act = 6, but the exact Hessian is 4x4 -> reject (not pos def).
    irc = _guard(np.zeros(6), np.eye(4), np.eye(6))
    assert irc._exact_endpoint_is_pos_def() is False


def test_legacy_predicate_never_evaluates_the_exact_hessian() -> None:
    # With require_pos_def_hessian disabled the run loop only ever consults the
    # quasi-Newton matrix; the exact getter must stay untouched.
    irc = _guard(np.zeros(6), -np.eye(6), np.eye(6))
    irc._mw_hessian_is_pos_def()
    assert irc.geometry.hessian_calls == 0


def test_coordinate_change_requests_a_fresh_exact_hessian() -> None:
    irc = _guard(np.zeros(6), 2.0 * np.eye(6), np.eye(6))

    irc._exact_endpoint_is_pos_def()
    irc._exact_endpoint_is_pos_def()  # same coords -> cached, no new eval
    assert irc.geometry.hessian_calls == 1

    irc.geometry._coords = np.ones(6)  # a later IRC step moved the endpoint
    irc._exact_endpoint_is_pos_def()
    assert irc.geometry.hessian_calls == 2


def test_exact_evaluation_does_not_mutate_real_geometry_result_state() -> None:
    class _Calculator:
        def __init__(self, exact_hessian):
            self.exact_hessian = exact_hessian
            self.seen_coords = []

        def get_hessian(self, _atoms, coords):
            self.seen_coords.append(np.asarray(coords, dtype=float).copy())
            return {
                "energy": 12.0,
                "forces": np.ones(6),
                "hessian": self.exact_hessian,
                "within_partial_hessian": {
                    "active_n_dof": 6,
                    "full_n_dof": 6,
                },
            }

    coords = np.arange(6, dtype=float) / 10.0
    exact_hessian = 3.0 * np.eye(6)
    calculator = _Calculator(exact_hessian)
    geom = Geometry(["H", "H"], coords, coord_type="cart")
    geom.set_calculator(calculator)

    prior_hessian = np.eye(6)
    prior_results = {"energy": -1.0, "tag": "accepted"}
    prior_partial = {"active_n_dof": 3, "full_n_dof": 6}
    prior_forces = np.arange(6, dtype=float)
    prior_true_forces = np.arange(6, dtype=float) + 10.0
    prior_true_hessian = 2.0 * np.eye(6)
    prior_all_energies = [-2.0, -1.0]
    geom._hessian = prior_hessian
    geom.results = prior_results
    geom.within_partial_hessian = prior_partial
    geom._energy = -1.0
    geom._forces = prior_forces
    geom.true_energy = -2.0
    geom.true_forces = prior_true_forces
    geom.true_hessian = prior_true_hessian
    geom._all_energies = prior_all_energies

    irc = IRC.__new__(IRC)
    irc.geometry = geom
    result = irc._exact_cart_hessian_at_current_coords()

    np.testing.assert_array_equal(result, exact_hessian)
    np.testing.assert_array_equal(calculator.seen_coords[0], coords)
    assert geom._hessian is prior_hessian
    assert geom.results is prior_results
    assert geom.within_partial_hessian is prior_partial
    assert geom._energy == -1.0
    assert geom._forces is prior_forces
    assert geom.true_energy == -2.0
    assert geom.true_forces is prior_true_forces
    assert geom.true_hessian is prior_true_hessian
    assert geom._all_energies is prior_all_energies
