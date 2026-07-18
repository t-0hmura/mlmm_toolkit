"""M55: the opt-in IRC endpoint guard must decide on an EXACT Hessian at the
exact current coordinates, not on the integrator's quasi-Newton ``mw_hessian``.

The guard is exercised at the method level with a stub geometry — the same
"direct current-source probe" shape the audit used — because a full EulerPC run
needs a real MLIP calculator.  ``IRC.__new__`` bypasses ``__init__`` so the
predicate/helper methods can be driven in isolation.
"""

from __future__ import annotations

import numpy as np

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

    @property
    def cart_hessian(self):
        self.hessian_calls += 1
        self.seen_coords.append(self._coords.copy())
        return self._exact_hessian


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
