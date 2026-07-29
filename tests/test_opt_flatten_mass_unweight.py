"""The ordinary-opt flatten displacement is mass-unweighted
exactly once.

``mlmm.workflows.opt._flatten_all_imag_modes_for_geom`` turns a mass-weighted
eigenvector q into the Cartesian normal-mode direction
u = M^(-1/2) q / ||M^(-1/2) q|| and displaces by ``flatten_amp / BOHR2ANG``
along u.  A second per-atom mass factor (the old ``mass_scale``) would rotate
the displacement toward M^(-1) q and change its amplitude.  This test binds to
the production helper (no reimplemented gate logic) and would fail under either
a double- or a missing mass-unweight.
"""

from __future__ import annotations

import numpy as np
import torch

from pysisyphus.constants import BOHR2ANG
import mlmm.workflows.opt as opt_mod
from mlmm.workflows.opt import _flatten_all_imag_modes_for_geom


class _Geometry:
    """Minimal geometry exposing the cart_coords / coords surface the flatten
    helper drives (mirrors the production Geometry contract used here)."""

    def __init__(self) -> None:
        self.cart_coords = np.zeros(6, dtype=float)

    @property
    def coords(self):
        return self.cart_coords

    @coords.setter
    def coords(self, value) -> None:
        self.cart_coords = np.asarray(value, dtype=float).copy()


# H then C so an erroneous second mass factor visibly rotates a mixed mode.
_MASSES = np.array([1.008, 12.011])
_AMP_ANG = 0.10

# Correct single-unweight ratio |u_H / u_C| = sqrt(m_C / m_H).
_EXPECTED_RATIO = float(np.sqrt(12.011 / 1.008))          # 3.451908834714...
# Double-scaled (buggy) ratio = m_C / m_H.
_DOUBLE_SCALED_RATIO = float(12.011 / 1.008)              # 11.915674603175...


def _run_flatten(monkeypatch, energies):
    geom = _Geometry()
    # One imaginary mode with equal-magnitude H-x and C-x components.
    modes = torch.zeros((1, 6), dtype=torch.float64)
    modes[0, 0] = 1.0  # H x
    modes[0, 3] = 1.0  # C x
    freqs_cm = np.array([-100.0])

    seq = iter(energies)
    monkeypatch.setattr(
        opt_mod,
        "_calc_energy",
        lambda geometry, kwargs, calc=None: next(seq),
    )

    flattened = _flatten_all_imag_modes_for_geom(
        geom,
        _MASSES,
        {},
        freqs_cm,
        modes,
        neg_freq_thresh_cm=5.0,
        flatten_amp_ang=_AMP_ANG,
        calculator=None,
    )
    assert flattened is True
    return geom


def _expected_unit():
    u = np.array([1.0, 0.0, 0.0, 1.0, 0.0, 0.0]) / np.sqrt(np.repeat(_MASSES, 3))
    return u / np.linalg.norm(u)


def test_opt_flatten_plus_branch_is_single_mass_unweight(monkeypatch) -> None:
    # E_ref, E_plus, E_minus -> keep the plus branch (E_plus <= E_minus).
    geom = _run_flatten(monkeypatch, (0.0, -1.0, -0.5))

    disp = geom.cart_coords.copy()
    # (6) amplitude preserved: ||delta_x|| == flatten_amp / BOHR2ANG.
    np.testing.assert_allclose(
        np.linalg.norm(disp), _AMP_ANG / BOHR2ANG, atol=1.0e-12
    )
    # (5) direction equals M^(-1/2) q normalized.
    np.testing.assert_allclose(
        disp / np.linalg.norm(disp), _expected_unit(), atol=1.0e-12
    )
    # (7) H/C component ratio is the single-unweight value, not double-scaled.
    ratio = abs(disp[0] / disp[3])
    np.testing.assert_allclose(ratio, _EXPECTED_RATIO, atol=1.0e-12)
    assert abs(ratio - _DOUBLE_SCALED_RATIO) > 1.0e-3


def test_opt_flatten_minus_branch_is_single_mass_unweight(monkeypatch) -> None:
    # E_ref, E_plus, E_minus -> keep the minus branch (E_minus < E_plus).
    geom = _run_flatten(monkeypatch, (0.0, -0.5, -1.0))

    disp = geom.cart_coords.copy()
    np.testing.assert_allclose(
        np.linalg.norm(disp), _AMP_ANG / BOHR2ANG, atol=1.0e-12
    )
    # Check collinearity magnitude first; the sign itself is pinned to the
    # deterministic minus branch by the assert below.
    unit = disp / np.linalg.norm(disp)
    np.testing.assert_allclose(np.abs(unit), np.abs(_expected_unit()), atol=1.0e-12)
    assert unit @ _expected_unit() < 0.0  # minus branch flips the sign
    ratio = abs(disp[0] / disp[3])
    np.testing.assert_allclose(ratio, _EXPECTED_RATIO, atol=1.0e-12)
    assert abs(ratio - _DOUBLE_SCALED_RATIO) > 1.0e-3
