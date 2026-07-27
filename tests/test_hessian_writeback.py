"""CHEMISTRY-RULE:7 — Hessian updates must be written back, not added in place to a copy.

`H[tensor_idx]` and `H[scalar, :, tensor_idx, :]` are *advanced* indexing: PyTorch returns a
COPY, so `H[...].add_(x)` mutates a temporary and the contribution is silently lost. Both
Hessian-assembly sites in this package once did exactly that:

* `mlmm/workflows/tsopt.py::_bofill_update_active` — the whole Bofill update was a no-op
* `mlmm/backends/mlmm_calc.py` — the link-atom <-> ML-region coupling blocks were discarded

The behavioural tests below pin the Bofill helper, whose inputs are small enough to construct
directly. The coupling blocks in `mlmm_calc.py` need a real ML/MM system with link atoms, so they
are NOT covered here and are not covered by the release smoke either — `tests/smoke/
backend_analytical_hessian.py` compares a bare two-atom `Atoms("H2")` Hessian between backends and
never exercises a link boundary. That gap is recorded as an open item; do not read this file as
coverage for it.
"""

from __future__ import annotations

import numpy as np
import torch

from mlmm.workflows.tsopt import _bofill_update_active


def _system(n: int = 9, seed: int = 3):
    """A representative non-degenerate (H, step, gradient-pair) triple."""
    g = torch.Generator().manual_seed(seed)
    A = torch.randn(n, n, generator=g, dtype=torch.float64)
    H = (A + A.T) * 0.5
    delta = np.linspace(0.01, 0.09, n)
    g_old = np.linspace(-0.04, 0.05, n)
    g_new = g_old + np.linspace(0.006, 0.03, n)
    return H, delta, g_old, g_new


def test_bofill_update_active_mutates_the_hessian_in_place():
    """The pre-fix `.add_()` form left H bit-identical; that must never pass again."""
    H, delta, g_old, g_new = _system()
    H0 = H.clone()

    returned = _bofill_update_active(H, delta, g_new, g_old)

    assert returned is H, "the helper documents an in-place update on the passed tensor"
    assert not torch.allclose(H, H0), (
        "Bofill update was discarded — advanced indexing returned a copy "
        "(write back with `H[idx] = H[idx] + inc`, never `H[idx].add_(inc)`)"
    )


def test_bofill_update_active_satisfies_the_secant_condition():
    """Independent invariant: both the SR1 and PSB parts of Bofill map d onto xi exactly,
    for any mixing factor, so (H_new - H_old) @ d == y - H_old @ d."""
    H, delta, g_old, g_new = _system()
    H0 = H.clone()

    d = torch.as_tensor(delta, dtype=H.dtype)
    y = torch.as_tensor(g_new - g_old, dtype=H.dtype)
    xi = y - H0 @ d
    # guard the guard: stay off the eps fallback branches so this exercises the real formula
    assert abs(float(torch.dot(d, xi))) > 1e-8
    assert float(torch.dot(d, d)) > 1e-8

    _bofill_update_active(H, delta, g_new, g_old)
    dH = H - H0

    assert torch.allclose(dH @ d, xi, rtol=1e-10, atol=1e-12)
    assert torch.allclose(dH, dH.T, rtol=1e-12, atol=1e-14), "the update must stay symmetric"


def test_bofill_update_active_is_not_the_identity_on_the_off_diagonal():
    """The upper-triangle path mirrors into the lower triangle; a write-back regression that
    only reached the diagonal would still fail here."""
    H, delta, g_old, g_new = _system()
    H0 = H.clone()

    _bofill_update_active(H, delta, g_new, g_old)
    dH = H - H0
    off = dH - torch.diag(torch.diagonal(dH))

    assert float(off.abs().max()) > 1e-12


# No source-level guard here on purpose. A text rule like "no `].add_(` in these files" also
# rejects the legitimate basic-index (view) form, and it misses the same defect written across a
# line break or through an alias — it constrains style without testing behaviour.
