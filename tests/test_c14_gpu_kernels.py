"""C14 numerical-parity falsifiers for the measured GPU kernels.

Each adopted candidate is pinned to its baseline to ordinary fp64 rounding and,
where the change touches ownership, to the trial-rejection / accepted-state
rollback invariant. Tests run on CPU tensors by default; CUDA variants are added
when a device is present. No MLIP model is needed.
"""
from __future__ import annotations

import inspect

import numpy as np
import pytest
import torch

from pysisyphus.optimizers import hessian_updates as HU
from pysisyphus.optimizers.RFOptimizer import RFOptimizer
from pysisyphus.optimizers.hessian_updates import (
    bofill_update,
    bofill_rank2_factors,
    _bofill_update_cpu_offload,
    bofill_cpu_offload_enabled,
)
from pysisyphus._array import active_square


DEVICES = ["cpu"] + (["cuda"] if torch.cuda.is_available() else [])
DTYPES = [torch.float64, torch.float32]


def _tol(dt):
    # dtype-specific tolerances; do not loosen to pass a faster kernel.
    return dict(rtol=1e-12, atol=1e-12) if dt == torch.float64 else dict(rtol=2e-5, atol=2e-6)


def _dense_bofill_ref_np(H, dx, dg):
    """Independent NumPy reference for mix*SR1 + (1-mix)*PSB."""
    H = np.asarray(H, dtype=np.float64)
    dx = np.asarray(dx, dtype=np.float64)
    dg = np.asarray(dg, dtype=np.float64)
    z = dg - H @ dx
    zdx = z @ dx
    dxdx = dx @ dx
    sr1 = np.outer(z, z) / zdx
    psb = (np.outer(dx, z) + np.outer(z, dx)) / dxdx - zdx * np.outer(dx, dx) / (dxdx ** 2)
    mix = (zdx ** 2) / ((z @ z) * dxdx)
    return mix * sr1 + (1 - mix) * psb


def _rand_system(n, dt, dev, seed=0):
    g = torch.Generator().manual_seed(seed)
    A = torch.randn(n, n, generator=g, dtype=dt)
    H = ((A + A.t()) * 0.5).to(dev)
    dx = torch.randn(n, generator=g, dtype=dt).to(dev)
    dg = torch.randn(n, generator=g, dtype=dt).to(dev)
    return H, dx, dg


@pytest.mark.parametrize("dev", DEVICES)
@pytest.mark.parametrize("dt", DTYPES)
def test_rfo_accelerated_step_preserves_reference_device(dev, dt):
    opt = RFOptimizer.__new__(RFOptimizer)
    opt.trust_radius = 0.1
    opt.log = lambda *_: None
    ref_step = torch.tensor([0.08, 0.0], dtype=dt, device=dev)
    ip_step = torch.tensor([0.03, 0.0], dtype=dt, device=dev)

    accepted = opt._accept_accelerated_step(
        np.array([0.04, 0.0]), ip_step, ref_step,
    )

    assert isinstance(accepted, torch.Tensor)
    assert accepted.dtype == dt
    assert accepted.device == ref_step.device
    torch.testing.assert_close(
        accepted, torch.tensor([0.07, 0.0], dtype=dt, device=dev), **_tol(dt)
    )


# ---------------------------------------------------------------------------
# Candidate 1 — rank-two Bofill parity
# ---------------------------------------------------------------------------
@pytest.mark.parametrize("dev", DEVICES)
@pytest.mark.parametrize("dt", DTYPES)
def test_bofill_rank2_matches_dense_reference(dev, dt):
    H, dx, dg = _rand_system(48, dt, dev, seed=11)
    ref = _dense_bofill_ref_np(H.cpu(), dx.cpu(), dg.cpu())
    dH, key = bofill_update(H, dx, dg)
    assert key == "Bofill"
    assert dH.device.type == H.device.type
    assert dH.dtype == H.dtype
    np.testing.assert_allclose(dH.detach().cpu().numpy(), ref, **_tol(dt))
    assert torch.allclose(dH, dH.t(), **_tol(dt))


@pytest.mark.parametrize("dev", DEVICES)
@pytest.mark.parametrize("dt", DTYPES)
def test_bofill_rank2_matches_cpu_offload_fallback(dev, dt):
    H, dx, dg = _rand_system(60, dt, dev, seed=12)
    dH_r2, _ = bofill_update(H, dx, dg)
    dH_cpu, _ = _bofill_update_cpu_offload(H, dx, dg)
    np.testing.assert_allclose(
        dH_r2.detach().cpu().numpy(), dH_cpu.detach().cpu().numpy(), **_tol(dt)
    )


def test_bofill_cpu_offload_flag_dispatch(monkeypatch):
    H, dx, dg = _rand_system(24, torch.float64, "cpu", seed=13)
    monkeypatch.setenv("PYSIS_BOFILL_CPU_OFFLOAD", "1")
    assert bofill_cpu_offload_enabled() is True
    dH_fb, _ = bofill_update(H, dx, dg)
    monkeypatch.setenv("PYSIS_BOFILL_CPU_OFFLOAD", "0")
    assert bofill_cpu_offload_enabled() is False
    dH_def, _ = bofill_update(H, dx, dg)
    np.testing.assert_allclose(dH_fb.numpy(), dH_def.numpy(), rtol=1e-12, atol=1e-12)


def test_bofill_numpy_path_unchanged():
    # The numpy dense path must still agree with the reference (regression guard).
    H, dx, dg = _rand_system(32, torch.float64, "cpu", seed=14)
    Hn, dxn, dgn = H.numpy(), dx.numpy(), dg.numpy()
    dH, key = bofill_update(Hn, dxn, dgn)
    assert key == "Bofill"
    ref = _dense_bofill_ref_np(Hn, dxn, dgn)
    np.testing.assert_allclose(dH, ref, rtol=1e-12, atol=1e-12)


@pytest.mark.parametrize("dev", DEVICES)
@pytest.mark.parametrize("dt", DTYPES)
def test_bofill_degenerate_secants_are_finite(dev, dt):
    H = torch.diag(torch.arange(1, 5, dtype=dt, device=dev))
    dx = torch.tensor([1.0, -0.5, 0.25, 0.75], dtype=dt, device=dev)

    # An exact secant already satisfies dg = H dx, so no Hessian update is due.
    exact_dg = H @ dx
    for update_func in (bofill_update, _bofill_update_cpu_offload):
        dH, _ = update_func(H, dx, exact_dg)
        assert torch.isfinite(dH).all()
        torch.testing.assert_close(dH, torch.zeros_like(H), **_tol(dt))

    # z orthogonal to dx is the pure-PSB limit; it must not evaluate 0/0.
    z = torch.tensor([0.5, 1.0, 0.0, 0.0], dtype=dt, device=dev)
    assert torch.allclose(torch.dot(z, dx), torch.zeros((), dtype=dt, device=dev))
    dxdx = torch.dot(dx, dx)
    expected = (torch.outer(dx, z) + torch.outer(z, dx)) / dxdx
    for update_func in (bofill_update, _bofill_update_cpu_offload):
        dH, _ = update_func(H, dx, exact_dg + z)
        assert torch.isfinite(dH).all()
        torch.testing.assert_close(dH, expected, **_tol(dt))

    # A zero step provides no secant information and is a no-op.
    for update_func in (bofill_update, _bofill_update_cpu_offload):
        dH, _ = update_func(H, torch.zeros_like(dx), torch.ones_like(dx))
        assert torch.isfinite(dH).all()
        torch.testing.assert_close(dH, torch.zeros_like(H), **_tol(dt))


@pytest.mark.parametrize("dt", (np.float32, np.float64))
def test_bofill_numpy_degenerate_secants_are_finite(dt):
    H = np.diag(np.arange(1, 5, dtype=dt))
    dx = np.array([1.0, -0.5, 0.25, 0.75], dtype=dt)
    exact_dg = H @ dx

    exact, _ = bofill_update(H, dx, exact_dg)
    zero_step, _ = bofill_update(H, np.zeros_like(dx), np.ones_like(dx))
    z = np.array([0.5, 1.0, 0.0, 0.0], dtype=dt)
    orthogonal, _ = bofill_update(H, dx, exact_dg + z)
    expected = (np.outer(dx, z) + np.outer(z, dx)) / np.dot(dx, dx)
    assert np.isfinite(exact).all()
    assert np.isfinite(zero_step).all()
    assert np.isfinite(orthogonal).all()
    np.testing.assert_array_equal(exact, np.zeros_like(H))
    np.testing.assert_array_equal(zero_step, np.zeros_like(H))
    torch_dtype = torch.float32 if dt == np.float32 else torch.float64
    np.testing.assert_allclose(orthogonal, expected, **_tol(torch_dtype))


def test_bofill_kernel_has_no_advanced_index_addto():
    # H01 negative source test: the low-rank kernel must not use triu_indices or
    # advanced-index .add_() (the no-op MLMM pattern) or a CPU round-trip default.
    src = inspect.getsource(bofill_rank2_factors) + inspect.getsource(bofill_update)
    assert "triu_indices" not in src
    assert ".add_(" not in src
    # the default torch branch must not detach the full Hessian to CPU
    default_branch = inspect.getsource(bofill_update).split("bofill_cpu_offload_enabled")[-1]
    assert "detach().cpu()" not in default_branch


# ---------------------------------------------------------------------------
# Candidate 2 — copy-on-write ownership / trial-rejection rollback
# ---------------------------------------------------------------------------
@pytest.mark.parametrize("dev", DEVICES)
@pytest.mark.parametrize("dt", DTYPES)
def test_cow_update_does_not_mutate_accepted_hessian(dev, dt):
    H, dx, dg = _rand_system(40, dt, dev, seed=21)
    accepted = H  # the "accepted" Hessian a trial snapshot would reference
    accepted_copy = accepted.detach().clone()

    # copy-on-write update as performed by HessianOptimizer.update_hessian
    U, C = bofill_rank2_factors(accepted, dx, dg)
    new_H = accepted.clone()
    new_H.addmm_(U, C @ U.t())

    # accepted state (what rollback restores) must be byte-identical afterwards
    assert torch.equal(accepted, accepted_copy)
    # new_H must equal accepted + dense update to roundoff
    ref = _dense_bofill_ref_np(accepted_copy.cpu(), dx.cpu(), dg.cpu())
    np.testing.assert_allclose(
        new_H.detach().cpu().numpy(), accepted_copy.cpu().numpy() + ref, **_tol(dt)
    )
    # and new_H is a distinct object (replacement, not in-place) -> rollback safe
    assert new_H.data_ptr() != accepted.data_ptr()


def test_cow_rollback_restores_prior_state():
    # Emulate accept-then-reject: replacement semantics let the old reference win.
    H, dx, dg = _rand_system(30, torch.float64, "cpu", seed=22)
    snapshot = H  # the H reference the trial snapshot holds and _restore_hessian_trial_state puts back
    before = snapshot.detach().clone()
    U, C = bofill_rank2_factors(H, dx, dg)
    trial = H.clone()
    trial.addmm_(U, C @ U.t())
    # reject: restore self.H = snapshot
    restored = snapshot
    assert torch.equal(restored, before)
    assert not torch.equal(trial, before)


# ---------------------------------------------------------------------------
# Candidate 3 — vector mass scaling parity
# ---------------------------------------------------------------------------
@pytest.mark.parametrize("dev", DEVICES)
@pytest.mark.parametrize("dt", DTYPES)
def test_vector_mass_scaling_matches_dense(dev, dt):
    g = torch.Generator().manual_seed(31)
    n = 36
    A = torch.randn(n, n, generator=g, dtype=dt)
    H = ((A + A.t()) * 0.5).to(dev)
    masses = (torch.rand(n, generator=g, dtype=dt).abs() + 0.5).to(dev)
    d = 1.0 / masses.sqrt()
    vec = d.unsqueeze(1) * H * d.unsqueeze(0)
    # fp64-exact reference for D @ H @ D (avoids fp32 TF32 matmul confounds; the
    # vector path is exact elementwise scaling, so it must match the true
    # product to storage-dtype rounding).
    Hd, dd = H.double(), d.double()
    ref = (torch.diag(dd) @ Hd @ torch.diag(dd)).to(dt)
    assert torch.allclose(vec, ref, **_tol(dt))
    # on CPU (no TF32) the plain dense product is bit-identical for a diagonal D
    if dev == "cpu":
        assert torch.equal(vec, torch.diag(d) @ H @ torch.diag(d))
    # mode-matrix row scaling
    M = torch.randn(n, 5, generator=g, dtype=dt).to(dev)
    ref_M = (torch.diag(dd) @ M.double()).to(dt)
    assert torch.allclose(d.unsqueeze(1) * M, ref_M, **_tol(dt))


def test_irc_mw_hessian_active_vector_path():
    # Drive IRC._mw_hessian_active directly with a diagonal mm_inv2.
    from pysisyphus.irc.IRC import IRC

    irc = IRC.__new__(IRC)
    g = torch.Generator().manual_seed(41)
    n = 12
    A = torch.randn(n, n, generator=g, dtype=torch.float64)
    H = (A + A.t()) * 0.5
    d = torch.rand(n, generator=g, dtype=torch.float64).abs() + 0.5
    irc.mm_inv2 = torch.diag(d)
    out = irc._mw_hessian_active(H)
    ref = torch.diag(d) @ H @ torch.diag(d)
    assert torch.allclose(out, ref, rtol=1e-12, atol=1e-12)
    irc.mm_inv2 = np.diag(d.numpy())
    out_np = irc._mw_hessian_active(H.numpy())
    np.testing.assert_allclose(out_np, ref.numpy(), rtol=1e-12, atol=1e-12)


# ---------------------------------------------------------------------------
# Candidate 4 — bounded row-chunk active square extraction parity
# ---------------------------------------------------------------------------
@pytest.mark.parametrize("dev", DEVICES)
@pytest.mark.parametrize("dt", [torch.float64, torch.float32])
@pytest.mark.parametrize("n,m", [(64, 40), (64, 64), (100, 1), (50, 0)])
def test_active_square_matches_chained(dev, dt, n, m):
    g = torch.Generator().manual_seed(51)
    H = torch.randn(n, n, generator=g, dtype=dt).to(dev)
    if m == 0:
        idx = torch.empty(0, dtype=torch.long, device=dev)
        out = active_square(H, idx)
        assert out.shape == (0, 0)
        return
    idx = torch.randperm(n, generator=g)[:m].to(dev)
    idx_sorted = idx.sort().values
    ref = H.index_select(0, idx_sorted).index_select(1, idx_sorted)
    got = active_square(H, idx_sorted)
    assert torch.equal(got, ref)
    # unordered / noncontiguous index set must also match the chained gather
    ref_u = H.index_select(0, idx).index_select(1, idx)
    got_u = active_square(H, idx)
    assert torch.equal(got_u, ref_u)
    # tiny row-chunk budget forces multiple chunks -> still identical
    got_tiny = active_square(H, idx, row_chunk_bytes=1)
    assert torch.equal(got_tiny, ref_u)


def test_active_square_numpy_matches_ix():
    g = torch.Generator().manual_seed(52)
    H = torch.randn(40, 40, generator=g, dtype=torch.float64).numpy()
    idx = np.array([3, 0, 17, 39, 5])
    got = active_square(H, idx)
    ref = H[np.ix_(idx, idx)]
    np.testing.assert_array_equal(got, ref)


# ---------------------------------------------------------------------------
# Candidate 5 — device-native FD force provider parity (synthetic provider)
# ---------------------------------------------------------------------------
@pytest.mark.parametrize("dt", DTYPES)
def test_device_fd_force_path_matches_numpy_roundtrip(dt):
    dev = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    n_atoms = 8
    dof = n_atoms * 3
    gg = torch.Generator().manual_seed(61)
    M = torch.randn(dof, dof, generator=gg, dtype=dt)
    K = ((M + M.t()) * 0.5).to(dev)
    x0 = torch.randn(dof, generator=gg, dtype=dt).to(dev)
    coord = np.zeros((n_atoms, 3))
    eps = 1e-3

    def forces_tensor(c):
        x = torch.as_tensor(c.reshape(-1), dtype=K.dtype, device=dev)
        return (-(K @ (x - x0))).reshape(-1, 3).detach()

    def assemble(native):
        # Mirror mlmm hessian_fd: both branches cast each displacement force to the
        # Hessian dtype/device before the central difference, so the assembly
        # device is held constant and this isolates the force SOURCE change
        # (device-resident tensor vs NumPy round-trip) from the compute device.
        H = torch.zeros((dof, dof), device=dev, dtype=dt)
        cp = coord.copy(); cm = coord.copy()
        for k in range(dof):
            a, c = k // 3, k % 3
            cp[a, c] = coord[a, c] + eps; cm[a, c] = coord[a, c] - eps
            if native:
                Fp = forces_tensor(cp).reshape(-1).to(dev, dtype=dt)
                Fm = forces_tensor(cm).reshape(-1).to(dev, dtype=dt)
            else:
                Fp = torch.from_numpy(forces_tensor(cp).reshape(-1).cpu().numpy()).to(dev, dtype=dt)
                Fm = torch.from_numpy(forces_tensor(cm).reshape(-1).cpu().numpy()).to(dev, dtype=dt)
            H[:, k] = -(Fp - Fm) / (2 * eps)
            cp[a, c] = coord[a, c]; cm[a, c] = coord[a, c]
        return H

    Hn = assemble(True)
    Hy = assemble(False)
    assert torch.equal(Hn, Hy)


def test_uma_backend_exposes_forces_tensor():
    # The native-tensor FD port must exist on the in-process UMA backend
    # (additive capability); mlmm's ML backends drive predictions through
    # ``eval``/``forces_tensor``.
    from mlmm.backends import mlmm_calc as mm

    backend_cls = None
    for _, obj in inspect.getmembers(mm, inspect.isclass):
        if hasattr(obj, "forces_tensor") and hasattr(obj, "eval"):
            backend_cls = obj
            break
    assert backend_cls is not None, "a UMA backend must expose forces_tensor for device-native FD"


def test_hessian_fd_gates_on_forces_tensor():
    # The generic FD assembler must take the device-native branch only when the
    # backend exposes ``forces_tensor`` and is not in the parallel-predict mode,
    # otherwise fall back to the NumPy round-trip. Guard the gate wiring so a
    # refactor cannot silently drop the fast path or force it on ASE/MM backends.
    from mlmm.backends import mlmm_calc as mm

    src = inspect.getsource(mm._MLBackend.hessian_fd)
    assert "forces_tensor" in src
    assert "parallel_predict" in src
