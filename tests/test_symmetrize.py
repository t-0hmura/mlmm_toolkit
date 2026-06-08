"""Regression tests for mlmm.core.utils.symmetrize_inplace.

Guards the chunked diagonal-block average against the strided
self-overlapping-view RuntimeError that aborted freq / TSopt mid-run when the
Hessian edge length straddled a chunk boundary (e.g. N=513 with chunk=512,
n_hess_active=171). The in-place `diag.add_(diag_tmp)` on a strided self-view
raised "some elements of the input tensor and the written-to tensor refer to a
single memory location"; the fix is an out-of-place add-then-assign.
"""

import pytest
import torch

from mlmm.core.utils import symmetrize_inplace

# N=512/513 bracket the default chunk=512 boundary; 170/171 mirror the
# observed n_hess_active that reproduced the crash; 1024 exercises a second
# full chunk; 1 is the degenerate scalar case.
_NS = [1, 170, 171, 512, 513, 600, 1024]
_DTYPES = [torch.float32, torch.float64]


@pytest.mark.parametrize("dtype", _DTYPES)
@pytest.mark.parametrize("N", _NS)
def test_symmetrize_inplace_matches_reference(N, dtype):
    """symmetrize_inplace(A) == 0.5 * (A + A.T) for every (N, dtype)."""
    torch.manual_seed(0)
    A = torch.randn(N, N, dtype=dtype)
    ref = (A + A.t()).mul(0.5)
    out = symmetrize_inplace(A.clone())
    assert out.shape == (N, N)
    assert out.dtype == dtype
    # fp32/fp64 averaging is exact here (sum of two values then *0.5), so an
    # exact-ish tolerance is appropriate; allow tiny rounding slack only.
    assert torch.allclose(out, ref, rtol=0.0, atol=1e-6 if dtype == torch.float32 else 1e-12)
    # Result must be symmetric.
    assert torch.allclose(out, out.t(), rtol=0.0, atol=1e-6 if dtype == torch.float32 else 1e-12)


@pytest.mark.parametrize("dtype", _DTYPES)
def test_symmetrize_inplace_n513_no_raise(dtype):
    """The N=513 strided self-overlapping-view case must not raise RuntimeError."""
    torch.manual_seed(0)
    A = torch.randn(513, 513, dtype=dtype)
    # Previously raised: "some elements of the input tensor and the written-to
    # tensor refer to a single memory location".
    out = symmetrize_inplace(A.clone())
    ref = (A + A.t()).mul(0.5)
    assert torch.allclose(out, ref, rtol=0.0, atol=1e-6 if dtype == torch.float32 else 1e-12)
