"""Ordinary Torch-fp32 secular solves use reliable working precision.

The NumPy reference sees exactly the same fp32-rounded input numbers, promoted
to float64. Return tensors must retain the input dtype/device; no model or
optimizer execution is involved.
"""
import numpy as np
import pytest
import torch

from pysisyphus.Geometry import Geometry
from pysisyphus.optimizers.RFOptimizer import RFOptimizer

DEVICES = ["cpu"] + (["cuda"] if torch.cuda.is_available() else [])


@pytest.fixture
def optimizer(tmp_path):
    geometry = Geometry(["H"], np.array([1., 0., 0.]), coord_type="cart")
    return RFOptimizer(geometry, hessian_init="unit", max_cycles=1, dump=False,
                       line_search=False, gdiis=False, out_dir=tmp_path)


@pytest.mark.parametrize("device", DEVICES)
@pytest.mark.parametrize("alpha", [1., 1e4])
@pytest.mark.parametrize("kind", ["min", "max"])
def test_ordinary_fp32_secular_route_and_return_contract(optimizer, device, alpha, kind):
    values = torch.tensor([-.4, .2, .7], dtype=torch.float32, device=device)
    gradient = torch.tensor([.11, -.08, .05], dtype=torch.float32, device=device)
    values64 = values.detach().cpu().numpy().astype(np.float64)
    gradient64 = gradient.detach().cpu().numpy().astype(np.float64)
    reference = optimizer.solve_rfo_secular(values64, gradient64, alpha, kind=kind)
    result = optimizer.solve_rfo_secular(values, gradient, alpha, kind=kind)
    assert reference is not None and result is not None, "Ordinary coupled inputs must retain the secular route"
    eps = torch.finfo(values.dtype).eps
    for index in (0, 3):  # Cartesian step and normalized tracking vector
        assert result[index].dtype == values.dtype and result[index].device == values.device
        actual = result[index].detach().cpu().numpy()
        assert np.isfinite(actual).all()
        np.testing.assert_allclose(actual, reference[index], rtol=8*eps, atol=0.)
    assert float(result[1]) == pytest.approx(float(reference[1]), rel=1e-12, abs=0.)
    assert np.isfinite(float(result[2])) and float(result[2]) != 0.
    np.testing.assert_allclose(float(torch.linalg.vector_norm(result[3])), 1., rtol=8*eps, atol=0.)
    np.testing.assert_array_equal(values.detach().cpu().numpy(), values64.astype(np.float32))
    np.testing.assert_array_equal(gradient.detach().cpu().numpy(), gradient64.astype(np.float32))
