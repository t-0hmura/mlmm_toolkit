"""Trust-radius growth requires a boundary-sized step at every radius."""
import numpy as np
import pytest

from pysisyphus.Geometry import Geometry
from pysisyphus.optimizers.RFOptimizer import RFOptimizer


def make_optimizer(tmp_path, radius, *, trust_min=1e-8, trust_band=False):
    # Use the actual concrete constructor, not an underinitialized __new__
    # object: this initializes logging, every trust-band field and both floors.
    geometry = Geometry(["H"], np.array([1.0, 0.0, 0.0]), coord_type="cart")
    return RFOptimizer(
        geometry, hessian_init="unit", max_cycles=1, dump=False,
        trust_radius=radius, trust_min=trust_min, trust_max=1.0,
        trust_band=trust_band, line_search=False, gdiis=False, out_dir=tmp_path,
    )


@pytest.mark.parametrize(
    "radius,step,rho,expected",
    [
        pytest.param(1e-4, 3.74316573217e-5, 1.032508995, 1e-4, id="small-interior-must-hold"),
        pytest.param(1e-6, .98e-6, 1.0, 1e-6, id="small-1e-6-outside-one-percent"),
        pytest.param(1e-5, .995e-5, 1.0, 2e-5, id="small-1e-5-inside-one-percent"),
        pytest.param(1e-3, .98e-3, 1.0, 1e-3, id="small-1e-3-outside-one-percent"),
        pytest.param(1e-2, .995e-2, 1.0, 2e-2, id="small-1e-2-inside-one-percent"),
        pytest.param(1e-4, 1e-4*(1.0+5.3e-11), 1.0, 2e-4, id="boundary-relative-roundoff"),
        pytest.param(.1, .0992, 1.0, .2, id="normal-radius-keeps-absolute-tolerance"),
        pytest.param(.3, .2992, 1.0, .6, id="large-radius-inside-absolute-tolerance"),
        pytest.param(.3, .2985, 1.0, .3, id="large-radius-absolute-cap-not-one-percent"),
        pytest.param(1e-4, 0.0, 1.0, 1e-4, id="zero-step-does-not-grow"),
    ],
)
def test_default_trust_growth_boundary(tmp_path, radius, step, rho, expected):
    optimizer = make_optimizer(tmp_path, radius)
    optimizer.set_new_trust_radius(rho, step)
    assert optimizer.trust_radius == pytest.approx(expected, rel=1e-13, abs=0.0), (
        f"radius={radius:.17g}, step={step:.17g}, rho={rho:.17g}, "
        f"measured={optimizer.trust_radius:.17g}, expected={expected:.17g}"
    )


@pytest.mark.parametrize(
    "radius,step,rho,min_radius,expected",
    [
        pytest.param(8e-4, 8e-4, .1, None, 2e-4, id="ordinary-quarter-shrink"),
        pytest.param(1e-5, 1e-5, -.5, None, 1e-5, id="below-ordinary-floor-does-not-increase"),
        pytest.param(1e-5, 1e-5, -.5, 1e-7, 2.5e-6, id="explicit-emergency-floor-shrink"),
    ],
)
def test_existing_shrink_and_emergency_floor(tmp_path, radius, step, rho, min_radius, expected):
    optimizer = make_optimizer(tmp_path, radius, trust_min=1e-4)
    optimizer.set_new_trust_radius(rho, step, min_radius=min_radius)
    assert optimizer.trust_radius == pytest.approx(expected, rel=1e-13, abs=0.0)


@pytest.mark.parametrize(
    "rho,step,expected",
    [
        pytest.param(1.0, .95e-4, 1.0925e-4, id="band-growth-still-uses-sigma-not-boundary"),
        pytest.param(10.0, .8e-4, .52e-4, id="band-shrink-still-uses-step-times-sigma"),
        pytest.param(2.0, 1e-4, 1e-4, id="band-intermediate-ratio-keeps-radius"),
    ],
)
def test_opt_in_trust_band_is_unchanged(tmp_path, rho, step, expected):
    optimizer = make_optimizer(tmp_path, 1e-4, trust_band=True)
    optimizer.set_new_trust_radius(rho, step)
    assert optimizer.trust_radius == pytest.approx(expected, rel=1e-13, abs=0.0)
