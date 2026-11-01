"""RS-P-RFO combined step against the trust radius after the basis change back to Cartesians."""
import numpy as np
import pytest

from test_rfo_symmetric import rs_optimizer


def model(scale, excess):
    """Micro cycles accept a step of length R * (1 + excess) in the eigenvector basis.

    Eigenvectors of norm ``scale`` stand for the round-off of the basis change.
    """
    opt = rs_optimizer()
    eigvals = np.array([-0.2, 0.3])
    gradient = np.array([0.1, 0.12])
    eigvecs = scale * np.eye(2)
    opt.housekeeping = lambda: (0.0, gradient, opt.cur_H, eigvals, eigvecs, False)

    def accepted_step(values, grad, alpha, kind, **kwargs):
        length = opt.trust_radius * (1 + excess) if kind == "max" else 0.0
        return np.array([length]), -0.01, 1.0, np.array([0.0, 1.0])

    opt.solve_rfo_secular = accepted_step
    return opt


@pytest.mark.parametrize("roundoff", [1e-15, 1e-10, 9e-10])
def test_roundoff_beyond_the_acceptance_edge_is_scaled_to_the_trust_radius(roundoff):
    # 1e-12 is the allowance of the micro-cycle check, so this step is accepted.
    opt = model(1.0 + roundoff, 1e-12)
    step = opt.optimize()
    assert np.linalg.norm(step) == pytest.approx(opt.trust_radius, rel=1e-14)
    assert step[1] == 0.0 and step[0] > 0.0
    assert len(opt.predicted_energy_changes) == 1


@pytest.mark.parametrize("roundoff", [1.1e-9, 1e-6])
def test_combined_step_beyond_1e9_of_the_trust_radius_still_raises(roundoff):
    opt = model(1.0 + roundoff, 1e-12)
    with pytest.raises(ValueError, match="combined step exceeds the trust radius"):
        opt.optimize()


def test_step_inside_the_trust_radius_is_not_rescaled():
    opt = model(1.0, -0.5)
    step = opt.optimize()
    assert np.linalg.norm(step) == pytest.approx(0.5 * opt.trust_radius, rel=1e-14)
