"""Regression coverage for Orb's saved-tensor mutation during double backward."""

from __future__ import annotations

import pytest

torch = pytest.importorskip("torch")

from mlmm.backends.mlmm_calc import (  # noqa: E402
    _autograd_hessian_with_mutation_guard,
)


def _mutating_energy(coords):
    intermediate = coords.clone()
    energy = (intermediate**2).sum()
    # Mutate a value that PowBackward saved.  This reproduces the same autograd
    # version-counter failure class seen in Orb message passing.
    intermediate.sin_()
    return energy


def test_orb_hessian_wraps_forward_and_backward_in_mutation_guard() -> None:
    coords = torch.tensor([1.0, 2.0], requires_grad=True)

    with pytest.raises(RuntimeError, match="modified by an inplace operation"):
        torch.autograd.functional.hessian(
            _mutating_energy,
            coords,
            vectorize=False,
            create_graph=False,
        )

    hessian = _autograd_hessian_with_mutation_guard(_mutating_energy, coords)
    torch.testing.assert_close(hessian, 2.0 * torch.eye(2))
