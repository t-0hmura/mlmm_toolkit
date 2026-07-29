"""Ordered active-atom validation shared by both public entry points.

An invalid active list must fail with one typed ``ValueError`` — the same error
vocabulary from the direct ``build_analytical_hessian`` path and the public
``workflows.torch_hessian`` path — and a valid list's order is preserved exactly.
"""

from __future__ import annotations

import numpy as np
import pytest
import torch

from hessian_ff.active_space import validate_active_atoms
from hessian_ff.analytical_hessian import build_analytical_hessian
from hessian_ff.workflows import _normalize_active_atoms


_INVALID = [
    [],                 # empty
    [1, 1],             # duplicate
    [-1],               # negative
    [3],                # == natom
    [5],                # > natom
    [1.0],              # float truncation
    [True],             # bool
    ["1"],              # numeric string
    [None],             # non-index
]


@pytest.mark.parametrize("bad", _INVALID)
def test_validate_active_atoms_rejects_invalid(bad) -> None:
    with pytest.raises(ValueError):
        validate_active_atoms(3, bad)


def test_validate_active_atoms_preserves_order() -> None:
    assert validate_active_atoms(5, [3, 1]) == [3, 1]
    assert validate_active_atoms(5, [4, 0, 2]) == [4, 0, 2]


def test_validate_active_atoms_accepts_numpy_integers() -> None:
    assert validate_active_atoms(5, [np.int64(3), np.int64(1)]) == [3, 1]


@pytest.mark.parametrize("bad", _INVALID)
def test_workflows_normalizer_uses_the_shared_validator(bad) -> None:
    # workflows.torch_hessian routes through _normalize_active_atoms; it must
    # raise the same typed error (no silent duplicate-drop) as the validator.
    with pytest.raises(ValueError):
        _normalize_active_atoms(3, bad)


class _StubSystem:
    """Minimal stand-in exposing only ``natom`` for the validation-first path."""

    natom = 3

    def to(self, *args, **kwargs):  # pragma: no cover - never reached
        return self


@pytest.mark.parametrize("bad", [[], [1, 1], [-1], [3], [1.0], [True]])
def test_build_analytical_hessian_validates_before_native(bad) -> None:
    # The direct build path validates *before* any native probing / allocation,
    # so the same invalid lists raise the same typed error as the workflow path.
    with pytest.raises(ValueError):
        build_analytical_hessian(_StubSystem(), torch.zeros((3, 3)), bad)
