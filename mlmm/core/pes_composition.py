"""Representation-safe composition of additive PES contributions.

Calculators expose one result mapping for energy, forces, Hessian, and optional
partial-Hessian metadata.  Restraint terms are added at that boundary so the
base mapping remains reusable and constraints are enforced after every
physical contribution.
"""

from __future__ import annotations

from typing import Any, Mapping, Optional, Sequence

import numpy as np
from pysisyphus._array import active_square


def _is_torch_tensor(value: Any) -> bool:
    return type(value).__module__.startswith("torch") and hasattr(value, "clone")


def _clone_array(value: Any) -> Any:
    if _is_torch_tensor(value):
        return value.clone()
    if isinstance(value, np.ndarray):
        return value.copy()
    return value


def clone_pes_result(result: Mapping[str, Any]) -> dict[str, Any]:
    """Clone top-level numerical buffers without rewriting metadata objects."""

    return {key: _clone_array(value) for key, value in result.items()}


def _add_inplace(target: Any, correction: Any) -> Any:
    """Add *correction* to a result-owned numerical buffer."""

    if _is_torch_tensor(target):
        import torch

        delta = torch.as_tensor(
            correction, dtype=target.dtype, device=target.device
        ).reshape(target.shape)
        target.add_(delta)
        return target

    out = np.asarray(target)
    if _is_torch_tensor(correction):
        correction = correction.detach().cpu().numpy()
    delta = np.asarray(correction, dtype=out.dtype).reshape(out.shape)
    np.add(out, delta, out=out)
    return target


def _zero_indices_inplace(
    value: Any, indices: Sequence[int], *, matrix: bool
) -> Any:
    if not indices:
        return value
    if _is_torch_tensor(value):
        import torch

        idx = torch.as_tensor(indices, dtype=torch.long, device=value.device)
        if matrix:
            n_dof = int(round(value.numel() ** 0.5))
            if n_dof * n_dof != value.numel():
                raise ValueError("Hessian buffer is not square.")
            square = value.reshape(n_dof, n_dof)
            square.index_fill_(0, idx, 0.0)
            square.index_fill_(1, idx, 0.0)
        else:
            value.reshape(-1).index_fill_(0, idx, 0.0)
        return value

    out = np.asarray(value)
    if matrix:
        n_dof = int(round(out.size ** 0.5))
        if n_dof * n_dof != out.size:
            raise ValueError("Hessian buffer is not square.")
        square = out.reshape(n_dof, n_dof)
        square[np.asarray(indices, dtype=int), :] = 0.0
        square[:, np.asarray(indices, dtype=int)] = 0.0
    else:
        out.reshape(-1)[np.asarray(indices, dtype=int)] = 0.0
    return out


def _array_size(value: Any) -> int:
    if _is_torch_tensor(value):
        return int(value.numel())
    return int(np.asarray(value).size)


def _mapped_hessian_delta(
    delta_full: Any,
    represented_dofs: np.ndarray,
    *,
    full_n_dof: int,
    is_full: bool,
) -> Any:
    """Return a full or compact correction without cloning the output buffer."""

    if _array_size(delta_full) != full_n_dof * full_n_dof:
        raise ValueError(
            f"Hessian correction must have shape ({full_n_dof}, {full_n_dof})."
        )
    if _is_torch_tensor(delta_full):
        square = delta_full.reshape(full_n_dof, full_n_dof)
        if is_full:
            return square
        return active_square(square, represented_dofs)

    square = np.asarray(delta_full).reshape(full_n_dof, full_n_dof)
    if is_full:
        return square
    return square[np.ix_(represented_dofs, represented_dofs)]


def _hessian_map(
    result: Mapping[str, Any], *, n_atoms: int
) -> tuple[np.ndarray, bool]:
    """Return global DOFs represented by H and whether H is full-system."""

    full_n_dof = 3 * int(n_atoms)
    hessian = result["hessian"]
    if int(np.prod(tuple(hessian.shape))) == full_n_dof * full_n_dof:
        return np.arange(full_n_dof, dtype=int), True

    metadata = result.get("within_partial_hessian")
    if not isinstance(metadata, Mapping) or "active_dofs" not in metadata:
        raise ValueError(
            "A compact Hessian requires within_partial_hessian.active_dofs "
            "to map additive PES terms."
        )
    active_dofs = np.asarray(metadata["active_dofs"], dtype=int).reshape(-1)
    if np.any(active_dofs < 0) or np.any(active_dofs >= full_n_dof):
        raise ValueError("within_partial_hessian.active_dofs contains invalid indices.")
    if np.unique(active_dofs).size != active_dofs.size:
        raise ValueError("within_partial_hessian.active_dofs contains duplicates.")
    if int(np.prod(tuple(hessian.shape))) != active_dofs.size * active_dofs.size:
        raise ValueError(
            "Hessian dimensions do not match within_partial_hessian.active_dofs."
        )
    return active_dofs, False


def _validated_atoms(
    atoms: Sequence[int], *, n_atoms: int, label: str
) -> list[int]:
    normalized = sorted({int(atom) for atom in atoms})
    if any(atom < 0 or atom >= n_atoms for atom in normalized):
        raise ValueError(
            f"{label} atom indices must be in [0, {n_atoms}); got {normalized}."
        )
    return normalized


def compose_additive_pes_result(
    base_result: Mapping[str, Any],
    *,
    n_atoms: int,
    energy_delta: Optional[float] = None,
    force_delta_full: Optional[Any] = None,
    hessian_delta_full: Optional[Any] = None,
    constrained_atoms: Sequence[int] = (),
    hessian_constrained_atoms: Optional[Sequence[int]] = None,
) -> dict[str, Any]:
    """Return ``base + delta`` with base representation and final constraints.

    Force corrections use full Cartesian order. Hessian corrections are full
    ``3N x 3N`` arrays and are mapped through the base result's
    ``within_partial_hessian.active_dofs`` metadata when needed. ML/MM may use
    a stricter Hessian-only constrained set than the force constrained set;
    other callers inherit ``constrained_atoms`` for both by default. The input
    mapping and its numerical buffers are never mutated.
    """

    n_atoms = int(n_atoms)
    constrained = _validated_atoms(
        constrained_atoms, n_atoms=n_atoms, label="Constrained"
    )
    hessian_constrained = _validated_atoms(
        constrained if hessian_constrained_atoms is None else hessian_constrained_atoms,
        n_atoms=n_atoms,
        label="Hessian-constrained",
    )

    result = clone_pes_result(base_result)
    if energy_delta is not None:
        result["energy"] = float(base_result["energy"]) + float(energy_delta)

    constrained_dofs = [
        3 * atom + axis for atom in constrained for axis in range(3)
    ]
    hessian_constrained_dofs = [
        3 * atom + axis for atom in hessian_constrained for axis in range(3)
    ]

    if force_delta_full is not None:
        expected = 3 * n_atoms
        if _array_size(force_delta_full) != expected:
            raise ValueError(
                f"Force correction must contain {expected} values for {n_atoms} atoms."
            )
        result["forces"] = _add_inplace(
            result["forces"], force_delta_full
        )
    if "forces" in result and constrained_dofs:
        result["forces"] = _zero_indices_inplace(
            result["forces"], constrained_dofs, matrix=False
        )

    if hessian_delta_full is not None:
        full_n_dof = 3 * n_atoms
        represented_dofs, is_full = _hessian_map(base_result, n_atoms=n_atoms)
        mapped = _mapped_hessian_delta(
            hessian_delta_full,
            represented_dofs,
            full_n_dof=full_n_dof,
            is_full=is_full,
        )
        result["hessian"] = _add_inplace(result["hessian"], mapped)
    elif "hessian" in result:
        represented_dofs, is_full = _hessian_map(base_result, n_atoms=n_atoms)
    else:
        represented_dofs = np.zeros(0, dtype=int)
        is_full = False

    if "hessian" in result and hessian_constrained_dofs:
        constrained_set = set(hessian_constrained_dofs)
        local_constrained = (
            hessian_constrained_dofs
            if is_full
            else [
                local
                for local, global_dof in enumerate(represented_dofs.tolist())
                if global_dof in constrained_set
            ]
        )
        if local_constrained:
            result["hessian"] = _zero_indices_inplace(
                result["hessian"], local_constrained, matrix=True
            )

    return result


__all__ = ["clone_pes_result", "compose_additive_pes_result"]
