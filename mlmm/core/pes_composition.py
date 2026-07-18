"""Representation-safe composition of additive PES contributions.

Calculators expose one result mapping for energy, forces, Hessian, and optional
partial-Hessian metadata.  Restraint terms are added at that boundary so the
base mapping remains reusable and constraints are enforced after every
physical contribution.
"""

from __future__ import annotations

from typing import Any, Mapping, Optional, Sequence

import numpy as np


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


def _add_like(reference: Any, correction: Any) -> Any:
    """Add *correction* in the reference array's representation."""

    if _is_torch_tensor(reference):
        import torch

        delta = torch.as_tensor(
            correction, dtype=reference.dtype, device=reference.device
        ).reshape(reference.shape)
        return reference.clone().add(delta)

    ref = np.asarray(reference)
    delta = np.asarray(correction, dtype=ref.dtype).reshape(ref.shape)
    return ref.copy() + delta


def _zero_indices(value: Any, indices: Sequence[int], *, matrix: bool) -> Any:
    if not indices:
        return _clone_array(value)
    if _is_torch_tensor(value):
        import torch

        out = value.clone()
        idx = torch.as_tensor(indices, dtype=torch.long, device=out.device)
        if matrix:
            n_dof = int(round(out.numel() ** 0.5))
            if n_dof * n_dof != out.numel():
                raise ValueError("Hessian buffer is not square.")
            square = out.reshape(n_dof, n_dof)
            square.index_fill_(0, idx, 0.0)
            square.index_fill_(1, idx, 0.0)
        else:
            out.reshape(-1).index_fill_(0, idx, 0.0)
        return out

    out = np.asarray(value).copy()
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
        if np.asarray(force_delta_full).size != expected:
            raise ValueError(
                f"Force correction must contain {expected} values for {n_atoms} atoms."
            )
        result["forces"] = _add_like(
            base_result["forces"], np.asarray(force_delta_full).reshape(-1)
        )
    if "forces" in result and constrained_dofs:
        result["forces"] = _zero_indices(
            result["forces"], constrained_dofs, matrix=False
        )

    if hessian_delta_full is not None:
        full_n_dof = 3 * n_atoms
        delta = np.asarray(hessian_delta_full)
        if delta.size != full_n_dof * full_n_dof:
            raise ValueError(
                f"Hessian correction must have shape ({full_n_dof}, {full_n_dof})."
            )
        delta = delta.reshape(full_n_dof, full_n_dof)
        represented_dofs, is_full = _hessian_map(base_result, n_atoms=n_atoms)
        mapped = delta[np.ix_(represented_dofs, represented_dofs)]
        result["hessian"] = _add_like(base_result["hessian"], mapped)
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
            result["hessian"] = _zero_indices(
                result["hessian"], local_constrained, matrix=True
            )

    return result


__all__ = ["clone_pes_result", "compose_additive_pes_result"]
