"""Geometry-identified Hessian files for ``freq`` -> ``irc`` handoff."""

from __future__ import annotations

from pathlib import Path
from typing import Any, Mapping, Optional

import numpy as np


SCHEMA_VERSION = 1


def save_hessian_file(
    path: Path | str,
    *,
    hessian: Any,
    energy_ha: float,
    cart_coords_bohr: Any,
    atomic_numbers: Any,
    partial_metadata: Optional[Mapping[str, Any]] = None,
) -> Path:
    """Atomically write a Hessian at the exact requested path and return it."""
    hess = np.asarray(hessian, dtype=float)
    coords = np.asarray(cart_coords_bohr, dtype=float).reshape(-1)
    numbers = np.asarray(atomic_numbers, dtype=np.int64).reshape(-1)
    if hess.ndim != 2 or hess.shape[0] != hess.shape[1]:
        raise ValueError(f"Hessian must be square, got shape {hess.shape}.")
    if coords.size != 3 * numbers.size:
        raise ValueError("Coordinate and atomic-number sizes are inconsistent.")
    if not np.isfinite(hess).all() or not np.isfinite(coords).all():
        raise ValueError("Hessian file data must be finite.")

    payload: dict[str, Any] = {
        "schema_version": np.int64(SCHEMA_VERSION),
        "hessian": hess,
        "energy_ha": np.float64(energy_ha),
        "cart_coords_bohr": coords,
        "atomic_numbers": numbers,
    }
    if isinstance(partial_metadata, Mapping) and partial_metadata.get("active_dofs") is not None:
        active_dofs = np.asarray(partial_metadata["active_dofs"], dtype=np.int64).reshape(-1)
        payload.update({
            "wph_active_dofs": active_dofs,
            "wph_active_n_dof": np.int64(
                int(partial_metadata.get("active_n_dof", active_dofs.size))
            ),
            "wph_full_n_dof": np.int64(
                int(partial_metadata.get("full_n_dof", coords.size))
            ),
        })
    from mlmm.core.result_commit import atomic_write_exact

    destination = Path(path)
    return atomic_write_exact(
        destination,
        lambda stream: np.savez_compressed(stream, **payload),
    )


def load_hessian_file(
    path: Path | str,
    *,
    cart_coords_bohr: Any,
    atomic_numbers: Any,
    expected_active_dofs: Any = None,
    atol_bohr: float = 1.1e-3,
) -> dict[str, Any]:
    """Load and validate a Hessian against the current geometry and DOF basis."""
    current_coords = np.asarray(cart_coords_bohr, dtype=float).reshape(-1)
    current_numbers = np.asarray(atomic_numbers, dtype=np.int64).reshape(-1)
    required = {
        "schema_version",
        "hessian",
        "energy_ha",
        "cart_coords_bohr",
        "atomic_numbers",
    }
    with np.load(path, allow_pickle=False) as data:
        missing = sorted(required.difference(data.files))
        if missing:
            raise ValueError(
                "Hessian file lacks geometry identity metadata "
                f"({', '.join(missing)}); regenerate it with `mlmm freq --dump-hess`."
            )
        version = int(np.asarray(data["schema_version"]).item())
        if version != SCHEMA_VERSION:
            raise ValueError(
                f"Unsupported Hessian file schema {version}; expected {SCHEMA_VERSION}."
            )
        saved_coords = np.asarray(data["cart_coords_bohr"], dtype=float).reshape(-1)
        saved_numbers = np.asarray(data["atomic_numbers"], dtype=np.int64).reshape(-1)
        if not np.array_equal(saved_numbers, current_numbers):
            raise ValueError("Hessian file atomic numbers/order do not match the IRC start geometry.")
        if saved_coords.shape != current_coords.shape or not np.allclose(
            saved_coords, current_coords, atol=float(atol_bohr), rtol=0.0
        ):
            raise ValueError("Hessian file coordinates do not match the IRC start geometry.")

        hess = np.asarray(data["hessian"], dtype=float)
        if hess.ndim != 2 or hess.shape[0] != hess.shape[1] or not np.isfinite(hess).all():
            raise ValueError(f"Hessian file contains an invalid matrix with shape {hess.shape}.")
        full_n_dof = int(current_coords.size)
        partial: Optional[dict[str, Any]] = None
        saved_active_dofs = np.arange(full_n_dof, dtype=np.int64)
        if "wph_active_dofs" in data.files:
            active_dofs = np.asarray(data["wph_active_dofs"], dtype=np.int64).reshape(-1)
            active_n_dof = int(data["wph_active_n_dof"]) if "wph_active_n_dof" in data.files else int(active_dofs.size)
            saved_full_n_dof = int(data["wph_full_n_dof"]) if "wph_full_n_dof" in data.files else full_n_dof
            if (
                active_n_dof != active_dofs.size
                or hess.shape[0] != active_dofs.size
                or saved_full_n_dof != full_n_dof
                or np.unique(active_dofs).size != active_dofs.size
                or np.any(active_dofs < 0)
                or np.any(active_dofs >= full_n_dof)
            ):
                raise ValueError("Hessian file contains inconsistent active-DOF metadata.")
            partial = {
                "active_n_dof": active_n_dof,
                "full_n_dof": saved_full_n_dof,
                "active_dofs": active_dofs.tolist(),
                "active_atoms": sorted(set((active_dofs // 3).tolist())),
            }
            saved_active_dofs = active_dofs
        elif hess.shape != (full_n_dof, full_n_dof):
            raise ValueError(
                "Partial Hessian file has no active-DOF metadata; regenerate it with `mlmm freq --dump-hess`."
            )

        if expected_active_dofs is not None:
            expected = np.asarray(expected_active_dofs, dtype=np.int64).reshape(-1)
            if (
                np.unique(expected).size != expected.size
                or np.any(expected < 0)
                or np.any(expected >= full_n_dof)
            ):
                raise ValueError("Current calculator supplied an invalid active-DOF basis.")
            if not np.array_equal(saved_active_dofs, expected):
                raise ValueError(
                    "Hessian file active-DOF basis does not match the current "
                    "ML/MM layer selection; regenerate it with the same layer/Hessian settings."
                )

        return {
            "hessian": hess.copy(),
            "energy_ha": float(np.asarray(data["energy_ha"]).item()),
            "partial_metadata": partial,
        }
