"""Geometry-identified Hessian files for ``freq`` -> ``irc`` handoff."""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any, Mapping, Optional

import numpy as np


SCHEMA_VERSION = 3
LEGACY_IDENTIFIED_SCHEMA_VERSION = 1
ELECTRONIC_STATE_SCHEMA_VERSION = 2


def _canonical_identity_json(identity: Mapping[str, Any]) -> str:
    """Serialize a persistent PES identity without pickle/object arrays."""

    if not isinstance(identity, Mapping):
        raise ValueError("potential_identity must be a mapping.")
    try:
        return json.dumps(
            dict(identity),
            sort_keys=True,
            separators=(",", ":"),
            allow_nan=False,
        )
    except (TypeError, ValueError) as exc:
        raise ValueError(
            "potential_identity must contain finite JSON-compatible values."
        ) from exc


def _validated_state_integer(value: Any, *, name: str, minimum: Optional[int] = None) -> int:
    """Return a scalar integer suitable for persistent electronic-state metadata."""
    array = np.asarray(value)
    if array.ndim != 0 or array.dtype.kind not in {"i", "u"}:
        raise ValueError(f"{name} must be a scalar integer.")
    result = int(array.item())
    if minimum is not None and result < minimum:
        raise ValueError(f"{name} must be >= {minimum}.")
    return result


def _validated_finite_scalar(value: Any, *, name: str) -> float:
    """Return a finite scalar floating-point value."""
    array = np.asarray(value)
    if array.ndim != 0:
        raise ValueError(f"{name} must be a finite scalar.")
    try:
        result = float(array.item())
    except (TypeError, ValueError) as exc:
        raise ValueError(f"{name} must be a finite scalar.") from exc
    if not np.isfinite(result):
        raise ValueError(f"{name} must be a finite scalar.")
    return result


def save_hessian_file(
    path: Path | str,
    *,
    hessian: Any,
    energy_ha: float,
    cart_coords_bohr: Any,
    atomic_numbers: Any,
    model_charge: int,
    model_mult: int,
    potential_identity: Mapping[str, Any],
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
    saved_energy = _validated_finite_scalar(energy_ha, name="energy_ha")
    saved_charge = _validated_state_integer(model_charge, name="model_charge")
    saved_mult = _validated_state_integer(model_mult, name="model_mult", minimum=1)

    payload: dict[str, Any] = {
        "schema_version": np.int64(SCHEMA_VERSION),
        "hessian": hess,
        "energy_ha": np.float64(saved_energy),
        "cart_coords_bohr": coords,
        "atomic_numbers": numbers,
        "model_charge": np.int64(saved_charge),
        "model_mult": np.int64(saved_mult),
        "potential_identity_json": np.str_(
            _canonical_identity_json(potential_identity)
        ),
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
    expected_model_charge: int,
    expected_model_mult: int,
    expected_potential_identity: Optional[Mapping[str, Any]] = None,
    expected_active_dofs: Any = None,
    allow_unverified_state: bool = False,
    allow_unverified_pes: bool = False,
    atol_bohr: float = 1.1e-3,
) -> dict[str, Any]:
    """Load and validate a Hessian against the current geometry and DOF basis."""
    current_coords = np.asarray(cart_coords_bohr, dtype=float).reshape(-1)
    current_numbers = np.asarray(atomic_numbers, dtype=np.int64).reshape(-1)
    current_charge = _validated_state_integer(
        expected_model_charge, name="expected_model_charge"
    )
    current_mult = _validated_state_integer(
        expected_model_mult, name="expected_model_mult", minimum=1
    )
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
        version = _validated_state_integer(
            data["schema_version"], name="schema_version", minimum=1
        )
        if version not in {
            LEGACY_IDENTIFIED_SCHEMA_VERSION,
            ELECTRONIC_STATE_SCHEMA_VERSION,
            SCHEMA_VERSION,
        }:
            raise ValueError(
                f"Unsupported Hessian file schema {version}; expected "
                f"{LEGACY_IDENTIFIED_SCHEMA_VERSION}, "
                f"{ELECTRONIC_STATE_SCHEMA_VERSION}, or {SCHEMA_VERSION}."
            )
        saved_charge: Optional[int]
        saved_mult: Optional[int]
        electronic_state_verified: bool
        state_fields = {"model_charge", "model_mult"}
        if version >= ELECTRONIC_STATE_SCHEMA_VERSION:
            missing_state = sorted(state_fields.difference(data.files))
            if missing_state:
                raise ValueError(
                    "Hessian file schema 2 lacks electronic-state metadata "
                    f"({', '.join(missing_state)}); regenerate it with "
                    "`mlmm freq --dump-hess`."
                )
            saved_charge = _validated_state_integer(
                data["model_charge"], name="saved model_charge"
            )
            saved_mult = _validated_state_integer(
                data["model_mult"], name="saved model_mult", minimum=1
            )
            if saved_charge != current_charge or saved_mult != current_mult:
                raise ValueError(
                    "Hessian file electronic state does not match the current "
                    "calculation: saved charge/multiplicity "
                    f"{saved_charge}/{saved_mult}, current "
                    f"{current_charge}/{current_mult}."
                )
            electronic_state_verified = True
        else:
            saved_charge = None
            saved_mult = None
            if not allow_unverified_state:
                raise ValueError(
                    "Hessian file schema 1 does not identify model charge and "
                    "multiplicity. Regenerate it with `mlmm freq --dump-hess`, "
                    "or independently verify the state and pass "
                    "`--allow-unverified-hess-state`."
                )
            electronic_state_verified = False

        saved_potential_identity: Optional[dict[str, Any]] = None
        potential_identity_verified = False
        if version == SCHEMA_VERSION:
            if "potential_identity_json" not in data.files:
                raise ValueError(
                    "Hessian file schema 3 lacks PES identity metadata; "
                    "regenerate it with `mlmm freq --dump-hess`."
                )
            try:
                saved_potential_identity = json.loads(
                    str(np.asarray(data["potential_identity_json"]).item())
                )
            except (TypeError, ValueError, json.JSONDecodeError) as exc:
                raise ValueError(
                    "Hessian file contains invalid PES identity metadata."
                ) from exc
            if not isinstance(saved_potential_identity, dict):
                raise ValueError(
                    "Hessian file PES identity metadata must be a mapping."
                )
            if expected_potential_identity is None:
                raise ValueError(
                    "Current calculation did not supply a PES identity for "
                    "Hessian verification."
                )
            expected_identity = json.loads(
                _canonical_identity_json(expected_potential_identity)
            )
            if saved_potential_identity != expected_identity:
                raise ValueError(
                    "Hessian file PES identity does not match the current "
                    "backend/model/precision/topology/layer/link/Hessian settings."
                )
            potential_identity_verified = True
        elif not allow_unverified_pes:
            raise ValueError(
                f"Hessian file schema {version} does not identify the generating "
                "PES. Regenerate it with `mlmm freq --dump-hess`, or independently "
                "verify the evaluator and pass `--allow-unverified-hess-state`."
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

        energy_ha = _validated_finite_scalar(data["energy_ha"], name="energy_ha")
        return {
            "hessian": hess.copy(),
            "energy_ha": energy_ha,
            "partial_metadata": partial,
            "schema_version": version,
            "model_charge": saved_charge,
            "model_mult": saved_mult,
            "electronic_state_verified": electronic_state_verified,
            "potential_identity": saved_potential_identity,
            "potential_identity_verified": potential_identity_verified,
        }
