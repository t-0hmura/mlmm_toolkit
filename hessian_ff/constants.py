"""Shared constants and validation utilities for hessian_ff."""

from __future__ import annotations

import math

import torch

# Coulomb constant for kcal/mol when charges are in elementary charge units
# and distance in Angstrom.
COULOMB_K = 332.0637132991921

TWO_PI = 2.0 * math.pi


def validate_coords(
    coords: torch.Tensor,
    natom: int,
    expected_dtype: torch.dtype,
    expected_device: torch.device,
    *,
    label: str = "coords",
) -> None:
    """Validate coordinate tensor shape, dtype, and device."""
    if not torch.is_tensor(coords):
        raise TypeError(f"{label} must be torch.Tensor, got {type(coords)!r}")
    if coords.ndim != 2:
        raise ValueError(f"{label} must have shape [N,3], got ndim={coords.ndim}")
    if int(coords.shape[0]) != int(natom) or int(coords.shape[1]) != 3:
        raise ValueError(
            f"{label} shape mismatch: expected [{natom},3], got {tuple(coords.shape)}"
        )
    if not coords.is_floating_point():
        raise TypeError(f"{label} must be floating tensor, got dtype={coords.dtype}")
    if coords.dtype != expected_dtype:
        raise ValueError(
            f"{label} dtype mismatch: "
            f"got {coords.dtype}, expected {expected_dtype}. "
            "Use matching precision for system and coordinates."
        )
    if coords.device != expected_device:
        raise ValueError(
            f"{label} device mismatch: "
            f"got {coords.device}, expected {expected_device}. "
            "Move coordinates to the same device as the force-field parameters."
        )
