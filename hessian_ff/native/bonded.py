from __future__ import annotations

import torch

from .loader import bonded_extension_status, get_bonded_extension


def bonded_energy_force_native(
    *,
    coords: torch.Tensor,
    bond_i: torch.Tensor,
    bond_j: torch.Tensor,
    bond_k: torch.Tensor,
    bond_r0: torch.Tensor,
    angle_i: torch.Tensor,
    angle_j: torch.Tensor,
    angle_k: torch.Tensor,
    angle_k0: torch.Tensor,
    angle_t0: torch.Tensor,
    dihed_i: torch.Tensor,
    dihed_j: torch.Tensor,
    dihed_k: torch.Tensor,
    dihed_l: torch.Tensor,
    dihed_force: torch.Tensor,
    dihed_period: torch.Tensor,
    dihed_phase: torch.Tensor,
    cmap_type: torch.Tensor,
    cmap_i: torch.Tensor,
    cmap_j: torch.Tensor,
    cmap_k: torch.Tensor,
    cmap_l: torch.Tensor,
    cmap_m: torch.Tensor,
    cmap_size: torch.Tensor,
    cmap_delta: torch.Tensor,
    cmap_offset: torch.Tensor,
    cmap_coeff: torch.Tensor,
    verbose: bool = False,
) -> tuple[torch.Tensor, torch.Tensor, torch.Tensor, torch.Tensor, torch.Tensor]:
    """Call C++ bonded extension.

    Returns
    -------
    (E_bond, E_angle, E_dihedral, E_cmap, force)
    """
    ext = get_bonded_extension(verbose=verbose, force_rebuild=False)
    if ext is None:
        status = bonded_extension_status()
        raise RuntimeError(
            "native bonded extension is unavailable. "
            f"detail: {status.get('note')}"
        )

    out = ext.bonded_energy_force(
        coords,
        bond_i,
        bond_j,
        bond_k,
        bond_r0,
        angle_i,
        angle_j,
        angle_k,
        angle_k0,
        angle_t0,
        dihed_i,
        dihed_j,
        dihed_k,
        dihed_l,
        dihed_force,
        dihed_period,
        dihed_phase,
        cmap_type,
        cmap_i,
        cmap_j,
        cmap_k,
        cmap_l,
        cmap_m,
        cmap_size,
        cmap_delta,
        cmap_offset,
        cmap_coeff,
    )
    return out[0], out[1], out[2], out[3], out[4]
