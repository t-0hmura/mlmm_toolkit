from __future__ import annotations

import torch

from .loader import get_nonbonded_extension, nonbonded_extension_status


def nonbonded_energy_force_native(
    *,
    coords: torch.Tensor,
    charge: torch.Tensor,
    atom_type: torch.Tensor,
    lj_acoef: torch.Tensor,
    lj_bcoef: torch.Tensor,
    hb_acoef: torch.Tensor,
    hb_bcoef: torch.Tensor,
    nb_index: torch.Tensor,
    pair_i: torch.Tensor,
    pair_j: torch.Tensor,
    pair14_i: torch.Tensor,
    pair14_j: torch.Tensor,
    pair14_inv_scee: torch.Tensor,
    pair14_inv_scnb: torch.Tensor,
    cpu_fast: bool = True,
    verbose: bool = False,
) -> tuple[torch.Tensor, torch.Tensor, torch.Tensor, torch.Tensor, torch.Tensor]:
    """Call C++ nonbonded extension.

    Returns
    -------
    (E_coul, E_lj, E_coul14, E_lj14, force)

    References:
    - native/nonbonded_ext.cpp in this repository.
    """
    ext = get_nonbonded_extension(verbose=verbose, force_rebuild=False)
    if ext is None:
        status = nonbonded_extension_status()
        raise RuntimeError(
            "native nonbonded extension is unavailable; torch fallback is disabled. "
            f"detail: {status.get('note')}"
        )

    chunk_size = max(int(pair_i.numel()), int(pair14_i.numel()), 1)
    out = ext.nonbonded_energy_force(
        coords,
        charge,
        atom_type,
        lj_acoef,
        lj_bcoef,
        hb_acoef,
        hb_bcoef,
        nb_index,
        pair_i,
        pair_j,
        pair14_i,
        pair14_j,
        pair14_inv_scee,
        pair14_inv_scnb,
        int(chunk_size),
        bool(cpu_fast),
        -1,
    )
    return out[0], out[1], out[2], out[3], out[4]


def nonbonded_energy_force_preparam_native(
    *,
    coords: torch.Tensor,
    pair_i: torch.Tensor,
    pair_j: torch.Tensor,
    pair_coul_coeff: torch.Tensor,
    pair_a12_coeff: torch.Tensor,
    pair_b6_coeff: torch.Tensor,
    pair_b10_coeff: torch.Tensor,
    pair14_i: torch.Tensor,
    pair14_j: torch.Tensor,
    pair14_coul_coeff: torch.Tensor,
    pair14_a12_coeff: torch.Tensor,
    pair14_b6_coeff: torch.Tensor,
    pair14_b10_coeff: torch.Tensor,
    cpu_fast: bool = True,
    verbose: bool = False,
) -> tuple[torch.Tensor, torch.Tensor, torch.Tensor, torch.Tensor, torch.Tensor]:
    """Call C++ nonbonded extension using precomputed pair coefficients.

    Returns
    -------
    (E_coul, E_lj, E_coul14, E_lj14, force)
    """
    ext = get_nonbonded_extension(verbose=verbose, force_rebuild=False)
    if ext is None:
        status = nonbonded_extension_status()
        raise RuntimeError(
            "native nonbonded extension is unavailable; torch fallback is disabled. "
            f"detail: {status.get('note')}"
        )

    chunk_size = max(int(pair_i.numel()), int(pair14_i.numel()), 1)
    out = ext.nonbonded_energy_force_preparam(
        coords,
        pair_i,
        pair_j,
        pair_coul_coeff,
        pair_a12_coeff,
        pair_b6_coeff,
        pair_b10_coeff,
        pair14_i,
        pair14_j,
        pair14_coul_coeff,
        pair14_a12_coeff,
        pair14_b6_coeff,
        pair14_b10_coeff,
        int(chunk_size),
        bool(cpu_fast),
        -1,
    )
    return out[0], out[1], out[2], out[3], out[4]
