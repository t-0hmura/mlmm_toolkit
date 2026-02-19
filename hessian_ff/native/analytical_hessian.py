from __future__ import annotations

import torch

from .loader import (
    analytical_hessian_extension_status,
    get_analytical_hessian_extension,
)


def _require_ext(verbose: bool = False):
    ext = get_analytical_hessian_extension(verbose=verbose, force_rebuild=False)
    if ext is None:
        st = analytical_hessian_extension_status()
        raise RuntimeError(
            "analytical Hessian native extension is unavailable. "
            f"detail: {st.get('note')}"
        )
    return ext


def bond_hessian_native(
    *,
    coords: torch.Tensor,
    bond_i: torch.Tensor,
    bond_j: torch.Tensor,
    bond_k: torch.Tensor,
    bond_r0: torch.Tensor,
    active_map: torch.Tensor,
    ndof: int,
) -> torch.Tensor:
    ext = _require_ext(verbose=False)
    return ext.bond_hessian(
        coords,
        bond_i,
        bond_j,
        bond_k,
        bond_r0,
        active_map,
        int(ndof),
    )


def nonbonded_hessian_native(
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
    inv_scee: torch.Tensor | None,
    inv_scnb: torch.Tensor | None,
    active_map: torch.Tensor,
    ndof: int,
) -> torch.Tensor:
    ext = _require_ext(verbose=False)
    # pybind optional Tensor args accept None.
    return ext.nonbonded_hessian(
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
        inv_scee,
        inv_scnb,
        active_map,
        int(ndof),
    )


def scatter_local_hessian_native(
    *,
    local_h: torch.Tensor,
    dof_idx: torch.Tensor,
    ndof: int,
) -> torch.Tensor:
    ext = _require_ext(verbose=False)
    return ext.scatter_local_hessian(local_h, dof_idx, int(ndof))
