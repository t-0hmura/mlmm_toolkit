from __future__ import annotations

from dataclasses import dataclass

import torch
from torch import nn

from ..constants import COULOMB_K, validate_coords
from ..native.nonbonded import (
    nonbonded_energy_force_native,
    nonbonded_energy_force_preparam_native,
)


@dataclass
class NonbondedEnergies:
    coulomb: torch.Tensor
    lj: torch.Tensor
    coulomb14: torch.Tensor
    lj14: torch.Tensor


class NonbondedTerm(nn.Module):
    """No-PBC AMBER nonbonded term (native C++ backend only)."""

    def __init__(
        self,
        natom: int,
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
        cpu_fast_kernel: bool = True,
    ):
        super().__init__()
        self.natom = int(natom)
        self.cpu_fast_kernel = bool(cpu_fast_kernel)

        self.register_buffer("charge", charge)
        self.register_buffer("atom_type", atom_type.long())
        self.register_buffer("lj_acoef", lj_acoef)
        self.register_buffer("lj_bcoef", lj_bcoef)
        self.register_buffer("hb_acoef", hb_acoef)
        self.register_buffer("hb_bcoef", hb_bcoef)
        self.register_buffer("nb_index", nb_index.long())

        self.register_buffer("pair_i", pair_i.long())
        self.register_buffer("pair_j", pair_j.long())
        self.register_buffer("pair14_i", pair14_i.long())
        self.register_buffer("pair14_j", pair14_j.long())
        self.register_buffer("pair14_inv_scee", pair14_inv_scee)
        self.register_buffer("pair14_inv_scnb", pair14_inv_scnb)

        # Precompute pair constants once so repeated coordinate updates on the
        # same topology skip atom-type/table lookups.
        p_coul, p_a12, p_b6, p_b10 = self._build_pair_params(
            self.pair_i,
            self.pair_j,
            inv_scee=None,
            inv_scnb=None,
        )
        p14_coul, p14_a12, p14_b6, p14_b10 = self._build_pair_params(
            self.pair14_i,
            self.pair14_j,
            inv_scee=self.pair14_inv_scee,
            inv_scnb=self.pair14_inv_scnb,
        )
        self.register_buffer("pair_coul_coeff", p_coul)
        self.register_buffer("pair_a12_coeff", p_a12)
        self.register_buffer("pair_b6_coeff", p_b6)
        self.register_buffer("pair_b10_coeff", p_b10)
        self.register_buffer("pair14_coul_coeff", p14_coul)
        self.register_buffer("pair14_a12_coeff", p14_a12)
        self.register_buffer("pair14_b6_coeff", p14_b6)
        self.register_buffer("pair14_b10_coeff", p14_b10)

    def _build_pair_params(
        self,
        pair_i: torch.Tensor,
        pair_j: torch.Tensor,
        inv_scee: torch.Tensor | None,
        inv_scnb: torch.Tensor | None,
    ) -> tuple[torch.Tensor, torch.Tensor, torch.Tensor, torch.Tensor]:
        dtype = self.charge.dtype
        device = self.charge.device
        n = int(pair_i.numel())
        if n == 0:
            z = torch.zeros((0,), dtype=dtype, device=device)
            return z, z, z, z

        ii = pair_i
        jj = pair_j
        qq = self.charge[ii] * self.charge[jj]
        coul = qq * float(COULOMB_K)
        if inv_scee is not None and inv_scee.numel() > 0:
            coul = coul * inv_scee.to(dtype=dtype)

        a12 = torch.zeros_like(coul)
        b6 = torch.zeros_like(coul)
        b10 = torch.zeros_like(coul)

        ti = self.atom_type[ii]
        tj = self.atom_type[jj]
        raw = self.nb_index[ti, tj]

        lj_mask = raw > 0
        if bool(torch.any(lj_mask)):
            lj_idx = raw[lj_mask] - 1
            a = self.lj_acoef.index_select(0, lj_idx)
            b = self.lj_bcoef.index_select(0, lj_idx)
            if inv_scnb is not None and inv_scnb.numel() > 0:
                s = inv_scnb[lj_mask].to(dtype=dtype)
                a = a * s
                b = b * s
            a12.index_put_((lj_mask,), a)
            b6.index_put_((lj_mask,), b)

        hb_mask = raw < 0
        if bool(torch.any(hb_mask)) and int(self.hb_acoef.numel()) > 0 and int(self.hb_bcoef.numel()) > 0:
            hb_idx_full = (-raw[hb_mask]) - 1
            valid = (
                (hb_idx_full >= 0)
                & (hb_idx_full < int(self.hb_acoef.numel()))
                & (hb_idx_full < int(self.hb_bcoef.numel()))
            )
            if bool(torch.any(valid)):
                hb_pos = torch.nonzero(hb_mask, as_tuple=False).reshape(-1)[valid]
                hb_idx = hb_idx_full[valid]
                a = self.hb_acoef.index_select(0, hb_idx)
                b = self.hb_bcoef.index_select(0, hb_idx)
                if inv_scnb is not None and inv_scnb.numel() > 0:
                    s = inv_scnb[hb_pos].to(dtype=dtype)
                    a = a * s
                    b = b * s
                a12.index_put_((hb_pos,), a)
                b10.index_put_((hb_pos,), b)

        return coul, a12, b6, b10

    def _native_energy_force(
        self,
        coords: torch.Tensor,
    ) -> tuple[NonbondedEnergies, torch.Tensor]:
        validate_coords(
            coords, self.natom, self.charge.dtype, self.charge.device,
            label="coords (nonbonded)",
        )
        if not self.cpu_fast_kernel:
            out = nonbonded_energy_force_preparam_native(
                coords=coords,
                pair_i=self.pair_i,
                pair_j=self.pair_j,
                pair_coul_coeff=self.pair_coul_coeff,
                pair_a12_coeff=self.pair_a12_coeff,
                pair_b6_coeff=self.pair_b6_coeff,
                pair_b10_coeff=self.pair_b10_coeff,
                pair14_i=self.pair14_i,
                pair14_j=self.pair14_j,
                pair14_coul_coeff=self.pair14_coul_coeff,
                pair14_a12_coeff=self.pair14_a12_coeff,
                pair14_b6_coeff=self.pair14_b6_coeff,
                pair14_b10_coeff=self.pair14_b10_coeff,
                cpu_fast=False,
            )
        else:
            out = nonbonded_energy_force_native(
                coords=coords,
                charge=self.charge,
                atom_type=self.atom_type,
                lj_acoef=self.lj_acoef,
                lj_bcoef=self.lj_bcoef,
                hb_acoef=self.hb_acoef,
                hb_bcoef=self.hb_bcoef,
                nb_index=self.nb_index,
                pair_i=self.pair_i,
                pair_j=self.pair_j,
                pair14_i=self.pair14_i,
                pair14_j=self.pair14_j,
                pair14_inv_scee=self.pair14_inv_scee,
                pair14_inv_scnb=self.pair14_inv_scnb,
                cpu_fast=self.cpu_fast_kernel,
            )

        e_coul, e_lj, e_coul14, e_lj14, force = out
        return (
            NonbondedEnergies(
                coulomb=e_coul,
                lj=e_lj,
                coulomb14=e_coul14,
                lj14=e_lj14,
            ),
            force,
        )

    def forward(self, coords: torch.Tensor) -> NonbondedEnergies:
        return self._native_energy_force(coords)[0]

    def energy_force(self, coords: torch.Tensor) -> tuple[NonbondedEnergies, torch.Tensor]:
        """Return nonbonded energies and analytical force (native backend)."""
        return self._native_energy_force(coords)
