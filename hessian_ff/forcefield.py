from __future__ import annotations

from typing import Dict, Optional, Union

import torch
from torch import nn

from .constants import validate_coords
from .system import AmberSystem
from .native.bonded import bonded_energy_force_native
from .terms.angle import AngleTerm
from .terms.bond import BondTerm
from .terms.cmap import CMapTerm
from .terms.dihedral import DihedralTerm
from .terms.nonbonded import NonbondedTerm


class ForceFieldTorch(nn.Module):
    """AMBER force field evaluator in PyTorch (no PBC, no Ewald)."""

    ENERGY_KEYS = (
        "E_total",
        "E_bond",
        "E_angle",
        "E_dihedral",
        "E_cmap",
        "E_coul",
        "E_lj",
        "E_coul14",
        "E_lj14",
        "E_nonbonded_total",
    )
    def __init__(
        self,
        system: AmberSystem,
        nonbonded_cpu_fast: bool = True,
    ):
        super().__init__()
        self.system = system

        self.bond = BondTerm(system.bond_i, system.bond_j, system.bond_k, system.bond_r0)
        self.angle = AngleTerm(
            system.angle_i, system.angle_j, system.angle_k, system.angle_k0, system.angle_t0
        )
        self.dihedral = DihedralTerm(
            system.dihed_i,
            system.dihed_j,
            system.dihed_k,
            system.dihed_l,
            system.dihed_force,
            system.dihed_period,
            system.dihed_phase,
        )
        self.cmap = CMapTerm(
            natom=system.natom,
            cmap_type=system.cmap_type,
            cmap_i=system.cmap_i,
            cmap_j=system.cmap_j,
            cmap_k=system.cmap_k,
            cmap_l=system.cmap_l,
            cmap_m=system.cmap_m,
            cmap_resolution=system.cmap_resolution,
            cmap_maps=system.cmap_maps,
            precomputed_size=system.cmap_size,
            precomputed_delta=system.cmap_delta,
            precomputed_offset=system.cmap_offset,
            precomputed_coeff=system.cmap_coeff,
        )
        self.nonbonded = NonbondedTerm(
            natom=system.natom,
            charge=system.charge,
            atom_type=system.atom_type,
            lj_acoef=system.lj_acoef,
            lj_bcoef=system.lj_bcoef,
            hb_acoef=system.hb_acoef,
            hb_bcoef=system.hb_bcoef,
            nb_index=system.nb_index,
            pair_i=system.pair_i,
            pair_j=system.pair_j,
            pair14_i=system.pair14_i,
            pair14_j=system.pair14_j,
            pair14_inv_scee=system.pair14_inv_scee,
            pair14_inv_scnb=system.pair14_inv_scnb,
            cpu_fast_kernel=nonbonded_cpu_fast,
        )

    def _validate_coords(self, coords: torch.Tensor) -> None:
        validate_coords(
            coords,
            self.system.natom,
            self.system.charge.dtype,
            self.system.charge.device,
        )

    def _bonded_energy_force_torch(
        self, coords: torch.Tensor
    ) -> tuple[torch.Tensor, torch.Tensor, torch.Tensor, torch.Tensor, torch.Tensor]:
        e_bond, f_bond = self.bond.energy_force(coords)
        e_angle, f_angle = self.angle.energy_force(coords)
        e_dihed, f_dihed = self.dihedral.energy_force(coords)
        e_cmap, f_cmap = self.cmap.energy_force(coords)
        force = f_bond + f_angle + f_dihed + f_cmap
        return e_bond, e_angle, e_dihed, e_cmap, force

    def _native_bonded_energy_force(
        self, coords: torch.Tensor
    ) -> tuple[torch.Tensor, torch.Tensor, torch.Tensor, torch.Tensor, torch.Tensor]:
        s = self.system
        return bonded_energy_force_native(
            coords=coords,
            bond_i=s.bond_i,
            bond_j=s.bond_j,
            bond_k=s.bond_k,
            bond_r0=s.bond_r0,
            angle_i=s.angle_i,
            angle_j=s.angle_j,
            angle_k=s.angle_k,
            angle_k0=s.angle_k0,
            angle_t0=s.angle_t0,
            dihed_i=s.dihed_i,
            dihed_j=s.dihed_j,
            dihed_k=s.dihed_k,
            dihed_l=s.dihed_l,
            dihed_force=s.dihed_force,
            dihed_period=s.dihed_period,
            dihed_phase=s.dihed_phase,
            cmap_type=s.cmap_type,
            cmap_i=s.cmap_i,
            cmap_j=s.cmap_j,
            cmap_k=s.cmap_k,
            cmap_l=s.cmap_l,
            cmap_m=s.cmap_m,
            cmap_size=s.cmap_size,
            cmap_delta=s.cmap_delta,
            cmap_offset=s.cmap_offset,
            cmap_coeff=s.cmap_coeff,
        )

    def _bonded_energy_force(
        self, coords: torch.Tensor
    ) -> tuple[torch.Tensor, torch.Tensor, torch.Tensor, torch.Tensor, torch.Tensor]:
        """Compute bonded E/F. Uses native on CPU+no-grad, otherwise torch path."""
        if (not bool(coords.requires_grad)) and coords.device.type == "cpu":
            return self._native_bonded_energy_force(coords)
        return self._bonded_energy_force_torch(coords)

    @staticmethod
    def _make_energy_dict(
        e_bond: torch.Tensor,
        e_angle: torch.Tensor,
        e_dihed: torch.Tensor,
        e_cmap: torch.Tensor,
        nb,
    ) -> Dict[str, torch.Tensor]:
        """Build the standard energy dict from bonded energies and nonbonded result."""
        e_nb_total = nb.coulomb + nb.lj + nb.coulomb14 + nb.lj14
        return {
            "E_total": e_bond + e_angle + e_dihed + e_cmap + e_nb_total,
            "E_bond": e_bond,
            "E_angle": e_angle,
            "E_dihedral": e_dihed,
            "E_cmap": e_cmap,
            "E_coul": nb.coulomb,
            "E_lj": nb.lj,
            "E_coul14": nb.coulomb14,
            "E_lj14": nb.lj14,
            "E_nonbonded_total": e_nb_total,
        }

    def forward(self, coords: torch.Tensor) -> Dict[str, torch.Tensor]:
        """Compute energies.

        Parameters
        ----------
        coords:
            Tensor of shape [N,3] in Angstrom. Should be float64 for best agreement.
        """
        self._validate_coords(coords)
        e_bond, e_angle, e_dihed, e_cmap, _ = self._bonded_energy_force(coords)
        nb = self.nonbonded(coords)
        return self._make_energy_dict(e_bond, e_angle, e_dihed, e_cmap, nb)

    def _validate_batch_coords(self, coords_batch: torch.Tensor) -> None:
        if coords_batch.ndim != 3:
            raise ValueError(
                f"coords_batch must have shape [B,N,3], got ndim={coords_batch.ndim}"
            )
        if int(coords_batch.shape[1]) != int(self.system.natom) or int(coords_batch.shape[2]) != 3:
            raise ValueError(
                f"coords_batch shape mismatch: expected [B,{self.system.natom},3], got {tuple(coords_batch.shape)}"
            )
        if not coords_batch.is_floating_point():
            raise TypeError(f"coords_batch must be floating tensor, got dtype={coords_batch.dtype}")
        if coords_batch.dtype != self.system.charge.dtype:
            raise ValueError(
                "coords_batch dtype mismatch: "
                f"got {coords_batch.dtype}, expected {self.system.charge.dtype}."
            )
        if coords_batch.device != self.system.charge.device:
            raise ValueError(
                "coords_batch device mismatch: "
                f"got {coords_batch.device}, expected {self.system.charge.device}."
            )

    def _forward_packed_single(self, coords: torch.Tensor) -> tuple[torch.Tensor, ...]:
        out = self.forward(coords)
        return tuple(out[k] for k in self.ENERGY_KEYS)

    def _energy_force_packed_single(
        self, coords: torch.Tensor, force_calc_mode: str
    ) -> tuple[torch.Tensor, ...]:
        out, force = self.energy_force(coords, force_calc_mode=force_calc_mode)
        return tuple(out[k] for k in self.ENERGY_KEYS) + (force,)

    def _microbatch_slices(self, total_batch: int, microbatch_size: Optional[int]) -> list[slice]:
        if total_batch == 0:
            return []
        if microbatch_size is None:
            return [slice(0, total_batch)]
        m = int(microbatch_size)
        if m < 1:
            raise ValueError(f"microbatch_size must be >=1, got {microbatch_size}")
        return [slice(i, min(i + m, total_batch)) for i in range(0, total_batch, m)]

    def forward_batch(
        self,
        coords_batch: torch.Tensor,
        batch_mode: str = "vmap",
        microbatch_size: Optional[int] = None,
    ) -> Dict[str, torch.Tensor]:
        """Compute energy terms for batched coordinates [B,N,3]."""
        self._validate_batch_coords(coords_batch)
        bsz = int(coords_batch.shape[0])
        if bsz == 0:
            z = coords_batch.new_zeros((0,))
            return {k: z for k in self.ENERGY_KEYS}

        mode = str(batch_mode).strip().lower()
        use_vmap = mode in {"vmap", "vectorized"}
        packed_chunks: list[tuple[torch.Tensor, ...]] = []
        for sl in self._microbatch_slices(bsz, microbatch_size):
            xb = coords_batch[sl]
            if use_vmap:
                packed = torch.vmap(self._forward_packed_single)(xb)
            else:
                outs = [self.forward(xb[i]) for i in range(int(xb.shape[0]))]
                packed = tuple(torch.stack([o[k] for o in outs], dim=0) for k in self.ENERGY_KEYS)
            packed_chunks.append(packed)

        merged = tuple(torch.cat([pc[i] for pc in packed_chunks], dim=0) for i in range(len(self.ENERGY_KEYS)))
        return {k: merged[i] for i, k in enumerate(self.ENERGY_KEYS)}

    def energy_force(
        self,
        coords: torch.Tensor,
        force_calc_mode: str = "Analytical",
    ) -> tuple[Dict[str, torch.Tensor], torch.Tensor]:
        """Compute energies and force.

        Parameters
        ----------
        force_calc_mode:
            - ``\"Analytical\"``: use analytical force expressions per term.
            - ``\"Autograd\"``: use Torch autograd on ``E_total``.
        """
        self._validate_coords(coords)

        mode = str(force_calc_mode).strip().lower()
        if mode == "autograd":
            x = coords.clone().detach().requires_grad_(True)
            out = self.forward(x)
            out["E_total"].backward()
            force = -x.grad.detach()
            out_detached = {k: v.detach() for k, v in out.items()}
            return out_detached, force

        if mode != "analytical":
            raise ValueError(f"Unknown force_calc_mode: {force_calc_mode!r}")

        e_bond, e_angle, e_dihed, e_cmap, force = self._bonded_energy_force(coords)
        nb_e, f_nb = self.nonbonded.energy_force(coords)
        force = force + f_nb
        out = self._make_energy_dict(e_bond, e_angle, e_dihed, e_cmap, nb_e)
        return out, force

    def energy_force_batch(
        self,
        coords_batch: torch.Tensor,
        force_calc_mode: str = "Analytical",
        batch_mode: str = "vmap",
        microbatch_size: Optional[int] = None,
    ) -> tuple[Dict[str, torch.Tensor], torch.Tensor]:
        """Compute energy and force for batched coordinates [B,N,3]."""
        self._validate_batch_coords(coords_batch)
        bsz = int(coords_batch.shape[0])
        if bsz == 0:
            z = coords_batch.new_zeros((0,))
            return {k: z for k in self.ENERGY_KEYS}, coords_batch.new_zeros(coords_batch.shape)

        mode = str(batch_mode).strip().lower()
        force_mode = str(force_calc_mode).strip()
        use_vmap = mode in {"vmap", "vectorized"} and force_mode.lower() == "analytical"

        e_chunks: list[Dict[str, torch.Tensor]] = []
        f_chunks: list[torch.Tensor] = []
        for sl in self._microbatch_slices(bsz, microbatch_size):
            xb = coords_batch[sl]
            if use_vmap:
                packed = torch.vmap(lambda x: self._energy_force_packed_single(x, force_mode))(xb)
                e_chunks.append({k: packed[i] for i, k in enumerate(self.ENERGY_KEYS)})
                f_chunks.append(packed[-1])
            else:
                e_list: list[Dict[str, torch.Tensor]] = []
                f_list: list[torch.Tensor] = []
                for i in range(int(xb.shape[0])):
                    out_i, force_i = self.energy_force(xb[i], force_calc_mode=force_mode)
                    e_list.append(out_i)
                    f_list.append(force_i)
                e_chunks.append({k: torch.stack([d[k] for d in e_list], dim=0) for k in self.ENERGY_KEYS})
                f_chunks.append(torch.stack(f_list, dim=0))

        e_out = {k: torch.cat([d[k] for d in e_chunks], dim=0) for k in self.ENERGY_KEYS}
        f_out = torch.cat(f_chunks, dim=0)
        return e_out, f_out

    @classmethod
    def from_prmtop(
        cls,
        prmtop_path: Union[str, "os.PathLike[str]"],
        device: Union[str, torch.device] = "cpu",
        nonbonded_cpu_fast: bool = True,
    ) -> "ForceFieldTorch":
        from .loaders import load_system

        system = load_system(prmtop_path, device=device)
        return cls(
            system,
            nonbonded_cpu_fast=nonbonded_cpu_fast,
        )
