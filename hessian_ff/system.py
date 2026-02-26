from __future__ import annotations

import dataclasses
from dataclasses import dataclass
from typing import Optional, Tuple

import torch


@dataclass(frozen=True)
class AmberSystem:
    """Tensorized representation of an AMBER prmtop system.

    Floating tensors are loaded as torch.float64 by default and can be cast to
    torch.float32 when needed.
    Coordinates are expected in Angstrom.

    Notes on charge units:
      - This package reads prmtop via ParmEd, which returns atomic charges in units
        of elementary charge (i.e., already de-scaled from raw %FLAG CHARGE).
      - Coulomb energy in kcal/mol is computed as:
            E = 332.0637132991921 * sum_{i<j} q_i q_j / r_ij
        where r_ij is in Angstrom.
    """

    # ---- atom-level ----
    natom: int
    charge: torch.Tensor  # [N] atomic charges in e
    atom_type: torch.Tensor  # [N] Lennard-Jones type index (0-based int64)

    # ---- bonded terms ----
    bond_i: torch.Tensor  # [Nb]
    bond_j: torch.Tensor  # [Nb]
    bond_k: torch.Tensor  # [Nb] force constant
    bond_r0: torch.Tensor  # [Nb] equilibrium distance

    angle_i: torch.Tensor  # [Na]
    angle_j: torch.Tensor  # [Na]
    angle_k: torch.Tensor  # [Na]
    angle_k0: torch.Tensor  # [Na] force constant
    angle_t0: torch.Tensor  # [Na] equilibrium angle (radians)

    dihed_i: torch.Tensor  # [Nd]
    dihed_j: torch.Tensor  # [Nd]
    dihed_k: torch.Tensor  # [Nd]
    dihed_l: torch.Tensor  # [Nd]
    dihed_force: torch.Tensor  # [Nd]
    dihed_period: torch.Tensor  # [Nd]
    dihed_phase: torch.Tensor  # [Nd] (radians)

    # ---- nonbonded parameter tables ----
    lj_acoef: torch.Tensor  # [Nlj] A in A/r^12
    lj_bcoef: torch.Tensor  # [Nlj] B in B/r^6
    nb_index: torch.Tensor  # [ntypes, ntypes] raw NONBONDED_PARM_INDEX (Fortran-style signed int)
    hb_acoef: torch.Tensor  # [Nhb] HBOND A in A/r^12 (used when nb_index<0)
    hb_bcoef: torch.Tensor  # [Nhb] HBOND B in B/r^10 (used when nb_index<0)

    # ---- pair lists (precomputed for no-PBC O(N^2) evaluation) ----
    pair_i: torch.Tensor  # [Np] nonbonded general pairs (i<j), excludes exclusions and 1-4
    pair_j: torch.Tensor  # [Np]

    pair14_i: torch.Tensor  # [N14] 1-4 pairs (i<j)
    pair14_j: torch.Tensor  # [N14]
    pair14_inv_scee: torch.Tensor  # [N14] multiply Coulomb by this
    pair14_inv_scnb: torch.Tensor  # [N14] multiply LJ by this

    # ---- CMAP term (optional; CHARMM-style correction map) ----
    cmap_type: torch.Tensor  # [Ncmap] map type index (0-based)
    cmap_i: torch.Tensor  # [Ncmap] first torsion atom i
    cmap_j: torch.Tensor  # [Ncmap] first torsion atom j
    cmap_k: torch.Tensor  # [Ncmap] first torsion atom k
    cmap_l: torch.Tensor  # [Ncmap] first torsion atom l
    cmap_m: torch.Tensor  # [Ncmap] second torsion terminal atom m
    cmap_resolution: Tuple[int, ...]  # one resolution per map type
    cmap_maps: Tuple[torch.Tensor, ...]  # flattened map values (kcal/mol) in OpenMM ordering
    cmap_size: torch.Tensor  # [Nmap] grid size per map type
    cmap_delta: torch.Tensor  # [Nmap] angular grid width (radian)
    cmap_offset: torch.Tensor  # [Nmap] flattened patch row offset
    cmap_coeff: torch.Tensor  # [sum(size*size),16] bicubic coefficients

    def to(
        self,
        device: torch.device | str | None = None,
        dtype: Optional[torch.dtype] = None,
    ) -> "AmberSystem":
        """Return a copy of this system moved to a device and/or floating dtype."""

        target_device = torch.device(device) if device is not None else None

        def mv(x: torch.Tensor) -> torch.Tensor:
            kw = {}
            if target_device is not None:
                kw["device"] = target_device
            if dtype is not None and x.is_floating_point():
                kw["dtype"] = dtype
            return x.to(**kw) if kw else x

        kwargs = {}
        for f in dataclasses.fields(self):
            val = getattr(self, f.name)
            if isinstance(val, torch.Tensor):
                kwargs[f.name] = mv(val)
            elif isinstance(val, tuple) and val and isinstance(val[0], torch.Tensor):
                kwargs[f.name] = tuple(mv(x) for x in val)
            else:
                kwargs[f.name] = val
        return AmberSystem(**kwargs)
