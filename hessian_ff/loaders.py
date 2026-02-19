from __future__ import annotations

import math
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional, Sequence, Tuple, Union

import torch

from .prmtop_parmed import read_prmtop_with_parmed
from .system import AmberSystem
from .terms.cmap import _calc_map_coefficients


def _chunk(seq: Sequence, n: int) -> Iterable[Tuple]:
    for i in range(0, len(seq), n):
        yield tuple(seq[i : i + n])


def _atom_index_from_coord_index(coord_index: int) -> int:
    """prmtop stores atom references as (atom_index * 3) (0-based) with sign flags."""
    return abs(coord_index) // 3


def _infer_ntypes(nonbonded_parm_index: Sequence[int]) -> int:
    n = len(nonbonded_parm_index)
    rt = int(round(math.sqrt(n)))
    if rt * rt != n:
        raise ValueError(f"Cannot infer NTYPES: NONBONDED_PARM_INDEX length {n} is not a square")
    return rt


def _build_excluded_pair_keys(
    natom: int, num_excluded: Sequence[int], excluded_list: Sequence[int]
) -> torch.Tensor:
    """Reconstruct excluded pair keys from NUMBER_EXCLUDED_ATOMS + EXCLUDED_ATOMS_LIST.

    Returns a 1D tensor of keys = i*natom + j for i<j.
    """
    keys: List[int] = []
    cursor = 0
    for i in range(natom):
        n = int(num_excluded[i])
        partners = excluded_list[cursor : cursor + n]
        cursor += n
        for j1 in partners:
            j1 = int(j1)
            if j1 <= 0:
                continue
            j = j1 - 1  # prmtop stores 1-based atom numbers
            if j == i:
                continue
            a, b = (i, j) if i < j else (j, i)
            keys.append(a * natom + b)
    if not keys:
        return torch.zeros((0,), dtype=torch.int64)
    # unique
    return torch.unique(torch.tensor(keys, dtype=torch.int64))


def _build_14_pairs(
    natom: int,
    dihedrals_inc_h: Sequence[int],
    dihedrals_wo_h: Sequence[int],
    scee: Sequence[float],
    scnb: Sequence[float],
) -> Tuple[torch.Tensor, torch.Tensor, torch.Tensor, torch.Tensor]:
    """Extract unique 1-4 pairs and their scaling from dihedral lists.

    Per AMBER prmtop conventions:
      - If atom3 index is negative => do not compute 1-4 for that torsion
      - If atom4 index is negative => improper torsion

    We include 1-4 pairs only for entries with atom3>0 and atom4>0.
    """

    def add_from(dihed_list: Sequence[int], out: dict):
        for a1, a2, a3, a4, itype in _chunk(dihed_list, 5):
            a1 = int(a1)
            a2 = int(a2)
            a3 = int(a3)
            a4 = int(a4)
            it = int(itype) - 1  # Fortran 1-based
            if a3 < 0:
                continue  # explicitly no 1-4
            if a4 < 0:
                continue  # improper; do not generate 1-4
            i = _atom_index_from_coord_index(a1)
            j = _atom_index_from_coord_index(a4)
            if i == j:
                continue
            a, b = (i, j) if i < j else (j, i)
            key = a * natom + b
            inv_scee = 1.0 / float(scee[it]) if float(scee[it]) != 0 else 1.0
            inv_scnb = 1.0 / float(scnb[it]) if float(scnb[it]) != 0 else 1.0
            if key in out:
                # If multiple terms point to the same 1-4 pair, they should agree.
                prev = out[key]
                if (abs(prev[0] - inv_scee) > 1e-12) or (abs(prev[1] - inv_scnb) > 1e-12):
                    raise ValueError(
                        f"Inconsistent 1-4 scaling for pair ({a},{b}): {prev} vs {(inv_scee, inv_scnb)}"
                    )
            else:
                out[key] = (inv_scee, inv_scnb)

    mapping: dict[int, Tuple[float, float]] = {}
    add_from(dihedrals_inc_h, mapping)
    add_from(dihedrals_wo_h, mapping)

    if not mapping:
        z = torch.zeros((0,), dtype=torch.int64)
        f = torch.zeros((0,), dtype=torch.float64)
        return z, z, f, f

    keys = sorted(mapping.keys())
    i = torch.tensor([k // natom for k in keys], dtype=torch.int64)
    j = torch.tensor([k % natom for k in keys], dtype=torch.int64)
    inv_scee = torch.tensor([mapping[k][0] for k in keys], dtype=torch.float64)
    inv_scnb = torch.tensor([mapping[k][1] for k in keys], dtype=torch.float64)
    return i, j, inv_scee, inv_scnb


def _build_general_pairs(
    natom: int, excluded_keys: torch.Tensor, keys14: torch.Tensor
) -> Tuple[torch.Tensor, torch.Tensor]:
    """Build all i<j pairs excluding excluded_keys and 1-4 keys."""
    ti, tj = torch.triu_indices(natom, natom, offset=1)
    keys = ti * natom + tj
    mask = torch.ones_like(keys, dtype=torch.bool)
    if excluded_keys.numel() > 0:
        mask &= ~torch.isin(keys, excluded_keys)
    if keys14.numel() > 0:
        mask &= ~torch.isin(keys, keys14)
    return ti[mask], tj[mask]


def _parse_cmap(
    raw: Dict[str, List[Union[int, float, str]]]
) -> Tuple[
    torch.Tensor,
    torch.Tensor,
    torch.Tensor,
    torch.Tensor,
    torch.Tensor,
    torch.Tensor,
    Tuple[int, ...],
    Tuple[torch.Tensor, ...],
]:
    """Parse CMAP flags and return tensors in OpenMM-compatible ordering.

    Returns
    -------
    cmap_type, cmap_i, cmap_j, cmap_k, cmap_l, cmap_m, cmap_resolution, cmap_maps
    """

    count_flag = "CMAP_COUNT" if "CMAP_COUNT" in raw else "CHARMM_CMAP_COUNT"
    if count_flag not in raw:
        z = torch.zeros((0,), dtype=torch.int64)
        return z, z, z, z, z, z, tuple(), tuple()

    count_vals = raw.get(count_flag, [])
    # AMBER stores two ints; OpenMM reads the second.
    nmap = int(count_vals[1]) if len(count_vals) > 1 else int(count_vals[0])
    if nmap <= 0:
        z = torch.zeros((0,), dtype=torch.int64)
        return z, z, z, z, z, z, tuple(), tuple()

    res_flag = "CMAP_RESOLUTION" if "CMAP_RESOLUTION" in raw else "CHARMM_CMAP_RESOLUTION"
    if res_flag not in raw:
        raise ValueError("CMAP_COUNT is present but CMAP_RESOLUTION is missing")
    res_vals = raw.get(res_flag, [])
    if len(res_vals) < nmap:
        raise ValueError(
            f"CMAP_RESOLUTION too short: need {nmap}, got {len(res_vals)}"
        )
    cmap_resolution = tuple(int(x) for x in res_vals[:nmap])

    cmap_maps: List[torch.Tensor] = []
    for midx, ngrid in enumerate(cmap_resolution, start=1):
        key = f"CMAP_PARAMETER_{midx:02d}"
        if key not in raw:
            key = f"CHARMM_CMAP_PARAMETER_{midx:02d}"
        if key not in raw:
            raise ValueError(f"Missing {key} for CMAP map {midx}")
        vals = [float(x) for x in raw[key]]
        need = ngrid * ngrid
        if len(vals) < need:
            raise ValueError(
                f"{key} too short: need {need} values for {ngrid}x{ngrid}, got {len(vals)}"
            )
        vals = vals[:need]

        # Match OpenMM amber_file_parser readAmberSystem() remapping.
        remapped: List[float] = []
        for i in range(ngrid):
            for j in range(ngrid):
                src = ngrid * ((j + ngrid // 2) % ngrid) + ((i + ngrid // 2) % ngrid)
                remapped.append(vals[src])
        cmap_maps.append(torch.tensor(remapped, dtype=torch.float64))

    index_flag = "CMAP_INDEX" if "CMAP_INDEX" in raw else "CHARMM_CMAP_INDEX"
    if index_flag not in raw:
        raise ValueError("CMAP_COUNT is present but CMAP_INDEX is missing")
    cmap_idx = raw.get(index_flag, [])
    if len(cmap_idx) % 6 != 0:
        raise ValueError(f"{index_flag} length must be multiple of 6, got {len(cmap_idx)}")

    cmap_type: List[int] = []
    cmap_i: List[int] = []
    cmap_j: List[int] = []
    cmap_k: List[int] = []
    cmap_l: List[int] = []
    cmap_m: List[int] = []
    for a1, a2, a3, a4, a5, itype in _chunk(cmap_idx, 6):
        a1 = int(a1)
        a2 = int(a2)
        a3 = int(a3)
        a4 = int(a4)
        a5 = int(a5)
        t = int(itype) - 1
        if min(a1, a2, a3, a4, a5) <= 0:
            raise ValueError(
                f"Invalid CMAP atom pointer ({a1}, {a2}, {a3}, {a4}, {a5}); expected 1-based positive"
            )
        if t < 0 or t >= nmap:
            raise ValueError(f"Invalid CMAP map type index {int(itype)} (1-based), nmap={nmap}")
        cmap_type.append(t)
        # CMAP_INDEX stores 1-based atom indices (not coord-index*3)
        cmap_i.append(a1 - 1)
        cmap_j.append(a2 - 1)
        cmap_k.append(a3 - 1)
        cmap_l.append(a4 - 1)
        cmap_m.append(a5 - 1)

    return (
        torch.tensor(cmap_type, dtype=torch.int64),
        torch.tensor(cmap_i, dtype=torch.int64),
        torch.tensor(cmap_j, dtype=torch.int64),
        torch.tensor(cmap_k, dtype=torch.int64),
        torch.tensor(cmap_l, dtype=torch.int64),
        torch.tensor(cmap_m, dtype=torch.int64),
        cmap_resolution,
        tuple(cmap_maps),
    )


def _build_cmap_coeff_data(
    cmap_resolution: Tuple[int, ...],
    cmap_maps: Tuple[torch.Tensor, ...],
) -> tuple[torch.Tensor, torch.Tensor, torch.Tensor, torch.Tensor]:
    """Precompute CMAP bicubic coefficient tables once at load time."""
    if len(cmap_resolution) == 0:
        z_i64 = torch.zeros((0,), dtype=torch.int64)
        z_f64 = torch.zeros((0,), dtype=torch.float64)
        z_coeff = torch.zeros((0, 16), dtype=torch.float64)
        return z_i64, z_f64, z_i64, z_coeff

    map_size = torch.tensor([int(x) for x in cmap_resolution], dtype=torch.int64)
    coeff_parts: List[torch.Tensor] = []
    offsets: List[int] = []
    cursor = 0
    for ngrid, table in zip(cmap_resolution, cmap_maps):
        n = int(ngrid)
        offsets.append(cursor)
        c = _calc_map_coefficients(table, size=n)
        coeff_parts.append(c)
        cursor += n * n

    map_delta = (torch.full_like(map_size, 2.0 * math.pi, dtype=torch.float64) / map_size.to(torch.float64))
    map_offset = torch.tensor(offsets, dtype=torch.int64)
    map_coeff = torch.cat(coeff_parts, dim=0) if coeff_parts else torch.zeros((0, 16), dtype=torch.float64)
    return map_size, map_delta, map_offset, map_coeff


def load_system(
    prmtop_path: Union[str, Path],
    device: Union[str, torch.device] = "cpu",
) -> AmberSystem:
    """Load an AMBER parm7/prmtop into an :class:`~hessian_ff.system.AmberSystem`.

    Parameters
    ----------
    prmtop_path:
        Path to .prmtop/.parm7
    device:
        Target tensor device.
    """
    prmtop_path = Path(prmtop_path)

    raw = read_prmtop_with_parmed(prmtop_path)

    # ---- required sections ----
    charge = raw["CHARGE"]
    natom = len(charge)
    atom_type_index = raw["ATOM_TYPE_INDEX"]

    nbpi = raw["NONBONDED_PARM_INDEX"]
    ntypes = _infer_ntypes(nbpi)
    # Keep raw signed NONBONDED_PARM_INDEX values:
    #   >0 : index into Lennard-Jones A/B tables (Fortran 1-based)
    #   <0 : index into HBOND A/B tables (Fortran -1-based)
    nb_index = torch.tensor(nbpi, dtype=torch.int64).reshape(ntypes, ntypes)

    lj_acoef = torch.tensor(raw["LENNARD_JONES_ACOEF"], dtype=torch.float64)
    lj_bcoef = torch.tensor(raw["LENNARD_JONES_BCOEF"], dtype=torch.float64)
    hb_acoef = torch.tensor(raw.get("HBOND_ACOEF", []), dtype=torch.float64)
    hb_bcoef = torch.tensor(raw.get("HBOND_BCOEF", []), dtype=torch.float64)
    atom_type_t = torch.tensor(atom_type_index, dtype=torch.int64) - 1

    # ---- bonds ----
    bonds = list(_chunk(raw.get("BONDS_INC_HYDROGEN", []), 3)) + list(
        _chunk(raw.get("BONDS_WITHOUT_HYDROGEN", []), 3)
    )
    bond_k_list = raw.get("BOND_FORCE_CONSTANT", [])
    bond_r0_list = raw.get("BOND_EQUIL_VALUE", [])
    bond_i = []
    bond_j = []
    bond_k = []
    bond_r0 = []
    for a1, a2, itype in bonds:
        i = _atom_index_from_coord_index(int(a1))
        j = _atom_index_from_coord_index(int(a2))
        t = int(itype) - 1
        bond_i.append(i)
        bond_j.append(j)
        bond_k.append(float(bond_k_list[t]))
        bond_r0.append(float(bond_r0_list[t]))

    # ---- angles ----
    angles = list(_chunk(raw.get("ANGLES_INC_HYDROGEN", []), 4)) + list(
        _chunk(raw.get("ANGLES_WITHOUT_HYDROGEN", []), 4)
    )
    ang_k_list = raw.get("ANGLE_FORCE_CONSTANT", [])
    ang_t0_list = raw.get("ANGLE_EQUIL_VALUE", [])
    angle_i = []
    angle_j = []
    angle_k = []
    angle_k0 = []
    angle_t0 = []
    for a1, a2, a3, itype in angles:
        i = _atom_index_from_coord_index(int(a1))
        j = _atom_index_from_coord_index(int(a2))
        k = _atom_index_from_coord_index(int(a3))
        t = int(itype) - 1
        angle_i.append(i)
        angle_j.append(j)
        angle_k.append(k)
        angle_k0.append(float(ang_k_list[t]))
        angle_t0.append(float(ang_t0_list[t]))

    # ---- dihedrals ----
    dihedrals = list(_chunk(raw.get("DIHEDRALS_INC_HYDROGEN", []), 5)) + list(
        _chunk(raw.get("DIHEDRALS_WITHOUT_HYDROGEN", []), 5)
    )
    dih_force_list = raw.get("DIHEDRAL_FORCE_CONSTANT", [])
    dih_phase_list = raw.get("DIHEDRAL_PHASE", [])
    dih_per_list = raw.get("DIHEDRAL_PERIODICITY", [])

    dihed_i = []
    dihed_j = []
    dihed_k = []
    dihed_l = []
    dihed_force = []
    dihed_period = []
    dihed_phase = []
    for a1, a2, a3, a4, itype in dihedrals:
        t = int(itype) - 1
        dihed_i.append(_atom_index_from_coord_index(int(a1)))
        dihed_j.append(_atom_index_from_coord_index(int(a2)))
        dihed_k.append(_atom_index_from_coord_index(int(a3)))
        dihed_l.append(_atom_index_from_coord_index(int(a4)))
        dihed_force.append(float(dih_force_list[t]))
        dihed_period.append(float(dih_per_list[t]))
        dihed_phase.append(float(dih_phase_list[t]))

    # ---- exclusions & 1-4 ----
    excluded_keys = _build_excluded_pair_keys(
        natom,
        raw.get("NUMBER_EXCLUDED_ATOMS", [0] * natom),
        raw.get("EXCLUDED_ATOMS_LIST", []),
    )

    scee = raw.get("SCEE_SCALE_FACTOR", [])
    scnb = raw.get("SCNB_SCALE_FACTOR", [])
    if not scee or not scnb:
        # Defaults: AMBER typically uses 1.2 and 2.0 (i.e. scale 1/1.2 and 1/2.0)
        # but in prmtop these should exist. If absent, fall back to "no scaling".
        scee = [1.0] * len(dih_force_list)
        scnb = [1.0] * len(dih_force_list)

    d_inc = raw.get("DIHEDRALS_INC_HYDROGEN", [])
    d_wo = raw.get("DIHEDRALS_WITHOUT_HYDROGEN", [])
    p14_i, p14_j, inv_scee, inv_scnb = _build_14_pairs(natom, d_inc, d_wo, scee, scnb)
    keys14 = p14_i * natom + p14_j

    pair_i, pair_j = _build_general_pairs(natom, excluded_keys, keys14)
    (
        cmap_type,
        cmap_i,
        cmap_j,
        cmap_k,
        cmap_l,
        cmap_m,
        cmap_resolution,
        cmap_maps,
    ) = _parse_cmap(raw)
    cmap_size, cmap_delta, cmap_offset, cmap_coeff = _build_cmap_coeff_data(
        cmap_resolution=cmap_resolution,
        cmap_maps=cmap_maps,
    )

    device = torch.device(device)

    return AmberSystem(
        natom=natom,
        charge=torch.tensor(charge, dtype=torch.float64, device=device),
        atom_type=atom_type_t.to(device),
        bond_i=torch.tensor(bond_i, dtype=torch.int64, device=device),
        bond_j=torch.tensor(bond_j, dtype=torch.int64, device=device),
        bond_k=torch.tensor(bond_k, dtype=torch.float64, device=device),
        bond_r0=torch.tensor(bond_r0, dtype=torch.float64, device=device),
        angle_i=torch.tensor(angle_i, dtype=torch.int64, device=device),
        angle_j=torch.tensor(angle_j, dtype=torch.int64, device=device),
        angle_k=torch.tensor(angle_k, dtype=torch.int64, device=device),
        angle_k0=torch.tensor(angle_k0, dtype=torch.float64, device=device),
        angle_t0=torch.tensor(angle_t0, dtype=torch.float64, device=device),
        dihed_i=torch.tensor(dihed_i, dtype=torch.int64, device=device),
        dihed_j=torch.tensor(dihed_j, dtype=torch.int64, device=device),
        dihed_k=torch.tensor(dihed_k, dtype=torch.int64, device=device),
        dihed_l=torch.tensor(dihed_l, dtype=torch.int64, device=device),
        dihed_force=torch.tensor(dihed_force, dtype=torch.float64, device=device),
        dihed_period=torch.tensor(dihed_period, dtype=torch.float64, device=device),
        dihed_phase=torch.tensor(dihed_phase, dtype=torch.float64, device=device),
        lj_acoef=lj_acoef.to(device),
        lj_bcoef=lj_bcoef.to(device),
        hb_acoef=hb_acoef.to(device),
        hb_bcoef=hb_bcoef.to(device),
        nb_index=nb_index.to(device),
        pair_i=pair_i.to(device),
        pair_j=pair_j.to(device),
        pair14_i=p14_i.to(device),
        pair14_j=p14_j.to(device),
        pair14_inv_scee=inv_scee.to(device),
        pair14_inv_scnb=inv_scnb.to(device),
        cmap_type=cmap_type.to(device),
        cmap_i=cmap_i.to(device),
        cmap_j=cmap_j.to(device),
        cmap_k=cmap_k.to(device),
        cmap_l=cmap_l.to(device),
        cmap_m=cmap_m.to(device),
        cmap_resolution=cmap_resolution,
        cmap_maps=tuple(x.to(device) for x in cmap_maps),
        cmap_size=cmap_size.to(device),
        cmap_delta=cmap_delta.to(device),
        cmap_offset=cmap_offset.to(device),
        cmap_coeff=cmap_coeff.to(device),
    )


def load_coords(
    coords: Union[str, Path, torch.Tensor, Any],
    natom: Optional[int] = None,
    device: Union[str, torch.device] = "cpu",
    dtype: torch.dtype = torch.float64,
) -> torch.Tensor:
    """Load coordinates (Angstrom).

    Supported formats:
      - Amber restart / inpcrd (ASCII)
      - PDB (ATOM/HETATM records; uses X,Y,Z columns)
      - XYZ
      - ASE Atoms object
      - torch.Tensor with shape [N,3]
    """
    dev = torch.device(device)

    if torch.is_tensor(coords):
        xyz_t = coords.detach().to(device=dev, dtype=dtype)
        if xyz_t.ndim != 2 or int(xyz_t.shape[1]) != 3:
            raise ValueError(f"coords tensor must have shape [N,3], got {tuple(xyz_t.shape)}")
        if natom is not None and int(xyz_t.shape[0]) != int(natom):
            raise ValueError(f"NATOM mismatch: coords has {int(xyz_t.shape[0])}, expected {natom}")
        return xyz_t

    # ASE Atoms (duck-typed): object exposing get_positions() -> [N,3]
    if hasattr(coords, "get_positions") and callable(getattr(coords, "get_positions")):
        pos = coords.get_positions()
        xyz_t = torch.as_tensor(pos, dtype=dtype, device=dev)
        if xyz_t.ndim != 2 or int(xyz_t.shape[1]) != 3:
            raise ValueError(f"ASE positions must have shape [N,3], got {tuple(xyz_t.shape)}")
        if natom is not None and int(xyz_t.shape[0]) != int(natom):
            raise ValueError(f"NATOM mismatch: coords has {int(xyz_t.shape[0])}, expected {natom}")
        return xyz_t

    coords_path = Path(coords)
    suffix = coords_path.suffix.lower()
    if suffix in {".rst7", ".inpcrd", ".crd", ".restrt"}:
        xyz = _read_amber_inpcrd(coords_path, natom=natom)
    elif suffix == ".pdb":
        xyz = _read_pdb(coords_path)
    elif suffix == ".xyz":
        xyz = _read_xyz(coords_path, natom=natom)
    else:
        raise ValueError(
            f"Unsupported coordinate input: {coords!r}. "
            "Supported: rst7/inpcrd/crd/restrt, pdb, xyz, ASE Atoms, torch.Tensor[N,3]."
        )
    if natom is not None and len(xyz) != int(natom):
        raise ValueError(f"NATOM mismatch: coords has {len(xyz)}, expected {natom}")
    return torch.tensor(xyz, dtype=dtype, device=dev)


def _read_amber_inpcrd(path: Path, natom: Optional[int] = None) -> List[List[float]]:
    """Very small reader for Amber ASCII restart/inpcrd coordinate files."""
    raw_lines = path.read_text(errors="replace").splitlines()
    if len(raw_lines) < 2:
        raise ValueError(f"File too short: {path}")
    # 2nd line begins with NATOM
    nat = int(raw_lines[1].split()[0])
    if natom is not None and nat != natom:
        raise ValueError(f"NATOM mismatch: coords has {nat}, expected {natom}")
    # Collect floats from remaining lines.
    # Amber ASCII inpcrd/rst7 commonly uses fixed-width (typically 6E12.7 per line).
    # Parsing as fixed-width is robust even if there are spaces.
    floats: List[float] = []
    for ln in raw_lines[2:]:
        for i in range(0, len(ln), 12):
            field = ln[i : i + 12].strip()
            if not field:
                continue
            try:
                floats.append(float(field))
            except ValueError:
                continue
    need = 3 * nat
    if len(floats) < need:
        raise ValueError(f"Not enough coordinate floats in {path}: got {len(floats)}, need {need}")
    floats = floats[:need]
    return [[floats[3 * i], floats[3 * i + 1], floats[3 * i + 2]] for i in range(nat)]


def _read_pdb(path: Path) -> List[List[float]]:
    xyz: List[List[float]] = []
    for ln in path.read_text(errors="replace").splitlines():
        if ln.startswith("ATOM") or ln.startswith("HETATM"):
            x = float(ln[30:38])
            y = float(ln[38:46])
            z = float(ln[46:54])
            xyz.append([x, y, z])
    if not xyz:
        raise ValueError(f"No ATOM/HETATM coordinates found in {path}")
    return xyz


def _read_xyz(path: Path, natom: Optional[int] = None) -> List[List[float]]:
    """Read first frame from XYZ.

    Supports standard XYZ (atom count + comment header) and simple
    whitespace-separated rows with either ``sym x y z`` or ``x y z``.
    """

    def _parse_xyz_row(raw_line: str) -> List[float]:
        toks = raw_line.split()
        if len(toks) < 3:
            raise ValueError(f"Invalid XYZ row: {raw_line!r}")
        for start in (1, 0):
            if len(toks) < start + 3:
                continue
            try:
                x = float(toks[start + 0])
                y = float(toks[start + 1])
                z = float(toks[start + 2])
                return [x, y, z]
            except ValueError:
                continue
        raise ValueError(f"Invalid XYZ row: {raw_line!r}")

    lines = [ln.rstrip("\n") for ln in path.read_text(errors="replace").splitlines()]
    if not lines:
        raise ValueError(f"Empty XYZ file: {path}")

    xyz: List[List[float]] = []
    header_natom: Optional[int] = None
    try:
        header_natom = int(lines[0].strip().split()[0])
    except Exception:
        header_natom = None

    if header_natom is not None:
        if len(lines) < 2 + header_natom:
            raise ValueError(
                f"XYZ header count mismatch in {path}: need {header_natom} atoms, "
                f"but file has {max(0, len(lines) - 2)} rows"
            )
        for ln in lines[2 : 2 + header_natom]:
            if not ln.strip():
                continue
            xyz.append(_parse_xyz_row(ln))
    else:
        for ln in lines:
            if not ln.strip():
                continue
            xyz.append(_parse_xyz_row(ln))

    if not xyz:
        raise ValueError(f"No coordinates found in XYZ file: {path}")
    if natom is not None and len(xyz) != int(natom):
        raise ValueError(f"NATOM mismatch: coords has {len(xyz)}, expected {natom}")
    return xyz
