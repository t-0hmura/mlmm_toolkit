# mlmm_toolkit/oniom_export.py

"""
Export ML/MM system to Gaussian/ORCA ONIOM input format from Amber parm7 topology.

Example:
    mlmm oniom-gaussian --parm real.parm7 -i pocket.pdb --model-pdb ml.pdb -o out.com

For detailed documentation, see: docs/oniom_export.md
"""

from __future__ import annotations

import re
import shlex
import shutil
import subprocess
from pathlib import Path
from typing import Any, Dict, List, Optional, Set, Tuple

import logging

import click
import numpy as np

from .add_elem_info import guess_element as _guess_element

logger = logging.getLogger(__name__)

try:
    import parmed as pmd
except ImportError:
    pmd = None

_GAUSSIAN_DEFAULT_METHOD = "wB97XD/def2-TZVPD"
_ORCA_DEFAULT_METHOD = "B3LYP D3BJ def2-SVP"


def _check_parmed() -> None:
    """Check if ParmEd is available."""
    if pmd is None:
        raise ImportError(
            "ParmEd is required for ONIOM export. Install with: pip install parmed"
        )


# -----------------------------------------------
# Coordinates / element handling
# -----------------------------------------------

# Periodic table symbols (1-indexed; index 0 is dummy)
_PERIODIC_TABLE: List[str] = [
    "",
    "H",
    "He",
    "Li",
    "Be",
    "B",
    "C",
    "N",
    "O",
    "F",
    "Ne",
    "Na",
    "Mg",
    "Al",
    "Si",
    "P",
    "S",
    "Cl",
    "Ar",
    "K",
    "Ca",
    "Sc",
    "Ti",
    "V",
    "Cr",
    "Mn",
    "Fe",
    "Co",
    "Ni",
    "Cu",
    "Zn",
    "Ga",
    "Ge",
    "As",
    "Se",
    "Br",
    "Kr",
    "Rb",
    "Sr",
    "Y",
    "Zr",
    "Nb",
    "Mo",
    "Tc",
    "Ru",
    "Rh",
    "Pd",
    "Ag",
    "Cd",
    "In",
    "Sn",
    "Sb",
    "Te",
    "I",
    "Xe",
    "Cs",
    "Ba",
    "La",
    "Ce",
    "Pr",
    "Nd",
    "Pm",
    "Sm",
    "Eu",
    "Gd",
    "Tb",
    "Dy",
    "Ho",
    "Er",
    "Tm",
    "Yb",
    "Lu",
    "Hf",
    "Ta",
    "W",
    "Re",
    "Os",
    "Ir",
    "Pt",
    "Au",
    "Hg",
    "Tl",
    "Pb",
    "Bi",
    "Po",
    "At",
    "Rn",
    "Fr",
    "Ra",
    "Ac",
    "Th",
    "Pa",
    "U",
    "Np",
    "Pu",
    "Am",
    "Cm",
    "Bk",
    "Cf",
    "Es",
    "Fm",
    "Md",
    "No",
    "Lr",
    "Rf",
    "Db",
    "Sg",
    "Bh",
    "Hs",
    "Mt",
    "Ds",
    "Rg",
    "Cn",
    "Nh",
    "Fl",
    "Mc",
    "Lv",
    "Ts",
    "Og",
]

_TWO_LETTER_ELEMENT_UPPER: Set[str] = {sym.upper() for sym in _PERIODIC_TABLE if len(sym) == 2}

# A small set of common atomic masses for robust element inference when atomic number is missing.
# (Atomic masses vary slightly by isotope; we only need a reasonable guess.)
_COMMON_MASS_TABLE: List[Tuple[str, float]] = [
    ("H", 1.008),
    ("C", 12.011),
    ("N", 14.007),
    ("O", 15.999),
    ("F", 18.998),
    ("Na", 22.990),
    ("Mg", 24.305),
    ("Al", 26.982),
    ("Si", 28.085),
    ("P", 30.974),
    ("S", 32.06),
    ("Cl", 35.45),
    ("K", 39.098),
    ("Ca", 40.078),
    ("Mn", 54.938),
    ("Fe", 55.845),
    ("Co", 58.933),
    ("Ni", 58.693),
    ("Cu", 63.546),
    ("Zn", 65.38),
    ("Se", 78.971),
    ("Br", 79.904),
    ("I", 126.90),
]


def _normalize_element_symbol(sym: str) -> str:
    """
    Normalize an element symbol to canonical form (e.g., 'cl'/'CL' -> 'Cl').

    Returns 'X' if empty / unknown.
    """
    s = (sym or "").strip()
    if not s:
        return "X"
    # Keep only first 2 characters (typical element symbols)
    s = re.sub(r"[^A-Za-z]", "", s)
    if not s:
        return "X"
    if len(s) == 1:
        return s.upper()
    return s[0].upper() + s[1].lower()


def _infer_element_from_pdb_atom_name(atom_name_field: str) -> str:
    """
    Best-effort element inference from the 4-char PDB atom-name field (columns 13-16).

    This uses PDB alignment conventions:
    - one-letter elements are right-justified (often leading space), e.g. " CA " (C-alpha) -> C
    - two-letter elements are left-justified, e.g. "CA  " (calcium) -> Ca
    """
    field = (atom_name_field or "")[:4]
    if len(field) < 2:
        return _normalize_element_symbol(field.strip())

    # Example: "1H  " -> element "H"
    if field[0].isdigit():
        return _normalize_element_symbol(field[1])

    # Right-justified => one-letter element in column 14
    if field[0] == " ":
        return _normalize_element_symbol(field[1])

    # Left-justified => likely 2-letter element (or 1-letter + suffix)
    cand2 = field[0:2].strip().upper()
    if cand2 in _TWO_LETTER_ELEMENT_UPPER:
        return _normalize_element_symbol(cand2)

    return _normalize_element_symbol(field[0])


def _read_pdb_geometry(pdb_path: Path) -> Tuple[np.ndarray, List[str]]:
    """Read coordinates and element symbols from a PDB file (ATOM/HETATM records)."""
    coords: List[List[float]] = []
    elements: List[str] = []
    with pdb_path.open("r") as f:
        for line in f:
            if not (line.startswith("ATOM") or line.startswith("HETATM")):
                continue
            try:
                x = float(line[30:38])
                y = float(line[38:46])
                z = float(line[46:54])
            except Exception as e:
                raise ValueError(f"Failed to parse coordinates from PDB line: {line.rstrip()}") from e

            # Element: prefer columns 77-78; fall back to residue-aware guess_element
            # (which correctly handles HG2->H in protein residues, etc.).
            elem_field = line[76:78].strip()
            if elem_field:
                elem = elem_field
                # Guard against 1-letter misalignment when atom name encodes a 2-letter element (e.g., MG)
                elem_inferred = _infer_element_from_pdb_atom_name(line[12:16])
                if len(elem_field) == 1 and elem_inferred and len(elem_inferred) == 2:
                    if elem_field.upper() != elem_inferred[0].upper():
                        elem = elem_inferred
            else:
                # No element column — use residue-aware inference (add_elem_info.guess_element)
                atom_name = line[12:16].strip()
                resname = line[17:20].strip()
                is_het = line.startswith("HETATM")
                guessed = _guess_element(atom_name, resname, is_het)
                if guessed:
                    elem = guessed
                else:
                    elem = _infer_element_from_pdb_atom_name(line[12:16])

            coords.append([x, y, z])
            elements.append(_normalize_element_symbol(elem))

    if not coords:
        raise ValueError(f"No ATOM/HETATM records found in PDB: {pdb_path}")
    return np.asarray(coords, dtype=float), elements


def _read_xyz_geometry(xyz_path: Path) -> Tuple[np.ndarray, List[str]]:
    """Read coordinates and element symbols from a single-frame XYZ file."""
    coords: List[List[float]] = []
    elements: List[str] = []
    with xyz_path.open("r") as f:
        lines = f.readlines()

    if len(lines) < 3:
        raise ValueError(f"XYZ file is too short: {xyz_path}")

    try:
        n_atoms = int(lines[0].strip())
    except Exception as e:
        raise ValueError(f"First line of XYZ must be an integer atom count: {xyz_path}") from e

    if len(lines) < 2 + n_atoms:
        raise ValueError(
            f"XYZ file atom count ({n_atoms}) exceeds available lines: {xyz_path}"
        )

    for i in range(n_atoms):
        raw = lines[2 + i].strip()
        if not raw:
            continue
        parts = raw.split()
        if len(parts) < 4:
            raise ValueError(f"Invalid XYZ atom line: '{raw}'")
        elem = _normalize_element_symbol(parts[0])
        try:
            x, y, z = float(parts[1]), float(parts[2]), float(parts[3])
        except Exception as e:
            raise ValueError(f"Invalid XYZ coordinates in line: '{raw}'") from e
        coords.append([x, y, z])
        elements.append(elem)

    if len(coords) != n_atoms:
        raise ValueError(
            f"XYZ parsing produced {len(coords)} atoms but header says {n_atoms}: {xyz_path}"
        )

    return np.asarray(coords, dtype=float), elements


def _read_input_geometry(input_path: Path) -> Tuple[np.ndarray, List[str]]:
    """
    Read coordinates + element list from an input coordinate file.

    Supported:
      - .pdb
      - .xyz
    """
    suffix = input_path.suffix.lower()
    if suffix == ".pdb" or suffix == ".ent":
        return _read_pdb_geometry(input_path)
    if suffix == ".xyz":
        return _read_xyz_geometry(input_path)
    raise ValueError(f"Unsupported input coordinate format: {input_path} (expected .pdb or .xyz)")


def _apply_coordinates_to_parm(parm, coords: np.ndarray) -> None:
    """Attach coordinates to a ParmEd structure, keeping atom order unchanged."""
    coords = np.asarray(coords, dtype=float)
    if coords.shape != (len(parm.atoms), 3):
        raise ValueError(
            f"Atom count mismatch: parm7 has {len(parm.atoms)} atoms, "
            f"but coordinate file has {coords.shape[0]} atoms"
        )

    # Prefer setting structure-level coordinates
    try:
        parm.coordinates = coords
    except Exception:
        logger.debug("Failed to set structure-level coordinates on parm", exc_info=True)

    # Ensure per-atom cached coordinates are also available
    for i, atom in enumerate(parm.atoms):
        x, y, z = coords[i]
        try:
            atom.xx = float(x)
            atom.xy = float(y)
            atom.xz = float(z)
        except Exception:
            logger.debug("Failed to set per-atom coords on atom %d", i, exc_info=True)


def _infer_element_from_mass(mass: float, tol: float = 1.5) -> str:
    """
    Infer element from atomic mass using a small common-element table.

    Returns 'X' if no close match found.
    """
    try:
        m = float(mass)
    except Exception:
        return "X"

    best_sym = "X"
    best_diff = 1e9
    for sym, ref_m in _COMMON_MASS_TABLE:
        diff = abs(m - ref_m)
        if diff < best_diff:
            best_diff = diff
            best_sym = sym

    if best_diff <= tol:
        return _normalize_element_symbol(best_sym)
    return "X"


def _get_parm_element(atom) -> str:
    """
    Best-effort element symbol for a ParmEd atom.

    Tries (in order):
      1) atom.element_name
      2) atom.atomic_number -> periodic table
      3) atom.mass -> common mass table
      4) atom.name PDB-style heuristic (very weak fallback)
    """
    # 1) element_name
    elem = getattr(atom, "element_name", None)
    if elem:
        norm = _normalize_element_symbol(str(elem))
        if norm != "X":
            return norm

    # 2) atomic_number
    z = getattr(atom, "atomic_number", None)
    if z is not None:
        try:
            zi = int(z)
            if 0 < zi < len(_PERIODIC_TABLE):
                return _PERIODIC_TABLE[zi]
        except Exception:
            logger.debug("Failed to infer element from atomic_number=%s", z, exc_info=True)

    # 3) mass
    mass = getattr(atom, "mass", None)
    if mass is not None:
        guess = _infer_element_from_mass(mass)
        if guess != "X":
            return guess

    # 4) fallback: first letter of atom name (can be wrong for metals)
    name = getattr(atom, "name", "")
    if name:
        return _normalize_element_symbol(str(name)[0])

    return "X"


def _get_parm_elements(parm) -> List[str]:
    """Element symbols for each atom in a ParmEd structure."""
    return [_get_parm_element(atom) for atom in parm.atoms]


def _validate_element_order(
    parm,
    input_elements: List[str],
    *,
    strict: bool = True,
) -> None:
    """
    Validate that element sequence in the coordinate file matches the parm7 topology.

    If `strict` is True, raise on the first detected mismatch where both elements are known.
    Unknown elements ('X') are ignored.
    """
    parm_elements = _get_parm_elements(parm)
    if len(parm_elements) != len(input_elements):
        raise ValueError(
            f"Atom count mismatch: parm7 has {len(parm_elements)} atoms, "
            f"but input has {len(input_elements)} atoms"
        )

    unknown_in = sum(1 for e in input_elements if e == "X")
    unknown_parm = sum(1 for e in parm_elements if e == "X")
    if unknown_in > 0 or unknown_parm > 0:
        click.echo(
            f"[oniom-export] WARNING: element check is partial "
            f"(unknown elements: input={unknown_in}, parm={unknown_parm})"
        )

    for i, (e_parm, e_in) in enumerate(zip(parm_elements, input_elements)):
        if e_parm == "X" or e_in == "X":
            continue
        if e_parm != e_in:
            msg = (
                f"Element sequence mismatch at atom index {i} (0-based): "
                f"parm7={e_parm}, input={e_in}. "
                f"Atom order likely differs between parm7 and the coordinate file."
            )
            if strict:
                raise ValueError(msg)
            click.echo(f"[oniom-export] WARNING: {msg}")


def _get_total_charge(parm) -> int:
    """
    Calculate total charge (integer) from atom partial charges.

    Notes
    -----
    - `atom.charge` is the most reliable source (units of electron charge).
    - Amber prmtop stores CHARGE values scaled by ~18.2223. If we ever fall back
      to parm_data["CHARGE"], we divide by 18.2223.
    """
    q: float
    try:
        q = float(sum(float(getattr(a, "charge", 0.0)) for a in parm.atoms))
    except Exception:
        # Fallback: try prmtop raw charges (often scaled)
        try:
            q_raw = float(np.sum(parm.parm_data["CHARGE"]))
            q = q_raw / 18.2223
        except Exception:
            q = 0.0

    q_int = int(round(q))
    if abs(q - q_int) > 1e-3:
        click.echo(
            f"[oniom-export] WARNING: total charge {q:.6f} is not close to an integer; "
            f"rounded to {q_int}"
        )
    return q_int


def _fix_atom_type(atom_type: str) -> str:
    """
    Fix atom types for Gaussian compatibility.

    - 2C, 3C -> C2C, C3C (numeric prefix)
    - C*, N* -> C9, N9 (asterisk)
    - lowercase (GAFF2) -> L{uppercase} (e.g., ca -> LCA)
    """
    atom_type = str(atom_type)
    if atom_type == "2C":
        return "C2C"
    elif atom_type == "3C":
        return "C3C"
    elif atom_type == "C*":
        return "C9"
    elif atom_type == "N*":
        return "N9"
    elif bool(re.match(r"^[a-z]+", atom_type)):
        return f"L{atom_type.upper()}"
    else:
        return atom_type


def _parse_pdb_atoms_with_meta(pdb_path: Path) -> List[Dict[str, Any]]:
    """
    Parse ATOM/HETATM records from a PDB file.

    Returns a list of dictionaries:
      - idx (0-based, sequential in file)
      - atom_name, res_name, chain_id, res_seq, icode
      - coord (np.ndarray shape (3,))
      - element (best-effort)
      - bfactor (float, defaults to 0.0)
    """
    atoms: List[Dict[str, Any]] = []
    with pdb_path.open("r") as f:
        atom_idx = 0
        for line in f:
            if not line.startswith(("ATOM", "HETATM")):
                continue

            atom_name = line[12:16].strip()
            res_name = line[17:20].strip()
            chain_id = line[21:22].strip()
            res_seq_str = line[22:26].strip()
            icode = line[26:27].strip()

            try:
                res_seq = int(res_seq_str)
            except Exception:
                res_seq = 0

            try:
                x = float(line[30:38])
                y = float(line[38:46])
                z = float(line[46:54])
            except Exception:
                x, y, z = 0.0, 0.0, 0.0

            # B-factor (tempFactor) is columns 61-66 in PDB v3.3 (0-based slice 60:66)
            try:
                bfac = float(line[60:66])
            except Exception:
                bfac = 0.0

            elem_field = line[76:78].strip() if len(line) >= 78 else ""
            if elem_field:
                elem = elem_field
            else:
                is_het = line.startswith("HETATM")
                guessed = _guess_element(atom_name, res_name, is_het)
                elem = guessed if guessed else _infer_element_from_pdb_atom_name(line[12:16])

            atoms.append(
                {
                    "idx": atom_idx,
                    "atom_name": atom_name,
                    "res_name": res_name,
                    "chain_id": chain_id,
                    "res_seq": res_seq,
                    "icode": icode,
                    "coord": np.array([x, y, z], dtype=float),
                    "element": _normalize_element_symbol(elem),
                    "bfactor": float(bfac),
                }
            )
            atom_idx += 1
    return atoms


def _read_qm_atoms_from_pdb(
    model_pdb: Path,
    *,
    input_pdb: Optional[Path] = None,
    system_coords: Optional[np.ndarray] = None,
    system_elements: Optional[List[str]] = None,
    match_tol: float = 0.2,
) -> Set[int]:
    """
    Determine QM-region atom indices (0-based, topology order).

    This function tries (in order):
      1) Match model_pdb atoms to atoms in `input_pdb` via (atom_name, res_name, res_seq),
         with coordinate-based disambiguation when that ID is not unique.
      2) If `input_pdb` is not available, match by nearest coordinates against `system_coords`
         (optionally requiring element agreement when available).

    Parameters
    ----------
    model_pdb
        PDB containing QM-region atoms (typically a subset PDB produced by `define-layer`
        or `build_model_pdb_from_bfactors`).
    input_pdb
        Full-system PDB whose atom order matches the Amber topology (recommended).
    system_coords
        Full-system coordinates (shape (N,3)) in topology order.
    system_elements
        Optional element list (len N). If provided, element mismatches are rejected.
    match_tol
        Maximum allowed distance (Å) for coordinate matching/disambiguation.

    Returns
    -------
    Set[int]
        QM atom indices in 0-based topology order.
    """
    qm_indices: Set[int] = set()

    model_atoms = _parse_pdb_atoms_with_meta(model_pdb)
    if not model_atoms:
        return set()

    # Path 1: match by PDB identifiers using full-system PDB
    if input_pdb is not None and input_pdb.suffix.lower() in {".pdb", ".ent"} and input_pdb.exists():
        input_atoms = _parse_pdb_atoms_with_meta(input_pdb)
        if not input_atoms:
            raise ValueError(f"Failed to read any atoms from input PDB: {input_pdb}")

        if system_coords is None:
            # If the caller didn't provide coords, use the PDB coords for disambiguation
            system_coords = np.asarray([a["coord"] for a in input_atoms], dtype=float)

        # Map from (atom_name,res_name,res_seq) -> list of candidate indices
        key_to_candidates: Dict[Tuple[str, str, int], List[int]] = {}
        for a in input_atoms:
            key = (a["atom_name"], a["res_name"], int(a["res_seq"]))
            key_to_candidates.setdefault(key, []).append(int(a["idx"]))

        used: Set[int] = set()
        sys_coords = np.asarray(system_coords, dtype=float)

        missing: int = 0
        for ma in model_atoms:
            key = (ma["atom_name"], ma["res_name"], int(ma["res_seq"]))
            cand = key_to_candidates.get(key, [])
            if not cand:
                missing += 1
                continue

            if len(cand) == 1:
                chosen = cand[0]
                if chosen in used:
                    # Already used (duplicate identifiers). Fall back to coordinate disambiguation.
                    cand = cand
                else:
                    qm_indices.add(chosen)
                    used.add(chosen)
                    continue

            # Disambiguate by nearest coordinate among candidates (and avoid already-used indices)
            cand_free = [i for i in cand if i not in used]
            if not cand_free:
                cand_free = cand  # allow reuse as a last resort

            cand_coords = sys_coords[np.asarray(cand_free, dtype=int)]
            dists = np.linalg.norm(cand_coords - ma["coord"][None, :], axis=1)
            j = int(np.argmin(dists))
            chosen = int(cand_free[j])

            if float(dists[j]) > match_tol:
                # ID match exists but coordinates are far apart -> likely inconsistent inputs
                click.echo(
                    f"[oniom-export] WARNING: matched ID {key} but nearest distance is {dists[j]:.3f} Å "
                    f"(> {match_tol} Å). Check that input_pdb and model_pdb come from the same structure."
                )

            qm_indices.add(chosen)
            used.add(chosen)

        if missing > 0:
            click.echo(
                f"[oniom-export] WARNING: {missing} atoms in model_pdb could not be matched by "
                f"(atom_name,res_name,res_seq) to input_pdb. "
                "If this is unexpected, verify residue numbering and naming."
            )

        return qm_indices

    # Path 2: coordinate-only matching against system_coords
    if system_coords is None:
        raise ValueError(
            "Cannot match model_pdb to topology atoms without either `input_pdb` (full-system PDB) "
            "or `system_coords`."
        )

    sys_coords = np.asarray(system_coords, dtype=float)
    try:
        from scipy.spatial import cKDTree

        tree = cKDTree(sys_coords)
        for ma in model_atoms:
            dist, idx = tree.query(ma["coord"], k=1)
            if float(dist) > match_tol:
                continue
            if system_elements is not None and 0 <= int(idx) < len(system_elements):
                e_sys = _normalize_element_symbol(system_elements[int(idx)])
                e_mod = _normalize_element_symbol(ma["element"])
                if e_sys != "X" and e_mod != "X" and e_sys != e_mod:
                    continue
            qm_indices.add(int(idx))
    except Exception:
        # Slow fallback
        for ma in model_atoms:
            d = np.linalg.norm(sys_coords - ma["coord"][None, :], axis=1)
            idx = int(np.argmin(d))
            if float(d[idx]) > match_tol:
                continue
            if system_elements is not None and 0 <= idx < len(system_elements):
                e_sys = _normalize_element_symbol(system_elements[idx])
                e_mod = _normalize_element_symbol(ma["element"])
                if e_sys != "X" and e_mod != "X" and e_sys != e_mod:
                    continue
            qm_indices.add(idx)

    return qm_indices

def _identify_qm_atoms_by_distance(
    parm,
    qm_residue_indices: List[int],
    near_cutoff: float,
) -> Tuple[Set[int], Set[int]]:
    """
    Identify QM and movable atoms based on residue indices and distance cutoff.

    Returns:
        (qm_atom_indices, movable_atom_indices) - both 0-based
    """
    from scipy import spatial

    # Get QM atom indices from specified residues
    qm_atom_indices: Set[int] = set()
    for resi in qm_residue_indices:
        if 0 <= resi < len(parm.residues):
            for atom in parm.residues[resi].atoms:
                qm_atom_indices.add(atom.idx)

    if not qm_atom_indices:
        return set(), set()

    # Find movable atoms (within near_cutoff of QM atoms)
    qm_list = sorted(qm_atom_indices)
    neighbor_mask = np.any(
        spatial.distance.cdist(parm.coordinates, parm.coordinates[qm_list]) <= near_cutoff,
        axis=1,
    )

    # Include entire residues if any atom is within cutoff
    movable_indices: Set[int] = set()
    neighbor_residues = set(parm.atoms[i].residue for i in np.where(neighbor_mask)[0])
    for residue in neighbor_residues:
        for atom in residue.atoms:
            movable_indices.add(atom.idx)

    return qm_atom_indices, movable_indices


# -----------------------------------------------
# QM/MM covalent-boundary link helpers
# -----------------------------------------------

_LINK_H_BOND_LENGTH = {
    "C": 1.09,
    "N": 1.01,
}

_LINK_H_FF_TYPE = {
    "C": "HC",
    "N": "H",
}


def _atom_xyz(parm, atom_idx: int) -> np.ndarray:
    """Return Cartesian coordinate (Å) for a topology atom index."""
    atom = parm.atoms[int(atom_idx)]
    try:
        return np.array([float(atom.xx), float(atom.xy), float(atom.xz)], dtype=float)
    except Exception:
        return np.asarray(parm.coordinates[int(atom_idx)], dtype=float)


def _find_qmmm_boundary_pairs(parm, qm_indices: Set[int]) -> List[Tuple[int, int]]:
    """
    Detect covalent QM/MM boundary bonds from topology bonds.

    Returns
    -------
    List[Tuple[int, int]]
        A list of (qm_idx, mm_idx) index pairs (0-based).
    """
    per_mm_candidates: Dict[int, List[int]] = {}

    for bond in getattr(parm, "bonds", []):
        i = int(bond.atom1.idx)
        j = int(bond.atom2.idx)
        i_qm = i in qm_indices
        j_qm = j in qm_indices
        if i_qm == j_qm:
            continue
        qm_idx, mm_idx = (i, j) if i_qm else (j, i)
        per_mm_candidates.setdefault(mm_idx, []).append(qm_idx)

    pairs: List[Tuple[int, int]] = []
    for mm_idx, cands_raw in sorted(per_mm_candidates.items()):
        cands = sorted(set(int(x) for x in cands_raw))
        if len(cands) == 1:
            pairs.append((cands[0], mm_idx))
            continue

        mm_xyz = _atom_xyz(parm, mm_idx)
        best_qm = min(
            cands,
            key=lambda q: float(np.linalg.norm(_atom_xyz(parm, int(q)) - mm_xyz)),
        )
        click.echo(
            f"[oniom-export] WARNING: MM atom {mm_idx} is bonded to multiple QM atoms "
            f"{cands}; using closest QM atom {best_qm} as the link parent."
        )
        pairs.append((int(best_qm), int(mm_idx)))

    return pairs


def _estimate_link_h_position(parm, qm_idx: int, mm_idx: int, bond_length: float) -> Optional[np.ndarray]:
    """Estimate link-H position using the MLMMCore rule: r_H = r_QM + u(QM->MM) * d."""
    qm_xyz = _atom_xyz(parm, qm_idx)
    mm_xyz = _atom_xyz(parm, mm_idx)
    vec = mm_xyz - qm_xyz
    norm = float(np.linalg.norm(vec))
    if norm < 1.0e-12:
        return None
    return qm_xyz + (vec / norm) * float(bond_length)


def _build_link_atom_specs(
    parm,
    qm_indices: Set[int],
    *,
    elements: Optional[List[str]] = None,
) -> Dict[int, Dict[str, Any]]:
    """
    Build link-atom specs keyed by boundary MM atom index.

    Each value contains:
      - qm_idx: 0-based QM parent index
      - ff_type: Gaussian Amber atom type for the link hydrogen
      - bond_length: QM-H bond length used for placement
      - position: estimated link-H coordinate (Å), when available
    """
    specs: Dict[int, Dict[str, Any]] = {}
    warned_elems: Set[str] = set()

    for qm_idx, mm_idx in _find_qmmm_boundary_pairs(parm, qm_indices):
        if elements is not None and 0 <= qm_idx < len(elements):
            qm_elem = _normalize_element_symbol(elements[qm_idx])
        else:
            qm_elem = _normalize_element_symbol(_get_parm_element(parm.atoms[qm_idx]))

        if qm_elem in _LINK_H_BOND_LENGTH:
            bond_len = float(_LINK_H_BOND_LENGTH[qm_elem])
        else:
            bond_len = float(_LINK_H_BOND_LENGTH["C"])
            if qm_elem not in warned_elems:
                click.echo(
                    f"[oniom-export] WARNING: unsupported QM parent element '{qm_elem}' for link-H "
                    f"distance; using C-like default ({bond_len:.2f} Å)."
                )
                warned_elems.add(qm_elem)

        ff_type = _LINK_H_FF_TYPE.get(qm_elem, _LINK_H_FF_TYPE["C"])
        link_pos = _estimate_link_h_position(parm, qm_idx=qm_idx, mm_idx=mm_idx, bond_length=bond_len)
        if link_pos is None:
            click.echo(
                f"[oniom-export] WARNING: failed to estimate link-H position for boundary "
                f"(QM={qm_idx}, MM={mm_idx}); skipping link annotation for this bond."
            )
            continue

        specs[int(mm_idx)] = {
            "qm_idx": int(qm_idx),
            "ff_type": str(ff_type),
            "bond_length": float(bond_len),
            "position": link_pos,
        }

    return specs


# -----------------------------------------------
# Gaussian ONIOM Export
# -----------------------------------------------

def _write_gaussian_header(
    parm,
    parm_path: str,
    output_name: str,
    method: str = "wB97XD/def2-TZVPD",
    nproc: int = 8,
    mem: str = "16GB",
    qm_charge: int = 0,
    qm_mult: int = 1,
    real_charge: Optional[int] = None,
    real_mult: Optional[int] = None,
) -> str:
    """
    Generate Gaussian ONIOM input header.

    Gaussian ONIOM uses *three* charge/multiplicity pairs for a 2-layer ONIOM job:
      (real system @ low level)  (model system @ high level)  (model system @ low level)

    We default the real-system charge to the total charge of the topology, and the real-system
    multiplicity to `qm_mult` (since the MM region is typically closed-shell).
    """
    total_charge = _get_total_charge(parm)

    if real_charge is None:
        real_charge = total_charge
    if real_mult is None:
        real_mult = qm_mult

    if int(real_charge) != int(total_charge):
        click.echo(
            f"[oniom-export] WARNING: real_charge={real_charge} differs from topology total charge={total_charge}. "
            "Proceeding as requested."
        )

    chk_name = Path(output_name).stem

    header = f"""%chk={chk_name}.chk
%mem={mem}
%nprocshared={nproc}
#p oniom({method}:amber=softonly)
scf=(xqc,intrep,maxconventionalcyc=80)
nosymm iop(2/15=3) geom=connectivity Amber=(FirstEquiv)

ONIOM inputfile generated by mlmm oniom-export from {parm_path}.

{real_charge} {real_mult} {qm_charge} {qm_mult} {qm_charge} {qm_mult}
"""
    return header

def _write_gaussian_coordinates(
    parm,
    qm_indices: Set[int],
    movable_indices: Set[int],
    elements: Optional[List[str]] = None,
    link_specs: Optional[Dict[int, Dict[str, Any]]] = None,
) -> Tuple[str, str]:
    """
    Generate Gaussian ONIOM coordinate section and connectivity.

    Returns:
        (coords_section, connectivity_section)
    """
    coords_lines: List[str] = []

    # Optional elements list (preferred when coming from input PDB/XYZ)
    elements_list: Optional[List[str]] = None
    if elements is not None and len(elements) == len(parm.atoms):
        elements_list = elements

    for atom in parm.atoms:
        idx = atom.idx
        layer = "H" if idx in qm_indices else "L"
        movable = 0 if idx in movable_indices else -1

        x, y, z = atom.xx, atom.xy, atom.xz
        ff_type = _fix_atom_type(atom.atom_type)
        charge = atom.charge

        if elements_list is not None:
            element = elements_list[idx]
        else:
            element = _get_parm_element(atom)

        atom_section = f"{element}-{ff_type}-{charge:.6f}"
        link_suffix = ""
        if layer == "L" and link_specs is not None and idx in link_specs:
            spec = link_specs[idx]
            qm_parent = int(spec["qm_idx"]) + 1  # Gaussian connectivity is 1-based
            link_ff_type = _fix_atom_type(str(spec["ff_type"]))
            link_suffix = f" H-{link_ff_type} {qm_parent}"
        coords_lines.append(
            f"{atom_section:<20} {movable:>2}   {x:12.6f} {y:12.6f} {z:12.6f} {layer}{link_suffix}"
        )

    # Connectivity section
    bond_dict: Dict[int, List[int]] = {}
    for bond in parm.bonds:
        i, j = bond.atom1.idx, bond.atom2.idx
        if i not in bond_dict:
            bond_dict[i] = []
        if j not in bond_dict:
            bond_dict[j] = []
        bond_dict[i].append(j)
        bond_dict[j].append(i)

    conn_lines: List[str] = []
    for i in range(len(parm.atoms)):
        neighbors = sorted([j for j in bond_dict.get(i, []) if j > i])
        if neighbors:
            neighbor_str = " ".join(f"{j+1} 1.0" for j in neighbors)
            conn_lines.append(f"{i+1} {neighbor_str}")
        else:
            conn_lines.append(f"{i+1}")

    return "\n".join(coords_lines), "\n".join(conn_lines)


def _write_gaussian_ff_params(parm) -> str:
    """
    Extract Amber-style force field parameters from parm7 for Gaussian (Amber=SoftOnly).

    This attempts to write a *self-contained* parameter section with the core Amber terms:
      - NonBon (Amber mixing rule + standard 1-4 scaling)
      - HrmStr1 (bonds)
      - HrmBnd1 (angles)
      - AmbTrs (proper torsions, periodicities 1-4)
      - ImpTrs (improper torsions)
      - VDW (per-atom-type LJ parameters; Radius = Rmin/2, Well-depth = epsilon)

    Limitations
    -----------
    - Amber torsions with periodicity > 4 are not representable by AmbTrs; they are skipped with a warning.
    - Per-dihedral 1-4 scaling factors (SCEE/SCNB) are not emitted (Gaussian uses the global NonBon scaling).
      This matches most standard Amber/GAFF workflows where SCEE/SCNB are uniform.
    """
    lines: List[str] = []

    # Non-bonded master function: Amber arithmetic mixing + standard exclusions and 1-4 scaling.
    # V-type=3 (Amber arithmetic), C-type=1 (Coulomb), cutoffs 0/0 (no explicit cutoffs here),
    # VScale: 1-2=0, 1-3=0, 1-4=0.5 ; CScale: 1-2=0, 1-3=0, 1-4=Amber default (1/1.2 via -1.2).
    lines.append("! Nonbonded master function (Amber defaults)")
    lines.append("NonBon 3 1 0 0 0.0 0.0 0.5 0.0 0.0 -1.2")

    # -------------------------
    # Bonds
    # -------------------------
    lines.append("")
    lines.append("! Bond parameters")
    bond_params: Set[Tuple[str, str, float, float]] = set()
    for bond in getattr(parm, "bonds", []):
        btype = getattr(bond, "type", None)
        if btype is None:
            continue
        try:
            k = float(getattr(btype, "k"))
            req = float(getattr(btype, "req"))
        except Exception:
            continue

        t1 = _fix_atom_type(getattr(bond.atom1, "atom_type", "X"))
        t2 = _fix_atom_type(getattr(bond.atom2, "atom_type", "X"))
        if t1 > t2:
            t1, t2 = t2, t1
        bond_params.add((t1, t2, k, req))

    for t1, t2, k, req in sorted(bond_params):
        lines.append(f"HrmStr1 {t1} {t2} {k:.6f} {req:.6f}")

    # -------------------------
    # Angles
    # -------------------------
    lines.append("")
    lines.append("! Angle parameters")
    angle_params: Set[Tuple[str, str, str, float, float]] = set()
    for angle in getattr(parm, "angles", []):
        atype = getattr(angle, "type", None)
        if atype is None:
            continue
        try:
            k = float(getattr(atype, "k"))
            theteq = float(getattr(atype, "theteq"))
        except Exception:
            continue

        t1 = _fix_atom_type(getattr(angle.atom1, "atom_type", "X"))
        t2 = _fix_atom_type(getattr(angle.atom2, "atom_type", "X"))
        t3 = _fix_atom_type(getattr(angle.atom3, "atom_type", "X"))

        # Sort endpoints for consistency
        if t1 > t3:
            t1, t3 = t3, t1
        angle_params.add((t1, t2, t3, k, theteq))

    for t1, t2, t3, k, theteq in sorted(angle_params):
        lines.append(f"HrmBnd1 {t1} {t2} {t3} {k:.6f} {theteq:.6f}")

    # -------------------------
    # Torsions (proper)
    # -------------------------
    def _as_term_list(dtype_obj: Any) -> List[Any]:
        if dtype_obj is None:
            return []
        # ParmEd uses DihedralTypeList for multi-term torsions (iterable)
        terms = getattr(dtype_obj, "terms", None)
        if terms is not None:
            try:
                return list(terms)
            except Exception:
                logger.debug("Failed to convert terms to list", exc_info=True)
        if isinstance(dtype_obj, (list, tuple)):
            return list(dtype_obj)
        # Try iteration (DihedralTypeList behaves like a list)
        try:
            if hasattr(dtype_obj, "__iter__") and not isinstance(dtype_obj, (str, bytes)):
                return list(dtype_obj)
        except Exception:
            logger.debug("Failed to iterate dtype_obj", exc_info=True)
        return [dtype_obj]

    def _get_attr(obj: Any, names: List[str], default: Any = None) -> Any:
        for n in names:
            if hasattr(obj, n):
                v = getattr(obj, n)
                if v is not None:
                    return v
        return default

    lines.append("")
    lines.append("! Proper torsions (AmbTrs)")
    # key -> (phase[4], mag[4])
    tors_params: Dict[Tuple[str, str, str, str], Tuple[List[float], List[float]]] = {}

    # Separate proper vs improper if ParmEd exposes `impropers`
    dihedrals_all = list(getattr(parm, "dihedrals", []) or [])
    impropers_from_attr = list(getattr(parm, "impropers", []) or [])
    if impropers_from_attr:
        proper_dihedrals = dihedrals_all
        improper_dihedrals = impropers_from_attr
    else:
        proper_dihedrals = [d for d in dihedrals_all if not bool(getattr(d, "improper", False))]
        improper_dihedrals = [d for d in dihedrals_all if bool(getattr(d, "improper", False))]

    for dih in proper_dihedrals:
        dtype = getattr(dih, "type", None)
        for term in _as_term_list(dtype):
            try:
                per = _get_attr(term, ["per", "periodicity", "period"], None)
                phase = float(_get_attr(term, ["phase", "phi", "phase_shift"], 0.0))
                mag = float(_get_attr(term, ["phi_k", "pk", "k", "barrier"], 0.0))
                div = float(_get_attr(term, ["div", "divider", "idivf", "npaths"], 1.0))
                if div == 0.0:
                    div = 1.0
            except Exception:
                continue

            try:
                n = int(round(abs(float(per))))
            except Exception:
                continue
            if n < 1:
                continue
            if n > 4:
                click.echo(
                    f"[oniom-export] WARNING: skipping Amber torsion with periodicity {n} (>4) "
                    f"for types {_fix_atom_type(dih.atom1.atom_type)}-{_fix_atom_type(dih.atom2.atom_type)}-"
                    f"{_fix_atom_type(dih.atom3.atom_type)}-{_fix_atom_type(dih.atom4.atom_type)}"
                )
                continue

            t1 = _fix_atom_type(getattr(dih.atom1, "atom_type", "X"))
            t2 = _fix_atom_type(getattr(dih.atom2, "atom_type", "X"))
            t3 = _fix_atom_type(getattr(dih.atom3, "atom_type", "X"))
            t4 = _fix_atom_type(getattr(dih.atom4, "atom_type", "X"))
            key = (t1, t2, t3, t4)

            if key not in tors_params:
                tors_params[key] = ([0.0, 0.0, 0.0, 0.0], [0.0, 0.0, 0.0, 0.0])

            phases, mags = tors_params[key]
            idx = n - 1

            # Amber divides each term by IDIVF; Gaussian's AmbTrs uses a single NPaths.
            # We fold the division into the magnitude and set NPaths=1.
            mag_eff = mag / div

            if mags[idx] != 0.0 and abs(phases[idx] - phase) > 1e-6:
                click.echo(
                    f"[oniom-export] WARNING: multiple torsion terms with the same periodicity {n} but "
                    f"different phases for key {key}. Keeping the first phase {phases[idx]:.3f} and "
                    f"adding magnitudes."
                )
            if mags[idx] == 0.0:
                phases[idx] = phase
            mags[idx] += mag_eff

    for (t1, t2, t3, t4) in sorted(tors_params.keys()):
        phases, mags = tors_params[(t1, t2, t3, t4)]
        if all(abs(m) < 1e-12 for m in mags):
            continue
        po1, po2, po3, po4 = phases
        m1, m2, m3, m4 = mags
        lines.append(
            f"AmbTrs {t1} {t2} {t3} {t4} "
            f"{po1:.6f} {po2:.6f} {po3:.6f} {po4:.6f} "
            f"{m1:.6f} {m2:.6f} {m3:.6f} {m4:.6f} 1"
        )

    # -------------------------
    # Improper torsions
    # -------------------------
    lines.append("")
    lines.append("! Improper torsions (ImpTrs)")
    improper_params: Set[Tuple[str, str, str, str, float, float, float]] = set()

    for imp in improper_dihedrals:
        dtype = getattr(imp, "type", None)
        for term in _as_term_list(dtype):
            try:
                per = float(abs(float(_get_attr(term, ["per", "periodicity", "period"], 2.0))))
                phase = float(_get_attr(term, ["phase", "phi", "phase_shift"], 0.0))
                mag = float(_get_attr(term, ["phi_k", "pk", "k", "barrier"], 0.0))
                div = float(_get_attr(term, ["div", "divider", "idivf", "npaths"], 1.0))
                if div == 0.0:
                    div = 1.0
            except Exception:
                continue

            t1 = _fix_atom_type(getattr(imp.atom1, "atom_type", "X"))
            t2 = _fix_atom_type(getattr(imp.atom2, "atom_type", "X"))
            t3 = _fix_atom_type(getattr(imp.atom3, "atom_type", "X"))
            t4 = _fix_atom_type(getattr(imp.atom4, "atom_type", "X"))

            mag_eff = mag / div
            improper_params.add((t1, t2, t3, t4, mag_eff, phase, per))

    for t1, t2, t3, t4, mag_eff, phase, per in sorted(improper_params):
        if abs(mag_eff) < 1e-12:
            continue
        lines.append(f"ImpTrs {t1} {t2} {t3} {t4} {mag_eff:.6f} {phase:.6f} {per:.6f}")

    # -------------------------
    # VDW parameters
    # -------------------------
    lines.append("")
    lines.append("! VDW parameters (Radius = Rmin/2 [Å], Well-depth = epsilon [kcal/mol])")
    vdw_params: Dict[str, Tuple[float, float]] = {}
    for atom in getattr(parm, "atoms", []):
        atype = _fix_atom_type(getattr(atom, "atom_type", "X"))
        if atype in vdw_params:
            continue

        # ParmEd often provides rmin (Amber Rmin/2), epsilon, and/or sigma
        radius: Optional[float] = None
        epsilon: Optional[float] = None

        try:
            if hasattr(atom, "rmin") and getattr(atom, "rmin") is not None:
                radius = float(getattr(atom, "rmin"))
            elif hasattr(atom, "rmin_half") and getattr(atom, "rmin_half") is not None:
                radius = float(getattr(atom, "rmin_half"))
            elif hasattr(atom, "sigma") and getattr(atom, "sigma") is not None:
                # Convert sigma -> Rmin/2 using rmin = 2^(1/6)*sigma, and radius = rmin/2
                radius = float(getattr(atom, "sigma")) * (2.0 ** (1.0 / 6.0)) / 2.0
        except Exception:
            radius = None

        try:
            if hasattr(atom, "epsilon") and getattr(atom, "epsilon") is not None:
                epsilon = float(getattr(atom, "epsilon"))
        except Exception:
            epsilon = None

        if radius is None or epsilon is None:
            continue

        vdw_params[atype] = (radius, epsilon)

    for atype, (radius, epsilon) in sorted(vdw_params.items()):
        lines.append(f"VDW {atype} {radius:.6f} {epsilon:.6f}")

    return "\n".join(lines)

def export_gaussian(
    parm7_path: Path,
    model_pdb: Optional[Path],
    output_path: Path,
    method: str = "wB97XD/def2-TZVPD",
    qm_charge: int = 0,
    qm_mult: int = 1,
    near_cutoff: float = 6.0,
    nproc: int = 8,
    mem: str = "16GB",
    qm_residues: Optional[List[int]] = None,
    input_path: Optional[Path] = None,
    element_check: bool = True,
) -> None:
    """
    Generate Gaussian ONIOM input file from parm7.

    Args:
        parm7_path: Path to Amber parm7 topology
        model_pdb: Path to PDB defining QM region atoms (optional)
        output_path: Output Gaussian input file
        method: QM method and basis set
        qm_charge: Charge of QM region
        qm_mult: Multiplicity of QM region
        near_cutoff: Distance cutoff for movable atoms (Angstrom)
        nproc: Number of processors
        mem: Memory allocation
        qm_residues: List of 0-based residue indices for QM region (alternative to model_pdb)
        input_path: Coordinate file (.pdb or .xyz). If omitted, uses coordinates stored in the ParmEd object.
        element_check: If True, validate element sequence between parm7 and input coordinates.
    """
    _check_parmed()

    parm = pmd.load_file(str(parm7_path))

    # Load / attach coordinates
    elements_for_output: Optional[List[str]] = None
    if input_path is not None:
        coords, input_elems = _read_input_geometry(input_path)
        _apply_coordinates_to_parm(parm, coords)
        elements_for_output = input_elems
        if element_check:
            _validate_element_order(parm, input_elems, strict=True)
    else:
        coords_attr = getattr(parm, "coordinates", None)
        n_coords = 0
        try:
            n_coords = len(coords_attr) if coords_attr is not None else 0
        except Exception:
            n_coords = 0
        if n_coords == 0:
            raise ValueError(
                "No coordinates found in the loaded parm7. "
                "Please provide a coordinate file with -i/--input (PDB or XYZ)."
            )
        elements_for_output = _get_parm_elements(parm)

    # Detect layer indices from B-factors if the input is a layered PDB produced by mlmm_toolkit.
    layer_info: Optional[Dict[str, List[int]]] = None
    if input_path is not None and input_path.suffix.lower() in {".pdb", ".ent"}:
        try:
            from .utils import (
                has_valid_layer_bfactors,
                parse_layer_indices_from_bfactors,
                read_bfactors_from_pdb,
            )

            bfactors = read_bfactors_from_pdb(input_path)
            if len(bfactors) == len(parm.atoms) and has_valid_layer_bfactors(bfactors):
                layer_info = parse_layer_indices_from_bfactors(bfactors)
                click.echo(
                    "[oniom-export] Detected ML/MM layer B-factors in the input PDB; "
                    "using them to decide movable/frozen atoms."
                )
        except Exception:
            layer_info = None

    # Determine QM region
    qm_indices: Set[int] = set()
    movable_indices: Set[int] = set()

    if model_pdb is not None:
        qm_indices = _read_qm_atoms_from_pdb(
            model_pdb,
            input_pdb=input_path
            if (input_path is not None and input_path.suffix.lower() in {".pdb", ".ent"})
            else None,
            system_coords=getattr(parm, "coordinates", None),
            system_elements=elements_for_output,
        )
    elif qm_residues:
        qm_indices, movable_indices = _identify_qm_atoms_by_distance(parm, qm_residues, near_cutoff)
    elif layer_info is not None and layer_info.get("ml_indices"):
        qm_indices = set(int(i) for i in layer_info["ml_indices"])
    else:
        raise ValueError(
            "No QM region specified. Provide --model-pdb, or supply a layered input PDB "
            "(B-factor=0 marks the ML/QM region), or use qm_residues in the Python API."
        )

    if not qm_indices:
        raise ValueError("No QM atoms identified")

    if max(qm_indices) >= len(parm.atoms):
        raise ValueError(
            f"QM index out of range: max(qm_indices)={max(qm_indices)} but topology has {len(parm.atoms)} atoms. "
            "Check that your model PDB / input PDB and parm7 have consistent atom ordering."
        )

    # Determine movable atoms for partial optimization.
    if layer_info is not None:
        frozen = set(int(i) for i in layer_info.get("frozen_indices", []))
        movable_indices = set(range(len(parm.atoms))) - frozen
    elif not movable_indices:
        # Distance-based selection: include all atoms in residues within `near_cutoff` of any QM atom.
        from scipy import spatial

        qm_list = sorted(qm_indices)
        neighbor_mask = np.any(
            spatial.distance.cdist(parm.coordinates, parm.coordinates[qm_list]) <= near_cutoff,
            axis=1,
        )
        neighbor_residues = set(parm.atoms[i].residue for i in np.where(neighbor_mask)[0])
        for residue in neighbor_residues:
            for atom in residue.atoms:
                movable_indices.add(atom.idx)

    movable_indices |= qm_indices  # QM atoms must always be movable

    # Detect covalent QM/MM boundaries and generate link-atom metadata.
    link_specs = _build_link_atom_specs(
        parm,
        qm_indices,
        elements=elements_for_output,
    )

    # Generate sections
    header = _write_gaussian_header(
        parm,
        str(parm7_path),
        str(output_path),
        method=method,
        nproc=nproc,
        mem=mem,
        qm_charge=qm_charge,
        qm_mult=qm_mult,
    )
    coords, connectivity = _write_gaussian_coordinates(
        parm,
        qm_indices,
        movable_indices,
        elements=elements_for_output,
        link_specs=link_specs,
    )
    ff_params = _write_gaussian_ff_params(parm)

    # Write output
    with output_path.open("w") as f:
        f.write(header)
        f.write(coords)
        f.write("\n\n")
        f.write(connectivity)
        f.write("\n\n")
        f.write(ff_params)
        f.write("\n\n")

    click.echo(f"[oniom-gaussian] Wrote '{output_path}'")
    click.echo(f"[oniom-gaussian] QM atoms: {len(qm_indices)}, Movable atoms: {len(movable_indices)}")
    click.echo(f"[oniom-gaussian] Link boundaries: {len(link_specs)}")


# -----------------------------------------------
# ORCA QM/MM Export
# -----------------------------------------------

# -----------------------------------------------
# ORCA QM/MM Export
# -----------------------------------------------

def _format_orca_index_set(indices: Set[int]) -> str:
    """
    Format a set of 0-based atom indices using ORCA's compact range syntax.

    Example:
        {0:3 7 10:12}
    """
    if not indices:
        return "{}"
    sorted_idx = sorted(int(i) for i in indices)
    parts: List[str] = []
    start = prev = sorted_idx[0]
    for i in sorted_idx[1:]:
        if i == prev + 1:
            prev = i
            continue
        parts.append(f"{start}:{prev}" if prev > start else f"{start}")
        start = prev = i
    parts.append(f"{start}:{prev}" if prev > start else f"{start}")
    return "{" + " ".join(parts) + "}"


def _manual_orcaff_command(parm7_path: Path, out_dir: Path) -> str:
    """Return a shell command users can run to generate ORCAFF.prms manually."""
    return (
        f"cd {shlex.quote(str(out_dir.resolve()))} && "
        f"orca_mm -convff -AMBER {shlex.quote(str(parm7_path.resolve()))}"
    )


def _resolve_oniom_mode(mode: Optional[str], output_path: Path) -> str:
    """Resolve export mode using explicit `--mode` first, then output suffix."""
    if mode is not None:
        return str(mode).strip().lower()

    suffix = output_path.suffix.lower()
    if suffix in {".gjf", ".com"}:
        return "g16"
    if suffix == ".inp":
        return "orca"

    raise ValueError(
        f"Could not infer export mode from -o/--output '{output_path}'. "
        "Specify --mode (g16/orca) or use an output suffix: .gjf/.com (g16), .inp (orca)."
    )


def export_orca(
    parm7_path: Path,
    model_pdb: Optional[Path],
    output_path: Path,
    method: str = "B3LYP D3BJ def2-SVP",
    qm_charge: int = 0,
    qm_mult: int = 1,
    total_charge: Optional[int] = None,
    total_mult: Optional[int] = None,
    nproc: int = 8,
    near_cutoff: float = 6.0,
    qm_residues: Optional[List[int]] = None,
    input_path: Optional[Path] = None,
    element_check: bool = True,
    orcaff_path: Optional[Path] = None,
    convert_orcaff: bool = True,
) -> None:
    """
    Generate an ORCA QM/MM input file.

    ORCA's QM/MM implementation expects an ORCA force-field parameter file (ORCAFF.prms).
    This can be generated from an Amber topology (prmtop/parm7) using ORCA's `orca_mm` utility:

        orca_mm -convff -AMBER <topology.prmtop>

    This exporter will try to generate the ORCAFF file automatically when:
      - `orcaff_path` is not provided, and
      - `convert_orcaff=True`, and
      - `orca_mm` is found in PATH.

    Args:
        parm7_path: Path to Amber parm7/prmtop topology file.
        model_pdb: PDB defining QM region (typically subset PDB).
        output_path: Output ORCA input file (.inp).
        method: ORCA QM method line (e.g., "B3LYP D3BJ def2-SVP").
        qm_charge: Charge of the QM region.
        qm_mult: Multiplicity of the QM region.
        total_charge: Charge of the full QM+MM system for Charge_Total in %qmmm.
            If None, uses topology total charge.
        total_mult: Multiplicity of the full QM+MM system for Mult_Total in %qmmm.
            If None, uses qm_mult.
        nproc: Number of processors.
        near_cutoff: Distance cutoff (Å) used to define ActiveAtoms when no layer B-factors exist.
        qm_residues: Alternative QM definition by 0-based residue indices in the ParmEd structure.
        input_path: Coordinate file (.pdb or .xyz). Atom order must match the topology.
        element_check: Validate element sequence between input and topology (best-effort).
        orcaff_path: Path to ORCAFF.prms file. If None, uses/creates <parm7_stem>.ORCAFF.prms in output dir.
        convert_orcaff: If True, try to run `orca_mm -convff -AMBER` when ORCAFF.prms is missing.
    """
    _check_parmed()

    parm = pmd.load_file(str(parm7_path))

    # Load / attach coordinates
    elements_for_output: Optional[List[str]] = None
    if input_path is not None:
        coords, input_elems = _read_input_geometry(input_path)
        _apply_coordinates_to_parm(parm, coords)
        elements_for_output = input_elems
        if element_check:
            _validate_element_order(parm, input_elems, strict=True)
    else:
        coords_attr = getattr(parm, "coordinates", None)
        n_coords = 0
        try:
            n_coords = len(coords_attr) if coords_attr is not None else 0
        except Exception:
            n_coords = 0
        if n_coords == 0:
            raise ValueError(
                "No coordinates found in the loaded topology/structure. "
                "Please provide a coordinate file with -i/--input (PDB or XYZ)."
            )
        elements_for_output = _get_parm_elements(parm)

    # Detect layer indices from B-factors if input is a layered PDB produced by mlmm_toolkit.
    layer_info: Optional[Dict[str, List[int]]] = None
    if input_path is not None and input_path.suffix.lower() in {".pdb", ".ent"}:
        try:
            from .utils import (
                has_valid_layer_bfactors,
                parse_layer_indices_from_bfactors,
                read_bfactors_from_pdb,
            )

            bfactors = read_bfactors_from_pdb(input_path)
            if len(bfactors) == len(parm.atoms) and has_valid_layer_bfactors(bfactors):
                layer_info = parse_layer_indices_from_bfactors(bfactors)
                click.echo(
                    "[oniom-export] Detected ML/MM layer B-factors in the input PDB; "
                    "using them to decide movable/frozen atoms."
                )
        except Exception:
            layer_info = None

    # Determine QM region
    qm_indices: Set[int] = set()
    movable_indices: Set[int] = set()

    if model_pdb is not None:
        qm_indices = _read_qm_atoms_from_pdb(
            model_pdb,
            input_pdb=input_path
            if (input_path is not None and input_path.suffix.lower() in {".pdb", ".ent"})
            else None,
            system_coords=getattr(parm, "coordinates", None),
            system_elements=elements_for_output,
        )
    elif qm_residues:
        qm_indices, movable_indices = _identify_qm_atoms_by_distance(parm, qm_residues, near_cutoff)
    elif layer_info is not None and layer_info.get("ml_indices"):
        qm_indices = set(int(i) for i in layer_info["ml_indices"])
    else:
        raise ValueError(
            "No QM region specified. Provide --model-pdb, or supply a layered input PDB "
            "(B-factor=0 marks the ML/QM region), or use qm_residues in the Python API."
        )

    if not qm_indices:
        raise ValueError("No QM atoms identified")

    if max(qm_indices) >= len(parm.atoms):
        raise ValueError(
            f"QM index out of range: max(qm_indices)={max(qm_indices)} but topology has {len(parm.atoms)} atoms. "
            "Check that your model PDB / input PDB and parm7 have consistent atom ordering."
        )

    # Determine ActiveAtoms (movable atoms)
    if layer_info is not None:
        frozen = set(int(i) for i in layer_info.get("frozen_indices", []))
        movable_indices = set(range(len(parm.atoms))) - frozen
    elif not movable_indices:
        # Distance-based selection: include all atoms in residues within `near_cutoff` of any QM atom.
        from scipy import spatial

        qm_list = sorted(qm_indices)
        neighbor_mask = np.any(
            spatial.distance.cdist(parm.coordinates, parm.coordinates[qm_list]) <= near_cutoff,
            axis=1,
        )
        neighbor_residues = set(parm.atoms[i].residue for i in np.where(neighbor_mask)[0])
        for residue in neighbor_residues:
            for atom in residue.atoms:
                movable_indices.add(atom.idx)

    movable_indices |= qm_indices

    if total_charge is None:
        total_charge = _get_total_charge(parm)
    if total_mult is None:
        total_mult = int(qm_mult)

    # ORCA generates link atoms automatically from QMAtoms and ORCAFF topology.
    # We still detect boundaries to report what was found and to keep behavior explicit.
    link_specs = _build_link_atom_specs(
        parm,
        qm_indices,
        elements=elements_for_output,
    )

    # Resolve/generate ORCAFF.prms
    out_dir = output_path.parent
    manual_orcaff_cmd = _manual_orcaff_command(parm7_path, out_dir)
    if orcaff_path is None:
        expected = out_dir / f"{parm7_path.stem}.ORCAFF.prms"
        orcaff_path = expected

        if not orcaff_path.exists() and convert_orcaff:
            orca_mm = shutil.which("orca_mm")
            if orca_mm is None:
                click.echo(
                    "[oniom-orca] WARNING: ORCAFF.prms not found and `orca_mm` is not available on PATH.\n"
                    f"[oniom-orca]          Run manually: {manual_orcaff_cmd}"
                )
            else:
                click.echo(
                    f"[oniom-orca] Generating ORCAFF.prms via: {manual_orcaff_cmd}"
                )
                proc = subprocess.run(
                    [orca_mm, "-convff", "-AMBER", str(parm7_path)],
                    cwd=str(out_dir),
                    stdout=subprocess.PIPE,
                    stderr=subprocess.STDOUT,
                    text=True,
                )
                if proc.returncode != 0:
                    raise RuntimeError(
                        "orca_mm failed "
                        f"(exit {proc.returncode}).\n"
                        f"Run manually: {manual_orcaff_cmd}\n"
                        f"Output:\n{proc.stdout}"
                    )

                # Try to locate the generated file (orca_mm typically writes <stem>.ORCAFF.prms)
                if not orcaff_path.exists():
                    candidates = sorted(out_dir.glob("*.ORCAFF.prms"), key=lambda p: p.stat().st_mtime, reverse=True)
                    if candidates:
                        orcaff_path = candidates[0]

    # ORCA input (use compact range syntax; indices are 0-based)
    qm_atoms_str = _format_orca_index_set(qm_indices)
    active_atoms_str = _format_orca_index_set(movable_indices)
    link_comment_block = ""
    if link_specs:
        link_lines = ["# Estimated link-H positions (QM/MM boundary caps; Angstrom)"]
        for mm_idx, spec in sorted(link_specs.items()):
            pos = np.asarray(spec["position"], dtype=float)
            link_lines.append(
                f"#   QM {int(spec['qm_idx']) + 1:>5d}  MM {int(mm_idx) + 1:>5d}  "
                f"Hcap ({pos[0]:10.6f}, {pos[1]:10.6f}, {pos[2]:10.6f})"
            )
        link_comment_block = "\n".join(link_lines) + "\n"

    # Prefer a relative filename when possible
    orcaff_ref = str(orcaff_path) if orcaff_path.is_absolute() else orcaff_path.name

    orca_input = f"""# ORCA QM/MM input generated by mlmm oniom-orca
# Amber topology: {parm7_path}
# ORCAFF parameters: {orcaff_ref}
# Coordinates: {input_path if input_path is not None else "(from topology/structure)"}
{link_comment_block}

%pal nprocs {nproc} end

! {method}
! QMMM

%qmmm
    ORCAFFFilename "{orcaff_ref}"
    QMAtoms {qm_atoms_str} end
    ActiveAtoms {active_atoms_str} end
    Charge_Total {int(total_charge)}
    Mult_Total {int(total_mult)}
end

* xyz {qm_charge} {qm_mult}
"""

    # Add coordinates
    elements_list: Optional[List[str]] = None
    if elements_for_output is not None and len(elements_for_output) == len(parm.atoms):
        elements_list = elements_for_output

    for atom in parm.atoms:
        x, y, z = atom.xx, atom.xy, atom.xz
        if elements_list is not None:
            element = elements_list[atom.idx]
        else:
            element = _get_parm_element(atom)
        orca_input += f"  {element:<2}  {x:12.6f} {y:12.6f} {z:12.6f}\n"

    orca_input += "*\n"

    # Write output
    with output_path.open("w") as f:
        f.write(orca_input)

    click.echo(f"[oniom-orca] Wrote '{output_path}'")
    click.echo(f"[oniom-orca] QM atoms: {len(qm_indices)}, Active atoms: {len(movable_indices)}")
    click.echo(f"[oniom-orca] Link boundaries (auto-capped by ORCA): {len(link_specs)}")
    if orcaff_path is not None:
        if orcaff_path.exists():
            click.echo(f"[oniom-orca] ORCAFF.prms: {orcaff_path}")
        else:
            click.echo(
                f"[oniom-orca] NOTE: ORCAFF.prms not found at '{orcaff_path}'. "
                f"Run manually: {manual_orcaff_cmd}"
            )



# -----------------------------------------------
# CLI Commands
# -----------------------------------------------

@click.command(
    name="oniom-export",
    help="Export ONIOM input from Amber parm7 topology (Gaussian g16 or ORCA).",
    context_settings={"help_option_names": ["-h", "--help"]},
)
@click.option(
    "--parm",
    "parm7",
    type=click.Path(path_type=Path, exists=True, dir_okay=False),
    required=True,
    help="Amber parm7 topology file.",
)
@click.option(
    "-i",
    "--input",
    "input_coords",
    type=click.Path(path_type=Path, exists=True, dir_okay=False),
    default=None,
    help="Coordinate file (.pdb or .xyz) for the current structure (atom order must match parm7).",
)
@click.option(
    "--element-check/--no-element-check",
    default=True,
    show_default=True,
    help="Validate that the element sequence in --input matches the parm7 topology.",
)
@click.option(
    "--model-pdb",
    type=click.Path(path_type=Path, exists=True, dir_okay=False),
    default=None,
    help="PDB file defining QM region atoms.",
)
@click.option(
    "-o",
    "--output",
    type=click.Path(path_type=Path, dir_okay=False),
    required=True,
    help="Output file path (.gjf/.com for g16, .inp for ORCA when --mode is omitted).",
)
@click.option(
    "--mode",
    type=click.Choice(["g16", "orca"], case_sensitive=False),
    default=None,
    help="Export mode. If omitted, inferred from -o suffix: .gjf/.com -> g16, .inp -> orca.",
)
@click.option(
    "--method",
    type=str,
    default=None,
    help="QM method and basis set. Defaults depend on mode.",
)
@click.option(
    "-q",
    "--charge",
    type=int,
    required=True,
    help="Charge of QM region.",
)
@click.option(
    "-m",
    "--multiplicity",
    type=int,
    default=1,
    show_default=True,
    help="Multiplicity of QM region.",
)
@click.option(
    "--near",
    type=float,
    default=6.0,
    show_default=True,
    help="Distance cutoff for movable/active atoms (Angstrom).",
)
@click.option(
    "--nproc",
    type=int,
    default=8,
    show_default=True,
    help="Number of processors.",
)
@click.option(
    "--mem",
    type=str,
    default="16GB",
    show_default=True,
    help="Memory allocation (g16 mode).",
)
@click.option(
    "--total-charge",
    type=int,
    default=None,
    help="Total charge of full QM+MM system for ORCA Charge_Total (orca mode).",
)
@click.option(
    "--total-mult",
    type=int,
    default=None,
    help="Total multiplicity of full QM+MM system for ORCA Mult_Total (orca mode).",
)
@click.option(
    "--orcaff",
    type=click.Path(exists=True, path_type=Path),
    default=None,
    help="Path to ORCAFF.prms (orca mode). If omitted, uses/creates <parm7_stem>.ORCAFF.prms in output directory.",
)
@click.option(
    "--convert-orcaff/--no-convert-orcaff",
    default=True,
    show_default=True,
    help="If ORCAFF.prms is missing, try `orca_mm -convff -AMBER` automatically (orca mode).",
)
def cli(
    parm7: Path,
    input_coords: Optional[Path],
    element_check: bool,
    model_pdb: Optional[Path],
    output: Path,
    mode: Optional[str],
    method: Optional[str],
    charge: int,
    multiplicity: int,
    near: float,
    nproc: int,
    mem: str,
    total_charge: Optional[int],
    total_mult: Optional[int],
    orcaff: Optional[Path],
    convert_orcaff: bool,
) -> None:
    """Export Gaussian/ORCA ONIOM input via a unified entrypoint."""
    try:
        resolved_mode = _resolve_oniom_mode(mode, output)

        if resolved_mode == "g16":
            export_gaussian(
                parm7_path=parm7,
                model_pdb=model_pdb,
                output_path=output,
                method=method or _GAUSSIAN_DEFAULT_METHOD,
                qm_charge=charge,
                qm_mult=multiplicity,
                near_cutoff=near,
                nproc=nproc,
                mem=mem,
                input_path=input_coords,
                element_check=element_check,
            )
            return

        export_orca(
            parm7_path=parm7,
            model_pdb=model_pdb,
            output_path=output,
            method=method or _ORCA_DEFAULT_METHOD,
            qm_charge=charge,
            qm_mult=multiplicity,
            total_charge=total_charge,
            total_mult=total_mult,
            nproc=nproc,
            near_cutoff=near,
            input_path=input_coords,
            element_check=element_check,
            orcaff_path=orcaff,
            convert_orcaff=convert_orcaff,
        )
    except Exception as e:
        click.echo(f"ERROR: {e}", err=True)
        raise SystemExit(1)
