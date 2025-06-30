"""mlmm/utils.py — Lightweight CLI utilities bundled with ML/MM calculator.

Three commands are installed via ``setup.py``::

    xyz_geom2pdb   – Convert XYZ geometry / trajectory to PDB.
    add_elem_info  – Add element symbols (cols 77–78) to a PDB using Biopython.
    get_freeze_indices – List atom indices to freeze based on distance from ML region.

Each command can also be accessed as a Python API (see the function names
below).
"""
from __future__ import annotations

import argparse
import os
from collections import defaultdict
from pathlib import Path
from typing import List, Set

import numpy as np
from ase.io import read, write
from Bio.PDB import PDBParser, PDBIO

# ---------------------------------------------------------------------
# xyz_geom2pdb                                                               
# ---------------------------------------------------------------------

def convert_xyz_to_pdb(xyz_path: Path, ref_pdb_path: Path, out_pdb_path: Path) -> None:
    """Overlay coordinates from *xyz_path* onto *ref_pdb_path* topology and write *out_pdb_path*.

    *xyz_path* may contain one or many frames.  For multi‑frame trajectories a
    MODEL/ENDMDL block is appended for each subsequent frame.
    """
    ref_atoms = read(ref_pdb_path)  # Reference topology (single frame)
    traj = read(xyz_path, index=":", format="xyz")  # All XYZ frames
    if not traj:
        raise ValueError(f"No frames found in {xyz_path}.")

    for step, frame in enumerate(traj):
        atoms = ref_atoms.copy()
        atoms.set_positions(frame.get_positions())
        if step == 0:
            write(out_pdb_path, atoms)  # overwrite / create
        else:
            write(out_pdb_path, atoms, append=True)


def _xyz_argparser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="xyz_geom2pdb",
        description="Convert an XYZ geometry/trajectory to PDB using a reference PDB topology.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument("-i", "--input", dest='input', required=True, metavar="FILE.xyz", help="Input XYZ file")
    p.add_argument("-r", "--ref", dest='ref', required=True, metavar="REF.pdb", help="Reference PDB file with topology")
    p.add_argument("-o", "--output", dest='output', metavar="OUT.pdb", help="Output PDB (default: INPUT with .pdb)")
    return p


def xyz_geom2pdb_console() -> None:
    """Entry point exposed via *xyz_geom2pdb* console script."""
    args = _xyz_argparser().parse_args()
    xyz_path = Path(args.input).expanduser().resolve()
    ref_path = Path(args.ref).expanduser().resolve()
    out_path = Path(args.output).expanduser().resolve() if args.output else xyz_path.with_suffix(".pdb")

    if not xyz_path.is_file():
        raise FileNotFoundError(xyz_path)
    if not ref_path.is_file():
        raise FileNotFoundError(ref_path)

    convert_xyz_to_pdb(xyz_path, ref_path, out_path)
    print(f"[xyz_geom2pdb] Wrote {out_path}")

# ---------------------------------------------------------------------
# add_elem_info                                                              
# ---------------------------------------------------------------------

def _infer_element(atom_name: str) -> str:
    """Infer element symbol from atom name (e.g. 'C', 'CA', 'O2')."""
    name = "".join(c for c in atom_name if c.isalpha())
    if not name:
        return "  "
    if name[0].isdigit():
        name = name[1:]
    elem = name[0].upper()
    if len(name) > 1 and name[1].islower():
        elem += name[1].lower()
    return elem.rjust(2)


def add_elem_info(input_file: str, output_file: str | None = None) -> None:
    """Appends element symbols (cols 77–78) to a PDB via Biopython."""
    if output_file is None:
        output_file = input_file

    parser = PDBParser(QUIET=True, PERMISSIVE=True)
    structure = parser.get_structure("model", input_file)

    for atom in structure.get_atoms():
        atom.element = _infer_element(atom.get_name()).strip()
        
    io = PDBIO()
    io.set_structure(structure)
    io.save(output_file)


def _elem_argparser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="add_elem_info",
        description="Add element symbols (PDB cols 77–78) using Biopython.",
    )
    p.add_argument("input", metavar="INPUT.pdb", help="Input PDB file")
    p.add_argument("output", nargs="?", default=None, help="Output file (in‑place if omitted)")
    return p


def add_elem_info_console() -> None:
    args = _elem_argparser().parse_args()
    add_elem_info(args.input, args.output)
    target = args.output if args.output else args.input
    print(f"[add_elem_info] Element symbols added → {target}")

# ---------------------------------------------------------------------
# get_freeze_indices                                                         
# ---------------------------------------------------------------------

def get_freeze_indices(real_pdb: str, model_pdb: str, relaxed_range: float) -> List[int]:
    """Return 0‑based ASE atom indices farther than *relaxed_range* Å from ML region."""
    # 1) Build ID set for model PDB
    model_ids: Set[str] = set()
    with open(model_pdb, "r", encoding="utf-8") as f:
        for line in f:
            if line.startswith(("ATOM", "HETATM")):
                atom_name = line[12:16].strip()
                res_name = line[17:20].strip()
                res_seq = line[22:26].strip()
                model_ids.add(f"{atom_name} {res_name} {res_seq}")

    # 2) Parse real PDB into atom records & residue groupings
    leap_atoms: List[dict] = []
    residues: dict[tuple[str, str], List[dict]] = defaultdict(list)
    with open(real_pdb, "r", encoding="utf-8") as f:
        for line in f:
            if not line.startswith(("ATOM", "HETATM")):
                continue
            idx = int(line[6:11])
            atom_name = line[12:16].strip()
            res_name = line[17:20].strip()
            chain_id = line[21]
            res_seq = line[22:26].strip()
            coord = np.array([
                float(line[30:38]),
                float(line[38:46]),
                float(line[46:54]),
            ])
            atom_id = f"{atom_name} {res_name} {res_seq}"
            atom = {"idx": idx, "id": atom_id, "coord": coord}
            leap_atoms.append(atom)
            residues[(chain_id, res_seq)].append(atom)

    # 3) Indices of ML‑region atoms in real PDB
    ml_idx_set = {a["idx"] for a in leap_atoms if a["id"] in model_ids}
    if not ml_idx_set:
        raise ValueError("No matching atom IDs between model_pdb and real_pdb.")
    ml_coords = np.array([a["coord"] for a in leap_atoms if a["idx"] in ml_idx_set])

    # 4) Residue‑wise min‑distance test
    freeze_indices: List[int] = []
    for atom_list in residues.values():
        res_coords = np.stack([a["coord"] for a in atom_list])
        dists = np.linalg.norm(ml_coords[:, None, :] - res_coords[None, :, :], axis=-1)
        if dists.min() > relaxed_range:
            freeze_indices.extend(a["idx"] - 1 for a in atom_list)  # 0‑based
    return freeze_indices


def _freeze_argparser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="get_freeze_indices",
        description="Identify atoms (> range Å from ML region) to freeze in ASE.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument("--real", required=True, metavar="REAL.pdb", help="Full complex PDB (e.g. leap.pdb)")
    p.add_argument("--model", required=True, metavar="MODEL.pdb", help="ML‑region PDB")
    p.add_argument("--range", type=float, required=True, metavar="Å", help="Relaxed range in Å")
    p.add_argument("--out", default=None, help="Output text file (indices)")
    return p


def get_freeze_indices_console() -> None:
    """Entry point exposed via *get_freeze_indices* console script."""
    args = _freeze_argparser().parse_args()
    indices = get_freeze_indices(args.real, args.model, args.range)
    if args.out:
        Path(args.out).write_text("[" + ",".join(map(str, indices)) + "]" + "\n", encoding="utf-8")
        print(f"[get_freeze_indices] {len(indices)} atoms frozen → {args.out}")
    else:
        print("[" + ",".join(map(str, indices)) + "]" + "\n")

# ---------------------------------------------------------------------
# Fallback dispatch when executed directly                                   
# ---------------------------------------------------------------------

def _dispatch() -> None:  # pragma: no cover
    if len(os.sys.argv) <= 1:
        print(__doc__)
        raise SystemExit(0)
    subcmd = Path(os.sys.argv[1]).name.lower()
    os.sys.argv.pop(1)  # remove subcmd so parser sees its own args
    if subcmd == "xyz_geom2pdb":
        xyz_geom2pdb_console()
    elif subcmd == "add_elem_info":
        add_elem_info_console()
    elif subcmd == "get_freeze_indices":
        get_freeze_indices_console()
    else:
        print(f"Unknown subcommand: {subcmd}")
        raise SystemExit(1)


if __name__ == "__main__":
    _dispatch()
