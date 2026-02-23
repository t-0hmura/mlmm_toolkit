# mlmm_toolkit/mm_parm.py

"""
AmberTools prmtop/rst7 builder with automatic GAFF2 ligand parameterization.

Example:
    mlmm mm-parm -i input.pdb --out-prefix complex --ligand-charge "GPP=-3"

For detailed documentation, see: docs/mm_parm.md
"""

from __future__ import annotations

import os
import re
import shutil
import subprocess
import sys
import tempfile
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Set, Tuple

import click

# ===================== User dictionaries & constants =====================

AMINO_ACIDS: Dict[str, int] = {
    # --- Standard 20 (L) ---
    "ALA":  0, "ARG": +1, "ASN":  0, "ASP": -1, "CYS":  0,
    "GLU": -1, "GLN":  0, "GLY":  0, "HIS":  0, "ILE":  0,
    "LEU":  0, "LYS": +1, "MET":  0, "PHE":  0, "PRO":  0,
    "SER":  0, "THR":  0, "TRP":  0, "TYR":  0, "VAL":  0,

    # --- Canonical extras ---
    "SEC":  0,   # selenocysteine
    "PYL": +1,   # pyrrolysine

    # --- Protonation / tautomers (Amber/CHARMM style) ---
    "HIP": +1,   # fully protonated His
    "HID":  0,   # Nδ-protonated His
    "HIE":  0,   # Nε-protonated His
    "ASH":  0,   # neutral Asp
    "GLH":  0,   # neutral Glu
    "LYN":  0,   # neutral Lys
    "ARN":  0,   # neutral Arg
    "TYM": -1,   # deprotonated Tyr (phenolate)

    # --- Phosphorylated residues ---
    "SEP": -2, "TPO": -2, "PTR": -2,

    # --- Cys family ---
    "CYX":  0,   # disulfide Cys
    "CSO":  0,   # Cys sulfenic acid
    "CSD": -1,   # Cys sulfinic acid
    "CSX":  0,   # generic Cys derivative
    "OCS": -1,   # cysteic acid
    "CYM": -1,   # deprotonated Cys

    # --- Lys variants / carboxylation ---
    "MLY": +1, "LLP": +1, "DLY": +1,
    "KCX": -1,   # Lysine Nz-Carboxylic Acid

    # --- D isomers (19 residues) ---
    "DAL":  0, "DAR": +1, "DSG": 0, "DAS": -1, "DCY": 0,
    "DGN":  0, "DGL": -1, "DHI": 0, "DIL":  0, "DLE": 0,
    "DLY": +1, "MED":  0, "DPN": 0, "DPR":  0, "DSN": 0,
    "DTH":  0, "DTR":  0, "DTY": 0, "DVA":  0,

    # --- Carboxylation / cyclization / others ---
    "CGU": -2,   # gamma-carboxy-glutamate
    "CGA": -1,   # carboxymethylated glutamate
    "PCA":  0,   # pyroglutamate
    "MSE":  0,   # selenomethionine
    "OMT":  0,   # methionine sulfone

    # --- Other modified residues possibly encountered ---
    "ASA": 0, "CIR": 0, "FOR": 0, "MVA": 0, "IIL": 0, "AIB": 0, "HTN": 0,
    "SAR": 0, "NMC": 0, "PFF": 0, "NFA": 0, "ALY": 0, "AZF": 0, "CNX": 0, "CYF": 0,

    # --- Hydroxyproline ---
    "HYP": 0,

    # --- All C-terminus forms ---
    "CALA": -1, "CARG":  0, "CASN": -1, "CASP": -2, "CCYS": -1,
    "CCYX": -1, "CGLN": -1, "CGLU": -2, "CGLY": -1, "CHID": -1,
    "CHIE": -1, "CHIP":  0, "CHYP": -1, "CILE": -1, "CLEU": -1,
    "CLYS":  0, "CMET": -1, "CPHE": -1, "CPRO": -1, "CSER": -1,
    "CTHR": -1, "CTRP": -1, "CTYR": -1, "CVAL": -1, "NHE": 0,
    "NME": 0,
    "CTER": -1,  # generic C-terminus

    # --- All N-terminus forms ---
    "NALA": +1, "NARG": +2, "NASN": +1, "NASP":  0, "NCYS": +1,
    "NCYX": +1, "NGLN": +1, "NGLU":  0, "NGLY": +1, "NHID": +1,
    "NHIE": +1, "NHIP": +2, "NILE": +1, "NLEU": +1, "NLYS": +2,
    "NMET": +1, "NPHE": +1, "NPRO": +1, "NSER": +1, "NTHR": +1,
    "NTRP": +1, "NTYR": +1, "NVAL": +1, "ACE": 0,
    "NTER": +1,  # generic N-terminus
}

# Common ions (by residue name) and their formal charges
ION: Dict[str, int] = {
    # +1
    "LI": +1, "NA": +1, "K": +1, "RB": +1, "CS": +1, "TL": +1, "AG": +1, "CU1": +1,
    "Ag": +1, "K+": +1, "Na+": +1, "NH4": +1, "H3O+": +1, "HE+": +1, "HZ+": +1, "Tl": +1,

    # +2
    "MG": +2, "CA": +2, "SR": +2, "BA": +2, "MN": +2, "FE2": +2, "CO": +2, "NI": +2,
    "CU": +2, "ZN": +2, "CD": +2, "HG": +2, "PB": +2, "Be": +2, "PD": +2, "PT": +2,
    "Sn": +2, "Ra": +2, "YB2": +2, "V2+": +2,

    # +3
    "FE": +3, "AU3": +3, "AL": +3, "GA": +3, "IN": +3,
    "CE": +3, "Ce": +3, "CR": +3, "Cr": +3, "Dy": +3, "EU": +3, "EU3": +3, "Er": +3,
    "GD3": +3, "LA": +3, "LU": +3, "Nd": +3, "PR": +3, "SM": +3, "Sm": +3, "TB": +3,
    "Tm": +3, "Y": +3, "Pu": +3,

    # +4
    "U4+": +4, "Th": +4, "Hf": +4, "Zr": +4,

    # -1
    "F": -1, "CL": -1, "BR": -1, "I": -1, "Cl-": -1, "IOD": -1,
}

# Water residue names considered as "water"
WATER_RES = {"HOH", "WAT", "H2O", "DOD", "TIP", "TIP3", "SOL"}

# Distance cutoff (Å) for disulfide detection (SG–SG)
DISULFIDE_CUTOFF = 2.3  # Å

# Hint message printed when the build fails
HINT_MESSAGE = (
    "[HINT] When the build fails, please check:\n"
    "  - TER records are present between protein chains in the input PDB.\n"
    "  - Ligand formal charges and spin multiplicities (defaults: 0 and 1) are set correctly via --ligand-charge/--ligand-mult.\n"
    "  - Hydrogens have been correctly added to the ligand (e.g. with --add-h/--ph).\n"
)

# ff19SB/OPC3 + nucleic/lipid/GLYCAM + GAFF2
LEAPRC_LINES = [
    "source leaprc.protein.ff19SB",
    "source leaprc.phosaa19SB",
    "source leaprc.protein.ff19SB_modAA",
    "source leaprc.lipid21",
    "source leaprc.RNA.OL3",
    "source leaprc.DNA.OL21",
    "source leaprc.GLYCAM_06j-1",
    "source leaprc.water.opc3",
    "source leaprc.gaff2",
    "loadamberparams frcmod.ionslm_126_opc3",
]

# AmberTools leaprc set for ff14SB/TIP3P
LEAPRC_LINES_OLD = [
    "source leaprc.protein.ff14SB",
    "source leaprc.phosaa14SB",
    "source leaprc.protein.ff14SB_modAA",
    "source leaprc.lipid21",
    "source leaprc.RNA.OL3",
    "source leaprc.DNA.OL21",
    "source leaprc.GLYCAM_06j-1",
    "source leaprc.water.tip3p",
    "source leaprc.gaff2",
    "loadamberparams frcmod.ionsjc_tip3p",
    "loadamberparams frcmod.ions1lm_126_tip3p",
    "loadamberparams frcmod.ions234lm_126_tip3p",
]

# ===================== Utilities =====================


def which(cmd: str) -> Optional[str]:
    """Return the path if command exists in PATH; otherwise None."""
    return shutil.which(cmd)


def ambertools_available() -> bool:
    """Return True if tleap, antechamber and parmchk2 are available on PATH."""
    return which("tleap") is not None and which("antechamber") is not None and which("parmchk2") is not None


def run(cmd: List[str], cwd: Optional[Path] = None, logfile: Optional[Path] = None) -> int:
    """Run a subprocess, capture stdout+stderr into a log file, and return the return code."""
    with subprocess.Popen(
        cmd,
        cwd=str(cwd) if cwd else None,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True,
        bufsize=1,
    ) as p:
        lines: List[str] = []
        for line in p.stdout:  # type: ignore
            lines.append(line)
        rc = p.wait()
    if logfile:
        logfile.write_text("".join(lines), encoding="utf-8", errors="ignore")
    return rc


def parse_ligand_charge(expr: Optional[str]) -> Dict[str, int]:
    """
    Parse '--ligand-charge' string into a dict.
    Accepts 'GPP=-3,MMT=-1' OR 'GPP:-3,MMT:-1' → {'GPP': -3, 'MMT': -1}.
    """
    if not expr:
        return {}
    out: Dict[str, int] = {}
    for tok in expr.split(","):
        tok = tok.strip()
        if not tok:
            continue
        if "=" in tok:
            k, v = tok.split("=", 1)
        elif ":" in tok:
            k, v = tok.split(":", 1)
        else:
            raise ValueError(f"Invalid format in --ligand-charge: {tok} (use RES=Q or RES:Q)")
        out[k.strip()] = int(v.strip())
    return out


def parse_ligand_mult(expr: Optional[str]) -> Dict[str, int]:
    """
    Parse '--ligand-mult' string into a dict of spin multiplicities.
    Accepts 'HEM=1,NO=2' OR 'HEM:1,NO:2' → {'HEM': 1, 'NO': 2}.
    """
    if not expr:
        return {}
    out: Dict[str, int] = {}
    for tok in expr.split(","):
        tok = tok.strip()
        if not tok:
            continue
        if "=" in tok:
            k, v = tok.split("=", 1)
        elif ":" in tok:
            k, v = tok.split(":", 1)
        else:
            raise ValueError(f"Invalid format in --ligand-mult: {tok} (use RES=M or RES:M)")
        out[k.strip()] = int(v.strip())
    return out


def copy_pdb_no_fix(pdb_path: Path, tmpdir: Path) -> Path:
    """Copy the input PDB verbatim to tmpdir/fixed.pdb (no structural fixing)."""
    fixed_pdb = tmpdir / "fixed.pdb"
    shutil.copy2(pdb_path, fixed_pdb)
    return fixed_pdb


def add_hydrogens_with_pdbfixer(pdb_in: Path, pdb_out: Path, ph: float) -> None:
    """
    Add hydrogens at the specified pH using PDBFixer, without adding missing heavy atoms/residues.
    """
    try:
        from pdbfixer import PDBFixer
        from pdbfixer import pdbfixer as _pdbfixer_mod
    except Exception as e:
        raise RuntimeError(
            "PDBFixer is required to use --add-h True, but it was not found."
        ) from e

    pdbfile_writer = getattr(getattr(_pdbfixer_mod, "app", None), "PDBFile", None)
    if pdbfile_writer is None:
        raise RuntimeError(
            "PDBFixer installation is incomplete: could not access PDB writer."
        )

    fixer = PDBFixer(filename=str(pdb_in))
    fixer.addMissingHydrogens(pH=ph)  # only Hs
    with open(pdb_out, "w") as f:
        pdbfile_writer.writeFile(fixer.topology, fixer.positions, f, keepIds=True)


def detect_disulfides_from_pdb(
    pdb_path: Path,
    cutoff: float = DISULFIDE_CUTOFF,
) -> List[Tuple[Tuple[str, int], Tuple[str, int]]]:
    """
    Extract SG (or S) atoms from CYS/CYM/CYX in a PDB and return residue-pairs
    with SG–SG distance ≤ cutoff Å.
    Return format: [((chainID, resSeq), (chainID, resSeq)), ...]
    """
    sg_sites: List[Tuple[str, int, float, float, float]] = []
    with open(pdb_path, "r") as f:
        for line in f:
            if not (line.startswith("ATOM") or line.startswith("HETATM")):
                continue
            resname = line[17:20].strip()
            # Restrict disulfide detection to residues defined in AMINO_ACIDS.
            if resname not in AMINO_ACIDS:
                continue
            if resname not in {"CYS", "CYM", "CYX"}:
                continue
            atom_name = line[12:16].strip()
            if atom_name not in {"SG", "S"}:
                continue
            # altLoc: only blank or 'A'
            altloc = line[16].strip()
            if altloc not in ("", "A"):
                continue
            chain = line[21]
            resseq_field = line[22:26]
            try:
                resseq = int(resseq_field)
            except Exception:
                try:
                    resseq = int(resseq_field.strip())
                except Exception:
                    continue
            try:
                x = float(line[30:38])
                y = float(line[38:46])
                z = float(line[46:54])
            except Exception:
                continue
            sg_sites.append((chain, resseq, x, y, z))
    pairs: List[Tuple[Tuple[str, int], Tuple[str, int]]] = []
    for i in range(len(sg_sites)):
        ci, ri, xi, yi, zi = sg_sites[i]
        for j in range(i + 1, len(sg_sites)):
            cj, rj, xj, yj, zj = sg_sites[j]
            dx = xi - xj
            dy = yi - yj
            dz = zi - zj
            dist = (dx * dx + dy * dy + dz * dz) ** 0.5
            if dist <= cutoff:
                pairs.append(((ci, ri), (cj, rj)))
    return pairs


def build_leap_residue_index(pdb_path: Path) -> Dict[Tuple[str, str], int]:
    """
    Build a mapping (chainID, resSeq as 4-char string) → LEaP 1-based residue index
    by scanning the PDB in order. LEaP numbers residues by appearance order, which can
    differ from RESSEQ integers; this avoids mismatches when bonding.
    """
    mapping: Dict[Tuple[str, str], int] = {}
    seen: Set[Tuple[str, str]] = set()
    idx = 0
    with open(pdb_path, "r") as f:
        for line in f:
            if not (line.startswith("ATOM") or line.startswith("HETATM")):
                continue
            chain = line[21]
            resseq = line[22:26]  # 4-character, right-justified
            key = (chain, resseq)
            if key not in seen:
                idx += 1
                seen.add(key)
                mapping[key] = idx
    return mapping


def parse_tleap_unknown_residues(leap_log: Path) -> Set[str]:
    """Parse a LEaP log and collect residue names reported as unknown/failed."""
    txt = leap_log.read_text(encoding="utf-8", errors="ignore")
    res: Set[str] = set()
    patterns = [
        r"Unknown residue:\s+([A-Za-z0-9\+\-]+)",
        r"Could not find in database the residue:\s+([A-Za-z0-9\+\-]+)",
        r"createAtomUnit:.*\bresidue\s+([A-Za-z0-9\+\-]+)\b",
        r"Creating new UNIT for residue:\s+([A-Za-z0-9_+\-]+)",
    ]
    for pat in patterns:
        for m in re.finditer(pat, txt):
            rn = m.group(1).strip()
            if rn in WATER_RES:
                continue
            res.add(rn)
    return res


# -------- helpers for TER insertion --------


def insert_ter_around_special_residues(pdb_in: Path, pdb_out: Path, special_resnames: Set[str]) -> None:
    """
    Insert TER records before/after contiguous blocks of residues whose names are in
    `special_resnames` (e.g., ligand names from --ligand-charge, WATER_RES, and ION).
    If such residues are consecutive, do not insert TER between them. Existing TER
    records are preserved and duplicate consecutive TERs are avoided.
    """

    def recname(line: str) -> str:
        return line[:6].strip()

    out_lines: List[str] = []
    prev_key: Optional[Tuple[str, str]] = None  # (chain, resseq[22:26])
    prev_special: Optional[bool] = None
    last_written_was_TER = False

    with open(pdb_in, "r") as f:
        for line in f:
            rn = recname(line)
            if rn in {"ATOM", "HETATM"}:
                chain = line[21]
                resseq = line[22:26]
                resname = line[17:20].strip()
                key = (chain, resseq)
                curr_special = resname in special_resnames

                if prev_key is None:
                    pass
                elif key != prev_key:
                    # boundary (previous residue → current residue)
                    if prev_special != curr_special and (prev_special or curr_special):
                        if not last_written_was_TER:
                            out_lines.append("TER\n")
                            last_written_was_TER = True
                out_lines.append(line)
                last_written_was_TER = False
                prev_key = key
                prev_special = curr_special

            elif rn == "TER":
                if not last_written_was_TER:
                    out_lines.append(line)
                last_written_was_TER = True
                prev_key = None
                prev_special = None
            else:
                out_lines.append(line)

    # Add trailing TER if the file ends with a "special" residue and no TER has been written
    if prev_key is not None and prev_special and not last_written_was_TER:
        out_lines.append("TER\n")

    with open(pdb_out, "w") as w:
        w.writelines(out_lines)


# -------- helper for amino-acid–like residue detection --------


def is_protein_like_residue_in_pdb(residue_pdb: Path) -> bool:
    """Return True if the PDB contains N, CA, and C atom names (amino-acid–like residue)."""
    need = {"N", "CA", "C"}
    present: Set[str] = set()
    with open(residue_pdb, "r") as f:
        for line in f:
            if not (line.startswith("ATOM") or line.startswith("HETATM")):
                continue
            name = line[12:16].strip()
            if name in need:
                present.add(name)
    return need.issubset(present)


# ===================== AmberTools route =====================


def extract_first_residue_pdb(src_pdb: Path, resname: str, dst_pdb: Path) -> bool:
    """
    Extract only the *first occurrence* of the specified residue name from src_pdb,
    write it as a standalone PDB to dst_pdb, and return True on success.
    """
    found = False
    out_lines: List[str] = []
    with open(src_pdb, "r") as f:
        for line in f:
            if not (line.startswith("ATOM") or line.startswith("HETATM")):
                continue
            rn = line[17:20].strip()
            if rn == resname and not found:
                resseq = line[22:26]
                chain = line[21]
                found = True
                break
    if not found:
        return False
    with open(src_pdb, "r") as f:
        for line in f:
            if not (line.startswith("ATOM") or line.startswith("HETATM")):
                continue
            rn = line[17:20].strip()
            if rn == resname and line[22:26] == resseq and line[21] == chain:
                out_lines.append(line)
    if not out_lines:
        return False
    with open(dst_pdb, "w") as w:
        for ln in out_lines:
            w.write(ln)
        w.write("END\n")
    return True


def antechamber_parametrize(resname: str, res_charge: int, res_mult: int, workdir: Path) -> Tuple[Path, Path]:
    """
    Run antechamber (GAFF2 + AM1-BCC) and parmchk2 to generate mol2/frcmod.
    Input: {resname}.pdb → Output: {resname}.mol2, {resname}.frcmod (paths returned).
    """
    pdb = workdir / f"{resname}.pdb"
    mol2 = workdir / f"{resname}.mol2"
    frcmod = workdir / f"{resname}.frcmod"
    cmd1 = [
        "antechamber",
        "-i",
        pdb.name,
        "-fi",
        "pdb",
        "-o",
        mol2.name,
        "-fo",
        "mol2",
        "-at",
        "gaff2",
        "-c",
        "bcc",
        "-nc",
        str(res_charge),
        "-m",
        str(res_mult),
        "-rn",
        resname,
        "-s",
        "2",
    ]
    rc1 = run(cmd1, cwd=workdir, logfile=workdir / f"{resname}.antechamber.log")
    if rc1 != 0 or not mol2.exists():
        raise RuntimeError(f"[{resname}] antechamber failed (see log).")
    cmd2 = ["parmchk2", "-i", mol2.name, "-f", "mol2", "-o", frcmod.name, "-s", "2"]
    rc2 = run(cmd2, cwd=workdir, logfile=workdir / f"{resname}.parmchk2.log")
    if rc2 != 0 or not frcmod.exists():
        raise RuntimeError(f"[{resname}] parmchk2 failed (see log).")
    click.echo(f"[mm-parm] Built GAFF2 parameter for [{resname}] (charge={res_charge}, mult={res_mult}).")
    return mol2, frcmod


def write_tleap_input(
    fixed_pdb: Path,
    lig_defs: List[Tuple[str, Path, Path]],
    ss_pairs: List[Tuple[Tuple[str, int], Tuple[str, int]]],
    out_prefix: str,
    tleap_in: Path,
    leaprc_lines: List[str],
) -> None:
    """
    Compose a LEaP input script.
    - lig_defs: list of (RESNAME, lib_or_mol2_path, frcmod_path)
        * .lib  → loadoff + loadamberparams frcmod
        * .mol2 → RES = loadmol2 + loadamberparams frcmod
    - ss_pairs: ((chainID, resSeq), (chainID, resSeq)) residue pairs to bond (S–S)
      (LEaP residue indices are resolved from PDB order via an internal mapping).
    """
    lines: List[str] = []
    lines += leaprc_lines

    # ligands / nonstandard residues
    for resname, param_file, frcmod in lig_defs:
        if param_file.suffix.lower() == ".lib":
            lines.append(f"loadoff {param_file.name}")
            lines.append(f"loadamberparams {frcmod.name}")
        else:
            lines.append(f"{resname} = loadmol2 {param_file.name}")
            lines.append(f"loadamberparams {frcmod.name}")

    # complex
    lines.append(f"complex = loadpdb {fixed_pdb.name}")

    # S–S bonds
    resnum_map = build_leap_residue_index(fixed_pdb)
    for (c1, r1), (c2, r2) in ss_pairs:
        key1 = (c1, f"{r1:>4}")
        key2 = (c2, f"{r2:>4}")
        if key1 in resnum_map and key2 in resnum_map:
            n1, n2 = resnum_map[key1], resnum_map[key2]
            lines.append(f"bond complex.{n1}.SG complex.{n2}.SG")
        else:
            lines.append(f"# WARN: could not resolve SS pair ({c1}{r1})-({c2}{r2})")

    # For logging: print charge in tleap output
    lines.append("charge complex")

    # outputs (parm7/inpcrd + pdb)
    lines.append(f"saveamberparm complex {out_prefix}.parm7 {out_prefix}.inpcrd")
    lines.append(f"savepdb       complex {out_prefix}.pdb")
    lines.append("quit")
    tleap_in.write_text("\n".join(lines) + "\n", encoding="utf-8")


def ambertools_route(
    pdb: Path,
    out_prefix: str,
    ligand_charge: Dict[str, int],
    ligand_mult: Dict[str, int],
    allow_nonstandard_aa: bool,
    keep_temp: bool,
    tmpdir: Path,
    ff_set: str,
    add_ter: bool,
) -> Tuple[Path, Path]:
    """
    AmberTools route:
      - (Optionally) add hydrogens beforehand (done in run_pipeline).
      - Use the input PDB as-is (copied to fixed.pdb). Optionally insert TERs.
      - Detect candidate S–S bonds by SG–SG geometry.
      - First, run LEaP without extra parameters. If unknown residues are reported,
        parameterize them with antechamber+parmchk2 (GAFF2/AM1‑BCC).
      - Nonstandard amino-acid–like residues (containing N/CA/C) are not treated
        automatically; the build aborts with an explanatory message.
      - Load parameters and re-run LEaP. LEaP writes complex.parm7/complex.inpcrd/
        complex.pdb; this function copies complex.parm7 and complex.inpcrd to
        <out_prefix>.parm7 and <out_prefix>.rst7. The caller handles the PDB export.
    """
    leaprc_lines = LEAPRC_LINES if ff_set == "ff19SB" else LEAPRC_LINES_OLD
    protein_ff = "ff19SB" if ff_set == "ff19SB" else "ff14SB"
    _ = protein_ff  # currently unused, kept for clarity/extension

    # PDB as-is, optional TER insertion
    fixed_pdb = copy_pdb_no_fix(pdb, tmpdir)
    if add_ter:
        special_resnames: Set[str] = set(ligand_charge.keys()) | set(WATER_RES) | set(ION.keys())
        fixed_pdb_with_ter = tmpdir / "fixed_withTER.pdb"
        insert_ter_around_special_residues(fixed_pdb, fixed_pdb_with_ter, special_resnames)
        fixed_pdb = fixed_pdb_with_ter

    # Detect S–S candidates
    ss_pairs = detect_disulfides_from_pdb(fixed_pdb, cutoff=DISULFIDE_CUTOFF)

    # Pass 1 (no extra params) -> will write complex.parm7/.inpcrd/.pdb
    leap_in = tmpdir / "tleap_1.in"
    write_tleap_input(
        fixed_pdb,
        lig_defs=[],
        ss_pairs=ss_pairs,
        out_prefix="complex",
        tleap_in=leap_in,
        leaprc_lines=leaprc_lines,
    )
    log1 = tmpdir / "tleap_1.log"
    run(["tleap", "-f", leap_in.name], cwd=tmpdir, logfile=log1)

    # Collect unknown residues
    need_params: Set[str] = parse_tleap_unknown_residues(log1)

    # Parameterize unknown residues
    lig_defs: List[Tuple[str, Path, Path]] = []
    for rn in sorted(need_params):
        charge = ligand_charge.get(rn, AMINO_ACIDS.get(rn, 0))
        mult = ligand_mult.get(rn, 1)
        lig_pdb = tmpdir / f"{rn}.pdb"
        ok = extract_first_residue_pdb(fixed_pdb, rn, lig_pdb)
        if not ok:
            raise RuntimeError(f"Failed to extract PDB for unknown residue {rn}")

        # Detect nonstandard amino-acid–like residues and abort with a message.
        if is_protein_like_residue_in_pdb(lig_pdb) and not allow_nonstandard_aa:
            raise RuntimeError(
                f"Nonstandard amino acid residue '{rn}' detected. "
                "mm_parm does not support nonstandard amino acids by default. "
                "Please parameterize this residue manually with AmberTools, "
                "edit the residue to a standard amino acid, or rerun with "
                "--allow-nonstandard-aa to force antechamber parameterization."
            )

        mol2, frcmod = antechamber_parametrize(rn, charge, mult, tmpdir)
        lig_defs.append((rn, mol2, frcmod))

    # Pass 2 (with generated parameters) -> will (re)write complex.* including PDB
    if need_params:
        leap_in2 = tmpdir / "tleap_2.in"
        write_tleap_input(
            fixed_pdb,
            lig_defs=lig_defs,
            ss_pairs=ss_pairs,
            out_prefix="complex",
            tleap_in=leap_in2,
            leaprc_lines=leaprc_lines,
        )
        log2 = tmpdir / "tleap_2.log"
        run(["tleap", "-f", leap_in2.name], cwd=tmpdir, logfile=log2)
        if not (tmpdir / "complex.parm7").exists():
            raise RuntimeError(f"tleap failed to produce parm7; see {log2.name}.")

    # Copy outputs (parm7, inpcrd) to final names
    src_parm = tmpdir / "complex.parm7"
    src_inp = tmpdir / "complex.inpcrd"
    if not (src_parm.exists() and src_inp.exists()):
        msg = f"LEaP outputs not found in {tmpdir}. Check logs: {tmpdir / 'tleap_1.log'}"
        if (tmpdir / "tleap_2.log").exists():
            msg += f" and {tmpdir / 'tleap_2.log'}"
        raise FileNotFoundError(msg)

    parm7 = Path(f"{out_prefix}.parm7").resolve()
    rst7 = Path(f"{out_prefix}.rst7").resolve()
    shutil.copy2(src_parm, parm7)
    shutil.copy2(src_inp, rst7)  # copy LEaP ASCII inpcrd as <prefix>.rst7

    # Return paths for prmtop/rst7; the caller will copy PDB using naming rule
    return parm7, rst7


# ===================== Main pipeline (library/CLI entry) =====================


@dataclass
class Args:
    pdb: Path
    out_prefix: str
    ligand_charge: Dict[str, int]
    ligand_mult: Dict[str, int]
    allow_nonstandard_aa: bool
    keep_temp: bool
    add_ter: bool
    add_h: bool
    ph: float
    ff_set: str  # "ff19SB" or "ff14SB"
    out_prefix_given: bool  # whether user explicitly provided --out-prefix


def run_pipeline(args: Args) -> None:
    if not args.pdb.exists():
        sys.exit(f"PDB not found: {args.pdb}")

    if not ambertools_available():
        sys.exit("AmberTools (tleap, antechamber, parmchk2) not found on PATH; this script requires AmberTools.")

    # Decide PDB filename to export/copy (used both on success and as H-added fallback)
    # When --out-prefix is omitted and --add-h False, do not write <input_stem>_parm.pdb.
    final_pdb_out: Optional[Path]
    if args.out_prefix_given:
        final_pdb_out = Path(f"{args.out_prefix}.pdb").resolve()
    else:
        if args.add_h:
            final_pdb_out = Path(f"{Path(args.pdb).stem}_parm.pdb").resolve()
        else:
            final_pdb_out = None

    # Prepare temporary working directory
    tmp_mgr: Optional[tempfile.TemporaryDirectory] = None
    if args.keep_temp:
        tmpdir_path = Path(tempfile.mkdtemp(prefix="parm7build_", dir=os.getcwd()))
    else:
        tmp_mgr = tempfile.TemporaryDirectory(prefix="parm7build_")
        tmpdir_path = Path(tmp_mgr.name)

    fixed_pdb_with_H: Optional[Path] = None  # for fallback export

    try:
        # Copy input PDB locally (avoid path/lock issues)
        local_pdb = tmpdir_path / "input.pdb"
        shutil.copy2(args.pdb, local_pdb)

        # Optional: add hydrogens via PDBFixer at specified pH
        prepared_pdb = local_pdb
        if args.add_h:
            fixed_pdb = tmpdir_path / "input_withH.pdb"
            click.echo(f"[mm-parm] Adding hydrogens with PDBFixer at pH={args.ph:.2f} ...")
            try:
                add_hydrogens_with_pdbfixer(local_pdb, fixed_pdb, args.ph)
            except Exception as e:
                if args.keep_temp:
                    click.echo(
                        f"[mm-parm] ERROR: PDBFixer hydrogen addition failed: {e}\n"
                        f"Temporary working directory kept at: {tmpdir_path}",
                        err=True,
                    )
                raise
            prepared_pdb = fixed_pdb
            fixed_pdb_with_H = fixed_pdb
            click.echo("[mm-parm] Hydrogens added (PDBFixer).")

        try:
            click.echo("[mm-parm] AmberTools detected. Using tleap + GAFF2 (AM1-BCC).")
            click.echo(
                f"[mm-parm] FF set: {args.ff_set} | add_ter: {args.add_ter} | "
                f"add_h: {args.add_h} (pH={args.ph:.2f})"
            )
            parm7, rst7 = ambertools_route(
                prepared_pdb,
                args.out_prefix,
                args.ligand_charge,
                args.ligand_mult,
                args.allow_nonstandard_aa,
                args.keep_temp,
                tmpdir_path,
                ff_set=args.ff_set,
                add_ter=args.add_ter,
            )
        except Exception as e:
            # Fallback export of H-added PDB on failure
            if fixed_pdb_with_H is not None and fixed_pdb_with_H.exists() and final_pdb_out is not None:
                try:
                    shutil.copy2(fixed_pdb_with_H, final_pdb_out)
                    click.echo(f"[mm-parm] Build failed, but wrote hydrogen-added PDB fallback: {final_pdb_out}")
                except Exception as copy_e:
                    click.echo(f"[mm-parm] WARNING: Failed to write fallback hydrogen-added PDB: {copy_e}", err=True)
            if args.keep_temp:
                click.echo(f"[mm-parm] ERROR: Failed: {e}\nTemporary working directory kept at: {tmpdir_path}", err=True)
            # Re-raise to preserve error behavior
            raise

        # Copy LEaP PDB (complex.pdb) to final name, if requested
        if final_pdb_out is not None:
            src_pdb = tmpdir_path / "complex.pdb"
            if src_pdb.exists():
                shutil.copy2(src_pdb, final_pdb_out)
                click.echo(f"[mm-parm] Wrote: {final_pdb_out}")
            else:
                click.echo("[mm-parm] WARNING: LEaP PDB (complex.pdb) was not found; skipping PDB export copy.", err=True)

        click.echo(f"[mm-parm] Wrote: {parm7}")
        click.echo(f"[mm-parm] Wrote: {rst7}")

        if args.keep_temp:
            click.echo(f"[mm-parm] Temporary directory kept: {tmpdir_path}")
            info = f"[mm-parm] LEaP logs: {tmpdir_path / 'tleap_1.log'}"
            if (tmpdir_path / "tleap_2.log").exists():
                info += f", {tmpdir_path / 'tleap_2.log'}"
            click.echo(info)
    except Exception:
        # Print a generic hint message on failure, then re-raise
        click.echo(HINT_MESSAGE, err=True)
        raise
    finally:
        if tmp_mgr is not None:
            try:
                tmp_mgr.cleanup()
            except Exception:
                pass


# ===================== Click CLI entry point =====================


@click.command(
    context_settings={"help_option_names": ["-h", "--help"]},
    help="Generate Amber parm7/rst7 (and a LEaP-exported PDB) from a PDB using AmberTools only.",
)
@click.option(
    "-i",
    "--input",
    "pdb",
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
    required=True,
    help="Input PDB file (used as-is; optional hydrogens via --add-h/--ph).",
)
@click.option(
    "--out-prefix",
    default=None,
    help=(
        "Output prefix (default: input PDB stem). For LEaP PDB: "
        "if omitted and --add-h True, <input_stem>_parm.pdb is used."
    ),
)
@click.option(
    "--ligand-charge",
    default=None,
    help=(
        'Comma-separated mapping of residue=charge or residue:charge '
        '(e.g., "GPP=-3,MMT=-1" or "GPP:-3,MMT:-1")'
    ),
)
@click.option(
    "--ligand-mult",
    default=None,
    help=(
        'Comma-separated mapping of residue=multiplicity or residue:multiplicity '
        '(e.g., "HEM=1,NO:2")'
    ),
)
@click.option(
    "--allow-nonstandard-aa",
    is_flag=True,
    help=(
        "Allow antechamber parameterization of residues that look amino-acid-like "
        "(contain N/CA/C). Use with care for modified amino acids."
    ),
)
@click.option(
    "--keep-temp",
    is_flag=True,
    help="Keep temporary working directory (in current dir) for debugging.",
)
@click.option(
    "--add-ter/--no-add-ter",
    "add_ter",
    default=True,
    show_default=True,
    help=(
        "Insert TER before/after target residues. "
        "When contiguous, TER is not inserted between them."
    ),
)
@click.option(
    "--add-h/--no-add-h",
    "add_h",
    default=False,
    show_default=True,
    help="Add hydrogens using PDBFixer at the specified --ph.",
)
@click.option(
    "--ph",
    "ph",
    type=float,
    default=7.0,
    help="pH used by PDBFixer when adding hydrogens (--add-h True). Default: 7.0",
)
@click.option(
    "--ff-set",
    type=click.Choice(["ff19SB", "ff14SB"]),
    default="ff19SB",
    help="Force-field set for proteins/backbone typing and water/ion parameters (default: ff19SB).",
)
def cli(
    pdb: Path,
    out_prefix: Optional[str],
    ligand_charge: Optional[str],
    ligand_mult: Optional[str],
    allow_nonstandard_aa: bool,
    keep_temp: bool,
    add_ter: bool,
    add_h: bool,
    ph: float,
    ff_set: str,
) -> None:
    """Click entry point that mirrors the documented CLI."""
    args = Args(
        pdb=pdb,
        out_prefix=out_prefix if out_prefix is not None else Path(pdb).stem,
        ligand_charge=parse_ligand_charge(ligand_charge),
        ligand_mult=parse_ligand_mult(ligand_mult),
        allow_nonstandard_aa=allow_nonstandard_aa,
        keep_temp=keep_temp,
        add_ter=bool(add_ter),
        add_h=bool(add_h),
        ph=ph,
        ff_set=ff_set,
        out_prefix_given=(out_prefix is not None),
    )
    run_pipeline(args)
