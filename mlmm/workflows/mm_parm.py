"""
AmberTools prmtop/rst7 builder with automatic GAFF2 ligand parameterization.

Example:
    mlmm mm-parm -i input.pdb --out-prefix complex -l "GPP=-3"

For detailed documentation, see: docs/mm_parm.md
"""

from __future__ import annotations

import logging
import os
import re
import shutil
import subprocess
import sys
import tempfile
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional, Set, Tuple

import click

from mlmm.core.result_commit import commit_payloads

logger = logging.getLogger(__name__)

# ===================== User dictionaries & constants =====================

from mlmm.core.residue_data import (  # single source of truth (was a drifted local copy)
    AMINO_ACIDS,
    ION,
    WATER_RES,
    DISULFIDE_CUTOFF,
)

# Distance cutoff (Å) for retaining a peptide C-N connection when adding TERs.
PEPTIDE_BOND_CUTOFF = 1.9

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


# DO NOT INLINE: HPC submission often `qsub /.../envs/mlmm/bin/mlmm all ...` without conda activate; bare shutil.which("tleap") fails because conda env bin not on $PATH. The sys.executable + sys.prefix dual fallback recovers automatically.
def which(cmd: str) -> Optional[str]:
    """Return the path if *cmd* is available; otherwise None.

    Searches ``$PATH`` first, then the directory of the running interpreter
    and ``<sys.prefix>/bin`` (and ``<sys.prefix>/Scripts`` on Windows). This
    makes the AmberTools preflight succeed when ``mlmm`` is invoked by an
    absolute path from a conda env (e.g. ``/.../envs/mlmm/bin/mlmm``) WITHOUT
    ``conda activate`` — tleap/antechamber/parmchk2 then sit next to the
    interpreter but are not on ``$PATH``.
    """
    found = shutil.which(cmd)
    if found:
        return found
    extra_dirs = [
        os.path.dirname(os.path.realpath(sys.executable)),
        os.path.join(sys.prefix, "bin"),
        os.path.join(sys.prefix, "Scripts"),
    ]
    seen = set()
    for d in extra_dirs:
        if not d or d in seen:
            continue
        seen.add(d)
        cand = shutil.which(cmd, path=d)
        if cand:
            return cand
    return None


_AMBERTOOLS_REQUIRED_COMMANDS: Tuple[str, ...] = ("tleap", "antechamber", "parmchk2")


def ambertools_command_paths() -> Dict[str, Optional[str]]:
    """Return resolved paths for required AmberTools executables."""
    return {cmd: which(cmd) for cmd in _AMBERTOOLS_REQUIRED_COMMANDS}


def missing_ambertools_commands(paths: Optional[Dict[str, Optional[str]]] = None) -> List[str]:
    """Return required AmberTools commands that are missing from PATH."""
    resolved = paths if paths is not None else ambertools_command_paths()
    return [cmd for cmd in _AMBERTOOLS_REQUIRED_COMMANDS if not resolved.get(cmd)]


def copy_pdb_with_element_fields(source: Path, destination: Path) -> Tuple[int, int]:
    """Copy a LEaP PDB while filling missing element columns in place.

    LEaP commonly leaves columns 77--78 blank. Preserve every record, line
    ending, and atom ordering; short atom records are padded only as needed to
    fill those columns. The exported PDB therefore remains topology-matched to
    the generated ``parm7``.

    Returns ``(assigned, unresolved)``.
    """
    from mlmm.domain.add_elem_info import guess_element

    assigned = 0
    unresolved = 0
    output_lines: List[str] = []
    with source.open(encoding="utf-8", errors="replace", newline="") as handle:
        for raw_line in handle:
            if raw_line.endswith("\r\n"):
                newline = "\r\n"
            elif raw_line.endswith("\n"):
                newline = "\n"
            elif raw_line.endswith("\r"):
                newline = "\r"
            else:
                newline = ""
            line = raw_line.rstrip("\r\n")
            if line.startswith(("ATOM  ", "HETATM")):
                padded = line.ljust(78)
                if not padded[76:78].strip():
                    element = guess_element(
                        padded[12:16].strip(),
                        padded[17:20].strip(),
                        padded.startswith("HETATM"),
                    )
                    if element:
                        line = padded[:76] + f"{element:>2}" + padded[78:]
                        assigned += 1
                    else:
                        line = padded
                        unresolved += 1
            output_lines.append(line + newline)

    destination = destination.resolve()
    destination.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile(
        mode="w", encoding="utf-8", newline="", dir=destination.parent,
        delete=False,
    ) as handle:
        temporary = Path(handle.name)
        handle.writelines(output_lines)
    try:
        shutil.copymode(source, temporary)
        os.replace(temporary, destination)
    finally:
        if temporary.exists():
            temporary.unlink()
    return assigned, unresolved


def ambertools_available() -> bool:
    """Return True if tleap, antechamber and parmchk2 are available on PATH."""
    return not missing_ambertools_commands()


# DO NOT INLINE: tleap/antechamber emit large stdout; subprocess.run(capture_output=True) buffers everything (memory) and only shows at exit (debugging pain). Line-by-line Popen gives real-time progress + bounded memory.
def run(cmd: List[str], cwd: Optional[Path] = None, logfile: Optional[Path] = None) -> int:
    """Run a subprocess, capture stdout+stderr into a log file, and return the return code."""
    if not cmd:
        raise ValueError("run() requires a non-empty command")
    executable = which(cmd[0])
    if executable is None:
        raise FileNotFoundError(f"Required executable '{cmd[0]}' was not found.")
    resolved_cmd = [executable, *cmd[1:]]
    with subprocess.Popen(
        resolved_cmd,
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


_TLEAP_COMPLEX_OUTPUTS = ("complex.parm7", "complex.inpcrd", "complex.pdb")
_TLEAP_REQUIRED_OUTPUTS = ("complex.parm7", "complex.inpcrd")


def _invalidate_tleap_complex_outputs(tmpdir: Path) -> None:
    """Remove the preceding LEaP pass before starting a new generation."""
    for name in _TLEAP_COMPLEX_OUTPUTS:
        (Path(tmpdir) / name).unlink(missing_ok=True)


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
            raise click.BadParameter(
                f"invalid --ligand-charge token '{tok}': use RES=Q or RES:Q "
                f"(e.g. 'SAM:1,GPP:-3')."
            )
        try:
            out[k.strip()] = int(v.strip())
        except ValueError:
            raise click.BadParameter(
                f"--ligand-charge value for '{k.strip()}' must be an integer, got '{v.strip()}'."
            )
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
            raise click.BadParameter(
                f"invalid --ligand-mult token '{tok}': use RES=M or RES:M "
                f"(e.g. 'HEM:1,NO:2')."
            )
        try:
            multiplicity = int(v.strip())
        except ValueError:
            raise click.BadParameter(
                f"--ligand-mult value for '{k.strip()}' must be an integer, "
                f"got '{v.strip()}'."
            )
        if multiplicity < 1:
            raise click.BadParameter(
                f"--ligand-mult value for '{k.strip()}' must be >= 1, "
                f"got {multiplicity}."
            )
        out[k.strip()] = multiplicity
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
            "PDBFixer is required to use --add-h, but it was not found."
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
    cyx_only: bool = False,
) -> List[Tuple[Tuple[str, int], Tuple[str, int]]]:
    """
    Extract SG (or S) atoms from CYS/CYM/CYX in a PDB and return residue-pairs
    with SG–SG distance ≤ cutoff Å.

    With ``cyx_only`` the scan is restricted to residues already named CYX, so a
    disulfide is formed only where the input says so explicitly and a CYS that
    merely happens to sit close to another one is left alone.

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
            if resname not in ({"CYX"} if cyx_only else {"CYS", "CYM", "CYX"}):
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


def rename_disulfide_cys_to_cyx(
    pdb_in: Path,
    pdb_out: Path,
    ss_pairs: List[Tuple[Tuple[str, int], Tuple[str, int]]],
) -> int:
    """
    Rename CYS -> CYX for every residue taking part in a detected disulfide.

    tleap's ``bond`` only adds the S-S connection; it does not strip the CYS
    template's HG. A disulfide cysteine left named CYS therefore keeps HG and
    its SG becomes hypervalent (CB + HG + SG). CYX is the Amber template for a
    disulfide-bonded cysteine (no HG), so the rename must happen before
    ``loadpdb``. CYM (thiolate, formal -1) is left untouched because renaming
    it would silently change the net charge; it is reported instead.

    Returns the number of renamed residues.
    """
    targets = {res for pair in ss_pairs for res in pair}
    renamed_res: Set[Tuple[str, int]] = set()
    cym_left: List[Tuple[str, int]] = []
    with open(pdb_in, "r") as fi, open(pdb_out, "w") as fo:
        for line in fi:
            if line.startswith(("ATOM", "HETATM")) and len(line) >= 26:
                resname = line[17:20].strip()
                chain = line[21]
                try:
                    resseq = int(line[22:26])
                except ValueError:
                    fo.write(line)
                    continue
                if (chain, resseq) in targets:
                    if resname == "CYS":
                        line = line[:17] + "CYX" + line[20:]
                        renamed_res.add((chain, resseq))
                    elif resname == "CYM" and (chain, resseq) not in cym_left:
                        cym_left.append((chain, resseq))
            fo.write(line)
    for chain, resseq in cym_left:
        print(
            f"[mm-parm] WARNING: {chain}{resseq} is CYM but takes part in a detected "
            "disulfide; leaving it as CYM (renaming would change the net charge). "
            "Rename it to CYX in the input if the disulfide is intended."
        )
    return len(renamed_res)


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
    If such residues are consecutive, do not insert TER between them. Also insert
    TERs between amino-acid residues that are adjacent in the file but not connected
    by a peptide C-N bond. Existing TER records are preserved and duplicate
    consecutive TERs are avoided.
    """

    def recname(line: str) -> str:
        return line[:6].strip()

    def residue_key(line: str) -> Tuple[str, str, str]:
        return (line[21], line[22:26], line[26:27])

    def atom_name(line: str) -> str:
        return line[12:16].strip()

    def xyz(line: str) -> Optional[Tuple[float, float, float]]:
        try:
            return (float(line[30:38]), float(line[38:46]), float(line[46:54]))
        except ValueError:
            return None

    def dist(a: Tuple[float, float, float], b: Tuple[float, float, float]) -> float:
        dx = a[0] - b[0]
        dy = a[1] - b[1]
        dz = a[2] - b[2]
        return (dx * dx + dy * dy + dz * dz) ** 0.5

    def residue_info(lines: List[str]) -> Tuple[str, Tuple[str, str, str], Dict[str, Tuple[float, float, float]]]:
        first = lines[0]
        coords: Dict[str, Tuple[float, float, float]] = {}
        for ln in lines:
            pos = xyz(ln)
            if pos is not None:
                coords[atom_name(ln)] = pos
        return first[17:20].strip(), residue_key(first), coords

    def peptide_break(prev_lines: List[str], curr_lines: List[str]) -> bool:
        prev_res, prev_key, prev_atoms = residue_info(prev_lines)
        curr_res, curr_key, curr_atoms = residue_info(curr_lines)
        if prev_res not in AMINO_ACIDS or curr_res not in AMINO_ACIDS:
            return False
        if prev_key[0] != curr_key[0]:
            return True
        prev_c = prev_atoms.get("C")
        curr_n = curr_atoms.get("N")
        if prev_c is None or curr_n is None:
            return False
        return dist(prev_c, curr_n) > PEPTIDE_BOND_CUTOFF

    def write_ter_if_needed() -> None:
        nonlocal last_written_was_TER
        if not last_written_was_TER:
            out_lines.append("TER\n")
            last_written_was_TER = True

    out_lines: List[str] = []
    prev_residue_lines: Optional[List[str]] = None
    current_residue_lines: Optional[List[str]] = None
    current_key: Optional[Tuple[str, str, str]] = None
    last_written_was_TER = False

    def flush_current() -> None:
        nonlocal prev_residue_lines, current_residue_lines, current_key, last_written_was_TER
        if current_residue_lines is None:
            return
        if prev_residue_lines is not None:
            prev_resname = prev_residue_lines[0][17:20].strip()
            curr_resname = current_residue_lines[0][17:20].strip()
            prev_special = prev_resname in special_resnames
            curr_special = curr_resname in special_resnames
            if (prev_special != curr_special and (prev_special or curr_special)) or peptide_break(
                prev_residue_lines,
                current_residue_lines,
            ):
                write_ter_if_needed()
        out_lines.extend(current_residue_lines)
        last_written_was_TER = False
        prev_residue_lines = current_residue_lines
        current_residue_lines = None
        current_key = None

    with open(pdb_in, "r") as f:
        for raw in f:
            line = raw if raw.endswith("\n") else raw + "\n"
            rn = recname(line)
            if rn in {"ATOM", "HETATM"}:
                key = residue_key(line)
                if current_key is None or key != current_key:
                    flush_current()
                    current_key = key
                    current_residue_lines = []
                current_residue_lines.append(line)

            elif rn == "TER":
                flush_current()
                write_ter_if_needed()
                prev_residue_lines = None
            else:
                flush_current()
                out_lines.append(line)

    flush_current()
    if prev_residue_lines is not None:
        prev_resname = prev_residue_lines[0][17:20].strip()
        if prev_resname in special_resnames:
            write_ter_if_needed()

    with open(pdb_out, "w") as w:
        w.writelines(out_lines)


# -------- helper for amino-acid–like residue detection --------


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
        for raw in f:
            if not (raw.startswith("ATOM") or raw.startswith("HETATM")):
                continue
            line = raw if raw.endswith("\n") else raw + "\n"
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

    # Electron-count preflight (C-003): catch an odd-electron ligand BEFORE
    # antechamber/sqm, which otherwise fails opaquely with the sub-log gone.
    # Most often this is a mis-specified ligand charge or protonation state
    # (e.g. SAM: 22 H = neutral / 23 H = +1).
    try:
        from mlmm.core.utils import validate_charge_spin
        from mlmm.domain.add_elem_info import guess_element

        _elems = []
        with open(pdb, encoding="utf-8", errors="ignore") as _fh:
            for _ln in _fh:
                if not _ln.startswith(("ATOM", "HETATM")):
                    continue
                _e = _ln[76:78].strip()
                if not _e:
                    _atname = _ln[12:16].strip()
                    _rn = _ln[17:20].strip()
                    _e = guess_element(_atname, _rn, _ln.startswith("HETATM"))
                if _e:
                    _elems.append(_e)
    except Exception:
        _elems = []
    if _elems:
        try:
            validate_charge_spin(_elems, res_charge, res_mult, source=str(pdb))
        except ValueError as _vc:
            raise RuntimeError(
                f"[{resname}] electron-count check failed before antechamber: "
                f"{_vc} Verify the charge/protonation of [{resname}] "
                f"(e.g. via --ligand-charge/--ligand-mult). For SAM, 22 H = "
                f"neutral (charge 0) and 23 H = +1 cation; the input PDB's "
                f"actual protonation must match the specified charge."
            ) from _vc

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
    keep_temp: bool,
    tmpdir: Path,
    ff_set: str,
    add_ter: bool,
    auto_disulfide: bool = True,
) -> Tuple[Path, Path]:
    """
    AmberTools route:
      - (Optionally) add hydrogens beforehand (done in run_pipeline).
      - Use the input PDB as-is (copied to fixed.pdb). Optionally insert TERs.
      - Detect candidate S–S bonds by SG–SG geometry.
      - First, run LEaP without extra parameters. If unknown residues are reported,
        parameterize them with antechamber+parmchk2 (GAFF2/AM1‑BCC).
      - Residues listed in ``AMINO_ACIDS`` that remain unknown to LEaP are not treated
        automatically; the build aborts with an explanatory message.
      - Load parameters and re-run LEaP. LEaP writes complex.parm7/complex.inpcrd/
        complex.pdb; this function copies complex.parm7 and complex.inpcrd to
        <out_prefix>.parm7 and <out_prefix>.rst7. The caller handles the PDB export.
    """
    leaprc_lines = LEAPRC_LINES if ff_set == "ff19SB" else LEAPRC_LINES_OLD

    # PDB as-is, optional TER insertion
    fixed_pdb = copy_pdb_no_fix(pdb, tmpdir)
    if add_ter:
        special_resnames: Set[str] = set(ligand_charge.keys()) | set(WATER_RES) | set(ION.keys())
        fixed_pdb_with_ter = tmpdir / "fixed_withTER.pdb"
        insert_ter_around_special_residues(fixed_pdb, fixed_pdb_with_ter, special_resnames)
        fixed_pdb = fixed_pdb_with_ter

    # Detect S–S candidates. With --no-auto-disulfide only residues already named
    # CYX are considered, so a CYS is never bonded on proximity alone.
    ss_pairs = detect_disulfides_from_pdb(
        fixed_pdb, cutoff=DISULFIDE_CUTOFF, cyx_only=not auto_disulfide
    )

    # A disulfide cysteine must be CYX *before* loadpdb: tleap's `bond` adds the
    # S-S connection but does not strip the CYS template's HG, which would leave
    # SG hypervalent (CB + HG + SG). Only reachable with auto-detection on; with
    # --no-auto-disulfide every pair is already CYX, so nothing is renamed.
    if auto_disulfide and ss_pairs:
        cyx_pdb = tmpdir / "fixed_cyx.pdb"
        n_cyx = rename_disulfide_cys_to_cyx(fixed_pdb, cyx_pdb, ss_pairs)
        if n_cyx:
            print(f"[mm-parm] Renamed {n_cyx} disulfide CYS -> CYX before tleap.")
        fixed_pdb = cyx_pdb

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
    rc1 = run(["tleap", "-f", leap_in.name], cwd=tmpdir, logfile=log1)

    # Collect unknown residues
    need_params: Set[str] = parse_tleap_unknown_residues(log1)
    if rc1 != 0 and not need_params:
        raise RuntimeError(f"tleap pass 1 exited with code {rc1}; see {log1.name}.")

    # Parameterize unknown residues
    lig_defs: List[Tuple[str, Path, Path]] = []
    for rn in sorted(need_params):
        charge = ligand_charge.get(rn, AMINO_ACIDS.get(rn, 0))
        mult = ligand_mult.get(rn, 1)
        lig_pdb = tmpdir / f"{rn}.pdb"
        ok = extract_first_residue_pdb(fixed_pdb, rn, lig_pdb)
        if not ok:
            raise RuntimeError(f"Failed to extract PDB for unknown residue {rn}")

        # Explicit ligand-charge mappings take highest priority and force GAFF2 parameterization.
        if rn in ligand_charge:
            mol2, frcmod = antechamber_parametrize(rn, charge, mult, tmpdir)
            lig_defs.append((rn, mol2, frcmod))
            continue

        # Amino-acid residues must be handled by the selected Amber protein force field.
        if rn in AMINO_ACIDS:
            raise RuntimeError(
                f"Nonstandard amino acid residue '{rn}' is not supported by mm_parm. "
                "This workflow does not auto-parameterize amino-acid residues. Options: "
                "(1) explicitly list the residue in --ligand-charge to force GAFF2 "
                "parameterization; (2) edit the input structure to use a supported residue; or "
                "(3) for a modified residue / covalent cross-link (which tleap cannot build "
                "automatically), prepare the topology yourself with tleap (lib/frcmod + a `bond` "
                "command for the cross-link) and supply it to the compute subcommands via --parm. "
                "If that parm uses 4-point water (OPC/TIP4P), run with --mm-backend openmm, or "
                "convert the water to 3-point OPC3 so mlmm-toolkit's default MM backend handles it."
            )

        mol2, frcmod = antechamber_parametrize(rn, charge, mult, tmpdir)
        lig_defs.append((rn, mol2, frcmod))

    # Pass 2 (with generated parameters) -> will (re)write complex.* including PDB
    if need_params:
        # Pass 2 is a distinct output generation.  Remove every pass-1
        # candidate before dispatch so a failed LEaP process cannot satisfy
        # the final existence checks with the preceding generation.
        _invalidate_tleap_complex_outputs(tmpdir)
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
        rc2 = run(["tleap", "-f", leap_in2.name], cwd=tmpdir, logfile=log2)
        if rc2 != 0:
            raise RuntimeError(f"tleap pass 2 exited with code {rc2}; see {log2.name}.")
        if not all((tmpdir / name).exists() for name in _TLEAP_REQUIRED_OUTPUTS):
            raise RuntimeError(
                f"tleap pass 2 did not produce a complete complex generation; "
                f"see {log2.name}."
            )

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
    commit_payloads(
        parm7,
        {
            parm7: src_parm.read_bytes(),
            rst7: src_inp.read_bytes(),
        },
    )

    # Return paths for prmtop/rst7; the caller will copy PDB using naming rule
    return parm7, rst7


# ===================== Main pipeline (library/CLI entry) =====================


@dataclass
class Args:
    pdb: Path
    out_prefix: str
    ligand_charge: Dict[str, int]
    ligand_mult: Dict[str, int]
    keep_temp: bool
    add_ter: bool
    auto_disulfide: bool
    add_h: bool
    ph: float
    ff_set: str  # "ff19SB" or "ff14SB"
    out_prefix_given: bool  # whether user explicitly provided --out-prefix


def run_pipeline(args: Args) -> None:
    if not args.pdb.exists():
        sys.exit(f"PDB not found: {args.pdb}")

    amber_paths = ambertools_command_paths()
    missing_cmds = missing_ambertools_commands(amber_paths)
    if missing_cmds:
        found_lines = [
            f"  {name}: {amber_paths[name]}"
            for name in _AMBERTOOLS_REQUIRED_COMMANDS
            if amber_paths.get(name)
        ]
        missing_text = ", ".join(missing_cmds)
        details = "\n".join(found_lines) if found_lines else "  (none found)"
        sys.exit(
            "AmberTools preflight failed.\n"
            f"Missing required command(s): {missing_text}\n"
            "Required: tleap, antechamber, parmchk2\n"
            "Detected command paths:\n"
            f"{details}"
        )

    # Decide PDB filename to export/copy (used both on success and as H-added fallback)
    # Without --out-prefix or --add-h, do not write <input_stem>_parm.pdb.
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
                args.keep_temp,
                tmpdir_path,
                ff_set=args.ff_set,
                add_ter=args.add_ter,
                auto_disulfide=args.auto_disulfide,
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
                assigned, unresolved = copy_pdb_with_element_fields(src_pdb, final_pdb_out)
                click.echo(f"[mm-parm] Wrote: {final_pdb_out}")
                if assigned:
                    click.echo(f"[mm-parm] Populated element columns for {assigned} atoms.")
                if unresolved:
                    click.echo(
                        f"[mm-parm] WARNING: Could not infer element columns for {unresolved} atoms.",
                        err=True,
                    )
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
                logger.debug("Failed to clean up temporary directory", exc_info=True)


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
    "-o",
    "--out-prefix",
    "out_prefix",
    default=None,
    help=(
        "Output prefix (default: input PDB stem). For LEaP PDB: "
        "if omitted with --add-h, <input_stem>_parm.pdb is used."
    ),
)
@click.option(
    "-l",
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
    "--keep-temp/--no-keep-temp",
    "keep_temp",
    default=False,
    show_default=True,
    help="Keep temporary working directory (in current dir) for debugging.",
)
@click.option(
    "--add-ter/--no-add-ter",
    "add_ter",
    default=True,
    show_default=True,
    help=(
        "Insert TER before/after target residues and disconnected peptide blocks. "
        "When target residues are contiguous, TER is not inserted between them."
    ),
)
@click.option(
    "--auto-disulfide/--no-auto-disulfide",
    "auto_disulfide",
    default=True,
    show_default=True,
    help=(
        "Detect disulfides from SG-SG geometry (<= 2.5 A) across CYS/CYM/CYX and bond "
        "them, renaming a bonded CYS to CYX so tleap drops its HG. With "
        "--no-auto-disulfide only residues already named CYX are bonded and CYS is "
        "left untouched."
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
    help="pH used by PDBFixer when adding hydrogens (--add-h). Default: 7.0",
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
    keep_temp: bool,
    add_ter: bool,
    auto_disulfide: bool,
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
        keep_temp=keep_temp,
        add_ter=bool(add_ter),
        auto_disulfide=bool(auto_disulfide),
        add_h=bool(add_h),
        ph=ph,
        ff_set=ff_set,
        out_prefix_given=(out_prefix is not None),
    )
    run_pipeline(args)
