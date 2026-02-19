# mlmm_toolkit/add_elem_info.py

"""
Add/repair PDB element symbols (columns 77-78) using Biopython inference.

Example:
    mlmm add-elem-info -i input.pdb -o fixed.pdb

For detailed documentation, see: docs/add_elem_info.md
"""

from __future__ import annotations

import argparse
import collections
import os
import re
import sys
from pathlib import Path
from typing import Optional, Set

import click
from Bio.PDB import PDBParser, PDBIO

# Reuse residue/ion dictionaries from extract.py to keep definitions in sync
from .extract import AMINO_ACIDS, ION, WATER_RES

# -----------------------------
# Element symbols (IUPAC, 1–118)
# -----------------------------
ELEMENTS: Set[str] = {
    "H","He","Li","Be","B","C","N","O","F","Ne","Na","Mg","Al","Si","P","S","Cl","Ar",
    "K","Ca","Sc","Ti","V","Cr","Mn","Fe","Co","Ni","Cu","Zn","Ga","Ge","As","Se","Br","Kr",
    "Rb","Sr","Y","Zr","Nb","Mo","Tc","Ru","Rh","Pd","Ag","Cd","In","Sn","Sb","Te","I","Xe",
    "Cs","Ba","La","Ce","Pr","Nd","Pm","Sm","Eu","Gd","Tb","Dy","Ho","Er","Tm","Yb","Lu",
    "Hf","Ta","W","Re","Os","Ir","Pt","Au","Hg","Tl","Pb","Bi","Po","At","Rn","Fr","Ra",
    "Ac","Th","Pa","U","Np","Pu","Am","Cm","Bk","Cf","Es","Fm","Md","No","Lr","Rf","Db",
    "Sg","Bh","Hs","Mt","Ds","Rg","Cn","Fl","Lv","Ts","Og"
}

# Common residue classes
PROTEIN_RES = set(AMINO_ACIDS.keys())
NUCLEIC_RES = {
    # DNA/RNA (minimum set)
    "DA","DT","DG","DC","DI",
    "A","U","G","C","I",
}

# -----------------------------
# Helper: normalize strings to element symbols
# -----------------------------
_re_letters = re.compile(r"[A-Za-z]+")

def _normalize_symbol(s: str) -> Optional[str]:
    """Remove non-letters; prefer a 2-letter match, then 1-letter, against known elements.
    Returns the correctly cased symbol if matched.
    Treat deuterium 'D' as hydrogen 'H' (PDB often uses D interchangeably with H).
    """
    if not s:
        return None
    m = _re_letters.findall(s)
    if not m:
        return None
    letters = "".join(m)
    if len(letters) >= 2:
        cand2 = (letters[:2][0].upper() + letters[:2][1].lower())
        if cand2 in ELEMENTS:
            return cand2
    cand1 = letters[0].upper()
    if cand1 in ELEMENTS:
        return cand1
    # Deuterium -> Hydrogen fallback
    if letters[0].upper() == "D":
        return "H"
    return None

def _symbol_from_resname(resname: str) -> Optional[str]:
    """
    Extract an element symbol from an ion residue name (e.g., CA, FE2, Cl-, YB2, IOD).
    """
    res = resname.strip()
    sym = _normalize_symbol(res)
    if sym is None and res.upper().startswith("IOD"):
        sym = "I"
    return sym

# -----------------------------
# Element inference (use residue to disambiguate)
# -----------------------------
def guess_element(atom_name: str, resname: str, is_het: bool) -> Optional[str]:
    """
    Infer the element from atom name + residue name.
    Priority:
      1) Ion residues: prefer the residue name (NH4 / H3O+ handled per-atom as H/N/O)
      2) Polymers (protein/nucleic acid) and water: follow convention (H/C/N/O/S/P/Se)
         - e.g., CA = Carbon (Cα), HG = Hydrogen, etc.
      3) Other ligands: use atom-name prefix; prioritize Carbon for C* (except CL) and P for P*
      4) Fallback to 2-letter then 1-letter normalization; return None if still ambiguous
    """
    name_u = atom_name.strip().upper()
    res_u = resname.strip().upper()

    # 1) Ion residues — strongly prefer the residue-derived element
    if res_u in {k.upper() for k in ION.keys()}:
        # Polyatomic ions (NH4, H3O+, …): decide per atom name (treat D* as H)
        if name_u.startswith(("H", "D")):
            return "H"
        if name_u.startswith("N"):
            return "N"
        if name_u.startswith("O"):
            return "O"
        # Monatomic metals/halogens: from residue name
        sym = _symbol_from_resname(res_u)
        if sym:
            return sym
        # If residue is atypical, allow atom-name halogens (CL/BR/I/F)
        if name_u.startswith("CL"):
            return "Cl"
        if name_u.startswith("BR"):
            return "Br"
        if name_u.startswith("I"):
            return "I"
        if name_u.startswith("F"):
            return "F"

    # 2) Polymers (protein/nucleic) / water
    is_protein = res_u in PROTEIN_RES
    is_nucl = res_u in NUCLEIC_RES
    is_water = res_u in WATER_RES
    if is_protein or is_nucl or is_water:
        # Water: only O and H (treat D* as H)
        if is_water:
            if name_u.startswith(("H", "D")):
                return "H"
            return "O"

        # Hydrogen (including D*)
        if name_u.startswith(("H", "D")):
            return "H"

        # Selenium (e.g., selenomethionine/selenocysteine)
        if name_u.startswith("SE"):
            return "Se"

        # P, N, O, S map directly by first letter
        if name_u.startswith("P"):
            return "P"
        if name_u.startswith("N"):
            return "N"
        if name_u.startswith("O"):
            return "O"
        if name_u.startswith("S"):
            return "S"

        # Carbon for Cα/sidechain labels (CA, CB, CG, CD, CE, CZ, CH*, etc.)
        if name_u.startswith("C"):
            return "C"

        # Rare halogens in polymers: final fallback to normalization
        sym = _normalize_symbol(name_u)
        if sym:
            return sym

    # 3) Non-polymers (ligands / cofactors)
    #    Hydrogen (including D*) for ligands/cofactors
    if name_u.startswith(("H", "D")):
        return "H"
    #    Carbon/Phosphorus-like labels (C*, P*) -> C/P (exclude CL)
    if name_u.startswith("C") and not name_u.startswith("CL"):
        return "C"
    if name_u.startswith("P"):
        return "P"

    # Metals and halogens often appear as the atom name (FE, ZN, MG, HG, CL, BR, I, F ...)
    sym = _normalize_symbol(name_u)
    if sym:
        return sym

    # 4) Unresolved
    return None

# -----------------------------
# Detect whether the input originally had element fields,
# keyed by atom serial number (columns 7–11)
# -----------------------------
def scan_existing_elements_by_serial(pdb_path: str) -> Set[int]:
    """
    Scan the raw PDB lines and return the serial numbers of ATOM/HETATM records whose
    element field (columns 77–78) was non-empty in the original file.
    This avoids Biopython side effects and reflects the true presence/absence in the input.
    """
    serials_with_elem: Set[int] = set()
    try:
        with open(pdb_path, "r", encoding="utf-8", errors="ignore") as fh:
            for line in fh:
                if not (line.startswith("ATOM") or line.startswith("HETATM")):
                    continue
                if len(line) < 78:
                    # No element field present
                    continue
                serial_str = line[6:11].strip()
                elem_raw = line[76:78].strip()
                if not serial_str:
                    continue
                try:
                    serial = int(serial_str)
                except ValueError:
                    continue
                # If non-empty, consider that the original file had an element entry
                # (keep isotopic labels like D as-is)
                if elem_raw:
                    serials_with_elem.add(serial)
    except Exception:
        # If the file can't be read, return empty (treat all as unset)
        pass
    return serials_with_elem

def _get_atom_serial(atom) -> Optional[int]:
    """
    Safely obtain the serial number from a Biopython Atom, handling version differences.
    """
    sn = getattr(atom, "serial_number", None)
    if sn is None and hasattr(atom, "get_serial_number"):
        try:
            sn = atom.get_serial_number()
        except Exception:
            sn = None
    return sn

# -----------------------------
# Main processing
# -----------------------------
def assign_elements(in_pdb: str, out_pdb: Optional[str], overwrite: bool = False) -> None:
    # Scan the input file for the original presence of element fields
    existing_by_serial = scan_existing_elements_by_serial(in_pdb)

    parser = PDBParser(QUIET=True)
    structure_id = os.path.splitext(os.path.basename(in_pdb))[0]
    structure = parser.get_structure(structure_id, in_pdb)

    total = 0
    assigned_new = 0          # newly set for atoms that lacked an element field
    overwritten = 0           # element existed originally but was re-inferred due to --overwrite
    kept_existing = 0         # element existed originally and was preserved (no --overwrite)
    unknown = []              # could not infer (left unchanged)

    by_element = collections.Counter()

    for model in structure:
        for chain in model:
            for residue in chain:
                hetflag = residue.id[0].strip()  # '' (empty) = standard; 'W' = water; 'H_' = HETATM
                is_het = (hetflag != "")
                resname = residue.get_resname()
                for atom in residue:
                    total += 1
                    name = atom.get_name()

                    serial = _get_atom_serial(atom)
                    had_element_in_input = (serial in existing_by_serial) if serial is not None else False

                    if had_element_in_input and not overwrite:
                        kept_existing += 1
                        continue  # Respect existing element: do not modify without --overwrite

                    sym = guess_element(name, resname, is_het)
                    if sym is None:
                        unknown.append((model.id, chain.id, residue.id, resname, name, serial))
                        # If inference failed: keep the previous value (if any), otherwise leave unset
                        continue

                    # Biopython uses atom.element to populate columns 77–78 on output
                    prev = getattr(atom, "element", None)
                    atom.element = sym
                    by_element[sym] += 1
                    if had_element_in_input:
                        if prev != sym:
                            overwritten += 1
                    else:
                        assigned_new += 1

    io = PDBIO()
    io.set_structure(structure)
    out_path = out_pdb if out_pdb else in_pdb  # overwrite input if not specified
    io.save(out_path)

    # Summary
    click.echo(f"[add-elem-info] Wrote: {out_path}")
    click.echo(f"  total atoms                 : {total}")
    click.echo(f"  newly assigned              : {assigned_new}")
    click.echo(f"  kept existing (no overwrite): {kept_existing}")
    click.echo(f"  overwritten (--overwrite)   : {overwritten}")
    if by_element:
        top = ", ".join(f"{k}:{v}" for k, v in by_element.most_common())
        click.echo(f"  assignment breakdown        : {top}")
    if unknown:
        click.echo(f"[add-elem-info] WARNING: Could not confidently assign {len(unknown)} atoms; left unchanged.")
        for (mid, chid, resid, resn, aname, serial) in unknown[:50]:
            if isinstance(resid, tuple):
                resseq = resid[1]
                icode = resid[2].strip()
            else:
                resseq, icode = "?", ""
            s_str = f" serial {serial}" if serial is not None else ""
            click.echo(f"    model {mid} chain {chid} {resn} {resseq}{icode} : {aname}{s_str}")
        if len(unknown) > 50:
            click.echo("    ... (truncated) ...")

def main():
    ap = argparse.ArgumentParser(
        description="Add/repair element columns (77–78) in a PDB using Biopython."
    )
    ap.add_argument("pdb", help="input PDB filepath")
    ap.add_argument("-o", "--out", help="output PDB filepath (omit to overwrite input)")
    ap.add_argument(
        "--overwrite",
        action="store_true",
        help="Re-infer and overwrite element fields even if present (by default, existing values are preserved).",
    )
    args = ap.parse_args()

    if not os.path.isfile(args.pdb):
        click.echo(f"[add-elem-info] ERROR: Input not found: {args.pdb}", err=True)
        sys.exit(1)

    try:
        assign_elements(args.pdb, args.out, overwrite=args.overwrite)
    except Exception as e:
        click.echo(f"[add-elem-info] ERROR: Failed: {e}", err=True)
        sys.exit(2)

# -----------------------------
# Click subcommand (mlmm add-elem-info)
# -----------------------------
@click.command(
    help="Add/repair element columns (77–78) in a PDB using Biopython.",
    context_settings={"help_option_names": ["-h", "--help"]},
)
@click.option(
    "-i", "--input",
    "in_pdb",
    type=click.Path(path_type=Path, exists=True, dir_okay=False),
    required=True,
    help="Input PDB filepath",
)
@click.option(
    "-o", "--out",
    "out_pdb",
    type=click.Path(path_type=Path, dir_okay=False),
    default=None,
    help="Output PDB filepath (omit to overwrite input)",
)
@click.option(
    "--overwrite",
    is_flag=True,
    help="Re-infer and overwrite element fields even if present (by default, existing values are preserved).",
)
def cli(in_pdb: Path, out_pdb: Optional[Path], overwrite: bool) -> None:
    """
    Click wrapper to run via the `mlmm add-elem-info` subcommand.
    """
    try:
        assign_elements(str(in_pdb), (str(out_pdb) if out_pdb else None), overwrite=overwrite)
    except SystemExit as e:
        # Match argparse-like behavior: propagate SystemExit as-is
        raise e
    except Exception as e:
        click.echo(f"[ERR] Failed: {e}", err=True)
        sys.exit(2)

if __name__ == "__main__":
    main()
