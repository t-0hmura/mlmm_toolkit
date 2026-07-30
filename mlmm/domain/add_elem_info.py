# mlmm/domain/add_elem_info.py

"""
Add/repair PDB element symbols (columns 77-78) without rewriting other records.

Example:
    mlmm add-elem-info -i input.pdb -o fixed.pdb

For detailed documentation, see: docs/add-elem-info.md
"""

from __future__ import annotations

import argparse
import collections
import os
import re
import sys
import time
from pathlib import Path
from typing import Optional, Set

import click

# Residue / ion / water tables live in the L5 foundation layer so the
# L3 domain module can consume them without re-importing the L2 workflows
# layer (which would invert the L1 -> L2 -> {L3, L4} -> L5 direction).
from mlmm.core.residue_data import AMINO_ACIDS, ION, WATER_RES

# Element symbols (IUPAC, 1–118)
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

# Helper: normalize strings to element symbols
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


def _symbol_from_aligned_atom_name(atom_name: str) -> Optional[str]:
    """Infer an element from the four-column PDB atom-name alignment."""
    if len(atom_name) < 4:
        return None
    raw = atom_name[:4]
    if raw[0].isspace():
        return _normalize_symbol(raw.lstrip()[:1])
    if raw[0].isdigit():
        return _normalize_symbol(raw.lstrip("0123456789")[:1])
    return _normalize_symbol(raw[:2])


# Element inference (use residue to disambiguate)
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
    is_protein = res_u in PROTEIN_RES
    is_nucl = res_u in NUCLEIC_RES
    is_water = res_u in WATER_RES

    if res_u in {k.upper() for k in ION.keys()} and not is_nucl:
        # Genuinely polyatomic ions (NH4, H3O+) contain more than one element,
        # so decide per atom name (treat D* as H). Monatomic metal/halogen ions
        # fall through to residue-name resolution below, so an ion whose symbol
        # starts with H/N/O (Na, Ni, Hg, Nd, Hf, He, …) is not mislabelled.
        if res_u in {"NH4", "H3O+", "H3O"}:
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

    if is_protein or is_nucl or is_water:
        # Water: only O and H (treat D* as H)
        if is_water:
            water_name = name_u.lstrip("0123456789")
            if water_name.startswith(("EP", "LP")) or water_name in {"M", "MW"}:
                return "EP"
            if water_name.startswith(("H", "D")):
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

    aligned = _symbol_from_aligned_atom_name(atom_name)
    if aligned is not None:
        return aligned

    # Unaligned programmatic inputs retain the historical prefix fallback.
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

    return None

def _replace_element_field(line: str, symbol: str) -> str:
    if line.endswith("\r\n"):
        content, ending = line[:-2], "\r\n"
    elif line.endswith(("\n", "\r")):
        content, ending = line[:-1], line[-1:]
    else:
        content, ending = line, ""
    content = content.ljust(78)
    return content[:76] + f"{symbol:>2}" + content[78:] + ending


def _default_out_pdb_path(in_pdb: str) -> str:
    path = Path(in_pdb)
    if path.suffix.lower() == ".pdb":
        return str(path.with_name(path.stem + "_add_elem.pdb"))
    return str(path) + "_add_elem.pdb"


def assign_elements(
    in_pdb: str,
    out_pdb: Optional[str],
    overwrite: bool = False,
    inplace: bool = False,
) -> None:
    total = 0
    assigned_new = 0
    overwritten = 0
    kept_existing = 0
    unknown = []
    by_element = collections.Counter()
    with open(in_pdb, "r", encoding="utf-8", errors="surrogateescape", newline="") as handle:
        lines = handle.readlines()

    model_id: object = 0
    rewritten = []
    for line in lines:
        if line.startswith("MODEL"):
            token = line[10:14].strip()
            model_id = int(token) if token.isdigit() else token or model_id
        if not line.startswith(("ATOM  ", "HETATM")):
            rewritten.append(line)
            continue

        total += 1
        previous = line[76:78].strip() if len(line.rstrip("\r\n")) >= 78 else ""
        if previous and not overwrite:
            kept_existing += 1
            rewritten.append(line)
            continue

        atom_name = line[12:16]
        resname = line[17:20]
        symbol = guess_element(atom_name, resname, line.startswith("HETATM"))
        serial_text = line[6:11].strip()
        serial = int(serial_text) if serial_text.isdigit() else None
        if symbol is None:
            unknown.append(
                (
                    model_id,
                    line[21:22].strip(),
                    resname.strip(),
                    line[22:26].strip(),
                    line[26:27].strip(),
                    atom_name.strip(),
                    serial,
                )
            )
            rewritten.append(line)
            continue

        by_element[symbol] += 1
        if previous:
            if previous != symbol:
                overwritten += 1
        else:
            assigned_new += 1
        rewritten.append(_replace_element_field(line, symbol))

    out_path = (
        out_pdb
        if out_pdb
        else (in_pdb if inplace else _default_out_pdb_path(in_pdb))
    )
    with open(out_path, "w", encoding="utf-8", errors="surrogateescape", newline="") as handle:
        handle.writelines(rewritten)

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
        click.echo(
            "[add-elem-info] WARNING: Could not confidently assign "
            f"{len(unknown)} atoms; left unchanged."
        )
        for mid, chid, resn, resseq, icode, aname, serial in unknown[:50]:
            s_str = f" serial {serial}" if serial is not None else ""
            click.echo(
                f"    model {mid} chain {chid} {resn} {resseq}{icode} : "
                f"{aname}{s_str}"
            )
        if len(unknown) > 50:
            click.echo("    ... (truncated) ...")


def main():
    ap = argparse.ArgumentParser(
        description="Add/repair element columns (77–78) in a PDB."
    )
    ap.add_argument("pdb", help="input PDB filepath")
    ap.add_argument(
        "-o",
        "--out",
        help="output PDB filepath (default: <input>_add_elem.pdb)",
    )
    ap.add_argument(
        "--inplace",
        action="store_true",
        help="replace the input file when --out is omitted",
    )
    ap.add_argument(
        "--overwrite",
        action="store_true",
        help=(
            "Re-infer and overwrite element fields even if present "
            "(by default, existing values are preserved)."
        ),
    )
    args = ap.parse_args()

    if not os.path.isfile(args.pdb):
        click.echo(f"[add-elem-info] ERROR: Input not found: {args.pdb}", err=True)
        sys.exit(1)

    try:
        assign_elements(args.pdb, args.out, overwrite=args.overwrite, inplace=args.inplace)
    except Exception as e:
        click.echo(f"[add-elem-info] ERROR: Failed: {e}", err=True)
        sys.exit(2)


# Click subcommand (mlmm add-elem-info)
@click.command(
    help="Add/repair element columns (77–78) in a PDB.",
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
    help="Output PDB filepath (default: <input>_add_elem.pdb; overrides --inplace)",
)
@click.option(
    "--inplace/--no-inplace",
    default=False,
    show_default=True,
    help="Replace the input file when -o/--out is omitted.",
)
@click.option(
    "--overwrite/--no-overwrite",
    "overwrite",
    default=False,
    show_default=True,
    help=(
        "Re-infer and overwrite element fields even if present "
        "(by default, existing values are preserved)."
    ),
)
def cli(in_pdb: Path, out_pdb: Optional[Path], inplace: bool, overwrite: bool) -> None:
    """
    Click wrapper to run via the `mlmm add-elem-info` subcommand.
    """
    time_start = time.perf_counter()
    try:
        assign_elements(
            str(in_pdb),
            (str(out_pdb) if out_pdb else None),
            overwrite=overwrite,
            inplace=inplace,
        )
    except SystemExit as e:
        # Match argparse-like behavior: propagate SystemExit as-is
        raise e
    except Exception as e:
        click.echo(f"[ERR] Failed: {e}", err=True)
        sys.exit(2)
    from mlmm.core.output import emit
    from mlmm.core.utils import format_elapsed

    emit(
        format_elapsed("[time] Elapsed Time for Add Element Info", time_start),
        narrative=True,
    )

if __name__ == "__main__":
    main()
