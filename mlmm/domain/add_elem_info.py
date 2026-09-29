# mlmm/domain/add_elem_info.py

"""
Add/repair PDB element symbols (columns 77-78) without rewriting other records.

Example:
    mlmm add-elem-info -i input.pdb -o fixed.pdb

For detailed documentation, see: docs/add-elem-info.md
"""

from __future__ import annotations

import collections
import os
import re
import sys
import time
from pathlib import Path
from typing import Optional

import click

# Residue / ion / water tables live in the L5 foundation layer so the
# L3 domain module can consume them without re-importing the L2 workflows
# layer (which would invert the L1 -> L2 -> {L3, L4} -> L5 direction).
from mlmm.core.residue_data import AMINO_ACIDS, ION, WATER_RES

# Element symbols (IUPAC, 1–118)
ELEMENTS: set[str] = {
    "H","He","Li","Be","B","C","N","O","F","Ne","Na","Mg","Al","Si","P","S","Cl","Ar",
    "K","Ca","Sc","Ti","V","Cr","Mn","Fe","Co","Ni","Cu","Zn","Ga","Ge","As","Se","Br","Kr",
    "Rb","Sr","Y","Zr","Nb","Mo","Tc","Ru","Rh","Pd","Ag","Cd","In","Sn","Sb","Te","I","Xe",
    "Cs","Ba","La","Ce","Pr","Nd","Pm","Sm","Eu","Gd","Tb","Dy","Ho","Er","Tm","Yb","Lu",
    "Hf","Ta","W","Re","Os","Ir","Pt","Au","Hg","Tl","Pb","Bi","Po","At","Rn","Fr","Ra",
    "Ac","Th","Pa","U","Np","Pu","Am","Cm","Bk","Cf","Es","Fm","Md","No","Lr","Rf","Db",
    "Sg","Bh","Hs","Mt","Ds","Rg","Cn","Fl","Lv","Ts","Og"
}

# Common residue classes
# SEC, PYL, CSO, CSX and LLP are amino acids whether or not the charge table lists them.
PROTEIN_RES = set(AMINO_ACIDS.keys()) | {"SEC", "PYL", "CSO", "CSX", "LLP"}
NUCLEIC_RES = {
    "DA","DT","DG","DC","DI",
    "A","U","G","C","I",
}


def _format_elapsed(start_time: float) -> str:
    elapsed = max(0.0, time.perf_counter() - start_time)
    hours, rem = divmod(elapsed, 3600)
    minutes, seconds = divmod(rem, 60)
    return (
        "[time] Elapsed Time for Add Element Info: "
        f"{int(hours):02d}:{int(minutes):02d}:{seconds:06.3f}"
    )


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
    if letters[0].upper() == "D":
        return "H"
    return None

def _symbol_from_resname(resname: str) -> Optional[str]:
    """Extract an element symbol from an ion residue name (e.g., CA, FE2, Cl-, YB2)."""
    res = resname.strip()
    return _normalize_symbol(res)


def _symbol_from_aligned_atom_name(atom_name: str) -> Optional[str]:
    """Infer an element from the four-column PDB atom-name alignment."""
    if len(atom_name) < 4:
        return None
    raw = atom_name[:4]
    if raw[0].isspace():
        # LEaP writes GAFF chlorine and bromine from column 14 (" CL1", " BR1").
        if raw[1:3].upper() in {"CL", "BR"}:
            return raw[1:3].capitalize()
        return _normalize_symbol(raw.lstrip()[:1])
    if raw[0].isdigit():
        return _normalize_symbol(raw.lstrip("0123456789")[:1])
    # Four-character hydrogen names start in column 13 (HG11 is H, not Hg).
    if raw[0] in "Hh" and not any(ch.isspace() for ch in raw):
        return "H"
    return _normalize_symbol(raw[:2])


# Element inference (use residue to disambiguate)
def guess_element(atom_name: str, resname: str, _is_het: bool = False) -> Optional[str]:
    """
    Infer the element from atom name + residue name.
    Pass the raw four-character PDB atom-name field (columns 13-16): its column
    alignment tells `` NA `` (N) from ``NA  `` (Na).
    Priority:
      1) Ion residues: prefer the residue name (NH4 / H3O+ handled per-atom as H/N/O)
      2) Polymers (protein/nucleic acid) and water: follow convention (H/C/N/O/S/P/Se)
      3) Other ligands: PDB column alignment; shorter names use H/C/P prefixes, then normalization
      4) Unresolved → None
    """
    name_u = atom_name.strip().upper()
    res_u = resname.strip().upper()
    is_protein = res_u in PROTEIN_RES
    is_nucl = res_u in NUCLEIC_RES
    is_water = res_u in WATER_RES

    # Residue I is inosine, not iodide.
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

    if is_protein or is_nucl or is_water:
        # Water: O and H (D* is H; leading digits as in 1HW are skipped); virtual sites are EP
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

        # Fallback to normalization (rare halogens or atypical labels)
        sym = _normalize_symbol(name_u)
        if sym:
            return sym

    aligned = _symbol_from_aligned_atom_name(atom_name)
    if aligned is not None:
        return aligned

    # Names shorter than four characters carry no column alignment.
    if name_u.startswith(("H", "D")):
        return "H"
    if name_u.startswith("C") and not name_u.startswith("CL"):
        return "C"
    if name_u.startswith("P"):
        return "P"

    # Metals and halogens often appear as the atom name (FE, ZN, MG, CL, BR, I, F, ...)
    sym = _normalize_symbol(name_u)
    if sym:
        return sym

    return None


def _replace_element_field(line: str, symbol: str, *, field_offset: int = 0) -> str:
    if line.endswith("\r\n"):
        content, ending = line[:-2], "\r\n"
    elif line.endswith(("\n", "\r")):
        content, ending = line[:-1], line[-1:]
    else:
        content, ending = line, ""
    start = 76 + int(field_offset)
    content = content.ljust(start + 2)
    return content[:start] + f"{symbol:>2}" + content[start + 2:] + ending


def pdb_decimal_overflow_shifts(line: str) -> tuple[int, int]:
    """Return atom-serial and residue-number decimal overflow widths."""

    serial_shift = 0
    while 11 + serial_shift < len(line) and line[11 + serial_shift].isdigit():
        serial_shift += 1
    residue_shift = 0
    overflow_end = 26 + serial_shift
    if overflow_end < len(line) and line[overflow_end].isdigit():
        while overflow_end < len(line) and line[overflow_end].isdigit():
            overflow_end += 1
        overflow_value = line[22 + serial_shift : overflow_end].strip()
        if re.fullmatch(r"[-+]?\d{5,}", overflow_value):
            residue_shift = overflow_end - (26 + serial_shift)
    return serial_shift, residue_shift


def _default_out_pdb_path(in_pdb: str) -> str:
    path = Path(in_pdb)
    if path.suffix.lower() == ".pdb":
        return str(path.with_name(path.stem + "_add_elem.pdb"))
    return str(path) + "_add_elem.pdb"


def assign_elements(in_pdb: str, out_pdb: Optional[str], overwrite: bool = False) -> None:
    # If an explicit output path is provided, never overwrite in-place even when --overwrite is
    # passed. This keeps -o/--out as the higher-priority choice.
    effective_overwrite = overwrite and out_pdb is None
    if effective_overwrite:
        out_path = in_pdb
    else:
        out_path = out_pdb if out_pdb else _default_out_pdb_path(in_pdb)
    input_path = Path(in_pdb).expanduser().resolve()
    output_path = Path(out_path).expanduser().resolve()
    aliases_input = input_path == output_path
    if not aliases_input and input_path.exists() and output_path.exists():
        aliases_input = os.path.samefile(input_path, output_path)
    if aliases_input and not overwrite:
        raise ValueError(
            "Output physically aliases the input; use --overwrite to replace it."
        )

    total = 0
    assigned_or_updated = 0
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
        serial_shift, residue_shift = pdb_decimal_overflow_shifts(line)
        field_offset = serial_shift + residue_shift
        atom_name = line[12 + serial_shift : 16 + serial_shift]
        resname = line[17 + serial_shift : 20 + serial_shift]
        symbol = guess_element(atom_name, resname, line.startswith("HETATM"))
        serial_text = line[6 : 11 + serial_shift].strip()
        serial = int(serial_text) if serial_text.isdigit() else None
        if symbol is None:
            unknown.append(
                (
                    model_id,
                    line[21 + serial_shift : 22 + serial_shift].strip(),
                    resname.strip(),
                    line[22 + serial_shift : 26 + field_offset].strip(),
                    line[26 + field_offset : 27 + field_offset].strip(),
                    atom_name.strip(),
                    serial,
                )
            )
            rewritten.append(line)
            continue

        previous = line[76 + field_offset : 78 + field_offset].strip()
        by_element[symbol] += 1
        if previous != symbol:
            assigned_or_updated += 1
        rewritten.append(_replace_element_field(line, symbol, field_offset=field_offset))

    with open(out_path, "w", encoding="utf-8", errors="surrogateescape", newline="") as handle:
        handle.writelines(rewritten)

    # Summary
    click.echo(f"[add-elem-info] Wrote: {out_path}")
    click.echo(f"  total atoms                 : {total}")
    click.echo(f"  assigned/updated            : {assigned_or_updated}")
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
    "-o", "--output", "--out",
    "out_pdb",
    type=click.Path(path_type=Path, dir_okay=False),
    default=None,
    help='Output PDB filepath (default: replace ".pdb" with "_add_elem.pdb"; when provided, --overwrite is ignored unless this is the input file, which requires it).',
)
@click.option(
    "--overwrite/--no-overwrite",
    default=False,
    show_default=True,
    help="Overwrite the input file in-place when -o/--out is omitted.",
)
@click.option(
    "--inplace/--no-inplace",
    "inplace",
    default=None,
    show_default="no-inplace",
    hidden=True,
    help="Alias of --overwrite.",
)
def cli(
    in_pdb: Path,
    out_pdb: Optional[Path],
    overwrite: bool,
    inplace: Optional[bool],
) -> None:
    """
    Click wrapper to run via the `mlmm add-elem-info` subcommand.
    """
    if inplace is not None:
        source = click.get_current_context().get_parameter_source("overwrite")
        if source is click.core.ParameterSource.COMMANDLINE and overwrite != inplace:
            raise click.UsageError(
                "Conflicting values were supplied through aliases: "
                "--overwrite, --inplace."
            )
        overwrite = inplace
    time_start = time.perf_counter()
    try:
        assign_elements(
            str(in_pdb),
            (str(out_pdb) if out_pdb else None),
            overwrite=overwrite,
        )
    except SystemExit as e:
        # Match argparse-like behavior: propagate SystemExit as-is
        raise e
    except Exception as e:
        click.echo(f"[ERR] Failed: {e}", err=True)
        sys.exit(2)
    from mlmm.core.output import emit

    emit(
        _format_elapsed(time_start),
        narrative=True,
    )
