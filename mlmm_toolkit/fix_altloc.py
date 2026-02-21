#!/usr/bin/env python3
"""
fix_altloc.py - Drop alternate locations from PDB files

What it does
------------
1) Blank the PDB altLoc column (column 17, 1-based) with a single space.
   - This is a 1-character replacement (no shifting / no reformatting).
2) If the same atom appears multiple times due to alternate locations
   (altLoc like A/B/... or custom labels like H/L),
   keep the "best" one by the default rule:
      - Highest occupancy first
      - If tied (or occupancy missing), keep the earliest one in the file

Handled records
---------------
- ATOM / HETATM
- ANISOU is also handled: ANISOU lines are kept only if the corresponding
  ATOM/HETATM line (same serial) is kept.

Notes
-----
- Atom serial numbers are NOT renumbered (gaps may remain).
- CONECT and other connectivity/annotation records are NOT updated.

Usage
-----
  mlmm fix-altloc -i input.pdb -o output.pdb
  mlmm fix-altloc -i ./dir -o ./dir_clean --recursive
  mlmm fix-altloc -i ./dir --inplace --recursive
"""

import argparse
import shutil
from pathlib import Path
from typing import Dict, Iterable, Iterator, List, Optional, Set, Tuple

import click
try:
    from .cli_utils import argparse_bool
except Exception:  # pragma: no cover - support direct execution
    from mlmm_toolkit.cli_utils import argparse_bool

COORD_RECORDS = ("ATOM  ", "HETATM")
ANISOU_RECORD = "ANISOU"

# PDB fixed columns (0-based Python indices)
ALTLOC_IDX = 16                 # column 17 (1-based)
SERIAL_SLICE = slice(6, 11)     # columns 7-11 (1-based), width 5
OCC_SLICE = slice(54, 60)       # columns 55-60 (1-based), width 6


def split_newline(line: str) -> Tuple[str, str]:
    """Split a line into (core, newline) while preserving the newline exactly."""
    if line.endswith("\r\n"):
        return line[:-2], "\r\n"
    if line.endswith("\n"):
        return line[:-1], "\n"
    if line.endswith("\r"):
        return line[:-1], "\r"
    return line, ""


def ensure_len(core: str, n: int) -> str:
    """Right-pad with spaces to guarantee at least n characters (no shifting)."""
    return core if len(core) >= n else core.ljust(n)


def blank_altloc(line: str) -> str:
    """
    Blank the altLoc field (column 17, 1-based) with a single space.

    IMPORTANT: This does NOT remove characters; it replaces exactly one character,
    so the fixed-width PDB formatting is preserved.
    """
    core, nl = split_newline(line)
    core = ensure_len(core, ALTLOC_IDX + 1)  # make sure core[ALTLOC_IDX] exists
    core = core[:ALTLOC_IDX] + " " + core[ALTLOC_IDX + 1:]
    return core + nl


def atom_serial_5(line: str) -> str:
    """Return the 5-character atom serial field exactly as it appears (cols 7-11)."""
    core, _ = split_newline(line)
    core = ensure_len(core, SERIAL_SLICE.stop)
    return core[SERIAL_SLICE]


def parse_occupancy(line: str) -> Optional[float]:
    """
    Parse occupancy from columns 55-60 (1-based).
    Returns None if missing/unparseable.
    """
    core, _ = split_newline(line)
    core = ensure_len(core, OCC_SLICE.stop)
    s = core[OCC_SLICE].strip()
    if not s:
        return None
    try:
        return float(s)
    except ValueError:
        return None


def atom_identity_key(line: str) -> Tuple[str, str, str, str, str, str, str]:
    """
    Build a key to identify the "same atom" while IGNORING altLoc.

    Fields used (fixed columns, classic PDB):
      - record name (ATOM/HETATM)     cols 1-6
      - atom name                      cols 13-16
      - residue name                   cols 18-20
      - chain ID                       col 22
      - residue sequence number        cols 23-26
      - insertion code                 col 27
      - segID (non-standard, common)   cols 73-76

    segID is included to reduce accidental merging in MD-style PDBs where chain ID may be blank.
    """
    core, _ = split_newline(line)
    core = ensure_len(core, 76)

    record = core[0:6]
    atom_name = core[12:16]
    res_name = core[17:20]
    chain_id = core[21:22]
    res_seq = core[22:26]
    i_code = core[26:27]
    seg_id = core[72:76]  # 73-76 (1-based), optional/non-standard

    return (record, atom_name, res_name, chain_id, res_seq, i_code, seg_id)


def process_block(lines: List[str]) -> List[str]:
    """
    Two-pass processing for a block (either the whole file if no MODEL,
    or the content between MODEL and ENDMDL):

    Pass 1: determine the best coordinate line per atom key
            by (occupancy desc, first-appearance asc).
    Pass 2: output only the chosen coordinate lines (altLoc blanked),
            and keep only ANISOU lines whose serial is chosen (altLoc blanked).
            All other records are passed through unchanged.

    Handling of different atom counts between altLoc states:
    --------------------------------------------------------
    When different altLoc states have different atoms (e.g., altLoc A has
    atoms N,CA,CB,CG while altLoc B has N,CA,CB,CD), this function:
    - For DUPLICATE atoms (same identity key, e.g., N,CA,CB): selects the best
      one based on occupancy
    - For UNIQUE atoms (only in one altLoc, e.g., CG in A, CD in B): keeps ALL
      of them in the output

    This ensures the output structure contains all unique atoms from all altLoc
    states, with duplicates resolved to the best conformer.
    """
    # key -> (occ_val_for_compare, line_index, serial5)
    best: Dict[Tuple[str, str, str, str, str, str, str], Tuple[float, int, str]] = {}

    for idx, line in enumerate(lines):
        if line.startswith(COORD_RECORDS):
            key = atom_identity_key(line)
            occ = parse_occupancy(line)
            occ_val = occ if occ is not None else float("-inf")
            serial = atom_serial_5(line)

            if key not in best:
                best[key] = (occ_val, idx, serial)
            else:
                best_occ, best_idx, _best_serial = best[key]
                # Prefer higher occupancy; if tied, prefer earlier line (smaller idx)
                if (occ_val > best_occ) or (occ_val == best_occ and idx < best_idx):
                    best[key] = (occ_val, idx, serial)

    chosen_serials: Set[str] = set(v[2] for v in best.values())

    out: List[str] = []
    for idx, line in enumerate(lines):
        if line.startswith(COORD_RECORDS):
            key = atom_identity_key(line)
            # Keep only the selected "best" line for this key
            if key in best and best[key][1] == idx:
                out.append(blank_altloc(line))
            continue

        if line.startswith(ANISOU_RECORD):
            serial = atom_serial_5(line)
            if serial in chosen_serials:
                out.append(blank_altloc(line))
            continue

        out.append(line)

    return out


def process_stream(lines: Iterable[str]) -> Iterator[str]:
    """
    Handle MODEL/ENDMDL blocks:
      - If MODEL records exist, apply the selection independently within each MODEL block.
      - Text outside MODEL blocks is processed as a single block.
    """
    buffer: List[str] = []
    in_model = False

    for line in lines:
        if line.startswith("MODEL "):
            # Flush anything accumulated before this MODEL
            if buffer:
                for x in process_block(buffer):
                    yield x
                buffer = []
            in_model = True
            yield line
            continue

        if in_model and line.startswith("ENDMDL"):
            # Process the model contents, then emit ENDMDL
            for x in process_block(buffer):
                yield x
            buffer = []
            in_model = False
            yield line
            continue

        buffer.append(line)

    # Flush remaining lines at EOF
    if buffer:
        for x in process_block(buffer):
            yield x


def clean_pdb_file(in_path: Path, out_path: Path) -> None:
    """Process a PDB file and write the cleaned output."""
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with in_path.open("r", newline="") as fin, out_path.open("w", newline="") as fout:
        for out_line in process_stream(fin):
            fout.write(out_line)


def collect_pdb_files(input_path: Path, recursive: bool) -> List[Path]:
    """Collect *.pdb files from a file or directory (optionally recursive)."""
    if input_path.is_file():
        return [input_path]
    pattern = "**/*.pdb" if recursive else "*.pdb"
    return sorted([p for p in input_path.glob(pattern) if p.is_file()])


# =============================================================================
# Public API for programmatic use
# =============================================================================

def has_altloc(pdb_path: Path) -> bool:
    """
    Check if a PDB file contains any non-blank altLoc characters (column 17, 1-based).

    Returns True if at least one ATOM/HETATM record has a non-space character
    in the altLoc column. Returns False if no altLoc is found.
    """
    try:
        with open(pdb_path, "r", encoding="utf-8", errors="ignore") as fh:
            for line in fh:
                if line.startswith(COORD_RECORDS):
                    # altLoc is at column 17 (1-based), which is index 16 (0-based)
                    if len(line) > ALTLOC_IDX:
                        altloc_char = line[ALTLOC_IDX]
                        if altloc_char != " " and altloc_char != "":
                            return True
        return False
    except Exception:
        return False


def fix_altloc_file(
    in_path: str | Path,
    out_path: str | Path,
    *,
    overwrite: bool = False,
    skip_if_no_altloc: bool = True,
) -> bool:
    """
    Fix alternate locations in a PDB file.

    Parameters
    ----------
    in_path : str | Path
        Input PDB file path.
    out_path : str | Path
        Output PDB file path.
    overwrite : bool
        If True, overwrite existing output file. Default False.
    skip_if_no_altloc : bool
        If True, skip processing if no altLoc is detected. Default True.

    Returns
    -------
    bool
        True if the file was processed (altloc found and fixed),
        False if skipped (no altloc detected).

    Raises
    ------
    FileExistsError
        If output file exists and overwrite=False.
    FileNotFoundError
        If input file does not exist.
    """
    in_path = Path(in_path)
    out_path = Path(out_path)

    if not in_path.exists():
        raise FileNotFoundError(f"Input file not found: {in_path}")

    if out_path.exists() and not overwrite:
        raise FileExistsError(f"Output file exists: {out_path}")

    if skip_if_no_altloc and not has_altloc(in_path):
        return False

    clean_pdb_file(in_path, out_path)
    return True


# =============================================================================
# CLI
# =============================================================================

def _fix_altloc_short_help() -> str:
    return "\n".join(
        [
            "Usage: mlmm fix-altloc [OPTIONS]",
            "",
            "Drop alternate locations from PDB files.",
            "",
            "Core options:",
            "  -i, --input PATH                 Input PDB file or directory.",
            "  -o, --out PATH                   Output file/directory.",
            "  --recursive BOOL                 Recurse when input is a directory.",
            "  --inplace BOOL                   Rewrite input file(s) in place.",
            "  --help-advanced                  Show full fix-altloc options and exit.",
            "",
            "Use '--help-advanced' to see all fix-altloc options.",
        ]
    )


def main(argv: Optional[List[str]] = None) -> int:
    def _parse_bool(value: str) -> bool:
        return argparse_bool(value)

    argv_list = list(argv) if argv is not None else None
    if argv_list is not None:
        wants_adv = "--help-advanced" in argv_list
        wants_help = ("--help" in argv_list) or ("-h" in argv_list)
        if wants_help and not wants_adv:
            click.echo(_fix_altloc_short_help())
            return 0
        if wants_adv:
            argv_list = [a for a in argv_list if a != "--help-advanced"]
            if ("--help" not in argv_list) and ("-h" not in argv_list):
                argv_list.append("--help")

    parser = argparse.ArgumentParser(
        description=(
            "Blank PDB altLoc column (col 17) without shifting, and keep one altLoc per atom "
            "by default rule: highest occupancy, then earliest appearance."
        )
    )
    parser.add_argument(
        "-i", "--input", dest="input", type=Path, required=True,
        help="Input PDB file or directory"
    )
    parser.add_argument(
        "-o", "--out", type=Path, default=None,
        help="Output file (if input is a file) or output directory (if input is a directory)."
    )
    parser.add_argument(
        "--recursive",
        type=_parse_bool,
        default=False,
        help="When input is a directory, process *.pdb recursively (including subdirectories). Use True/False.",
    )
    parser.add_argument(
        "--inplace",
        type=_parse_bool,
        default=False,
        help="Overwrite input file(s) in place (creates .bak next to each file). Use True/False.",
    )
    parser.add_argument(
        "--overwrite",
        type=_parse_bool,
        default=False,
        help="Allow overwriting existing output files. Use True/False.",
    )
    parser.add_argument(
        "--force",
        type=_parse_bool,
        default=False,
        help="Process files even if no altLoc is detected (default: skip files without altLoc). Use True/False.",
    )
    args = parser.parse_args(argv_list)

    pdb_files = collect_pdb_files(args.input, args.recursive)
    if not pdb_files:
        parser.error(f"No .pdb files found in: {args.input}")

    skip_if_no_altloc = not args.force
    processed_count = 0
    skipped_count = 0

    # In-place mode
    if args.inplace:
        for in_path in pdb_files:
            if skip_if_no_altloc and not has_altloc(in_path):
                skipped_count += 1
                continue

            bak_path = in_path.with_suffix(in_path.suffix + ".bak")
            if not bak_path.exists():
                shutil.copy2(in_path, bak_path)

            tmp_path = in_path.with_suffix(in_path.suffix + ".tmp")
            clean_pdb_file(in_path, tmp_path)
            tmp_path.replace(in_path)
            processed_count += 1

        if processed_count > 0:
            click.echo(f"[fix-altloc] Processed {processed_count} file(s) in-place.")
        if skipped_count > 0:
            click.echo(f"[fix-altloc] Skipped {skipped_count} file(s) (no altLoc detected).")
        return 0

    # File input
    if args.input.is_file():
        in_path = args.input

        if skip_if_no_altloc and not has_altloc(in_path):
            click.echo(f"[fix-altloc] Skipped {in_path} (no altLoc detected).")
            return 0

        if args.out is None:
            out_path = in_path.with_name(in_path.stem + "_clean.pdb")
        else:
            # If --out ends with .pdb, treat it as a file; otherwise treat it as a directory.
            if args.out.suffix.lower() == ".pdb":
                out_path = args.out
            else:
                args.out.mkdir(parents=True, exist_ok=True)
                out_path = args.out / in_path.name

        if out_path.exists() and not args.overwrite:
            parser.error(f"Output exists: {out_path} (use --overwrite to overwrite)")

        clean_pdb_file(in_path, out_path)
        click.echo(f"[fix-altloc] Fixed altLoc → {out_path}")
        return 0

    # Directory input
    in_dir = args.input
    out_dir = args.out if args.out is not None else in_dir.with_name(in_dir.name + "_clean")
    out_dir.mkdir(parents=True, exist_ok=True)

    for in_path in pdb_files:
        if skip_if_no_altloc and not has_altloc(in_path):
            skipped_count += 1
            continue

        rel = in_path.relative_to(in_dir)
        out_path = out_dir / rel

        if out_path.exists() and not args.overwrite:
            parser.error(f"Output exists: {out_path} (use --overwrite to overwrite)")

        clean_pdb_file(in_path, out_path)
        processed_count += 1

    if processed_count > 0:
        click.echo(f"[fix-altloc] Processed {processed_count} file(s) → {out_dir}")
    if skipped_count > 0:
        click.echo(f"[fix-altloc] Skipped {skipped_count} file(s) (no altLoc detected).")

    return 0


@click.command(
    name="fix-altloc",
    help="Drop alternate locations from PDB, keeping best by occupancy.",
    context_settings={
        "ignore_unknown_options": True,
        "allow_extra_args": True,
        "help_option_names": [],
    },
)
@click.pass_context
def cli(ctx: click.Context) -> None:
    main(list(ctx.args))


if __name__ == "__main__":
    raise SystemExit(main())
