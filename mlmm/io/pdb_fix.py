#!/usr/bin/env python3
"""
pdb_fix.py - Drop alternate locations from PDB files

What it does
------------
1) Blank the PDB altLoc column (column 17, 1-based) with a single space.
   - This is a 1-character replacement (no shifting / no reformatting).
2) Select one coherent non-blank altLoc label per residue (A/B/... or custom
   labels such as H/L):
      - Highest mean occupancy across that residue's labelled atoms
      - A label with no parsed occupancies ranks below any parsed mean
      - Break equal scores by earliest appearance (including all-missing cases)
   Blank/shared atoms are retained except when the selected label supplies the
   same atom identity. Atoms that exist only in an unselected conformer are
   dropped instead of being merged into a hybrid residue.

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

import shutil
from pathlib import Path
from typing import Dict, Iterable, Iterator, List, Optional, Set, Tuple

import click

from mlmm.io.altloc import choose_altloc_label, occupancy_choice_key

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


def residue_identity_key(line: str) -> Tuple[str, str, str, str, str]:
    """Return the residue identity used for coherent altLoc selection."""
    core, _ = split_newline(line)
    core = ensure_len(core, 76)
    return (
        core[17:20],
        core[21:22],
        core[22:26],
        core[26:27],
        core[72:76],
    )


def altloc_label(line: str) -> str:
    """Return the stripped one-character altLoc label (empty for blank)."""
    core, _ = split_newline(line)
    core = ensure_len(core, ALTLOC_IDX + 1)
    return core[ALTLOC_IDX].strip()


def process_block(lines: List[str]) -> List[str]:
    """
    Multi-pass processing for a block (either the whole file if no MODEL,
    or the content between MODEL and ENDMDL):

    Pass 1: choose one non-blank altLoc label per residue by mean occupancy,
            breaking ties by first appearance.
    Pass 2: retain blank/shared records and records carrying the selected label,
            except when the selected label supplies the same atom identity;
            choose at most one coordinate line per identity.
    Pass 3: output the selected coordinates with altLoc blanked and keep only
            matching ANISOU records. Other record types pass through unchanged.

    Selecting at residue scope prevents combining coordinates from different
    deposited conformers into a hybrid residue.
    """
    label_observations: Dict[
        Tuple[str, str, str, str, str],
        List[Tuple[str, Optional[float], int]],
    ] = {}

    for idx, line in enumerate(lines):
        if not line.startswith(COORD_RECORDS):
            continue
        label = altloc_label(line)
        if not label:
            continue
        residue = residue_identity_key(line)
        label_observations.setdefault(residue, []).append(
            (label, parse_occupancy(line), idx)
        )

    selected_labels: Dict[Tuple[str, str, str, str, str], str] = {}
    for residue, observations in label_observations.items():
        selected_labels[residue] = choose_altloc_label(observations)

    # A labelled version of an atom supersedes a blank/shared record, matching
    # the typed structure reader.  Other blank/shared atoms remain untouched.
    chosen_label_keys: Set[Tuple[str, str, str, str, str, str, str]] = set()
    for line in lines:
        if not line.startswith(COORD_RECORDS):
            continue
        label = altloc_label(line)
        selected = selected_labels.get(residue_identity_key(line))
        if label and label == selected:
            chosen_label_keys.add(atom_identity_key(line))

    # atom key -> (occupancy rank, line index, serial field)
    best: Dict[
        Tuple[str, str, str, str, str, str, str],
        Tuple[Tuple[int, float, int], int, str],
    ] = {}

    for idx, line in enumerate(lines):
        if line.startswith(COORD_RECORDS):
            label = altloc_label(line)
            selected = selected_labels.get(residue_identity_key(line))
            if label and label != selected:
                continue
            key = atom_identity_key(line)
            if not label and key in chosen_label_keys:
                continue
            occ = parse_occupancy(line)
            serial = atom_serial_5(line)

            if key not in best:
                rank = occupancy_choice_key(occ, idx)
                best[key] = (rank, idx, serial)
            else:
                best_rank, best_idx, _best_serial = best[key]
                # Resolve malformed duplicates after coherent label selection.
                rank = occupancy_choice_key(occ, idx)
                if rank > best_rank:
                    best[key] = (rank, idx, serial)

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


def _run_fix_altloc(
    input_path: Path,
    out: Optional[Path],
    recursive: bool,
    inplace: bool,
    overwrite: bool,
    force: bool,
) -> None:
    """Core business logic for fix-altloc (called from Click CLI)."""
    pdb_files = collect_pdb_files(input_path, recursive)
    if not pdb_files:
        raise click.ClickException(f"No .pdb files found in: {input_path}")

    skip_if_no_altloc = not force
    processed_count = 0
    skipped_count = 0

    # In-place mode
    if inplace:
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
        return

    # File input
    if input_path.is_file():
        in_path = input_path

        if skip_if_no_altloc and not has_altloc(in_path):
            click.echo(f"[fix-altloc] Skipped {in_path} (no altLoc detected).")
            return

        if out is None:
            out_path = in_path.with_name(in_path.stem + "_clean.pdb")
        else:
            if out.suffix.lower() == ".pdb":
                out_path = out
            else:
                out.mkdir(parents=True, exist_ok=True)
                out_path = out / in_path.name

        if out_path.exists() and not overwrite:
            raise click.ClickException(f"Output exists: {out_path} (use --overwrite to overwrite)")

        clean_pdb_file(in_path, out_path)
        click.echo(f"[fix-altloc] Fixed altLoc → {out_path}")
        return

    # Directory input
    in_dir = input_path
    out_dir = out if out is not None else in_dir.with_name(in_dir.name + "_clean")
    out_dir.mkdir(parents=True, exist_ok=True)

    for in_path in pdb_files:
        if skip_if_no_altloc and not has_altloc(in_path):
            skipped_count += 1
            continue

        rel = in_path.relative_to(in_dir)
        out_path = out_dir / rel

        if out_path.exists() and not overwrite:
            raise click.ClickException(f"Output exists: {out_path} (use --overwrite to overwrite)")

        clean_pdb_file(in_path, out_path)
        processed_count += 1

    if processed_count > 0:
        click.echo(f"[fix-altloc] Processed {processed_count} file(s) → {out_dir}")
    if skipped_count > 0:
        click.echo(f"[fix-altloc] Skipped {skipped_count} file(s) (no altLoc detected).")


@click.command(
    name="fix-altloc",
    help=(
        "Blank PDB altLoc column (col 17) without shifting, and select one "
        "coherent label per residue by highest mean occupancy, then earliest appearance."
    ),
    context_settings={"help_option_names": ["-h", "--help"]},
)
@click.option(
    "-i", "--input", "input_path",
    type=click.Path(exists=True, path_type=Path),
    required=True,
    help="Input PDB file or directory.",
)
@click.option(
    "-o", "--out",
    type=click.Path(path_type=Path),
    default=None,
    help="Output file (if input is a file) or output directory (if input is a directory).",
)
@click.option(
    "--recursive/--no-recursive",
    default=False, show_default=True,
    help="When input is a directory, process *.pdb recursively (including subdirectories).",
)
@click.option(
    "--inplace/--no-inplace",
    default=False, show_default=True,
    help="Overwrite input file(s) in place (creates .bak next to each file).",
)
@click.option(
    "--overwrite/--no-overwrite",
    default=False, show_default=True,
    help="Allow overwriting existing output files.",
)
@click.option(
    "--force/--no-force",
    default=False, show_default=True,
    help="Process files even if no altLoc is detected (default: skip files without altLoc).",
)
def cli(
    input_path: Path,
    out: Optional[Path],
    recursive: bool,
    inplace: bool,
    overwrite: bool,
    force: bool,
) -> None:
    _run_fix_altloc(input_path, out, recursive, inplace, overwrite, force)


if __name__ == "__main__":
    cli()
