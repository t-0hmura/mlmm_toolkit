"""Import ONIOM input (Gaussian/ORCA) and reconstruct XYZ + layered PDB.

Outputs:
  - <out_prefix>.xyz
  - <out_prefix>_layered.pdb

Layer encoding in output PDB B-factor:
  - ML(QM):       0.00
  - Movable MM:  10.00
  - Frozen MM:   20.00
"""

from __future__ import annotations

import re
from pathlib import Path
from typing import List, Optional, Sequence, Set, Tuple

import click
import numpy as np

from .defaults import BFACTOR_ML, BFACTOR_MOVABLE_MM, BFACTOR_FROZEN


_FLOAT_RE = r"[-+]?\d*\.?\d+(?:[eE][-+]?\d+)?"
_G16_COORD_RE = re.compile(
    rf"^\s*(\S+)\s+(-?\d+)\s+({_FLOAT_RE})\s+({_FLOAT_RE})\s+({_FLOAT_RE})\s+([HL])(?:\s+.*)?$"
)
_ORCA_XYZ_RE = re.compile(rf"^\s*([A-Za-z][A-Za-z]?)\s+({_FLOAT_RE})\s+({_FLOAT_RE})\s+({_FLOAT_RE})\s*$")
_SIX_INT_RE = re.compile(r"^\s*[-+]?\d+(?:\s+[-+]?\d+){5}\s*$")


def _normalize_element_symbol(sym: str) -> str:
    s = re.sub(r"[^A-Za-z]", "", (sym or "").strip())
    if not s:
        return "X"
    if len(s) == 1:
        return s.upper()
    return s[0].upper() + s[1].lower()


def _resolve_mode(mode: Optional[str], input_path: Path) -> str:
    if mode is not None:
        m = str(mode).strip().lower()
        if m in {"g16", "orca"}:
            return m
        raise click.BadParameter("--mode must be one of: g16, orca")

    suf = input_path.suffix.lower()
    if suf in {".gjf", ".com"}:
        return "g16"
    if suf == ".inp":
        return "orca"
    raise click.ClickException(
        f"Could not infer mode from '{input_path.name}'. Use --mode g16|orca "
        "or input extension .gjf/.com/.inp."
    )


def _parse_orca_index_set(raw: str) -> Set[int]:
    """Parse ORCA compact index set (0-based), e.g. "{0:3 7 10:12}"."""
    txt = (raw or "").strip()
    if not txt:
        return set()
    if txt.startswith("{") and txt.endswith("}"):
        txt = txt[1:-1].strip()
    if not txt:
        return set()

    out: Set[int] = set()
    for tok in txt.split():
        if ":" in tok:
            parts = tok.split(":", 1)
            if len(parts) != 2:
                raise ValueError(f"Invalid ORCA range token: '{tok}'")
            a = int(parts[0])
            b = int(parts[1])
            if b < a:
                raise ValueError(f"Invalid descending ORCA range token: '{tok}'")
            out.update(range(a, b + 1))
        else:
            out.add(int(tok))
    return out


def _parse_gaussian_oniom(path: Path) -> Tuple[np.ndarray, List[str], Set[int], Set[int], int, int]:
    """Parse Gaussian ONIOM input produced by mlmm oniom-export.

    Returns:
      coords, elements, qm_indices, movable_indices, qm_charge, qm_mult
    """
    lines = path.read_text(encoding="utf-8", errors="replace").splitlines()

    charge_line_idx = None
    for i, line in enumerate(lines):
        if _SIX_INT_RE.match(line):
            charge_line_idx = i
            break
    if charge_line_idx is None:
        raise click.ClickException("Failed to locate Gaussian ONIOM charge/multiplicity line.")

    parts = lines[charge_line_idx].split()
    if len(parts) != 6:
        raise click.ClickException("Invalid Gaussian ONIOM charge/multiplicity line.")

    try:
        qm_charge = int(parts[2])
        qm_mult = int(parts[3])
    except Exception as exc:
        raise click.ClickException("Failed to parse Gaussian QM charge/multiplicity.") from exc

    # Coordinate block starts after the charge line and optional blanks.
    i = charge_line_idx + 1
    while i < len(lines) and (not lines[i].strip()):
        i += 1

    coords: List[List[float]] = []
    elems: List[str] = []
    qm_indices: Set[int] = set()
    movable_indices: Set[int] = set()

    idx = 0
    while i < len(lines):
        line = lines[i]
        if not line.strip():
            break

        m = _G16_COORD_RE.match(line)
        if m is None:
            # We intentionally support only mlmm-export style rows.
            raise click.ClickException(
                f"Unsupported Gaussian coordinate line format at line {i + 1}: {line!r}"
            )

        atom_token = m.group(1)
        movable = int(m.group(2))
        x = float(m.group(3))
        y = float(m.group(4))
        z = float(m.group(5))
        layer = m.group(6)

        elem_raw = atom_token.split("-", 1)[0]
        elem = _normalize_element_symbol(elem_raw)

        coords.append([x, y, z])
        elems.append(elem)

        if layer == "H":
            qm_indices.add(idx)
            movable_indices.add(idx)
        else:
            if movable == 0:
                movable_indices.add(idx)

        idx += 1
        i += 1

    if not coords:
        raise click.ClickException("No coordinate rows found in Gaussian ONIOM input.")

    return np.asarray(coords, dtype=float), elems, qm_indices, movable_indices, qm_charge, qm_mult


def _parse_orca_qmmm(path: Path) -> Tuple[np.ndarray, List[str], Set[int], Set[int], int, int]:
    """Parse ORCA QM/MM input produced by mlmm oniom-export.

    Returns:
      coords, elements, qm_indices, movable_indices, qm_charge, qm_mult
    """
    lines = path.read_text(encoding="utf-8", errors="replace").splitlines()

    qmmm_start = None
    qmmm_end = None
    for i, line in enumerate(lines):
        if re.match(r"^\s*%qmmm\b", line, flags=re.IGNORECASE):
            qmmm_start = i
            break
    if qmmm_start is None:
        raise click.ClickException("Failed to find %qmmm block in ORCA input.")

    for i in range(qmmm_start + 1, len(lines)):
        if re.match(r"^\s*end\s*$", lines[i], flags=re.IGNORECASE):
            qmmm_end = i
            break
    if qmmm_end is None:
        raise click.ClickException("Failed to find end of %qmmm block in ORCA input.")

    qmmm_lines = lines[qmmm_start : qmmm_end + 1]
    qm_indices: Optional[Set[int]] = None
    active_indices: Optional[Set[int]] = None

    qm_pat = re.compile(r"QMAtoms\s+(\{.*\})\s+end", flags=re.IGNORECASE)
    act_pat = re.compile(r"ActiveAtoms\s+(\{.*\})\s+end", flags=re.IGNORECASE)

    for raw in qmmm_lines:
        m_qm = qm_pat.search(raw)
        if m_qm:
            qm_indices = _parse_orca_index_set(m_qm.group(1))
        m_act = act_pat.search(raw)
        if m_act:
            active_indices = _parse_orca_index_set(m_act.group(1))

    if qm_indices is None:
        raise click.ClickException("Failed to parse QMAtoms from ORCA %qmmm block.")
    if active_indices is None:
        raise click.ClickException("Failed to parse ActiveAtoms from ORCA %qmmm block.")

    xyz_start = None
    qm_charge = 0
    qm_mult = 1
    xyz_header_re = re.compile(r"^\s*\*\s*xyz\s+([-+]?\d+)\s+([-+]?\d+)\s*$", flags=re.IGNORECASE)
    for i, line in enumerate(lines):
        m = xyz_header_re.match(line)
        if m:
            xyz_start = i
            qm_charge = int(m.group(1))
            qm_mult = int(m.group(2))
            break
    if xyz_start is None:
        raise click.ClickException("Failed to find '* xyz <charge> <mult>' block in ORCA input.")

    coords: List[List[float]] = []
    elems: List[str] = []
    i = xyz_start + 1
    while i < len(lines):
        line = lines[i]
        if re.match(r"^\s*\*\s*$", line):
            break
        if not line.strip():
            i += 1
            continue

        m = _ORCA_XYZ_RE.match(line)
        if m is None:
            raise click.ClickException(
                f"Unsupported ORCA xyz line format at line {i + 1}: {line!r}"
            )

        elem = _normalize_element_symbol(m.group(1))
        x = float(m.group(2))
        y = float(m.group(3))
        z = float(m.group(4))

        elems.append(elem)
        coords.append([x, y, z])
        i += 1

    if not coords:
        raise click.ClickException("No coordinates found in ORCA xyz block.")

    n_atoms = len(coords)
    out_of_range = [i for i in qm_indices | active_indices if i < 0 or i >= n_atoms]
    if out_of_range:
        raise click.ClickException(
            f"ORCA QMAtoms/ActiveAtoms contain out-of-range indices for {n_atoms} atoms."
        )

    movable_indices = set(active_indices)
    movable_indices |= set(qm_indices)

    return np.asarray(coords, dtype=float), elems, set(qm_indices), movable_indices, qm_charge, qm_mult


def _bfactor_for_atom(idx: int, qm_indices: Set[int], movable_indices: Set[int]) -> float:
    if idx in qm_indices:
        return float(BFACTOR_ML)
    if idx in movable_indices:
        return float(BFACTOR_MOVABLE_MM)
    return float(BFACTOR_FROZEN)


def _write_xyz(path: Path, coords: np.ndarray, elements: Sequence[str], comment: str = "") -> None:
    n_atoms = int(coords.shape[0])
    lines = [str(n_atoms), comment]
    for i in range(n_atoms):
        e = _normalize_element_symbol(elements[i] if i < len(elements) else "X")
        x, y, z = coords[i]
        lines.append(f"{e:>2s} {x: .8f} {y: .8f} {z: .8f}")
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def _format_pdb_atom_line(
    serial: int,
    atom_name: str,
    res_name: str,
    chain_id: str,
    res_seq: int,
    x: float,
    y: float,
    z: float,
    bfac: float,
    element: str,
) -> str:
    an = (atom_name or "X")[:4]
    rn = (res_name or "MOL")[:3]
    ch = (chain_id or "A")[:1]
    rs = int(res_seq)
    el = _normalize_element_symbol(element)
    return (
        f"ATOM  {serial:5d} {an:>4s} {rn:>3s} {ch:1s}{rs:4d}    "
        f"{x:8.3f}{y:8.3f}{z:8.3f}"
        f"{1.00:6.2f}{bfac:6.2f}          {el:>2s}\n"
    )


def _patch_ref_pdb_line(line: str, x: float, y: float, z: float, bfac: float) -> str:
    s = line.rstrip("\n")
    if len(s) < 80:
        s = s.ljust(80)
    # columns: x[30:38], y[38:46], z[46:54], b[60:66]
    s = s[:30] + f"{x:8.3f}{y:8.3f}{z:8.3f}" + s[54:60] + f"{bfac:6.2f}" + s[66:]
    return s + "\n"


def _write_layered_pdb_without_ref(
    path: Path,
    coords: np.ndarray,
    elements: Sequence[str],
    qm_indices: Set[int],
    movable_indices: Set[int],
) -> None:
    lines: List[str] = []
    for i in range(int(coords.shape[0])):
        x, y, z = coords[i]
        elem = _normalize_element_symbol(elements[i] if i < len(elements) else "X")
        atom_name = elem if len(elem) <= 2 else elem[:2]
        bfac = _bfactor_for_atom(i, qm_indices, movable_indices)
        lines.append(
            _format_pdb_atom_line(
                serial=i + 1,
                atom_name=atom_name,
                res_name="MOL",
                chain_id="A",
                res_seq=1,
                x=float(x),
                y=float(y),
                z=float(z),
                bfac=bfac,
                element=elem,
            )
        )
    lines.append("END\n")
    path.write_text("".join(lines), encoding="utf-8")


def _write_layered_pdb_with_ref(
    path: Path,
    ref_pdb: Path,
    coords: np.ndarray,
    qm_indices: Set[int],
    movable_indices: Set[int],
) -> None:
    ref_lines = ref_pdb.read_text(encoding="utf-8", errors="replace").splitlines(keepends=True)
    atom_line_indices: List[int] = [
        i for i, line in enumerate(ref_lines) if line.startswith(("ATOM  ", "HETATM"))
    ]

    n_atoms = int(coords.shape[0])
    if len(atom_line_indices) != n_atoms:
        raise click.ClickException(
            f"--ref-pdb atom count mismatch: ref has {len(atom_line_indices)} ATOM/HETATM rows, "
            f"but ONIOM input has {n_atoms} atoms."
        )

    out_lines = list(ref_lines)
    for idx, line_idx in enumerate(atom_line_indices):
        x, y, z = coords[idx]
        bfac = _bfactor_for_atom(idx, qm_indices, movable_indices)
        out_lines[line_idx] = _patch_ref_pdb_line(out_lines[line_idx], float(x), float(y), float(z), bfac)

    path.write_text("".join(out_lines), encoding="utf-8")


@click.command(
    name="oniom-import",
    help=(
        "Import ONIOM input (Gaussian g16 or ORCA) and reconstruct XYZ + B-factor layered PDB."
    ),
    context_settings={"help_option_names": ["-h", "--help"]},
)
@click.option(
    "-i",
    "--input",
    "input_path",
    type=click.Path(path_type=Path, exists=True, dir_okay=False),
    required=True,
    help="Input ONIOM file (.gjf/.com for g16, .inp for ORCA).",
)
@click.option(
    "--mode",
    type=click.Choice(["g16", "orca"], case_sensitive=False),
    default=None,
    help="Input mode. If omitted, inferred from input suffix.",
)
@click.option(
    "-o",
    "--out-prefix",
    "out_prefix",
    type=click.Path(path_type=Path),
    default=None,
    help="Output prefix. Defaults to input stem in the current working directory.",
)
@click.option(
    "--ref-pdb",
    type=click.Path(path_type=Path, exists=True, dir_okay=False),
    default=None,
    help="Reference PDB to preserve atom naming/residue metadata (atom count must match).",
)
def cli(
    input_path: Path,
    mode: Optional[str],
    out_prefix: Optional[Path],
    ref_pdb: Optional[Path],
) -> None:
    mode_resolved = _resolve_mode(mode, input_path)

    if out_prefix is None:
        prefix = Path.cwd() / input_path.stem
    else:
        prefix = Path(out_prefix)
    prefix = prefix.resolve()
    prefix.parent.mkdir(parents=True, exist_ok=True)

    if mode_resolved == "g16":
        coords, elements, qm_indices, movable_indices, qm_charge, qm_mult = _parse_gaussian_oniom(input_path)
    else:
        coords, elements, qm_indices, movable_indices, qm_charge, qm_mult = _parse_orca_qmmm(input_path)

    n_atoms = int(coords.shape[0])
    if n_atoms <= 0:
        raise click.ClickException("No atoms parsed from ONIOM input.")

    xyz_path = prefix.with_suffix(".xyz")
    pdb_path = prefix.parent / f"{prefix.name}_layered.pdb"

    _write_xyz(
        xyz_path,
        coords,
        elements,
        comment=(
            f"mode={mode_resolved} atoms={n_atoms} qm={len(qm_indices)} movable={len(movable_indices)} "
            f"q={qm_charge} m={qm_mult}"
        ),
    )

    if ref_pdb is not None:
        _write_layered_pdb_with_ref(pdb_path, ref_pdb, coords, qm_indices, movable_indices)
    else:
        _write_layered_pdb_without_ref(pdb_path, coords, elements, qm_indices, movable_indices)

    click.echo(f"[oniom-import] mode={mode_resolved}")
    click.echo(
        f"[oniom-import] atoms={n_atoms}, qm={len(qm_indices)}, movable={len(movable_indices)}, "
        f"frozen={n_atoms - len(set(movable_indices) | set(qm_indices))}"
    )
    click.echo(f"[oniom-import] wrote: {xyz_path}")
    click.echo(f"[oniom-import] wrote: {pdb_path}")
