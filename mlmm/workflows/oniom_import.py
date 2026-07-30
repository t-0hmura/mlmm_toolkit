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
import tempfile
import time
from pathlib import Path
from typing import List, Optional, Sequence, Set, Tuple

import click
import numpy as np

from mlmm.core.defaults import BFACTOR_ML, BFACTOR_MOVABLE_MM, BFACTOR_FROZEN
from mlmm.core.result_commit import commit_payloads


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
    elements: List[str],
    qm_indices: Set[int],
    movable_indices: Set[int],
    *,
    expected_ref_order_digest: Optional[str] = None,
    allow_unverified_ref_order: bool = False,
) -> str:
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
    try:
        from mlmm.io.structure_formats import read_pdb_atom_sites

        ref_elements = [
            _normalize_element_symbol(record.element)
            for record in read_pdb_atom_sites(ref_pdb, warn_altloc=False)[0]
        ]
    except Exception as exc:
        raise click.ClickException(
            f"Failed to read ordered elements from --ref-pdb: {exc}"
        ) from exc
    oniom_elements = [
        _normalize_element_symbol(symbol) for symbol in elements
    ]
    if ref_elements != oniom_elements:
        mismatch = next(
            (
                index
                for index, (ref_element, oniom_element) in enumerate(
                    zip(ref_elements, oniom_elements)
                )
                if ref_element != oniom_element
            ),
            0,
        )
        raise click.ClickException(
            "--ref-pdb atom-order element mismatch at 1-based atom "
            f"{mismatch + 1}: reference={ref_elements[mismatch]}, "
            f"ONIOM={oniom_elements[mismatch]}."
        )
    from mlmm.io.oniom_identity import elements_are_unique, pdb_order_digest

    if expected_ref_order_digest is not None:
        current_digest = pdb_order_digest(ref_pdb)
        if current_digest != expected_ref_order_digest:
            raise click.ClickException(
                "--ref-pdb atom identity/order does not match the reference "
                "embedded by `mlmm oniom-export`."
            )
        verification = "identity-verified"
    elif elements_are_unique(oniom_elements):
        verification = "element-verified"
    elif not allow_unverified_ref_order:
        raise click.ClickException(
            "The ONIOM input has no embedded reference-order identity and "
            "contains repeated elements, so --ref-pdb atom order cannot be "
            "verified. Independently verify the order, then pass "
            "--allow-unverified-ref-order to use legacy positional mapping."
        )
    else:
        verification = "unverified-opt-in"
        click.echo(
            "[oniom-import] WARNING: --ref-pdb atom order is not identity-verified; "
            "using explicitly requested positional mapping.",
            err=True,
        )

    out_lines = list(ref_lines)
    for idx, line_idx in enumerate(atom_line_indices):
        x, y, z = coords[idx]
        bfac = _bfactor_for_atom(idx, qm_indices, movable_indices)
        out_lines[line_idx] = _patch_ref_pdb_line(out_lines[line_idx], float(x), float(y), float(z), bfac)

    path.write_text("".join(out_lines), encoding="utf-8")
    return verification


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
@click.option(
    "--allow-unverified-ref-order/--no-allow-unverified-ref-order",
    default=False,
    show_default=True,
    help="Allow legacy positional --ref-pdb mapping when repeated elements make "
         "atom identity unverifiable. Use only after independently checking order.",
)
def cli(
    input_path: Path,
    mode: Optional[str],
    out_prefix: Optional[Path],
    ref_pdb: Optional[Path],
    allow_unverified_ref_order: bool,
) -> None:
    time_start = time.perf_counter()
    if allow_unverified_ref_order and ref_pdb is None:
        raise click.UsageError("--allow-unverified-ref-order requires --ref-pdb.")
    mode_resolved = _resolve_mode(mode, input_path)
    from mlmm.io.oniom_identity import extract_embedded_order_digest

    try:
        expected_ref_order_digest = extract_embedded_order_digest(input_path)
    except ValueError as exc:
        raise click.ClickException(str(exc)) from exc

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

    with tempfile.TemporaryDirectory(prefix="mlmm-oniom-import-") as staging_dir:
        staging_root = Path(staging_dir)
        staged_xyz = staging_root / xyz_path.name
        staged_pdb = staging_root / pdb_path.name
        _write_xyz(
            staged_xyz,
            coords,
            elements,
            comment=(
                f"mode={mode_resolved} atoms={n_atoms} qm={len(qm_indices)} "
                f"movable={len(movable_indices)} q={qm_charge} m={qm_mult}"
            ),
        )

        if ref_pdb is not None:
            ref_verification = _write_layered_pdb_with_ref(
                staged_pdb,
                ref_pdb,
                coords,
                elements,
                qm_indices,
                movable_indices,
                expected_ref_order_digest=expected_ref_order_digest,
                allow_unverified_ref_order=allow_unverified_ref_order,
            )
        else:
            _write_layered_pdb_without_ref(
                staged_pdb, coords, elements, qm_indices, movable_indices
            )

        commit_payloads(
            xyz_path,
            {
                xyz_path: staged_xyz.read_bytes(),
                pdb_path: staged_pdb.read_bytes(),
            },
        )

    click.echo(f"[oniom-import] mode={mode_resolved}")
    if ref_pdb is not None:
        click.echo(f"[oniom-import] ref_order={ref_verification}")
    click.echo(
        f"[oniom-import] atoms={n_atoms}, qm={len(qm_indices)}, movable={len(movable_indices)}, "
        f"frozen={n_atoms - len(set(movable_indices) | set(qm_indices))}"
    )
    click.echo(f"[oniom-import] wrote: {xyz_path}")
    click.echo(f"[oniom-import] wrote: {pdb_path}")
    from mlmm.core.output import emit
    from mlmm.core.utils import format_elapsed

    emit(
        format_elapsed("[time] Elapsed Time for ONIOM Import", time_start),
        narrative=True,
    )
