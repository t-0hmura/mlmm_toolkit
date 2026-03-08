# mlmm_toolkit/utils.py

"""
utils — concise utilities for configuration, plotting, and coordinates
====================================================================

Usage (API)
-----
    from mlmm_toolkit.utils import (
        build_energy_diagram,
        convert_xyz_to_pdb,
        merge_freeze_atom_indices,
        pretty_block,
    )

Examples::
    >>> from pathlib import Path
    >>> block = pretty_block("Geometry", {"freeze_atoms": [0, 1, 5]})
    >>> diagram = build_energy_diagram([0.0, 12.3, 5.4], ["R", "TS", "P"])

Description
-----
- **Generic helpers**
  - `pretty_block(title, content)`: Return a YAML-formatted block with an underlined title. Uses `yaml.safe_dump` with `allow_unicode=True`, `sort_keys=False`. Renders `{}` when `content` is empty.
  - `format_freeze_atoms_for_echo(cfg, key="freeze_atoms")`: Normalize geometry configuration for CLI echo. If the key is an iterable (but not a string), summarize to a compact single-line form like `"11073 atoms [0,1,2,3,4,...,12981,12982,12983,12984,12985]"`.
  - `format_elapsed(prefix, start_time, end_time=None)`: Format a wall-clock duration (HH:MM:SS.sss) given a start time and optional end time, using `time.perf_counter()` when the end time is omitted.
  - `merge_freeze_atom_indices(geom_cfg, *indices)`: Merge one or more iterables of atom indices into `geom_cfg["freeze_atoms"]`. Preserve existing entries, de-duplicate, sort numerically, and return the updated list (in place).
  - `apply_layer_freeze_constraints(geom_cfg, calc_cfg, layer_info, echo_fn=None)`: Merge layer-detected frozen indices (`layer_info["frozen_indices"]`) into both `geom_cfg["freeze_atoms"]` and `calc_cfg["freeze_atoms"]`, then optionally emit a concise summary line.
  - `deep_update(dst, src)`: Recursively update mapping `dst` with `src`. Nested dicts are merged, non-dicts overwrite; returns `dst`.
  - `_get_mapping_section(cfg, path)`: Internal helper to resolve a nested mapping section. Returns a `dict` or `None`.
  - `apply_yaml_overrides(yaml_cfg, overrides)`: For each target dictionary and its candidate key paths, find the first existing path in `yaml_cfg` and apply it via `deep_update`. Centralizes repeated `yaml_cfg.get(...)`-style merging.
  - `load_yaml_dict(path)`: Load a YAML file whose root must be a mapping. Returns `{}` when `path` is `None`. Raises `ValueError` if the YAML root is not a mapping.

- **Plotly: Energy diagram builder**
  - `build_energy_diagram(energies, labels, ylabel="ΔE", baseline=False, showgrid=False)`:
    Render an energy diagram where each state is a thick horizontal segment and adjacent states are connected by dotted diagonals (right end of left state → left end of right state). Segment length shrinks as the number of states grows to keep gaps readable. X ticks are centered on states and labeled by `labels`. Optional dotted baseline at the first state’s energy; optional grid. Energies are plotted as provided (no unit conversion). Returns a `plotly.graph_objs.Figure`. Validates equal lengths for `energies`/`labels` and non-empty input.

- **Coordinate conversion utilities**
  - `convert_xyz_to_pdb(xyz_path, ref_pdb_path, out_pdb_path)`:
    Overlay coordinates from an XYZ file (single or multi-frame) onto the atom ordering/topology of a reference PDB and write to `out_pdb_path`. The first frame creates/overwrites; subsequent frames append using `MODEL`/`ENDMDL`. Implemented with ASE (`ase.io.read`/`write`). Raises `ValueError` if no frames are found in the XYZ.

Outputs (& Directory Layout)
-----
- This module does not create directories.
- Functions primarily return Python objects or mutate dictionaries in place.
- On-disk output occurs only when explicitly requested by the caller:
  - `convert_xyz_to_pdb` writes a PDB file to `out_pdb_path` (first frame create/overwrite; subsequent frames append with `MODEL`/`ENDMDL` blocks).
  - `build_energy_diagram` returns a Plotly `Figure`; it does not write files unless the caller saves/exports the figure.

Notes:
-----
- Energy units in `build_energy_diagram` are passed through unchanged; ensure consistent units across states.
- Axis/line styling in `build_energy_diagram` is fixed-width with automatic padding; segment length adapts to the number of states.
- `load_yaml_dict` uses `yaml.safe_load` and enforces a mapping at the YAML root; empty files yield `{}`.
- `apply_yaml_overrides` tries candidate key paths in order and applies only the first existing mapping section per target.
- Dependencies: PyYAML, ASE (`ase.io.read`/`write`), Plotly (graph objects).
"""

import ast
import logging
import math
import os
import re
import time
import tempfile
from collections.abc import Iterable as _Iterable, Mapping, Sequence as _Sequence
from dataclasses import dataclass
from numbers import Real, Integral
from pathlib import Path
from typing import Any, Callable, Dict, Optional, Sequence, List, Tuple

import click
import numpy as np
import yaml
from ase.io import read, write
import plotly.graph_objs as go

from pysisyphus.helpers import geom_loader
from pysisyphus.constants import ANG2BOHR

from .add_elem_info import guess_element

logger = logging.getLogger(__name__)

# =============================================================================
# Generic helpers
# =============================================================================


def ensure_dir(path: Path) -> None:
    """Create a directory (parents ok); noop if it already exists."""
    path.mkdir(parents=True, exist_ok=True)


def read_xyz_as_blocks(path: Path, *, strict: bool = False) -> List[List[str]]:
    """Read an XYZ-style trajectory into blocks of lines.

    When *strict* is True, malformed headers or truncated frames raise a ClickException.
    """
    try:
        lines = path.read_text(encoding="utf-8").splitlines()
    except Exception as e:
        import click
        raise click.ClickException(f"Failed to read {path}: {e}")

    blocks: List[List[str]] = []
    i = 0
    while i < len(lines):
        if not lines[i].strip():
            i += 1
            continue
        try:
            n_atoms = int(lines[i].strip().split()[0])
        except Exception:
            if strict:
                import click
                raise click.ClickException(f"[xyz] Malformed XYZ/TRJ header at line {i+1} of {path}")
            break
        end = i + n_atoms + 2
        if end > len(lines):
            if strict:
                import click
                raise click.ClickException(f"[xyz] Incomplete XYZ frame at line {i+1} of {path}")
            break
        blocks.append(lines[i:end])
        i = end
    return blocks


def parse_xyz_block(
    block: Sequence[str],
    *,
    path: Path,
    frame_idx: int,
) -> Tuple[List[str], "np.ndarray"]:
    """Parse a single XYZ frame block into (elements, coords_angstrom)."""
    import click

    if not block:
        raise click.ClickException(f"[xyz] Empty XYZ frame in {path}")
    try:
        nat = int(block[0].strip().split()[0])
    except Exception:
        raise click.ClickException(
            f"[xyz] Malformed XYZ/TRJ header in frame {frame_idx} of {path}"
        )
    if len(block) < 2 + nat:
        raise click.ClickException(
            f"[xyz] Incomplete XYZ frame {frame_idx} in {path} (expected {nat} atoms)."
        )
    elems: List[str] = []
    coords: List[List[float]] = []
    for k in range(nat):
        parts = block[2 + k].split()
        if len(parts) < 4:
            raise click.ClickException(
                f"[xyz] Malformed atom line in frame {frame_idx} of {path}"
            )
        elems.append(parts[0])
        coords.append([float(parts[1]), float(parts[2]), float(parts[3])])
    return elems, np.array(coords, dtype=float)


def xyz_blocks_first_last(
    blocks: Sequence[Sequence[str]],
    *,
    path: Path,
) -> Tuple[List[str], "np.ndarray", "np.ndarray"]:
    """Return (elements, first_coords_ang, last_coords_ang) from pre-parsed XYZ blocks."""
    import click

    if not blocks:
        raise click.ClickException(f"[xyz] No frames found in {path}")
    first_elems, first_coords = parse_xyz_block(blocks[0], path=path, frame_idx=1)
    last_elems, last_coords = parse_xyz_block(blocks[-1], path=path, frame_idx=len(blocks))
    if first_elems != last_elems:
        raise click.ClickException(f"[xyz] Element list changed across frames in {path}")
    return first_elems, first_coords, last_coords


def read_xyz_first_last(trj_path: Path) -> Tuple[List[str], "np.ndarray", "np.ndarray"]:
    """Lightweight XYZ trajectory reader: return (elements, first_coords[Å], last_coords[Å])."""
    blocks = read_xyz_as_blocks(trj_path, strict=True)
    return xyz_blocks_first_last(blocks, path=trj_path)


def close_matplotlib_figures() -> None:
    """Best-effort cleanup for matplotlib figures to avoid open-figure warnings."""
    try:
        import matplotlib.pyplot as plt
        plt.close("all")
    except Exception:
        pass


def distance_A_from_coords(coords_bohr: "np.ndarray", i: int, j: int) -> float:
    """Return interatomic distance in Å given coords in Bohr."""
    diff = coords_bohr[i] - coords_bohr[j]
    return float(np.linalg.norm(diff) / ANG2BOHR)


def distance_tag(value_A: float, *, digits: int = 2, pad: int = 3) -> str:
    """Format a distance in Å as a zero-padded integer tag (default: ×10^2)."""
    scale = 10 ** digits
    return f"{int(round(value_A * scale)):0{pad}d}"


def values_from_bounds(low: float, high: float, h: float) -> "np.ndarray":
    """Return evenly spaced values from low→high with step cap h (inclusive)."""
    if h <= 0.0:
        raise click.BadParameter("--max-step-size must be > 0.")
    delta = abs(high - low)
    if delta < 1e-12:
        return np.array([low], dtype=float)
    N = int(math.ceil(delta / h))
    return np.linspace(low, high, N + 1, dtype=float)


def geom_from_xyz_string(
    xyz_text: str,
    *,
    coord_type: str,
    freeze_atoms: Optional[Sequence[int]] = None,
) -> Any:
    """Load a pysisyphus Geometry from an XYZ text string (tempfile-backed)."""
    s = xyz_text if xyz_text.endswith("\n") else (xyz_text + "\n")
    freeze_atoms = list(freeze_atoms) if freeze_atoms is not None else []
    tmp = tempfile.NamedTemporaryFile("w+", suffix=".xyz", delete=False)
    try:
        tmp.write(s)
        tmp.flush()
        tmp.close()

        g = geom_loader(
            Path(tmp.name),
            coord_type=coord_type,
            freeze_atoms=freeze_atoms,
        )
        try:
            g.freeze_atoms = np.array(sorted(set(map(int, freeze_atoms))), dtype=int)
        except Exception:
            click.echo(
                "[geom] WARNING: Failed to attach freeze_atoms to geometry.",
                err=True,
            )
        return g
    finally:
        try:
            os.unlink(tmp.name)
        except Exception:
            logger.debug("Failed to unlink temp file %s", tmp.name, exc_info=True)


def append_xyz_trajectory(dst_path: Path, src_path: Path, *, reset: bool = False) -> bool:
    """Append an XYZ trajectory segment to a concatenated trajectory file."""
    if not src_path.exists():
        return False
    mode = "w" if reset else "a"
    with src_path.open("r", encoding="utf-8") as src, dst_path.open(mode, encoding="utf-8") as dst:
        while True:
            chunk = src.read(1024 * 1024)
            if not chunk:
                break
            dst.write(chunk)
    return True


def snapshot_geometry(geom: Any, *, coord_type_default: str) -> Any:
    """Create an independent pysisyphus Geometry snapshot from the given Geometry."""
    s = geom.as_xyz()
    return geom_from_xyz_string(
        s,
        coord_type=getattr(geom, "coord_type", coord_type_default),
        freeze_atoms=getattr(geom, "freeze_atoms", []),
    )


def unbiased_energy_hartree(geom, base_calc) -> float:
    """Evaluate UMA energy (Hartree) without harmonic bias."""
    coords_bohr = np.asarray(geom.coords)
    elems = getattr(geom, "atoms", None)
    if elems is None:
        return float("nan")
    try:
        return float(base_calc.get_energy(elems, coords_bohr)["energy"])
    except Exception:
        return float("nan")


def pretty_block(title: str, content: Dict[str, Any]) -> str:
    """
    Return a YAML-formatted block with an underlined title.
    """
    body = yaml.safe_dump(_to_yaml_safe(content), sort_keys=False, allow_unicode=True).strip()
    return f"{title}\n" + "-" * len(title) + "\n" + (body if body else "(empty)") + "\n"


def _to_yaml_safe(value: Any) -> Any:
    """Recursively convert NumPy values/containers into YAML-safe builtins."""
    if isinstance(value, np.generic):
        return value.item()
    if isinstance(value, np.ndarray):
        return [_to_yaml_safe(v) for v in value.tolist()]
    if isinstance(value, Mapping):
        out: Dict[Any, Any] = {}
        for k, v in value.items():
            nk = _to_yaml_safe(k)
            if isinstance(nk, (list, tuple, set, dict)):
                nk = str(nk)
            out[nk] = _to_yaml_safe(v)
        return out
    if isinstance(value, tuple):
        return [_to_yaml_safe(v) for v in value]
    if isinstance(value, list):
        return [_to_yaml_safe(v) for v in value]
    if isinstance(value, set):
        return [_to_yaml_safe(v) for v in sorted(value, key=lambda x: str(x))]
    return value


# Backend-specific key prefixes in MLMM_CALC_KW.
# Keys with these prefixes are only relevant when the corresponding backend is active.
_BACKEND_KEY_PREFIXES: Dict[str, tuple] = {
    "uma": ("uma_model", "uma_task_name"),
    "orb": ("orb_model", "orb_precision"),
    "mace": ("mace_model", "mace_dtype"),
    "aimnet2": ("aimnet2_model",),
}


def filter_calc_for_echo(calc_cfg: Dict[str, Any]) -> Dict[str, Any]:
    """Remove backend-specific keys that are irrelevant for the active backend.

    Also hides xTB/embedcharge keys when embedcharge is disabled.
    """
    cfg = dict(calc_cfg)
    active = cfg.get("backend", "uma")

    # Remove keys belonging to inactive ML backends
    for backend, keys in _BACKEND_KEY_PREFIXES.items():
        if backend != active:
            for k in keys:
                cfg.pop(k, None)

    # Hide xTB-specific keys when embedcharge is disabled
    if not cfg.get("embedcharge"):
        for k in list(cfg):
            if k.startswith("xtb_"):
                cfg.pop(k)
        cfg.pop("embedcharge_step", None)

    return cfg


def strip_inherited_keys(
    child_cfg: Dict[str, Any],
    base_cfg: Dict[str, Any],
    *,
    mode: str = "present",
) -> Dict[str, Any]:
    """Return child_cfg without inherited keys (for concise logs).

    Parameters
    ----------
    child_cfg : Dict[str, Any]
        The child configuration dictionary to trim.
    base_cfg : Dict[str, Any]
        The base configuration dictionary to compare against.
    mode : str
        - "present": Remove keys that exist in base_cfg regardless of value.
        - "same": Remove keys only when the value matches base_cfg.

    Returns
    -------
    Dict[str, Any]
        A new dictionary with inherited keys removed.
    """
    if mode not in {"present", "same"}:
        raise ValueError(f"Unknown strip_inherited_keys mode: {mode}")
    trimmed: Dict[str, Any] = {}
    for key, value in child_cfg.items():
        if key in base_cfg:
            if mode == "present":
                continue
            if base_cfg.get(key) == value:
                continue
        trimmed[key] = value
    return trimmed


def _summarize_atom_indices(items: Sequence[Any]) -> str:
    """Return a compact single-line summary for atom indices."""
    if not items:
        return ""

    count = len(items)
    if count <= 64:
        return f"{count} atoms [{','.join(map(str, items))}]"

    head = ",".join(map(str, items[:5]))
    tail = ",".join(map(str, items[-5:]))
    return f"{count} atoms [{head},...,{tail}]"


def format_freeze_atoms_for_echo(
    cfg: Dict[str, Any],
    *,
    key: str = "freeze_atoms",
) -> Dict[str, Any]:
    """
    Normalize freeze-atoms fields for concise CLI echo output.
    """
    g = dict(cfg)
    freeze_atoms = g.get(key)
    if freeze_atoms is None:
        return g

    if isinstance(freeze_atoms, str):
        return g

    try:
        items = list(freeze_atoms)
    except TypeError:
        return g

    g[key] = _summarize_atom_indices(items)
    return g


def format_elapsed(prefix: str, start_time: float, end_time: Optional[float] = None) -> str:
    """Return a formatted elapsed-time string with the provided ``prefix`` label."""
    finish = end_time if end_time is not None else time.perf_counter()
    elapsed = max(0.0, finish - start_time)
    hours, rem = divmod(elapsed, 3600)
    minutes, seconds = divmod(rem, 60)
    return f"{prefix}: {int(hours):02d}:{int(minutes):02d}:{seconds:06.3f}"


def normalize_freeze_atoms(raw: Any) -> List[int]:
    """Normalize freeze_atoms values (string/list/iterable) into a list of integers.

    Parameters
    ----------
    raw : Any
        Input value that can be a string (e.g., "1,2,3" or "1 2 3"),
        a list of integers, or any iterable of numeric values.

    Returns
    -------
    List[int]
        List of integer indices.

    Examples
    --------
    >>> normalize_freeze_atoms("1, 2, 3")
    [1, 2, 3]
    >>> normalize_freeze_atoms([1, 2, 3])
    [1, 2, 3]
    >>> normalize_freeze_atoms(None)
    []
    """
    import re

    if raw is None:
        return []
    if isinstance(raw, str):
        tokens = re.findall(r"-?\d+", raw)
        return [int(tok) for tok in tokens]
    try:
        return [int(i) for i in raw]
    except Exception:
        return []


def merge_freeze_atom_indices(
    geom_cfg: Dict[str, Any],
    *indices: _Iterable[int],
) -> List[int]:
    """Merge one or more iterables of indices into ``geom_cfg['freeze_atoms']``.

    Existing entries are preserved, duplicates removed, and the result sorted.
    The updated list is returned.
    """
    merged: set[int] = set()
    base = geom_cfg.get("freeze_atoms", None)
    merged.update(normalize_freeze_atoms(base))
    for seq in indices:
        merged.update(normalize_freeze_atoms(seq))
    result = sorted(merged)
    geom_cfg["freeze_atoms"] = result
    return result


# =============================================================================
# Link-freezing helpers
# =============================================================================


def parse_pdb_coords(pdb_path):
    """Parse ATOM/HETATM records from *pdb_path* and separate link hydrogen (HL) atoms.

    Returns:
        A tuple (others, lkhs) where:
            - others: list of tuples (index, x, y, z, line) for all atoms except the
              'HL' atom of residue 'LKH'. ``index`` is the 0-based position in the
              atom sequence as loaded from the *first* MODEL (or the full file if no
              MODEL records are present).
            - lkhs: list of tuples (x, y, z, line) for atoms where residue name is
              'LKH' and atom name is 'HL' in the same MODEL selection.

    Notes
    -----
        - Coordinates are read from standard PDB columns:
          X: columns 31-38, Y: 39-46, Z: 47-54 (1-based indexing).
        - If multiple MODEL blocks are present, only the first model is considered,
          matching typical geom_loader behavior.
    """
    with open(pdb_path, "r") as f:
        lines = f.readlines()

    others = []
    lkhs = []
    model_seen = False
    in_first_model = True
    atom_index = 0
    for line in lines:
        if line.startswith("MODEL"):
            if not model_seen:
                model_seen = True
                in_first_model = True
            else:
                in_first_model = False
            continue
        if line.startswith("ENDMDL"):
            if model_seen and in_first_model:
                break
            continue
        if model_seen and not in_first_model:
            continue
        if not (line.startswith("ATOM") or line.startswith("HETATM")):
            continue

        current_index = atom_index
        atom_index += 1

        name = line[12:16].strip()
        resname = line[17:20].strip()
        try:
            x = float(line[30:38])
            y = float(line[38:46])
            z = float(line[46:54])
        except ValueError:
            continue

        if resname == "LKH" and name == "HL":
            lkhs.append((x, y, z, line))
        else:
            others.append((current_index, x, y, z, line))
    return others, lkhs


def nearest_index(point, pool):
    """Find the nearest point in *pool* to *point* using Euclidean distance.

    Args:
        point: Tuple (x, y, z) representing the query coordinate.
        pool: Iterable of tuples (index, x, y, z, line) to search.

    Returns:
        A tuple (index, distance) where:
            - index is the 0-based index of the nearest entry in *pool* (or -1 if *pool* is empty).
            - distance is the Euclidean distance to that entry (``inf`` if *pool* is empty).
    """
    x, y, z = point
    best_i = -1
    best_d2 = float("inf")
    for atom_index, a, b, c, _ in pool:
        d2 = (a - x) ** 2 + (b - y) ** 2 + (c - z) ** 2
        if d2 < best_d2:
            best_d2 = d2
            best_i = atom_index
    return best_i, math.sqrt(best_d2)


def detect_freeze_links(pdb_path):
    """Identify link-parent atom indices for 'LKH'/'HL' link hydrogens.

    For each 'HL' atom in residue 'LKH', find the nearest atom among all other
    ATOM/HETATM records and return the indices of those nearest neighbors in the
    same atom ordering used by geometry loading (first MODEL if present).

    Args:
        pdb_path: Path to the input PDB file.

    Returns:
        List of 0-based indices into the full atom sequence (including any link H atoms)
        corresponding to the nearest neighbors (link parents). Returns an empty list if
        no LKH/HL atoms are present or if link hydrogens exist without any other atoms.
    """
    others, lkhs = parse_pdb_coords(pdb_path)

    if not lkhs or not others:
        return []

    indices = []
    for (x, y, z, line) in lkhs:
        idx, dist = nearest_index((x, y, z), others)
        if idx >= 0:
            indices.append(idx)
    return indices


def detect_freeze_links_logged(pdb_path: Path) -> List[int]:
    """Return link-parent indices and raise a user-facing error on failure."""
    try:
        return list(detect_freeze_links(pdb_path))
    except Exception as e:
        raise click.ClickException(
            f"[freeze-links] Failed to detect link parents for '{pdb_path.name}': {e}"
        ) from e


def merge_detected_freeze_links(
    geom_cfg: Dict[str, Any],
    pdb_path: Path,
    *,
    prefix: str = "[freeze-links]",
) -> List[int]:
    """Detect link-parent atoms and merge them into ``geom_cfg['freeze_atoms']``."""
    detected = detect_freeze_links_logged(pdb_path)
    merged = merge_freeze_atom_indices(geom_cfg, detected)
    if merged:
        click.echo(f"{prefix} Freeze atoms (0-based): {','.join(map(str, merged))}")
    return merged


def apply_layer_freeze_constraints(
    geom_cfg: Dict[str, Any],
    calc_cfg: Dict[str, Any],
    layer_info: Optional[Dict[str, Sequence[int]]],
    *,
    echo_fn: Optional[Callable[[str], None]] = None,
) -> List[int]:
    """Merge frozen-layer atoms into geometry/calculator freeze lists."""
    base_freeze = normalize_freeze_atoms(geom_cfg.get("freeze_atoms"))
    frozen_from_layer = normalize_freeze_atoms((layer_info or {}).get("frozen_indices", []))

    if frozen_from_layer:
        before = set(base_freeze)
        merged = sorted(before | set(frozen_from_layer))
        added = len(set(merged) - before)
        if echo_fn is not None:
            echo_fn(
                f"[layer] Applied freeze constraints from frozen layer: "
                f"total={len(merged)} (added_from_layer={added})"
            )
    else:
        merged = sorted(set(base_freeze))

    geom_cfg["freeze_atoms"] = merged
    calc_cfg["freeze_atoms"] = merged
    return merged


def deep_update(dst: Dict[str, Any], src: Optional[Dict[str, Any]]) -> Dict[str, Any]:
    """
    Recursively update mapping *dst* with *src*, returning *dst*.
    """
    for k, v in (src or {}).items():
        if isinstance(v, dict) and isinstance(dst.get(k), dict):
            deep_update(dst[k], v)
        else:
            dst[k] = v
    return dst


def collect_single_option_values(
    argv: _Sequence[str],
    names: _Sequence[str],
    label: str,
) -> List[str]:
    """Collect values following a flag that must appear at most once."""
    vals: List[str] = []
    seen = 0
    i = 0
    while i < len(argv):
        tok = argv[i]
        if tok in names:
            seen += 1
            j = i + 1
            while j < len(argv) and not argv[j].startswith("-"):
                vals.append(argv[j])
                j += 1
            i = j
        else:
            i += 1
    if seen > 1:
        raise click.BadParameter(
            f"Use a single {label} followed by multiple values; repeated flags are not accepted."
        )
    return vals


def load_pdb_atom_metadata(pdb_path: Path) -> List[Dict[str, Any]]:
    """Return per-atom metadata (serial, name, resname, resseq, element) in file order."""
    atoms: List[Dict[str, Any]] = []
    with open(pdb_path, "r") as f:
        for line in f:
            if not (line.startswith("ATOM") or line.startswith("HETATM")):
                continue

            serial_txt = line[6:11].strip()
            resseq_txt = line[22:26].strip()
            atom_name = line[12:16].strip()
            res_name = line[17:20].strip()
            element_txt = line[76:78].strip()
            is_hetatm = line.startswith("HETATM")

            try:
                serial = int(serial_txt) if serial_txt else None
            except ValueError:
                serial = None
            try:
                resseq = int(resseq_txt) if resseq_txt else None
            except ValueError:
                resseq = None

            if not element_txt:
                inferred = guess_element(atom_name, res_name, is_hetatm)
                element_txt = inferred or ""

            atoms.append(
                {
                    "serial": serial,
                    "name": atom_name,
                    "resname": res_name,
                    "resseq": resseq,
                    "element": element_txt,
                }
            )
    return atoms


def resolve_atom_spec_index(spec: str, atom_meta: _Sequence[Dict[str, Any]]) -> int:
    """Resolve an atom selector string into a 0-based atom index using PDB metadata."""
    tokens = [t for t in re.split(r"[\s/`,\\]+", spec.strip().replace(" ", ",")) if t]
    if len(tokens) != 3:
        raise ValueError(
            f"Atom spec '{spec}' must have exactly 3 fields (resname, resseq, atomname)."
        )

    tokens_upper = [t.upper() for t in tokens]
    matches: List[int] = []
    for idx, meta in enumerate(atom_meta):
        resname = (meta.get("resname") or "").strip().upper()
        resseq = meta.get("resseq")
        atom = (meta.get("name") or "").strip().upper()
        if resseq is None:
            continue
        fields = {resname, str(resseq), atom}
        if all(tok in fields for tok in tokens_upper):
            matches.append(idx)

    if len(matches) == 1:
        return matches[0]
    if len(matches) > 1:
        raise ValueError(
            f"Atom spec '{spec}' matches {len(matches)} atoms; use an explicit atom index."
        )

    resname, resseq_str, atom = tokens_upper
    if not resseq_str.isdigit():
        raise ValueError(
            f"Atom spec '{spec}' could not be resolved and residue number '{tokens[1]}' is not numeric."
        )
    resseq_int = int(resseq_str)
    ordered_matches = [
        idx
        for idx, meta in enumerate(atom_meta)
        if (meta.get("resname") or "").strip().upper() == resname
        and meta.get("resseq") == resseq_int
        and (meta.get("name") or "").strip().upper() == atom
    ]
    if len(ordered_matches) == 1:
        return ordered_matches[0]
    if len(ordered_matches) > 1:
        raise ValueError(
            f"Atom spec '{spec}' matches {len(ordered_matches)} atoms after ordered fallback; "
            "use an explicit atom index."
        )

    raise ValueError(f"Atom spec '{spec}' did not match any atom.")


def atom_label_from_meta(atom_meta: _Sequence[Dict[str, Any]], index: int) -> str:
    if index < 0 or index >= len(atom_meta):
        return f"idx{index}"
    meta = atom_meta[index]
    resname = (meta.get("resname") or "?").strip() or "?"
    resseq = meta.get("resseq")
    resseq_txt = "?" if resseq is None else str(resseq)
    atom = (meta.get("name") or "?").strip() or "?"
    return f"{resname}-{resseq_txt}-{atom}"


def axis_label_csv(
    axis_name: str,
    i_idx: int,
    j_idx: int,
    one_based: bool,
    atom_meta: Optional[_Sequence[Dict[str, Any]]] = None,
    pair_raw: Optional[Tuple[Any, Any, float, float]] = None,
) -> str:
    if pair_raw and (isinstance(pair_raw[0], str) or isinstance(pair_raw[1], str)) and atom_meta:
        i_label = atom_label_from_meta(atom_meta, i_idx)
        j_label = atom_label_from_meta(atom_meta, j_idx)
        return f"{axis_name}_{i_label}_{j_label}_A"
    i_disp = i_idx + 1 if one_based else i_idx
    j_disp = j_idx + 1 if one_based else j_idx
    return f"{axis_name}_{i_disp}_{j_disp}_A"


def axis_label_html(label: str) -> str:
    parts = label.split("_")
    if len(parts) >= 4 and parts[-1] == "A":
        axis = parts[0]
        i_disp = parts[1]
        j_disp = parts[2]
        return f"{axis} ({i_disp},{j_disp}) (Å)"
    return label


def resolve_scan_index(
    value: Any,
    *,
    one_based: bool,
    atom_meta: Optional[_Sequence[Dict[str, Any]]],
    context: str,
) -> int:
    """Resolve an index or atom-spec string for scan lists with consistent errors."""
    if isinstance(value, Integral):
        idx_val = int(value)
        if one_based:
            idx_val -= 1
        if idx_val < 0:
            raise click.BadParameter(
                f"Negative atom index after base conversion in {context}: {idx_val} (0-based expected)."
            )
        return idx_val
    if isinstance(value, str):
        if not atom_meta:
            raise click.BadParameter(
                f"{context} uses a string atom spec, but no PDB metadata is available."
            )
        try:
            return resolve_atom_spec_index(value, atom_meta)
        except ValueError as exc:
            raise click.BadParameter(f"{context} {exc}")
    raise click.BadParameter(f"{context} must be an int index or atom spec string.")


def parse_scan_list_triples(
    raw: str,
    *,
    one_based: bool,
    atom_meta: Optional[_Sequence[Dict[str, Any]]],
    option_name: str,
    return_one_based: bool = False,
) -> Tuple[List[Tuple[int, int, float]], List[Tuple[Any, Any, float]]]:
    """Parse --scan-lists triples into indices (0-based by default)."""
    try:
        obj = ast.literal_eval(raw)
    except Exception as e:
        raise click.BadParameter(f"Invalid literal for {option_name}: {e}")

    if not isinstance(obj, (list, tuple)):
        raise click.BadParameter(f"{option_name} must be a list/tuple of (i,j,target).")

    parsed: List[Tuple[int, int, float]] = []
    for entry_idx, t in enumerate(obj, start=1):
        if not (
            isinstance(t, (list, tuple))
            and len(t) == 3
            and isinstance(t[2], Real)
        ):
            raise click.BadParameter(
                f"{option_name} entry {entry_idx} must be (i,j,target): got {t}"
            )

        i = resolve_scan_index(
            t[0],
            one_based=one_based,
            atom_meta=atom_meta,
            context=f"{option_name} entry {entry_idx} (i)",
        )
        j = resolve_scan_index(
            t[1],
            one_based=one_based,
            atom_meta=atom_meta,
            context=f"{option_name} entry {entry_idx} (j)",
        )
        if return_one_based:
            i += 1
            j += 1
        parsed.append((i, j, float(t[2])))

    return parsed, list(obj)


def parse_scan_list_quads(
    raw: str,
    *,
    expected_len: int,
    one_based: bool,
    atom_meta: Optional[_Sequence[Dict[str, Any]]],
    option_name: str,
) -> Tuple[List[Tuple[int, int, float, float]], List[Tuple[Any, Any, float, float]]]:
    """Parse --scan-lists quadruples into 0-based indices."""
    try:
        obj = ast.literal_eval(raw)
    except Exception as e:
        raise click.BadParameter(f"Invalid literal for {option_name}: {e}")

    if not (isinstance(obj, (list, tuple)) and len(obj) == expected_len):
        quads = ",".join([f"(i{n},j{n},low{n},high{n})" for n in range(1, expected_len + 1)])
        raise click.BadParameter(
            f"{option_name} must contain exactly {expected_len} quadruples: [{quads}]"
        )

    parsed: List[Tuple[int, int, float, float]] = []
    for entry_idx, q in enumerate(obj, start=1):
        if not (
            isinstance(q, (list, tuple))
            and len(q) == 4
            and isinstance(q[2], Real)
            and isinstance(q[3], Real)
        ):
            raise click.BadParameter(f"{option_name} entry must be (i,j,low,high): got {q}")

        i = resolve_scan_index(
            q[0],
            one_based=one_based,
            atom_meta=atom_meta,
            context=f"{option_name} entry {entry_idx} (i)",
        )
        j = resolve_scan_index(
            q[1],
            one_based=one_based,
            atom_meta=atom_meta,
            context=f"{option_name} entry {entry_idx} (j)",
        )
        parsed.append((i, j, float(q[2]), float(q[3])))

    for i, j, low, high in parsed:
        if low <= 0.0 or high <= 0.0:
            raise click.BadParameter(f"Distances must be positive: {(i, j, low, high)}")

    return parsed, list(obj)


def _load_scan_spec_root(
    spec_path: Path,
    *,
    option_name: str = "--spec",
) -> Mapping[str, Any]:
    """Load a scan spec (YAML/JSON) and ensure mapping root."""
    try:
        with open(spec_path, "r", encoding="utf-8") as handle:
            data = yaml.safe_load(handle)
    except Exception as exc:
        raise click.BadParameter(
            f"Failed to parse {option_name} file '{spec_path}': {exc}"
        )

    if data is None:
        raise click.BadParameter(f"{option_name} file '{spec_path}' is empty.")
    if not isinstance(data, Mapping):
        raise click.BadParameter(
            f"{option_name} file '{spec_path}' must have a mapping at the YAML/JSON root."
        )
    return data


def _spec_one_based(
    value: Any,
    *,
    default: bool,
    option_name: str = "--spec",
) -> bool:
    """Resolve one_based value from spec with CLI fallback."""
    if value is None:
        return bool(default)
    if isinstance(value, bool):
        return value
    if isinstance(value, str):
        key = value.strip().lower()
        if key in {"1", "true", "yes", "y", "on"}:
            return True
        if key in {"0", "false", "no", "n", "off"}:
            return False
    raise click.BadParameter(
        f"{option_name} field 'one_based' must be a boolean (true/false)."
    )


def _first_spec_field(
    spec_cfg: Mapping[str, Any],
    candidates: _Sequence[str],
) -> Tuple[Optional[str], Any]:
    for key in candidates:
        if key in spec_cfg:
            return key, spec_cfg[key]
    return None, None


def parse_scan_spec_stages(
    spec_path: Path,
    *,
    one_based_default: bool,
    atom_meta: Optional[_Sequence[Dict[str, Any]]],
    option_name: str = "--spec",
) -> Tuple[List[List[Tuple[int, int, float]]], bool]:
    """Parse staged 1D scan spec into 0-based stage triples."""
    spec_cfg = _load_scan_spec_root(spec_path, option_name=option_name)
    stages_key, stages_raw = _first_spec_field(spec_cfg, ("stages",))
    if stages_key is None:
        raise click.BadParameter(f"{option_name} must define 'stages'.")
    if not isinstance(stages_raw, (list, tuple)) or len(stages_raw) == 0:
        raise click.BadParameter(f"{option_name} field '{stages_key}' must be a non-empty list.")

    one_based = _spec_one_based(
        spec_cfg.get("one_based"), default=one_based_default, option_name=option_name
    )
    stages: List[List[Tuple[int, int, float]]] = []
    for stage_idx, stage_raw in enumerate(stages_raw, start=1):
        if not isinstance(stage_raw, (list, tuple)):
            raise click.BadParameter(
                f"{option_name} {stages_key}[{stage_idx}] must be a list of (i,j,target) entries."
            )
        parsed, _ = parse_scan_list_triples(
            repr(list(stage_raw)),
            one_based=one_based,
            atom_meta=atom_meta,
            option_name=f"{option_name} {stages_key}[{stage_idx}]",
        )
        if not parsed:
            raise click.BadParameter(
                f"{option_name} {stages_key}[{stage_idx}] must contain at least one (i,j,target) triple."
            )
        for i, j, target in parsed:
            if target <= 0.0:
                raise click.BadParameter(
                    f"Non-positive target distance in {option_name} {stages_key}[{stage_idx}]: {(i, j, target)}."
                )
        stages.append(parsed)
    return stages, one_based


def parse_scan_spec_quads(
    spec_path: Path,
    *,
    expected_len: int,
    one_based_default: bool,
    atom_meta: Optional[_Sequence[Dict[str, Any]]],
    option_name: str = "--spec",
) -> Tuple[List[Tuple[int, int, float, float]], List[Tuple[Any, Any, float, float]], bool]:
    """Parse 2D/3D scan spec into 0-based quad tuples."""
    spec_cfg = _load_scan_spec_root(spec_path, option_name=option_name)
    pairs_key, pairs_raw = _first_spec_field(spec_cfg, ("pairs",))
    if pairs_key is None:
        raise click.BadParameter(f"{option_name} must define 'pairs'.")
    if not isinstance(pairs_raw, (list, tuple)):
        raise click.BadParameter(f"{option_name} field '{pairs_key}' must be a list.")

    one_based = _spec_one_based(
        spec_cfg.get("one_based"), default=one_based_default, option_name=option_name
    )
    parsed, raw_pairs = parse_scan_list_quads(
        repr(list(pairs_raw)),
        expected_len=expected_len,
        one_based=one_based,
        atom_meta=atom_meta,
        option_name=f"{option_name} {pairs_key}",
    )
    return parsed, raw_pairs, one_based


PDB_ATOM_META_HEADER = f"{'id':>5} {'atom':<4} {'res':<4} {'resid':>4} {'el':<2}"


def format_pdb_atom_metadata(atom_meta: _Sequence[Dict[str, Any]], index: int) -> str:
    """Format metadata for atom *index* as aligned text: serial name resname resseq element."""
    fallback_serial = index + 1
    if index < 0 or index >= len(atom_meta):
        return f"{fallback_serial:>5} {'?':<4} {'?':<4} {'?':>4} {'?':<2}"

    meta = atom_meta[index]
    serial = meta.get("serial") or fallback_serial
    name = meta.get("name") or "?"
    resname = meta.get("resname") or "?"
    resseq = meta.get("resseq")
    resseq_str = "?" if resseq is None else str(resseq)
    element = (meta.get("element") or "?").strip() or "?"

    return f"{serial:>5} {name:<4} {resname:<4} {resseq_str:>4} {element:<2}"


def normalize_choice(
    value: str,
    *,
    param: str,
    alias_groups: Sequence[Tuple[Sequence[str], str]],
    allowed_hint: str,
) -> str:
    """Normalize a mode choice using alias groups and raise error on failure.

    Parameters
    ----------
    value : str
        The value to normalize.
    param : str
        Parameter name for error messages.
    alias_groups : Sequence[Tuple[Sequence[str], str]]
        Sequence of (aliases, canonical) pairs where aliases is a sequence of strings.
    allowed_hint : str
        Description of allowed values for error messages.

    Returns
    -------
    str
        The canonical value corresponding to the matched alias.

    Raises
    ------
    click.BadParameter
        If the value does not match any alias.
    """
    key = (value or "").strip().lower()
    for aliases, canonical in alias_groups:
        if any(key == alias.lower() for alias in aliases):
            return canonical

    hint = allowed_hint.strip()
    detail = f" Allowed: {hint}." if hint else ""
    raise click.BadParameter(f"Unknown value for {param} '{value}'.{detail}")


def _get_mapping_section(cfg: Mapping[str, Any], path: _Sequence[str]) -> Optional[Dict[str, Any]]:
    cur: Any = cfg
    for key in path:
        if not isinstance(cur, Mapping):
            return None
        cur = cur.get(key)
        if cur is None:
            return None
    return cur if isinstance(cur, dict) else None


def apply_yaml_overrides(
    yaml_cfg: Mapping[str, Any],
    overrides: _Sequence[Tuple[Dict[str, Any], _Sequence[_Sequence[str]]]],
) -> None:
    """Apply YAML overrides to multiple target dictionaries.

    Parameters
    ----------
    yaml_cfg : Mapping[str, Any]
        Parsed YAML configuration (root-level mapping).
    overrides : Sequence[Tuple[Dict[str, Any], Sequence[Sequence[str]]]]
        Each entry consists of the target dictionary to update followed by one or
        more candidate key paths. The first existing path is used. For example::

            apply_yaml_overrides(
                yaml_cfg,
                [
                    (geom_cfg, (("geom",),)),
                    (lbfgs_cfg, (("stopt", "lbfgs"), ("lbfgs",))),
                ],
            )

        This mirrors the previous ``deep_update(..., yaml_cfg.get(...))`` pattern
        while centralizing the shared logic.
    """
    for target, paths in overrides:
        for path in paths:
            norm_path = tuple(path)
            section = _get_mapping_section(yaml_cfg, norm_path)
            if section is not None:
                deep_update(target, section)
                break


def yaml_section_has_key(
    yaml_cfg: Mapping[str, Any],
    paths: _Sequence[_Sequence[str]],
    key: str,
) -> bool:
    """Return True when any candidate YAML section explicitly defines ``key``."""
    for path in paths:
        section = _get_mapping_section(yaml_cfg, tuple(path))
        if isinstance(section, Mapping) and (key in section):
            return True
    return False


def load_yaml_dict(path: Optional[Path]) -> Dict[str, Any]:
    """
    Load a YAML file whose root must be a mapping. Return an empty dict if *path* is None.
    """
    if not path:
        return {}

    with open(path, "r") as f:
        data = yaml.safe_load(f) or {}

    if not isinstance(data, dict):
        raise ValueError(f"YAML root must be a mapping, got: {type(data)}")

    return data


# =============================================================================
# Plotly: Energy diagram builder
# =============================================================================
def build_energy_diagram(
    energies: Sequence[float],
    labels: Sequence[str],
    ylabel: str = "ΔE",
    baseline: bool = False,
    showgrid: bool = False,
) -> go.Figure:
    """
    Plot an energy diagram using Plotly.

    Parameters
    ----------
    energies : Sequence[float]
        Energies for each state (same unit). Values are plotted without conversion.
    labels : Sequence[str]
        Labels corresponding to each state (for example, ["R", "TS1", "IM1", "TS2", "P"]).
        Must be the same length as ``energies``.
    ylabel : str, optional
        Y-axis label (for example, "ΔE" or "ΔG"). Defaults to ``"ΔE"``.
    baseline : bool, optional
        If ``True``, draw a dotted baseline at the energy of the first state across the plot.
    showgrid : bool, optional
        If ``True``, show grid lines on both axes. Defaults to ``False``.

    Returns
    -------
    plotly.graph_objs.Figure
        Figure containing the energy diagram.

    Notes
    -----
    - Each state is rendered as a thick horizontal segment (width ``HLINE_WIDTH``).
    - Adjacent states are connected by dotted diagonal segments from the right end of
      the left state to the left end of the right state.
    - Segment length automatically shrinks with additional states so that gaps remain
      between neighbors.
    - X-axis ticks are centered on each state and labeled using ``labels``.
    """
    if len(energies) == 0:
        raise ValueError("`energies` must contain at least one value.")
    if len(energies) != len(labels):
        raise ValueError("`energies` and `labels` must have the same length.")

    n = len(energies)
    energies = [float(e) for e in energies]

    # -----------------------------
    # Layout/style constants
    # -----------------------------
    AXIS_WIDTH = 3
    FONT_SIZE = 18
    AXIS_TITLE_SIZE = 20
    HLINE_WIDTH = 6           # Width of the horizontal state segments
    CONNECTOR_WIDTH = 2       # Width of the dotted connectors
    LINE_COLOR = "#1C1C1C"
    GRID_COLOR = "lightgrey"

    # -----------------------------
    # Geometry along the X axis (centers and segment lengths)
    # -----------------------------
    # Place segment centers at 0.5, 1.5, 2.5, ... (equally spaced)
    centers = [i + 0.5 for i in range(n)]

    # Shorten the segment as n grows (min 0.35, max 0.85)
    # Examples: n=5 -> 0.7, n=10 -> 0.5, n>=20 -> 0.35
    seg_width = min(0.85, max(0.35, 0.90 - 0.04 * n))
    half = seg_width / 2.0

    lefts = [c - half for c in centers]
    rights = [c + half for c in centers]

    # -----------------------------
    # Assemble the figure
    # -----------------------------
    fig = go.Figure()

    # Baseline (dotted line at the first energy level)
    if baseline:
        fig.add_trace(
            go.Scatter(
                x=[lefts[0], rights[-1]],
                y=[energies[0], energies[0]],
                mode="lines",
                line=dict(color=GRID_COLOR, dash="dot", width=2),
                hoverinfo="skip",
                showlegend=False,
            )
        )

    # Horizontal segments for each state
    for i, (e, lab) in enumerate(zip(energies, labels)):
        fig.add_trace(
            go.Scatter(
                x=[lefts[i], rights[i]],
                y=[e, e],
                mode="lines",
                line=dict(color=LINE_COLOR, width=HLINE_WIDTH),
                hovertemplate=f"{lab}: %{{y:.6f}}<extra></extra>",
                showlegend=False,
            )
        )

    # Dotted diagonals between adjacent states (right end -> left end)
    for i in range(n - 1):
        fig.add_trace(
            go.Scatter(
                x=[rights[i], lefts[i + 1]],
                y=[energies[i], energies[i + 1]],
                mode="lines",
                line=dict(color=LINE_COLOR, width=CONNECTOR_WIDTH, dash="dot"),
                hoverinfo="skip",
                showlegend=False,
            )
        )

    # -----------------------------
    # Axis ranges and styling
    # -----------------------------
    # Add a small margin beyond the first/last segments on X
    xpad = max(0.08, 0.15 * (1.0 - seg_width))
    x_min = lefts[0] - xpad
    x_max = rights[-1] + xpad

    # Add vertical padding above and below
    y_min = min(energies)
    y_max = max(energies)
    span = max(1e-6, y_max - y_min)  # Avoid zero span even if all values match
    ypad_low = 0.10 * span
    ypad_high = 0.20 * span
    y_range = [y_min - ypad_low, y_max + ypad_high]

    xaxis_config = dict(
        range=[x_min, x_max],
        showline=True,
        linewidth=AXIS_WIDTH,
        linecolor=LINE_COLOR,
        mirror=True,
        ticks="inside",
        tickwidth=AXIS_WIDTH,
        tickcolor=LINE_COLOR,
        tickfont=dict(size=FONT_SIZE, color=LINE_COLOR),
        showgrid=showgrid,
        gridcolor=GRID_COLOR,
        gridwidth=0.5,
        zeroline=False,
        tickmode="array",
        tickvals=centers,
        ticktext=list(labels),
        title=dict(text="", font=dict(size=AXIS_TITLE_SIZE, color=LINE_COLOR)),
    )

    yaxis_config = dict(
        range=y_range,
        showline=True,
        linewidth=AXIS_WIDTH,
        linecolor=LINE_COLOR,
        mirror=True,
        ticks="inside",
        tickwidth=AXIS_WIDTH,
        tickcolor=LINE_COLOR,
        tickfont=dict(size=FONT_SIZE, color=LINE_COLOR),
        showgrid=showgrid,
        gridcolor=GRID_COLOR,
        gridwidth=0.5,
        zeroline=False,
        title=dict(text=ylabel, font=dict(size=AXIS_TITLE_SIZE, color=LINE_COLOR)),
    )

    fig.update_layout(
        xaxis=xaxis_config,
        yaxis=yaxis_config,
        plot_bgcolor="white",
        paper_bgcolor="white",
        margin=dict(l=80, r=40, t=40, b=80),
    )

    return fig


# =============================================================================
# Coordinate conversion utilities
# =============================================================================
def convert_xyz_to_pdb(xyz_path: Path, ref_pdb_path: Path, out_pdb_path: Path) -> None:
    """Overlay coordinates from *xyz_path* onto the topology of *ref_pdb_path* and write to *out_pdb_path*.

    Notes:
        - *xyz_path* may contain one or many frames. For multi‑frame trajectories,
          a MODEL/ENDMDL block is appended for each subsequent frame in the output PDB.
        - On the first frame the output file is created/overwritten; subsequent frames are appended.
        - Validates that atom ordering (element symbols) matches between XYZ and PDB.

    Args:
        xyz_path: Path to an XYZ file (single or multi-frame).
        ref_pdb_path: Path to a reference PDB providing atom ordering/topology.
        out_pdb_path: Destination PDB file to write.

    Raises:
        ValueError: If no frames found in XYZ file or atom ordering mismatch.
    """
    ref_atoms = read(ref_pdb_path)  # Reference topology/ordering (single frame)
    traj = read(xyz_path, index=":", format="xyz")  # Load all frames from the XYZ
    if not traj:
        raise ValueError(f"No frames found in {xyz_path}.")

    ref_symbols = ref_atoms.get_chemical_symbols()

    for step, frame in enumerate(traj):
        xyz_symbols = frame.get_chemical_symbols()
        xyz_positions = frame.get_positions()

        if xyz_symbols != ref_symbols:
            # If atom counts match, the PDB likely has missing/wrong element columns.
            # Trust XYZ symbols (from pysisyphus which uses proper element detection)
            # and patch the reference atoms so PDB topology is preserved.
            if len(xyz_symbols) == len(ref_symbols):
                ref_atoms.set_chemical_symbols(xyz_symbols)
                ref_symbols = xyz_symbols
            else:
                raise ValueError(
                    "Atom count mismatch between XYZ and PDB; "
                    f"XYZ has {len(xyz_symbols)} atoms, PDB has {len(ref_symbols)}."
                )

        atoms = ref_atoms.copy()
        atoms.set_positions(xyz_positions)
        if step == 0:
            write(out_pdb_path, atoms)  # Create/overwrite on the first frame
        else:
            write(out_pdb_path, atoms, append=True)  # Append subsequent frames using MODEL/ENDMDL


# =============================================================================
# Global toggle for XYZ/TRJ → PDB conversion
# =============================================================================
_CONVERT_FILES_ENABLED: bool = True


def set_convert_file_enabled(enabled: bool) -> None:
    """Globally enable or disable XYZ/TRJ conversions to PDB outputs."""
    global _CONVERT_FILES_ENABLED
    _CONVERT_FILES_ENABLED = bool(enabled)


def is_convert_file_enabled() -> bool:
    """Check if convert-files is globally enabled."""
    return _CONVERT_FILES_ENABLED


def convert_xyz_like_outputs(
    xyz_path: Path,
    ref_pdb_path: Optional[Path],
    out_pdb_path: Optional[Path] = None,
    *,
    context: str = "outputs",
    on_error: str = "raise",
) -> bool:
    """Convert an XYZ file to PDB output using ref topology.

    Respects the global _CONVERT_FILES_ENABLED toggle.
    Returns True when conversion succeeded; False otherwise.
    """
    if not _CONVERT_FILES_ENABLED:
        return False
    if ref_pdb_path is None or out_pdb_path is None:
        return False
    try:
        convert_xyz_to_pdb(xyz_path, ref_pdb_path, out_pdb_path)
        return True
    except Exception as e:
        if on_error == "warn":
            click.echo(f"[convert] WARNING: Failed to convert {context}: {e}", err=True)
            return False
        raise click.ClickException(f"[convert] Failed to convert {context}: {e}") from e


def pdb_keys_from_line(line: str) -> Tuple[Tuple, Tuple]:
    """Extract robust keys from a PDB ATOM/HETATM record.

    Returns:
        key_full: (chain, resseq, icode, resname, atomname, altloc)
        key_simple: (chain, resseq, icode, atomname)
    """
    atomname = line[12:16].strip()
    altloc = line[16:17].strip()
    resname = line[17:20].strip()
    chain = line[21:22].strip()
    resseq_str = line[22:26].strip()
    try:
        resseq = int(resseq_str)
    except ValueError:
        resseq = -10**9  # unlikely sentinel when missing
    icode = line[26:27].strip()
    key_full = (chain, resseq, icode, resname, atomname, altloc)
    key_simple = (chain, resseq, icode, atomname)
    return key_full, key_simple


def collect_ml_atom_keys(model_pdb: Path) -> Tuple[set, set]:
    """Collect ML-region atom keys from model_pdb.

    Returns:
        keys_full: Set of (chain, resseq, icode, resname, atomname, altloc)
        keys_simple: Set of (chain, resseq, icode, atomname)
    """
    from typing import Set as SetType
    keys_full: SetType[Tuple] = set()
    keys_simple: SetType[Tuple] = set()
    try:
        with model_pdb.open("r") as fh:
            for line in fh:
                if line.startswith("ATOM") or line.startswith("HETATM"):
                    kf, ks = pdb_keys_from_line(line)
                    keys_full.add(kf)
                    keys_simple.add(ks)
    except Exception:
        # If anything goes wrong, leave sets empty; caller will handle gracefully.
        pass
    return keys_full, keys_simple


def format_pdb_with_bfactor(line: str, b: float) -> str:
    """Return PDB line with B-factor field (cols 61-66) set to b (6.2f)."""
    if len(line) < 66:
        line = line.rstrip("\n")
        line = line + " " * max(0, 66 - len(line))
        line = line + "\n"
    bf_str = f"{b:6.2f}"
    # Preserve occupancy (cols 55-60), overwrite tempFactor (61-66).
    new_line = line[:60] + bf_str + line[66:]
    return new_line


def annotate_pdb_bfactors_inplace(
    pdb_path: Path,
    model_pdb: Path,
    freeze_indices_0based: Sequence[int],
    beta_ml: float = 100.0,
    beta_frz: float = 50.0,
    beta_both: float = 150.0,
) -> None:
    """Overwrite B-factors in-place for ML/frozen atoms visualization.

    - ML-region atoms: beta_ml (default 100.00)
    - frozen atoms: beta_frz (default 50.00)
    - ML ∩ frozen: beta_both (default 150.00)

    Indexing for 'frozen' is 0-based and resets at each MODEL.
    """
    ml_full, ml_simple = collect_ml_atom_keys(model_pdb)
    frozen_set = set(int(i) for i in (freeze_indices_0based or []))

    try:
        lines = pdb_path.read_text().splitlines(keepends=True)
    except Exception:
        return

    out_lines: List[str] = []
    atom_idx = 0  # resets per MODEL

    for line in lines:
        rec = line[:6]
        if rec.startswith("MODEL"):
            # reset atom counter for each model
            atom_idx = 0
            out_lines.append(line)
            continue
        if rec.startswith("ATOM  ") or rec.startswith("HETATM"):
            kf, ks = pdb_keys_from_line(line)
            is_ml = (kf in ml_full) or (ks in ml_simple)
            is_frz = (atom_idx in frozen_set)
            if is_ml and is_frz:
                out_lines.append(format_pdb_with_bfactor(line, beta_both))
            elif is_ml:
                out_lines.append(format_pdb_with_bfactor(line, beta_ml))
            elif is_frz:
                out_lines.append(format_pdb_with_bfactor(line, beta_frz))
            else:
                out_lines.append(line)
            atom_idx += 1
        else:
            out_lines.append(line)

    try:
        pdb_path.write_text("".join(out_lines))
    except Exception:
        # Silently ignore if we cannot write; conversion outputs are still present.
        pass


def convert_and_annotate_xyz_to_pdb(
    src_xyz_or_trj: Path,
    ref_pdb: Path,
    dst_pdb: Path,
    model_pdb: Path,
    freeze_indices_0based: Sequence[int],
) -> None:
    """Convert an XYZ/TRJ file to PDB and annotate B-factors to highlight ML and frozen atoms.

    This mirrors the behaviour of the `opt` tool:
      - ML-region atoms: 100.00
      - frozen atoms: 50.00
      - ML ∩ frozen: 150.00
    """
    try:
        convert_xyz_to_pdb(src_xyz_or_trj, ref_pdb, dst_pdb)
        annotate_pdb_bfactors_inplace(
            dst_pdb,
            model_pdb=model_pdb,
            freeze_indices_0based=freeze_indices_0based,
        )
    except Exception as exc:
        click.echo(
            f"[convert] WARNING: Failed to convert '{src_xyz_or_trj}' to PDB: {exc}",
            err=True,
        )


# =============================================================================
# Input preparation helpers
# =============================================================================


@dataclass
class PreparedInputStructure:
    source_path: Path
    geom_path: Path

    def cleanup(self) -> None:
        """No-op: no temporary files are created."""
        return None

    def __enter__(self) -> "PreparedInputStructure":
        return self

    def __exit__(self, exc_type, exc, tb) -> None:
        self.cleanup()


def prepare_input_structure(path: Path) -> PreparedInputStructure:
    """Return a lightweight wrapper for the provided structure path."""
    return PreparedInputStructure(source_path=path, geom_path=path)


def _count_atoms_in_file(path: Path) -> int:
    """Count atoms in a structure file (PDB or XYZ)."""
    suffix = path.suffix.lower()
    if suffix == ".pdb":
        count = 0
        with open(path, "r") as f:
            for line in f:
                if line.startswith(("ATOM  ", "HETATM")):
                    count += 1
        return count
    elif suffix == ".xyz":
        # XYZ format: first line is atom count
        with open(path, "r") as f:
            first_line = f.readline().strip()
            try:
                return int(first_line)
            except ValueError:
                return 0
    return 0


def apply_ref_pdb_override(
    prepared_input: PreparedInputStructure,
    ref_pdb: Optional[Path],
) -> Optional[Path]:
    """Use a reference PDB topology while keeping XYZ coordinates for geometry loading.

    When --ref-pdb is provided:
    - geom_path remains the original input (xyz) for high-precision coordinates
    - source_path is updated to ref_pdb for topology/residue information
    """
    import click
    if ref_pdb is None:
        return None
    ref_pdb = Path(ref_pdb).resolve()
    if ref_pdb.suffix.lower() != ".pdb":
        raise click.BadParameter("--ref-pdb must be a .pdb file.")
    geom_count = _count_atoms_in_file(prepared_input.geom_path)
    ref_count = _count_atoms_in_file(ref_pdb)
    if geom_count != ref_count:
        raise click.BadParameter(
            f"Atom count mismatch: {prepared_input.geom_path.name} has {geom_count} atoms, "
            f"but --ref-pdb {ref_pdb.name} has {ref_count} atoms."
        )
    prepared_input.source_path = ref_pdb
    return ref_pdb


def resolve_charge_spin_or_raise(
    prepared: PreparedInputStructure,
    charge: Optional[int],
    spin: Optional[int],
    *,
    spin_default: int = 1,
    charge_default: Optional[int] = None,
) -> Tuple[int, int]:
    """Resolve charge/spin from inputs.

    The ``prepared`` argument is accepted for API stability and is not used.
    Charge defaults to unresolved; callers must provide an explicit value or
    pass ``charge_default`` intentionally.
    """
    if charge is None:
        if charge_default is None:
            raise click.ClickException(
                "Total charge is unresolved. Provide -q/--charge."
            )
        charge = charge_default
    if spin is None:
        spin = spin_default
    return int(charge), int(spin)


# -----------------------------------------------
# B-factor based 3-layer ML/MM system utilities
# -----------------------------------------------

def read_bfactors_from_pdb(pdb_path: Path) -> List[float]:
    """
    Read B-factor (temperature factor) values from a PDB file.

    Returns a list of B-factors in atom order (0-indexed).
    Only ATOM and HETATM records are processed.
    """
    bfactors: List[float] = []
    with open(pdb_path, "r") as f:
        for line in f:
            if line.startswith(("ATOM  ", "HETATM")):
                # B-factor is at columns 61-66 (1-indexed), i.e., [60:66]
                try:
                    bfac = float(line[60:66].strip())
                except (ValueError, IndexError):
                    bfac = 0.0
                bfactors.append(bfac)
    return bfactors


def parse_layer_indices_from_bfactors(
    bfactors: List[float],
    tolerance: float = 1.0,
) -> Dict[str, List[int]]:
    """
    Parse B-factor values into layer indices for 3-layer ML/MM system.

    B-factor encoding:
        0.0 (±tolerance): ML atoms
        10.0 (±tolerance): Movable MM atoms
        20.0 (±tolerance): Frozen MM atoms

    Parameters
    ----------
    bfactors : List[float]
        B-factor values for each atom (0-indexed).
    tolerance : float
        Tolerance for B-factor matching (default: 1.0).

    Returns
    -------
    Dict[str, List[int]]
        Dictionary with keys:
        - "ml_indices": ML region atoms
        - "hess_mm_indices": Compatibility key (empty in 3-layer encoding)
        - "movable_mm_indices": Movable MM atoms
        - "frozen_indices": Frozen atoms
        - "unassigned_indices": Atoms with B-factors not matching any layer
    """
    from .defaults import BFACTOR_ML, BFACTOR_HESS_MM, BFACTOR_MOVABLE_MM, BFACTOR_FROZEN

    ml_indices: List[int] = []
    hess_mm_indices: List[int] = []
    movable_mm_indices: List[int] = []
    frozen_indices: List[int] = []
    unassigned_indices: List[int] = []

    for i, bfac in enumerate(bfactors):
        if abs(bfac - BFACTOR_ML) <= tolerance:
            ml_indices.append(i)
        elif abs(bfac - BFACTOR_FROZEN) <= tolerance:
            frozen_indices.append(i)
        elif abs(bfac - BFACTOR_MOVABLE_MM) <= tolerance:
            movable_mm_indices.append(i)
        elif (
            BFACTOR_HESS_MM != BFACTOR_MOVABLE_MM
            and abs(bfac - BFACTOR_HESS_MM) <= tolerance
        ):
            hess_mm_indices.append(i)
        else:
            unassigned_indices.append(i)

    return {
        "ml_indices": ml_indices,
        "hess_mm_indices": hess_mm_indices,
        "movable_mm_indices": movable_mm_indices,
        "frozen_indices": frozen_indices,
        "unassigned_indices": unassigned_indices,
    }


def has_valid_layer_bfactors(bfactors: List[float], tolerance: float = 1.0) -> bool:
    """
    Check if PDB B-factors contain valid 3-layer encoding.

    Returns True if at least one atom has ML B-factor and the B-factors are
    predominantly in the expected range (0, 10, 20).
    """
    from .defaults import BFACTOR_ML, BFACTOR_HESS_MM, BFACTOR_MOVABLE_MM, BFACTOR_FROZEN

    valid_bfactors = {BFACTOR_ML, BFACTOR_MOVABLE_MM, BFACTOR_FROZEN, BFACTOR_HESS_MM}
    has_ml = False
    valid_count = 0

    for bfac in bfactors:
        for valid in valid_bfactors:
            if abs(bfac - valid) <= tolerance:
                valid_count += 1
                if abs(bfac - BFACTOR_ML) <= tolerance:
                    has_ml = True
                break

    # Consider valid if:
    # 1. Has at least one ML atom
    # 2. At least 80% of atoms have valid B-factors
    return has_ml and (valid_count / max(len(bfactors), 1) >= 0.8)


def parse_indices_string(indices_str: str, one_based: bool = True) -> List[int]:
    """
    Parse a comma-separated index string into a sorted list of 0-based ints.

    Supports ranges like "1-5" (inclusive). By default, inputs are 1-based.
    """
    import click
    if indices_str is None:
        return []
    tokens = [tok.strip() for tok in str(indices_str).replace(" ", ",").split(",") if tok.strip()]
    indices: List[int] = []
    for token in tokens:
        if "-" in token and not token.startswith("-"):
            parts = token.split("-")
            if len(parts) == 2 and parts[0] and parts[1]:
                try:
                    start = int(parts[0])
                    end = int(parts[1])
                except ValueError as exc:
                    raise click.BadParameter(f"Invalid range token in --model-indices: '{token}'") from exc
                if one_based:
                    start -= 1
                    end -= 1
                if start < 0 or end < 0 or start > end:
                    raise click.BadParameter(f"Invalid range in --model-indices: '{token}'")
                indices.extend(range(start, end + 1))
                continue
        try:
            value = int(token)
        except ValueError as exc:
            raise click.BadParameter(f"Invalid index in --model-indices: '{token}'") from exc
        if one_based:
            value -= 1
        if value < 0:
            raise click.BadParameter(f"--model-indices expects positive indices; got {value + (1 if one_based else 0)}")
        indices.append(value)
    return sorted(set(indices))


def write_model_pdb_from_indices(
    input_pdb_path: Path,
    output_pdb_path: Path,
    indices: Sequence[int],
) -> None:
    """
    Write a model PDB containing only atoms at the specified 0-based indices.
    """
    import click
    if not indices:
        raise ValueError("No indices provided to build model PDB.")
    n_atoms = _count_atoms_in_file(input_pdb_path)
    if n_atoms <= 0:
        raise ValueError(f"No atoms found in input PDB: {input_pdb_path}")
    for idx in indices:
        if idx < 0 or idx >= n_atoms:
            raise click.BadParameter(
                f"model index out of range: {idx} (valid: 0 <= idx < {n_atoms})"
            )

    keep = set(int(i) for i in indices)
    lines_out: List[str] = []
    atom_idx = 0
    with open(input_pdb_path, "r") as f:
        for line in f:
            if line.startswith(("ATOM  ", "HETATM")):
                if atom_idx in keep:
                    # Auto-fill element column (77-78) if missing
                    raw = line.rstrip("\n")
                    elem_field = raw[76:78].strip() if len(raw) >= 78 else ""
                    if not elem_field:
                        atom_name = raw[12:16].strip()
                        res_name = raw[17:20].strip()
                        is_hetatm = raw.startswith("HETATM")
                        elem = guess_element(atom_name, res_name, is_hetatm)
                        if elem:
                            padded = raw.ljust(76) + f"{elem:>2}" + "\n"
                            lines_out.append(padded)
                        else:
                            lines_out.append(line)
                    else:
                        lines_out.append(line)
                atom_idx += 1
    if not lines_out:
        raise ValueError("Model PDB would be empty; check indices and input PDB.")
    if not lines_out[-1].endswith("\n"):
        lines_out[-1] = lines_out[-1] + "\n"
    lines_out.append("END\n")
    with open(output_pdb_path, "w") as f:
        f.writelines(lines_out)


def build_model_pdb_from_indices(
    input_pdb_path: Path,
    out_dir: Path,
    indices: Sequence[int],
    *,
    label: str = "model_from_indices",
) -> Path:
    """
    Create a temporary model PDB under out_dir using explicit indices.
    """
    out_dir.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile(
        mode="w",
        suffix=".pdb",
        prefix=f"{label}_",
        dir=out_dir,
        delete=False,
    ) as tmp:
        tmp_path = Path(tmp.name)
    write_model_pdb_from_indices(input_pdb_path, tmp_path, indices)
    return tmp_path


def build_model_pdb_from_bfactors(
    input_pdb_path: Path,
    out_dir: Path,
    *,
    tolerance: float = None,
    label: str = "model_from_bfactor",
) -> Tuple[Path, Dict[str, List[int]]]:
    """
    Create a model PDB using ML indices derived from B-factors.

    Returns (model_pdb_path, layer_info).
    """
    from .defaults import BFACTOR_TOLERANCE
    tol = BFACTOR_TOLERANCE if tolerance is None else float(tolerance)
    bfactors = read_bfactors_from_pdb(input_pdb_path)
    if not bfactors:
        raise ValueError(f"No ATOM/HETATM records found in {input_pdb_path}.")
    if not has_valid_layer_bfactors(bfactors, tolerance=tol):
        raise ValueError(
            "Invalid or missing layer B-factors (expected ~0/10/20). "
            "Provide --no-detect-layer with --model-pdb/--model-indices."
        )
    layer_info = parse_layer_indices_from_bfactors(bfactors, tolerance=tol)
    ml_indices = layer_info.get("ml_indices") or []
    if not ml_indices:
        raise ValueError("No ML atoms detected from B-factors (value ~0).")
    out_dir.mkdir(parents=True, exist_ok=True)
    tmp_path = out_dir / f"{label}.pdb"
    write_model_pdb_from_indices(input_pdb_path, tmp_path, ml_indices)
    return tmp_path, layer_info


def write_layer_bfactors_to_pdb(
    input_pdb_path: Path,
    output_pdb_path: Path,
    ml_indices: List[int],
    hess_mm_indices: Optional[List[int]] = None,
    movable_mm_indices: Optional[List[int]] = None,
    frozen_indices: Optional[List[int]] = None,
) -> None:
    """
    Write a PDB file with B-factors set according to 3-layer assignments.

    B-factor encoding:
        ML atoms: 0.0
        Movable MM atoms: 10.0
        Frozen MM atoms: 20.0
        Hessian MM atoms: encoded with the same B-factor as movable MM

    Parameters
    ----------
    input_pdb_path : Path
        Source PDB file to read atom records from.
    output_pdb_path : Path
        Output PDB file path.
    ml_indices : List[int]
        0-based indices of ML region atoms.
    hess_mm_indices : Optional[List[int]]
        0-based indices of MM atoms with Hessian (written as movable B-factor).
    movable_mm_indices : Optional[List[int]]
        0-based indices of movable MM atoms without Hessian.
    frozen_indices : Optional[List[int]]
        0-based indices of frozen atoms.

    Notes
    -----
    Supports multi-MODEL PDB files (e.g., trajectories): atom index resets
    at each MODEL record.
    """
    from .defaults import BFACTOR_ML, BFACTOR_HESS_MM, BFACTOR_MOVABLE_MM, BFACTOR_FROZEN

    ml_set = set(ml_indices or [])
    hess_mm_set = set(hess_mm_indices or [])
    movable_mm_set = set(movable_mm_indices or [])
    frozen_set = set(frozen_indices or [])

    lines_out: List[str] = []
    atom_idx = 0

    with open(input_pdb_path, "r") as f:
        for line in f:
            rec = line[:6]
            # Reset atom counter at each MODEL record (for trajectory files)
            if rec.startswith("MODEL"):
                atom_idx = 0
                lines_out.append(line)
                continue

            if line.startswith(("ATOM  ", "HETATM")):
                # Determine B-factor for this atom
                if atom_idx in ml_set:
                    bfac = BFACTOR_ML
                elif atom_idx in hess_mm_set:
                    bfac = BFACTOR_HESS_MM
                elif atom_idx in movable_mm_set:
                    bfac = BFACTOR_MOVABLE_MM
                elif atom_idx in frozen_set:
                    bfac = BFACTOR_FROZEN
                else:
                    # Default: treat as movable MM (layer 3)
                    bfac = BFACTOR_MOVABLE_MM

                # Replace B-factor (columns 61-66, 1-indexed)
                # PDB format: columns 61-66 are B-factor with format %6.2f
                # Ensure line is long enough before modifying
                if len(line) >= 66:
                    new_line = line[:60] + f"{bfac:6.2f}" + line[66:]
                else:
                    # Pad line if too short
                    padded = line.rstrip("\n").ljust(66)
                    new_line = padded[:60] + f"{bfac:6.2f}" + "\n"
                lines_out.append(new_line)
                atom_idx += 1
            else:
                lines_out.append(line)

    with open(output_pdb_path, "w") as f:
        f.writelines(lines_out)


def update_pdb_bfactors_from_layers(
    pdb_path: Path,
    ml_indices: List[int],
    hess_mm_indices: Optional[List[int]] = None,
    movable_mm_indices: Optional[List[int]] = None,
    frozen_indices: Optional[List[int]] = None,
) -> None:
    """
    Update B-factors in a PDB file in-place based on layer assignments.

    This is a convenience wrapper that reads and writes to the same file.
    """
    import tempfile
    import shutil

    with tempfile.NamedTemporaryFile(mode="w", suffix=".pdb", delete=False) as tmp:
        tmp_path = Path(tmp.name)

    try:
        write_layer_bfactors_to_pdb(
            pdb_path,
            tmp_path,
            ml_indices,
            hess_mm_indices,
            movable_mm_indices,
            frozen_indices,
        )
        shutil.move(str(tmp_path), str(pdb_path))
    finally:
        if tmp_path.exists():
            tmp_path.unlink()
