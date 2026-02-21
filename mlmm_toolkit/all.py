# mlmm_toolkit/all.py

"""
End-to-end enzymatic reaction workflow: extract, scan, MEP, TS, IRC, freq, DFT.

Example:
    mlmm all -i R.pdb P.pdb -c 'GPP,MMT' --ligand-charge 'GPP:-3,MMT:-1'

For detailed documentation, see: docs/all.md
"""

from __future__ import annotations

import ast
from collections import defaultdict
from pathlib import Path
from typing import List, Sequence, Optional, Tuple, Dict, Any
import shutil
import tempfile
import os

import sys
import math
import click
from click.core import ParameterSource
import time  # timing
import yaml
import numpy as np
import torch

# Biopython for PDB parsing (post-processing helpers)
from Bio import PDB

# pysisyphus helpers/constants
from pysisyphus.helpers import geom_loader
from pysisyphus.constants import BOHR2ANG, AU2KCALPERMOL

# Local imports from the package
from .extract import extract_api
from . import path_search as _path_search
from . import tsopt as _ts_opt
from . import freq as _freq_cli
from . import dft as _dft_cli

from .trj2fig import run_trj2fig
from .summary_log import write_summary_log
from .utils import (
    build_energy_diagram,
    deep_update,
    ensure_dir,
    format_elapsed,
    prepare_input_structure,
    collect_single_option_values,
    load_yaml_dict,
    load_pdb_atom_metadata,
    parse_scan_list_triples,
    read_xyz_as_blocks,
)
from .preflight import validate_existing_files, ensure_commands_available
from . import scan as _scan_cli
from .add_elem_info import assign_elements as _assign_elem_info
from .define_layer import define_layers as _define_layers
from .mlmm_calc import mlmm as _mlmm_calc
from .mm_parm import (
    Args as _AutoMMArgs,
    parse_ligand_charge as _mm_parse_ligand_charge,
    parse_ligand_mult as _mm_parse_ligand_mult,
    run_pipeline as _mm_run,
)

AtomKey = Tuple[str, str, str, str, str, str]

_log_started = False


def _echo(*args, **kwargs) -> None:
    """Echo with local output tracking for section spacing."""
    global _log_started
    click.echo(*args, **kwargs)
    _log_started = True


def _echo_section(message: str, **kwargs) -> None:
    """Echo a section header with a leading blank line unless it's the first log."""
    global _log_started
    if _log_started:
        click.echo()
    click.echo(message, **kwargs)
    _log_started = True


def _first_existing_artifact(out_dir: Path, patterns: Sequence[str]) -> Optional[Path]:
    """Resolve the first existing file for a list of relative patterns."""
    for pattern in patterns:
        if any(ch in pattern for ch in "*?[]"):
            for candidate in sorted(out_dir.glob(pattern)):
                if candidate.is_file():
                    return candidate.resolve()
            continue
        candidate = out_dir / pattern
        if candidate.is_file():
            return candidate.resolve()
    return None


def _link_or_copy_file(src: Path, dst: Path) -> bool:
    """Create a symlink when possible; fall back to copy."""
    try:
        if dst.exists() or dst.is_symlink():
            if dst.is_dir():
                return False
            dst.unlink()
        rel = os.path.relpath(src, start=dst.parent)
        dst.symlink_to(rel)
        return True
    except Exception:
        try:
            shutil.copy2(src, dst)
            return True
        except Exception:
            return False


def _write_output_summary_md(out_dir: Path) -> None:
    """Write summary.md and expose key outputs at out_dir root."""
    try:
        out_dir = out_dir.resolve()
        if not out_dir.exists():
            return

        root_specs: List[Tuple[str, str]] = [
            ("summary.yaml", "YAML summary"),
            ("summary.log", "Text summary"),
            ("mep.trj", "MEP trajectory"),
            ("mep.pdb", "MEP trajectory (PDB)"),
            ("mep_w_ref.pdb", "MEP merged with reference"),
            ("mep_plot.png", "MEP profile plot"),
            ("energy_diagram_MEP.png", "State energy diagram"),
            ("energy_diagram_UMA_all.png", "All-segment UMA diagram"),
            ("irc_plot_all.png", "Aggregated IRC plot"),
        ]
        root_lines: List[str] = []
        for rel, label in root_specs:
            if (out_dir / rel).is_file():
                root_lines.append(f"- {label}: [`{rel}`]({rel})")

        shortcut_specs: List[Tuple[str, str, Sequence[str]]] = [
            ("key_mep.trj", "Primary MEP trajectory", ["mep.trj", "path_search/mep.trj", "path_opt/final_geometries.trj"]),
            ("key_mep.pdb", "Primary MEP PDB", ["mep.pdb", "path_search/mep.pdb", "path_opt/final_geometries.pdb"]),
            ("key_ts.pdb", "TS structure (PDB)", ["ts_seg_01.pdb", "path_search/post_seg_*/ts/final_geometry.pdb", "path_opt/post_seg_*/ts/final_geometry.pdb", "tsopt_single/ts/final_geometry.pdb"]),
            ("key_ts.xyz", "TS structure (XYZ)", ["ts_seg_01.xyz", "path_search/post_seg_*/ts/final_geometry.xyz", "path_opt/post_seg_*/ts/final_geometry.xyz", "tsopt_single/ts/final_geometry.xyz"]),
            ("key_freq_TS.csv", "TS frequencies", ["path_search/post_seg_*/freq/TS/frequencies.csv", "path_opt/post_seg_*/freq/TS/frequencies.csv", "tsopt_single/freq/TS/frequencies.csv"]),
            ("key_dft_TS.yaml", "TS DFT result", ["path_search/post_seg_*/dft/TS/result.yaml", "path_opt/post_seg_*/dft/TS/result.yaml", "tsopt_single/dft/TS/result.yaml"]),
            ("key_irc_plot.png", "IRC plot", ["irc_plot_all.png", "path_search/post_seg_*/irc/irc_plot.png", "path_opt/post_seg_*/irc/irc_plot.png", "tsopt_single/irc/irc_plot.png"]),
        ]
        shortcut_lines: List[str] = []
        for name, label, patterns in shortcut_specs:
            src = _first_existing_artifact(out_dir, patterns)
            if src is None:
                continue
            dst = out_dir / name
            try:
                same = src.resolve() == dst.resolve()
            except Exception:
                same = False
            if not same and not _link_or_copy_file(src, dst):
                continue
            src_rel = os.path.relpath(src, start=out_dir)
            shortcut_lines.append(
                f"- {label}: [`{name}`]({name}) (source: `{src_rel}`)"
            )

        lines: List[str] = [
            "# Run Summary",
            "",
            f"- Generated: `{time.strftime('%Y-%m-%d %H:%M:%S %Z')}`",
            f"- Output directory: `{out_dir}`",
            "",
            "## Primary Artifacts",
        ]
        if root_lines:
            lines.extend(root_lines)
        else:
            lines.append("- No primary artifacts detected yet.")
        lines.extend(["", "## Root Shortcuts", "- `key_*` files are symlinks when possible, otherwise copied files."])
        if shortcut_lines:
            lines.extend(shortcut_lines)
        else:
            lines.append("- No shortcuts were generated.")
        lines.extend(["", "## Notes", "- Start from `summary.log` for a concise narrative and from `summary.yaml` for machine-readable values."])

        summary_md = out_dir / "summary.md"
        summary_md.write_text("\n".join(lines).rstrip() + "\n", encoding="utf-8")
        _echo(f"[write] Wrote '{summary_md}'.")
    except Exception as e:
        _echo(f"[write] WARNING: Failed to write summary.md: {e}", err=True)


def _run_cli_main(
    cmd_name: str,
    cli_obj,
    args: Sequence[str],
    *,
    on_nonzero: str = "warn",
    on_exception: str = "raise",
    prefix: Optional[str] = None,
) -> None:
    """Run a Click command with temporary argv and consistent error handling."""
    saved = list(sys.argv)
    label = prefix or cmd_name
    try:
        sys.argv = ["mlmm", cmd_name] + list(args)
        _echo("\n")
        cli_obj.main(args=list(args), standalone_mode=False)
    except SystemExit as e:
        code = getattr(e, "code", 1)
        if code not in (None, 0):
            if on_nonzero == "raise":
                raise click.ClickException(f"[{label}] {cmd_name} exit code {code}.")
            _echo(f"[{label}] WARNING: {cmd_name} exited with code {code}")
    except Exception as e:
        if on_exception == "raise":
            raise click.ClickException(f"[{label}] {cmd_name} failed: {e}")
        _echo(f"[{label}] WARNING: {cmd_name} failed: {e}")
    finally:
        sys.argv = saved
        _echo("\n")


# -----------------------------
# Helpers
# -----------------------------

def _append_cli_arg(args: List[str], flag: str, value: Any | None) -> None:
    """Append ``flag`` and ``value`` (converted to string) to ``args`` when ``value`` is not ``None``."""
    if value is None:
        return
    if isinstance(value, bool):
        args.extend([flag, "True" if value else "False"])
    else:
        args.extend([flag, str(value)])


def _resolve_override_dir(default: Path, override: Path | None) -> Path:
    """Return ``override`` when provided (respecting absolute paths); otherwise ``default``."""
    if override is None:
        return default
    if override.is_absolute():
        return override
    return default.parent / override


def _resolve_yaml_sources(
    config_yaml: Optional[Path],
    override_yaml: Optional[Path],
    args_yaml_legacy: Optional[Path],
) -> Tuple[Optional[Path], Optional[Path], bool]:
    """Resolve config/override YAML inputs and legacy alias usage."""
    if override_yaml is not None and args_yaml_legacy is not None:
        raise click.BadParameter(
            "Use either --override-yaml or --args-yaml (legacy alias), not both."
        )

    used_legacy_alias = False
    if args_yaml_legacy is not None:
        used_legacy_alias = True
        override_yaml = args_yaml_legacy

    return config_yaml, override_yaml, used_legacy_alias


def _build_effective_args_yaml(
    config_yaml: Optional[Path],
    override_yaml: Optional[Path],
    *,
    tmp_prefix: str,
) -> Tuple[Optional[Path], Dict[str, Any]]:
    """
    Build an effective args-yaml file path.

    Precedence for file layering:
      config_yaml < override_yaml
    """
    base_cfg = load_yaml_dict(config_yaml)
    override_cfg = load_yaml_dict(override_yaml)

    if config_yaml is None and override_yaml is None:
        return None, {}
    if config_yaml is None:
        return override_yaml, override_cfg
    if override_yaml is None:
        return config_yaml, base_cfg

    merged: Dict[str, Any] = {}
    deep_update(merged, base_cfg)
    deep_update(merged, override_cfg)

    with tempfile.NamedTemporaryFile(
        mode="w",
        encoding="utf-8",
        suffix=".yaml",
        prefix=tmp_prefix,
        delete=False,
    ) as tf:
        yaml.safe_dump(merged, tf, sort_keys=False, allow_unicode=True)
        effective = Path(tf.name).resolve()

    return effective, merged


def _write_ml_region_definition(pocket_pdb: Path, dest: Path) -> Path:
    """
    Copy ``pocket_pdb`` to ``dest`` for downstream ML/MM commands.

    The copy preserves whatever link-hydrogen policy was used during extraction; set ``--add-linkH False``
    if you need a link-free ML-region definition.
    """
    dest.parent.mkdir(parents=True, exist_ok=True)
    try:
        shutil.copyfile(pocket_pdb, dest)
    except FileNotFoundError:
        raise click.ClickException(f"[all] Pocket PDB not found while building ML region: {pocket_pdb}")
    return dest.resolve()


def _mm_charge_mapping(expr: Optional[str]) -> Dict[str, int]:
    """Return a ligand-charge mapping for mm_parm when ``expr`` uses RES=Q or RES:Q syntax."""
    if not expr:
        return {}
    if ("=" not in expr) and (":" not in expr):
        return {}
    try:
        return _mm_parse_ligand_charge(expr)
    except Exception as exc:  # pragma: no cover - defensive
        raise click.ClickException(f"[all] Invalid --ligand-charge mapping for mm_parm: {exc}")


def _mm_mult_mapping(expr: Optional[str]) -> Dict[str, int]:
    """Return a ligand-multiplicity mapping for mm_parm when ``expr`` uses RES=M or RES:M syntax."""
    if not expr:
        return {}
    try:
        return _mm_parse_ligand_mult(expr)
    except Exception as exc:  # pragma: no cover - defensive
        raise click.ClickException(f"[all] Invalid --auto-mm-ligand-mult mapping for mm_parm: {exc}")


def _build_mm_parm7(
    pdb: Path,
    ligand_charge_expr: Optional[str],
    ligand_mult_expr: Optional[str],
    out_dir: Path,
    ff_set: str,
    add_TER: bool,
    keep_temp: bool,
    allow_nonstandard_aa: bool,
) -> Tuple[Path, Path]:
    """Run mm_parm on ``pdb`` and return (parm7, rst7)."""
    out_dir.mkdir(parents=True, exist_ok=True)
    out_prefix = (out_dir / pdb.stem).resolve()
    out_prefix.parent.mkdir(parents=True, exist_ok=True)
    args = _AutoMMArgs(
        pdb=pdb.resolve(),
        out_prefix=str(out_prefix),
        ligand_charge=_mm_charge_mapping(ligand_charge_expr),
        ligand_mult=_mm_mult_mapping(ligand_mult_expr),
        allow_nonstandard_aa=bool(allow_nonstandard_aa),
        keep_temp=bool(keep_temp),
        add_TER=bool(add_TER),
        add_H=False,
        ph=7.0,
        ff_set=str(ff_set),
        out_prefix_given=True,
    )
    try:
        _mm_run(args)
    except SystemExit as exc:  # pragma: no cover - click exit translation
        code = getattr(exc, "code", 1)
        raise click.ClickException(f"[all] mm_parm exited with code {code}.")
    except Exception as exc:
        raise click.ClickException(f"[all] mm_parm failed: {exc}")

    parm7 = Path(f"{args.out_prefix}.parm7").resolve()
    rst7 = Path(f"{args.out_prefix}.rst7").resolve()
    if not parm7.exists():
        raise click.ClickException(f"[all] mm_parm did not produce parm7 at {parm7}")
    if not rst7.exists():
        raise click.ClickException(f"[all] mm_parm did not produce rst7 at {rst7}")
    return parm7, rst7


def _parse_atom_key_from_line(line: str) -> Optional[AtomKey]:
    """Extract a structural identity key from a PDB ATOM/HETATM record."""
    if not (line.startswith("ATOM") or line.startswith("HETATM")):
        return None
    atomname = line[12:16].strip()
    altloc = (line[16] if len(line) > 16 else " ").strip()
    resname = line[17:20].strip()
    chain = (line[21] if len(line) > 21 else " ").strip()
    resseq = line[22:26].strip()
    icode = (line[26] if len(line) > 26 else " ").strip()
    return (chain, resname, resseq, icode, atomname, altloc)


def _key_variants(key: AtomKey) -> List[AtomKey]:
    """Return key variants with progressively relaxed identity fields (deduplicated)."""
    chain, resn, resseq, icode, atom, alt = key
    raw_variants = [
        (chain, resn, resseq, icode, atom, alt),
        (chain, resn, resseq, icode, atom, ""),
        (chain, resn, resseq, "", atom, alt),
        (chain, resn, resseq, "", atom, ""),
    ]
    seen: set[AtomKey] = set()
    variants: List[AtomKey] = []
    for variant in raw_variants:
        if variant in seen:
            continue
        seen.add(variant)
        variants.append(variant)
    return variants


def _build_variant_occurrence_table(keys: Sequence[AtomKey]) -> List[Dict[AtomKey, int]]:
    """Track how many times each relaxed key variant has appeared up to each atom index."""
    counts: Dict[AtomKey, int] = defaultdict(int)
    per_atom: List[Dict[AtomKey, int]] = []
    for key in keys:
        current: Dict[AtomKey, int] = {}
        for variant in _key_variants(key):
            counts[variant] += 1
            current[variant] = counts[variant]
        per_atom.append(current)
    return per_atom


def _pocket_key_to_index(pocket_pdb: Path) -> Dict[AtomKey, List[int]]:
    """Build mapping: structural atom key -> list of pocket indices (1-based by file order)."""
    key2idx: Dict[AtomKey, List[int]] = defaultdict(list)
    idx = 0
    try:
        with open(pocket_pdb, "r", encoding="utf-8", errors="ignore") as fh:
            for line in fh:
                if line.startswith("ATOM") or line.startswith("HETATM"):
                    key = _parse_atom_key_from_line(line)
                    if key is None:
                        continue
                    idx += 1
                    for variant in _key_variants(key):
                        key2idx[variant].append(idx)
    except FileNotFoundError:
        raise click.ClickException(f"[all] Pocket PDB not found: {pocket_pdb}")
    if not key2idx:
        raise click.ClickException(f"[all] Pocket PDB {pocket_pdb} has no ATOM/HETATM records.")
    return dict(key2idx)


def _read_full_atom_keys_in_file_order(full_pdb: Path) -> List[AtomKey]:
    """Read ATOM/HETATM lines and return keys in the original file order."""
    keys: List[AtomKey] = []
    try:
        with open(full_pdb, "r", encoding="utf-8", errors="ignore") as fh:
            for line in fh:
                if line.startswith("ATOM") or line.startswith("HETATM"):
                    key = _parse_atom_key_from_line(line)
                    if key is not None:
                        keys.append(key)
    except FileNotFoundError:
        raise click.ClickException(f"[all] File not found while parsing PDB: {full_pdb}")
    if not keys:
        raise click.ClickException(f"[all] No ATOM/HETATM records detected in {full_pdb}.")
    return keys


def _format_atom_key_for_msg(key: AtomKey) -> str:
    """Pretty string for diagnostics."""
    chain, resn, resseq, icode, atom, alt = key
    res = f"{chain}:{resn}{resseq}{(icode if icode else '')}"
    alt_sfx = f",alt={alt}" if alt else ""
    return f"{res}:{atom}{alt_sfx}"


def _parse_scan_lists_literals(
    scan_lists_raw: Sequence[str],
    atom_meta: Optional[Sequence[Dict[str, Any]]] = None,
) -> List[List[Tuple[int, int, float]]]:
    """Parse ``--scan-lists`` literals without re-basing atom indices."""
    stages: List[List[Tuple[int, int, float]]] = []
    for idx_stage, literal in enumerate(scan_lists_raw, start=1):
        tuples, _ = parse_scan_list_triples(
            literal,
            one_based=True,
            atom_meta=atom_meta,
            option_name=f"--scan-lists #{idx_stage}",
            return_one_based=True,
        )
        if not tuples:
            raise click.BadParameter(
                f"--scan-lists #{idx_stage} must contain at least one (i,j,target) triple."
            )
        stages.append(tuples)
    return stages


def _format_scan_stage(stage: List[Tuple[int, int, float]]) -> str:
    """Serialize a scan stage back into a Python-like literal string."""
    return "[" + ", ".join(f"({i},{j},{target})" for (i, j, target) in stage) + "]"


def _convert_scan_lists_to_pocket_indices(
    scan_lists_raw: Sequence[str],
    full_input_pdb: Path,
    pocket_pdb: Path,
) -> List[List[Tuple[int, int, float]]]:
    """
    Convert user-provided atom indices (based on the full input PDB) to pocket indices.
    Returns the converted stages as lists of (i,j,target) with 1-based pocket indices.
    """
    if not scan_lists_raw:
        return []

    full_atom_meta = load_pdb_atom_metadata(full_input_pdb)
    stages = _parse_scan_lists_literals(scan_lists_raw, atom_meta=full_atom_meta)

    orig_keys_in_order = _read_full_atom_keys_in_file_order(full_input_pdb)
    key_to_pocket_idx = _pocket_key_to_index(pocket_pdb)
    variant_occ_table = _build_variant_occurrence_table(orig_keys_in_order)

    n_atoms_full = len(orig_keys_in_order)

    def _map_full_index_to_pocket(idx_one_based: int, stage_idx: int, tuple_idx: int, side_label: str) -> int:
        key = orig_keys_in_order[idx_one_based - 1]
        variant_occ = variant_occ_table[idx_one_based - 1]
        for variant in _key_variants(key):
            occurrence = variant_occ.get(variant)
            indices = key_to_pocket_idx.get(variant)
            if occurrence is None or not indices:
                continue
            if occurrence <= len(indices):
                return indices[occurrence - 1]

        msg_key = _format_atom_key_for_msg(key)
        raise click.BadParameter(
            f"--scan-lists #{stage_idx} tuple #{tuple_idx} ({side_label}) references atom index {idx_one_based} "
            f"(key {msg_key}) which is not present in the pocket after extraction. "
            "Increase extraction coverage (e.g., --radius/--radius-het2het, --selected_resn, or set --exclude-backbone False), "
            "or choose atoms that survive in the pocket."
        )

    converted: List[List[Tuple[int, int, float]]] = []
    for stage_idx, stage in enumerate(stages, start=1):
        stage_converted: List[Tuple[int, int, float]] = []
        for tuple_idx, (idx_i, idx_j, target) in enumerate(stage, start=1):
            if idx_i <= 0 or idx_j <= 0:
                raise click.BadParameter(
                    f"--scan-lists #{stage_idx} tuple #{tuple_idx} must use 1-based atom indices."
                )
            if idx_i > n_atoms_full or idx_j > n_atoms_full:
                raise click.BadParameter(
                    f"--scan-lists #{stage_idx} tuple #{tuple_idx} references an atom index "
                    f"beyond the input PDB atom count ({n_atoms_full})."
                )

            stage_converted.append(
                (
                    _map_full_index_to_pocket(idx_i, stage_idx, tuple_idx, "i"),
                    _map_full_index_to_pocket(idx_j, stage_idx, tuple_idx, "j"),
                    target,
                )
            )
        converted.append(stage_converted)
    return converted


def _round_charge_with_note(q: float) -> int:
    """
    Cast the extractor's total charge (float) to an integer suitable for the path search.
    If it is not already an integer within 1e-6, round to the nearest integer with a console note.
    """
    q_rounded = int(round(float(q)))
    if not math.isfinite(q):
        raise click.BadParameter(f"Computed total charge is non-finite: {q!r}")
    if abs(float(q) - q_rounded) > 1e-6:
        click.echo(f"[all] NOTE: extractor total charge = {q:g} → rounded to integer {q_rounded} for the path search.")
    return q_rounded


def _pdb_needs_elem_fix(p: Path) -> bool:
    """
    Return True if the PDB has at least one ATOM/HETATM record whose element field (cols 77–78) is empty.
    This is a light-weight check to decide whether to run add_elem_info.
    """
    try:
        with p.open("r", encoding="utf-8", errors="ignore") as fh:
            saw_atom = False
            for line in fh:
                if line.startswith("ATOM") or line.startswith("HETATM"):
                    saw_atom = True
                    if len(line) < 78 or not line[76:78].strip():
                        return True
        # If no ATOM/HETATM was seen, fall back to "no fix"
        return False
    except Exception:
        # On I/O errors, skip fixing (use original)
        return False


# ---------- Post-processing helpers (minimal, reuse internals) ----------

def _read_summary(summary_yaml: Path) -> List[Dict[str, Any]]:
    """
    Read path_search/summary.yaml and return segments list (empty if not found).
    """
    try:
        if not summary_yaml.exists():
            return []
        data = yaml.safe_load(summary_yaml.read_text(encoding="utf-8")) or {}
        segs = data.get("segments", []) or []
        if not isinstance(segs, list):
            return []
        return segs
    except Exception:
        return []


def _pdb_models_to_coords_and_elems(pdb_path: Path) -> Tuple[List[np.ndarray], List[str]]:
    """
    Return ([coords_model1, coords_model2, ...] in Å), [elements] from a multi-model PDB.
    """
    parser = PDB.PDBParser(QUIET=True)
    st = parser.get_structure("seg", str(pdb_path))
    models = list(st.get_models())
    if not models:
        raise click.ClickException(f"[post] No MODEL found in PDB: {pdb_path}")
    # atom order taken from first model
    atoms0 = [a for a in models[0].get_atoms()]
    elems: List[str] = []
    for a in atoms0:
        el = (a.element or "").strip()
        if not el:
            # fall back: derive from atom name
            nm = a.get_name().strip()
            el = "".join([c for c in nm if c.isalpha()])[:2].title() or "C"
        elems.append(el)
    coords_list: List[np.ndarray] = []
    for m in models:
        atoms = [a for a in m.get_atoms()]
        if len(atoms) != len(atoms0):
            raise click.ClickException(f"[post] Atom count mismatch across models in {pdb_path}")
        coords = np.array([a.get_coord() for a in atoms], dtype=float)
        coords_list.append(coords)
    return coords_list, elems


def _geom_from_angstrom(elems: Sequence[str],
                        coords_ang: np.ndarray,
                        freeze_atoms: Sequence[int]) -> Any:
    """
    Create a Geometry from Å coordinates using _path_search._new_geom_from_coords (expects Bohr).
    """
    coords_bohr = np.asarray(coords_ang, dtype=float) / BOHR2ANG
    return _path_search._new_geom_from_coords(elems, coords_bohr, coord_type="cart", freeze_atoms=freeze_atoms)


def _load_segment_end_geoms(seg_pdb: Path, freeze_atoms: Sequence[int]) -> Tuple[Any, Any]:
    """
    Load first/last model as Geometries from a per-segment pocket PDB.
    """
    coords_list, elems = _pdb_models_to_coords_and_elems(seg_pdb)
    gL = _geom_from_angstrom(elems, coords_list[0], freeze_atoms)
    gR = _geom_from_angstrom(elems, coords_list[-1], freeze_atoms)
    return gL, gR


def _compute_imag_mode_direction(ts_geom: Any,
                                 calc_kwargs: Dict[str, Any],
                                 freeze_atoms: Sequence[int],
                                 uma_kwargs: Optional[Dict[str, Any]] = None) -> np.ndarray:
    """
    Compute imaginary mode direction (N×3, unit vector in Cartesian space) at TS geometry.
    Uses tsopt internal helpers to minimize new code.
    """
    kw = calc_kwargs if calc_kwargs else uma_kwargs or {}
    device_str = kw.get("ml_device", kw.get("device", "auto"))
    if device_str == "auto":
        device_str = "cuda" if torch.cuda.is_available() else "cpu"
    # full analytic Hessian (torch tensor) — may be partial if B-factor layers are used
    H_t = _ts_opt._calc_full_hessian_torch(ts_geom, calc_kwargs=kw,
                                           device=torch.device(device_str))

    n_atoms_total = len(ts_geom.atoms)
    n_hess_atoms = H_t.shape[0] // 3  # partial Hessian: fewer than total atoms

    from ase.data import atomic_masses
    masses_amu = np.array([atomic_masses[z] for z in ts_geom.atomic_numbers])

    if n_hess_atoms < n_atoms_total:
        # Partial Hessian: restrict coords/masses to hess_active_atoms
        active_atoms = getattr(ts_geom, "_hess_active_atoms_last", None)
        if active_atoms is None or len(active_atoms) != n_hess_atoms:
            # Fallback: assume first n_hess_atoms (ML + Hessian-target MM) are active.
            active_atoms = np.arange(n_hess_atoms)
        coords_bohr_sub = ts_geom.coords.reshape(-1, 3)[active_atoms]
        masses_amu_sub = masses_amu[active_atoms]
        coords_bohr_t = torch.as_tensor(coords_bohr_sub, dtype=H_t.dtype, device=H_t.device)
        masses_au_t = torch.as_tensor(masses_amu_sub * _ts_opt.AMU2AU, dtype=H_t.dtype, device=H_t.device)
    else:
        coords_bohr_t = torch.as_tensor(ts_geom.coords.reshape(-1, 3), dtype=H_t.dtype, device=H_t.device)
        masses_au_t = torch.as_tensor(masses_amu * _ts_opt.AMU2AU, dtype=H_t.dtype, device=H_t.device)

    mode_sub = _ts_opt._mode_direction_by_root(H_t, coords_bohr_t, masses_au_t,
                                               root=0,
                                               freeze_idx=list(freeze_atoms) if len(freeze_atoms) > 0 else None)

    # If partial, embed mode back into full-atom space (non-active atoms get zero displacement)
    if n_hess_atoms < n_atoms_total:
        mode_full = np.zeros((n_atoms_total, 3), dtype=float)
        active_atoms_np = np.asarray(active_atoms)
        mode_full[active_atoms_np] = mode_sub.reshape(-1, 3)
        mode = mode_full
    else:
        mode = mode_sub

    # ensure unit length
    norm = float(np.linalg.norm(mode.reshape(-1)))
    if norm <= 0:
        raise click.ClickException("[post] Imaginary mode direction has zero norm.")
    return (mode / norm)


def _displaced_geometry_along_mode(geom: Any,
                                   mode_xyz: np.ndarray,
                                   amplitude_ang: float,
                                   freeze_atoms: Sequence[int]) -> Any:
    """
    Displace geometry along mode by ± amplitude (Å). Returns new Geometry.
    """
    coords_bohr = np.asarray(geom.coords3d, dtype=float)  # Bohr
    disp_bohr = (amplitude_ang / BOHR2ANG) * np.asarray(mode_xyz, dtype=float)  # (N,3)
    new_coords_bohr = coords_bohr + disp_bohr
    return _path_search._new_geom_from_coords(geom.atoms, new_coords_bohr, coord_type=geom.coord_type, freeze_atoms=freeze_atoms)


def _save_single_geom_as_pdb_for_tools(g: Any, ref_pdb: Path, out_dir: Path, name: str) -> Path:
    """
    Write a single-geometry XYZ/TRJ with energy and convert to PDB using the pocket ref (for downstream CLI tools).
    Returns PDB path.
    """
    out_dir.mkdir(parents=True, exist_ok=True)
    xyz_trj = out_dir / f"{name}.trj"
    _path_search._write_xyz_trj_with_energy([g], [float(g.energy)], xyz_trj)
    pdb_out = out_dir / f"{name}.pdb"
    _path_search._maybe_convert_to_pdb(xyz_trj, ref_pdb_path=ref_pdb, out_path=pdb_out)
    return pdb_out


def _run_tsopt_on_hei(hei_pdb: Path,
                      charge: int,
                      spin: int,
                      real_parm7: Path,
                      model_pdb: Path,
                      detect_layer: bool,
                      args_yaml: Optional[Path],
                      out_dir: Path,
                      opt_mode_default: str,
                      overrides: Optional[Dict[str, Any]] = None) -> Tuple[Path, Any]:
    """
    Run tsopt CLI on a HEI pocket PDB; return (final_ts_pdb_path, ts_geom)
    """
    overrides = overrides or {}
    prepared_input = prepare_input_structure(hei_pdb)
    ts_dir = _resolve_override_dir(out_dir / "ts", overrides.get("out_dir"))
    ensure_dir(ts_dir)

    opt_mode = overrides.get("opt_mode", opt_mode_default)

    ts_args: List[str] = [
        "-i", str(prepared_input.source_path),
        "--real-parm7", str(real_parm7),
        "--model-pdb", str(model_pdb),
        "-q", str(int(charge)),
        "-m", str(int(spin)),
        "--out-dir", str(ts_dir),
    ]
    ts_args.append("--detect-layer" if detect_layer else "--no-detect-layer")

    if opt_mode is not None:
        ts_args.extend(["--opt-mode", str(opt_mode)])

    _append_cli_arg(ts_args, "--max-cycles", overrides.get("max_cycles"))
    _append_cli_arg(ts_args, "--dump", overrides.get("dump"))

    hess_mode = overrides.get("hessian_calc_mode")
    if hess_mode:
        ts_args.extend(["--hessian-calc-mode", str(hess_mode)])

    if args_yaml is not None:
        ts_args.extend(["--args-yaml", str(args_yaml)])

    _echo(f"[tsopt] Running tsopt on HEI → out={ts_dir}")
    _run_cli_main("tsopt", _ts_opt.cli, ts_args, on_nonzero="raise", prefix="tsopt")

    # Prefer PDB (tsopt converts when input is PDB)
    ts_pdb = ts_dir / "final_geometry.pdb"
    final_xyz = ts_dir / "final_geometry.xyz"
    if not ts_pdb.exists():
        # fallback: use final .xyz and convert
        if not final_xyz.exists():
            prepared_input.cleanup()
            raise click.ClickException("[tsopt] TS outputs not found.")
        _path_search._maybe_convert_to_pdb(final_xyz, hei_pdb, ts_dir / "final_geometry.pdb")
    ts_pdb = ts_dir / "final_geometry.pdb"
    g_ts = geom_loader(ts_pdb, coord_type="cart")

    # Ensure calculator to have energy on g_ts
    calc = _mlmm_calc(
        model_charge=int(charge),
        model_mult=int(spin),
        input_pdb=str(ts_pdb),
        real_parm7=str(real_parm7),
        model_pdb=str(model_pdb),
        use_bfactor_layers=detect_layer,
    )
    g_ts.set_calculator(calc)
    _ = float(g_ts.energy)

    prepared_input.cleanup()
    return ts_pdb, g_ts


def _pseudo_irc_and_match(seg_idx: int,
                          seg_dir: Path,
                          ref_pdb_for_seg: Path,
                          seg_pocket_pdb: Path,
                          g_ts: Any,
                          q_int: int,
                          spin: int,
                          real_parm7: Optional[Path] = None,
                          model_pdb: Optional[Path] = None,
                          detect_layer: bool = False) -> Dict[str, Any]:
    """
    From a TS pocket geometry, perform pseudo-IRC:
      - compute imag. mode
      - displace ± (0.25 Å) and optimize both to minima (LBFGS)
      - map each min to left/right segment endpoint by bond-change check (if segment endpoints exist)
        *If no segment endpoints are available (single-structure TSOPT-only mode), fall back to (plus, minus).*
    Returns dict with paths/energies/geoms for {left, ts, right}, and small IRC plots.
    """
    # Mode direction
    calc_kwargs = dict(
        model_charge=int(q_int),
        model_mult=int(spin),
        input_pdb=str(ref_pdb_for_seg),
        real_parm7=str(real_parm7) if real_parm7 else None,
        model_pdb=str(model_pdb) if model_pdb else None,
        use_bfactor_layers=detect_layer,
    )
    mode_xyz = _compute_imag_mode_direction(g_ts, calc_kwargs=calc_kwargs, freeze_atoms=[])

    # Displace ± and optimize
    irc_dir = seg_dir / "irc"
    ensure_dir(irc_dir)
    amp = 0.25  # Å; small stable displacement
    g_plus0 = _displaced_geometry_along_mode(g_ts,  mode_xyz, +amp, [])
    g_minus0 = _displaced_geometry_along_mode(g_ts, mode_xyz, -amp, [])

    # Shared ML/MM calc
    shared_calc = _mlmm_calc(**calc_kwargs)
    # LBFGS settings (reuse defaults)
    sopt_cfg = dict(_path_search.LBFGS_KW)
    sopt_cfg["dump"] = True
    sopt_cfg["out_dir"] = str(irc_dir)

    # Optimize
    g_plus  = _path_search._optimize_single(g_plus0, shared_calc, sopt_cfg, irc_dir, tag=f"seg_{seg_idx:02d}_irc_plus",  ref_pdb_path=seg_pocket_pdb)
    g_minus = _path_search._optimize_single(g_minus0, shared_calc, sopt_cfg, irc_dir, tag=f"seg_{seg_idx:02d}_irc_minus", ref_pdb_path=seg_pocket_pdb)

    # IRC mini plots (TS→min)
    try:
        trj_plus  = irc_dir / f"seg_{seg_idx:02d}_irc_plus_opt/optimization.trj"
        trj_minus = irc_dir / f"seg_{seg_idx:02d}_irc_minus_opt/optimization.trj"
        if trj_plus.exists():
            run_trj2fig(trj_plus, [irc_dir / f"irc_plus_plot.png"], unit="kcal", reference="init", reverse_x=False)
        if trj_minus.exists():
            run_trj2fig(trj_minus, [irc_dir / f"irc_minus_plot.png"], unit="kcal", reference="init", reverse_x=False)
    except Exception as e:
        click.echo(f"[irc] WARNING: failed to plot IRC mini plots: {e}", err=True)

    # Try to load segment endpoints (pocket-only) — available only after path_search
    gL_end = None
    gR_end = None
    seg_pocket_path = seg_dir.parent / f"mep_seg_{seg_idx:02d}.pdb"
    if seg_pocket_path.exists():
        try:
            gL_end, gR_end = _load_segment_end_geoms(seg_pocket_path, [])
        except Exception as e:
            click.echo(f"[post] WARNING: failed to load segment endpoints: {e}", err=True)

    # Decide mapping
    bond_cfg = dict(_path_search.BOND_KW)

    def _rmsd_cart(g1, g2) -> float:
        """Cartesian RMSD (Bohr) between two geometries."""
        c1 = np.asarray(g1.coords).reshape(-1, 3)
        c2 = np.asarray(g2.coords).reshape(-1, 3)
        n = min(len(c1), len(c2))
        return float(np.sqrt(np.mean((c1[:n] - c2[:n]) ** 2)))

    def _matches(x, y) -> bool:
        try:
            chg, _ = _path_search._has_bond_change(x, y, bond_cfg)
            return (not chg)
        except Exception:
            # fallback: small RMSD threshold
            return (_rmsd_cart(x, y) < 1e-3)

    candidates = [("plus", g_plus), ("minus", g_minus)]
    mapping: Dict[str, Any] = {"left": None, "right": None}

    if (gL_end is not None) and (gR_end is not None):
        # First pass: exact match on bond changes
        for tag, g in candidates:
            if _matches(g, gL_end) and not _matches(g, gR_end):
                mapping["left"] = (tag, g)
            elif _matches(g, gR_end) and not _matches(g, gL_end):
                mapping["right"] = (tag, g)
        # Second pass: fill missing by RMSD
        for side, g_end in (("left", gL_end), ("right", gR_end)):
            if mapping[side] is None:
                remain = [(t, gg) for (t, gg) in candidates if mapping.get("left", (None, None))[0] != t and mapping.get("right", (None, None))[0] != t]
                if not remain:
                    remain = candidates
                best = min(remain, key=lambda p: _rmsd_cart(p[1], g_end))
                mapping[side] = best
    else:
        # Fallback (single-structure TSOPT-only mode): keep a deterministic assignment
        mapping["left"] = ("minus", g_minus)
        mapping["right"] = ("plus", g_plus)

    # Energies (ensure calculator is attached)
    for _, g in candidates:
        if g.calculator is None:
            g.set_calculator(shared_calc)
        _ = float(g.energy)
    if g_ts.calculator is None:
        g_ts.set_calculator(shared_calc)
    _ = float(g_ts.energy)

    # Dump tiny TS↔min trj for each direction
    try:
        for side in ("left", "right"):
            tag, gmin = mapping[side]
            trj = irc_dir / f"irc_{side}.trj"
            _path_search._write_xyz_trj_with_energy([g_ts, gmin], [float(g_ts.energy), float(gmin.energy)], trj)
            run_trj2fig(trj, [irc_dir / f"irc_{side}_plot.png"], unit="kcal", reference="init", reverse_x=False)
    except Exception:
        pass

    irc_trj = None
    irc_plot = None
    try:
        g_left = mapping["left"][1]
        g_right = mapping["right"][1]
        irc_trj = irc_dir / "finished_irc.trj"
        _path_search._write_xyz_trj_with_energy(
            [g_left, g_ts, g_right],
            [float(g_left.energy), float(g_ts.energy), float(g_right.energy)],
            irc_trj,
        )
        irc_plot = irc_dir / "irc_plot.png"
        run_trj2fig(irc_trj, [irc_plot], unit="kcal", reference="init", reverse_x=False)
    except Exception as e:
        click.echo(f"[irc] WARNING: failed to build combined IRC plot: {e}", err=True)

    return {
        "left_min_geom": mapping["left"][1],
        "right_min_geom": mapping["right"][1],
        "ts_geom": g_ts,
        "left_tag": mapping["left"][0],
        "right_tag": mapping["right"][0],
        "irc_trj": str(irc_trj) if irc_trj else None,
        "irc_plot": str(irc_plot) if irc_plot else None,
    }


def _write_segment_energy_diagram(
    prefix: Path,
    labels: List[str],
    energies_eh: List[float],
    title_note: str,
    ylabel: str = "ΔE (kcal/mol)",
    write_html: bool = False,
) -> Optional[Dict[str, Any]]:
    """
    Write energy diagram (PNG only) using utils.build_energy_diagram.
    """
    if not energies_eh:
        return None
    e0 = energies_eh[0]
    energies_kcal = [(e - e0) * AU2KCALPERMOL for e in energies_eh]
    fig = build_energy_diagram(
        energies=energies_kcal,
        labels=labels,
        ylabel=ylabel,
        baseline=True,
        showgrid=False,
    )
    if title_note:
        fig.update_layout(title=title_note)
    png = prefix.with_suffix(".png")
    try:
        fig.write_image(str(png), scale=2)
    except Exception as e:
        click.echo(f"[diagram] NOTE: PNG export skipped (install 'kaleido' to enable): {e}", err=True)
    else:
            click.echo(f"[diagram] Wrote energy diagram → {png.name}")

    payload: Dict[str, Any] = {
        "name": prefix.stem,
        "labels": labels,
        "energies_kcal": energies_kcal,
        "ylabel": ylabel,
        "energies_au": list(energies_eh),
        "image": str(png),
    }
    if title_note:
        payload["title"] = title_note
    return payload


def _build_global_segment_labels(n_segments: int) -> List[str]:
    """
    Build GSM-like labels for aggregated R/TS/P diagrams over multiple segments.

    Pattern:
      - n = 1: ["R", "TS1", "P"]
      - n >= 2: R, TS1, IM1_1, IM1_2, TS2, IM2_1, IM2_2, ..., TSN, P
    """
    if n_segments <= 0:
        return []
    if n_segments == 1:
        return ["R", "TS1", "P"]
    labels = ["R"]
    for seg_idx in range(1, n_segments + 1):
        if seg_idx == 1:
            labels.extend([f"TS{seg_idx}", f"IM{seg_idx}_1", f"IM{seg_idx}_2"])
        elif seg_idx == n_segments:
            labels.extend([f"IM{seg_idx - 1}_2", f"TS{seg_idx}", "P"])
        else:
            labels.extend(
                [f"IM{seg_idx - 1}_2", f"TS{seg_idx}", f"IM{seg_idx}_1"]
            )
    return labels


def _merge_irc_trajectories_to_single_plot(
    trj_and_flags: Sequence[Tuple[Path, bool]],
    out_png: Path,
) -> None:
    """
    Build a single IRC plot over all reactive segments using trj2fig.
    """
    all_blocks: List[str] = []
    for trj_path, reverse in trj_and_flags:
        if not isinstance(trj_path, Path) or not trj_path.exists():
            continue
        try:
            blocks = read_xyz_as_blocks(trj_path)
        except click.ClickException as e:
            click.echo(str(e), err=True)
            continue
        if not blocks:
            continue
        if reverse:
            blocks = list(reversed(blocks))
        all_blocks.extend("\n".join(b) for b in blocks)

    if not all_blocks:
        return

    tmp_trj = out_png.with_suffix(".trj")
    ensure_dir(tmp_trj.parent)
    try:
        tmp_trj.write_text("\n".join(all_blocks) + "\n", encoding="utf-8")
    except Exception as e:
        click.echo(f"[irc_all] WARNING: Failed to write concatenated IRC trajectory: {e}", err=True)
        return

    try:
        run_trj2fig(tmp_trj, [out_png], unit="kcal", reference="init", reverse_x=False)
        click.echo(f"[irc_all] Wrote aggregated IRC plot → {out_png}")
    except Exception as e:
        click.echo(f"[irc_all] WARNING: failed to plot concatenated IRC trajectory: {e}", err=True)
    finally:
        try:
            tmp_trj.unlink()
        except Exception:
            pass


def _run_freq_for_state(pdb_path: Path,
                        q_int: int,
                        spin: int,
                        real_parm7: Path,
                        model_pdb: Path,
                        detect_layer: bool,
                        out_dir: Path,
                        args_yaml: Optional[Path],
                        overrides: Optional[Dict[str, Any]] = None) -> Dict[str, Any]:
    """
    Run freq CLI; return parsed thermo dict (may be empty).
    """
    fdir = out_dir
    ensure_dir(fdir)
    overrides = overrides or {}

    dump_use = overrides.get("dump")
    if dump_use is None:
        dump_use = True

    args = [
        "-i", str(pdb_path),
        "--real-parm7", str(real_parm7),
        "--model-pdb", str(model_pdb),
        "-q", str(int(q_int)),
        "-m", str(int(spin)),
        "--out-dir", str(fdir),
    ]
    args.append("--detect-layer" if detect_layer else "--no-detect-layer")

    _append_cli_arg(args, "--max-write", overrides.get("max_write"))
    _append_cli_arg(args, "--amplitude-ang", overrides.get("amplitude_ang"))
    _append_cli_arg(args, "--n-frames", overrides.get("n_frames"))
    if overrides.get("sort") is not None:
        args.extend(["--sort", str(overrides.get("sort"))])
    _append_cli_arg(args, "--temperature", overrides.get("temperature"))
    _append_cli_arg(args, "--pressure", overrides.get("pressure"))
    _append_cli_arg(args, "--dump", dump_use)

    hess_mode = overrides.get("hessian_calc_mode")
    if hess_mode:
        args.extend(["--hessian-calc-mode", str(hess_mode)])

    if args_yaml is not None:
        args.extend(["--args-yaml", str(args_yaml)])
    _run_cli_main("freq", _freq_cli.cli, args, on_nonzero="warn", on_exception="raise", prefix="freq")
    # parse thermoanalysis.yaml if any
    y = fdir / "thermoanalysis.yaml"
    if y.exists():
        try:
            return yaml.safe_load(y.read_text(encoding="utf-8")) or {}
        except Exception:
            return {}
    return {}


def _run_dft_for_state(pdb_path: Path,
                       q_int: int,
                       spin: int,
                       real_parm7: Path,
                       model_pdb: Path,
                       detect_layer: bool,
                       out_dir: Path,
                       args_yaml: Optional[Path],
                       func_basis: str = "wb97x-v/def2-tzvp",
                       overrides: Optional[Dict[str, Any]] = None) -> Dict[str, Any]:
    """
    Run dft CLI; return parsed result.yaml dict (may be empty).
    """
    ddir = out_dir
    ensure_dir(ddir)
    overrides = overrides or {}

    func_basis_use = overrides.get("func_basis", func_basis)

    args = [
        "-i", str(pdb_path),
        "--real-parm7", str(real_parm7),
        "--model-pdb", str(model_pdb),
        "-q", str(int(q_int)),
        "-m", str(int(spin)),
        "--func-basis", str(func_basis_use),
        "--out-dir", str(ddir),
    ]
    args.append("--detect-layer" if detect_layer else "--no-detect-layer")

    _append_cli_arg(args, "--max-cycle", overrides.get("max_cycle"))
    _append_cli_arg(args, "--conv-tol", overrides.get("conv_tol"))
    _append_cli_arg(args, "--grid-level", overrides.get("grid_level"))

    if args_yaml is not None:
        args.extend(["--args-yaml", str(args_yaml)])
    _run_cli_main("dft", _dft_cli.cli, args, on_nonzero="warn", on_exception="raise", prefix="dft")
    y = out_dir / "result.yaml"
    if y.exists():
        try:
            return yaml.safe_load(y.read_text(encoding="utf-8")) or {}
        except Exception:
            return {}
    return {}


# -----------------------------
# CLI
# -----------------------------

_ALL_PRIMARY_HELP_OPTIONS = frozenset(
    {
        "-i",
        "--input",
        "-c",
        "--center",
        "--ligand-charge",
        "--out-dir",
        "--tsopt",
        "--thermo",
        "--dft",
        "--config",
        "--dry-run",
        "--help-advanced",
    }
)


def _show_advanced_help(
    ctx: click.Context, _param: click.Parameter, value: bool
) -> None:
    """Print full option help (including hidden advanced options) and exit."""
    if not value or ctx.resilient_parsing:
        return

    hidden = getattr(ctx.command, "_advanced_hidden_options", ())
    restored: list[click.Option] = []
    for opt in hidden:
        if opt.hidden:
            opt.hidden = False
            restored.append(opt)
    try:
        click.echo(ctx.command.get_help(ctx))
    finally:
        for opt in restored:
            opt.hidden = True
    ctx.exit()


def _configure_all_help_visibility(command: click.Command) -> None:
    """Hide advanced options from default --help while keeping them functional."""
    hidden_options: list[click.Option] = []
    for param in command.params:
        if not isinstance(param, click.Option):
            continue
        names = set(param.opts + param.secondary_opts)
        if names & _ALL_PRIMARY_HELP_OPTIONS:
            continue
        if param.hidden:
            continue
        param.hidden = True
        hidden_options.append(param)
    setattr(command, "_advanced_hidden_options", tuple(hidden_options))


@click.command(
    help="Run pocket extraction → (optional single-structure staged scan) → MEP search → merge to full PDBs in one shot.\n"
         "If exactly one input is provided: (a) with --scan-lists, stage results feed into path_search; "
         "(b) with --tsopt True and no --scan-lists, run TSOPT-only mode.",
    context_settings={
        "help_option_names": ["-h", "--help"],
        "ignore_unknown_options": True,
        "allow_extra_args": True,
    },
)
@click.option(
    "--help-advanced",
    is_flag=True,
    is_eager=True,
    expose_value=False,
    callback=_show_advanced_help,
    help="Show all options (including advanced settings) and exit.",
)
# ===== Inputs =====
@click.option(
    "-i", "--input", "input_paths",
    type=click.Path(path_type=Path, exists=True, dir_okay=False),
    multiple=True, required=True,
    help=("Two or more **full** PDBs in reaction order (reactant [intermediates ...] product), "
          "or a single **full** PDB (with --scan-lists or with --tsopt True). "
          "You may pass a single '-i' followed by multiple space-separated files (e.g., '-i A.pdb B.pdb C.pdb').")
)
@click.option(
    "-c", "--center", "center_spec",
    type=str, required=True,
    help=("Substrate specification for the extractor: "
          "a PDB path, a residue-ID list like '123,124' or 'A:123,B:456' "
          "(insertion codes OK: '123A' / 'A:123A'), "
          "or a residue-name list like 'GPP,MMT'.")
)
@click.option(
    "--out-dir", "out_dir",
    type=click.Path(path_type=Path, file_okay=False),
    default=Path("./result_all/"), show_default=True,
    help="Top-level output directory for the pipeline."
)
# ===== Extractor knobs (subset of extract.parse_args) =====
@click.option("-r", "--radius", type=float, default=2.6, show_default=True,
              help="Inclusion cutoff (Å) around substrate atoms.")
@click.option("--radius-het2het", type=float, default=0.0, show_default=True,
              help="Independent hetero–hetero cutoff (Å) for non‑C/H pairs.")
@click.option("--include-H2O", "--include-h2o", "include_h2o", type=click.BOOL, default=True, show_default=True,
              help="Include waters (HOH/WAT/TIP3/SOL) in the pocket.")
@click.option("--exclude-backbone", "exclude_backbone", type=click.BOOL, default=True, show_default=True,
              help="Remove backbone atoms on non‑substrate amino acids (with PRO/HYP safeguards).")
@click.option("--add-linkH", "add_linkh", type=click.BOOL, default=False, show_default=True,
              help="Add link hydrogens for severed bonds (carbon-only) in pockets.")
@click.option("--selected_resn", type=str, default="", show_default=True,
              help="Force-include residues (comma/space separated; chain/insertion codes allowed).")
@click.option("--ligand-charge", type=str, default=None,
              help=("Either a total charge (number) to distribute across unknown residues "
                    "or a mapping like 'GPP:-3,MMT:-1'."))
@click.option("--auto-mm-ff-set", "mm_ff_set",
              type=click.Choice(["ff19SB", "ff14SB"], case_sensitive=False),
              default="ff19SB", show_default=True,
              help="Force-field set forwarded to mm_parm (ff19SB uses OPC3; ff14SB uses TIP3P).")
@click.option("--auto-mm-add-ter/--auto-mm-no-add-ter", "mm_add_ter",
              default=True, show_default=True,
              help="Control mm_parm TER insertion around ligand/water/ion blocks.")
@click.option("--auto-mm-keep-temp", "mm_keep_temp", is_flag=True, default=False, show_default=True,
              help="Keep the mm_parm temporary working directory (for debugging).")
@click.option(
    "--auto-mm-ligand-mult",
    "mm_ligand_mult",
    type=str,
    default=None,
    help=("Spin multiplicity mapping forwarded to mm_parm (e.g., 'GPP:2,SAM:1'). "
          "If omitted, mm_parm defaults to 1 for all ligands.")
)
@click.option(
    "--auto-mm-allow-nonstandard-aa",
    "mm_allow_nonstandard_aa",
    is_flag=True,
    default=False,
    show_default=True,
    help="Allow mm_parm to parameterize amino-acid-like residues via antechamber.",
)
@click.option("--verbose", type=click.BOOL, default=True, show_default=True, help="Enable INFO-level logging inside extractor.")
# ===== Path search knobs (subset of path_search.cli) =====
@click.option("-m", "--multiplicity", "spin", type=int, default=1, show_default=True, help="Multiplicity (2S+1).")
@click.option("--max-nodes", type=int, default=10, show_default=True,
              help="Max internal nodes for **segment** GSM (String has max_nodes+2 images including endpoints).")
@click.option("--max-cycles", type=int, default=300, show_default=True, help="Maximum GSM optimization cycles.")
@click.option("--climb", type=click.BOOL, default=True, show_default=True,
              help="Enable transition-state climbing after growth for the **first** segment in each pair.")
@click.option("--sopt-mode", type=click.Choice(["lbfgs", "rfo", "light", "heavy"], case_sensitive=False),
              default="lbfgs", show_default=True,
              help="Single-structure optimizer kind for HEI±1 and kink nodes.")
@click.option("--opt-mode", type=str, default=None,
              help="Common optimizer mode forwarded to scan/tsopt (--opt-mode). When unset, tools use their defaults.")
@click.option("--dump", type=click.BOOL, default=False, show_default=True,
              help="Dump GSM / single-structure trajectories during the run, forwarding the same flag to scan/tsopt/freq.")
@click.option("--config", "config_yaml", type=click.Path(path_type=Path, exists=True, dir_okay=False),
              default=None, help="Base YAML configuration file applied before explicit CLI options.")
@click.option("--override-yaml", type=click.Path(path_type=Path, exists=True, dir_okay=False),
              default=None, help="Final YAML override file (highest priority YAML layer).")
@click.option("--args-yaml", "args_yaml_legacy", type=click.Path(path_type=Path, exists=True, dir_okay=False),
              default=None, help="[legacy] Alias of --override-yaml; kept for backward compatibility.")
@click.option("--show-config/--no-show-config", "show_config", default=False, show_default=True,
              help="Print resolved configuration and continue execution.")
@click.option("--dry-run/--no-dry-run", "dry_run", default=False, show_default=True,
              help="Validate options and print the execution plan without running any stage.")
@click.option("--pre-opt", "pre_opt", type=click.BOOL, default=True, show_default=True,
              help="If False, skip initial single-structure optimizations of the pocket inputs.")
@click.option("--hessian-calc-mode",
              type=click.Choice(["Analytical", "FiniteDifference"], case_sensitive=False),
              default=None,
              help="Common UMA Hessian calculation mode forwarded to tsopt and freq.")
@click.option(
    "--detect-layer/--no-detect-layer",
    "detect_layer",
    default=True,
    show_default=True,
    help="Detect ML/MM layers from input PDB B-factors (B=0/10/20) in downstream tools. "
         "If disabled, downstream tools require --model-pdb or --model-indices.",
)
# ===== Post-processing toggles =====
@click.option("--tsopt", "do_tsopt", type=click.BOOL, default=False, show_default=True,
              help="TS optimization + pseudo-IRC per reactive segment (or TSOPT-only mode for single-structure), and build energy diagrams.")
@click.option("--thermo", "do_thermo", type=click.BOOL, default=False, show_default=True,
              help="Run freq on (R,TS,P) per reactive segment (or TSOPT-only mode) and build Gibbs free-energy diagram (UMA).")
@click.option("--dft", "do_dft", type=click.BOOL, default=False, show_default=True,
              help="Run DFT single-point on (R,TS,P) and build DFT energy diagram. With --thermo True, also generate a DFT//UMA Gibbs diagram.")
@click.option("--tsopt-mode", type=click.Choice(["light", "heavy"], case_sensitive=False),
              default="light", show_default=True,
              help="TS optimizer mode: light (Dimer) or heavy (RS-I-RFO).")
@click.option("--tsopt-max-cycles", type=int, default=None,
              help="Override tsopt --max-cycles value.")
@click.option("--tsopt-out-dir", type=click.Path(path_type=Path, file_okay=False), default=None,
              help="Override tsopt output subdirectory (relative paths are resolved against the default).")
@click.option("--freq-out-dir", type=click.Path(path_type=Path, file_okay=False), default=None,
              help="Override freq output base directory (relative paths resolved against the default).")
@click.option("--freq-max-write", type=int, default=None,
              help="Override freq --max-write value.")
@click.option("--freq-amplitude-ang", type=float, default=None,
              help="Override freq --amplitude-ang (Å).")
@click.option("--freq-n-frames", type=int, default=None,
              help="Override freq --n-frames value.")
@click.option("--freq-sort", type=click.Choice(["value", "abs"], case_sensitive=False), default=None,
              help="Override freq mode sorting.")
@click.option("--freq-temperature", type=float, default=None,
              help="Override freq thermochemistry temperature (K).")
@click.option("--freq-pressure", type=float, default=None,
              help="Override freq thermochemistry pressure (atm).")
@click.option("--dft-out-dir", type=click.Path(path_type=Path, file_okay=False), default=None,
              help="Override dft output base directory (relative paths resolved against the default).")
@click.option("--dft-func-basis", type=str, default=None,
              help="Override dft --func-basis value.")
@click.option("--dft-max-cycle", type=int, default=None,
              help="Override dft --max-cycle value.")
@click.option("--dft-conv-tol", type=float, default=None,
              help="Override dft --conv-tol value.")
@click.option("--dft-grid-level", type=int, default=None,
              help="Override dft --grid-level value.")
# ===== NEW: staged scan specification for single-structure route =====
@click.option(
    "--scan-lists",
    "--scan-list",
    "scan_lists_raw",
    type=str, multiple=True, required=False,
    help='Python-like list of (i,j,target_Å) per stage for **single-structure** scan. '
         'Pass a single --scan-list(s) followed by multiple literals, e.g. '
         '"[(12,45,1.35)]" "[(10,55,2.20),(23,34,1.80)]". '
         'Indices refer to the original full PDB (1-based) or PDB atom selectors like "TYR,285,CA"; '
         'they are auto-mapped to the pocket after extraction. Stage results feed into path_search.',
)
@click.option("--scan-out-dir", type=click.Path(path_type=Path, file_okay=False), default=None,
              help="Override the scan output directory (default: <out-dir>/scan/). Relative paths are resolved against the default parent.")
@click.option("--scan-one-based", type=click.BOOL, default=None,
              help="Override scan indexing interpretation (True = 1-based, False = 0-based).")
@click.option("--scan-max-step-size", type=float, default=None,
              help="Override scan --max-step-size (Å).")
@click.option("--scan-bias-k", type=float, default=None,
              help="Override scan harmonic bias strength k (eV/Å^2).")
@click.option("--scan-relax-max-cycles", type=int, default=None,
              help="Override scan relaxation max cycles per step.")
@click.option("--scan-preopt", "scan_preopt_override", type=click.BOOL, default=None,
              help="Override scan --preopt flag.")
@click.option("--scan-endopt", "scan_endopt_override", type=click.BOOL, default=None,
              help="Override scan --endopt flag.")
@click.pass_context
def cli(
    ctx: click.Context,
    input_paths: Sequence[Path],
    center_spec: str,
    out_dir: Path,
    radius: float,
    radius_het2het: float,
    include_h2o: bool,
    exclude_backbone: bool,
    add_linkh: bool,
    selected_resn: str,
    ligand_charge: Optional[str],
    mm_ff_set: str,
    mm_add_ter: bool,
    mm_keep_temp: bool,
    mm_ligand_mult: Optional[str],
    mm_allow_nonstandard_aa: bool,
    verbose: bool,
    spin: int,
    max_nodes: int,
    max_cycles: int,
    climb: bool,
    sopt_mode: str,
    opt_mode: Optional[str],
    dump: bool,
    config_yaml: Optional[Path],
    override_yaml: Optional[Path],
    args_yaml_legacy: Optional[Path],
    show_config: bool,
    dry_run: bool,
    pre_opt: bool,
    hessian_calc_mode: Optional[str],
    detect_layer: bool,
    do_tsopt: bool,
    do_thermo: bool,
    do_dft: bool,
    scan_lists_raw: Sequence[str],
    scan_out_dir: Optional[Path],
    scan_one_based: Optional[bool],
    scan_max_step_size: Optional[float],
    scan_bias_k: Optional[float],
    scan_relax_max_cycles: Optional[int],
    scan_preopt_override: Optional[bool],
    scan_endopt_override: Optional[bool],
    tsopt_mode: str,
    tsopt_max_cycles: Optional[int],
    tsopt_out_dir: Optional[Path],
    freq_out_dir: Optional[Path],
    freq_max_write: Optional[int],
    freq_amplitude_ang: Optional[float],
    freq_n_frames: Optional[int],
    freq_sort: Optional[str],
    freq_temperature: Optional[float],
    freq_pressure: Optional[float],
    dft_out_dir: Optional[Path],
    dft_func_basis: Optional[str],
    dft_max_cycle: Optional[int],
    dft_conv_tol: Optional[float],
    dft_grid_level: Optional[int],
) -> None:
    """
    The **all** command composes `extract` → (optional `scan` on pocket) → `path_search` and hides ref-template bookkeeping.
    It also accepts the sloppy `-i A B C` style like `path_search` does. With single input:
      - with --scan-lists: run staged scan on the pocket and use stage results as inputs for path_search,
      - with --tsopt True and no --scan-lists: run TSOPT-only mode (no path_search).
    """
    global _log_started
    _log_started = False

    time_start = time.perf_counter()
    command_str = "mlmm all " + " ".join(sys.argv[1:])

    dump_override_requested = False
    try:
        dump_source = ctx.get_parameter_source("dump")
        dump_override_requested = dump_source not in (None, ParameterSource.DEFAULT)
    except Exception:
        dump_override_requested = False

    config_yaml, override_yaml, used_legacy_yaml = _resolve_yaml_sources(
        config_yaml=config_yaml,
        override_yaml=override_yaml,
        args_yaml_legacy=args_yaml_legacy,
    )
    if used_legacy_yaml:
        click.echo(
            "[deprecation] --args-yaml is deprecated; use --override-yaml.",
            err=True,
        )
    args_yaml, merged_yaml_cfg = _build_effective_args_yaml(
        config_yaml=config_yaml,
        override_yaml=override_yaml,
        tmp_prefix="mlmm_all_merged_",
    )

    mm_ff_set = "ff14SB" if str(mm_ff_set).lower().startswith("ff14") else "ff19SB"

    # --- Robustly accept a single "-i" followed by multiple paths (like path_search.cli) ---
    argv_all = sys.argv[1:]
    i_vals = collect_single_option_values(argv_all, ("-i", "--input"), label="-i/--input")
    if i_vals:
        i_parsed = validate_existing_files(
            i_vals,
            option_name="-i/--input",
            hint="When using '-i', list only existing file paths (multiple paths may follow a single '-i').",
        )
        input_paths = tuple(i_parsed)

    scan_vals = collect_single_option_values(
        argv_all, ("--scan-lists", "--scan-list"), "--scan-list(s)"
    )
    if scan_vals:
        scan_lists_raw = tuple(scan_vals)

    # --------------------------
    # Validate input count / single-structure modes
    # --------------------------
    is_single = (len(input_paths) == 1)
    has_scan = bool(scan_lists_raw)
    single_tsopt_mode = (is_single and (not has_scan) and do_tsopt)

    if (len(input_paths) < 2) and (not (is_single and (has_scan or do_tsopt))):
        raise click.BadParameter(
            "Provide at least two PDBs with -i/--input in reaction order, "
            "or use a single PDB with --scan-lists, or a single PDB with --tsopt True."
        )

    tsopt_opt_mode_default = opt_mode.lower() if opt_mode else "light"
    tsopt_overrides: Dict[str, Any] = {}
    if tsopt_max_cycles is not None:
        tsopt_overrides["max_cycles"] = int(tsopt_max_cycles)
    if dump_override_requested:
        tsopt_overrides["dump"] = bool(dump)
    if tsopt_out_dir is not None:
        tsopt_overrides["out_dir"] = tsopt_out_dir
    if hessian_calc_mode is not None:
        tsopt_overrides["hessian_calc_mode"] = hessian_calc_mode
    if opt_mode is not None:
        tsopt_overrides["opt_mode"] = tsopt_opt_mode_default

    freq_overrides: Dict[str, Any] = {}
    if freq_max_write is not None:
        freq_overrides["max_write"] = int(freq_max_write)
    if freq_amplitude_ang is not None:
        freq_overrides["amplitude_ang"] = float(freq_amplitude_ang)
    if freq_n_frames is not None:
        freq_overrides["n_frames"] = int(freq_n_frames)
    if freq_sort is not None:
        freq_overrides["sort"] = freq_sort.lower()
    if freq_temperature is not None:
        freq_overrides["temperature"] = float(freq_temperature)
    if freq_pressure is not None:
        freq_overrides["pressure"] = float(freq_pressure)
    if dump_override_requested:
        freq_overrides["dump"] = bool(dump)
    if hessian_calc_mode is not None:
        freq_overrides["hessian_calc_mode"] = hessian_calc_mode

    dft_overrides: Dict[str, Any] = {}
    if dft_max_cycle is not None:
        dft_overrides["max_cycle"] = int(dft_max_cycle)
    if dft_conv_tol is not None:
        dft_overrides["conv_tol"] = float(dft_conv_tol)
    if dft_grid_level is not None:
        dft_overrides["grid_level"] = int(dft_grid_level)

    dft_func_basis_use = dft_func_basis or "wb97m-v/6-31g**"
    dft_method_fallback = str(dft_overrides.get("func_basis", dft_func_basis_use))

    if show_config or dry_run:
        config_payload: Dict[str, Any] = {
            "yaml": {
                "config": str(config_yaml) if config_yaml else None,
                "override_yaml": str(override_yaml) if override_yaml else None,
                "effective_args_yaml": str(args_yaml) if args_yaml else None,
            },
            "all": {
                "inputs": [str(p) for p in input_paths],
                "center": center_spec,
                "out_dir": str(out_dir),
                "spin": int(spin),
                "max_nodes": int(max_nodes),
                "max_cycles": int(max_cycles),
                "climb": bool(climb),
                "sopt_mode": str(sopt_mode),
                "opt_mode": (None if opt_mode is None else str(opt_mode)),
                "dump": bool(dump),
                "pre_opt": bool(pre_opt),
                "detect_layer": bool(detect_layer),
                "tsopt": bool(do_tsopt),
                "thermo": bool(do_thermo),
                "dft": bool(do_dft),
            },
            "overrides": {
                "tsopt": tsopt_overrides,
                "freq": freq_overrides,
                "dft": dft_overrides,
            },
        }
        if merged_yaml_cfg:
            config_payload["effective_yaml"] = merged_yaml_cfg
        _echo_section("=== [all] Effective configuration ===")
        click.echo(
            yaml.safe_dump(config_payload, sort_keys=False, allow_unicode=True).rstrip()
        )

    if dry_run:
        _echo("[all] Dry-run mode: no extraction/search/post-processing was executed.")
        _echo(
            "[all] Planned stages: extract -> mm_parm -> optional scan -> path_search -> optional tsopt/freq/dft."
        )
        _echo(format_elapsed("[all] Elapsed for Whole Pipeline", time_start))
        return

    # --------------------------
    # Prepare directories
    # --------------------------
    out_dir = out_dir.resolve()
    pockets_dir = out_dir / "pockets"
    path_dir = out_dir / "path_search"
    scan_dir = _resolve_override_dir(out_dir / "scan", scan_out_dir)  # for single-structure scan mode
    ensure_dir(out_dir)
    ensure_dir(pockets_dir)
    if not single_tsopt_mode:
        ensure_dir(path_dir)  # path_search might be skipped only in tsopt-only mode

    # --------------------------
    # Preflight: add_elem_info only for inputs lacking element fields
    # → Create fixed copies under a temporary folder inside out_dir (used ONLY for extraction)
    # --------------------------
    elem_tmp_dir = out_dir / "add_elem_info"
    inputs_for_extract: List[Path] = []
    elem_fix_echo=False
    for p in input_paths:
        if _pdb_needs_elem_fix(p):
            if elem_fix_echo==False:
                _echo_section("=== [all] Preflight — add_elem_info (only when element fields are missing) ===")
                elem_fix_echo=True
            ensure_dir(elem_tmp_dir)
            out_p = (elem_tmp_dir / p.name).resolve()
            try:
                _assign_elem_info(str(p), str(out_p), overwrite=False)
                _echo(f"[all] add_elem_info: fixed elements → {out_p}")
                inputs_for_extract.append(out_p)
            except SystemExit as e:
                code = getattr(e, "code", 1)
                _echo(f"[all] WARNING: add_elem_info exited with code {code} for {p}; using original.", err=True)
                inputs_for_extract.append(p.resolve())
            except Exception as e:
                _echo(f"[all] WARNING: add_elem_info failed for {p}: {e} — using original file.", err=True)
                inputs_for_extract.append(p.resolve())
        else:
            inputs_for_extract.append(p.resolve())

    extract_inputs = tuple(inputs_for_extract)

    _echo_section("=== [all] Stage 1/3 — Active-site pocket extraction (multi-structure union when applicable) ===")

    # Build per-structure pocket output file list (one per input full PDB for extraction)
    pocket_outputs: List[Path] = []
    for p in extract_inputs:
        pocket_outputs.append((pockets_dir / f"pocket_{p.stem}.pdb").resolve())

    # Run extractor via its public API (multi-structure union mode)
    try:
        ex_res = extract_api(
            complex_pdb=[str(p) for p in extract_inputs],
            center=center_spec,
            output=[str(p) for p in pocket_outputs],
            radius=float(radius),
            radius_het2het=float(radius_het2het),
            include_H2O=bool(include_h2o),
            exclude_backbone=bool(exclude_backbone),
            add_linkH=bool(add_linkh),
            selected_resn=selected_resn or "",
            ligand_charge=ligand_charge,
            verbose=bool(verbose),
        )
    except Exception as e:
        raise click.ClickException(f"[all] Extractor failed: {e}")

    # Report extractor outputs and charge breakdown
    _echo("[all] Pocket files:")
    for op in pocket_outputs:
        _echo(f"  - {op}")

    try:
        cs = ex_res.get("charge_summary", {})
        q_total = float(cs.get("total_charge", 0.0))
        q_prot = float(cs.get("protein_charge", 0.0))
        q_lig = float(cs.get("ligand_total_charge", 0.0))
        q_ion = float(cs.get("ion_total_charge", 0.0))
        _echo("\n[all] Charge summary from extractor (model #1):")
        _echo(f"  Protein: {q_prot:+g},  Ligand: {q_lig:+g},  Ions: {q_ion:+g},  Total: {q_total:+g}")
        q_int = _round_charge_with_note(q_total)
    except Exception as e:
        raise click.ClickException(f"[all] Could not obtain total charge from extractor: {e}")

    # --------------------------
    # Stage 1b: ML-region definition (copy first pocket) and mm_parm on the first full input
    # --------------------------
    _echo_section("=== [all] ML/MM preparation — ML region + parm7 via mm_parm ===")
    first_pocket = pocket_outputs[0]
    first_full_input = extract_inputs[0]
    ml_region_pdb = _write_ml_region_definition(first_pocket, out_dir / "ml_region.pdb")
    _echo(f"[all] ML region definition → {ml_region_pdb}")
    _echo(f"[all] mm_parm source PDB → {first_full_input}")

    mm_dir = out_dir / "mm_parm"
    ensure_commands_available(
        ("tleap", "antechamber", "parmchk2"),
        context="mm_parm (AmberTools)",
    )
    real_parm7_path, real_rst7_path = _build_mm_parm7(
        pdb=first_full_input,
        ligand_charge_expr=ligand_charge,
        ligand_mult_expr=mm_ligand_mult,
        out_dir=mm_dir,
        ff_set=mm_ff_set,
        add_TER=mm_add_ter,
        keep_temp=mm_keep_temp,
        allow_nonstandard_aa=mm_allow_nonstandard_aa,
    )
    _echo(f"[all] mm_parm outputs → parm7: {real_parm7_path.name}, rst7: {real_rst7_path.name}")

    # --------------------------
    # define-layer: assign 3-layer B-factors to each full-system PDB
    # --------------------------
    _echo_section("=== [all] define-layer — assign 3-layer B-factors to full-system PDBs ===")
    layered_dir = out_dir / "layered"
    ensure_dir(layered_dir)
    layered_inputs: List[Path] = []
    for idx, full_pdb in enumerate(extract_inputs):
        out_layered = layered_dir / f"{full_pdb.stem}_layered.pdb"
        try:
            layer_info = _define_layers(
                input_pdb=full_pdb,
                output_pdb=out_layered,
                model_pdb=ml_region_pdb,
            )
            _echo(f"[all] define-layer [{idx}]: {full_pdb.name} → {out_layered.name}  "
                  f"(ML={len(layer_info.get('ml_indices', []))}, "
                  f"MovableMM={len(layer_info.get('movable_mm_indices', []))}, "
                  f"FrozenMM={len(layer_info.get('frozen_indices', []))})")
            layered_inputs.append(out_layered)
        except Exception as e:
            _echo(f"[all] WARNING: define-layer failed for {full_pdb.name}: {e}", err=True)
            _echo(f"[all] Falling back to original PDB (no B-factor layers).", err=True)
            layered_inputs.append(full_pdb)

    # --------------------------
    # Other path: single-structure + --tsopt True (and NO scan-lists) → TSOPT-only mode
    # --------------------------
    if single_tsopt_mode:
        _echo_section("=== [all] TSOPT-only single-structure mode ===")
        tsroot = out_dir / "tsopt_single"
        ensure_dir(tsroot)

        # Use the layered full-system PDB as TS initial guess
        layered_pdb = layered_inputs[0]
        # TS optimization
        ts_pdb, g_ts = _run_tsopt_on_hei(
            layered_pdb,
            q_int,
            spin,
            real_parm7_path,
            ml_region_pdb,
            detect_layer,
            args_yaml,
            tsroot,
            tsopt_opt_mode_default,
            overrides=tsopt_overrides,
        )

        # Pseudo-IRC & minimize both ends (no segment endpoints exist → fallback mapping in helper)
        irc_res = _pseudo_irc_and_match(seg_idx=1,
                                        seg_dir=tsroot,
                                        ref_pdb_for_seg=ts_pdb,
                                        seg_pocket_pdb=first_pocket,
                                        g_ts=g_ts,
                                        q_int=q_int,
                                        spin=spin,
                                        real_parm7=real_parm7_path,
                                        model_pdb=ml_region_pdb,
                                        detect_layer=detect_layer)
        gL = irc_res["left_min_geom"]
        gR = irc_res["right_min_geom"]
        gT = irc_res["ts_geom"]
        irc_plot_path = irc_res.get("irc_plot")
        irc_trj_path = irc_res.get("irc_trj")
        if irc_trj_path:
            try:
                irc_trj_for_all.append((Path(irc_trj_path), False))
            except Exception:
                pass

        # Ensure UMA energies
        eL = float(gL.energy)
        eT = float(gT.energy)
        eR = float(gR.energy)

        # In this mode ONLY: assign Reactant/Product so that higher-energy end is the Reactant
        if eL >= eR:
            g_react, e_react = gL, eL
            g_prod,  e_prod  = gR, eR
        else:
            g_react, e_react = gR, eR
            g_prod,  e_prod  = gL, eL

        # Save standardized PDBs
        struct_dir = tsroot / "structures"
        ensure_dir(struct_dir)
        pocket_ref = first_pocket
        pR = _save_single_geom_as_pdb_for_tools(g_react, pocket_ref, struct_dir, "reactant")
        pT = _save_single_geom_as_pdb_for_tools(gT,       pocket_ref, struct_dir, "ts")
        pP = _save_single_geom_as_pdb_for_tools(g_prod,   pocket_ref, struct_dir, "product")

        # UMA energy diagram (R, TS, P)
        uma_prefix = tsroot / "energy_diagram_UMA"
        uma_diag = _write_segment_energy_diagram(
            uma_prefix,
            labels=["R", "TS", "P"],
            energies_eh=[e_react, eT, e_prod],
            title_note="(UMA, TSOPT/IRC)",
        )
        g_uma_diag = None
        dft_diag = None
        g_dft_diag = None

        # Thermochemistry (UMA) Gibbs
        thermo_payloads: Dict[str, Dict[str, Any]] = {}
        GR = GT = GP = None
        eR_dft = eT_dft = eP_dft = None
        GR_dftUMA = GT_dftUMA = GP_dftUMA = None
        freq_root = _resolve_override_dir(tsroot / "freq", freq_out_dir)
        dft_root = _resolve_override_dir(tsroot / "dft", dft_out_dir)

        if do_thermo:
            _echo(f"[thermo] Single TSOPT: freq on R/TS/P")
            tR = _run_freq_for_state(pR, q_int, spin, real_parm7_path, ml_region_pdb, detect_layer,
                                     freq_root / "R", args_yaml, overrides=freq_overrides)
            tT = _run_freq_for_state(pT, q_int, spin, real_parm7_path, ml_region_pdb, detect_layer,
                                     freq_root / "TS", args_yaml, overrides=freq_overrides)
            tP = _run_freq_for_state(pP, q_int, spin, real_parm7_path, ml_region_pdb, detect_layer,
                                     freq_root / "P", args_yaml, overrides=freq_overrides)
            thermo_payloads = {"R": tR, "TS": tT, "P": tP}
            try:
                GR = float(tR.get("sum_EE_and_thermal_free_energy_ha", e_react))
                GT = float(tT.get("sum_EE_and_thermal_free_energy_ha", eT))
                GP = float(tP.get("sum_EE_and_thermal_free_energy_ha", e_prod))
                g_uma_diag = _write_segment_energy_diagram(
                    tsroot / "energy_diagram_G_UMA",
                    labels=["R", "TS", "P"],
                    energies_eh=[GR, GT, GP],
                    title_note="(Gibbs, UMA)",
                    ylabel="ΔG (kcal/mol)",
                )
            except Exception as e:
                _echo(f"[thermo] WARNING: failed to build Gibbs diagram: {e}", err=True)

        # DFT & DFT//UMA
        if do_dft:
            _echo(f"[dft] Single TSOPT: DFT on R/TS/P")
            dR = _run_dft_for_state(pR, q_int, spin, real_parm7_path, ml_region_pdb, detect_layer,
                                     dft_root / "R", args_yaml, func_basis=dft_func_basis_use, overrides=dft_overrides)
            dT = _run_dft_for_state(pT, q_int, spin, real_parm7_path, ml_region_pdb, detect_layer,
                                     dft_root / "TS", args_yaml, func_basis=dft_func_basis_use, overrides=dft_overrides)
            dP = _run_dft_for_state(pP, q_int, spin, real_parm7_path, ml_region_pdb, detect_layer,
                                     dft_root / "P", args_yaml, func_basis=dft_func_basis_use, overrides=dft_overrides)
            try:
                eR_dft = float(dR.get("energy", {}).get("hartree", e_react) if dR else e_react)
                eT_dft = float(dT.get("energy", {}).get("hartree", eT) if dT else eT)
                eP_dft = float(dP.get("energy", {}).get("hartree", e_prod) if dP else e_prod)
                dft_diag = _write_segment_energy_diagram(
                    tsroot / "energy_diagram_DFT",
                    labels=["R", "TS", "P"],
                    energies_eh=[eR_dft, eT_dft, eP_dft],
                    title_note=f"({dft_method_fallback})",
                )
            except Exception as e:
                _echo(f"[dft] WARNING: failed to build DFT diagram: {e}", err=True)

            if do_thermo:
                try:
                    dG_R = float(thermo_payloads.get("R", {}).get("thermal_correction_free_energy_ha", 0.0))
                    dG_T = float(thermo_payloads.get("TS", {}).get("thermal_correction_free_energy_ha", 0.0))
                    dG_P = float(thermo_payloads.get("P", {}).get("thermal_correction_free_energy_ha", 0.0))
                    GR_dftUMA = eR_dft + dG_R
                    GT_dftUMA = eT_dft + dG_T
                    GP_dftUMA = eP_dft + dG_P
                    g_dft_diag = _write_segment_energy_diagram(
                        tsroot / "energy_diagram_G_DFT_plus_UMA",
                        labels=["R", "TS", "P"],
                        energies_eh=[GR_dftUMA, GT_dftUMA, GP_dftUMA],
                        title_note="(Gibbs, DFT//UMA)",
                        ylabel="ΔG (kcal/mol)",
                    )
                except Exception as e:
                    _echo(f"[dft//uma] WARNING: failed to build DFT//UMA Gibbs diagram: {e}", err=True)

        # Summary.yaml / summary.log for TSOPT-only mode
        bond_cfg = dict(_path_search.BOND_KW)
        bond_summary = ""
        try:
            changed, bond_summary = _path_search._has_bond_change(g_react, g_prod, bond_cfg)
            if not changed:
                bond_summary = "(no covalent changes detected)"
        except Exception:
            bond_summary = "(no covalent changes detected)"

        barrier = (eT - e_react) * AU2KCALPERMOL
        delta = (e_prod - e_react) * AU2KCALPERMOL

        energy_diagrams: List[Dict[str, Any]] = []
        def _promote_all(diag: Optional[Dict[str, Any]], stem: str) -> None:
            if not diag:
                return
            diag_all = dict(diag)
            diag_all["name"] = f"{stem}_all"
            diag_all["image"] = str(out_dir / f"{stem}_all.png")
            energy_diagrams.append(diag_all)

        _promote_all(uma_diag, "energy_diagram_UMA")
        _promote_all(g_uma_diag, "energy_diagram_G_UMA")
        _promote_all(dft_diag, "energy_diagram_DFT")
        _promote_all(g_dft_diag, "energy_diagram_G_DFT_plus_UMA")

        summary = {
            "out_dir": str(tsroot),
            "n_images": 0,
            "n_segments": 1,
            "segments": [
                {
                    "index": 1,
                    "tag": "seg_01",
                    "kind": "tsopt",
                    "barrier_kcal": float(barrier),
                    "delta_kcal": float(delta),
                    "bond_changes": bond_summary,
                }
            ],
            "energy_diagrams": list(energy_diagrams),
        }
        try:
            with open(tsroot / "summary.yaml", "w") as f:
                yaml.safe_dump(summary, f, sort_keys=False, allow_unicode=True)
            shutil.copy2(tsroot / "summary.yaml", out_dir / "summary.yaml")
        except Exception as e:
            _echo(f"[write] WARNING: failed to write summary.yaml: {e}", err=True)

        segment_log: Dict[str, Any] = {
            "index": 1,
            "tag": "seg_01",
            "kind": "tsopt",
            "bond_changes": bond_summary,
            "post_dir": str(tsroot),
            "mep_barrier_kcal": barrier,
            "mep_delta_kcal": delta,
        }
        if irc_plot_path:
            segment_log["irc_plot"] = str(irc_plot_path)
        if irc_trj_path:
            segment_log["irc_traj"] = str(irc_trj_path)
        if do_thermo:
            n_imag = None
            try:
                n_imag = int(thermo_payloads.get("TS", {}).get("num_imag_freq"))
            except Exception:
                n_imag = None
            if n_imag is not None:
                segment_log["ts_imag"] = {"n_imag": n_imag}

        segment_log["uma"] = {
            "labels": ["R", "TS", "P"],
            "energies_au": [e_react, eT, e_prod],
            "energies_kcal": [0.0, barrier, delta],
            "diagram": str((tsroot / "energy_diagram_UMA").with_suffix(".png")),
            "structures": {"R": pR, "TS": pT, "P": pP},
            "barrier_kcal": barrier,
            "delta_kcal": delta,
        }
        if GR is not None and GT is not None and GP is not None:
            segment_log["gibbs_uma"] = {
                "labels": ["R", "TS", "P"],
                "energies_au": [GR, GT, GP],
                "energies_kcal": [
                    0.0,
                    (GT - GR) * AU2KCALPERMOL,
                    (GP - GR) * AU2KCALPERMOL,
                ],
                "diagram": str((tsroot / "energy_diagram_G_UMA").with_suffix(".png")),
                "structures": {"R": pR, "TS": pT, "P": pP},
                "barrier_kcal": (GT - GR) * AU2KCALPERMOL,
                "delta_kcal": (GP - GR) * AU2KCALPERMOL,
            }
        if eR_dft is not None and eT_dft is not None and eP_dft is not None:
            segment_log["dft"] = {
                "labels": ["R", "TS", "P"],
                "energies_au": [eR_dft, eT_dft, eP_dft],
                "energies_kcal": [
                    0.0,
                    (eT_dft - eR_dft) * AU2KCALPERMOL,
                    (eP_dft - eR_dft) * AU2KCALPERMOL,
                ],
                "diagram": str((tsroot / "energy_diagram_DFT").with_suffix(".png")),
                "structures": {"R": pR, "TS": pT, "P": pP},
                "barrier_kcal": (eT_dft - eR_dft) * AU2KCALPERMOL,
                "delta_kcal": (eP_dft - eR_dft) * AU2KCALPERMOL,
            }
        if GR_dftUMA is not None and GT_dftUMA is not None and GP_dftUMA is not None:
            segment_log["gibbs_dft_uma"] = {
                "labels": ["R", "TS", "P"],
                "energies_au": [GR_dftUMA, GT_dftUMA, GP_dftUMA],
                "energies_kcal": [
                    0.0,
                    (GT_dftUMA - GR_dftUMA) * AU2KCALPERMOL,
                    (GP_dftUMA - GR_dftUMA) * AU2KCALPERMOL,
                ],
                "diagram": str((tsroot / "energy_diagram_G_DFT_plus_UMA").with_suffix(".png")),
                "structures": {"R": pR, "TS": pT, "P": pP},
                "barrier_kcal": (GT_dftUMA - GR_dftUMA) * AU2KCALPERMOL,
                "delta_kcal": (GP_dftUMA - GR_dftUMA) * AU2KCALPERMOL,
            }

        summary_payload = {
            "root_out_dir": str(out_dir),
            "path_dir": str(tsroot),
            "path_module_dir": "tsopt_single",
            "pipeline_mode": "tsopt-only",
            "refine_path": True,
            "tsopt": do_tsopt,
            "thermo": do_thermo,
            "dft": do_dft,
            "opt_mode": tsopt_opt_mode_default,
            "mep_mode": "tsopt-only",
            "uma_model": None,
            "command": command_str,
            "charge": q_int,
            "spin": spin,
            "mep": {"n_images": 0, "n_segments": 1},
            "segments": summary.get("segments", []),
            "energy_diagrams": list(energy_diagrams),
            "post_segments": [segment_log],
            "key_files": {
                "summary.yaml": "YAML-format summary",
                "summary.log": "This summary",
                "energy_diagram_UMA_all.png": "UMA R-TS-P energies",
                "energy_diagram_G_UMA_all.png": "UMA Gibbs R-TS-P",
                "energy_diagram_DFT_all.png": "DFT R-TS-P",
                "energy_diagram_G_DFT_plus_UMA_all.png": "DFT//UMA Gibbs R-TS-P",
                "irc_plot_all.png": "Aggregated IRC plot",
            },
        }
        try:
            write_summary_log(tsroot / "summary.log", summary_payload)
            shutil.copy2(tsroot / "summary.log", out_dir / "summary.log")
        except Exception as e:
            _echo(f"[write] WARNING: failed to write summary.log: {e}", err=True)

        try:
            for stem in (
                "energy_diagram_UMA",
                "energy_diagram_G_UMA",
                "energy_diagram_DFT",
                "energy_diagram_G_DFT_plus_UMA",
            ):
                src = tsroot / f"{stem}.png"
                if src.exists():
                    shutil.copy2(src, out_dir / f"{stem}_all.png")
        except Exception as e:
            _echo(f"[all] WARNING: failed to copy *_all diagrams: {e}", err=True)

        try:
            if irc_plot_path:
                irc_plot_src = Path(irc_plot_path)
                if irc_plot_src.exists():
                    shutil.copy2(irc_plot_src, out_dir / "irc_plot_all.png")
        except Exception as e:
            _echo(f"[all] WARNING: failed to copy irc_plot_all.png: {e}", err=True)

        _write_output_summary_md(out_dir)
        _echo_section("=== [all] TSOPT-only pipeline finished successfully ===")
        _echo(format_elapsed("[all] Elapsed for Whole Pipeline", time_start))
        return

    # --------------------------
    # Stage 1c: Optional scan (single-structure only) to build ordered pocket inputs
    # --------------------------
    pockets_for_path: List[Path]
    if is_single and has_scan:
        _echo_section("=== [all] Stage 1c — Staged scan on layered full-system PDB (single-structure mode) ===")
        ensure_dir(scan_dir)
        layered_pdb = Path(layered_inputs[0]).resolve()
        full_input_pdb = Path(input_paths[0]).resolve()
        # Use the layered full-system PDB for scan (no pocket index remapping needed)
        full_atom_meta = load_pdb_atom_metadata(full_input_pdb)
        converted_scan_stages = _parse_scan_lists_literals(scan_lists_raw, atom_meta=full_atom_meta)
        scan_one_based_effective = True if scan_one_based is None else bool(scan_one_based)
        scan_stage_literals: List[str] = []
        for stage in converted_scan_stages:
            if scan_one_based_effective:
                stage_use = stage
            else:
                stage_use = [(i - 1, j - 1, target) for (i, j, target) in stage]
            scan_stage_literals.append(_format_scan_stage(stage_use))
        _echo("[all] Remapped --scan-lists indices from the full PDB to the pocket ordering.")
        scan_preopt_use = pre_opt if scan_preopt_override is None else bool(scan_preopt_override)
        scan_endopt_use = True if scan_endopt_override is None else bool(scan_endopt_override)
        scan_opt_mode_use = (opt_mode.lower() if opt_mode else
                              ("lbfgs" if sopt_mode.lower() in ("lbfgs", "light") else "rfo"))

        scan_args: List[str] = [
            "-i", str(layered_pdb),
            "--real-parm7", str(real_parm7_path),
            "-q", str(int(q_int)),
            "-m", str(int(spin)),
            "--out-dir", str(scan_dir),
            "--preopt", "True" if scan_preopt_use else "False",
            "--endopt", "True" if scan_endopt_use else "False",
            "--opt-mode", str(scan_opt_mode_use),
        ]
        scan_args.append("--detect-layer" if detect_layer else "--no-detect-layer")

        if dump_override_requested:
            scan_args.extend(["--dump", "True" if dump else "False"])

        if scan_one_based is not None:
            scan_args.append("--one-based" if scan_one_based else "--zero-based")

        _append_cli_arg(scan_args, "--max-step-size", scan_max_step_size)
        _append_cli_arg(scan_args, "--bias-k", scan_bias_k)
        _append_cli_arg(scan_args, "--relax-max-cycles", scan_relax_max_cycles)
        if args_yaml is not None:
            scan_args.extend(["--args-yaml", str(args_yaml)])
        # Forward all converted --scan-lists (aligned to the pocket atom order)
        for literal in scan_stage_literals:
            scan_args.extend(["--scan-lists", literal])

        _echo("[all] Invoking scan with arguments:")
        _echo("  " + " ".join(scan_args))

        _run_cli_main("scan", _scan_cli.cli, scan_args, on_nonzero="raise", on_exception="raise", prefix="all")

        # Collect stage results as pocket inputs
        stage_pdbs = sorted((scan_dir).glob("stage_*/*"))
        stage_results: List[Path] = []
        for p in stage_pdbs:
            if p.name == "result.pdb":
                stage_results.append(p.resolve())
        if not stage_results:
            raise click.ClickException("[all] No stage result PDBs found under scan/.")
        _echo("[all] Collected scan stage pocket files:")
        for p in stage_results:
            _echo(f"  - {p}")

        # Input series to path_search: [initial layered PDB, scan stage results ...]
        # Note: scan stage results are also full-system layered PDBs since scan ran on layered_pdb
        pockets_for_path = [layered_pdb] + stage_results
    else:
        # Multi-structure standard route: use layered full-system PDBs
        pockets_for_path = list(layered_inputs)

    # --------------------------
    # Stage 2: Path search on full-system layered PDBs
    # --------------------------
    _echo_section("=== [all] Stage 2/3 — MEP search on full-system layered PDBs (recursive GSM) ===")

    # Build path_search CLI args using *repeated* options (robust for Click)
    ps_args: List[str] = []

    # Inputs: single -i followed by all layered full-system PDBs
    ps_args.append("-i")
    for p in pockets_for_path:
        ps_args.append(str(p))

    # Charge & spin
    ps_args.extend(["-q", str(q_int)])
    ps_args.extend(["-m", str(int(spin))])
    ps_args.extend(["--real-parm7", str(real_parm7_path)])
    # Layered PDBs have B-factors → detect-layer will auto-identify layers
    ps_args.append("--detect-layer")

    # Nodes, cycles, climb, optimizer, dump, out-dir, pre-opt, args-yaml
    ps_args.extend(["--max-nodes", str(int(max_nodes))])
    ps_args.extend(["--max-cycles", str(int(max_cycles))])
    ps_args.extend(["--climb", "True" if climb else "False"])
    ps_sopt = "light" if sopt_mode.lower() in ("lbfgs", "light") else "heavy"
    ps_args.extend(["--sopt-mode", ps_sopt])
    ps_args.extend(["--dump", "True" if dump else "False"])
    ps_args.extend(["--out-dir", str(path_dir)])
    ps_args.extend(["--pre-opt", "True" if pre_opt else "False"])
    if args_yaml is not None:
        ps_args.extend(["--args-yaml", str(args_yaml)])

    # Auto-provide ref templates (original full PDBs) for full-system merge.
    # Multi-structure: one ref per original input; single+scan: reuse the single template for all pockets.
    ps_args.append("--ref-pdb")
    if is_single and has_scan:
        for _ in pockets_for_path:
            ps_args.append(str(input_paths[0]))
    else:
        for p in input_paths:
            ps_args.append(str(p))

    _echo("[all] Invoking path_search with arguments:")
    _echo("  " + " ".join(ps_args))

    _run_cli_main("path_search", _path_search.cli, ps_args, on_nonzero="raise", on_exception="raise", prefix="all")

    # --------------------------
    # Stage 3: Merge (performed by path_search when --ref-pdb was supplied)
    # --------------------------
    _echo_section("=== [all] Stage 3/3 — Merge into full-system templates ===")
    _echo("[all] Merging was carried out by path_search using the original inputs as templates.")
    _echo(f"[all] Final products can be found under: {path_dir}")
    _echo("  - mep_w_ref.pdb            (full-system merged trajectory)")
    _echo("  - mep_w_ref_seg_XX.pdb     (per-segment merged trajectories for covalent-change segments)")
    _echo("  - summary.yaml             (segment barriers, ΔE, labels)")
    _echo("  - mep_plot.png / energy_diagram_MEP.png / summary.log")
    _echo_section("=== [all] Pipeline finished successfully (core path) ===")

    summary_yaml = path_dir / "summary.yaml"
    summary_loaded = {}
    if summary_yaml.exists():
        try:
            summary_loaded = yaml.safe_load(summary_yaml.read_text(encoding="utf-8")) or {}
        except Exception:
            summary_loaded = {}
    summary: Dict[str, Any] = summary_loaded if isinstance(summary_loaded, dict) else {}
    segments = _read_summary(summary_yaml)
    energy_diagrams: List[Dict[str, Any]] = []
    existing_diagrams = summary.get("energy_diagrams", [])
    if isinstance(existing_diagrams, list):
        energy_diagrams.extend(existing_diagrams)

    def _copy_path_outputs_to_root() -> None:
        try:
            for name in (
                "mep_plot.png",
                "energy_diagram_MEP.png",
                "mep.pdb",
                "mep_w_ref.pdb",
                "summary.yaml",
                "summary.log",
            ):
                src = path_dir / name
                if src.exists():
                    shutil.copy2(src, out_dir / name)
            for stem in ("mep", "mep_w_ref"):
                for ext in (".trj", ".xyz"):
                    src = path_dir / f"{stem}{ext}"
                    if src.exists():
                        shutil.copy2(src, out_dir / src.name)
        except Exception as e:
            _echo(f"[all] WARNING: Failed to copy path_search outputs: {e}", err=True)

    def _write_pipeline_summary_log(post_segment_logs: Sequence[Dict[str, Any]]) -> None:
        try:
            diag_for_log: Dict[str, Any] = {}
            for diag in summary.get("energy_diagrams", []) or []:
                if isinstance(diag, dict) and str(diag.get("name", "")).lower().endswith("mep"):
                    diag_for_log = diag
                    break
            mep_info = {
                "n_images": summary.get("n_images"),
                "n_segments": summary.get("n_segments"),
                "traj_pdb": str(path_dir / "mep.pdb") if (path_dir / "mep.pdb").exists() else None,
                "mep_plot": str(path_dir / "mep_plot.png") if (path_dir / "mep_plot.png").exists() else None,
                "diagram": diag_for_log,
            }
            summary_payload = {
                "root_out_dir": str(out_dir),
                "path_dir": str(path_dir),
                "path_module_dir": path_dir.name,
                "pipeline_mode": "path-search",
                "refine_path": True,
                "tsopt": do_tsopt,
                "thermo": do_thermo,
                "dft": do_dft,
                "opt_mode": opt_mode.lower() if opt_mode else None,
                "mep_mode": "path-search",
                "uma_model": None,
                "command": command_str,
                "charge": q_int,
                "spin": spin,
                "mep": mep_info,
                "segments": summary.get("segments", []),
                "energy_diagrams": summary.get("energy_diagrams", []),
                "post_segments": list(post_segment_logs),
                "key_files": {
                    "summary.yaml": "YAML-format summary",
                    "summary.log": "This summary",
                    "mep_plot.png": "ML/MM MEP energy plot",
                    "energy_diagram_MEP.png": "Compressed MEP diagram",
                    "energy_diagram_UMA_all.png": "UMA R-TS-P energies (all segments)",
                    "energy_diagram_G_UMA_all.png": "UMA Gibbs R-TS-P (all segments)",
                    "energy_diagram_DFT_all.png": "DFT R-TS-P (all segments)",
                    "energy_diagram_G_DFT_plus_UMA_all.png": "DFT//UMA Gibbs R-TS-P (all segments)",
                    "irc_plot_all.png": "Aggregated IRC plot",
                },
            }
            write_summary_log(path_dir / "summary.log", summary_payload)
            _copy_path_outputs_to_root()
        except Exception as e:
            _echo(f"[write] WARNING: Failed to write summary.log: {e}", err=True)

    # --------------------------
    # Optional Stage 4: TSOPT / THERMO / DFT (per reactive segment)
    # --------------------------
    if not (do_tsopt or do_thermo or do_dft):
        _write_pipeline_summary_log([])
        _write_output_summary_md(out_dir)
        # Elapsed time
        _echo(format_elapsed("[all] Elapsed for Whole Pipeline", time_start))
        return

    _echo_section("=== [all] Stage 4 — Post-processing per reactive segment ===")

    # Use segment summary from path_search
    if not segments:
        _echo("[post] No segments found in summary; nothing to do.")
        _write_pipeline_summary_log([])
        _write_output_summary_md(out_dir)
        _echo(format_elapsed("[all] Elapsed for Whole Pipeline", time_start))
        return

    # Iterate only bond-change segments (kind='seg' and bond_changes not empty and not '(no covalent...)')
    reactive = [s for s in segments if (s.get("kind", "seg") == "seg" and str(s.get("bond_changes", "")).strip() and str(s.get("bond_changes", "")).strip() != "(no covalent changes detected)")]
    if not reactive:
        _echo("[post] No bond-change segments. Skipping TS/thermo/DFT.")
        _write_pipeline_summary_log([])
        _write_output_summary_md(out_dir)
        _echo(format_elapsed("[all] Elapsed for Whole Pipeline", time_start))
        return

    post_segment_logs: List[Dict[str, Any]] = []
    tsopt_seg_energies: List[Tuple[float, float, float]] = []
    g_uma_seg_energies: List[Tuple[float, float, float]] = []
    dft_seg_energies: List[Tuple[float, float, float]] = []
    g_dftuma_seg_energies: List[Tuple[float, float, float]] = []
    irc_trj_for_all: List[Tuple[Path, bool]] = []

    # For each reactive segment
    for s in reactive:
        seg_idx = int(s.get("index", 0) or 0)
        seg_tag = s.get("tag", f"seg_{seg_idx:02d}")
        _echo_section(f"--- [post] Segment {seg_idx:02d} ({seg_tag}) ---")

        seg_root = path_dir  # base
        seg_dir = seg_root / f"post_seg_{seg_idx:02d}"
        ensure_dir(seg_dir)

        # HEI pocket file prepared by path_search (only for bond-change segments)
        hei_pocket_pdb = seg_root / f"hei_seg_{seg_idx:02d}.pdb"
        if not hei_pocket_pdb.exists():
            _echo(f"[post] WARNING: HEI pocket PDB not found for segment {seg_idx:02d}; skipping TSOPT.", err=True)
            continue

        # 4.1 TS optimization (optional; still needed to drive IRC & diagrams)
        if do_tsopt:
            ts_pdb, g_ts = _run_tsopt_on_hei(
                hei_pocket_pdb,
                q_int,
                spin,
                real_parm7_path,
                ml_region_pdb,
                detect_layer,
                args_yaml,
                seg_dir,
                tsopt_opt_mode_default,
                overrides=tsopt_overrides,
            )
        else:
            # If TSOPT off: use the GSM HEI (pocket) as TS geometry
            ts_pdb = hei_pocket_pdb
            g_ts = geom_loader(ts_pdb, coord_type="cart")
            calc = _mlmm_calc(
                model_charge=int(q_int),
                model_mult=int(spin),
                input_pdb=str(ts_pdb),
                real_parm7=str(real_parm7_path),
                model_pdb=str(ml_region_pdb),
                use_bfactor_layers=detect_layer,
            )
            g_ts.set_calculator(calc); _ = float(g_ts.energy)

        # 4.2 Pseudo-IRC & mapping to (left,right)
        irc_plot_path = None
        irc_trj_path = None
        irc_res = _pseudo_irc_and_match(seg_idx=seg_idx,
                                        seg_dir=seg_dir,
                                        ref_pdb_for_seg=ts_pdb,
                                        seg_pocket_pdb=hei_pocket_pdb,
                                        g_ts=g_ts,
                                        q_int=q_int,
                                        spin=spin,
                                        real_parm7=real_parm7_path,
                                        model_pdb=ml_region_pdb,
                                        detect_layer=detect_layer)
        irc_plot_path = irc_res.get("irc_plot")
        irc_trj_path = irc_res.get("irc_trj")
        if irc_trj_path:
            try:
                irc_trj_for_all.append((Path(irc_trj_path), False))
            except Exception:
                pass

        gL = irc_res["left_min_geom"]
        gR = irc_res["right_min_geom"]
        gT = irc_res["ts_geom"]
        # Save standardized PDBs for tools
        struct_dir = seg_dir / "structures"
        ensure_dir(struct_dir)
        pL = _save_single_geom_as_pdb_for_tools(gL, hei_pocket_pdb, struct_dir, "reactant")
        pT = _save_single_geom_as_pdb_for_tools(gT, hei_pocket_pdb, struct_dir, "ts")
        pR = _save_single_geom_as_pdb_for_tools(gR, hei_pocket_pdb, struct_dir, "product")

        # 4.3 Segment-level energy diagram from UMA (R,TS,P)
        eR = float(gL.energy)
        eT = float(gT.energy)
        eP = float(gR.energy)
        tsopt_seg_energies.append((eR, eT, eP))
        uma_prefix = seg_dir / "energy_diagram_UMA"
        _write_segment_energy_diagram(
            uma_prefix,
            labels=["R", f"TS{seg_idx}", "P"],
            energies_eh=[eR, eT, eP],
            title_note="(UMA, TSOPT/IRC)",
        )

        # 4.4 Thermochemistry (UMA freq) and Gibbs diagram
        thermo_payloads: Dict[str, Dict[str, Any]] = {}
        GR = GT = GP = None
        freq_seg_root = _resolve_override_dir(seg_dir / "freq", freq_out_dir)
        dft_seg_root = _resolve_override_dir(seg_dir / "dft", dft_out_dir)

        if do_thermo:
            _echo(f"[thermo] Segment {seg_idx:02d}: freq on R/TS/P")
            tR = _run_freq_for_state(pL, q_int, spin, real_parm7_path, ml_region_pdb, detect_layer,
                                     freq_seg_root / "R", args_yaml, overrides=freq_overrides)
            tT = _run_freq_for_state(pT, q_int, spin, real_parm7_path, ml_region_pdb, detect_layer,
                                     freq_seg_root / "TS", args_yaml, overrides=freq_overrides)
            tP = _run_freq_for_state(pR, q_int, spin, real_parm7_path, ml_region_pdb, detect_layer,
                                     freq_seg_root / "P", args_yaml, overrides=freq_overrides)
            thermo_payloads = {"R": tR, "TS": tT, "P": tP}
            try:
                GR = float(tR.get("sum_EE_and_thermal_free_energy_ha", eR))
                GT = float(tT.get("sum_EE_and_thermal_free_energy_ha", eT))
                GP = float(tP.get("sum_EE_and_thermal_free_energy_ha", eP))
                g_uma_seg_energies.append((GR, GT, GP))
                _write_segment_energy_diagram(
                    seg_dir / "energy_diagram_G_UMA",
                    labels=["R", f"TS{seg_idx}", "P"],
                    energies_eh=[GR, GT, GP],
                    title_note="(Gibbs, UMA)",
                    ylabel="ΔG (kcal/mol)",
                )
            except Exception as e:
                _echo(f"[thermo] WARNING: failed to build Gibbs diagram: {e}", err=True)

        # 4.5 DFT single-point and (optionally) DFT//UMA Gibbs
        eR_dft = eT_dft = eP_dft = None
        GR_dftUMA = GT_dftUMA = GP_dftUMA = None
        if do_dft:
            _echo(f"[dft] Segment {seg_idx:02d}: DFT on R/TS/P")
            dR = _run_dft_for_state(pL, q_int, spin, real_parm7_path, ml_region_pdb, detect_layer,
                                     dft_seg_root / "R", args_yaml, func_basis=dft_func_basis_use, overrides=dft_overrides)
            dT = _run_dft_for_state(pT, q_int, spin, real_parm7_path, ml_region_pdb, detect_layer,
                                     dft_seg_root / "TS", args_yaml, func_basis=dft_func_basis_use, overrides=dft_overrides)
            dP = _run_dft_for_state(pR, q_int, spin, real_parm7_path, ml_region_pdb, detect_layer,
                                     dft_seg_root / "P", args_yaml, func_basis=dft_func_basis_use, overrides=dft_overrides)
            try:
                eR_dft = float(dR.get("energy", {}).get("hartree", np.nan) if dR else np.nan)
                eT_dft = float(dT.get("energy", {}).get("hartree", np.nan) if dT else np.nan)
                eP_dft = float(dP.get("energy", {}).get("hartree", np.nan) if dP else np.nan)
                if all(map(np.isfinite, [eR_dft, eT_dft, eP_dft])):
                    dft_seg_energies.append((eR_dft, eT_dft, eP_dft))
                    _write_segment_energy_diagram(seg_dir / "energy_diagram_DFT",
                                                  labels=["R", f"TS{seg_idx}", "P"],
                                                  energies_eh=[eR_dft, eT_dft, eP_dft],
                                                  title_note=f"({dft_method_fallback})")
                else:
                    _echo("[dft] WARNING: some DFT energies missing; diagram skipped.", err=True)
            except Exception as e:
                _echo(f"[dft] WARNING: failed to build DFT diagram: {e}", err=True)

            # DFT//UMA thermal Gibbs (E_DFT + ΔG_therm(UMA))
            if do_thermo:
                try:
                    dG_R = float(thermo_payloads.get("R", {}).get("thermal_correction_free_energy_ha", 0.0))
                    dG_T = float(thermo_payloads.get("TS", {}).get("thermal_correction_free_energy_ha", 0.0))
                    dG_P = float(thermo_payloads.get("P", {}).get("thermal_correction_free_energy_ha", 0.0))
                    eR_dft = float(dR.get("energy", {}).get("hartree", eR) if dR else eR)
                    eT_dft = float(dT.get("energy", {}).get("hartree", eT) if dT else eT)
                    eP_dft = float(dP.get("energy", {}).get("hartree", eP) if dP else eP)
                    GR_dftUMA = eR_dft + dG_R
                    GT_dftUMA = eT_dft + dG_T
                    GP_dftUMA = eP_dft + dG_P
                    g_dftuma_seg_energies.append((GR_dftUMA, GT_dftUMA, GP_dftUMA))
                    _write_segment_energy_diagram(
                        seg_dir / "energy_diagram_G_DFT_plus_UMA",
                        labels=["R", f"TS{seg_idx}", "P"],
                        energies_eh=[GR_dftUMA, GT_dftUMA, GP_dftUMA],
                        title_note="(Gibbs, DFT//UMA)",
                        ylabel="ΔG (kcal/mol)",
                    )
                except Exception as e:
                    _echo(f"[dft//uma] WARNING: failed to build DFT//UMA Gibbs diagram: {e}", err=True)

        segment_log: Dict[str, Any] = {
            "index": seg_idx,
            "tag": seg_tag,
            "kind": s.get("kind", "seg"),
            "bond_changes": s.get("bond_changes", ""),
            "mep_barrier_kcal": s.get("barrier_kcal"),
            "mep_delta_kcal": s.get("delta_kcal"),
            "post_dir": str(seg_dir),
        }
        if irc_plot_path:
            segment_log["irc_plot"] = str(irc_plot_path)
        if irc_trj_path:
            segment_log["irc_traj"] = str(irc_trj_path)
        if do_thermo:
            n_imag = None
            try:
                n_imag = int(thermo_payloads.get("TS", {}).get("num_imag_freq"))
            except Exception:
                n_imag = None
            if n_imag is not None:
                segment_log["ts_imag"] = {"n_imag": n_imag}
        segment_log["uma"] = {
            "labels": ["R", "TS", "P"],
            "energies_au": [eR, eT, eP],
            "energies_kcal": [
                0.0,
                (eT - eR) * AU2KCALPERMOL,
                (eP - eR) * AU2KCALPERMOL,
            ],
            "diagram": str((seg_dir / "energy_diagram_UMA").with_suffix(".png")),
            "structures": {"R": pL, "TS": pT, "P": pR},
            "barrier_kcal": (eT - eR) * AU2KCALPERMOL,
            "delta_kcal": (eP - eR) * AU2KCALPERMOL,
        }
        if GR is not None and GT is not None and GP is not None:
            segment_log["gibbs_uma"] = {
                "labels": ["R", "TS", "P"],
                "energies_au": [GR, GT, GP],
                "energies_kcal": [
                    0.0,
                    (GT - GR) * AU2KCALPERMOL,
                    (GP - GR) * AU2KCALPERMOL,
                ],
                "diagram": str((seg_dir / "energy_diagram_G_UMA").with_suffix(".png")),
                "structures": {"R": pL, "TS": pT, "P": pR},
                "barrier_kcal": (GT - GR) * AU2KCALPERMOL,
                "delta_kcal": (GP - GR) * AU2KCALPERMOL,
            }
        if eR_dft is not None and eT_dft is not None and eP_dft is not None and all(
            map(np.isfinite, [eR_dft, eT_dft, eP_dft])
        ):
            segment_log["dft"] = {
                "labels": ["R", "TS", "P"],
                "energies_au": [eR_dft, eT_dft, eP_dft],
                "energies_kcal": [
                    0.0,
                    (eT_dft - eR_dft) * AU2KCALPERMOL,
                    (eP_dft - eR_dft) * AU2KCALPERMOL,
                ],
                "diagram": str((seg_dir / "energy_diagram_DFT").with_suffix(".png")),
                "structures": {"R": pL, "TS": pT, "P": pR},
                "barrier_kcal": (eT_dft - eR_dft) * AU2KCALPERMOL,
                "delta_kcal": (eP_dft - eR_dft) * AU2KCALPERMOL,
            }
        if GR_dftUMA is not None and GT_dftUMA is not None and GP_dftUMA is not None:
            segment_log["gibbs_dft_uma"] = {
                "labels": ["R", "TS", "P"],
                "energies_au": [GR_dftUMA, GT_dftUMA, GP_dftUMA],
                "energies_kcal": [
                    0.0,
                    (GT_dftUMA - GR_dftUMA) * AU2KCALPERMOL,
                    (GP_dftUMA - GR_dftUMA) * AU2KCALPERMOL,
                ],
                "diagram": str((seg_dir / "energy_diagram_G_DFT_plus_UMA").with_suffix(".png")),
                "structures": {"R": pL, "TS": pT, "P": pR},
                "barrier_kcal": (GT_dftUMA - GR_dftUMA) * AU2KCALPERMOL,
                "delta_kcal": (GP_dftUMA - GR_dftUMA) * AU2KCALPERMOL,
            }

        post_segment_logs.append(segment_log)

    # -------------------------------------------------------------------------
    # Aggregate diagrams over all reactive segments
    # -------------------------------------------------------------------------
    if tsopt_seg_energies:
        tsopt_all_energies = [e for triple in tsopt_seg_energies for e in triple]
        tsopt_all_labels = _build_global_segment_labels(len(tsopt_seg_energies))
        if tsopt_all_labels and len(tsopt_all_labels) == len(tsopt_all_energies):
            diag_payload = _write_segment_energy_diagram(
                out_dir / "energy_diagram_UMA_all",
                labels=tsopt_all_labels,
                energies_eh=tsopt_all_energies,
                title_note="(UMA, TSOPT + IRC; all segments)",
                write_html=False,
            )
            if diag_payload:
                energy_diagrams.append(diag_payload)

    if do_thermo and g_uma_seg_energies:
        g_uma_all_energies = [e for triple in g_uma_seg_energies for e in triple]
        g_uma_all_labels = _build_global_segment_labels(len(g_uma_seg_energies))
        if g_uma_all_labels and len(g_uma_all_labels) == len(g_uma_all_energies):
            diag_payload = _write_segment_energy_diagram(
                out_dir / "energy_diagram_G_UMA_all",
                labels=g_uma_all_labels,
                energies_eh=g_uma_all_energies,
                title_note="(UMA + Thermal Correction; all segments)",
                ylabel="ΔG (kcal/mol)",
                write_html=False,
            )
            if diag_payload:
                energy_diagrams.append(diag_payload)

    if do_dft and dft_seg_energies:
        dft_all_energies = [e for triple in dft_seg_energies for e in triple]
        dft_all_labels = _build_global_segment_labels(len(dft_seg_energies))
        if dft_all_labels and len(dft_all_labels) == len(dft_all_energies):
            diag_payload = _write_segment_energy_diagram(
                out_dir / "energy_diagram_DFT_all",
                labels=dft_all_labels,
                energies_eh=dft_all_energies,
                title_note=f"({dft_method_fallback}; all segments)",
                write_html=False,
            )
            if diag_payload:
                energy_diagrams.append(diag_payload)

    if do_dft and do_thermo and g_dftuma_seg_energies:
        g_dftuma_all_energies = [e for triple in g_dftuma_seg_energies for e in triple]
        g_dftuma_all_labels = _build_global_segment_labels(len(g_dftuma_seg_energies))
        if g_dftuma_all_labels and len(g_dftuma_all_labels) == len(g_dftuma_all_energies):
            diag_payload = _write_segment_energy_diagram(
                out_dir / "energy_diagram_G_DFT_plus_UMA_all",
                labels=g_dftuma_all_labels,
                energies_eh=g_dftuma_all_energies,
                title_note=f"({dft_method_fallback} // UMA + Thermal Correction; all segments)",
                ylabel="ΔG (kcal/mol)",
                write_html=False,
            )
            if diag_payload:
                energy_diagrams.append(diag_payload)

    # -------------------------------------------------------------------------
    # Aggregated IRC plot over all reactive segments
    # -------------------------------------------------------------------------
    if irc_trj_for_all:
        _merge_irc_trajectories_to_single_plot(
            irc_trj_for_all, out_dir / "irc_plot_all.png"
        )

    # Refresh summary.yaml with final energy diagram metadata
    try:
        summary["energy_diagrams"] = list(energy_diagrams)
        with open(path_dir / "summary.yaml", "w") as f:
            yaml.safe_dump(summary, f, sort_keys=False, allow_unicode=True)
        try:
            shutil.copy2(path_dir / "summary.yaml", out_dir / "summary.yaml")
        except Exception as e:
            _echo(f"[all] WARNING: Failed to mirror summary.yaml to {out_dir}: {e}", err=True)
    except Exception as e:
        _echo(f"[write] WARNING: Failed to refresh summary.yaml with energy diagram metadata: {e}", err=True)

    _write_pipeline_summary_log(post_segment_logs)
    _write_output_summary_md(out_dir)
    _echo(format_elapsed("[all] Elapsed for Whole Pipeline", time_start))


_configure_all_help_visibility(cli)


if __name__ == "__main__":
    cli()
