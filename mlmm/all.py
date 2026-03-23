# mlmm/all.py

"""
End-to-end enzymatic reaction workflow: extract, scan, MEP, TS, IRC, freq, DFT.

Example:
    mlmm all -i R.pdb P.pdb -c 'GPP,MMT' -l 'GPP:-3,MMT:-1'

For detailed documentation, see: docs/all.md
"""

from __future__ import annotations

from collections import defaultdict
from pathlib import Path
from typing import List, Sequence, Optional, Tuple, Dict, Any
import shutil
import tempfile

import gc
import logging
import sys
import math
import click
from .cli_utils import make_is_param_explicit
import time  # timing
import yaml
import numpy as np
import torch

logger = logging.getLogger(__name__)

# Biopython for PDB parsing (post-processing helpers)
from Bio import PDB

# pysisyphus helpers/constants
from pysisyphus.helpers import geom_loader
from pysisyphus.constants import BOHR2ANG, AU2KCALPERMOL

# Local imports from the package
from .extract import extract_api, compute_charge_summary, log_charge_summary
from . import path_search as _path_search
from . import path_opt as _path_opt
from . import opt as _opt_cli
from . import tsopt as _ts_opt
from . import freq as _freq_cli
from . import dft as _dft_cli
from . import irc as _irc_cli

from .trj2fig import run_trj2fig
from .summary_log import write_summary_log
from .utils import (
    apply_ref_pdb_override,
    build_energy_diagram,
    close_matplotlib_figures,
    deep_update,
    ensure_dir,
    format_elapsed,
    parse_xyz_block,
    prepare_input_structure,
    collect_single_option_values,
    load_yaml_dict,
    load_pdb_atom_metadata,
    parse_scan_list_triples,
    read_xyz_as_blocks,
    read_xyz_first_last,
    xyz_blocks_first_last,
)
from .cli_utils import resolve_yaml_sources, load_merged_yaml_cfg
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

class _EchoState:
    """Encapsulate CLI output state for section-spacing logic."""

    def __init__(self) -> None:
        self._started = False

    def reset(self) -> None:
        self._started = False

    def echo(self, *args, **kwargs) -> None:
        click.echo(*args, **kwargs)
        self._started = True

    def section(self, message: str, **kwargs) -> None:
        if self._started:
            click.echo()
        click.echo(message, **kwargs)
        self._started = True


_echo_state = _EchoState()


def _echo(*args, **kwargs) -> None:
    """Echo with local output tracking for section spacing."""
    _echo_state.echo(*args, **kwargs)


def _echo_section(message: str, **kwargs) -> None:
    """Echo a section header with a leading blank line unless it's the first log."""
    _echo_state.section(message, **kwargs)


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
        _echo("")
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
        # Release GPU memory between pipeline stages to prevent OOM.
        # Subcommand finally blocks unbind their heavy locals (= None).
        # gc.collect() is needed to break cyclic refs inside torch.nn.Module,
        # then empty_cache() reclaims the CUDA allocator cache.
        gc.collect()
        if torch.cuda.is_available():
            torch.cuda.empty_cache()
        _echo("")


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


def _append_toggle_arg(args: List[str], flag: str, value: Any | None) -> None:
    """Append Click bool-toggle option as ``--flag`` / ``--no-flag`` when value is not ``None``."""
    if value is None:
        return
    if not isinstance(value, bool):
        raise TypeError(f"Toggle flag '{flag}' requires bool value, got {type(value).__name__}.")
    base = flag if not flag.startswith("--no-") else f"--{flag[5:]}"
    neg = f"--no-{base[2:]}"
    args.append(base if value else neg)


def _resolve_override_dir(default: Path, override: Path | None) -> Path:
    """Return ``override`` when provided (respecting absolute paths); otherwise ``default``."""
    if override is None:
        return default
    if override.is_absolute():
        return override
    return default.parent / override


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
    merged, base_cfg, override_cfg = load_merged_yaml_cfg(config_yaml, override_yaml)

    if config_yaml is None and override_yaml is None:
        return None, {}
    if config_yaml is None:
        return override_yaml, override_cfg
    if override_yaml is None:
        return config_yaml, base_cfg

    with tempfile.NamedTemporaryFile(
        mode="w",
        encoding="utf-8",
        suffix=".yaml",
        prefix=tmp_prefix,
        delete=False,
    ) as tf:
        yaml.safe_dump(merged, tf, sort_keys=False, allow_unicode=True)
        effective = Path(tf.name).resolve()

    # Register cleanup so the temp file is removed when the process exits.
    import atexit
    atexit.register(lambda p=effective: p.unlink(missing_ok=True))

    return effective, merged


def _write_ml_region_definition(pocket_pdb: Path, dest: Path) -> Path:
    """
    Copy ``pocket_pdb`` to ``dest`` for downstream ML/MM commands.

    The copy preserves whatever link-hydrogen policy was used during extraction; set ``--add-linkh False``
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
    add_ter: bool,
    keep_temp: bool,
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
        keep_temp=bool(keep_temp),
        add_ter=bool(add_ter),
        add_h=False,
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


def _derive_charge_from_ligand_charge_when_extract_skipped(
    pdb_path: Path,
    ligand_charge: Optional[str],
) -> Optional[int]:
    """Derive total charge from a PDB using extract-style charge summary.

    *pdb_path* may be a full-complex PDB or a --model-pdb pocket.
    """
    if ligand_charge is None:
        return None
    try:
        parser = PDB.PDBParser(QUIET=True)
        complex_struct = parser.get_structure("complex", str(pdb_path))
        selected_ids = {res.get_full_id() for res in complex_struct.get_residues()}
        summary = compute_charge_summary(complex_struct, selected_ids, set(), ligand_charge)
        log_charge_summary("[all]", summary)
        q_total = float(summary.get("total_charge", 0.0))
        click.echo(f"[all] Charge summary from {pdb_path.name} (--ligand-charge without extraction):")
        click.echo(
            f"  Protein: {summary.get('protein_charge', 0.0):+g},  "
            f"Ligand: {summary.get('ligand_total_charge', 0.0):+g},  "
            f"Ions: {summary.get('ion_total_charge', 0.0):+g},  "
            f"Total: {q_total:+g}"
        )
        return _round_charge_with_note(q_total)
    except Exception as e:
        click.echo(
            f"[all] NOTE: failed to derive total charge from --ligand-charge: {e}",
            err=True,
        )
        return None


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


def _irc_and_match(seg_idx: int,
                   seg_dir: Path,
                   ref_pdb_for_seg: Path,
                   seg_pocket_pdb: Path,
                   g_ts: Any,
                   q_int: int,
                   spin: int,
                   real_parm7: Optional[Path] = None,
                   model_pdb: Optional[Path] = None,
                   detect_layer: bool = False,
                   backend: Optional[str] = None,
                   embedcharge: bool = False,
                   embedcharge_cutoff: Optional[float] = None,
                   link_atom_method: Optional[str] = None,
                   mm_backend: Optional[str] = None,
                   use_cmap: Optional[bool] = None,
                   args_yaml: Optional[Path] = None) -> Dict[str, Any]:
    """
    Run EulerPC IRC from a TS geometry, then map the IRC endpoints to (left, right)
    by comparing bond states with the GSM segment endpoints (when available).
    Falls back to raw IRC orientation in TSOPT-only mode.
    """
    irc_dir = seg_dir / "irc"
    ensure_dir(irc_dir)

    # Build irc CLI arguments
    irc_args: List[str] = [
        "-i", str(ref_pdb_for_seg),
        "--parm", str(real_parm7),
        "--model-pdb", str(model_pdb),
        "-q", str(int(q_int)),
        "-m", str(int(spin)),
        "--out-dir", str(irc_dir),
    ]
    irc_args.append("--detect-layer" if detect_layer else "--no-detect-layer")
    if backend is not None:
        irc_args.extend(["--backend", str(backend)])
    if embedcharge:
        irc_args.append("--embedcharge")
        if embedcharge_cutoff is not None:
            irc_args.extend(["--embedcharge-cutoff", str(embedcharge_cutoff)])
    else:
        irc_args.append("--no-embedcharge")
    if link_atom_method is not None:
        irc_args.extend(["--link-atom-method", str(link_atom_method)])
    if mm_backend is not None:
        irc_args.extend(["--mm-backend", str(mm_backend)])
    if use_cmap is not None:
        irc_args.extend(["--cmap" if use_cmap else "--no-cmap"])
    if args_yaml is not None:
        irc_args.extend(["--config", str(args_yaml)])

    _echo(f"[irc] Running EulerPC IRC → out={irc_dir}")
    _run_cli_main("irc", _irc_cli.cli, irc_args, on_nonzero="raise", prefix="irc")

    # Read IRC endpoints
    finished_trj = irc_dir / "finished_irc_trj.xyz"
    finished_pdb = irc_dir / "finished_irc.pdb"
    irc_plot = irc_dir / "irc_plot.png"

    if not finished_trj.exists():
        raise click.ClickException(f"[irc] IRC trajectory not found: {finished_trj}")

    # Convert to PDB if not already done
    if not finished_pdb.exists():
        _path_search._maybe_convert_to_pdb(finished_trj, ref_pdb_path=seg_pocket_pdb, out_path=finished_pdb)

    elems, c_first, c_last = read_xyz_first_last(finished_trj)

    # Create geometries from IRC endpoints
    _irc_calc_kwargs = dict(
        model_charge=int(q_int),
        model_mult=int(spin),
        input_pdb=str(ref_pdb_for_seg),
        real_parm7=str(real_parm7) if real_parm7 else None,
        model_pdb=str(model_pdb) if model_pdb else None,
        use_bfactor_layers=detect_layer,
        backend=backend,
        embedcharge=embedcharge,
    )
    if link_atom_method is not None:
        _irc_calc_kwargs["link_atom_method"] = link_atom_method
    if mm_backend is not None:
        _irc_calc_kwargs["mm_backend"] = mm_backend
    if use_cmap is not None:
        _irc_calc_kwargs["use_cmap"] = use_cmap
    calc = _mlmm_calc(**_irc_calc_kwargs)

    g_left = _path_search._new_geom_from_coords(
        elems, c_first / BOHR2ANG, coord_type="cart", freeze_atoms=[])
    g_right = _path_search._new_geom_from_coords(
        elems, c_last / BOHR2ANG, coord_type="cart", freeze_atoms=[])
    g_left.set_calculator(calc)
    g_right.set_calculator(calc)
    _ = float(g_left.energy)
    _ = float(g_right.energy)

    # Reload TS geometry with energy
    if g_ts.calculator is None:
        g_ts.set_calculator(calc)
    _ = float(g_ts.energy)

    left_tag = "backward"
    right_tag = "forward"
    reverse_irc = False

    # Try to load segment endpoints for mapping
    gL_end = None
    gR_end = None
    seg_pocket_path = seg_dir.parent / f"mep_seg_{seg_idx:02d}.pdb"
    if seg_pocket_path.exists():
        try:
            gL_end, gR_end = _load_segment_end_geoms(seg_pocket_path, [])
        except Exception as e:
            click.echo(f"[post] WARNING: failed to load segment endpoints: {e}", err=True)

    # Map IRC endpoints to left/right using bond-change analysis
    if gL_end is not None and gR_end is not None:
        bond_cfg = dict(_path_search.BOND_KW)

        def _matches(x, y) -> bool:
            try:
                chg, _ = _path_search._has_bond_change(x, y, bond_cfg)
                return not chg
            except Exception:
                return False

        def _rmsd_cart(g1, g2) -> float:
            c1 = np.asarray(g1.coords).reshape(-1, 3)
            c2 = np.asarray(g2.coords).reshape(-1, 3)
            n = min(len(c1), len(c2))
            return float(np.sqrt(np.mean((c1[:n] - c2[:n]) ** 2)))

        # Check if IRC endpoints need swapping
        match_LL = _matches(g_left, gL_end)
        match_LR = _matches(g_left, gR_end)
        match_RL = _matches(g_right, gL_end)
        match_RR = _matches(g_right, gR_end)

        if match_LR and match_RL and not (match_LL and match_RR):
            # Swap: IRC backward→right, forward→left
            g_left, g_right = g_right, g_left
            left_tag, right_tag = right_tag, left_tag
            reverse_irc = True
        elif not (match_LL and match_RR):
            # RMSD-based fallback
            d_LL = _rmsd_cart(g_left, gL_end)
            d_LR = _rmsd_cart(g_left, gR_end)
            d_RL = _rmsd_cart(g_right, gL_end)
            d_RR = _rmsd_cart(g_right, gR_end)
            if (d_LR + d_RL) < (d_LL + d_RR):
                g_left, g_right = g_right, g_left
                left_tag, right_tag = right_tag, left_tag
                reverse_irc = True

    return {
        "left_min_geom": g_left,
        "right_min_geom": g_right,
        "ts_geom": g_ts,
        "left_tag": left_tag,
        "right_tag": right_tag,
        "irc_trj": str(finished_trj) if finished_trj.exists() else None,
        "irc_plot": str(irc_plot) if irc_plot.exists() else None,
        "reverse_irc": reverse_irc,
    }


def _save_single_geom_for_tools(g: Any, ref_pdb: Path, out_dir: Path, name: str) -> Tuple[Path, Path]:
    """
    Write XYZ (primary, full precision) + PDB (companion) for a single geometry.
    Returns (xyz_path, pdb_path).
    """
    out_dir.mkdir(parents=True, exist_ok=True)
    # XYZ — full precision
    xyz_out = out_dir / f"{name}.xyz"
    with open(xyz_out, "w") as f:
        f.write(g.as_xyz() + "\n")
    # TRJ with energy (for PDB conversion and trajectory viewers)
    xyz_trj = out_dir / f"{name}_trj.xyz"
    _path_search._write_xyz_trj_with_energy([g], [float(g.energy)], xyz_trj)
    # PDB companion
    pdb_out = out_dir / f"{name}.pdb"
    _path_search._maybe_convert_to_pdb(xyz_trj, ref_pdb_path=ref_pdb, out_path=pdb_out)
    return xyz_out, pdb_out


def _save_single_geom_as_pdb_for_tools(g: Any, ref_pdb: Path, out_dir: Path, name: str) -> Path:
    """Backward-compatible wrapper: returns PDB path only."""
    _, pdb = _save_single_geom_for_tools(g, ref_pdb, out_dir, name)
    return pdb


def _run_tsopt_on_hei(hei_pdb: Path,
                      charge: int,
                      spin: int,
                      real_parm7: Path,
                      model_pdb: Path,
                      detect_layer: bool,
                      args_yaml: Optional[Path],
                      out_dir: Path,
                      opt_mode_default: str,
                      overrides: Optional[Dict[str, Any]] = None,
                      backend: Optional[str] = None,
                      embedcharge: bool = False,
                      embedcharge_cutoff: Optional[float] = None,
                      link_atom_method: Optional[str] = None,
                      mm_backend: Optional[str] = None,
                      use_cmap: Optional[bool] = None,
                      ref_pdb: Optional[Path] = None) -> Tuple[Path, Any]:
    """
    Run tsopt CLI on a HEI structure; return (final_ts_pdb_path, ts_geom).

    When *ref_pdb* (layered PDB with B-factor layer info) is given, the HEI XYZ
    is used as input and *ref_pdb* is passed via ``--ref-pdb`` so that the
    calculator correctly detects ML/MM layers from B-factors.
    """
    overrides = overrides or {}
    # Prefer HEI XYZ (full precision) + layered ref-pdb (B-factor layer info)
    hei_xyz = hei_pdb.with_suffix(".xyz")
    if ref_pdb is not None and hei_xyz.exists():
        input_file = hei_xyz
        topology_pdb = ref_pdb
    else:
        input_file = hei_pdb
        topology_pdb = hei_pdb
    prepared_input = prepare_input_structure(input_file)
    if input_file.suffix.lower() == ".xyz" and ref_pdb is not None:
        apply_ref_pdb_override(prepared_input, ref_pdb)
    try:
        ts_dir = _resolve_override_dir(out_dir / "ts", overrides.get("out_dir"))
        ensure_dir(ts_dir)

        opt_mode = overrides.get("opt_mode", opt_mode_default)

        ts_args: List[str] = [
            "-i", str(prepared_input.geom_path),
        ]
        if input_file.suffix.lower() == ".xyz" and ref_pdb is not None:
            ts_args.extend(["--ref-pdb", str(ref_pdb)])
        ts_args.extend([
            "--parm", str(real_parm7),
            "--model-pdb", str(model_pdb),
            "-q", str(int(charge)),
            "-m", str(int(spin)),
            "--out-dir", str(ts_dir),
        ])
        ts_args.append("--detect-layer" if detect_layer else "--no-detect-layer")

        if opt_mode is not None:
            ts_args.extend(["--opt-mode", str(opt_mode)])

        _append_cli_arg(ts_args, "--max-cycles", overrides.get("max_cycles"))
        _append_toggle_arg(ts_args, "--dump", overrides.get("dump"))
        _append_toggle_arg(ts_args, "--convert-files", overrides.get("convert_files"))
        _append_cli_arg(ts_args, "--thresh", overrides.get("thresh"))
        _append_toggle_arg(ts_args, "--flatten", overrides.get("flatten"))

        hess_mode = overrides.get("hessian_calc_mode")
        if hess_mode:
            ts_args.extend(["--hessian-calc-mode", str(hess_mode)])

        if args_yaml is not None:
            ts_args.extend(["--config", str(args_yaml)])

        if backend is not None:
            ts_args.extend(["--backend", str(backend)])
        if embedcharge:
            ts_args.append("--embedcharge")
            if embedcharge_cutoff is not None:
                ts_args.extend(["--embedcharge-cutoff", str(embedcharge_cutoff)])
        else:
            ts_args.append("--no-embedcharge")
        if link_atom_method is not None:
            ts_args.extend(["--link-atom-method", str(link_atom_method)])
        if mm_backend is not None:
            ts_args.extend(["--mm-backend", str(mm_backend)])
        if use_cmap is not None:
            ts_args.extend(["--cmap" if use_cmap else "--no-cmap"])
        if overrides.get("skip_final_freq"):
            ts_args.append("--skip-final-freq")

        _echo(f"[tsopt] Running tsopt on HEI → out={ts_dir}")
        _run_cli_main("tsopt", _ts_opt.cli, ts_args, on_nonzero="raise", prefix="tsopt")

        # Prefer XYZ (full precision) for geometry loading; PDB for topology
        final_xyz = ts_dir / "final_geometry.xyz"
        ts_pdb = ts_dir / "final_geometry.pdb"
        if not ts_pdb.exists() and final_xyz.exists():
            _path_search._maybe_convert_to_pdb(final_xyz, topology_pdb, ts_pdb)
        if not final_xyz.exists() and not ts_pdb.exists():
            raise click.ClickException("[tsopt] TS outputs not found.")
        geom_src = final_xyz if final_xyz.exists() else ts_pdb
        g_ts = geom_loader(geom_src, coord_type="cart")

        # Ensure calculator to have energy on g_ts
        _ts_calc_kwargs = dict(
            model_charge=int(charge),
            model_mult=int(spin),
            input_pdb=str(topology_pdb),
            real_parm7=str(real_parm7),
            model_pdb=str(model_pdb),
            use_bfactor_layers=detect_layer,
            backend=backend,
            embedcharge=embedcharge,
        )
        if link_atom_method is not None:
            _ts_calc_kwargs["link_atom_method"] = link_atom_method
        if mm_backend is not None:
            _ts_calc_kwargs["mm_backend"] = mm_backend
        if use_cmap is not None:
            _ts_calc_kwargs["use_cmap"] = use_cmap
        calc = _mlmm_calc(**_ts_calc_kwargs)
        g_ts.set_calculator(calc)
        _ = float(g_ts.energy)

        return ts_pdb, g_ts
    finally:
        prepared_input.cleanup()


def _pseudo_irc_and_match(seg_idx: int,
                          seg_dir: Path,
                          ref_pdb_for_seg: Path,
                          seg_pocket_pdb: Path,
                          g_ts: Any,
                          q_int: int,
                          spin: int,
                          real_parm7: Optional[Path] = None,
                          model_pdb: Optional[Path] = None,
                          detect_layer: bool = False,
                          backend: Optional[str] = None,
                          embedcharge: bool = False,
                          embedcharge_cutoff: Optional[float] = None,
                          link_atom_method: Optional[str] = None,
                          mm_backend: Optional[str] = None,
                          use_cmap: Optional[bool] = None) -> Dict[str, Any]:
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
    if backend is not None:
        calc_kwargs["backend"] = backend
    calc_kwargs["embedcharge"] = embedcharge
    if link_atom_method is not None:
        calc_kwargs["link_atom_method"] = link_atom_method
    if mm_backend is not None:
        calc_kwargs["mm_backend"] = mm_backend
    if use_cmap is not None:
        calc_kwargs["use_cmap"] = use_cmap
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
        trj_plus  = irc_dir / f"seg_{seg_idx:02d}_irc_plus_opt/optimization_trj.xyz"
        trj_minus = irc_dir / f"seg_{seg_idx:02d}_irc_minus_opt/optimization_trj.xyz"
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
                remain = [(t, gg) for (t, gg) in candidates if (mapping["left"] is None or mapping["left"][0] != t) and (mapping["right"] is None or mapping["right"][0] != t)]
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
            trj = irc_dir / f"irc_{side}_trj.xyz"
            _path_search._write_xyz_trj_with_energy([g_ts, gmin], [float(g_ts.energy), float(gmin.energy)], trj)
            run_trj2fig(trj, [irc_dir / f"irc_{side}_plot.png"], unit="kcal", reference="init", reverse_x=False)
    except Exception:
        logger.debug("Failed to write IRC mini-trajectory plot", exc_info=True)

    irc_trj = None
    irc_plot = None
    try:
        g_left = mapping["left"][1]
        g_right = mapping["right"][1]
        irc_trj = irc_dir / "finished_irc_trj.xyz"
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

    labels: List[str] = []
    for seg_idx in range(1, n_segments + 1):
        if seg_idx == 1:
            labels.extend(["R", "TS1", "IM1_1"])
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

    tmp_trj = out_png.with_name(f"{out_png.stem}_trj.xyz")
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
            logger.debug("Failed to unlink temp trajectory file", exc_info=True)


def _run_freq_for_state(pdb_path: Path,
                        q_int: int,
                        spin: int,
                        real_parm7: Path,
                        model_pdb: Path,
                        detect_layer: bool,
                        out_dir: Path,
                        args_yaml: Optional[Path],
                        overrides: Optional[Dict[str, Any]] = None,
                        backend: Optional[str] = None,
                        embedcharge: bool = False,
                        embedcharge_cutoff: Optional[float] = None,
                        link_atom_method: Optional[str] = None,
                        mm_backend: Optional[str] = None,
                        use_cmap: Optional[bool] = None,
                        xyz_path: Optional[Path] = None) -> Dict[str, Any]:
    """
    Run freq CLI; return parsed thermo dict (may be empty).
    When *xyz_path* is given, use it for full-precision coordinates with
    *pdb_path* as topology reference (--ref-pdb).
    """
    fdir = out_dir
    ensure_dir(fdir)
    overrides = overrides or {}

    dump_use = overrides.get("dump")
    if dump_use is None:
        dump_use = True

    # Prefer XYZ (full precision) with --ref-pdb for topology
    if xyz_path is not None and xyz_path.exists():
        args = ["-i", str(xyz_path), "--ref-pdb", str(pdb_path)]
    else:
        args = ["-i", str(pdb_path)]
    args.extend([
        "--parm", str(real_parm7),
        "--model-pdb", str(model_pdb),
        "-q", str(int(q_int)),
        "-m", str(int(spin)),
        "--out-dir", str(fdir),
    ])
    args.append("--detect-layer" if detect_layer else "--no-detect-layer")

    _append_cli_arg(args, "--max-write", overrides.get("max_write"))
    _append_cli_arg(args, "--amplitude-ang", overrides.get("amplitude_ang"))
    _append_cli_arg(args, "--n-frames", overrides.get("n_frames"))
    if overrides.get("sort") is not None:
        args.extend(["--sort", str(overrides.get("sort"))])
    _append_cli_arg(args, "--temperature", overrides.get("temperature"))
    _append_cli_arg(args, "--pressure", overrides.get("pressure"))
    _append_toggle_arg(args, "--dump", dump_use)
    _append_toggle_arg(args, "--convert-files", overrides.get("convert_files"))

    hess_mode = overrides.get("hessian_calc_mode")
    if hess_mode:
        args.extend(["--hessian-calc-mode", str(hess_mode)])

    if args_yaml is not None:
        args.extend(["--config", str(args_yaml)])
    if backend is not None:
        args.extend(["--backend", str(backend)])
    if embedcharge:
        args.append("--embedcharge")
        if embedcharge_cutoff is not None:
            args.extend(["--embedcharge-cutoff", str(embedcharge_cutoff)])
    else:
        args.append("--no-embedcharge")
    if link_atom_method is not None:
        args.extend(["--link-atom-method", str(link_atom_method)])
    if mm_backend is not None:
        args.extend(["--mm-backend", str(mm_backend)])
    if use_cmap is not None:
        args.extend(["--cmap" if use_cmap else "--no-cmap"])
    _run_cli_main("freq", _freq_cli.cli, args, on_nonzero="warn", on_exception="raise", prefix="freq")
    # parse thermoanalysis.yaml if any
    y = fdir / "thermoanalysis.yaml"
    if y.exists():
        try:
            return yaml.safe_load(y.read_text(encoding="utf-8")) or {}
        except Exception:
            return {}
    return {}


def _run_opt_for_state(
    pdb_path: Path,
    q_int: int,
    spin: int,
    real_parm7: Path,
    model_pdb: Path,
    detect_layer: bool,
    out_dir: Path,
    args_yaml: Optional[Path],
    opt_mode_default: str,
    convert_files: Optional[bool] = None,
    backend: Optional[str] = None,
    embedcharge: bool = False,
    embedcharge_cutoff: Optional[float] = None,
    link_atom_method: Optional[str] = None,
    mm_backend: Optional[str] = None,
    use_cmap: Optional[bool] = None,
    thresh: Optional[str] = None,
    xyz_path: Optional[Path] = None,
) -> Tuple[Any, Path]:
    """
    Run opt CLI for a single endpoint and return (optimized Geometry, final geometry path).
    When *xyz_path* is given, pass it as ``-i`` with ``--ref-pdb pdb_path`` to
    preserve full coordinate precision.
    """
    opt_dir = out_dir
    ensure_dir(opt_dir)

    # Use XYZ (full precision) when available; fall back to PDB
    if xyz_path is not None and xyz_path.exists():
        prepared_input = prepare_input_structure(xyz_path)
        apply_ref_pdb_override(prepared_input, pdb_path)
        input_label = xyz_path.name
    else:
        prepared_input = prepare_input_structure(pdb_path)
        input_label = pdb_path.name
    try:
        opt_mode = str(opt_mode_default or "heavy").lower()
        args = [
            "-i", str(prepared_input.geom_path),
        ]
        # Add --ref-pdb when input is XYZ
        if prepared_input.geom_path.suffix.lower() == ".xyz":
            args.extend(["--ref-pdb", str(prepared_input.source_path)])
        args.extend([
            "--parm", str(real_parm7),
            "--model-pdb", str(model_pdb),
            "-q", str(int(q_int)),
            "-m", str(int(spin)),
            "--out-dir", str(opt_dir),
            "--opt-mode", opt_mode,
        ])
        args.append("--detect-layer" if detect_layer else "--no-detect-layer")
        _append_toggle_arg(args, "--convert-files", convert_files)
        _append_cli_arg(args, "--thresh", thresh)

        if args_yaml is not None:
            args.extend(["--config", str(args_yaml)])

        if backend is not None:
            args.extend(["--backend", str(backend)])
        if embedcharge:
            args.append("--embedcharge")
            if embedcharge_cutoff is not None:
                args.extend(["--embedcharge-cutoff", str(embedcharge_cutoff)])
        else:
            args.append("--no-embedcharge")
        if link_atom_method is not None:
            args.extend(["--link-atom-method", str(link_atom_method)])
        if mm_backend is not None:
            args.extend(["--mm-backend", str(mm_backend)])
        if use_cmap is not None:
            args.extend(["--cmap" if use_cmap else "--no-cmap"])

        _echo(f"[endpoint-opt] Running opt on {input_label} (mode={opt_mode}) → out={opt_dir}")
        _run_cli_main("opt", _opt_cli.cli, args, on_nonzero="raise", on_exception="raise", prefix="endpoint-opt")

        final_pdb = opt_dir / "final_geometry.pdb"
        final_xyz = opt_dir / "final_geometry.xyz"
        # Prefer XYZ (full precision) for geometry loading
        if final_xyz.exists():
            final_geom_path = final_xyz
        elif final_pdb.exists():
            final_geom_path = final_pdb
        else:
            raise click.ClickException(f"[endpoint-opt] opt outputs not found under {opt_dir}")

        g_opt = geom_loader(final_geom_path, coord_type="cart")
        calc_input_pdb = final_pdb if final_pdb.exists() else pdb_path
        _opt_calc_kwargs = dict(
            model_charge=int(q_int),
            model_mult=int(spin),
            input_pdb=str(calc_input_pdb),
            real_parm7=str(real_parm7),
            model_pdb=str(model_pdb),
            use_bfactor_layers=detect_layer,
            backend=backend,
            embedcharge=embedcharge,
        )
        if link_atom_method is not None:
            _opt_calc_kwargs["link_atom_method"] = link_atom_method
        if mm_backend is not None:
            _opt_calc_kwargs["mm_backend"] = mm_backend
        if use_cmap is not None:
            _opt_calc_kwargs["use_cmap"] = use_cmap
        calc = _mlmm_calc(**_opt_calc_kwargs)
        g_opt.set_calculator(calc)
        _ = float(g_opt.energy)

        return g_opt, final_geom_path
    finally:
        prepared_input.cleanup()


def _run_dft_for_state(pdb_path: Path,
                       q_int: int,
                       spin: int,
                       real_parm7: Path,
                       model_pdb: Path,
                       detect_layer: bool,
                       out_dir: Path,
                       args_yaml: Optional[Path],
                       func_basis: str = "wb97m-v/def2-tzvpd",
                       overrides: Optional[Dict[str, Any]] = None,
                       backend: Optional[str] = None,
                       embedcharge: bool = False,
                       embedcharge_cutoff: Optional[float] = None,
                       link_atom_method: Optional[str] = None,
                       mm_backend: Optional[str] = None,
                       use_cmap: Optional[bool] = None,
                       xyz_path: Optional[Path] = None) -> Dict[str, Any]:
    """
    Run dft CLI; return parsed result.yaml dict (may be empty).
    When *xyz_path* is given, use it for full-precision coordinates with
    *pdb_path* as topology reference (--ref-pdb).
    """
    ddir = out_dir
    ensure_dir(ddir)
    overrides = overrides or {}

    func_basis_use = overrides.get("func_basis", func_basis)

    # Prefer XYZ (full precision) with --ref-pdb for topology
    if xyz_path is not None and xyz_path.exists():
        args = ["-i", str(xyz_path), "--ref-pdb", str(pdb_path)]
    else:
        args = ["-i", str(pdb_path)]
    args.extend([
        "--parm", str(real_parm7),
        "--model-pdb", str(model_pdb),
        "-q", str(int(q_int)),
        "-m", str(int(spin)),
        "--func-basis", str(func_basis_use),
        "--out-dir", str(ddir),
    ])
    args.append("--detect-layer" if detect_layer else "--no-detect-layer")

    _append_cli_arg(args, "--max-cycle", overrides.get("max_cycle"))
    _append_cli_arg(args, "--conv-tol", overrides.get("conv_tol"))
    _append_cli_arg(args, "--grid-level", overrides.get("grid_level"))
    _append_toggle_arg(args, "--convert-files", overrides.get("convert_files"))

    if args_yaml is not None:
        args.extend(["--config", str(args_yaml)])
    if backend is not None:
        args.extend(["--backend", str(backend)])
    if embedcharge:
        args.append("--embedcharge")
        if embedcharge_cutoff is not None:
            args.extend(["--embedcharge-cutoff", str(embedcharge_cutoff)])
    else:
        args.append("--no-embedcharge")
    if link_atom_method is not None:
        args.extend(["--link-atom-method", str(link_atom_method)])
    if mm_backend is not None:
        args.extend(["--mm-backend", str(mm_backend)])
    if use_cmap is not None:
        args.extend(["--cmap" if use_cmap else "--no-cmap"])
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
        "-l",
        "--ligand-charge",
        "-q",
        "--charge",
        "--out-dir",
        "--tsopt",
        "--thermo",
        "--dft",
        "--config",
        "--dry-run",
        "--embedcharge",
        "-s",
        "--scan-lists",
        "-b",
        "--backend",
        "-o",
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
    help="Run pocket extraction → (optional single-structure staged scan) → MEP search in one shot.\n"
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
    type=str, required=False, default=None,
    help=("Substrate specification for the extractor: "
          "a PDB path, a residue-ID list like '123,124' or 'A:123,B:456' "
          "(insertion codes OK: '123A' / 'A:123A'), "
          "or a residue-name list like 'GPP,MMT'. "
          "When omitted, extraction is skipped and full structures are used directly.")
)
@click.option(
    "-o", "--out-dir", "out_dir",
    type=click.Path(path_type=Path, file_okay=False),
    default=Path("./result_all/"), show_default=True,
    help="Top-level output directory for the pipeline."
)
# ===== Extractor knobs (subset of extract.parse_args) =====
@click.option("-r", "--radius", type=float, default=2.6, show_default=True,
              help="Inclusion cutoff (Å) around substrate atoms.")
@click.option("--radius-het2het", type=float, default=0.0, show_default=True,
              help="Independent hetero–hetero cutoff (Å) for non‑C/H pairs.")
@click.option("--include-h2o", "include_h2o", type=click.BOOL, default=True, show_default=True,
              help="Include waters (HOH/WAT/TIP3/SOL) in the pocket.")
@click.option("--exclude-backbone", "exclude_backbone", type=click.BOOL, default=False, show_default=True,
              help="Remove backbone atoms on non‑substrate amino acids (with PRO/HYP safeguards).")
@click.option("--add-linkh", "add_linkh", type=click.BOOL, default=False, show_default=True,
              help="Add link hydrogens for severed bonds (carbon-only) in pockets.")
@click.option("--selected-resn", type=str, default="", show_default=True,
              help="Force-include residues (comma/space separated; chain/insertion codes allowed).")
@click.option("-l", "--ligand-charge", type=str, default=None,
              help=("Either a total charge (number) to distribute across unknown residues "
                    "or a mapping like 'GPP:-3,MMT:-1'."))
@click.option(
    "-q",
    "--charge",
    "charge_override",
    type=int,
    default=None,
    help="Force total system charge. Highest priority over derived charges.",
)
@click.option(
    "--parm",
    "parm7_override",
    type=click.Path(path_type=Path, exists=True, dir_okay=False),
    default=None,
    help="Pre-built AMBER parm7 topology file. When provided, mm_parm generation is skipped.",
)
@click.option(
    "--model-pdb",
    "model_pdb_override",
    type=click.Path(path_type=Path, exists=True, dir_okay=False),
    default=None,
    help="Pre-built ML-region PDB (with B-factor layer info). When provided, ml_region generation is skipped.",
)
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
@click.option("--verbose", type=click.BOOL, default=True, show_default=True, help="Enable INFO-level logging inside extractor.")
# ===== Path search knobs (subset of path_search.cli) =====
@click.option("-m", "--multiplicity", "spin", type=int, default=1, show_default=True, help="Multiplicity (2S+1).")
@click.option("--max-nodes", type=int, default=_path_opt.GS_KW["max_nodes"], show_default=True,
              help="Max internal nodes for **segment** GSM (String has max_nodes+2 images including endpoints).")
@click.option("--max-cycles", type=int, default=300, show_default=True, help="Maximum GSM optimization cycles.")
@click.option("--climb", type=click.BOOL, default=True, show_default=True,
              help="Enable transition-state climbing after growth for the **first** segment in each pair.")
@click.option(
    "--opt-mode",
    type=click.Choice(["grad", "hess"], case_sensitive=False),
    default="grad",
    show_default=True,
    help=(
        "Optimizer mode forwarded to scan/path-search and used for single optimizations: "
        "grad (=LBFGS/Dimer) or hess (=RFO/RSIRFO)."
    ),
)
@click.option(
    "--opt-mode-post",
    type=click.Choice(["grad", "hess"], case_sensitive=False),
    default="hess",
    show_default=True,
    help=(
        "Optimizer mode for TSOPT and post-IRC endpoint optimizations. "
        "Takes precedence over --opt-mode for these stages."
    ),
)
@click.option("--dump", type=click.BOOL, default=False, show_default=True,
              help="Dump GSM / single-structure trajectories during the run, forwarding the same flag to scan/tsopt/freq.")
@click.option(
    "--refine-path/--no-refine-path",
    "refine_path",
    default=True,
    show_default=True,
    help=(
        "If True, run recursive path_search on the full ordered series; if False, run a single-pass "
        "path-opt GSM between each adjacent pair and concatenate the segments (no path_search)."
    ),
)
@click.option(
    "--thresh",
    type=str,
    default=None,
    show_default=False,
    help=(
        "Convergence preset (gau_loose|gau|gau_tight|gau_vtight|baker|never). "
        "Defaults to 'gau' when not provided."
    ),
)
@click.option(
    "--thresh-post",
    type=str,
    default="baker",
    show_default=True,
    help=(
        "Convergence preset for post-IRC endpoint optimizations "
        "(gau_loose|gau|gau_tight|gau_vtight|baker|never)."
    ),
)
@click.option("--config", "config_yaml", type=click.Path(path_type=Path, exists=True, dir_okay=False),
              default=None, help="Base YAML configuration file applied before explicit CLI options.")
@click.option("--show-config/--no-show-config", "show_config", default=False, show_default=True,
              help="Print resolved configuration and continue execution.")
@click.option("--dry-run/--no-dry-run", "dry_run", default=False, show_default=True,
              help="Validate options and print the execution plan without running any stage.")
@click.option("--preopt", "pre_opt", type=click.BOOL, default=True, show_default=True,
              help="If True, run initial single-structure optimizations of the pocket inputs.")
@click.option("--hessian-calc-mode",
              type=click.Choice(["Analytical", "FiniteDifference"], case_sensitive=False),
              default=None,
              help="Common MLIP Hessian calculation mode forwarded to tsopt and freq. Default: 'FiniteDifference'. Use 'Analytical' when VRAM is sufficient.")
@click.option(
    "--detect-layer/--no-detect-layer",
    "detect_layer",
    default=True,
    show_default=True,
    help="Detect ML/MM layers from input PDB B-factors (ML=0, MovableMM=10, FrozenMM=20) in downstream tools. "
         "If disabled, downstream tools require --model-pdb or --model-indices.",
)
# ===== Post-processing toggles =====
@click.option("--tsopt", "do_tsopt", type=click.BOOL, default=False, show_default=True,
              help="TS optimization + EulerPC IRC per reactive segment (or TSOPT-only mode for single-structure), and build energy diagrams.")
@click.option("--thermo", "do_thermo", type=click.BOOL, default=False, show_default=True,
              help="Run freq on (R,TS,P) per reactive segment (or TSOPT-only mode) and build Gibbs free-energy diagram (MLIP).")
@click.option("--dft", "do_dft", type=click.BOOL, default=False, show_default=True,
              help="Run DFT single-point on (R,TS,P) and build DFT energy diagram. With --thermo True, also generate a DFT//MLIP Gibbs diagram.")
@click.option("--tsopt-max-cycles", type=int, default=None,
              help="Override tsopt --max-cycles value.")
@click.option(
    "--flatten/--no-flatten",
    "flatten",
    default=False,
    show_default=True,
    help="Enable the extra-imaginary-mode flattening loop in tsopt (grad: dimer loop, hess: post-RSIRFO); --no-flatten forces flatten_max_iter=0.",
)
@click.option(
    "--skip-final-freq/--no-skip-final-freq",
    "skip_final_freq",
    default=False,
    show_default=True,
    help="Skip post-convergence frequency analysis in tsopt. Useful for large unfrozen systems.",
)
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
    "-s", "--scan-lists",
    "scan_lists_raw",
    type=str, multiple=True, required=False,
    help='Scan targets: inline Python literal or a YAML/JSON spec file path. '
         'Multiple inline literals define sequential stages, e.g. '
         '"[(12,45,1.35)]" "[(10,55,2.20),(23,34,1.80)]". '
         'Indices refer to the original full PDB (1-based) or PDB atom selectors like "TYR,285,CA"; '
         'they are auto-mapped to the pocket after extraction.',
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
@click.option("--convert-files/--no-convert-files", "convert_files", default=True, show_default=True,
              help="Convert XYZ/TRJ outputs to PDB format using reference topology; forwarded to all subcommands.")
@click.option(
    "--ref-pdb",
    "ref_pdb_cli",
    type=click.Path(path_type=Path, exists=True, dir_okay=False),
    default=None,
    help=(
        "Reference PDB for topology/B-factor layer information when -i provides XYZ inputs. "
        "Used for define-layer, mm_parm, ml_region, and forwarded to downstream tools "
        "(tsopt, irc, freq, path_search) as --ref-pdb."
    ),
)
@click.option(
    "-b", "--backend",
    type=click.Choice(["uma", "orb", "mace", "aimnet2"], case_sensitive=False),
    default=None,
    show_default=False,
    help="ML backend for the ONIOM high-level region (default: uma).",
)
@click.option(
    "--embedcharge/--no-embedcharge",
    "embedcharge",
    default=False,
    show_default=True,
    help="Enable xTB point-charge embedding correction for MM→ML environmental effects.",
)
@click.option(
    "--embedcharge-cutoff",
    "embedcharge_cutoff",
    type=float,
    default=None,
    show_default=False,
    help="Distance cutoff (Å) from ML region for MM point charges in xTB embedding. "
         "Default: 12.0 Å. Only used when --embedcharge is enabled.",
)
@click.option(
    "--link-atom-method",
    "link_atom_method",
    type=click.Choice(["scaled", "fixed"], case_sensitive=False),
    default=None,
    show_default=False,
    help="Link-atom position mode: scaled (g-factor, default) or fixed (legacy 1.09/1.01 Å).",
)
@click.option(
    "--mm-backend",
    "mm_backend",
    type=click.Choice(["hessian_ff", "openmm"], case_sensitive=False),
    default=None,
    show_default=False,
    help="MM backend: hessian_ff (analytical Hessian, default) or openmm (finite-difference Hessian, slower).",
)
@click.option(
    "--cmap/--no-cmap",
    "use_cmap",
    default=None,
    show_default=False,
    help="Enable CMAP (backbone cross-map) terms in model parm7. Default: disabled (Gaussian ONIOM-compatible).",
)
@click.pass_context
def cli(
    ctx: click.Context,
    input_paths: Sequence[Path],
    center_spec: Optional[str],
    out_dir: Path,
    radius: float,
    radius_het2het: float,
    include_h2o: bool,
    exclude_backbone: bool,
    add_linkh: bool,
    selected_resn: str,
    ligand_charge: Optional[str],
    charge_override: Optional[int],
    parm7_override: Optional[Path],
    model_pdb_override: Optional[Path],
    mm_ff_set: str,
    mm_add_ter: bool,
    mm_keep_temp: bool,
    mm_ligand_mult: Optional[str],
    verbose: bool,
    spin: int,
    max_nodes: int,
    max_cycles: int,
    climb: bool,
    opt_mode: str,
    opt_mode_post: Optional[str],
    dump: bool,
    refine_path: bool,
    thresh: Optional[str],
    thresh_post: str,
    config_yaml: Optional[Path],
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
    convert_files: bool,
    ref_pdb_cli: Optional[Path],
    backend: Optional[str],
    embedcharge: bool,
    embedcharge_cutoff: Optional[float],
    link_atom_method: Optional[str],
    mm_backend: Optional[str],
    use_cmap: Optional[bool],
    tsopt_max_cycles: Optional[int],
    flatten: bool,
    skip_final_freq: bool,
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
    _echo_state.reset()

    time_start = time.perf_counter()
    command_str = "mlmm all " + " ".join(sys.argv[1:])

    _is_param_explicit = make_is_param_explicit(ctx)
    dump_override_requested = _is_param_explicit("dump")
    opt_mode_set = _is_param_explicit("opt_mode")
    opt_mode_post_set = _is_param_explicit("opt_mode_post")

    config_yaml, override_yaml, _ = resolve_yaml_sources(config_yaml, None, None)
    args_yaml, merged_yaml_cfg = _build_effective_args_yaml(
        config_yaml=config_yaml,
        override_yaml=None,
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

    scan_vals = collect_single_option_values(argv_all, ("-s", "--scan-lists"), "--scan-lists")
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

    _mode_alias = {
        "grad": "grad",
        "hess": "hess",
        "light": "grad",
        "heavy": "hess",
    }
    opt_mode_norm = _mode_alias.get(str(opt_mode).strip().lower(), "grad")
    path_search_opt_mode = opt_mode_norm
    opt_mode_post_norm = (
        None
        if opt_mode_post is None
        else _mode_alias.get(str(opt_mode_post).strip().lower(), "hess")
    )
    endpoint_opt_mode_default = (
        opt_mode_post_norm if (opt_mode_post_set and opt_mode_post_norm is not None)
        else (opt_mode_norm if opt_mode_set else "hess")
    )
    if opt_mode_post_norm in {"grad", "hess"}:
        tsopt_opt_mode_default = opt_mode_post_norm
    elif opt_mode_set:
        tsopt_opt_mode_default = opt_mode_norm
    else:
        tsopt_opt_mode_default = "hess"
    tsopt_overrides: Dict[str, Any] = {}
    if tsopt_max_cycles is not None:
        tsopt_overrides["max_cycles"] = int(tsopt_max_cycles)
    if dump_override_requested:
        tsopt_overrides["dump"] = bool(dump)
    if tsopt_out_dir is not None:
        tsopt_overrides["out_dir"] = tsopt_out_dir
    if hessian_calc_mode is not None:
        tsopt_overrides["hessian_calc_mode"] = hessian_calc_mode
    if opt_mode_post_norm in {"grad", "hess"}:
        tsopt_overrides["opt_mode"] = opt_mode_post_norm
    elif opt_mode_set:
        tsopt_overrides["opt_mode"] = tsopt_opt_mode_default
    tsopt_overrides["convert_files"] = bool(convert_files)
    if thresh_post is not None:
        tsopt_overrides["thresh"] = str(thresh_post)
    if _is_param_explicit("flatten"):
        tsopt_overrides["flatten"] = bool(flatten)
    if skip_final_freq:
        tsopt_overrides["skip_final_freq"] = True

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
    freq_overrides["convert_files"] = bool(convert_files)

    dft_overrides: Dict[str, Any] = {}
    if dft_max_cycle is not None:
        dft_overrides["max_cycle"] = int(dft_max_cycle)
    if dft_conv_tol is not None:
        dft_overrides["conv_tol"] = float(dft_conv_tol)
    if dft_grid_level is not None:
        dft_overrides["grid_level"] = int(dft_grid_level)
    dft_overrides["convert_files"] = bool(convert_files)

    dft_func_basis_use = dft_func_basis or "wb97m-v/def2-tzvpd"
    dft_method_fallback = dft_func_basis_use

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
                "charge_override": charge_override,
                "skip_extract": bool(center_spec is None or str(center_spec).strip() == ""),
                "out_dir": str(out_dir),
                "spin": int(spin),
                "max_nodes": int(max_nodes),
                "max_cycles": int(max_cycles),
                "climb": bool(climb),
                "opt_mode": str(opt_mode),
                "opt_mode_post": (None if opt_mode_post is None else str(opt_mode_post)),
                "path_search_opt_mode": str(path_search_opt_mode),
                "endpoint_opt_mode": str(endpoint_opt_mode_default),
                "dump": bool(dump),
                "refine_path": bool(refine_path),
                "thresh": thresh,
                "thresh_post": thresh_post,
                "flatten": bool(flatten),
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
    path_dir = out_dir / ("path_search" if refine_path else "path_opt")
    scan_dir = _resolve_override_dir(out_dir / "scan", scan_out_dir)  # for single-structure scan mode
    ensure_dir(out_dir)
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
    skip_extract = center_spec is None or str(center_spec).strip() == ""

    # When inputs are XYZ and --ref-pdb is provided, use it for topology-requiring steps
    ref_pdb_for_topology: Optional[Path] = None
    if ref_pdb_cli is not None:
        ref_pdb_for_topology = ref_pdb_cli.resolve()
        _echo(f"[all] --ref-pdb provided: {ref_pdb_for_topology}")

    resolved_charge: Optional[int] = None
    pocket_outputs: List[Path] = []

    if skip_extract:
        _echo_section(
            "=== [all] Stage 1/3 — Extraction skipped (no -c/--center); using full structures as pockets ==="
        )
        pocket_outputs = [p.resolve() for p in extract_inputs]
        _echo("[all] Pocket inputs (full structures):")
        for op in pocket_outputs:
            _echo(f"  - {op}")
        # Use --model-pdb for charge derivation when provided (pocket charge),
        # otherwise fall back to full input PDB.
        charge_source_pdb = model_pdb_override if model_pdb_override is not None else extract_inputs[0]
        resolved_charge = _derive_charge_from_ligand_charge_when_extract_skipped(
            charge_source_pdb, ligand_charge
        )
    else:
        _echo_section(
            "=== [all] Stage 1/3 — Active-site pocket extraction (multi-structure union when applicable) ==="
        )
        ensure_dir(pockets_dir)
        for p in extract_inputs:
            pocket_outputs.append((pockets_dir / f"pocket_{p.stem}.pdb").resolve())

        try:
            ex_res = extract_api(
                complex_pdb=[str(p) for p in extract_inputs],
                center=center_spec,
                output=[str(p) for p in pocket_outputs],
                radius=float(radius),
                radius_het2het=float(radius_het2het),
                include_h2o=bool(include_h2o),
                exclude_backbone=bool(exclude_backbone),
                add_linkh=bool(add_linkh),
                selected_resn=selected_resn or "",
                ligand_charge=ligand_charge,
                verbose=bool(verbose),
            )
        except Exception as e:
            raise click.ClickException(f"[all] Extractor failed: {e}")

        _echo("[all] Pocket files:")
        for op in pocket_outputs:
            _echo(f"  - {op}")

        try:
            cs = ex_res.get("charge_summary", {})
            q_total = float(cs.get("total_charge", 0.0))
            q_prot = float(cs.get("protein_charge", 0.0))
            q_lig = float(cs.get("ligand_total_charge", 0.0))
            q_ion = float(cs.get("ion_total_charge", 0.0))
            _echo("")
            _echo("[all] Charge summary from extractor (model #1):")
            _echo(
                f"  Protein: {q_prot:+g},  Ligand: {q_lig:+g},  Ions: {q_ion:+g},  Total: {q_total:+g}"
            )
            resolved_charge = _round_charge_with_note(q_total)
        except Exception as e:
            raise click.ClickException(f"[all] Could not obtain total charge from extractor: {e}")

    if charge_override is not None:
        q_int = int(charge_override)
        override_msg = f"[all] WARNING: -q/--charge override supplied; forcing TOTAL system charge to {q_int:+d}"
        if resolved_charge is not None:
            override_msg += f" (would otherwise use {int(resolved_charge):+d} from workflow)"
        _echo(override_msg)
    else:
        if resolved_charge is None:
            raise click.ClickException(
                "[all] Total charge could not be resolved. Provide -q/--charge, "
                "or provide --ligand-charge when extraction is skipped."
            )
        q_int = int(resolved_charge)

    # --------------------------
    # Stage 1b: ML-region definition (copy first pocket) and mm_parm on the first full input
    # --------------------------
    _echo_section("=== [all] ML/MM preparation — ML region + parm7 ===")
    first_pocket = pocket_outputs[0]
    first_full_input = extract_inputs[0]
    # When --ref-pdb is provided, use it for PDB-requiring topology operations
    pocket_for_ml_region = ref_pdb_for_topology if ref_pdb_for_topology is not None else first_pocket
    pdb_for_mm_parm = ref_pdb_for_topology if ref_pdb_for_topology is not None else first_full_input

    # ML region definition: use --model-pdb if provided, otherwise generate from pocket
    if model_pdb_override is not None:
        ml_region_pdb = model_pdb_override.resolve()
        _echo(f"[all] ML region definition (--model-pdb override) → {ml_region_pdb}")
    else:
        ml_region_pdb = _write_ml_region_definition(pocket_for_ml_region, out_dir / "ml_region.pdb")
        _echo(f"[all] ML region definition → {ml_region_pdb}")

    # mm_parm: use --parm if provided, otherwise run tleap
    if parm7_override is not None:
        real_parm7_path = parm7_override.resolve()
        _echo(f"[all] parm7 (--parm override) → {real_parm7_path}")
    else:
        _echo(f"[all] mm_parm source PDB → {pdb_for_mm_parm}")
        mm_dir = out_dir / "mm_parm"
        ensure_commands_available(
            ("tleap", "antechamber", "parmchk2"),
            context="mm_parm (AmberTools)",
        )
        real_parm7_path, real_rst7_path = _build_mm_parm7(
            pdb=pdb_for_mm_parm,
            ligand_charge_expr=ligand_charge,
            ligand_mult_expr=mm_ligand_mult,
            out_dir=mm_dir,
            ff_set=mm_ff_set,
            add_ter=mm_add_ter,
            keep_temp=mm_keep_temp,
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
        # When --ref-pdb is given and input is not PDB, use ref_pdb for define-layer
        pdb_for_layer = full_pdb
        if ref_pdb_for_topology is not None and full_pdb.suffix.lower() != ".pdb":
            pdb_for_layer = ref_pdb_for_topology
        out_layered = layered_dir / f"{pdb_for_layer.stem}_layered.pdb"
        try:
            layer_info = _define_layers(
                input_pdb=pdb_for_layer,
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
    irc_trj_for_all: List[Tuple[Path, bool]] = []

    if single_tsopt_mode:
        _echo_section("=== [all] TSOPT-only single-structure mode ===")
        tsroot = out_dir / "tsopt_single"
        ensure_dir(tsroot)

        # Use the layered full-system PDB as TS initial guess
        layered_pdb = layered_inputs[0]
        # When --ref-pdb is given and input is XYZ, copy the XYZ next to the layered PDB
        # so that _run_tsopt_on_hei can use XYZ (full precision) + layered PDB (topology)
        if ref_pdb_for_topology is not None and extract_inputs[0].suffix.lower() != ".pdb":
            xyz_companion = layered_pdb.with_suffix(".xyz")
            if not xyz_companion.exists():
                shutil.copy2(extract_inputs[0], xyz_companion)
                _echo(f"[all] Copied XYZ input → {xyz_companion} (full precision for tsopt)")
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
            backend=backend,
            embedcharge=embedcharge,
            embedcharge_cutoff=embedcharge_cutoff,
            link_atom_method=link_atom_method,
            mm_backend=mm_backend,
            use_cmap=use_cmap,
            ref_pdb=layered_pdb,
        )

        # EulerPC IRC & map endpoints (no segment endpoints exist → fallback mapping)
        irc_pocket_ref = ref_pdb_for_topology if ref_pdb_for_topology is not None else first_pocket
        irc_res = _irc_and_match(seg_idx=1,
                                 seg_dir=tsroot,
                                 ref_pdb_for_seg=ts_pdb,
                                 seg_pocket_pdb=irc_pocket_ref,
                                 g_ts=g_ts,
                                 q_int=q_int,
                                 spin=spin,
                                 real_parm7=real_parm7_path,
                                 model_pdb=ml_region_pdb,
                                 detect_layer=detect_layer,
                                 backend=backend,
                                 embedcharge=embedcharge,
                                 embedcharge_cutoff=embedcharge_cutoff,
                                 link_atom_method=link_atom_method,
                                 mm_backend=mm_backend,
                                 use_cmap=use_cmap,
                                 args_yaml=args_yaml)
        gL = irc_res["left_min_geom"]
        gR = irc_res["right_min_geom"]
        gT = irc_res["ts_geom"]
        irc_plot_path = irc_res.get("irc_plot")
        irc_trj_path = irc_res.get("irc_trj")
        if irc_trj_path:
            try:
                irc_trj_for_all.append((Path(irc_trj_path), bool(irc_res.get("reverse_irc", False))))
            except Exception:
                logger.debug("Failed to append IRC trajectory path", exc_info=True)

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

        # Save XYZ (full precision) + PDB (companion) and run endpoint-opt
        struct_dir = tsroot / "structures"
        ensure_dir(struct_dir)
        pocket_ref = ref_pdb_for_topology if ref_pdb_for_topology is not None else first_pocket
        xR_irc, pR_irc = _save_single_geom_for_tools(g_react, pocket_ref, struct_dir, "reactant_irc")
        xT, pT         = _save_single_geom_for_tools(gT,       pocket_ref, struct_dir, "ts")
        xP_irc, pP_irc = _save_single_geom_for_tools(g_prod,   pocket_ref, struct_dir, "product_irc")

        endpoint_opt_dir = tsroot / "endpoint_opt"
        ensure_dir(endpoint_opt_dir)

        # Map IRC left/right Hessians → R/P endpoint (left=forward, right=backward)
        from .hessian_cache import load as _hess_load, store as _hess_store, clear as _clear_hess_cache
        _react_hk = "irc_left" if eL >= eR else "irc_right"
        _prod_hk  = "irc_right" if eL >= eR else "irc_left"

        _c = _hess_load(_react_hk)
        if _c:
            _hess_store("irc_endpoint", _c["hessian"], active_dofs=_c.get("active_dofs"), meta=_c.get("meta"))
        try:
            g_react, _ = _run_opt_for_state(
                pR_irc, q_int, spin, real_parm7_path, ml_region_pdb, detect_layer,
                endpoint_opt_dir / "R", args_yaml, endpoint_opt_mode_default,
                convert_files=convert_files,
                backend=backend,
                embedcharge=embedcharge,
                embedcharge_cutoff=embedcharge_cutoff,
                link_atom_method=link_atom_method,
                mm_backend=mm_backend,
                use_cmap=use_cmap,
                thresh=thresh_post,
                xyz_path=xR_irc,
            )
        except Exception as e:
            _echo(
                f"[post] WARNING: Reactant endpoint optimization failed in TSOPT-only mode: {e}",
                err=True,
            )

        _c = _hess_load(_prod_hk)
        if _c:
            _hess_store("irc_endpoint", _c["hessian"], active_dofs=_c.get("active_dofs"), meta=_c.get("meta"))
        try:
            g_prod, _ = _run_opt_for_state(
                pP_irc, q_int, spin, real_parm7_path, ml_region_pdb, detect_layer,
                endpoint_opt_dir / "P", args_yaml, endpoint_opt_mode_default,
                convert_files=convert_files,
                backend=backend,
                embedcharge=embedcharge,
                embedcharge_cutoff=embedcharge_cutoff,
                link_atom_method=link_atom_method,
                mm_backend=mm_backend,
                use_cmap=use_cmap,
                thresh=thresh_post,
                xyz_path=xP_irc,
            )
        except Exception as e:
            _echo(
                f"[post] WARNING: Product endpoint optimization failed in TSOPT-only mode: {e}",
                err=True,
            )
        shutil.rmtree(endpoint_opt_dir, ignore_errors=True)
        _echo("[endpoint-opt] Clean endpoint-opt working dir.")

        xR, pR = _save_single_geom_for_tools(g_react, pocket_ref, struct_dir, "reactant")
        xP, pP = _save_single_geom_for_tools(g_prod,   pocket_ref, struct_dir, "product")
        e_react = float(g_react.energy)
        e_prod = float(g_prod.energy)

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

        # ── Release GPU memory before freq/thermo/DFT ──
        for _g in (gL, gR, gT, g_react, g_prod):
            if _g is not None and hasattr(_g, "calculator"):
                _g.calculator = None
        gc.collect()
        if torch.cuda.is_available():
            torch.cuda.empty_cache()

        # Thermochemistry (UMA) Gibbs
        thermo_payloads: Dict[str, Dict[str, Any]] = {}
        GR = GT = GP = None
        eR_dft = eT_dft = eP_dft = None
        GR_dftUMA = GT_dftUMA = GP_dftUMA = None
        freq_root = _resolve_override_dir(tsroot / "freq", freq_out_dir)
        dft_root = _resolve_override_dir(tsroot / "dft", dft_out_dir)

        if do_thermo:
            _echo(f"[thermo] Single TSOPT: freq on TS/R/P")
            tT = _run_freq_for_state(pT, q_int, spin, real_parm7_path, ml_region_pdb, detect_layer,
                                     freq_root / "TS", args_yaml, overrides=freq_overrides,
                                     backend=backend, embedcharge=embedcharge, embedcharge_cutoff=embedcharge_cutoff,
                                     link_atom_method=link_atom_method, mm_backend=mm_backend, use_cmap=use_cmap, xyz_path=xT)
            _clear_hess_cache()  # TS Hessian consumed; R/P need exact computation
            tR = _run_freq_for_state(pR, q_int, spin, real_parm7_path, ml_region_pdb, detect_layer,
                                     freq_root / "R", args_yaml, overrides=freq_overrides,
                                     backend=backend, embedcharge=embedcharge, embedcharge_cutoff=embedcharge_cutoff,
                                     link_atom_method=link_atom_method, mm_backend=mm_backend, use_cmap=use_cmap, xyz_path=xR)
            tP = _run_freq_for_state(pP, q_int, spin, real_parm7_path, ml_region_pdb, detect_layer,
                                     freq_root / "P", args_yaml, overrides=freq_overrides,
                                     backend=backend, embedcharge=embedcharge, embedcharge_cutoff=embedcharge_cutoff,
                                     link_atom_method=link_atom_method, mm_backend=mm_backend, use_cmap=use_cmap, xyz_path=xP)
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
                                     dft_root / "R", args_yaml, func_basis=dft_func_basis_use, overrides=dft_overrides,
                                     backend=backend, embedcharge=embedcharge, embedcharge_cutoff=embedcharge_cutoff,
                                     link_atom_method=link_atom_method, mm_backend=mm_backend, use_cmap=use_cmap, xyz_path=xR)
            dT = _run_dft_for_state(pT, q_int, spin, real_parm7_path, ml_region_pdb, detect_layer,
                                     dft_root / "TS", args_yaml, func_basis=dft_func_basis_use, overrides=dft_overrides,
                                     backend=backend, embedcharge=embedcharge, embedcharge_cutoff=embedcharge_cutoff,
                                     link_atom_method=link_atom_method, mm_backend=mm_backend, use_cmap=use_cmap, xyz_path=xT)
            dP = _run_dft_for_state(pP, q_int, spin, real_parm7_path, ml_region_pdb, detect_layer,
                                     dft_root / "P", args_yaml, func_basis=dft_func_basis_use, overrides=dft_overrides,
                                     backend=backend, embedcharge=embedcharge, embedcharge_cutoff=embedcharge_cutoff,
                                     link_atom_method=link_atom_method, mm_backend=mm_backend, use_cmap=use_cmap, xyz_path=xP)
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
            "refine_path": bool(refine_path),
            "thresh": thresh,
            "thresh_post": thresh_post,
            "flatten": bool(flatten),
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
            "key_files": {},
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

        # summary.md and key_* outputs are disabled.
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
        scan_endopt_use = False if scan_endopt_override is None else bool(scan_endopt_override)
        scan_opt_mode_use = path_search_opt_mode

        scan_args: List[str] = [
            "-i", str(layered_pdb),
            "--parm", str(real_parm7_path),
            "-q", str(int(q_int)),
            "-m", str(int(spin)),
            "--out-dir", str(scan_dir),
            "--preopt" if scan_preopt_use else "--no-preopt",
            "--endopt" if scan_endopt_use else "--no-endopt",
            "--opt-mode", str(scan_opt_mode_use),
        ]
        scan_args.append("--detect-layer" if detect_layer else "--no-detect-layer")

        if dump_override_requested:
            scan_args.append("--dump" if dump else "--no-dump")

        if scan_one_based is not None:
            scan_args.append("--one-based" if scan_one_based else "--zero-based")

        _append_cli_arg(scan_args, "--max-step-size", scan_max_step_size)
        _append_cli_arg(scan_args, "--bias-k", scan_bias_k)
        _append_cli_arg(scan_args, "--relax-max-cycles", scan_relax_max_cycles)
        scan_args.append("--convert-files" if convert_files else "--no-convert-files")
        if thresh is not None:
            scan_args.extend(["--thresh", str(thresh)])
        if args_yaml is not None:
            scan_args.extend(["--config", str(args_yaml)])
        # Forward all converted --scan-lists (aligned to the pocket atom order)
        if scan_stage_literals:
            scan_args.append("--scan-lists")
            scan_args.extend(scan_stage_literals)

        if backend is not None:
            scan_args.extend(["--backend", str(backend)])
        if embedcharge:
            scan_args.append("--embedcharge")
            if embedcharge_cutoff is not None:
                scan_args.extend(["--embedcharge-cutoff", str(embedcharge_cutoff)])
        else:
            scan_args.append("--no-embedcharge")
        if link_atom_method is not None:
            scan_args.extend(["--link-atom-method", str(link_atom_method)])
        if mm_backend is not None:
            scan_args.extend(["--mm-backend", str(mm_backend)])
        if use_cmap is not None:
            scan_args.extend(["--cmap" if use_cmap else "--no-cmap"])

        _echo("[all] mlmm scan " + " ".join(scan_args))

        _run_cli_main("scan", _scan_cli.cli, scan_args, on_nonzero="raise", on_exception="raise", prefix="all")

        # Collect stage results — prefer XYZ (full precision), keep PDB as ref for topology
        stage_results: List[Path] = []
        stage_refs: List[Path] = []
        for st in sorted(scan_dir.glob("stage_*")):
            if not st.is_dir():
                continue
            xyz = st / "result.xyz"
            pdb = st / "result.pdb"
            if xyz.exists():
                stage_results.append(xyz.resolve())
                stage_refs.append(pdb.resolve() if pdb.exists() else layered_pdb)
            elif pdb.exists():
                stage_results.append(pdb.resolve())
                stage_refs.append(pdb.resolve())
        if not stage_results:
            raise click.ClickException("[all] No stage result files found under scan/.")
        _echo("[all] Collected scan stage files:")
        for p in stage_results:
            _echo(f"  - {p}")

        # Input series to path_search: [preopt result (if available), scan stage results ...]
        # When scan ran with --preopt, its optimised reactant geometry lives in
        # scan/preopt/result.xyz (full precision) or result.pdb.  Using this
        # avoids a redundant ~2000-cycle re-optimisation inside path_search.
        preopt_xyz = scan_dir / "preopt" / "result.xyz"
        preopt_pdb = scan_dir / "preopt" / "result.pdb"
        if preopt_xyz.exists():
            init0_geom = preopt_xyz.resolve()
            init0_ref = layered_pdb          # layered PDB has authoritative B-factor layers
        elif preopt_pdb.exists():
            init0_geom = preopt_pdb.resolve()
            init0_ref = layered_pdb
        else:
            # No preopt output — fall back to original layered PDB
            init0_geom = layered_pdb
            init0_ref = layered_pdb
        pockets_for_path = [init0_geom] + stage_results
        refs_for_path = [init0_ref] + stage_refs
    else:
        # Multi-structure standard route: use layered full-system PDBs
        pockets_for_path = list(layered_inputs)

    # --------------------------
    # Stage 2: Path search on full-system layered PDBs
    # --------------------------
    if refine_path:
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
        ps_args.extend(["--parm", str(real_parm7_path)])
        # Layered PDBs have B-factors → detect-layer will auto-identify layers
        ps_args.append("--detect-layer")

        # Nodes, cycles, climb, optimizer, dump, out-dir, preopt, args-yaml
        ps_args.extend(["--max-nodes", str(int(max_nodes))])
        ps_args.extend(["--max-cycles", str(int(max_cycles))])
        ps_args.append("--climb" if climb else "--no-climb")
        ps_args.extend(["--opt-mode", str(path_search_opt_mode)])
        ps_args.append("--dump" if dump else "--no-dump")
        ps_args.extend(["--out-dir", str(path_dir)])
        ps_args.append("--preopt" if pre_opt else "--no-preopt")
        ps_args.append("--convert-files" if convert_files else "--no-convert-files")
        if thresh is not None:
            ps_args.extend(["--thresh", str(thresh)])
        if args_yaml is not None:
            ps_args.extend(["--config", str(args_yaml)])

        # Provide --ref-pdb for topology/B-factor info (one per input)
        # MUST use layered PDBs (with B-factor layer info) so that downstream
        # PDB conversion preserves ML/MovableMM/FrozenMM layer encoding.
        ps_args.append("--ref-pdb")
        if is_single and has_scan:
            # single+scan: use refs_for_path which maps to each pocket (XYZ→PDB ref)
            for ref in refs_for_path:
                ps_args.append(str(ref))
        else:
            for lp in layered_inputs:
                ps_args.append(str(lp))

        if backend is not None:
            ps_args.extend(["--backend", str(backend)])
        if embedcharge:
            ps_args.append("--embedcharge")
            if embedcharge_cutoff is not None:
                ps_args.extend(["--embedcharge-cutoff", str(embedcharge_cutoff)])
        else:
            ps_args.append("--no-embedcharge")
        if link_atom_method is not None:
            ps_args.extend(["--link-atom-method", str(link_atom_method)])
        if mm_backend is not None:
            ps_args.extend(["--mm-backend", str(mm_backend)])
        if use_cmap is not None:
            ps_args.extend(["--cmap" if use_cmap else "--no-cmap"])

        _echo("[all] mlmm path-search " + " ".join(ps_args))

        _run_cli_main("path_search", _path_search.cli, ps_args, on_nonzero="raise", on_exception="raise", prefix="all")
    else:
        # --no-refine-path: run path-opt GSM between each adjacent pair and concatenate
        _echo_section("=== [all] Stage 2/3 — MEP path-opt on full-system layered PDBs (single-pass GSM per pair) ===")

        if len(pockets_for_path) < 2:
            raise click.ClickException("[all] Need at least two structures for path-opt MEP concatenation.")

        ensure_dir(path_dir)
        combined_blocks: List[str] = []
        path_opt_segments: List[Dict[str, Any]] = []

        for pair_idx in range(len(pockets_for_path) - 1):
            p_left = pockets_for_path[pair_idx]
            p_right = pockets_for_path[pair_idx + 1]
            seg_tag = f"seg_{pair_idx:02d}"
            seg_out = path_dir / f"{seg_tag}_mep"
            ensure_dir(seg_out)

            po_args: List[str] = [
                "-i", str(p_left), str(p_right),
                "-q", str(q_int),
                "-m", str(int(spin)),
                "--parm", str(real_parm7_path),
                "--detect-layer",
                "--max-nodes", str(int(max_nodes)),
                "--max-cycles", str(int(max_cycles)),
            ]
            po_args.append("--climb" if climb else "--no-climb")
            po_args.append("--dump" if dump else "--no-dump")
            po_args.extend(["--out-dir", str(seg_out)])
            po_args.append("--preopt" if pre_opt else "--no-preopt")
            po_args.append("--convert-files" if convert_files else "--no-convert-files")
            if thresh is not None:
                po_args.extend(["--thresh", str(thresh)])
            if args_yaml is not None:
                po_args.extend(["--config", str(args_yaml)])
            if backend is not None:
                po_args.extend(["--backend", str(backend)])
            if embedcharge:
                po_args.append("--embedcharge")
                if embedcharge_cutoff is not None:
                    po_args.extend(["--embedcharge-cutoff", str(embedcharge_cutoff)])
            else:
                po_args.append("--no-embedcharge")
            if link_atom_method is not None:
                po_args.extend(["--link-atom-method", str(link_atom_method)])
            if mm_backend is not None:
                po_args.extend(["--mm-backend", str(mm_backend)])
            if use_cmap is not None:
                po_args.extend(["--cmap" if use_cmap else "--no-cmap"])

            _echo(f"[all] mlmm path-opt " + " ".join(po_args))
            _run_cli_main("path_opt", _path_opt.cli, po_args, on_nonzero="raise", on_exception="raise", prefix="all")

            # --- Post-processing per segment ---
            seg_trj = seg_out / "final_geometries_trj.xyz"
            if not seg_trj.exists():
                raise click.ClickException(
                    f"[all] path-opt segment {pair_idx} did not produce final_geometries_trj.xyz"
                )

            # Copy per-segment trajectory to path_dir
            try:
                seg_mep_trj = path_dir / f"mep_seg_{pair_idx:02d}_trj.xyz"
                shutil.copy2(seg_trj, seg_mep_trj)
                if pockets_for_path[0].suffix.lower() == ".pdb":
                    _path_search._maybe_convert_to_pdb(
                        seg_mep_trj,
                        ref_pdb_path=pockets_for_path[0],
                        out_path=path_dir / f"mep_seg_{pair_idx:02d}.pdb",
                    )
            except Exception as e:
                _echo(
                    f"[all] WARNING: failed to emit per-segment trajectory copies for segment {pair_idx:02d}: {e}",
                    err=True,
                )

            # Mirror HEI artifacts
            hei_src = seg_out / "hei.xyz"
            if hei_src.exists():
                try:
                    shutil.copy2(hei_src, path_dir / f"hei_seg_{pair_idx:02d}.xyz")
                    hei_pdb_src = seg_out / "hei.pdb"
                    if hei_pdb_src.exists():
                        shutil.copy2(hei_pdb_src, path_dir / f"hei_seg_{pair_idx:02d}.pdb")
                except Exception as e:
                    _echo(
                        f"[all] WARNING: failed to prepare HEI artifacts for segment {pair_idx:02d}: {e}",
                        err=True,
                    )

            # Parse trajectory blocks for concatenation and energy extraction
            raw_blocks = read_xyz_as_blocks(seg_trj, strict=True)
            blocks = ["\n".join(b) + "\n" for b in raw_blocks]
            if not blocks:
                raise click.ClickException(
                    f"[all] No frames read from path-opt segment {pair_idx} trajectory: {seg_trj}"
                )
            # Skip duplicate first frame for subsequent segments
            if pair_idx > 0:
                blocks = blocks[1:]
            combined_blocks.extend(blocks)

            # Extract energies from trajectory comment lines
            energies_seg: List[float] = []
            for blk in raw_blocks:
                E = np.nan
                if len(blk) >= 2:
                    try:
                        E = float(blk[1].split()[0])
                    except Exception:
                        E = np.nan
                energies_seg.append(E)

            # Parse first/last frame coordinates for bond-change detection
            first_last = None
            try:
                first_last = xyz_blocks_first_last(raw_blocks, path=seg_trj)
            except Exception as e:
                _echo(
                    f"[all] WARNING: failed to parse first/last frames for segment {pair_idx:02d}: {e}",
                    err=True,
                )

            path_opt_segments.append(
                {
                    "tag": seg_tag,
                    "energies": energies_seg,
                    "traj": seg_trj,
                    "inputs": (p_left, p_right),
                    "first_last": first_last,
                }
            )

        # --- Concatenated MEP trajectory ---
        final_trj = path_dir / "mep_trj.xyz"
        try:
            final_trj.write_text("".join(combined_blocks), encoding="utf-8")
            _echo(f"[all] Wrote concatenated MEP trajectory: {final_trj}")
        except Exception as e:
            raise click.ClickException(f"[all] Failed to write concatenated MEP: {e}")

        # Energy plot for concatenated trajectory
        try:
            run_trj2fig(final_trj, [path_dir / "mep_plot.png"], unit="kcal", reference="init", reverse_x=False)
            close_matplotlib_figures()
            _echo(f"[plot] Saved energy plot → '{path_dir / 'mep_plot.png'}'")
        except Exception as e:
            _echo(f"[plot] WARNING: Failed to plot concatenated MEP: {e}", err=True)

        # PDB conversion of concatenated trajectory
        try:
            if pockets_for_path[0].suffix.lower() == ".pdb":
                mep_pdb_path = path_dir / "mep.pdb"
                _path_search._maybe_convert_to_pdb(
                    final_trj, ref_pdb_path=pockets_for_path[0], out_path=mep_pdb_path
                )
                if mep_pdb_path.exists():
                    shutil.copy2(mep_pdb_path, out_dir / mep_pdb_path.name)
                    _echo(f"[all] Copied concatenated MEP PDB → {out_dir / mep_pdb_path.name}")
        except Exception as e:
            _echo(
                f"[all] WARNING: Failed to convert/copy concatenated MEP to PDB: {e}",
                err=True,
            )

        # --- Energy diagram ---
        energy_diagrams_po: List[Dict[str, Any]] = []
        try:
            labels = _build_global_segment_labels(len(path_opt_segments))
            energies_chain: List[float] = []
            for si, seg_info in enumerate(path_opt_segments):
                Es = [float(x) for x in seg_info.get("energies", [])]
                if not Es:
                    continue
                if si == 0:
                    energies_chain.append(Es[0])
                energies_chain.append(float(np.nanmax(Es)))
                energies_chain.append(Es[-1])
            if labels and energies_chain and len(labels) == len(energies_chain):
                title_note = "(GSM; all segments)" if len(path_opt_segments) > 1 else "(GSM)"
                diag_payload = _write_segment_energy_diagram(
                    path_dir / "energy_diagram_mep",
                    labels=labels,
                    energies_eh=energies_chain,
                    title_note=title_note,
                )
                if diag_payload:
                    energy_diagrams_po.append(diag_payload)
        except Exception as e:
            _echo(f"[diagram] WARNING: Failed to build GSM diagram for path-opt branch: {e}", err=True)

        # --- Bond change detection and summary.yaml ---
        segments_summary: List[Dict[str, Any]] = []
        bond_cfg = dict(_path_search.BOND_KW)
        for seg_idx, info in enumerate(path_opt_segments):
            Es = [float(x) for x in info.get("energies", []) if np.isfinite(x)]
            if not Es:
                continue
            barrier = (max(Es) - Es[0]) * AU2KCALPERMOL
            delta = (Es[-1] - Es[0]) * AU2KCALPERMOL
            bond_summary = ""
            try:
                first_last = info.get("first_last")
                if first_last:
                    elems, c_first, c_last = first_last
                else:
                    elems, c_first, c_last = read_xyz_first_last(Path(info["traj"]))
                gL = _geom_from_angstrom(elems, c_first, [])
                gR = _geom_from_angstrom(elems, c_last, [])
                changed, bond_summary = _path_search._has_bond_change(gL, gR, bond_cfg)
                if not changed:
                    bond_summary = "(no covalent changes detected)"
            except Exception as e:
                _echo(
                    f"[all] WARNING: Failed to detect bond changes for segment {seg_idx:02d}: {e}",
                    err=True,
                )
                bond_summary = "(no covalent changes detected)"

            segments_summary.append(
                {
                    "index": seg_idx,
                    "tag": info.get("tag", f"seg_{seg_idx:02d}"),
                    "kind": "seg",
                    "barrier_kcal": float(barrier),
                    "delta_kcal": float(delta),
                    "bond_changes": bond_summary,
                }
            )

        po_summary: Dict[str, Any] = {
            "out_dir": str(path_dir),
            "n_images": len(read_xyz_as_blocks(final_trj)),
            "n_segments": len(segments_summary),
            "segments": segments_summary,
        }
        if energy_diagrams_po:
            po_summary["energy_diagrams"] = list(energy_diagrams_po)
        try:
            with open(path_dir / "summary.yaml", "w") as f:
                yaml.safe_dump(po_summary, f, sort_keys=False, allow_unicode=True)
            _echo(f"[write] Wrote '{path_dir / 'summary.yaml'}'.")
        except Exception as e:
            _echo(f"[write] WARNING: Failed to write summary.yaml for path-opt branch: {e}", err=True)

        # Copy key outputs to out_dir root
        try:
            for name in ("mep_plot.png", "energy_diagram_mep.png", "summary.yaml"):
                src = path_dir / name
                if src.exists():
                    shutil.copy2(src, out_dir / name)
            for ext in ("_trj.xyz", ".xyz"):
                src = path_dir / f"mep{ext}"
                if src.exists():
                    shutil.copy2(src, out_dir / src.name)
        except Exception as e:
            _echo(f"[all] WARNING: Failed to relocate path-opt summary files: {e}", err=True)

    # --------------------------
    # Stage 3: Merge (performed by path_search when --ref-pdb was supplied)
    # --------------------------
    _echo_section("=== [all] Stage 3/3 — Final outputs ===")
    _echo(f"[all] Final products can be found under: {path_dir}")
    _echo("  - mep_trj.xyz              (concatenated MEP trajectory)")
    _echo("  - mep.pdb                  (PDB conversion, if input was .pdb)")
    _echo("  - mep_seg_XX_trj.xyz       (per-segment trajectories)")
    _echo("  - hei_seg_XX.xyz/.pdb      (HEI per segment)")
    _echo("  - summary.yaml             (segment barriers, ΔE, bond changes)")
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
                "summary.yaml",
                "summary.log",
            ):
                src = path_dir / name
                if src.exists():
                    shutil.copy2(src, out_dir / name)
            for stem in ("mep",):
                for ext in ("_trj.xyz", ".xyz"):
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
                "pipeline_mode": "path-search" if refine_path else "path-opt",
                "refine_path": bool(refine_path),
                "thresh": thresh,
                "thresh_post": thresh_post,
                "flatten": bool(flatten),
                "tsopt": do_tsopt,
                "thermo": do_thermo,
                "dft": do_dft,
                "opt_mode": opt_mode_norm,
                "opt_mode_post": opt_mode_post.lower() if opt_mode_post else None,
                "mep_mode": "path-search" if refine_path else "path-opt",
                "uma_model": None,
                "command": command_str,
                "charge": q_int,
                "spin": spin,
                "mep": mep_info,
                "segments": summary.get("segments", []),
                "energy_diagrams": summary.get("energy_diagrams", []),
                "post_segments": list(post_segment_logs),
                "key_files": {},
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
        # summary.md and key_* outputs are disabled.
        # Elapsed time
        _echo(format_elapsed("[all] Elapsed for Whole Pipeline", time_start))
        return

    _echo_section("=== [all] Stage 4 — Post-processing per reactive segment ===")

    # Use segment summary from path_search / path-opt
    if not segments:
        _echo("[post] No segments found in summary; nothing to do.")
        _write_pipeline_summary_log([])
        # summary.md and key_* outputs are disabled.
        _echo(format_elapsed("[all] Elapsed for Whole Pipeline", time_start))
        return

    # Iterate only bond-change segments (kind='seg' and bond_changes not empty and not '(no covalent...)')
    reactive = [s for s in segments if (s.get("kind", "seg") == "seg" and str(s.get("bond_changes", "")).strip() and str(s.get("bond_changes", "")).strip() != "(no covalent changes detected)")]
    if not reactive:
        _echo("[post] No bond-change segments. Skipping TS/thermo/DFT.")
        _write_pipeline_summary_log([])
        # summary.md and key_* outputs are disabled.
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
                backend=backend,
                embedcharge=embedcharge,
                embedcharge_cutoff=embedcharge_cutoff,
                link_atom_method=link_atom_method,
                mm_backend=mm_backend,
                use_cmap=use_cmap,
                ref_pdb=layered_inputs[0] if layered_inputs else None,
            )
        else:
            # If TSOPT off: use the GSM HEI (pocket) as TS geometry
            ts_pdb = hei_pocket_pdb
            g_ts = geom_loader(ts_pdb, coord_type="cart")
            _hei_calc_kwargs = dict(
                model_charge=int(q_int),
                model_mult=int(spin),
                input_pdb=str(ts_pdb),
                real_parm7=str(real_parm7_path),
                model_pdb=str(ml_region_pdb),
                use_bfactor_layers=detect_layer,
                backend=backend,
                embedcharge=embedcharge,
            )
            if link_atom_method is not None:
                _hei_calc_kwargs["link_atom_method"] = link_atom_method
            if mm_backend is not None:
                _hei_calc_kwargs["mm_backend"] = mm_backend
            if use_cmap is not None:
                _hei_calc_kwargs["use_cmap"] = use_cmap
            calc = _mlmm_calc(**_hei_calc_kwargs)
            g_ts.set_calculator(calc); _ = float(g_ts.energy)

        # 4.2 EulerPC IRC & mapping to (left,right)
        irc_plot_path = None
        irc_trj_path = None
        irc_res = _irc_and_match(seg_idx=seg_idx,
                                 seg_dir=seg_dir,
                                 ref_pdb_for_seg=ts_pdb,
                                 seg_pocket_pdb=hei_pocket_pdb,
                                 g_ts=g_ts,
                                 q_int=q_int,
                                 spin=spin,
                                 real_parm7=real_parm7_path,
                                 model_pdb=ml_region_pdb,
                                 detect_layer=detect_layer,
                                 backend=backend,
                                 embedcharge=embedcharge,
                                 embedcharge_cutoff=embedcharge_cutoff,
                                 link_atom_method=link_atom_method,
                                 mm_backend=mm_backend,
                                 use_cmap=use_cmap,
                                 args_yaml=args_yaml)
        irc_plot_path = irc_res.get("irc_plot")
        irc_trj_path = irc_res.get("irc_trj")
        if irc_trj_path:
            try:
                irc_trj_for_all.append((Path(irc_trj_path), bool(irc_res.get("reverse_irc", False))))
            except Exception:
                logger.debug("Failed to append IRC trajectory path", exc_info=True)

        gL = irc_res["left_min_geom"]
        gR = irc_res["right_min_geom"]
        gT = irc_res["ts_geom"]
        # Save IRC endpoints (XYZ primary), run endpoint-opt, then save optimized structures
        struct_dir = seg_dir / "structures"
        ensure_dir(struct_dir)
        xL_irc, pL_irc = _save_single_geom_for_tools(gL, hei_pocket_pdb, struct_dir, "reactant_irc")
        xT, pT         = _save_single_geom_for_tools(gT, hei_pocket_pdb, struct_dir, "ts")
        xR_irc, pR_irc = _save_single_geom_for_tools(gR, hei_pocket_pdb, struct_dir, "product_irc")

        endpoint_opt_dir = seg_dir / "endpoint_opt"
        ensure_dir(endpoint_opt_dir)

        # Map IRC left/right Hessians → R/P endpoint
        # When reverse_irc is True, _irc_and_match swapped left/right to match GSM endpoints,
        # so "irc_left" (=forward) now corresponds to gR and "irc_right" (=backward) to gL.
        from .hessian_cache import load as _hess_load, store as _hess_store, clear as _clear_hess_cache
        _reversed = bool(irc_res.get("reverse_irc", False))
        _left_hk  = "irc_right" if _reversed else "irc_left"
        _right_hk = "irc_left"  if _reversed else "irc_right"

        _c = _hess_load(_left_hk)
        if _c:
            _hess_store("irc_endpoint", _c["hessian"], active_dofs=_c.get("active_dofs"), meta=_c.get("meta"))
        try:
            gL, _ = _run_opt_for_state(
                pL_irc, q_int, spin, real_parm7_path, ml_region_pdb, detect_layer,
                endpoint_opt_dir / "R", args_yaml, endpoint_opt_mode_default,
                convert_files=convert_files,
                backend=backend,
                embedcharge=embedcharge,
                embedcharge_cutoff=embedcharge_cutoff,
                link_atom_method=link_atom_method,
                mm_backend=mm_backend,
                use_cmap=use_cmap,
                thresh=thresh_post,
                xyz_path=xL_irc,
            )
        except Exception as e:
            _echo(
                f"[post] WARNING: Reactant endpoint optimization failed for segment {seg_idx:02d}: {e}",
                err=True,
            )

        _c = _hess_load(_right_hk)
        if _c:
            _hess_store("irc_endpoint", _c["hessian"], active_dofs=_c.get("active_dofs"), meta=_c.get("meta"))
        try:
            gR, _ = _run_opt_for_state(
                pR_irc, q_int, spin, real_parm7_path, ml_region_pdb, detect_layer,
                endpoint_opt_dir / "P", args_yaml, endpoint_opt_mode_default,
                convert_files=convert_files,
                backend=backend,
                embedcharge=embedcharge,
                embedcharge_cutoff=embedcharge_cutoff,
                link_atom_method=link_atom_method,
                mm_backend=mm_backend,
                use_cmap=use_cmap,
                thresh=thresh_post,
                xyz_path=xR_irc,
            )
        except Exception as e:
            _echo(
                f"[post] WARNING: Product endpoint optimization failed for segment {seg_idx:02d}: {e}",
                err=True,
            )
        shutil.rmtree(endpoint_opt_dir, ignore_errors=True)
        _echo("[endpoint-opt] Clean endpoint-opt working dir.")

        xL, pL = _save_single_geom_for_tools(gL, hei_pocket_pdb, struct_dir, "reactant")
        xR, pR = _save_single_geom_for_tools(gR, hei_pocket_pdb, struct_dir, "product")

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

        # ── Release GPU memory before freq/thermo/DFT ──
        for _g in (gL, gR, gT):
            if _g is not None and hasattr(_g, "calculator"):
                _g.calculator = None
        gc.collect()
        if torch.cuda.is_available():
            torch.cuda.empty_cache()

        # 4.4 Thermochemistry (UMA freq) and Gibbs diagram
        thermo_payloads: Dict[str, Dict[str, Any]] = {}
        GR = GT = GP = None
        freq_seg_root = _resolve_override_dir(seg_dir / "freq", freq_out_dir)
        dft_seg_root = _resolve_override_dir(seg_dir / "dft", dft_out_dir)

        if do_thermo:
            _echo(f"[thermo] Segment {seg_idx:02d}: freq on TS/R/P")
            tT = _run_freq_for_state(pT, q_int, spin, real_parm7_path, ml_region_pdb, detect_layer,
                                     freq_seg_root / "TS", args_yaml, overrides=freq_overrides,
                                     backend=backend, embedcharge=embedcharge, embedcharge_cutoff=embedcharge_cutoff,
                                     link_atom_method=link_atom_method, mm_backend=mm_backend, use_cmap=use_cmap, xyz_path=xT)
            _clear_hess_cache()  # TS Hessian consumed; R/P need exact computation
            tR = _run_freq_for_state(pL, q_int, spin, real_parm7_path, ml_region_pdb, detect_layer,
                                     freq_seg_root / "R", args_yaml, overrides=freq_overrides,
                                     backend=backend, embedcharge=embedcharge, embedcharge_cutoff=embedcharge_cutoff,
                                     link_atom_method=link_atom_method, mm_backend=mm_backend, use_cmap=use_cmap, xyz_path=xL)
            tP = _run_freq_for_state(pR, q_int, spin, real_parm7_path, ml_region_pdb, detect_layer,
                                     freq_seg_root / "P", args_yaml, overrides=freq_overrides,
                                     backend=backend, embedcharge=embedcharge, embedcharge_cutoff=embedcharge_cutoff,
                                     link_atom_method=link_atom_method, mm_backend=mm_backend, use_cmap=use_cmap, xyz_path=xR)
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
                                     dft_seg_root / "R", args_yaml, func_basis=dft_func_basis_use, overrides=dft_overrides,
                                     backend=backend, embedcharge=embedcharge, embedcharge_cutoff=embedcharge_cutoff,
                                     link_atom_method=link_atom_method, mm_backend=mm_backend, use_cmap=use_cmap, xyz_path=xL)
            dT = _run_dft_for_state(pT, q_int, spin, real_parm7_path, ml_region_pdb, detect_layer,
                                     dft_seg_root / "TS", args_yaml, func_basis=dft_func_basis_use, overrides=dft_overrides,
                                     backend=backend, embedcharge=embedcharge, embedcharge_cutoff=embedcharge_cutoff,
                                     link_atom_method=link_atom_method, mm_backend=mm_backend, use_cmap=use_cmap, xyz_path=xT)
            dP = _run_dft_for_state(pR, q_int, spin, real_parm7_path, ml_region_pdb, detect_layer,
                                     dft_seg_root / "P", args_yaml, func_basis=dft_func_basis_use, overrides=dft_overrides,
                                     backend=backend, embedcharge=embedcharge, embedcharge_cutoff=embedcharge_cutoff,
                                     link_atom_method=link_atom_method, mm_backend=mm_backend, use_cmap=use_cmap, xyz_path=xR)
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
    # summary.md and key_* outputs are disabled.
    _echo(format_elapsed("[all] Elapsed for Whole Pipeline", time_start))


_configure_all_help_visibility(cli)


if __name__ == "__main__":
    cli()
