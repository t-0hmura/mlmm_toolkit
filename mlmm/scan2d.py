"""ML/MM two-distance (d1, d2) grid scan with harmonic restraints.

Example:
    mlmm scan2d -i input.pdb --parm real.parm7 --model-pdb ml_region.pdb \
        -q 0 --scan-lists "[(12,45,1.30,3.10),(10,55,1.20,3.20)]"

For detailed documentation, see: docs/scan2d.md
"""

from __future__ import annotations

import functools
from copy import deepcopy
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence

import gc
import logging
import math
import shutil
import sys
import textwrap
import traceback
import tempfile
import time

logger = logging.getLogger(__name__)

import click
import numpy as np
import torch
import pandas as pd
from scipy.interpolate import Rbf
import plotly.graph_objects as go

from pysisyphus.helpers import geom_loader
from pysisyphus.optimizers.LBFGS import LBFGS
from pysisyphus.optimizers.exceptions import OptimizationError, ZeroStepLength
from pysisyphus.constants import ANG2BOHR, AU2KCALPERMOL

from .mlmm_calc import mlmm
from .defaults import BIAS_KW as _BIAS_KW_DEFAULT, MLMM_CALC_KW
from .opt import (
    GEOM_KW as _OPT_GEOM_KW,
    CALC_KW as _OPT_CALC_KW,
    OPT_BASE_KW as _OPT_BASE_KW,
    LBFGS_KW as _OPT_LBFGS_KW,
    HarmonicBiasCalculator,
    _parse_freeze_atoms,
    _normalize_geom_freeze,
)
from .opt import _convert_yaml_layer_atoms_1to0
from .utils import (
    apply_ref_pdb_override,
    apply_layer_freeze_constraints,
    set_convert_file_enabled,
    deep_update,
    load_yaml_dict,
    apply_yaml_overrides,
    pretty_block,
    strip_inherited_keys,
    filter_calc_for_echo,
    format_freeze_atoms_for_echo,
    format_elapsed,
    merge_freeze_atom_indices,
    prepare_input_structure,
    resolve_charge_spin_or_raise,
    convert_xyz_to_pdb,
    load_pdb_atom_metadata,
    parse_scan_list_quads,
    parse_scan_spec_quads,
    is_scan_spec_file,
    axis_label_csv,
    axis_label_html,
    PDB_ATOM_META_HEADER,
    format_pdb_atom_metadata,
    parse_indices_string,
    build_model_pdb_from_bfactors,
    build_model_pdb_from_indices,
    ensure_dir,
    distance_A_from_coords,
    distance_tag,
    values_from_bounds,
    unbiased_energy_hartree,
    snapshot_geometry,
    convert_and_annotate_xyz_to_pdb,
)
from .cli_utils import resolve_yaml_sources, load_merged_yaml_cfg, make_is_param_explicit

# Shared defaults (copied from opt.py to keep ML/MM behaviour consistent)
GEOM_KW: Dict[str, Any] = deepcopy(_OPT_GEOM_KW)
CALC_KW: Dict[str, Any] = deepcopy(_OPT_CALC_KW)
OPT_BASE_KW: Dict[str, Any] = deepcopy(_OPT_BASE_KW)
OPT_BASE_KW.update(
    {
        "out_dir": "./result_scan2d/",
        "dump": False,        # Keep LBFGS runs light; per-grid TRJs are handled separately via --dump
        "max_cycles": 10000,  # Overridden per relaxation through --relax-max-cycles
    }
)
LBFGS_KW: Dict[str, Any] = deepcopy(_OPT_LBFGS_KW)
LBFGS_KW.update({"out_dir": "./result_scan2d/"})
BIAS_KW: Dict[str, Any] = deepcopy(_BIAS_KW_DEFAULT)


_snapshot_geometry = functools.partial(snapshot_geometry, coord_type_default="cart")


def _select_closest_state(
    states: Sequence[Dict[str, Any]],
    d1_target: float,
    d2_target: float,
):
    """
    Return the Geometry from `states` whose (d1_A, d2_A) is closest to the target.

    Each state is a dict with at least the keys "d1_A", "d2_A" and "geom".
    """
    if not states:
        return None
    best_state = min(
        states,
        key=lambda st: math.hypot(st["d1_A"] - d1_target, st["d2_A"] - d2_target),
    )
    return best_state.get("geom")


def _select_closest_state_1d(
    states: Sequence[Dict[str, Any]],
    d1_target: float,
):
    """
    Return the Geometry from `states` whose d1_A is closest to the target.

    Used for choosing the initial structure for the d1-only biased relaxation.
    """
    if not states:
        return None
    best_state = min(
        states,
        key=lambda st: abs(st["d1_A"] - d1_target),
    )
    return best_state.get("geom")


def _make_lbfgs(
    geom,
    lbfgs_cfg: Dict[str, Any],
    opt_cfg: Dict[str, Any],
    *,
    max_step_bohr: float,
    relax_max_cycles: int,
    out_dir: Path,
    prefix: str,
) -> LBFGS:
    common = dict(opt_cfg)
    common["out_dir"] = str(out_dir)
    common["prefix"] = prefix
    args = {**lbfgs_cfg, **common}
    args["max_step"] = min(float(lbfgs_cfg.get("max_step", 0.30)), max_step_bohr)
    args["max_cycles"] = int(relax_max_cycles)
    return LBFGS(geom, **args)


@click.command(
    help="2D distance scan with harmonic restraints using the ML/MM calculator.",
    context_settings={"help_option_names": ["-h", "--help"]},
)
@click.option(
    "-i",
    "--input",
    "input_path",
    type=click.Path(path_type=Path, exists=True, dir_okay=False),
    required=True,
    help="Input enzyme complex PDB (required).",
)
@click.option(
    "--parm",
    "real_parm7",
    type=click.Path(path_type=Path, exists=True, dir_okay=False),
    required=True,
    help="Amber parm7 topology for the enzyme (required).",
)
@click.option(
    "--model-pdb",
    type=click.Path(path_type=Path, exists=True, dir_okay=False),
    required=False,
    help="PDB defining the ML region. Optional when --detect-layer is enabled.",
)
@click.option(
    "--model-indices",
    "model_indices_str",
    type=str,
    default=None,
    show_default=False,
    help="Comma-separated atom indices for the ML region (ranges allowed like 1-5). "
         "Used when --model-pdb is omitted.",
)
@click.option(
    "--model-indices-one-based/--model-indices-zero-based",
    "model_indices_one_based",
    default=True,
    show_default=True,
    help="Interpret --model-indices as 1-based (default) or 0-based.",
)
@click.option(
    "--detect-layer/--no-detect-layer",
    "detect_layer",
    default=True,
    show_default=True,
    help="Detect ML/MM layers from input PDB B-factors (ML=0, MovableMM=10, FrozenMM=20). "
         "If disabled, you must provide --model-pdb or --model-indices.",
)
@click.option("-q", "--charge", type=int, required=False,
              help="ML-region total charge. Required unless --ligand-charge is provided.")
@click.option("-l", "--ligand-charge", type=str, default=None, show_default=False,
              help="Total charge or per-resname mapping (e.g., GPP:-3,SAM:1) used to derive "
                   "charge when -q is omitted (requires PDB input or --ref-pdb).")
@click.option(
    "-m",
    "--multiplicity",
    "spin",
    type=int,
    default=None,
    show_default=False,
    help="Spin multiplicity (2S+1) for the ML region. Defaults to 1 when omitted.",
)
@click.option(
    "--freeze-atoms",
    "freeze_atoms_cli",
    type=str,
    default=None,
    show_default=False,
    help='Comma-separated 1-based atom indices to freeze (e.g., "1,3,5").',
)
@click.option(
    "--hess-cutoff",
    "hess_cutoff",
    type=float,
    default=None,
    show_default=False,
    help="Distance cutoff (Å) from ML region for MM atoms to include in Hessian calculation. "
         "Applied to movable MM atoms and can be combined with --detect-layer.",
)
@click.option(
    "--movable-cutoff",
    "movable_cutoff",
    type=float,
    default=None,
    show_default=False,
    help="Distance cutoff (Å) from ML region for movable MM atoms. MM atoms beyond this are frozen. "
         "Providing --movable-cutoff disables --detect-layer.",
)
@click.option(
    "-s", "--scan-lists",
    "scan_list_raw",
    type=str,
    required=False,
    help="Scan targets: inline Python literal or a YAML/JSON spec file path.",
)
@click.option(
    "--one-based/--zero-based",
    "one_based",
    default=True,
    show_default=True,
    help="Interpret (i,j) indices in --scan-lists as 1-based (default) or 0-based.",
)
@click.option(
    "--print-parsed/--no-print-parsed",
    "print_parsed",
    default=False,
    show_default=True,
    help="Print parsed scan targets after resolving --scan-lists.",
)
@click.option(
    "--max-step-size",
    type=float,
    default=0.20,
    show_default=True,
    help="Maximum spacing between successive distance targets [Å].",
)
@click.option("--bias-k", type=float, default=300.0, show_default=True, help="Harmonic well strength k [eV/Å^2].")
@click.option(
    "--relax-max-cycles",
    type=int,
    default=10000,
    show_default=True,
    help="Maximum LBFGS cycles per biased relaxation (also used for preopt).",
)
@click.option(
    "--dump/--no-dump",
    "dump",
    default=False,
    show_default=True,
    help="Write inner d2 scan TRJs per d1 slice.",
)
@click.option(
    "-o", "--out-dir",
    type=str,
    default="./result_scan2d/",
    show_default=True,
    help="Base output directory.",
)
@click.option(
    "--thresh",
    type=click.Choice(["gau_loose", "gau", "gau_tight", "gau_vtight", "baker", "never"], case_sensitive=False),
    default="baker",
    show_default=True,
    help="Convergence preset.",
)
@click.option(
    "--config",
    "config_yaml",
    type=click.Path(path_type=Path, exists=True, dir_okay=False),
    default=None,
    show_default=False,
    help="Base YAML configuration file applied before explicit CLI options.",
)
@click.option(
    "--ref-pdb",
    type=click.Path(path_type=Path, exists=True, dir_okay=False),
    default=None,
    help="Reference PDB topology to use when --input is XYZ (keeps XYZ coordinates).",
)
@click.option(
    "--preopt/--no-preopt",
    "preopt",
    default=False,
    show_default=True,
    help="Run an unbiased pre-optimization.",
)
@click.option(
    "--baseline",
    type=click.Choice(["min", "first"]),
    default="min",
    show_default=True,
    help="Reference for relative energy (kcal/mol): 'min' or 'first' (i=0,j=0).",
)
@click.option(
    "--zmin",
    type=float,
    default=None,
    show_default=False,
    help="Lower bound of the contour color scale (kcal/mol).",
)
@click.option(
    "--zmax",
    type=float,
    default=None,
    show_default=False,
    help="Upper bound of the contour color scale (kcal/mol).",
)
@click.option(
    "--convert-files/--no-convert-files",
    "convert_files",
    default=True,
    show_default=True,
    help="Convert XYZ/TRJ outputs into PDB companions based on the input format.",
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
@click.option(
    "--out-json/--no-out-json",
    "out_json",
    default=False,
    show_default=True,
    help="Write machine-readable result.json to out_dir.",
)
@click.pass_context
def cli(
    ctx: click.Context,
    input_path: Path,
    real_parm7: Path,
    model_pdb: Optional[Path],
    model_indices_str: Optional[str],
    model_indices_one_based: bool,
    detect_layer: bool,
    charge: Optional[int],
    ligand_charge: Optional[str],
    spin: Optional[int],
    freeze_atoms_cli: Optional[str],
    hess_cutoff: Optional[float],
    movable_cutoff: Optional[float],
    scan_list_raw: Optional[str],
    one_based: bool,
    print_parsed: bool,
    max_step_size: float,
    bias_k: float,
    relax_max_cycles: int,
    dump: bool,
    out_dir: str,
    thresh: Optional[str],
    config_yaml: Optional[Path],
    ref_pdb: Optional[Path],
    preopt: bool,
    baseline: str,
    zmin: Optional[float],
    zmax: Optional[float],
    convert_files: bool,
    backend: Optional[str],
    embedcharge: bool,
    embedcharge_cutoff: Optional[float],
    link_atom_method: Optional[str],
    mm_backend: Optional[str],
    use_cmap: Optional[bool],
    out_json: bool,
) -> None:
    _is_param_explicit = make_is_param_explicit(ctx)

    set_convert_file_enabled(convert_files)
    time_start = time.perf_counter()
    config_yaml, override_yaml, used_legacy_yaml = resolve_yaml_sources(
        config_yaml=config_yaml,
        override_yaml=None,
        args_yaml_legacy=None,
    )

    # Validate input format: PDB directly, or XYZ with --ref-pdb
    suffix = input_path.suffix.lower()
    if suffix not in (".pdb", ".xyz"):
        click.echo("ERROR: --input must be a PDB or XYZ file.", err=True)
        sys.exit(1)
    if suffix == ".xyz" and ref_pdb is None:
        click.echo("ERROR: --ref-pdb is required when --input is an XYZ file.", err=True)
        sys.exit(1)

    tmp_root = None
    try:
        with prepare_input_structure(input_path) as prepared_input:
            try:
                apply_ref_pdb_override(prepared_input, ref_pdb)
            except click.BadParameter as e:
                click.echo(f"ERROR: {e}", err=True)
                sys.exit(1)
            geom_input_path = prepared_input.geom_path
            source_path = prepared_input.source_path
            charge, spin = resolve_charge_spin_or_raise(
                prepared_input, charge, spin,
                ligand_charge=ligand_charge, prefix="[scan2d]",
            )

            try:
                freeze_atoms_list = _parse_freeze_atoms(freeze_atoms_cli)
            except click.BadParameter as exc:
                click.echo(f"ERROR: {exc}", err=True)
                sys.exit(1)

            model_indices: Optional[List[int]] = None
            if model_indices_str:
                try:
                    model_indices = parse_indices_string(model_indices_str, one_based=model_indices_one_based)
                except click.BadParameter as exc:
                    click.echo(f"ERROR: {exc}", err=True)
                    sys.exit(1)

            yaml_cfg, _, _ = load_merged_yaml_cfg(
                config_yaml=config_yaml,
                override_yaml=None,
            )

            geom_cfg = dict(GEOM_KW)
            calc_cfg = dict(CALC_KW)
            opt_cfg = dict(OPT_BASE_KW)
            lbfgs_cfg = dict(LBFGS_KW)
            bias_cfg = dict(BIAS_KW)

            apply_yaml_overrides(
                yaml_cfg,
                [
                    (geom_cfg, (("geom",),)),
                    (calc_cfg, (("calc",), ("mlmm",))),
                    (opt_cfg, (("opt",),)),
                    (lbfgs_cfg, (("lbfgs",), ("opt", "lbfgs"))),
                    (bias_cfg, (("bias",),)),
                ],
            )

            try:
                geom_freeze = _normalize_geom_freeze(geom_cfg.get("freeze_atoms"))
            except click.BadParameter as exc:
                click.echo(f"ERROR: {exc}", err=True)
                sys.exit(1)
            geom_cfg["freeze_atoms"] = geom_freeze
            _convert_yaml_layer_atoms_1to0(calc_cfg)
            if freeze_atoms_list:
                merge_freeze_atom_indices(geom_cfg, freeze_atoms_list)
            freeze_atoms_final = list(geom_cfg.get("freeze_atoms") or [])
            calc_cfg["freeze_atoms"] = freeze_atoms_final

            opt_cfg["out_dir"] = out_dir
            opt_cfg["dump"] = False
            opt_cfg["max_cycles"] = int(relax_max_cycles)
            if thresh is not None:
                opt_cfg["thresh"] = str(thresh)
            lbfgs_cfg["max_cycles"] = int(relax_max_cycles)
            if bias_k is not None:
                bias_cfg["k"] = float(bias_k)

            out_dir_path = Path(opt_cfg["out_dir"]).resolve()

            calc_cfg["model_charge"] = int(charge)
            calc_cfg["model_mult"] = int(spin)
            calc_cfg["input_pdb"] = str(source_path)
            calc_cfg["real_parm7"] = str(real_parm7)
            if backend is not None:
                calc_cfg["backend"] = str(backend).lower()
            if _is_param_explicit("embedcharge"):
                calc_cfg["embedcharge"] = bool(embedcharge)
            if _is_param_explicit("embedcharge_cutoff"):
                calc_cfg["embedcharge_cutoff"] = embedcharge_cutoff
            if link_atom_method is not None:
                calc_cfg["link_atom_method"] = str(link_atom_method).lower()
            if mm_backend is not None:
                calc_cfg["mm_backend"] = str(mm_backend).lower()
            if use_cmap is not None:
                calc_cfg["use_cmap"] = use_cmap

            # movable_cutoff implies full distance-based layer assignment.
            # hess_cutoff alone can be combined with --detect-layer.
            if movable_cutoff is not None:
                if detect_layer:
                    click.echo("[layer] --movable-cutoff provided; disabling --detect-layer.", err=True)
                detect_layer = False

            layer_source_pdb = source_path
            if detect_layer and layer_source_pdb.suffix.lower() != ".pdb":
                click.echo("ERROR: --detect-layer requires a PDB input (or --ref-pdb).", err=True)
                sys.exit(1)

            model_pdb_path: Optional[Path] = None
            layer_info: Optional[Dict[str, List[int]]] = None

            if detect_layer:
                try:
                    model_pdb_path, layer_info = build_model_pdb_from_bfactors(layer_source_pdb, out_dir_path)
                    calc_cfg["use_bfactor_layers"] = True
                    click.echo(
                        f"[layer] Detected B-factor layers: ML={len(layer_info.get('ml_indices', []))}, "
                        f"MovableMM={len(layer_info.get('movable_mm_indices', []))}, "
                        f"FrozenMM={len(layer_info.get('frozen_indices', []))}"
                    )
                except Exception as e:
                    if model_pdb is None and not model_indices:
                        click.echo(f"ERROR: {e}", err=True)
                        sys.exit(1)
                    click.echo(f"[layer] WARNING: {e} Falling back to explicit ML region.", err=True)
                    detect_layer = False

            if not detect_layer:
                if model_pdb is None and not model_indices:
                    click.echo("ERROR: Provide --model-pdb or --model-indices when --no-detect-layer.", err=True)
                    sys.exit(1)
                if model_pdb is not None:
                    model_pdb_path = Path(model_pdb)
                else:
                    if layer_source_pdb.suffix.lower() != ".pdb":
                        click.echo("ERROR: --model-indices requires a PDB input (or --ref-pdb).", err=True)
                        sys.exit(1)
                    try:
                        model_pdb_path = build_model_pdb_from_indices(layer_source_pdb, out_dir_path, model_indices or [])
                    except Exception as e:
                        click.echo(f"ERROR: {e}", err=True)
                        sys.exit(1)
                calc_cfg["use_bfactor_layers"] = False

            if model_pdb_path is None:
                click.echo("ERROR: Failed to resolve model PDB for the ML region.", err=True)
                sys.exit(1)

            calc_cfg["model_pdb"] = str(model_pdb_path)
            freeze_atoms_final = apply_layer_freeze_constraints(
                geom_cfg,
                calc_cfg,
                layer_info,
                echo_fn=click.echo,
            )

            # Distance-based overrides for Hessian-target and movable MM selection.
            if hess_cutoff is not None:
                calc_cfg["hess_cutoff"] = hess_cutoff
            if movable_cutoff is not None:
                calc_cfg["movable_cutoff"] = movable_cutoff
                calc_cfg["use_bfactor_layers"] = False

            for key in ("input_pdb", "real_parm7", "model_pdb", "mm_fd_dir"):
                val = calc_cfg.get(key)
                if val:
                    calc_cfg[key] = str(Path(val).expanduser().resolve())
            ensure_dir(out_dir_path)

            ref_pdb_resolve = source_path.resolve()

            click.echo(pretty_block("geom", format_freeze_atoms_for_echo(geom_cfg, key="freeze_atoms")))
            echo_calc = format_freeze_atoms_for_echo(filter_calc_for_echo(calc_cfg), key="freeze_atoms")
            click.echo(pretty_block("calc", echo_calc))
            echo_opt = strip_inherited_keys({**opt_cfg, "out_dir": str(out_dir_path)}, OPT_BASE_KW, mode="same")
            click.echo(pretty_block("opt", echo_opt))
            # Show only lbfgs-specific settings, not inherited from opt_cfg
            echo_lbfgs = strip_inherited_keys(lbfgs_cfg, opt_cfg)
            click.echo(pretty_block("lbfgs", echo_lbfgs))
            click.echo(pretty_block("bias", bias_cfg))

            pdb_atom_meta: List[Dict[str, Any]] = []
            if source_path.suffix.lower() == ".pdb":
                pdb_atom_meta = load_pdb_atom_metadata(source_path)

            if scan_list_raw is None:
                raise click.BadParameter("--scan-lists is required.")
            scan_one_based = bool(one_based)
            scan_source = "--scan-lists"
            if is_scan_spec_file(scan_list_raw):
                spec_path = Path(scan_list_raw)
                parsed, raw_pairs, scan_one_based = parse_scan_spec_quads(
                    spec_path,
                    expected_len=2,
                    one_based_default=one_based,
                    atom_meta=pdb_atom_meta,
                    option_name="--scan-lists",
                )
                scan_source = f"--scan-lists ({spec_path})"
            else:
                parsed, raw_pairs = parse_scan_list_quads(
                    scan_list_raw,
                    expected_len=2,
                    one_based=scan_one_based,
                    atom_meta=pdb_atom_meta,
                    option_name="--scan-lists",
                )
            (i1, j1, low1, high1), (i2, j2, low2, high2) = parsed
            d1_label_csv = axis_label_csv("d1", i1, j1, scan_one_based, pdb_atom_meta, raw_pairs[0])
            d2_label_csv = axis_label_csv("d2", i2, j2, scan_one_based, pdb_atom_meta, raw_pairs[1])
            d1_label_html = axis_label_html(d1_label_csv)
            d2_label_html = axis_label_html(d2_label_csv)
            if print_parsed:
                click.echo(
                    pretty_block(
                        "scan-parsed",
                        {
                            "source": scan_source,
                            "one_based": bool(scan_one_based),
                            "pairs_0based": parsed,
                        },
                    )
                )
            click.echo(
                pretty_block(
                    "scan-list",
                    {"d1": (i1 + 1, j1 + 1, low1, high1), "d2": (i2 + 1, j2 + 1, low2, high2)},
                )
            )
            if pdb_atom_meta:
                click.echo("[scan2d] PDB atom details for scanned pairs:")
                legend = PDB_ATOM_META_HEADER
                click.echo(f"        legend: {legend}")
                click.echo(f"  d1 i: {format_pdb_atom_metadata(pdb_atom_meta, i1)}")
                click.echo(f"     j: {format_pdb_atom_metadata(pdb_atom_meta, j1)}")
                click.echo(f"  d2 i: {format_pdb_atom_metadata(pdb_atom_meta, i2)}")
                click.echo(f"     j: {format_pdb_atom_metadata(pdb_atom_meta, j2)}")

            # Directory layout: final outputs under out_dir/, optimizer scratch in a temporary directory
            tmp_root = Path(tempfile.mkdtemp(prefix="scan2d_tmp_"))
            grid_dir = out_dir_path / "grid"
            tmp_opt_dir = tmp_root / "opt"
            ensure_dir(grid_dir)
            ensure_dir(tmp_opt_dir)
            final_dir = out_dir_path

            coord_type = geom_cfg.get("coord_type", "cart")
            geom_outer = geom_loader(geom_input_path, coord_type=coord_type)
            freeze = list(geom_cfg.get("freeze_atoms") or [])
            if freeze:
                try:
                    geom_outer.freeze_atoms = np.array(freeze, dtype=int)
                except Exception:
                    logger.debug("Failed to set freeze_atoms on geometry", exc_info=True)

            base_calc = mlmm(**calc_cfg)
            biased = HarmonicBiasCalculator(base_calc, k=float(bias_cfg["k"]))

            try:
                import torch as _torch
                _resolved_dev = "cuda" if _torch.cuda.is_available() else "cpu"
                click.echo(f"[calc] Resolved device: {_resolved_dev}")
            except Exception:
                pass

            if preopt:
                click.echo("[preopt] Unbiased relaxation of the initial structure ...")
                geom_outer.set_calculator(base_calc)
                optimizer0 = _make_lbfgs(
                    geom_outer,
                    lbfgs_cfg,
                    opt_cfg,
                    max_step_bohr=float(max_step_size) * ANG2BOHR,
                    relax_max_cycles=relax_max_cycles,
                    out_dir=tmp_opt_dir,
                    prefix="preopt",
                )
                try:
                    optimizer0.run()
                except ZeroStepLength:
                    click.echo("[preopt] ZeroStepLength — continuing.", err=True)
                except OptimizationError as exc:
                    click.echo(f"[preopt] OptimizationError — {exc}", err=True)

            records: List[Dict[str, Any]] = []
            # Keep track of previously visited structures (preopt + biased scans)
            # so that we can always start each new scan from the closest one in
            # terms of the scanned distances.
            grid_states: List[Dict[str, Any]] = []

            # Measure reference distances on the (pre)optimized structure
            d1_ref = distance_A_from_coords(np.asarray(geom_outer.coords).reshape(-1, 3), i1, j1)
            d2_ref = distance_A_from_coords(np.asarray(geom_outer.coords).reshape(-1, 3), i2, j2)
            if math.isfinite(d1_ref) and math.isfinite(d2_ref):
                click.echo(
                    f"[center] reference distances from (pre)optimized structure: "
                    f"d1 = {d1_ref:.3f} Å, d2 = {d2_ref:.3f} Å"
                )

                # Write preoptimized structure into the grid directory with distance-based name
                d1_ref_tag = distance_tag(d1_ref)
                d2_ref_tag = distance_tag(d2_ref)
                preopt_xyz_path = grid_dir / f"preopt_i{d1_ref_tag}_j{d2_ref_tag}.xyz"
                try:
                    xyz_pre = geom_outer.as_xyz()
                    if not xyz_pre.endswith("\n"):
                        xyz_pre += "\n"
                    with open(preopt_xyz_path, "w") as handle:
                        handle.write(xyz_pre)

                    convert_and_annotate_xyz_to_pdb(
                        preopt_xyz_path,
                        ref_pdb_resolve,
                        preopt_xyz_path.with_suffix(".pdb"),
                        model_pdb_path,
                        freeze_atoms_final,
                    )
                except Exception as exc:
                    click.echo(
                        f"[write] WARNING: failed to write or convert {preopt_xyz_path.name}: {exc}",
                        err=True,
                    )

                # Unbiased energy of the (pre)optimized structure for inclusion in the PES
                preopt_energy_h = unbiased_energy_hartree(geom_outer, base_calc)
                records.append(
                    {
                        "i": -1,
                        "j": -1,
                        "d1_A": float(d1_ref),
                        "d2_A": float(d2_ref),
                        "energy_hartree": preopt_energy_h,
                        "bias_converged": True,
                        "is_preopt": True,
                    }
                )
                # Also store a snapshot of this structure as the first candidate
                # starting point for subsequent biased scans.
                grid_states.append(
                    {
                        "d1_A": float(d1_ref),
                        "d2_A": float(d2_ref),
                        "geom": _snapshot_geometry(geom_outer),
                    }
                )
            else:
                click.echo(
                    "[center] WARNING: failed to determine reference distances; using grid order as-is.",
                    err=True,
                )
                d1_ref_tag = None
                d2_ref_tag = None

            # Build distance grids and reorder so that scanning starts near the reference structure
            d1_values = values_from_bounds(low1, high1, float(max_step_size))
            d2_values = values_from_bounds(low2, high2, float(max_step_size))

            if math.isfinite(d1_ref):
                d1_values = np.array(
                    sorted(d1_values, key=lambda v: abs(v - d1_ref)),
                    dtype=float,
                )
            if math.isfinite(d2_ref):
                d2_values = np.array(
                    sorted(d2_values, key=lambda v: abs(v - d2_ref)),
                    dtype=float,
                )

            N1, N2 = len(d1_values), len(d2_values)
            click.echo(f"[grid] d1 steps = {N1}  values(A)={list(map(lambda x: f'{x:.3f}', d1_values))}")
            click.echo(f"[grid] d2 steps = {N2}  values(A)={list(map(lambda x: f'{x:.3f}', d2_values))}")
            click.echo(f"[grid] total grid points = {N1 * N2}")

            max_step_bohr = float(max_step_size) * ANG2BOHR

            for i_idx, d1_target in enumerate(d1_values):
                d1_tag = distance_tag(d1_target)
                click.echo(f"\n--- d1 step {i_idx + 1}/{N1} : target = {d1_target:.3f} Å ---")

                # Choose the closest previously visited structure (in d1) as the
                # starting point for the d1-biased relaxation.
                start_outer = _select_closest_state_1d(grid_states, float(d1_target))
                if start_outer is None:
                    start_outer = geom_outer
                geom_outer = _snapshot_geometry(start_outer)

                geom_outer.set_calculator(biased)
                biased.set_pairs([(i1, j1, float(d1_target))])
                geom_outer.set_calculator(biased)

                opt1 = _make_lbfgs(
                    geom_outer,
                    lbfgs_cfg,
                    opt_cfg,
                    max_step_bohr=max_step_bohr,
                    relax_max_cycles=relax_max_cycles,
                    out_dir=tmp_opt_dir,
                    prefix=f"d1_{d1_tag}",
                )
                try:
                    opt1.run()
                except ZeroStepLength:
                    click.echo(f"[d1 {i_idx}] ZeroStepLength — continuing to d2 scan.", err=True)
                except OptimizationError as exc:
                    click.echo(f"[d1 {i_idx}] OptimizationError — {exc}", err=True)

                # Record the relaxed (d1-biased) structure as another candidate
                # starting point for subsequent grid points.
                d1_cur_outer = distance_A_from_coords(np.asarray(geom_outer.coords).reshape(-1, 3), i1, j1)
                d2_cur_outer = distance_A_from_coords(np.asarray(geom_outer.coords).reshape(-1, 3), i2, j2)
                if math.isfinite(d1_cur_outer) and math.isfinite(d2_cur_outer):
                    grid_states.append(
                        {
                            "d1_A": float(d1_cur_outer),
                            "d2_A": float(d2_cur_outer),
                            "geom": _snapshot_geometry(geom_outer),
                        }
                    )

                trj_blocks = [] if dump else None

                for j_idx, d2_target in enumerate(d2_values):
                    d2_tag = distance_tag(d2_target)

                    # For each (d1, d2) grid point, choose as initial structure the
                    # previously visited geometry whose scanned distances are closest
                    # to the current targets.
                    start_inner = _select_closest_state(grid_states, float(d1_target), float(d2_target))
                    if start_inner is None:
                        # Fallback: use the d1-relaxed structure for this slice.
                        start_inner = geom_outer
                    geom_inner = _snapshot_geometry(start_inner)
                    geom_inner.set_calculator(biased)

                    biased.set_pairs([(i1, j1, float(d1_target)), (i2, j2, float(d2_target))])

                    opt2 = _make_lbfgs(
                        geom_inner,
                        lbfgs_cfg,
                        opt_cfg,
                        max_step_bohr=max_step_bohr,
                        relax_max_cycles=relax_max_cycles,
                        out_dir=tmp_opt_dir,
                        prefix=f"d1_{d1_tag}_d2_{d2_tag}",
                    )
                    try:
                        opt2.run()
                        converged = True
                    except ZeroStepLength:
                        click.echo(
                            f"[d1 {i_idx}, d2 {j_idx}] ZeroStepLength — recorded anyway.",
                            err=True,
                        )
                        converged = False
                    except OptimizationError as exc:
                        click.echo(f"[d1 {i_idx}, d2 {j_idx}] OptimizationError — {exc}", err=True)
                        converged = False

                    energy_h = unbiased_energy_hartree(geom_inner, base_calc)

                    # Record this grid point as a new candidate starting structure
                    # for subsequent scans.
                    d1_cur = distance_A_from_coords(np.asarray(geom_inner.coords).reshape(-1, 3), i1, j1)
                    d2_cur = distance_A_from_coords(np.asarray(geom_inner.coords).reshape(-1, 3), i2, j2)
                    if math.isfinite(d1_cur) and math.isfinite(d2_cur):
                        grid_states.append(
                            {
                                "d1_A": float(d1_cur),
                                "d2_A": float(d2_cur),
                                "geom": _snapshot_geometry(geom_inner),
                            }
                        )

                    # Distance-based filenames: e.g., point_i125_j324.xyz for d1=1.25 Å, d2=3.24 Å
                    xyz_path = grid_dir / f"point_i{d1_tag}_j{d2_tag}.xyz"
                    try:
                        xyz = geom_inner.as_xyz()
                        if not xyz.endswith("\n"):
                            xyz += "\n"
                        with open(xyz_path, "w") as handle:
                            handle.write(xyz)

                        # Convert grid-point XYZ to PDB (with B-factor annotation)
                        convert_and_annotate_xyz_to_pdb(
                            xyz_path,
                            ref_pdb_resolve,
                            xyz_path.with_suffix(".pdb"),
                            model_pdb_path,
                            freeze_atoms_final,
                        )
                    except Exception as exc:
                        click.echo(
                            f"[write] WARNING: failed to write or convert {xyz_path.name}: {exc}",
                            err=True,
                        )

                    if dump and trj_blocks is not None:
                        block = geom_inner.as_xyz()
                        if not block.endswith("\n"):
                            block += "\n"
                        trj_blocks.append(block)

                    records.append(
                        {
                            "i": int(i_idx),
                            "j": int(j_idx),
                            "d1_A": float(d1_target),
                            "d2_A": float(d2_target),
                            "energy_hartree": energy_h,
                            "bias_converged": bool(converged),
                            "is_preopt": False,
                        }
                    )

                if dump and trj_blocks:
                    # Distance-based filename for inner path as well
                    trj_path = grid_dir / f"inner_path_d1_{d1_tag}_trj.xyz"
                    try:
                        with open(trj_path, "w") as handle:
                            handle.write("".join(trj_blocks))
                        click.echo(f"[write] Wrote '{trj_path}'.")

                        # Convert inner-path TRJ to multi-model PDB (with B-factor annotation)
                        convert_and_annotate_xyz_to_pdb(
                            trj_path,
                            ref_pdb_resolve,
                            trj_path.with_suffix(".pdb"),
                            model_pdb_path,
                            freeze_atoms_final,
                        )
                    except Exception as exc:
                        click.echo(
                            f"[write] WARNING: failed to write or convert '{trj_path}': {exc}",
                            err=True,
                        )

            df = pd.DataFrame.from_records(records)
            if df.empty:
                click.echo("No grid records produced; aborting.", err=True)
                sys.exit(1)

            if baseline == "first":
                mask = (df["i"] == 0) & (df["j"] == 0)
                if mask.sum() == 0:
                    click.echo("WARNING: baseline='first' but grid point (0,0) not found; falling back to min.", err=True)
                    ref = float(df["energy_hartree"].min())
                else:
                    ref = float(df.loc[mask, "energy_hartree"].iloc[0])
            else:
                ref = float(df["energy_hartree"].min())
            df["energy_kcal"] = (df["energy_hartree"] - ref) * AU2KCALPERMOL
            df["d1_label"] = d1_label_csv
            df["d2_label"] = d2_label_csv

            surface_csv = final_dir / "surface.csv"
            df.to_csv(surface_csv, index=False)
            click.echo(f"[write] Wrote '{surface_csv}'.")

            # ===== Plots (RBF on a fixed 50×50 grid, unified layout, placed under final_dir) =====
            d1_points = df["d1_A"].to_numpy(dtype=float)
            d2_points = df["d2_A"].to_numpy(dtype=float)
            z_points = df["energy_kcal"].to_numpy(dtype=float)
            mask = np.isfinite(d1_points) & np.isfinite(d2_points) & np.isfinite(z_points)
            if not np.any(mask):
                click.echo("[plot] No finite data for plotting.", err=True)
                sys.exit(1)

            x_min, x_max = float(np.min(d1_points[mask])), float(np.max(d1_points[mask]))
            y_min, y_max = float(np.min(d2_points[mask])), float(np.max(d2_points[mask]))

            xi = np.linspace(x_min, x_max, 50)
            yi = np.linspace(y_min, y_max, 50)
            XI, YI = np.meshgrid(xi, yi)

            rbf = Rbf(d1_points[mask], d2_points[mask], z_points[mask], function="multiquadric")
            ZI = rbf(XI, YI)

            vmin = float(np.nanmin(ZI)) if zmin is None else float(zmin)
            vmax = float(np.nanmax(ZI)) if zmax is None else float(zmax)
            if not np.isfinite(vmin) or not np.isfinite(vmax) or vmax <= vmin:
                vmin, vmax = float(np.nanmin(ZI)), float(np.nanmax(ZI))

            # Choose neat contour/tick steps
            def _nice_step(span: float) -> float:
                if span <= 0:
                    return 1.0
                raw = span / 6.0
                mag = 10 ** math.floor(math.log10(raw))
                candidates = (0.5, 1, 2, 5, 10, 20)
                # start with the first candidate
                best = candidates[0] * mag
                best_err = abs(best - raw)
                for m in candidates[1:]:
                    s = m * mag
                    err = abs(s - raw)
                    if err < best_err:
                        best, best_err = s, err
                return best

            c_step = _nice_step(vmax - vmin)
            c_start = math.floor(vmin / c_step) * c_step
            c_end = math.ceil(vmax / c_step) * c_step

            # ---- 2D contour plot (PNG with explicit size) ----
            fig2d = go.Figure(
                data=go.Contour(
                    z=ZI,
                    x=xi,
                    y=yi,
                    contours=dict(start=c_start, end=c_end, size=c_step),
                    zmin=vmin,
                    zmax=vmax,
                    contours_coloring="heatmap",
                    colorscale="plasma",
                    colorbar=dict(
                        title=dict(text="(kcal/mol)", side="top", font=dict(size=16, color="#1C1C1C")),
                        tickfont=dict(size=14, color="#1C1C1C"),
                        ticks="inside",
                        ticklen=10,
                        tickcolor="#1C1C1C",
                        outlinecolor="#1C1C1C",
                        outlinewidth=2,
                        lenmode="fraction",
                        len=1.11,
                        x=1.05,
                        y=0.53,
                        xanchor="left",
                        yanchor="middle",
                    ),
                )
            )
            fig2d.update_layout(
                width=640,
                height=600,
                xaxis_title=d1_label_html,
                yaxis_title=d2_label_html,
                plot_bgcolor="white",
                xaxis=dict(
                    range=[x_min, x_max],
                    showline=True,
                    linewidth=3,
                    linecolor="#1C1C1C",
                    mirror=True,
                    tickson="boundaries",
                    ticks="inside",
                    tickwidth=3,
                    tickcolor="#1C1C1C",
                    title_font=dict(size=18, color="#1C1C1C"),
                    tickfont=dict(size=18, color="#1C1C1C"),
                    tickvals=list(np.linspace(x_min, x_max, 6)),
                    tickformat=".2f",
                ),
                yaxis=dict(
                    range=[y_min, y_max],
                    showline=True,
                    linewidth=3,
                    linecolor="#1C1C1C",
                    mirror=True,
                    tickson="boundaries",
                    ticks="inside",
                    tickwidth=3,
                    tickcolor="#1C1C1C",
                    title_font=dict(size=18, color="#1C1C1C"),
                    tickfont=dict(size=18, color="#1C1C1C"),
                    tickvals=list(np.linspace(y_min, y_max, 6)),
                    tickformat=".2f",
                ),
                margin=dict(l=10, r=10, b=10, t=40),
            )
            png2d = final_dir / "scan2d_map.png"
            fig2d.write_image(str(png2d), scale=2, engine="kaleido", width=680, height=600)
            click.echo(f"[plot] Wrote '{png2d}'.")

            # ---- 3D surface plus base-plane projection ----
            spread = vmax - vmin if (vmax > vmin) else 1.0
            z_bottom = vmin - spread
            z_top = vmax

            # Avoid ticks below zmin (= vmin) and snap to sensible values
            z_step = _nice_step(vmax - vmin)
            z_start_tick = math.ceil(vmin / z_step) * z_step  # First tick must be ≥ vmin
            z_ticks = np.arange(z_start_tick, z_top + 0.5 * z_step, z_step).tolist()

            surface3d = go.Surface(
                x=XI,
                y=YI,
                z=ZI,
                colorscale="plasma",
                cmin=vmin,
                cmax=vmax,
                colorbar=dict(
                    title=dict(text="(kcal/mol)", side="top", font=dict(size=16, color="#1C1C1C")),
                    tickfont=dict(size=14, color="#1C1C1C"),
                    ticks="inside",
                    ticklen=10,
                    tickcolor="#1C1C1C",
                    outlinecolor="#1C1C1C",
                    outlinewidth=2,
                    lenmode="fraction",
                    len=1.11,
                    x=1.05,
                    y=0.53,
                    xanchor="left",
                    yanchor="middle",
                ),
                contours={
                    "z": {
                        "show": True,
                        "start": c_start,
                        "end": c_end,
                        "size": c_step,
                        "color": "black",
                        "project": {"z": True},
                    }
                },
                name="3D Surface",
            )

            plane_proj = go.Surface(
                x=XI,
                y=YI,
                z=np.full_like(ZI, z_bottom),
                surfacecolor=ZI,
                colorscale="plasma",
                cmin=vmin,
                cmax=vmax,
                showscale=False,
                opacity=1.0,
                name="2D Contour Projection (Bottom)",
            )

            fig3d = go.Figure(data=[surface3d, plane_proj])
            fig3d.update_layout(
                title="Energy Landscape with 2D PES Scan",
                width=800,
                height=700,
                scene=dict(
                    bgcolor="rgba(0,0,0,0)",
                    xaxis=dict(
                        title=d1_label_html,
                        range=[x_min, x_max],
                        showline=True,
                        linewidth=4,
                        linecolor="#1C1C1C",
                        mirror=True,
                        ticks="inside",
                        tickwidth=4,
                        tickcolor="#1C1C1C",
                        gridcolor="rgba(0,0,0,0.1)",
                        zerolinecolor="rgba(0,0,0,0.1)",
                        showbackground=False,
                    ),
                    yaxis=dict(
                        title=d2_label_html,
                        range=[y_min, y_max],
                        showline=True,
                        linewidth=4,
                        linecolor="#1C1C1C",
                        mirror=True,
                        ticks="inside",
                        tickwidth=4,
                        tickcolor="#1C1C1C",
                        gridcolor="rgba(0,0,0,0.1)",
                        zerolinecolor="rgba(0,0,0,0.1)",
                        showbackground=False,
                    ),
                    zaxis=dict(
                        title="Potential Energy (kcal/mol)",
                        range=[z_bottom, z_top],
                        tickmode="array",
                        tickvals=z_ticks,
                        showline=True,
                        linewidth=4,
                        linecolor="#1C1C1C",
                        mirror=True,
                        ticks="inside",
                        tickwidth=4,
                        tickcolor="#1C1C1C",
                        showgrid=True,
                        gridcolor="rgba(0,0,0,0.1)",
                        zerolinecolor="rgba(0,0,0,0.1)",
                        showbackground=False,
                    ),
                ),
                margin=dict(l=10, r=20, b=10, t=40),
                paper_bgcolor="white",
            )

            html3d = final_dir / "scan2d_landscape.html"
            fig3d.write_html(str(html3d))
            click.echo(f"[plot] Wrote '{html3d}'.")

            click.echo("=== 2D Scan finished ===\n")
            click.echo(format_elapsed("[time] Elapsed Time for 2D Scan", time_start))

            if out_json:
                from .utils import write_result_json
                min_energy = float(df["energy_hartree"].min()) if not df.empty else None
                result_data: Dict[str, Any] = {
                    "status": "completed",
                    "n_grid_points": len(df),
                    "pair1": {"i": int(i1 + 1), "j": int(j1 + 1), "low": float(low1), "high": float(high1)},
                    "pair2": {"i": int(i2 + 1), "j": int(j2 + 1), "low": float(low2), "high": float(high2)},
                    "backend": calc_cfg.get("backend", "uma"),
                    "charge": calc_cfg.get("model_charge"),
                    "spin": calc_cfg.get("model_mult"),
                    "min_energy_hartree": min_energy,
                    "files": {
                        "surface_csv": "surface.csv",
                        "scan2d_map_png": "scan2d_map.png",
                        "scan2d_landscape_html": "scan2d_landscape.html",
                    },
                }
                write_result_json(
                    final_dir, result_data,
                    command="scan2d",
                    elapsed_seconds=time.perf_counter() - time_start,
                )

    except KeyboardInterrupt:
        click.echo("\nInterrupted by user.", err=True)
        sys.exit(130)
    except Exception as exc:
        tb = "".join(traceback.format_exception(type(exc), exc, exc.__traceback__))
        click.echo("Unhandled exception during 2D scan:\n" + textwrap.indent(tb, "  "), err=True)
        sys.exit(1)
    finally:
        if tmp_root is not None:
            shutil.rmtree(tmp_root, ignore_errors=True)
        # Release GPU memory so subsequent pipeline stages don't OOM
        base_calc = geom_outer = geom_inner = optimizer0 = None
        gc.collect()  # break cyclic refs inside torch.nn.Module
        if torch.cuda.is_available():
            torch.cuda.empty_cache()
