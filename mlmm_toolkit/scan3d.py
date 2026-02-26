"""
ML/MM three-distance scan with harmonic restraints (d1, d2, d3 grid).

Example:
    mlmm scan3d -i input.pdb --real-parm7 real.parm7 --model-pdb ml_region.pdb -q 0 \
        --scan-lists "[(12,45,1.30,3.10),(10,55,1.20,3.20),(15,60,1.10,3.00)]"

For detailed documentation, see: docs/scan3d.md
"""

from __future__ import annotations

import functools
from copy import deepcopy
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence, Tuple, Set

import ast
import math
import sys
import textwrap
import traceback
import tempfile
import os
import time

import click
import numpy as np
import pandas as pd
from scipy.interpolate import Rbf
import plotly.graph_objects as go

from pysisyphus.helpers import geom_loader
from pysisyphus.optimizers.LBFGS import LBFGS
from pysisyphus.optimizers.exceptions import OptimizationError, ZeroStepLength
from pysisyphus.constants import ANG2BOHR, AU2KCALPERMOL

from .mlmm_calc import mlmm
from .opt import (
    GEOM_KW as _OPT_GEOM_KW,
    CALC_KW as _OPT_CALC_KW,
    OPT_BASE_KW as _OPT_BASE_KW,
    LBFGS_KW as _OPT_LBFGS_KW,
    HarmonicBiasCalculator,
    _parse_freeze_atoms,
    _normalize_geom_freeze,
)
from .utils import (
    apply_ref_pdb_override,
    apply_layer_freeze_constraints,
    deep_update,
    load_yaml_dict,
    apply_yaml_overrides,
    pretty_block,
    strip_inherited_keys,
    format_freeze_atoms_for_echo,
    format_elapsed,
    merge_freeze_atom_indices,
    prepare_input_structure,
    resolve_charge_spin_or_raise,
    convert_xyz_to_pdb,
    load_pdb_atom_metadata,
    parse_scan_list_quads,
    parse_scan_spec_quads,
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
from .cli_utils import resolve_yaml_sources, load_merged_yaml_cfg

# Shared defaults (copied from opt.py to keep ML/MM behaviour consistent)
GEOM_KW: Dict[str, Any] = deepcopy(_OPT_GEOM_KW)
CALC_KW: Dict[str, Any] = deepcopy(_OPT_CALC_KW)
OPT_BASE_KW: Dict[str, Any] = deepcopy(_OPT_BASE_KW)
OPT_BASE_KW.update(
    {
        "out_dir": "./result_scan3d/",
        "dump": False,
        "max_cycles": 10000,
    }
)
LBFGS_KW: Dict[str, Any] = deepcopy(_OPT_LBFGS_KW)
LBFGS_KW.update({"out_dir": "./result_scan3d/"})
BIAS_KW: Dict[str, Any] = {"k": 100.0}  # eV/Å^2

HARTREE_TO_KCAL_MOL = 627.50961
_VOLUME_GRID_N = 50  # 50×50×50 RBF interpolation grid

_snapshot_geometry = functools.partial(snapshot_geometry, coord_type_default="cart")


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
    help="3D distance scan with harmonic restraints using the ML/MM calculator.",
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
    "--real-parm7",
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
    help="Detect ML/MM layers from input PDB B-factors (B=0/10/20). "
         "If disabled, you must provide --model-pdb or --model-indices.",
)
@click.option("-q", "--charge", type=int, required=True, help="ML-region total charge.")
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
    "--scan-lists",
    "scan_list_raw",
    type=str,
    required=False,
    help='Python-like list with three quadruples: "[(i1,j1,low1,high1),(i2,j2,low2,high2),(i3,j3,low3,high3)]".',
)
@click.option(
    "--spec",
    "spec_path",
    type=click.Path(path_type=Path, exists=True, dir_okay=False),
    required=False,
    help="YAML/JSON scan spec file (recommended). Use this instead of --scan-lists.",
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
    help="Print parsed scan targets after resolving --spec/--scan-lists.",
)
@click.option(
    "--max-step-size",
    type=float,
    default=0.20,
    show_default=True,
    help="Maximum spacing between successive distance targets [Å].",
)
@click.option("--bias-k", type=float, default=100.0, show_default=True, help="Harmonic well strength k [eV/Å^2].")
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
    help="Write inner d3 scan TRJs per (d1,d2) slice.",
)
@click.option(
    "--out-dir",
    type=str,
    default="./result_scan3d/",
    show_default=True,
    help="Base output directory.",
)
@click.option(
    "--thresh",
    type=str,
    default=None,
    show_default=False,
    help="Convergence preset (gau_loose|gau|gau_tight|gau_vtight|baker|never).",
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
    default=True,
    show_default=True,
    help="Run an unbiased pre-optimization.",
)
@click.option(
    "--baseline",
    type=click.Choice(["min", "first"]),
    default="min",
    show_default=True,
    help="Reference for relative energy (kcal/mol): 'min' or 'first' (i=0,j=0,k=0).",
)
@click.option(
    "--zmin",
    type=float,
    default=None,
    show_default=False,
    help="Lower bound of the color scale (kcal/mol).",
)
@click.option(
    "--zmax",
    type=float,
    default=None,
    show_default=False,
    help="Upper bound of the color scale (kcal/mol).",
)
def cli(
    input_path: Path,
    real_parm7: Path,
    model_pdb: Optional[Path],
    model_indices_str: Optional[str],
    model_indices_one_based: bool,
    detect_layer: bool,
    charge: Optional[int],
    spin: Optional[int],
    freeze_atoms_cli: Optional[str],
    hess_cutoff: Optional[float],
    movable_cutoff: Optional[float],
    scan_list_raw: Optional[str],
    spec_path: Optional[Path],
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
) -> None:
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

    try:
        with prepare_input_structure(input_path) as prepared_input:
            try:
                apply_ref_pdb_override(prepared_input, ref_pdb)
            except click.BadParameter as e:
                click.echo(f"ERROR: {e}", err=True)
                sys.exit(1)
            geom_input_path = prepared_input.geom_path
            source_path = prepared_input.source_path
            charge, spin = resolve_charge_spin_or_raise(prepared_input, charge, spin)

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
            echo_calc = strip_inherited_keys(calc_cfg, CALC_KW, mode="same")
            echo_calc = format_freeze_atoms_for_echo(echo_calc, key="freeze_atoms")
            click.echo(pretty_block("calc", echo_calc))
            echo_opt = strip_inherited_keys({**opt_cfg, "out_dir": str(out_dir_path)}, OPT_BASE_KW, mode="same")
            click.echo(pretty_block("opt", echo_opt))
            echo_lbfgs = strip_inherited_keys(lbfgs_cfg, opt_cfg)
            click.echo(pretty_block("lbfgs", echo_lbfgs))
            click.echo(pretty_block("bias", bias_cfg))

            pdb_atom_meta: List[Dict[str, Any]] = []
            if source_path.suffix.lower() == ".pdb":
                pdb_atom_meta = load_pdb_atom_metadata(source_path)

            if spec_path is not None and scan_list_raw is not None:
                raise click.BadParameter("Use either --spec or --scan-lists, not both.")
            scan_one_based = bool(one_based)
            scan_source = "--scan-lists"
            if spec_path is not None:
                parsed, raw_pairs, scan_one_based = parse_scan_spec_quads(
                    spec_path,
                    expected_len=3,
                    one_based_default=one_based,
                    atom_meta=pdb_atom_meta,
                    option_name="--spec",
                )
                scan_source = f"--spec ({spec_path})"
            else:
                if scan_list_raw is None:
                    raise click.BadParameter("Provide either --spec or --scan-lists.")
                parsed, raw_pairs = parse_scan_list_quads(
                    scan_list_raw,
                    expected_len=3,
                    one_based=scan_one_based,
                    atom_meta=pdb_atom_meta,
                    option_name="--scan-lists",
                )
            (i1, j1, low1, high1), (i2, j2, low2, high2), (i3, j3, low3, high3) = parsed
            d1_label_csv = axis_label_csv("d1", i1, j1, scan_one_based, pdb_atom_meta, raw_pairs[0])
            d2_label_csv = axis_label_csv("d2", i2, j2, scan_one_based, pdb_atom_meta, raw_pairs[1])
            d3_label_csv = axis_label_csv("d3", i3, j3, scan_one_based, pdb_atom_meta, raw_pairs[2])
            d1_label_html = axis_label_html(d1_label_csv)
            d2_label_html = axis_label_html(d2_label_csv)
            d3_label_html = axis_label_html(d3_label_csv)
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
                    "scan-list (0-based)",
                    {
                        "d1": (i1, j1, low1, high1),
                        "d2": (i2, j2, low2, high2),
                        "d3": (i3, j3, low3, high3),
                    },
                )
            )
            if pdb_atom_meta:
                click.echo("[scan3d] PDB atom details for scanned pairs:")
                legend = PDB_ATOM_META_HEADER
                click.echo(f"        legend: {legend}")
                click.echo(f"  d1 i: {format_pdb_atom_metadata(pdb_atom_meta, i1)}")
                click.echo(f"     j: {format_pdb_atom_metadata(pdb_atom_meta, j1)}")
                click.echo(f"  d2 i: {format_pdb_atom_metadata(pdb_atom_meta, i2)}")
                click.echo(f"     j: {format_pdb_atom_metadata(pdb_atom_meta, j2)}")
                click.echo(f"  d3 i: {format_pdb_atom_metadata(pdb_atom_meta, i3)}")
                click.echo(f"     j: {format_pdb_atom_metadata(pdb_atom_meta, j3)}")

            # Directory layout
            tmp_root = Path(tempfile.mkdtemp(prefix="scan3d_tmp_"))
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
                    pass

            base_calc = mlmm(**calc_cfg)
            biased = HarmonicBiasCalculator(base_calc, k=float(bias_cfg["k"]))

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
                    prefix="preopt_",
                )
                try:
                    optimizer0.run()
                except ZeroStepLength:
                    click.echo("[preopt] ZeroStepLength — continuing.", err=True)
                except OptimizationError as exc:
                    click.echo(f"[preopt] OptimizationError — {exc}", err=True)

            records: List[Dict[str, Any]] = []

            # Measure reference distances on the (pre)optimized structure
            coords_outer = np.asarray(geom_outer.coords).reshape(-1, 3)
            d1_ref = distance_A_from_coords(coords_outer, i1, j1)
            d2_ref = distance_A_from_coords(coords_outer, i2, j2)
            d3_ref = distance_A_from_coords(coords_outer, i3, j3)

            if math.isfinite(d1_ref) and math.isfinite(d2_ref) and math.isfinite(d3_ref):
                click.echo(
                    f"[center] reference distances from (pre)optimized structure: "
                    f"d1 = {d1_ref:.3f} Å, d2 = {d2_ref:.3f} Å, d3 = {d3_ref:.3f} Å"
                )

                # Write preoptimized structure
                d1_ref_tag = distance_tag(d1_ref)
                d2_ref_tag = distance_tag(d2_ref)
                d3_ref_tag = distance_tag(d3_ref)
                preopt_xyz_path = grid_dir / f"preopt_i{d1_ref_tag}_j{d2_ref_tag}_k{d3_ref_tag}.xyz"
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

                preopt_energy_h = unbiased_energy_hartree(geom_outer, base_calc)
                records.append(
                    {
                        "i": -1,
                        "j": -1,
                        "k": -1,
                        "d1_A": float(d1_ref),
                        "d2_A": float(d2_ref),
                        "d3_A": float(d3_ref),
                        "energy_hartree": preopt_energy_h,
                        "bias_converged": True,
                        "is_preopt": True,
                    }
                )
            else:
                click.echo(
                    "[center] WARNING: failed to determine reference distances; using grid order as-is.",
                    err=True,
                )

            # Build distance grids and reorder so that scanning starts near the reference structure
            d1_values = values_from_bounds(low1, high1, float(max_step_size))
            d2_values = values_from_bounds(low2, high2, float(max_step_size))
            d3_values = values_from_bounds(low3, high3, float(max_step_size))

            if math.isfinite(d1_ref):
                d1_values = np.array(sorted(d1_values, key=lambda v: abs(v - d1_ref)), dtype=float)
            if math.isfinite(d2_ref):
                d2_values = np.array(sorted(d2_values, key=lambda v: abs(v - d2_ref)), dtype=float)
            if math.isfinite(d3_ref):
                d3_values = np.array(sorted(d3_values, key=lambda v: abs(v - d3_ref)), dtype=float)

            N1, N2, N3 = len(d1_values), len(d2_values), len(d3_values)
            click.echo(f"[grid] d1 steps = {N1}  values(A)={list(map(lambda x: f'{x:.3f}', d1_values))}")
            click.echo(f"[grid] d2 steps = {N2}  values(A)={list(map(lambda x: f'{x:.3f}', d2_values))}")
            click.echo(f"[grid] d3 steps = {N3}  values(A)={list(map(lambda x: f'{x:.3f}', d3_values))}")
            click.echo(f"[grid] total grid points = {N1 * N2 * N3}")

            max_step_bohr = float(max_step_size) * ANG2BOHR

            # Caches for nearest-neighbor starting geometries
            d1_geoms: Dict[int, Any] = {}
            d2_geoms: Dict[int, Dict[int, Any]] = {}
            d3_geoms: Dict[Tuple[int, int], Dict[int, Any]] = {}

            geom_outer_initial = _snapshot_geometry(geom_outer)

            # ===== 3D nested scan: d1 (outer) → d2 (middle) → d3 (inner) =====
            for i_idx, d1_target in enumerate(d1_values):
                d1_tag = distance_tag(d1_target)
                click.echo(f"\n--- d1 step {i_idx + 1}/{N1} : target = {d1_target:.3f} Å ---")

                # Choose initial geometry for this d1
                if not d1_geoms:
                    geom_outer_i = _snapshot_geometry(geom_outer_initial)
                else:
                    nearest_i = min(d1_geoms.keys(), key=lambda p: abs(d1_values[p] - d1_target))
                    geom_outer_i = _snapshot_geometry(d1_geoms[nearest_i])

                biased.set_pairs([(i1, j1, float(d1_target))])
                geom_outer_i.set_calculator(biased)

                opt1 = _make_lbfgs(
                    geom_outer_i,
                    lbfgs_cfg,
                    opt_cfg,
                    max_step_bohr=max_step_bohr,
                    relax_max_cycles=relax_max_cycles,
                    out_dir=tmp_opt_dir,
                    prefix=f"d1_{d1_tag}_",
                )
                try:
                    opt1.run()
                except ZeroStepLength:
                    click.echo(f"[d1 {i_idx}] ZeroStepLength — continuing to d2/d3 scan.", err=True)
                except OptimizationError as exc:
                    click.echo(f"[d1 {i_idx}] OptimizationError — {exc}", err=True)

                geom_after_d1 = _snapshot_geometry(geom_outer_i)
                d1_geoms[i_idx] = geom_after_d1

                if i_idx not in d2_geoms:
                    d2_geoms[i_idx] = {}

                for j_idx, d2_target in enumerate(d2_values):
                    d2_tag = distance_tag(d2_target)
                    click.echo(
                        f"  [stage] d1/d2 step ({i_idx + 1}/{N1}, {j_idx + 1}/{N2}): "
                        f"targets = ({d1_target:.3f}, {d2_target:.3f}) Å"
                    )

                    # Choose initial geometry for this (d1,d2)
                    d2_store = d2_geoms[i_idx]
                    if not d2_store:
                        geom_mid = _snapshot_geometry(geom_after_d1)
                    else:
                        nearest_j = min(d2_store.keys(), key=lambda p: abs(d2_values[p] - d2_target))
                        geom_mid = _snapshot_geometry(d2_store[nearest_j])

                    biased.set_pairs([
                        (i1, j1, float(d1_target)),
                        (i2, j2, float(d2_target)),
                    ])
                    geom_mid.set_calculator(biased)

                    opt2 = _make_lbfgs(
                        geom_mid,
                        lbfgs_cfg,
                        opt_cfg,
                        max_step_bohr=max_step_bohr,
                        relax_max_cycles=relax_max_cycles,
                        out_dir=tmp_opt_dir,
                        prefix=f"d1_{d1_tag}_d2_{d2_tag}_",
                    )
                    try:
                        opt2.run()
                    except ZeroStepLength:
                        click.echo(f"[d1 {i_idx}, d2 {j_idx}] ZeroStepLength — continuing to d3 scan.", err=True)
                    except OptimizationError as exc:
                        click.echo(f"[d1 {i_idx}, d2 {j_idx}] OptimizationError — {exc}", err=True)

                    geom_after_d2 = _snapshot_geometry(geom_mid)
                    d2_store[j_idx] = geom_after_d2

                    key_ij = (i_idx, j_idx)
                    if key_ij not in d3_geoms:
                        d3_geoms[key_ij] = {}
                    d3_store = d3_geoms[key_ij]

                    trj_blocks = [] if dump else None

                    for k_idx, d3_target in enumerate(d3_values):
                        d3_tag = distance_tag(d3_target)

                        # Choose initial geometry for this (d1,d2,d3)
                        if not d3_store:
                            geom_inner = _snapshot_geometry(geom_after_d2)
                        else:
                            nearest_k = min(d3_store.keys(), key=lambda p: abs(d3_values[p] - d3_target))
                            geom_inner = _snapshot_geometry(d3_store[nearest_k])

                        biased.set_pairs([
                            (i1, j1, float(d1_target)),
                            (i2, j2, float(d2_target)),
                            (i3, j3, float(d3_target)),
                        ])
                        geom_inner.set_calculator(biased)

                        opt3 = _make_lbfgs(
                            geom_inner,
                            lbfgs_cfg,
                            opt_cfg,
                            max_step_bohr=max_step_bohr,
                            relax_max_cycles=relax_max_cycles,
                            out_dir=tmp_opt_dir,
                            prefix=f"d1_{d1_tag}_d2_{d2_tag}_d3_{d3_tag}_",
                        )
                        try:
                            opt3.run()
                            converged = True
                        except ZeroStepLength:
                            click.echo(
                                f"    [d1 {i_idx}, d2 {j_idx}, d3 {k_idx}] ZeroStepLength — recorded anyway.",
                                err=True,
                            )
                            converged = False
                        except OptimizationError as exc:
                            click.echo(
                                f"    [d1 {i_idx}, d2 {j_idx}, d3 {k_idx}] OptimizationError — {exc}",
                                err=True,
                            )
                            converged = False

                        # Cache final geometry for nearest-neighbor reuse
                        d3_store[k_idx] = _snapshot_geometry(geom_inner)

                        energy_h = unbiased_energy_hartree(geom_inner, base_calc)

                        xyz_path = grid_dir / f"point_i{d1_tag}_j{d2_tag}_k{d3_tag}.xyz"
                        try:
                            xyz = geom_inner.as_xyz()
                            if not xyz.endswith("\n"):
                                xyz += "\n"
                            with open(xyz_path, "w") as handle:
                                handle.write(xyz)

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
                                "k": int(k_idx),
                                "d1_A": float(d1_target),
                                "d2_A": float(d2_target),
                                "d3_A": float(d3_target),
                                "energy_hartree": energy_h,
                                "bias_converged": bool(converged),
                                "is_preopt": False,
                            }
                        )

                    if dump and trj_blocks:
                        trj_path = grid_dir / f"inner_path_d1_{d1_tag}_d2_{d2_tag}_trj.xyz"
                        try:
                            with open(trj_path, "w") as handle:
                                handle.write("".join(trj_blocks))
                            click.echo(f"[write] Wrote '{trj_path}'.")

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

            # Baseline energy
            if baseline == "first":
                mask = (df["i"] == 0) & (df["j"] == 0) & (df["k"] == 0)
                if not mask.any():
                    click.echo(
                        "[baseline] 'first' requested but (i=0,j=0,k=0) missing; using global minimum instead.",
                        err=True,
                    )
                    ref_energy = float(df["energy_hartree"].min())
                else:
                    ref_energy = float(df.loc[mask, "energy_hartree"].iloc[0])
            else:
                ref_energy = float(df["energy_hartree"].min())
            df["energy_kcal"] = (df["energy_hartree"] - ref_energy) * HARTREE_TO_KCAL_MOL
            df["d1_label"] = d1_label_csv
            df["d2_label"] = d2_label_csv
            df["d3_label"] = d3_label_csv

            surface_csv = final_dir / "surface.csv"
            df.to_csv(surface_csv, index=False)
            click.echo(f"[write] Wrote '{surface_csv}'.")

            # ===== 3D RBF interpolation & visualization (isosurface) =====
            d1_points = df["d1_A"].to_numpy(dtype=float)
            d2_points = df["d2_A"].to_numpy(dtype=float)
            d3_points = df["d3_A"].to_numpy(dtype=float)
            z_points = df["energy_kcal"].to_numpy(dtype=float)

            mask = (
                np.isfinite(d1_points)
                & np.isfinite(d2_points)
                & np.isfinite(d3_points)
                & np.isfinite(z_points)
            )
            if not np.any(mask):
                click.echo("[plot] No finite data for plotting.", err=True)
                sys.exit(1)

            x_min, x_max = float(np.min(d1_points[mask])), float(np.max(d1_points[mask]))
            y_min, y_max = float(np.min(d2_points[mask])), float(np.max(d2_points[mask]))
            z_min_val, z_max_val = float(np.min(d3_points[mask])), float(np.max(d3_points[mask]))

            xi = np.linspace(x_min, x_max, _VOLUME_GRID_N)
            yi = np.linspace(y_min, y_max, _VOLUME_GRID_N)
            zi = np.linspace(z_min_val, z_max_val, _VOLUME_GRID_N)

            click.echo("[plot] 3D RBF interpolation on a 50×50×50 grid ...")
            rbf3d = Rbf(
                d1_points[mask],
                d2_points[mask],
                d3_points[mask],
                z_points[mask],
                function="multiquadric",
            )

            XI, YI, ZI = np.meshgrid(xi, yi, zi, indexing="xy")
            X_flat = XI.flatten()
            Y_flat = YI.flatten()
            Z_flat = ZI.flatten()
            E_flat = rbf3d(X_flat, Y_flat, Z_flat)

            vmin = float(np.nanmin(E_flat)) if zmin is None else float(zmin)
            vmax = float(np.nanmax(E_flat)) if zmax is None else float(zmax)
            if not np.isfinite(vmin) or not np.isfinite(vmax) or vmax <= vmin:
                vmin, vmax = float(np.nanmin(E_flat)), float(np.nanmax(E_flat))

            # Discrete isosurfaces with banded colors
            n_levels = 8
            level_values = np.linspace(vmin, vmax, n_levels + 2)[1:-1]
            level_colors = [
                "#0d0887",
                "#5b02a3",
                "#9c179e",
                "#cb4679",
                "#ed7953",
                "#fb9f3a",
                "#fdca26",
                "#f0f921",
            ]
            level_opacity = [
                1.000, 0.667, 0.444, 0.296, 0.198, 0.132, 0.088, 0.059,
            ]

            isosurfaces = []
            for lvl, color, opacity_lvl in zip(level_values, level_colors, level_opacity):
                trace = go.Isosurface(
                    x=X_flat,
                    y=Y_flat,
                    z=Z_flat,
                    value=E_flat,
                    isomin=lvl,
                    isomax=lvl,
                    surface_count=1,
                    opacity=opacity_lvl,
                    showscale=False,
                    colorscale=[[0.0, color], [1.0, color]],
                    caps=dict(x_show=False, y_show=False, z_show=False),
                    name=f"{lvl:.1f} kcal/mol",
                )
                isosurfaces.append(trace)

            # Add a dummy scatter trace to host a global colorbar
            colorbar_colorscale = [
                [idx / (len(level_colors) - 1), col]
                for idx, col in enumerate(level_colors)
            ]
            cb_tickvals = [float(v) for v in level_values]
            cb_ticktext = [f"{v:.1f}" for v in level_values]

            colorbar_trace = go.Scatter3d(
                x=[x_min],
                y=[y_min],
                z=[z_min_val],
                mode="markers",
                marker=dict(
                    size=0,
                    opacity=0.0,
                    color=[vmin, vmax],
                    colorscale=colorbar_colorscale,
                    showscale=True,
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
                        tickvals=cb_tickvals,
                        ticktext=cb_ticktext,
                    ),
                ),
                hoverinfo="none",
                showlegend=False,
            )

            # Axis labels
            d1_label = d1_label_html
            d2_label = d2_label_html
            d3_label = d3_label_html

            fig3d = go.Figure(data=isosurfaces + [colorbar_trace])
            fig3d.update_layout(
                title="3D Energy Landscape (ML/MM)",
                width=900,
                height=800,
                scene=dict(
                    bgcolor="rgba(0,0,0,0)",
                    xaxis=dict(
                        title=d1_label,
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
                        title=d2_label,
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
                        title=d3_label,
                        range=[z_min_val, z_max_val],
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
                    aspectmode="cube",
                ),
                margin=dict(l=10, r=20, b=10, t=40),
                paper_bgcolor="white",
            )

            html3d = final_dir / "scan3d_density.html"
            fig3d.write_html(str(html3d))
            click.echo(f"[plot] Wrote '{html3d}'.")

            click.echo("\n=== 3D Scan finished ===\n")
            click.echo(format_elapsed("[time] Elapsed Time for 3D Scan", time_start))

    except KeyboardInterrupt:
        click.echo("\nInterrupted by user.", err=True)
        sys.exit(130)
    except Exception as exc:
        tb = "".join(traceback.format_exception(type(exc), exc, exc.__traceback__))
        click.echo("Unhandled exception during 3D scan:\n" + textwrap.indent(tb, "  "), err=True)
        sys.exit(1)
