"""
ML/MM three-distance scan with harmonic restraints (d1, d2, d3 grid).

Example:
    mlmm scan3d -i input.pdb --parm real.parm7 --model-pdb ml_region.pdb -q 0 \
        --scan-lists "[(12,45,1.30,3.10),(10,55,1.20,3.20),(15,60,1.10,3.00)]"

For detailed documentation, see: docs/scan3d.md
"""

from __future__ import annotations

import functools
from copy import deepcopy
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

import gc
import logging
import math
import shutil
import sys
import tempfile
import time

logger = logging.getLogger(__name__)

import click
import torch
import numpy as np
import pandas as pd
from scipy.interpolate import Rbf
import plotly.graph_objects as go

from pysisyphus.helpers import geom_loader
from pysisyphus.optimizers.exceptions import OptimizationError, ZeroStepLength
from pysisyphus.constants import ANG2BOHR, AU2KCALPERMOL

from mlmm.backends.mlmm_calc import mlmm
from mlmm.core.defaults import BIAS_KW as _BIAS_KW_DEFAULT, GEOM_KW_DEFAULT, OUT_DIR_SCAN3D
from mlmm.workflows.opt import (
    GEOM_KW as _OPT_GEOM_KW,
    CALC_KW as _OPT_CALC_KW,
    OPT_BASE_KW as _OPT_BASE_KW,
    LBFGS_KW as _OPT_LBFGS_KW,
    _parse_freeze_atoms,
    _normalize_geom_freeze,
)
from mlmm.workflows.restraints import HarmonicBiasCalculator
from mlmm.workflows.opt import _convert_yaml_layer_atoms_1to0
from mlmm.core.utils import (
    apply_ref_pdb_override,
    apply_layer_freeze_constraints,
    set_convert_file_enabled,
    apply_yaml_overrides,
    pretty_block,
    strip_inherited_keys,
    filter_calc_for_echo,
    format_freeze_atoms_for_echo,
    format_elapsed,
    merge_freeze_atom_indices,
    prepare_input_structure,
    resolve_charge_spin_or_raise,
    load_pdb_atom_metadata,
    parse_scan_list_quads,
    parse_scan_spec_quads,
    is_scan_spec_file,
    axis_label_csv,
    axis_label_html,
    PDB_ATOM_META_HEADER,
    format_pdb_atom_metadata,
    parse_indices_string,
    resolve_ml_layer_assignment,
    ensure_dir,
    distance_A_from_coords,
    distance_tag,
    values_from_bounds,
    unbiased_energy_hartree,
    snapshot_geometry,
    convert_and_annotate_xyz_to_pdb,
    echo_resolved_device,
)
from mlmm.workflows.scan_common import add_scan_common_options, make_scan_lbfgs as _make_lbfgs
from mlmm.cli.common_options import (
    add_ml_layer_detection_options,
    add_coord_type_option,
    add_print_every_option,
    add_precision_option, add_backend_model_option,
    add_deterministic_option, add_allow_charge_mult_mismatch_option,
)
from mlmm.cli.decorators import resolve_yaml_sources, load_merged_yaml_cfg, make_is_param_explicit, render_cli_exception

# Shared defaults (copied from opt.py to keep ML/MM behavior consistent)
GEOM_KW: Dict[str, Any] = deepcopy(_OPT_GEOM_KW)
CALC_KW: Dict[str, Any] = deepcopy(_OPT_CALC_KW)
OPT_BASE_KW: Dict[str, Any] = deepcopy(_OPT_BASE_KW)
OPT_BASE_KW.update(
    {
        "out_dir": OUT_DIR_SCAN3D,
        "dump": False,
        "max_cycles": 10000,
    }
)
LBFGS_KW: Dict[str, Any] = deepcopy(_OPT_LBFGS_KW)
LBFGS_KW.update({"out_dir": OUT_DIR_SCAN3D})
BIAS_KW: Dict[str, Any] = deepcopy(_BIAS_KW_DEFAULT)

_VOLUME_GRID_N = 50  # 50×50×50 RBF interpolation grid

_snapshot_geometry = functools.partial(snapshot_geometry, coord_type_default="cart")


def _extract_axis_label(df: pd.DataFrame, column: str, fallback: Optional[str]) -> Optional[str]:
    if column not in df.columns:
        return fallback
    values = df[column].dropna()
    if values.empty:
        return fallback
    return str(values.iloc[0])


def _finalize_surface_and_plot(
    *,
    df: pd.DataFrame,
    final_dir: Path,
    baseline: str,
    zmin: Optional[float],
    zmax: Optional[float],
    d1_label_csv: Optional[str],
    d2_label_csv: Optional[str],
    d3_label_csv: Optional[str],
    write_surface_csv: bool,
    time_start: float,
) -> None:
    if df.empty:
        click.echo("No grid records produced; aborting.", err=True)
        sys.exit(1)

    d1_label_csv = _extract_axis_label(df, "d1_label", d1_label_csv)
    d2_label_csv = _extract_axis_label(df, "d2_label", d2_label_csv)
    d3_label_csv = _extract_axis_label(df, "d3_label", d3_label_csv)
    if d1_label_csv is None or d2_label_csv is None or d3_label_csv is None:
        click.echo(
            "[plot] WARNING: axis label metadata is missing in CSV; using generic labels.",
            err=True,
        )

    d1_label_html = axis_label_html(d1_label_csv) if d1_label_csv else "d1 (Å)"
    d2_label_html = axis_label_html(d2_label_csv) if d2_label_csv else "d2 (Å)"
    d3_label_html = axis_label_html(d3_label_csv) if d3_label_csv else "d3 (Å)"

    if "energy_kcal" not in df.columns:
        if "energy_hartree" not in df.columns:
            click.echo(
                "[baseline] energy_kcal is missing and energy_hartree is not available in CSV; aborting.",
                err=True,
            )
            sys.exit(1)

        if baseline == "first":
            ref_mask = (df["i"] == 0) & (df["j"] == 0) & (df["k"] == 0)
            if not ref_mask.any():
                click.echo(
                    "[baseline] 'first' requested but (i=0,j=0,k=0) missing; using global minimum instead.",
                    err=True,
                )
                ref_energy = float(df["energy_hartree"].min())
            else:
                ref_energy = float(df.loc[ref_mask, "energy_hartree"].iloc[0])
        else:
            ref_energy = float(df["energy_hartree"].min())
        df["energy_kcal"] = (df["energy_hartree"] - ref_energy) * AU2KCALPERMOL

    if write_surface_csv:
        surface_csv = final_dir / "surface.csv"
        df["d1_label"] = d1_label_csv
        df["d2_label"] = d2_label_csv
        df["d3_label"] = d3_label_csv
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

    fig3d = go.Figure(data=isosurfaces + [colorbar_trace])
    fig3d.update_layout(
        title="3D Energy Landscape (ML/MM)",
        width=900,
        height=800,
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
                title=d3_label_html,
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

    click.echo("\n====== 3D Scan finished ======\n", narrative=True)
    click.echo(format_elapsed("[time] Elapsed Time for 3D Scan", time_start), narrative=True)


@click.command(
    help="3D distance scan with harmonic restraints using the ML/MM calculator.",
    context_settings={"help_option_names": ["-h", "--help"]},
)
@click.option(
    "-i",
    "--input",
    "input_path",
    type=click.Path(path_type=Path, exists=True, dir_okay=False),
    required=False,
    help="Input structure file (.pdb/.xyz). Required unless --csv is provided.",
)
@click.option(
    "--parm",
    "real_parm7",
    type=click.Path(path_type=Path, exists=True, dir_okay=False),
    required=False,
    help="Amber parm7 topology for the enzyme. Required unless --csv is provided.",
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
    "--csv",
    "csv_path",
    type=click.Path(path_type=Path, exists=True, dir_okay=False),
    required=False,
    help="Plot-only mode: load a precomputed surface.csv and skip the 3D scan.",
)
@click.option(
    "--print-parsed/--no-print-parsed",
    "print_parsed",
    default=False,
    show_default=True,
    help="Print parsed scan targets after resolving --scan-lists.",
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
    help="Enable xTB point-charge embedding correction for MM→ML environmental effects (experimental).",
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
@add_scan_common_options(
    out_dir_default=OUT_DIR_SCAN3D,
    baseline_help="Reference for relative energy (kcal/mol): 'min' or 'first' (i=0,j=0,k=0).",
    dump_help="Write inner d3 scan TRJs per (d1,d2) slice.",
)
@add_ml_layer_detection_options()
@add_print_every_option()
@add_precision_option()
@add_backend_model_option()
@add_deterministic_option()
@add_allow_charge_mult_mismatch_option()
@click.pass_context
def cli(
    ctx: click.Context,
    input_path: Optional[Path],
    real_parm7: Optional[Path],
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
    csv_path: Optional[Path],
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
    print_every: Optional[int],
    precision: Optional[str],
    backend_model: Optional[str],
) -> None:
    _is_param_explicit = make_is_param_explicit(ctx)

    set_convert_file_enabled(convert_files)
    time_start = time.perf_counter()
    config_yaml, override_yaml, used_legacy_yaml = resolve_yaml_sources(
        config_yaml=config_yaml,
        override_yaml=None,
        args_yaml_legacy=None,
    )

    if csv_path is not None:
        final_dir = Path(out_dir).resolve()
        ensure_dir(final_dir)
        resolved_csv = Path(csv_path).resolve()
        try:
            df = pd.read_csv(resolved_csv)
        except Exception as exc:
            click.echo(f"[read] Failed to read CSV '{resolved_csv}': {exc}", err=True)
            sys.exit(1)
        click.echo(f"[read] Loaded precomputed grid from '{resolved_csv}'.")
        _finalize_surface_and_plot(
            df=df,
            final_dir=final_dir,
            baseline=baseline,
            zmin=zmin,
            zmax=zmax,
            d1_label_csv=None,
            d2_label_csv=None,
            d3_label_csv=None,
            write_surface_csv=False,
            time_start=time_start,
        )
        if out_json:
            from mlmm.core.utils import write_result_json
            min_energy = float(df["energy_hartree"].min()) if (not df.empty and "energy_hartree" in df.columns) else None
            result_data: Dict[str, Any] = {
                "status": "completed",
                "n_grid_points": len(df),
                "backend": None,
                "min_energy_hartree": min_energy,
                "files": {
                    "scan3d_density_html": "scan3d_density.html",
                },
            }
            write_result_json(
                final_dir, result_data,
                command="scan3d",
                elapsed_seconds=time.perf_counter() - time_start,
            )
        return

    if input_path is None:
        click.echo("ERROR: -i/--input is required unless --csv is provided.", err=True)
        sys.exit(1)
    if real_parm7 is None:
        click.echo("ERROR: --parm is required unless --csv is provided.", err=True)
        sys.exit(1)

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
                ligand_charge=ligand_charge, prefix="[scan3d]",
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
            # Honor the documented precedence defaults < --config YAML < CLI:
            # only override (already-YAML-merged) max_cycles when the user
            # explicitly passed --relax-max-cycles. Mirrors scan.py:554's
            # gate; an unconditional assignment silently clobbered the YAML
            # tier when the CLI default (10000) won.
            if _is_param_explicit("relax_max_cycles"):
                opt_cfg["max_cycles"] = int(relax_max_cycles)
                lbfgs_cfg["max_cycles"] = int(relax_max_cycles)
            if thresh is not None:
                opt_cfg["thresh"] = str(thresh)
            if bias_k is not None:
                bias_cfg["k"] = float(bias_k)

            out_dir_path = Path(opt_cfg["out_dir"]).resolve()

            calc_cfg["model_charge"] = int(charge)
            calc_cfg["model_mult"] = int(spin)
            calc_cfg["input_pdb"] = str(source_path)
            calc_cfg["real_parm7"] = str(real_parm7)
            if backend is not None:
                calc_cfg["backend"] = str(backend).lower()
            if precision is not None:
                from mlmm.backends import apply_precision_to_calc_cfg
                apply_precision_to_calc_cfg(calc_cfg, precision)
            if backend_model is not None:
                from mlmm.backends import apply_backend_model_to_calc_cfg
                apply_backend_model_to_calc_cfg(calc_cfg, backend_model)
            if _is_param_explicit("embedcharge"):
                calc_cfg["embedcharge"] = bool(embedcharge)
            if _is_param_explicit("embedcharge_cutoff"):
                calc_cfg["embedcharge_cutoff"] = embedcharge_cutoff
            if _is_param_explicit("print_every") and print_every is not None:
                opt_cfg["print_every"] = int(print_every)
            if link_atom_method is not None:
                calc_cfg["link_atom_method"] = str(link_atom_method).lower()
            if mm_backend is not None:
                calc_cfg["mm_backend"] = str(mm_backend).lower()
            if use_cmap is not None:
                calc_cfg["use_cmap"] = use_cmap

            try:
                model_pdb_path, layer_info = resolve_ml_layer_assignment(
                    source_path=source_path,
                    out_dir_path=out_dir_path,
                    model_pdb=model_pdb,
                    model_indices=model_indices,
                    detect_layer=detect_layer,
                    hess_cutoff=hess_cutoff,
                    movable_cutoff=movable_cutoff,
                    calc_cfg=calc_cfg,
                )
            except click.ClickException as e:
                click.echo(f"ERROR: {e.message}", err=True)
                sys.exit(1)

            freeze_atoms_final = apply_layer_freeze_constraints(
                geom_cfg,
                calc_cfg,
                layer_info,
                echo_fn=click.echo,
            )

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
                    expected_len=3,
                    one_based_default=one_based,
                    atom_meta=pdb_atom_meta,
                    option_name="--scan-lists",
                )
                scan_source = f"--scan-lists ({spec_path})"
            else:
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
                        {
                            "d1": (i1 + 1, j1 + 1, low1, high1),
                            "d2": (i2 + 1, j2 + 1, low2, high2),
                            "d3": (i3 + 1, j3 + 1, low3, high3),
                        },
                    )
                )
                # --print-parsed = "just show the parsed spec": exit before
                # any GPU calculation. scan3d has no --dry-run, so this is
                # also its only GPU-free spec-validation path.
                sys.exit(0)
            click.echo(
                pretty_block(
                    "scan-list",
                    {
                        "d1": (i1 + 1, j1 + 1, low1, high1),
                        "d2": (i2 + 1, j2 + 1, low2, high2),
                        "d3": (i3 + 1, j3 + 1, low3, high3),
                    },
                )
            )
            if pdb_atom_meta:
                click.echo("[scan3d] PDB atom details for scanned pairs:", detail=True)
                legend = PDB_ATOM_META_HEADER
                click.echo(f"        legend: {legend}", detail=True)
                click.echo(f"  d1 i: {format_pdb_atom_metadata(pdb_atom_meta, i1)}", detail=True)
                click.echo(f"     j: {format_pdb_atom_metadata(pdb_atom_meta, j1)}", detail=True)
                click.echo(f"  d2 i: {format_pdb_atom_metadata(pdb_atom_meta, i2)}", detail=True)
                click.echo(f"     j: {format_pdb_atom_metadata(pdb_atom_meta, j2)}", detail=True)
                click.echo(f"  d3 i: {format_pdb_atom_metadata(pdb_atom_meta, i3)}", detail=True)
                click.echo(f"     j: {format_pdb_atom_metadata(pdb_atom_meta, j3)}", detail=True)

            # Directory layout
            tmp_root = Path(tempfile.mkdtemp(prefix="scan3d_tmp_"))
            grid_dir = out_dir_path / "grid"
            tmp_opt_dir = tmp_root / "opt"
            ensure_dir(grid_dir)
            ensure_dir(tmp_opt_dir)
            final_dir = out_dir_path

            freeze = list(geom_cfg.get("freeze_atoms") or [])
            coord_type = geom_cfg.get("coord_type", GEOM_KW_DEFAULT["coord_type"])
            geom_outer = geom_loader(
                geom_input_path,
                coord_type=coord_type,
                freeze_atoms=freeze,
            )
            if freeze:
                try:
                    geom_outer.freeze_atoms = np.array(freeze, dtype=int)
                except Exception:
                    logger.debug("Failed to set freeze_atoms on geometry", exc_info=True)

            base_calc = mlmm(**calc_cfg)
            biased = HarmonicBiasCalculator(base_calc, k=float(bias_cfg["k"]))

            echo_resolved_device()

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

            # Measure reference distances on the (pre)optimized structure
            coords_outer = np.asarray(geom_outer.coords3d)
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
            click.echo(f"[grid] d1 steps = {N1}  values(A)={list(map(lambda x: f'{x:.3f}', d1_values))}", narrative=True)
            click.echo(f"[grid] d2 steps = {N2}  values(A)={list(map(lambda x: f'{x:.3f}', d2_values))}", narrative=True)
            click.echo(f"[grid] d3 steps = {N3}  values(A)={list(map(lambda x: f'{x:.3f}', d3_values))}", narrative=True)
            click.echo(f"[grid] total grid points = {N1 * N2 * N3}", narrative=True)

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
                    prefix=f"d1_{d1_tag}",
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
                        prefix=f"d1_{d1_tag}_d2_{d2_tag}",
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
                            prefix=f"d1_{d1_tag}_d2_{d2_tag}_d3_{d3_tag}",
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
            _finalize_surface_and_plot(
                df=df,
                final_dir=final_dir,
                baseline=baseline,
                zmin=zmin,
                zmax=zmax,
                d1_label_csv=d1_label_csv,
                d2_label_csv=d2_label_csv,
                d3_label_csv=d3_label_csv,
                write_surface_csv=True,
                time_start=time_start,
            )

            if out_json:
                from mlmm.core.utils import write_result_json
                min_energy = float(df["energy_hartree"].min()) if (not df.empty and "energy_hartree" in df.columns) else None
                result_data_main: Dict[str, Any] = {
                    "status": "completed",
                    "n_grid_points": len(df),
                    "pair1": {"i": int(i1 + 1), "j": int(j1 + 1), "low": float(low1), "high": float(high1)},
                    "pair2": {"i": int(i2 + 1), "j": int(j2 + 1), "low": float(low2), "high": float(high2)},
                    "pair3": {"i": int(i3 + 1), "j": int(j3 + 1), "low": float(low3), "high": float(high3)},
                    "backend": calc_cfg.get("backend", "uma"),
                    "charge": calc_cfg.get("model_charge"),
                    "spin": calc_cfg.get("model_mult"),
                    "min_energy_hartree": min_energy,
                    "files": {
                        "surface_csv": "surface.csv",
                        "scan3d_density_html": "scan3d_density.html",
                    },
                }
                write_result_json(
                    final_dir, result_data_main,
                    command="scan3d",
                    elapsed_seconds=time.perf_counter() - time_start,
                )

    except KeyboardInterrupt:
        click.echo("\nInterrupted by user.", err=True)
        sys.exit(130)
    except Exception as exc:
        render_cli_exception(exc, label="3D scan", out_dir=out_dir, command="scan3d", time_start=time_start)
    finally:
        if tmp_root is not None:
            shutil.rmtree(tmp_root, ignore_errors=True)
        # Release GPU memory so subsequent pipeline stages don't OOM
        base_calc = geom_outer = optimizer0 = None
        gc.collect()  # break cyclic refs inside torch.nn.Module
        if torch.cuda.is_available():
            torch.cuda.empty_cache()
