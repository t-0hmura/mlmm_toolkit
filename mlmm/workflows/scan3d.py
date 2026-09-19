"""
ML/MM three-coordinate grid scan with harmonic restraints.

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


def _result_calculator_fields(
    calc_cfg: Optional[Dict[str, Any]],
) -> Dict[str, Any]:
    """Return one stable calculator/charge schema for fresh and CSV-only runs."""
    if calc_cfg is None:
        return {
            "mlip_backend": None,
            "mlip_model": None,
            "mlip_model_label": None,
            "mlip_task": None,
            "mlip_precision": None,
            "mm_backend": None,
            "link_atom_method": None,
            "use_cmap": None,
            "charge": None,
            "spin": None,
        }

    from mlmm.core.utils import calculator_provenance

    return {
        **calculator_provenance(calc_cfg),
        "charge": calc_cfg.get("model_charge"),
        "spin": calc_cfg.get("model_mult"),
    }

import click
from mlmm.core.output import emit
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
from mlmm.workflows._outcomes import (
    attach_outcomes,
    is_preoptimization_record,
    make_scan_point,
    optimizer_converged_bit,
    scan_scientific_status,
    seed_eligible_mask,
)
from mlmm.workflows.charge_prep import resolve_charge_spin_or_raise
from mlmm.core.utils import (
    apply_ref_pdb_override,
    apply_layer_freeze_constraints,
    set_convert_file_enabled,
    apply_yaml_overrides,
    pretty_block,
    strip_inherited_keys,
    filter_calc_for_echo,
    format_freeze_atoms_for_echo,
    emit_dry_run_complete,
    format_elapsed,
    merge_freeze_atom_indices,
    prepare_input_structure,
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
    unique_tag_digits,
    values_from_bounds,
    unbiased_energy_hartree,
    snapshot_geometry,
    convert_and_annotate_xyz_to_pdb,
    echo_resolved_device,
)
from mlmm.workflows.scan_common import (
    add_scan_common_options,
    make_scan_lbfgs as _make_lbfgs,
    OutputCollisionError,
    prepare_grid_scan_output,
    prepare_scan_fixed_outputs,
    resolve_scan_optimizer_configs,
)
from mlmm.domain.scan_coordinates import (
    coordinate_atoms,
    coordinate_bounds,
    coordinate_kind,
    coordinate_step_cap,
    coordinate_unit,
    coordinate_value,
    format_coordinate,
)
from mlmm.cli.common_options import (
    add_ml_layer_detection_options,
    add_print_every_option,
    add_precision_option, add_backend_model_option, add_calc_file_option,
    add_workers_options,
    add_deterministic_option, add_allow_charge_mult_mismatch_option,
)
from mlmm.cli.common_options import add_dft_calculator_options
from mlmm.cli.decorators import resolve_yaml_sources, load_merged_yaml_cfg, make_is_param_explicit, render_cli_exception

# Shared defaults (copied from opt.py to keep ML/MM behavior consistent)
GEOM_KW: Dict[str, Any] = deepcopy(_OPT_GEOM_KW)
CALC_KW: Dict[str, Any] = deepcopy(_OPT_CALC_KW)
OPT_BASE_KW: Dict[str, Any] = deepcopy(_OPT_BASE_KW)
OPT_BASE_KW.update(
    {
        "out_dir": OUT_DIR_SCAN3D,
        "dump": False,
        "max_cycles": None,
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


def _explicit_true_series(series: pd.Series) -> pd.Series:
    return series.map(
        lambda value: value is True or str(value).strip().lower() == "true"
    )


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
) -> Dict[str, Any]:
    if df.empty:
        raise ValueError("No grid records were produced.")

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

    required_coords = ("d1_A", "d2_A", "d3_A")
    missing_coords = [name for name in required_coords if name not in df.columns]
    if missing_coords:
        raise ValueError(
            "surface.csv is missing coordinate column(s): "
            + ", ".join(missing_coords)
        )
    if "energy_hartree" not in df.columns and "energy_kcal" not in df.columns:
        raise ValueError(
            "surface.csv requires energy_hartree or energy_kcal."
        )

    grid_mask = pd.Series(
        [
            not is_preoptimization_record(record)
            for record in df.to_dict(orient="records")
        ],
        index=df.index,
    )

    coordinate_finite = np.ones(len(df), dtype=bool)
    for column in required_coords:
        coordinate_finite &= np.isfinite(df[column].to_numpy(dtype=float))
    energy_column = (
        "energy_hartree" if "energy_hartree" in df.columns else "energy_kcal"
    )
    finite_energy = np.isfinite(df[energy_column].to_numpy(dtype=float))
    complete_provenance = all(
        column in df.columns
        for column in ("bias_converged", "artifact_written")
    )

    usable_mask = grid_mask & coordinate_finite & finite_energy
    if "bias_converged" in df.columns:
        usable_mask &= _explicit_true_series(df["bias_converged"])
    if "artifact_written" in df.columns:
        usable_mask &= _explicit_true_series(df["artifact_written"])
    if not complete_provenance:
        click.echo(
            "[plot] WARNING: CSV lacks complete point provenance; using "
            "finite explicitly-converged rows in legacy plot-only mode.",
            err=True,
        )

    if not bool(usable_mask.any()):
        raise ValueError("No usable finite non-preoptimization grid point.")

    if baseline == "first":
        first_mask = (
            usable_mask
            & (df["i"] == 0)
            & (df["j"] == 0)
            & (df["k"] == 0)
        )
        if not bool(first_mask.any()):
            click.echo(
                "[baseline] 'first' requested but usable (i=0,j=0,k=0) "
                "is missing; using the usable minimum instead.",
                err=True,
            )
            ref_index = df.loc[usable_mask, energy_column].idxmin()
        else:
            ref_index = df.index[first_mask][0]
    else:
        ref_index = df.loc[usable_mask, energy_column].idxmin()

    if "energy_hartree" in df.columns:
        ref_energy = float(df.loc[ref_index, "energy_hartree"])
        df["energy_kcal"] = (
            df["energy_hartree"] - ref_energy
        ) * AU2KCALPERMOL
    else:
        ref_energy = float(df.loc[ref_index, "energy_kcal"])
        df["energy_kcal"] = df["energy_kcal"] - ref_energy

    if write_surface_csv:
        surface_csv = final_dir / "surface.csv"
        df["d1_label"] = d1_label_csv
        df["d2_label"] = d2_label_csv
        df["d3_label"] = d3_label_csv
        _csv_drop3 = [
            c
            for c in ("seed_eligible", "artifact_written", "geometry_file")
            if c in df.columns
        ]
        df.drop(columns=_csv_drop3).to_csv(surface_csv, index=False)
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
        & usable_mask.to_numpy(dtype=bool)
    )
    if not np.any(mask):
        raise ValueError("No finite data are available for plotting.")

    points = np.column_stack((d1_points[mask], d2_points[mask], d3_points[mask]))
    unique_points = np.unique(points, axis=0)
    if len(unique_points) != len(points):
        raise ValueError(
            "3D interpolation requires unique coordinate triples; found "
            f"{len(points) - len(unique_points)} duplicate row(s)."
        )
    support_rank = (
        int(np.linalg.matrix_rank(unique_points - unique_points[0]))
        if len(unique_points) > 1
        else 0
    )
    axis_spans = np.ptp(unique_points, axis=0) if len(unique_points) else np.zeros(3)
    if len(unique_points) < 4 or support_rank < 3 or np.any(axis_spans <= 0.0):
        raise ValueError(
            "3D interpolation requires at least four non-coplanar points "
            "spanning every axis; found "
            f"{len(unique_points)} unique point(s), rank {support_rank}, "
            f"axis spans {axis_spans.tolist()}."
        )

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
        autosize=True,
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
    fig3d.write_html(
        str(html3d),
        config={"responsive": True, "displaylogo": False},
        default_width="100%",
        default_height="100%",
    )
    click.echo(f"[plot] Wrote '{html3d}'.")

    emit("\n====== 3D Scan finished ======\n", narrative=True)
    min_energy_hartree = (
        float(df.loc[usable_mask, "energy_hartree"].min())
        if "energy_hartree" in df.columns
        else None
    )
    return {
        "complete_provenance": bool(complete_provenance),
        "n_grid_points": int(np.count_nonzero(grid_mask)),
        "n_points_usable": int(np.count_nonzero(usable_mask)),
        "min_energy_hartree": min_energy_hartree,
    }


@click.command(
    help="3D internal-coordinate scan with harmonic restraints using ML/MM.",
    context_settings={"help_option_names": ["-h", "--help"]},
)
@click.option(
    "-i",
    "--input",
    "input_path",
    type=click.Path(path_type=Path, exists=True, dir_okay=False),
    required=False,
    help=("Input PDB/mmCIF, or XYZ with --ref-pdb. Required unless --csv is "
          "provided."),
)
@click.option(
    "--parm7",
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
    help="ML-only, link-H-free PDB subset; atom identity/order must match the "
         "full PDB/parm7. When provided, it defines ML membership; "
         "--detect-layer still reads valid movable/frozen MM B-factors.",
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
              help="ML-region total charge. Required unless --ligand-charge or plot-only --csv is provided.")
@click.option("-l", "--ligand-charge", type=str, default=None, show_default=False,
              help="Total charge for unknown ligand residues or a per-resname mapping "
                   "(e.g., GPP:-3,SAM:1), used to derive the ML-region charge when -q "
                   "is omitted (requires PDB input or --ref-pdb).")
@click.option(
    "-m",
    "--multiplicity",
    "spin",
    type=int,
    default=None,
    show_default="1",
    help="Spin multiplicity (2S+1) for the ML region.",
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
    "--movable-cutoff",
    "movable_cutoff",
    type=float,
    default=None,
    show_default="use freeze_atoms",
    help="Distance cutoff (Å) from ML region for movable MM atoms. MM atoms beyond this are frozen. "
         "Providing --movable-cutoff disables --detect-layer.",
)
@click.option(
    "-s", "--scan-lists",
    "scan_list_raw",
    type=str,
    required=False,
    help=(
        "Three scan ranges as an inline literal or YAML/JSON file: distance "
        "(i,j,low,high), angle (i,j,k,low,high), or dihedral "
        "(i,j,k,l,low,high). Distances use Å; angles and dihedrals use degrees."
    ),
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
    "--dry-run/--no-dry-run",
    "dry_run",
    default=False,
    show_default=True,
    help="Validate options and print the execution plan without running the scan.",
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
    type=click.Choice(["uma", "orb", "mace", "aimnet2", "dft"], case_sensitive=False),
    default=None,
    show_default="uma",
    help="High-level backend for the ONIOM model region.",
)
@click.option(
    "--embedcharge/--no-embedcharge",
    "embedcharge",
    default=False,
    show_default=True,
    help="Enable electrostatic embedding. MLIP backends use the experimental, computationally expensive xTB point-charge delta correction; dft uses native PySCF MM point charges.",
)
@click.option(
    "--embedcharge-cutoff",
    "embedcharge_cutoff",
    type=float,
    default=None,
    show_default="12.0",
    help="Distance cutoff (Å) from the ML region for MM point charges used by embedding.",
)
@click.option(
    "--link-atom-method",
    "link_atom_method",
    type=click.Choice(["scaled", "fixed"], case_sensitive=False),
    default=None,
    show_default="scaled",
    help="Link-atom position mode: scaled (g-factor) or fixed (1.09/1.01 Å).",
)
@click.option(
    "--mm-backend",
    "mm_backend",
    type=click.Choice(["hessian_ff", "openmm"], case_sensitive=False),
    default=None,
    show_default="hessian_ff",
    help="MM backend. MM Hessians use finite differences by default; set calc.mm_fd: false for the hessian_ff analytical path.",
)
@click.option(
    "--cmap/--no-cmap",
    "use_cmap",
    default=None,
    show_default="cmap",
    help="Preserve CMAP terms in both real and model MM layers when present in parm7.",
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
@add_workers_options()
@add_backend_model_option()
@add_calc_file_option()
@add_deterministic_option()
@add_allow_charge_mult_mismatch_option()
@add_dft_calculator_options()
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
    movable_cutoff: Optional[float],
    scan_list_raw: Optional[str],
    csv_path: Optional[Path],
    one_based: bool,
    print_parsed: bool,
    dry_run: bool,
    max_step_size: float,
    max_angle_step_size: float,
    max_dihedral_step_size: float,
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
    workers: Optional[int],
    workers_per_node: Optional[int],
    backend_model: Optional[str],
    calc_file: Optional[str],
    calc_factory: Optional[str],
) -> None:
    _is_param_explicit = make_is_param_explicit(ctx)

    set_convert_file_enabled(convert_files)
    time_start = time.perf_counter()
    config_yaml, override_yaml, used_legacy_yaml = resolve_yaml_sources(
        config_yaml=config_yaml,
        override_yaml=None,
        args_yaml_legacy=None,
    )
    yaml_cfg, _, _ = load_merged_yaml_cfg(
        config_yaml=config_yaml,
        override_yaml=None,
    )
    from mlmm.cli.decorators import resolve_model_indices_setting
    model_indices_str, model_indices_one_based = resolve_model_indices_setting(
        ctx, yaml_cfg, model_indices_str, model_indices_one_based
    )

    if csv_path is not None:
        final_dir = Path(out_dir).resolve()
        resolved_csv = Path(csv_path).resolve()
        try:
            final_dir = prepare_scan_fixed_outputs(
                final_dir,
                fixed_names=(
                    "scan3d_density.html",
                    "result.json",
                    "summary.json",
                ),
                protected_inputs=(
                    resolved_csv,
                    config_yaml,
                    override_yaml,
                ),
            )
            df = pd.read_csv(resolved_csv)
            click.echo(f"[read] Loaded precomputed grid from '{resolved_csv}'.")
            surface_stats = _finalize_surface_and_plot(
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
                result_data: Dict[str, Any] = {
                    "status": "completed",
                    "energy_reference": "bare_mlmm_pes",
                    "n_grid_points": surface_stats["n_grid_points"],
                    **_result_calculator_fields(None),
                    "min_energy_hartree": surface_stats["min_energy_hartree"],
                    "files": {
                        "scan3d_density_html": "scan3d_density.html",
                    },
                    "current_output_paths": ["scan3d_density.html"],
                }
                if surface_stats["complete_provenance"]:
                    result_data["n_points_usable"] = surface_stats[
                        "n_points_usable"
                    ]
                write_result_json(
                    final_dir, result_data,
                    command="scan3d",
                    elapsed_seconds=time.perf_counter() - time_start,
                )
            emit(
                format_elapsed("[time] Elapsed Time for 3D Scan", time_start),
                narrative=True,
            )
        except KeyboardInterrupt:
            click.echo("\nInterrupted by user.", err=True)
            sys.exit(130)
        except OutputCollisionError:
            raise
        except Exception as exc:
            render_cli_exception(
                exc,
                label="3D scan",
                out_dir=final_dir,
                command="scan3d",
                time_start=time_start,
            )
        return

    if input_path is None:
        click.echo("ERROR: -i/--input is required unless --csv is provided.", err=True)
        sys.exit(1)
    if real_parm7 is None:
        click.echo("ERROR: --parm is required unless --csv is provided.", err=True)
        sys.exit(1)

    suffix = input_path.suffix.lower()
    if suffix not in (".pdb", ".cif", ".mmcif", ".xyz"):
        click.echo("ERROR: --input must be a PDB, mmCIF, or XYZ file.", err=True)
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

            geom_cfg = dict(GEOM_KW)
            calc_cfg = dict(CALC_KW)
            bias_cfg = dict(BIAS_KW)

            opt_cfg, lbfgs_cfg = resolve_scan_optimizer_configs(
                yaml_cfg,
                opt_defaults=OPT_BASE_KW,
                lbfgs_defaults=LBFGS_KW,
                thresh=thresh,
                relax_max_cycles=relax_max_cycles,
                is_param_explicit=_is_param_explicit,
            )

            apply_yaml_overrides(
                yaml_cfg,
                [
                    (geom_cfg, (("geom",),)),
                    (calc_cfg, (("calc",), ("mlmm",))),
                    (bias_cfg, (("bias",),)),
                ],
            )
            detect_layer_effective = (
                bool(detect_layer)
                if _is_param_explicit("detect_layer")
                else bool(calc_cfg.get("use_bfactor_layers", True))
            )
            charge, spin = resolve_charge_spin_or_raise(
                prepared_input, charge, spin,
                ligand_charge=ligand_charge, prefix="[scan3d]",
                model_pdb=model_pdb,
                model_indices_spec=model_indices_str,
                detect_layer=detect_layer_effective,
                yaml_cfg=yaml_cfg,
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

            if _is_param_explicit("out_dir"):
                opt_cfg["out_dir"] = out_dir
            opt_cfg["dump"] = False
            if bias_k is not None:
                bias_cfg["k"] = float(bias_k)

            out_dir_path = Path(opt_cfg["out_dir"]).resolve()
            spec_path = (
                Path(scan_list_raw)
                if scan_list_raw is not None
                and is_scan_spec_file(scan_list_raw)
                else None
            )

            calc_cfg["model_charge"] = int(charge)
            calc_cfg["model_mult"] = int(spin)
            calc_cfg["input_pdb"] = str(source_path)
            calc_cfg["real_parm7"] = str(real_parm7)
            if backend is not None:
                calc_cfg["backend"] = str(backend).lower()
            from mlmm.backends import apply_precision_to_calc_cfg
            # Always run so a YAML-set workers>1 also gets the analytical-Hessian guard.
            from mlmm.backends import apply_workers_to_calc_cfg
            apply_workers_to_calc_cfg(calc_cfg, workers, workers_per_node)
            from mlmm.backends import apply_backend_model_to_calc_cfg
            # Unconditional: also pops a raw backend_model token from a --config YAML.
            apply_backend_model_to_calc_cfg(calc_cfg, backend_model)
            # --calc-file overrides --backend with a user ASE Calculator (custom backend).
            from mlmm.backends import apply_calc_file_to_calc_cfg
            apply_calc_file_to_calc_cfg(calc_cfg, calc_file, calc_factory)
            # Unconditional: also dispatches a --config YAML calc.precision
            # after the final backend has been selected.
            apply_precision_to_calc_cfg(calc_cfg, precision)
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
            from mlmm.core.dft_settings import finalize_dft_calculator_config
            finalize_dft_calculator_config(
                ctx, calc_cfg, output_dir=out_dir_path
            )

            try:
                model_pdb_path, layer_info = resolve_ml_layer_assignment(
                    source_path=source_path,
                    out_dir_path=out_dir_path,
                    model_pdb=model_pdb,
                    model_indices=model_indices,
                    detect_layer=detect_layer_effective,
                    hess_cutoff=calc_cfg.get("hess_cutoff"),
                    movable_cutoff=movable_cutoff,
                    calc_cfg=calc_cfg,
                    protected_inputs=(
                        input_path,
                        source_path,
                        geom_input_path,
                        real_parm7,
                        model_pdb,
                        ref_pdb,
                        spec_path,
                        config_yaml,
                        override_yaml,
                        (
                            Path(calc_cfg["calc_file"])
                            if calc_cfg.get("calc_file")
                            else None
                        ),
                    ),
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
            if spec_path is not None:
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
            axis1, axis2, axis3 = parsed
            kinds = [coordinate_kind(axis, is_range=True) for axis in parsed]
            units = ["Å" if kind == "distance" else "deg" for kind in kinds]
            atoms1, atoms2, atoms3 = [coordinate_atoms(axis, is_range=True) for axis in parsed]
            (low1, high1), (low2, high2), (low3, high3) = [coordinate_bounds(axis) for axis in parsed]
            def _axis_payload(kind_axis, atoms_axis, low, high):
                payload = {
                    "kind": kind_axis,
                    "atoms_1based": [int(i) + 1 for i in atoms_axis],
                    "unit": coordinate_unit(kind_axis),
                    "low": float(low), "high": float(high),
                }
                if kind_axis == "distance":
                    payload.update(i=int(atoms_axis[0] + 1), j=int(atoms_axis[1] + 1))
                return payload
            frozen_set = set(map(int, freeze_atoms_final))
            for axis, entry in enumerate(parsed, start=1):
                atoms = coordinate_atoms(entry, is_range=True)
                if all(int(atom) in frozen_set for atom in atoms):
                    raise click.BadParameter(
                        "A scan restraint cannot contain only frozen atoms: "
                        f"axis d{axis}, atoms {[int(atom) + 1 for atom in atoms]}."
                    )
            labels = []
            for name, kind_axis, atoms_axis, raw_axis in zip(("d1", "d2", "d3"), kinds, (atoms1, atoms2, atoms3), raw_pairs):
                labels.append(axis_label_csv(name, *atoms_axis, scan_one_based, pdb_atom_meta, raw_axis)
                              if kind_axis == "distance" else f"{name}_{kind_axis}_{'_'.join(str(i + 1) for i in atoms_axis)}_deg")
            d1_label_csv, d2_label_csv, d3_label_csv = labels
            if print_parsed:
                click.echo(
                    pretty_block(
                        "scan-parsed",
                        {
                            "source": scan_source,
                            "one_based": bool(scan_one_based),
                            "pairs_0based": parsed,
                        },
                    force=True)
                )
                click.echo(
                    pretty_block(
                        "scan-list",
                        {
                            "d1": format_coordinate(axis1, is_range=True),
                            "d2": format_coordinate(axis2, is_range=True),
                            "d3": format_coordinate(axis3, is_range=True),
                        },
                    force=True)
                )
                # --print-parsed = "just show the parsed spec": exit before
                # any GPU calculation.
                if dry_run:
                    emit_dry_run_complete()
                else:
                    emit(
                        format_elapsed("[time] Elapsed Time for 3D Scan", time_start),
                        narrative=True,
                    )
                sys.exit(0)
            if dry_run:
                click.echo(
                    pretty_block(
                        "dry_run_plan",
                        {
                            "input_geometry": str(source_path),
                            "output_dir": str(out_dir_path),
                            "scan_source": scan_source,
                            "one_based": bool(scan_one_based),
                            "d1_0based": tuple(axis1),
                            "d2_0based": tuple(axis2),
                            "d3_0based": tuple(axis3),
                            "charge": int(charge),
                            "spin": int(spin),
                            "detect_layer": bool(detect_layer_effective),
                            "backend": calc_cfg.get("backend", "uma"),
                            "embedcharge": bool(calc_cfg.get("embedcharge", False)),
                        },
                        force=True,
                    )
                )
                emit_dry_run_complete()
                return
            click.echo(
                pretty_block(
                    "scan-list",
                    {
                        "d1": format_coordinate(axis1, is_range=True),
                        "d2": format_coordinate(axis2, is_range=True),
                        "d3": format_coordinate(axis3, is_range=True),
                    },
                )
            )
            if pdb_atom_meta:
                emit("[scan3d] PDB atom details for scanned pairs:", detail=True)
                legend = PDB_ATOM_META_HEADER
                emit(f"        legend: {legend}", detail=True)
                for label, atoms in zip(("d1", "d2", "d3"), (atoms1, atoms2, atoms3)):
                    for pos, atom in enumerate(atoms, start=1):
                        emit(f"  {label} atom {pos}: {format_pdb_atom_metadata(pdb_atom_meta, atom)}", detail=True)

            # Directory layout
            tmp_root = Path(tempfile.mkdtemp(prefix="scan3d_tmp_"))
            final_dir, grid_dir = prepare_grid_scan_output(
                out_dir_path,
                fixed_names=(
                    "surface.csv",
                    "scan3d_density.html",
                    "result.json",
                    "summary.json",
                ),
                protected_inputs=(
                    input_path,
                    source_path,
                    geom_input_path,
                    real_parm7,
                    model_pdb,
                    model_pdb_path,
                    ref_pdb,
                    spec_path,
                    config_yaml,
                    override_yaml,
                    (
                        Path(calc_cfg["calc_file"])
                        if calc_cfg.get("calc_file")
                        else None
                    ),
                ),
            )
            tmp_opt_dir = tmp_root / "opt"
            ensure_dir(tmp_opt_dir)

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

            # The reference/anchor structure is usable-by-default when no preopt
            # is requested; when preopt runs, its reported convergence bit
            # replaces the default.
            _preopt_conv: Optional[bool] = True
            if preopt:
                preopt_input = _snapshot_geometry(geom_outer)
                click.echo("[preopt] Unbiased relaxation of the initial structure ...")
                geom_outer.set_calculator(base_calc)
                optimizer0 = _make_lbfgs(
                    geom_outer,
                    lbfgs_cfg,
                    opt_cfg,
                    max_step_bohr=float(max_step_size) * ANG2BOHR,
                    out_dir=tmp_opt_dir,
                    prefix="preopt",
                )
                try:
                    optimizer0.run()
                    _preopt_conv = optimizer_converged_bit(optimizer0)
                except ZeroStepLength:
                    click.echo("[preopt] ZeroStepLength — continuing.", err=True)
                    _preopt_conv = optimizer_converged_bit(optimizer0)
                except OptimizationError as exc:
                    click.echo(f"[preopt] OptimizationError — {exc}", err=True)
                    _preopt_conv = False
                try:
                    preopt_energy_check = unbiased_energy_hartree(
                        geom_outer, base_calc
                    )
                    preopt_state_finite = bool(
                        np.all(np.isfinite(np.asarray(geom_outer.coords3d)))
                        and np.isfinite(float(preopt_energy_check))
                    )
                except Exception:
                    preopt_state_finite = False
                if _preopt_conv is not True or not preopt_state_finite:
                    click.echo(
                        "[preopt] Preoptimization was not a finite converged "
                        "result; restoring the input geometry for the scan.",
                        err=True,
                    )
                    geom_outer = _snapshot_geometry(preopt_input)

            records: List[Dict[str, Any]] = []

            # Measure reference distances on the (pre)optimized structure
            coords_outer = np.asarray(geom_outer.coords3d)
            d1_ref = coordinate_value(coords_outer, axis1, is_range=True)
            d2_ref = coordinate_value(coords_outer, axis2, is_range=True)
            d3_ref = coordinate_value(coords_outer, axis3, is_range=True)

            if math.isfinite(d1_ref) and math.isfinite(d2_ref) and math.isfinite(d3_ref):
                click.echo(
                    f"[center] reference coordinate values: "
                    f"d1 = {d1_ref:.3f} {units[0]}, d2 = {d2_ref:.3f} {units[1]}, "
                    f"d3 = {d3_ref:.3f} {units[2]}"
                )

                # Write preoptimized structure
                d1_ref_tag = distance_tag(d1_ref)
                d2_ref_tag = distance_tag(d2_ref)
                d3_ref_tag = distance_tag(d3_ref)
                preopt_xyz_path = grid_dir / f"preopt_i{d1_ref_tag}_j{d2_ref_tag}_k{d3_ref_tag}.xyz"
                preopt_artifact_written = False
                try:
                    xyz_pre = geom_outer.as_xyz()
                    if not xyz_pre.endswith("\n"):
                        xyz_pre += "\n"
                    with open(preopt_xyz_path, "w") as handle:
                        handle.write(xyz_pre)
                    preopt_artifact_written = True

                    if convert_files:
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
                        "bias_converged": _preopt_conv,
                        "artifact_written": preopt_artifact_written,
                        "is_preopt": True,
                    }
                )
            else:
                click.echo(
                    "[center] WARNING: failed to determine reference distances; using grid order as-is.",
                    err=True,
                )

            # Build distance grids and reorder so that scanning starts near the reference structure
            d1_values = values_from_bounds(low1, high1, coordinate_step_cap(kinds[0], max_step_size, max_angle_step_size, max_dihedral_step_size))
            d2_values = values_from_bounds(low2, high2, coordinate_step_cap(kinds[1], max_step_size, max_angle_step_size, max_dihedral_step_size))
            d3_values = values_from_bounds(low3, high3, coordinate_step_cap(kinds[2], max_step_size, max_angle_step_size, max_dihedral_step_size))

            # One tag precision per axis, so a fine grid cannot map two targets
            # onto the same point tag and truncate the earlier artifact.
            d1_digits = unique_tag_digits(d1_values)
            d2_digits = unique_tag_digits(d2_values)
            d3_digits = unique_tag_digits(d3_values)

            def _d1_tag(value: float) -> str:
                return distance_tag(value, digits=d1_digits, pad=d1_digits + 1)

            def _d2_tag(value: float) -> str:
                return distance_tag(value, digits=d2_digits, pad=d2_digits + 1)

            def _d3_tag(value: float) -> str:
                return distance_tag(value, digits=d3_digits, pad=d3_digits + 1)

            if math.isfinite(d1_ref):
                d1_values = np.array(sorted(d1_values, key=lambda v: abs(v - d1_ref)), dtype=float)
            if math.isfinite(d2_ref):
                d2_values = np.array(sorted(d2_values, key=lambda v: abs(v - d2_ref)), dtype=float)
            if math.isfinite(d3_ref):
                d3_values = np.array(sorted(d3_values, key=lambda v: abs(v - d3_ref)), dtype=float)

            N1, N2, N3 = len(d1_values), len(d2_values), len(d3_values)
            emit(f"[grid] d1 steps = {N1}  values({units[0]})={list(map(lambda x: f'{x:.3f}', d1_values))}", narrative=True)
            emit(f"[grid] d2 steps = {N2}  values({units[1]})={list(map(lambda x: f'{x:.3f}', d2_values))}", narrative=True)
            emit(f"[grid] d3 steps = {N3}  values({units[2]})={list(map(lambda x: f'{x:.3f}', d3_values))}", narrative=True)
            emit(f"[grid] total grid points = {N1 * N2 * N3}", narrative=True)

            max_step_bohr = float(max_step_size) * ANG2BOHR

            # Caches for nearest-neighbor starting geometries
            d1_geoms: Dict[int, Any] = {}
            d2_geoms: Dict[int, Dict[int, Any]] = {}
            d3_geoms: Dict[Tuple[int, int], Dict[int, Any]] = {}

            geom_outer_initial = _snapshot_geometry(geom_outer)

            # ===== 3D nested scan: d1 (outer) → d2 (middle) → d3 (inner) =====
            for i_idx, d1_target in enumerate(d1_values):
                d1_tag = _d1_tag(d1_target)
                click.echo(f"\n--- d1 step {i_idx + 1}/{N1} : target = {d1_target:.3f} {units[0]} ---")

                # Choose initial geometry for this d1
                if not d1_geoms:
                    geom_outer_i = _snapshot_geometry(geom_outer_initial)
                else:
                    nearest_i = min(d1_geoms.keys(), key=lambda p: abs(d1_values[p] - d1_target))
                    geom_outer_i = _snapshot_geometry(d1_geoms[nearest_i])

                biased.set_restraints([(*atoms1, float(d1_target))])
                geom_outer_i.set_calculator(biased)
                geom_outer_start = _snapshot_geometry(geom_outer_i)

                opt1 = _make_lbfgs(
                    geom_outer_i,
                    lbfgs_cfg,
                    opt_cfg,
                    max_step_bohr=max_step_bohr,
                    out_dir=tmp_opt_dir,
                    prefix=f"d1_{d1_tag}",
                )
                d1_converged = None
                try:
                    opt1.run()
                    d1_converged = optimizer_converged_bit(opt1)
                except ZeroStepLength:
                    click.echo(f"[d1 {i_idx}] ZeroStepLength — continuing to d2/d3 scan.", err=True)
                    d1_converged = optimizer_converged_bit(opt1)
                except OptimizationError as exc:
                    click.echo(f"[d1 {i_idx}] OptimizationError — {exc}", err=True)
                    d1_converged = False

                if d1_converged is True and np.isfinite(
                    np.asarray(geom_outer_i.coords3d, dtype=float)
                ).all():
                    geom_after_d1 = _snapshot_geometry(geom_outer_i)
                    d1_geoms[i_idx] = geom_after_d1
                else:
                    geom_after_d1 = _snapshot_geometry(geom_outer_start)

                if i_idx not in d2_geoms:
                    d2_geoms[i_idx] = {}

                for j_idx, d2_target in enumerate(d2_values):
                    d2_tag = _d2_tag(d2_target)
                    click.echo(
                        f"  [stage] d1/d2 step ({i_idx + 1}/{N1}, {j_idx + 1}/{N2}): "
                        f"targets = ({d1_target:.3f} {units[0]}, {d2_target:.3f} {units[1]})"
                    )

                    # Choose initial geometry for this (d1,d2)
                    d2_store = d2_geoms[i_idx]
                    if not d2_store:
                        geom_mid = _snapshot_geometry(geom_after_d1)
                    else:
                        nearest_j = min(d2_store.keys(), key=lambda p: abs(d2_values[p] - d2_target))
                        geom_mid = _snapshot_geometry(d2_store[nearest_j])

                    biased.set_restraints([
                        (*atoms1, float(d1_target)),
                        (*atoms2, float(d2_target)),
                    ])
                    geom_mid.set_calculator(biased)
                    geom_mid_start = _snapshot_geometry(geom_mid)

                    opt2 = _make_lbfgs(
                        geom_mid,
                        lbfgs_cfg,
                        opt_cfg,
                        max_step_bohr=max_step_bohr,
                        out_dir=tmp_opt_dir,
                        prefix=f"d1_{d1_tag}_d2_{d2_tag}",
                    )
                    d2_converged = None
                    try:
                        opt2.run()
                        d2_converged = optimizer_converged_bit(opt2)
                    except ZeroStepLength:
                        click.echo(f"[d1 {i_idx}, d2 {j_idx}] ZeroStepLength — continuing to d3 scan.", err=True)
                        d2_converged = optimizer_converged_bit(opt2)
                    except OptimizationError as exc:
                        click.echo(f"[d1 {i_idx}, d2 {j_idx}] OptimizationError — {exc}", err=True)
                        d2_converged = False

                    if d2_converged is True and np.isfinite(
                        np.asarray(geom_mid.coords3d, dtype=float)
                    ).all():
                        geom_after_d2 = _snapshot_geometry(geom_mid)
                        d2_store[j_idx] = geom_after_d2
                    else:
                        geom_after_d2 = _snapshot_geometry(geom_mid_start)

                    key_ij = (i_idx, j_idx)
                    if key_ij not in d3_geoms:
                        d3_geoms[key_ij] = {}
                    d3_store = d3_geoms[key_ij]

                    trj_blocks = [] if dump else None

                    for k_idx, d3_target in enumerate(d3_values):
                        d3_tag = _d3_tag(d3_target)

                        # Choose initial geometry for this (d1,d2,d3)
                        if not d3_store:
                            geom_inner = _snapshot_geometry(geom_after_d2)
                        else:
                            nearest_k = min(d3_store.keys(), key=lambda p: abs(d3_values[p] - d3_target))
                            geom_inner = _snapshot_geometry(d3_store[nearest_k])

                        biased.set_restraints([
                            (*atoms1, float(d1_target)),
                            (*atoms2, float(d2_target)),
                            (*atoms3, float(d3_target)),
                        ])
                        geom_inner.set_calculator(biased)

                        opt3 = _make_lbfgs(
                            geom_inner,
                            lbfgs_cfg,
                            opt_cfg,
                            max_step_bohr=max_step_bohr,
                            out_dir=tmp_opt_dir,
                            prefix=f"d1_{d1_tag}_d2_{d2_tag}_d3_{d3_tag}",
                        )
                        # a normal (non-raising) run is NOT convergence —
                        # read the optimizer's explicit tri-state bit.
                        converged: Optional[bool] = None
                        try:
                            opt3.run()
                            converged = optimizer_converged_bit(opt3)
                        except ZeroStepLength:
                            click.echo(
                                f"    [d1 {i_idx}, d2 {j_idx}, d3 {k_idx}] ZeroStepLength — recorded anyway.",
                                err=True,
                            )
                            converged = optimizer_converged_bit(opt3)
                        except OptimizationError as exc:
                            click.echo(
                                f"    [d1 {i_idx}, d2 {j_idx}, d3 {k_idx}] OptimizationError — {exc}",
                                err=True,
                            )
                            converged = False

                        energy_h = unbiased_energy_hartree(geom_inner, base_calc)

                        xyz_path = grid_dir / f"point_i{d1_tag}_j{d2_tag}_k{d3_tag}.xyz"
                        _artifact_written = False
                        try:
                            xyz = geom_inner.as_xyz()
                            if not xyz.endswith("\n"):
                                xyz += "\n"
                            with open(xyz_path, "w") as handle:
                                handle.write(xyz)
                            _artifact_written = True

                            if convert_files:
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

                        # Reuse only scientifically usable points. A converged
                        # state with a non-finite energy or no geometry artifact
                        # must not seed a later grid point.
                        if (
                            converged is True
                            and math.isfinite(energy_h)
                            and np.isfinite(
                                np.asarray(geom_inner.coords3d, dtype=float)
                            ).all()
                            and _artifact_written
                        ):
                            d3_store[k_idx] = _snapshot_geometry(geom_inner)

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
                                "bias_converged": converged,
                                "artifact_written": bool(_artifact_written),
                                "geometry_file": (
                                    str(Path("grid") / xyz_path.name)
                                    if _artifact_written
                                    else None
                                ),
                                "is_preopt": False,
                            }
                        )

                    if dump and trj_blocks:
                        trj_path = grid_dir / f"inner_path_d1_{d1_tag}_d2_{d2_tag}_trj.xyz"
                        try:
                            with open(trj_path, "w") as handle:
                                handle.write("".join(trj_blocks))
                            click.echo(f"[write] Wrote '{trj_path}'.")

                            if convert_files:
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
            for axis_number, axis_kind in enumerate(kinds, start=1):
                df[f"q{axis_number}"] = df[f"d{axis_number}_A"]
                df[f"q{axis_number}_unit"] = coordinate_unit(axis_kind)
            surface_stats = _finalize_surface_and_plot(
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
                grid_records = (
                    [
                        rec
                        for rec in records
                        if not bool(rec.get("is_preopt", False))
                    ]
                    if csv_path is None
                    else []
                )
                result_data_main: Dict[str, Any] = {
                    "status": "completed",
                    "energy_reference": "bare_mlmm_pes",
                    "n_grid_points": surface_stats["n_grid_points"],
                    "pair1": _axis_payload(kinds[0], atoms1, low1, high1),
                    "pair2": _axis_payload(kinds[1], atoms2, low2, high2),
                    "pair3": _axis_payload(kinds[2], atoms3, low3, high3),
                    **_result_calculator_fields(calc_cfg),
                    "min_energy_hartree": surface_stats["min_energy_hartree"],
                    "files": {
                        "surface_csv": "surface.csv",
                        "scan3d_density_html": "scan3d_density.html",
                    },
                }
                grid_geometry_files = [
                    str(rec["geometry_file"])
                    for rec in grid_records
                    if rec.get("geometry_file")
                ]
                result_data_main["grid_points"] = []
                for rec in grid_records:
                    point = {
                        "index": [int(rec["i"]), int(rec["j"]), int(rec["k"])],
                        "coordinate_values": [
                            float(rec["d1_A"]),
                            float(rec["d2_A"]),
                            float(rec["d3_A"]),
                        ],
                        "coordinate_targets": [
                            float(rec["d1_A"]),
                            float(rec["d2_A"]),
                            float(rec["d3_A"]),
                        ],
                        "coordinate_units": [coordinate_unit(kind) for kind in kinds],
                        "energy_hartree": rec.get("energy_hartree"),
                        "converged": rec.get("bias_converged"),
                        "geometry_file": rec.get("geometry_file"),
                    }
                    if all(kind == "distance" for kind in kinds):
                        point["distances_angstrom"] = list(point["coordinate_values"])
                        point["targets_angstrom"] = list(point["coordinate_targets"])
                    result_data_main["grid_points"].append(point)
                result_data_main["current_output_paths"] = [
                    "surface.csv",
                    "scan3d_density.html",
                    *grid_geometry_files,
                ]
                # Additive outcome fields: every attempted point and aggregate
                # scientific_status. Legacy ``status`` stays "completed".
                _point_outcomes3 = [
                    make_scan_point(
                        f"i{rec.get('i')}_j{rec.get('j')}_k{rec.get('k')}",
                        executed=True,
                        converged=rec.get("bias_converged"),
                        energy=rec.get("energy_hartree"),
                        artifact_written=bool(rec.get("artifact_written", False)),
                    )
                    for rec in grid_records
                ]
                _sci3, _sci3_reasons = scan_scientific_status(_point_outcomes3)
                result_data_main["execution_status"] = "completed"
                result_data_main["n_points_attempted"] = len(_point_outcomes3)
                result_data_main["n_points_usable"] = sum(
                    1 for p in _point_outcomes3 if p.seed_eligible
                )
                attach_outcomes(
                    result_data_main,
                    point_outcomes=_point_outcomes3,
                    scientific_status=_sci3,
                    scientific_status_reasons=_sci3_reasons,
                )
                write_result_json(
                    final_dir, result_data_main,
                    command="scan3d",
                    elapsed_seconds=time.perf_counter() - time_start,
                )
            emit(
                format_elapsed("[time] Elapsed Time for 3D Scan", time_start),
                narrative=True,
            )

    except KeyboardInterrupt:
        click.echo("\nInterrupted by user.", err=True)
        sys.exit(130)
    except OutputCollisionError:
        raise
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
