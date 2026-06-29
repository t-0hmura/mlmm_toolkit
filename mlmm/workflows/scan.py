"""
ML/MM staged bond-length scan with harmonic restraints.

Example:
    mlmm scan -i pocket.pdb --parm real.parm7 --model-pdb ml_region.pdb -q 0 --scan-lists "[(12,45,2.20)]"

For detailed documentation, see: docs/scan.md
"""

from __future__ import annotations

from copy import deepcopy
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional, Sequence, Tuple

import gc
import logging
import math
import sys
import textwrap

logger = logging.getLogger(__name__)

import click
import numpy as np
import time
import torch

from pysisyphus.helpers import geom_loader
from pysisyphus.optimizers.LBFGS import LBFGS
from pysisyphus.optimizers.exceptions import OptimizationError, ZeroStepLength
from pysisyphus.constants import BOHR2ANG, ANG2BOHR

from mlmm.backends.mlmm_calc import mlmm
from mlmm.core.defaults import (
    BIAS_KW as _BIAS_KW_DEFAULT,
    BOND_KW as _BOND_KW_DEFAULT,
    GEOM_KW_DEFAULT,
    OUT_DIR_SCAN,
    THRESH_CHOICES,
)
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
    convert_xyz_to_pdb,
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
    collect_single_option_values,
    load_pdb_atom_metadata,
    parse_scan_list_triples,
    parse_scan_spec_stages,
    is_scan_spec_file,
    PDB_ATOM_META_HEADER,
    format_pdb_atom_metadata,
    parse_indices_string,
    resolve_ml_layer_assignment,
    snapshot_geometry,
    echo_resolved_device,
)
from mlmm.domain.bond_changes import compare_structures, summarize_changes
from mlmm.cli.common_options import (
    add_ml_charge_spin_options,
    add_ml_layer_detection_options,
    add_coord_type_option,
    add_print_every_option,
    add_precision_option, add_backend_model_option, add_calc_file_option,
    add_deterministic_option, add_allow_charge_mult_mismatch_option,
)
from mlmm.cli.decorators import resolve_yaml_sources, load_merged_yaml_cfg, make_is_param_explicit, render_cli_exception


# Defaults (merge order: defaults ← YAML ← CLI)

# Geometry handling (Cartesian recommended for scans)
GEOM_KW: Dict[str, Any] = deepcopy(_OPT_GEOM_KW)

# ML/MM calculator defaults (shared with opt/path_*)
CALC_KW: Dict[str, Any] = deepcopy(_OPT_CALC_KW)

# Optimizer base (convergence, dumping, etc.)
OPT_BASE_KW: Dict[str, Any] = deepcopy(_OPT_BASE_KW)
OPT_BASE_KW.update({
    "out_dir": OUT_DIR_SCAN,
})

# LBFGS specifics
LBFGS_KW: Dict[str, Any] = deepcopy(_OPT_LBFGS_KW)
LBFGS_KW.update({
    "out_dir": OUT_DIR_SCAN,
})

# Bias (harmonic well) defaults; can be overridden via YAML: section "bias"
BIAS_KW: Dict[str, Any] = deepcopy(_BIAS_KW_DEFAULT)

# Bond-change detection (as in path_search)
BOND_KW: Dict[str, Any] = deepcopy(_BOND_KW_DEFAULT)


def _coords3d_to_xyz_string(geom, energy: Optional[float] = None) -> str:
    s = geom.as_xyz()
    lines = s.splitlines()
    if energy is not None and len(lines) >= 2 and lines[0].strip().isdigit():
        lines[1] = f"{energy:.12f}"
        s = "\n".join(lines)
    if not s.endswith("\n"):
        s += "\n"
    return s


def _pair_distances(coords_ang: np.ndarray, pairs: Iterable[Tuple[int, int]]) -> List[float]:
    """
    coords_ang: (N,3) in Å; returns a list of distances (Å) for the given pairs.
    """
    dists: List[float] = []
    for i, j in pairs:
        v = coords_ang[i] - coords_ang[j]
        d = float(np.linalg.norm(v))
        dists.append(d)
    return dists


def _schedule_for_stage(
    coords_ang: np.ndarray,
    tuples: List[Tuple[int, int, float]],
    max_step_size_ang: float,
) -> Tuple[int, List[float], List[float], List[float]]:
    """
    Given current *Å* coords and stage tuples, compute:
      N: number of steps
      r0: initial distances per tuple (Å)
      rT: target distances per tuple (Å)
      step_widths: δ_k per tuple (Å, signed)
    """
    pairs = [(i, j) for (i, j, _) in tuples]
    r0 = _pair_distances(coords_ang, pairs)
    rT = [t for (_, _, t) in tuples]
    deltas = [RT - R0 for (R0, RT) in zip(r0, rT)]
    d_max = max((abs(d) for d in deltas), default=0.0)
    if d_max <= 0.0:
        return 0, r0, rT, [0.0] * len(tuples)
    if max_step_size_ang <= 0.0:
        raise click.BadParameter("--max-step-size must be > 0.")
    N = int(math.ceil(d_max / max_step_size_ang))
    step_widths = [d / N for d in deltas]
    return N, r0, rT, step_widths



def _has_bond_change(x, y, bond_cfg: Dict[str, Any]) -> Tuple[bool, str]:
    """
    Return for covalent bonds forming/breaking between `x` and `y`.
    """
    res = compare_structures(
        x, y,
        device=bond_cfg.get("device", "cuda"),
        bond_factor=float(bond_cfg.get("bond_factor", 1.20)),
        margin_fraction=float(bond_cfg.get("margin_fraction", 0.05)),
        delta_fraction=float(bond_cfg.get("delta_fraction", 0.05)),
    )
    formed = len(getattr(res, "formed_covalent", [])) > 0
    broken = len(getattr(res, "broken_covalent", [])) > 0
    summary = summarize_changes(x, res, one_based=True)
    return (formed or broken), summary


def _snapshot_geometry(g) -> Any:
    """Create an independent pysisyphus Geometry snapshot from the given Geometry."""
    return snapshot_geometry(g, coord_type_default="cart")


@click.command(
    help="Bond-length driven scan with staged harmonic restraints and relaxation (ML/MM).",
    context_settings={
        "help_option_names": ["-h", "--help"],
        "ignore_unknown_options": True,
        "allow_extra_args": True,
    },
)
@click.option(
    "-i", "--input",
    "input_path",
    type=click.Path(path_type=Path, exists=True, dir_okay=False),
    required=True,
    help="Full-enzyme PDB used by the ML/MM calculator and as reference for conversions.",
)
@click.option(
    "--parm",
    "real_parm7",
    type=click.Path(path_type=Path, exists=True, dir_okay=False),
    required=True,
    help="Amber parm7 topology covering the entire enzyme complex.",
)
@click.option(
    "--model-pdb",
    type=click.Path(path_type=Path, exists=True, dir_okay=False),
    required=False,
    help="PDB defining the ML-region atoms for ML/MM. Optional when --detect-layer is enabled.",
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
    "--freeze-atoms",
    "freeze_atoms_cli",
    type=str,
    default=None,
    show_default=False,
    help="Comma-separated 1-based atom indices to freeze (e.g., '1,3,5').",
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
    "scan_lists_raw",
    type=str,
    multiple=True,
    required=False,
    help="Scan targets: inline Python literal (e.g. '[(1,5,1.4)]') or a YAML/JSON spec file path. "
         "Multiple inline literals define sequential stages.",
)
@click.option("--one-based/--zero-based", "one_based", default=True, show_default=True,
              help="Interpret (i,j) indices in --scan-lists as 1-based (default) or 0-based.")
@click.option(
    "--print-parsed/--no-print-parsed",
    "print_parsed",
    default=False,
    show_default=True,
    help="Print parsed scan targets after resolving -s/--scan-lists.",
)
@click.option("--max-step-size", type=float, default=0.20, show_default=True,
              help="Maximum change in any scanned bond length per step [Å].")
@click.option("--bias-k", type=float, default=None, show_default=False,
              help=(
                  "Harmonic well strength k [eV/Å^2]. "
                  "Defaults to YAML bias.k (BIAS_KW['k']=300 in defaults.py) when omitted; "
                  "explicit CLI value overrides YAML."
              ))
@click.option(
    "--opt-mode",
    type=click.Choice(["grad", "hess", "lbfgs", "rfo", "light", "heavy"], case_sensitive=False),
    default=None,
    show_default=False,
    help="Compatibility option for mlmm all forwarding. "
         "Scan relaxations always use LBFGS; values other than grad/lbfgs/light emit a warning.",
)
@click.option(
    "--max-cycles",
    type=int,
    default=10000,
    show_default=True,
    help="Maximum LBFGS cycles per biased step and per (pre|end)opt stage.",
)
@click.option(
    "--relax-max-cycles",
    type=int,
    default=None,
    show_default=False,
    help="Compatibility alias of --max-cycles (overrides it when provided).",
)
@click.option(
    "--dump/--no-dump",
    "dump",
    default=False,
    show_default=True,
    help="Write per-step optimizer trajectory files. "
         "scan_trj.xyz and scan.pdb are always written per-stage and as a combined file in out-dir, regardless of this flag.",
)
@click.option("-o", "--out-dir", type=str, default=OUT_DIR_SCAN, show_default=True,
              help="Base output directory.")
@click.option(
    "--thresh",
    type=click.Choice(THRESH_CHOICES, case_sensitive=False),
    default=None,
    help="Convergence preset for relaxations.",
)
@click.option(
    "--config",
    "config_yaml",
    type=click.Path(path_type=Path, exists=True, dir_okay=False),
    default=None,
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
    help="Pre-optimize initial structure without bias before the scan.",
)
@click.option(
    "--endopt/--no-endopt",
    "endopt",
    default=False,
    show_default=True,
    help="After each stage, run an additional unbiased optimization of the stage result.",
)
@click.option(
    "--dry-run/--no-dry-run",
    "dry_run",
    default=False,
    show_default=True,
    help="Validate options and print the execution plan without running the scan.",
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
@add_ml_layer_detection_options()
@add_ml_charge_spin_options()
@add_print_every_option()
@add_precision_option()
@add_backend_model_option()
@add_calc_file_option()
@add_deterministic_option()
@add_allow_charge_mult_mismatch_option()
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
    scan_lists_raw: Sequence[str],
    one_based: bool,
    print_parsed: bool,
    max_step_size: float,
    bias_k: Optional[float],
    opt_mode: Optional[str],
    max_cycles: int,
    relax_max_cycles: Optional[int],
    dump: bool,
    out_dir: str,
    thresh: Optional[str],
    config_yaml: Optional[Path],
    ref_pdb: Optional[Path],
    preopt: bool,
    endopt: bool,
    dry_run: bool,
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
    calc_file: Optional[str],
    calc_factory: str,
) -> None:
    _is_param_explicit = make_is_param_explicit(ctx)

    set_convert_file_enabled(convert_files)
    time_start = time.perf_counter()
    config_yaml, override_yaml, used_legacy_yaml = resolve_yaml_sources(
        config_yaml=config_yaml,
        override_yaml=None,
        args_yaml_legacy=None,
    )

    if relax_max_cycles is not None:
        max_cycles = int(relax_max_cycles)
    if max_cycles <= 0:
        raise click.BadParameter("--max-cycles must be > 0.")
    if opt_mode is not None and str(opt_mode).lower() not in {"lbfgs", "light", "grad"}:
        click.echo(
            f"[scan] NOTE: --opt-mode={opt_mode} is accepted for compatibility, "
            "but scan relaxations use LBFGS.",
            err=True,
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
            charge, spin = resolve_charge_spin_or_raise(
                prepared_input, charge, spin,
                ligand_charge=ligand_charge, prefix="[scan]",
            )

            try:
                freeze_atoms_list = _parse_freeze_atoms(freeze_atoms_cli)
            except click.BadParameter as e:
                click.echo(f"ERROR: {e}", err=True)
                sys.exit(1)

            model_indices: Optional[List[int]] = None
            if model_indices_str:
                try:
                    model_indices = parse_indices_string(model_indices_str, one_based=model_indices_one_based)
                except click.BadParameter as e:
                    click.echo(f"ERROR: {e}", err=True)
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
            bond_cfg = dict(BOND_KW)

            apply_yaml_overrides(
                yaml_cfg,
                [
                    (geom_cfg, (("geom",),)),
                    (calc_cfg, (("calc",), ("mlmm",))),
                    (opt_cfg, (("opt",),)),
                    (lbfgs_cfg, (("lbfgs",), ("opt", "lbfgs"))),
                    (bias_cfg, (("bias",),)),
                    (bond_cfg, (("bond",),)),
                ],
            )

            # Staged scans run restrained L-BFGS with no microiteration, so DLC
            # over the ML/MM system is meaningless (it crashes poly_line_search
            # with a Cartesian/internal dimension mismatch); force Cartesian,
            # matching path-opt / path-search.
            geom_cfg["coord_type"] = "cart"

            try:
                geom_freeze = _normalize_geom_freeze(geom_cfg.get("freeze_atoms"))
            except click.BadParameter as e:
                click.echo(f"ERROR: {e}", err=True)
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
            # only override the (already YAML-merged) max_cycles when the user
            # explicitly passed --max-cycles or --relax-max-cycles. Mirrors
            # opt.py's _is_param_explicit gating; an unconditional assignment
            # here silently clobbered the --config YAML tier.
            if _is_param_explicit("max_cycles") or relax_max_cycles is not None:
                opt_cfg["max_cycles"] = int(max_cycles)
                lbfgs_cfg["max_cycles"] = int(max_cycles)
            if thresh is not None:
                opt_cfg["thresh"] = str(thresh)

            if bias_k is not None:
                bias_cfg["k"] = float(bias_k)

            out_dir_path = Path(out_dir).resolve()

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
            # --calc-file overrides --backend with a user ASE Calculator (custom backend).
            from mlmm.backends import apply_calc_file_to_calc_cfg
            apply_calc_file_to_calc_cfg(calc_cfg, calc_file, calc_factory)
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
            echo_geom = format_freeze_atoms_for_echo(geom_cfg, key="freeze_atoms")
            echo_calc = format_freeze_atoms_for_echo(filter_calc_for_echo(calc_cfg), key="freeze_atoms")
            echo_opt = strip_inherited_keys({**opt_cfg, "out_dir": str(out_dir_path)}, OPT_BASE_KW, mode="same")
            # Show only lbfgs-specific settings, not inherited from opt_cfg
            echo_lbfgs = strip_inherited_keys(lbfgs_cfg, opt_cfg)
            click.echo(pretty_block("geom", echo_geom))
            click.echo(pretty_block("calc", echo_calc))
            click.echo(pretty_block("opt", echo_opt))
            click.echo(pretty_block("lbfgs", echo_lbfgs))
            click.echo(pretty_block("bias", bias_cfg))
            click.echo(pretty_block("bond", bond_cfg))

            pdb_atom_meta: List[Dict[str, Any]] = []
            if source_path.suffix.lower() == ".pdb":
                pdb_atom_meta = load_pdb_atom_metadata(source_path)

            cli_scan_values = collect_single_option_values(
                sys.argv[1:], ("-s", "--scan-lists"), "--scan-lists"
            )
            if not cli_scan_values:
                raise click.BadParameter("--scan-lists is required.")

            stages: List[List[Tuple[int, int, float]]]
            scan_one_based = bool(one_based)
            scan_source = "--scan-lists"
            # Bidirectional scan support (4-tuple): track which stages
            # need geometry snapshot/reset.
            _bidir_reset_before: set = set()
            _bidir_snapshot_before: set = set()
            # Auto-detect: single value that is a YAML/JSON file → spec mode
            if len(cli_scan_values) == 1 and is_scan_spec_file(cli_scan_values[0]):
                spec_path = Path(cli_scan_values[0])
                stages, scan_one_based = parse_scan_spec_stages(
                    spec_path,
                    one_based_default=one_based,
                    atom_meta=pdb_atom_meta,
                    option_name="--scan-lists",
                )
                scan_source = f"--scan-lists ({spec_path})"
            else:
                stages = []
                for idx, raw in enumerate(cli_scan_values, start=1):
                    parsed, _ = parse_scan_list_triples(
                        raw,
                        one_based=scan_one_based,
                        atom_meta=pdb_atom_meta,
                        option_name=f"--scan-lists #{idx}",
                    )
                    for t in parsed:
                        for dist in t[2:]:
                            if dist <= 0.0:
                                raise click.BadParameter(
                                    f"Non-positive target length in --scan-lists #{idx}: {t}."
                                )
                    # Expand 4-tuples into two stages with reset marker
                    has_4tuple = any(len(t) == 4 for t in parsed)
                    if has_4tuple:
                        for t in parsed:
                            if len(t) == 4:
                                i, j, start, end = t
                                stage_a_idx = len(stages)
                                stages.append([(i, j, start)])
                                _bidir_snapshot_before.add(stage_a_idx)
                                _bidir_reset_before.add(stage_a_idx + 1)
                                stages.append([(i, j, end)])
                            else:
                                stages.append([t])
                    else:
                        stages.append(parsed)
            K = len(stages)
            click.echo(f"[scan] Received {K} stage(s).", narrative=True)
            if print_parsed:
                click.echo(
                    pretty_block(
                        "scan-parsed",
                        {
                            "source": scan_source,
                            "one_based": bool(scan_one_based),
                            "stages_0based": stages,
                        },
                    )
                )
                # --print-parsed means "just show the parsed spec": exit
                # before any GPU calculation. (Also gives scan a GPU-free
                # spec-validation path.)
                sys.exit(0)

            if dry_run:
                model_region_source = "bfactor"
                if not detect_layer:
                    if model_pdb is not None:
                        model_region_source = "model_pdb"
                    elif model_indices:
                        model_region_source = "model_indices"
                click.echo(
                    pretty_block(
                        "dry_run_plan",
                        {
                            "input_geometry": str(geom_input_path),
                            "output_dir": str(out_dir_path),
                            "detect_layer": bool(detect_layer),
                            "model_region_source": model_region_source,
                            "num_stages": len(stages),
                            "stages_0based": stages,
                            "preopt": bool(preopt),
                            "endopt": bool(endopt),
                            "bias_k": float(bias_cfg["k"]),
                            "max_step_size": float(max_step_size),
                            "max_cycles": int(max_cycles),
                            "backend": calc_cfg.get("backend", "uma"),
                            "embedcharge": bool(calc_cfg.get("embedcharge", False)),
                        },
                    )
                )
                click.echo("[dry-run] Validation complete. Scan execution was skipped.")
                return

            if pdb_atom_meta:
                click.echo("[scan] PDB atom details for scanned pairs:", detail=True)
                legend = PDB_ATOM_META_HEADER
                click.echo(f"        legend: {legend}", detail=True)
                for stage_idx, tuples in enumerate(stages, start=1):
                    click.echo(f"  Stage {stage_idx}:", detail=True)
                    for pair_idx, (i, j, _) in enumerate(tuples, start=1):
                        click.echo(
                            f"    pair {pair_idx} i: {format_pdb_atom_metadata(pdb_atom_meta, i)}",
                            detail=True,
                        )
                        click.echo(
                            f"           j: {format_pdb_atom_metadata(pdb_atom_meta, j)}",
                            detail=True,
                        )
            stages_summary: List[Dict[str, Any]] = []

            out_dir_path.mkdir(parents=True, exist_ok=True)
            freeze = list(geom_cfg.get("freeze_atoms") or [])
            coord_type = geom_cfg.get("coord_type", GEOM_KW_DEFAULT["coord_type"])
            geom = geom_loader(
                geom_input_path,
                coord_type=coord_type,
                freeze_atoms=freeze,
            )
            if freeze:
                try:
                    geom.freeze_atoms = np.array(freeze, dtype=int)
                except Exception:
                    logger.debug("Failed to set freeze_atoms on geometry", exc_info=True)

            base_calc = mlmm(**calc_cfg)

            echo_resolved_device()

            max_step_bohr = float(max_step_size) * ANG2BOHR

            def _make_lbfgs(_out_dir: Path, _prefix: str) -> LBFGS:
                common = dict(opt_cfg)
                common["out_dir"] = str(_out_dir)
                common["prefix"] = _prefix
                args = {**lbfgs_cfg, **common}
                args["max_step"] = min(float(lbfgs_cfg.get("max_step", 0.30)), max_step_bohr)
                return LBFGS(geom, **args)

            if preopt:
                pre_dir = out_dir_path / "preopt"
                pre_dir.mkdir(parents=True, exist_ok=True)
                geom.set_calculator(base_calc)
                click.echo("[preopt] Unbiased relaxation (LBFGS) ...")
                optimizer0 = _make_lbfgs(pre_dir, "preopt")
                try:
                    optimizer0.run()
                except ZeroStepLength:
                    click.echo("[preopt] ZeroStepLength — continuing.", err=True)
                except OptimizationError as e:
                    click.echo(f"[preopt] OptimizationError — {e}", err=True)

                pre_xyz = pre_dir / "result.xyz"
                with open(pre_xyz, "w") as f:
                    f.write(_coords3d_to_xyz_string(geom))
                click.echo(f"[write] Wrote '{pre_xyz}'.")
                try:
                    convert_xyz_to_pdb(pre_xyz, source_path.resolve(), pre_dir / "result.pdb")
                    click.echo(f"[convert] Wrote '{pre_dir / 'result.pdb'}'.")
                except Exception as e:
                    click.echo(f"[convert] WARNING: Failed to convert preopt result to PDB: {e}", err=True)

            biased = HarmonicBiasCalculator(base_calc, k=float(bias_cfg["k"]))
            geom.set_calculator(biased)

            all_trj_blocks: List[str] = []
            # For bidirectional 4-tuple scans: save geometry before pass 1,
            # restore before pass 2, and reverse pass 1 trajectory.
            _bidir_saved_geom = None
            _bidir_pass1_trj: List[str] = []

            for k, tuples in enumerate(stages, start=1):
                # Bidirectional support: snapshot before pass 1
                stage_idx_0 = k - 1  # 0-based
                if stage_idx_0 in _bidir_snapshot_before:
                    _bidir_saved_geom = _snapshot_geometry(geom)
                    _bidir_pass1_trj = []
                # Bidirectional support: restore geometry before pass 2
                if stage_idx_0 in _bidir_reset_before and _bidir_saved_geom is not None:
                    click.echo("[bidir] Restoring initial geometry for reverse-direction pass.", narrative=True)
                    geom.coords = _bidir_saved_geom.coords.copy()

                stage_dir = out_dir_path / f"stage_{k:02d}"
                stage_dir.mkdir(parents=True, exist_ok=True)
                click.echo(f"\n--- Stage {k}/{K} ---")
                click.echo(f"Targets (i,j,target Å, 1-based): {[(i + 1, j + 1, t) for (i, j, t) in tuples]}", narrative=True)

                start_geom_for_stage = _snapshot_geometry(geom)

                R_bohr = np.array(geom.coords3d, dtype=float)
                R_ang = R_bohr * BOHR2ANG
                Nsteps, r0, rT, step_widths = _schedule_for_stage(R_ang, tuples, float(max_step_size))
                click.echo(f"[stage {k}] initial distances (Å) = {['{:.3f}'.format(x) for x in r0]}", narrative=True)
                click.echo(f"[stage {k}] target distances  (Å) = {['{:.3f}'.format(x) for x in rT]}", narrative=True)
                click.echo(f"[stage {k}] steps N = {Nsteps}", narrative=True)

                srec: Dict[str, Any] = {
                    "index": int(k),
                    "pairs_1based": [(int(i) + 1, int(j) + 1) for (i, j, _) in tuples],
                    "initial_distances_A": [float(f"{x:.3f}") for x in r0],
                    "target_distances_A": [float(f"{x:.3f}") for x in rT],
                    "per_pair_step_A": [float(f"{x:.3f}") for x in step_widths],
                    "num_steps": int(Nsteps),
                    "bond_change": {"changed": None, "summary": ""},
                }
                stages_summary.append(srec)

                trj_blocks: List[str] = []
                stage_energies: List[Optional[float]] = []
                stage_trj_path = stage_dir / "scan_trj.xyz"
                stage_trj_path.write_text("")
                pairs = [(i, j) for (i, j, _) in tuples]

                if Nsteps == 0:
                    if endopt:
                        geom.set_calculator(base_calc)
                        click.echo(f"[stage {k}] endopt (unbiased) ...", narrative=True)
                        end_optimizer = None
                        try:
                            end_optimizer = _make_lbfgs(stage_dir, "endopt")
                            end_optimizer.run()
                        except ZeroStepLength:
                            click.echo(f"[stage {k}] endopt ZeroStepLength — continuing.", err=True)
                        except OptimizationError as e:
                            click.echo(f"[stage {k}] endopt OptimizationError — {e}", err=True)
                        finally:
                            geom.set_calculator(biased)
                        srec["converged"] = getattr(end_optimizer, 'is_converged', None) if end_optimizer is not None else None

                    # No scan steps: empty energy trajectory
                    srec["energies_hartree"] = []

                    try:
                        changed, summary = _has_bond_change(start_geom_for_stage, geom, bond_cfg)
                        click.echo(f"[stage {k}] Covalent-bond changes (start vs final): {'Yes' if changed else 'No'}", narrative=True)
                        if changed and summary and summary.strip():
                            click.echo(textwrap.indent(summary.strip(), prefix="  "))
                        if not changed:
                            click.echo("  (no covalent changes detected)")
                        try:
                            srec["bond_change"]["changed"] = bool(changed)
                            srec["bond_change"]["summary"] = (summary.strip() if (summary and summary.strip()) else "")
                        except Exception:
                            logger.debug("Failed to store bond_change record", exc_info=True)
                    except Exception as e:
                        click.echo(f"[stage {k}] WARNING: Failed to evaluate bond changes: {e}", err=True)

                    final_xyz = stage_dir / "result.xyz"
                    with open(final_xyz, "w") as f:
                        f.write(_coords3d_to_xyz_string(geom))
                    click.echo(f"[write] Wrote '{final_xyz}'.")
                    # Capture final energy directly from geometry object
                    try:
                        srec["final_energy_hartree"] = float(geom.energy) if geom.energy is not None else None
                    except Exception:
                        srec["final_energy_hartree"] = None
                    try:
                        convert_xyz_to_pdb(final_xyz, source_path.resolve(), stage_dir / "result.pdb")
                        click.echo(f"[convert] Wrote '{stage_dir / 'result.pdb'}'.")
                    except Exception as e:
                        click.echo(f"[convert] WARNING: Failed to convert stage result to PDB: {e}", err=True)
                    continue

                for s in range(1, Nsteps + 1):
                    step_targets = [r0_i + s * dw for (r0_i, dw) in zip(r0, step_widths)]
                    biased.set_pairs([(i, j, t) for ((i, j), t) in zip(pairs, step_targets)])
                    geom.set_calculator(biased)

                    prefix = f"scan_s{s:04d}"
                    optimizer = _make_lbfgs(stage_dir, prefix)
                    click.echo(f"\n[stage {k}] step {s}/{Nsteps}: relaxation (LBFGS) ...", narrative=True)
                    try:
                        optimizer.run()
                    except ZeroStepLength:
                        click.echo(f"[stage {k}] step {s}: ZeroStepLength — continuing to next step.", err=True)
                    except OptimizationError as e:
                        click.echo(f"[stage {k}] step {s}: OptimizationError — {e}", err=True)

                    trj_blocks.append(_coords3d_to_xyz_string(geom))
                    stage_energies.append(float(geom.energy) if geom.energy is not None else None)
                    with open(stage_trj_path, "a") as _tf:
                        _tf.write(trj_blocks[-1])

                # Record convergence of the last optimizer (endopt if used, else last scan step)
                _last_opt = None
                if endopt:
                    geom.set_calculator(base_calc)
                    click.echo(f"[stage {k}] endopt (unbiased) ...", narrative=True)
                    try:
                        end_optimizer = _make_lbfgs(stage_dir, "endopt")
                        end_optimizer.run()
                    except ZeroStepLength:
                        click.echo(f"[stage {k}] endopt ZeroStepLength — continuing.", err=True)
                    except OptimizationError as e:
                        click.echo(f"[stage {k}] endopt OptimizationError — {e}", err=True)
                    finally:
                        geom.set_calculator(biased)
                    _last_opt = end_optimizer
                else:
                    _last_opt = optimizer
                srec["converged"] = getattr(_last_opt, 'is_converged', None) if _last_opt is not None else None

                # Store per-step energies in stage record
                srec["energies_hartree"] = stage_energies

                try:
                    changed, summary = _has_bond_change(start_geom_for_stage, geom, bond_cfg)
                    click.echo(f"[stage {k}] Covalent-bond changes (start vs final): {'Yes' if changed else 'No'}", narrative=True)
                    if changed and summary and summary.strip():
                        click.echo(textwrap.indent(summary.strip(), prefix="  "))
                    if not changed:
                        click.echo("  (no covalent changes detected)")
                    try:
                        srec["bond_change"]["changed"] = bool(changed)
                        srec["bond_change"]["summary"] = (summary.strip() if (summary and summary.strip()) else "")
                    except Exception:
                        logger.debug("Failed to store bond_change record", exc_info=True)
                except Exception as e:
                    click.echo(f"[stage {k}] WARNING: Failed to evaluate bond changes: {e}", err=True)

                if trj_blocks:
                    click.echo(f"[write] Wrote '{stage_trj_path}'.")
                    # Bidirectional trajectory assembly:
                    # pass 1 (initial→start) is saved; pass 2 (initial→end)
                    # triggers assembly: reversed(pass1) + pass2 → start→initial→end
                    if stage_idx_0 in _bidir_snapshot_before:
                        _bidir_pass1_trj = list(trj_blocks)
                    elif stage_idx_0 in _bidir_reset_before:
                        all_trj_blocks.extend(reversed(_bidir_pass1_trj))
                        all_trj_blocks.extend(trj_blocks)
                        _bidir_pass1_trj = []
                    else:
                        all_trj_blocks.extend(trj_blocks)
                    try:
                        convert_xyz_to_pdb(stage_trj_path, source_path.resolve(), stage_dir / "scan.pdb")
                        click.echo(f"[convert] Wrote '{stage_dir / 'scan.pdb'}'.")
                    except Exception as e:
                        click.echo(f"[convert] WARNING: Failed to convert stage trajectory to PDB: {e}", err=True)

                final_xyz = stage_dir / "result.xyz"
                with open(final_xyz, "w") as f:
                    f.write(_coords3d_to_xyz_string(geom))
                click.echo(f"[write] Wrote '{final_xyz}'.")
                # Capture final energy directly from geometry object
                try:
                    srec["final_energy_hartree"] = float(geom.energy) if geom.energy is not None else None
                except Exception:
                    srec["final_energy_hartree"] = None
                try:
                    convert_xyz_to_pdb(final_xyz, source_path.resolve(), stage_dir / "result.pdb")
                    click.echo(f"[convert] Wrote '{stage_dir / 'result.pdb'}'.")
                except Exception as e:
                    click.echo(f"[convert] WARNING: Failed to convert stage result to PDB: {e}", err=True)

        # Write combined scan_trj.xyz + scan.pdb to out_dir
        if all_trj_blocks:
            combined_trj = out_dir_path / "scan_trj.xyz"
            with open(combined_trj, "w") as f:
                f.write("".join(all_trj_blocks))
            click.echo(f"[write] Wrote '{combined_trj}'.")
            try:
                convert_xyz_to_pdb(combined_trj, source_path.resolve(), out_dir_path / "scan.pdb")
                click.echo(f"[convert] Wrote '{out_dir_path / 'scan.pdb'}'.")
            except Exception as e:
                click.echo(f"[convert] WARNING: Failed to convert combined trajectory to PDB: {e}", err=True)

        def _echo_human_summary(_stages: List[Dict[str, Any]], _max_step_size: float) -> None:
            """
            Print a readable end-of-run summary like the requested example.
            """
            def _fmt_target_value(x: float) -> str:
                # 2.600 -> "2.6", 1.500 -> "1.5"
                s = f"{x:.3f}".rstrip("0").rstrip(".")
                return s

            def _targets_triplet_str(pairs_1based: List[Tuple[int, int]], targets: List[float]) -> str:
                triples = [f"({i}, {j}, {_fmt_target_value(t)})" for (i, j), t in zip(pairs_1based, targets)]
                return "[" + ", ".join(triples) + "]"

            def _list_of_str_3f(values: List[float]) -> str:
                return "[" + ", ".join(f"'{v:.3f}'" for v in values) + "]"

            click.echo("\nSummary")
            click.echo("------------------")
            for s in _stages:
                idx = int(s.get("index", 0))
                pairs_1b = list(s.get("pairs_1based", []))
                r0 = list(s.get("initial_distances_A", []))
                rT = list(s.get("target_distances_A", []))
                dA = list(s.get("per_pair_step_A", []))
                N = int(s.get("num_steps", 0))
                bchg = s.get("bond_change", {}) or {}
                changed = bool(bchg.get("changed"))
                summary_txt = (bchg.get("summary") or "").strip()

                click.echo(f"[stage {idx}] Targets (i,j,target Å, 1-based): { _targets_triplet_str(pairs_1b, rT) }", narrative=True)
                click.echo(f"[stage {idx}] initial distances (Å) = { _list_of_str_3f(r0) }", narrative=True)
                click.echo(f"[stage {idx}] target distances  (Å) = { _list_of_str_3f(rT) }", narrative=True)
                click.echo(f"[stage {idx}] per_pair_step     (Å) = { _list_of_str_3f(dA) }", narrative=True)
                click.echo(f"[stage {idx}] steps N = {N}", narrative=True)
                click.echo(f"[stage {idx}] Covalent-bond changes (start vs final): {'Yes' if changed else 'No'}", narrative=True)
                if changed and summary_txt:
                    click.echo(textwrap.indent(summary_txt, prefix="  "))
                if not changed:
                    click.echo("  (no covalent changes detected)")
                click.echo("")  # blank line between stages

        _echo_human_summary(stages_summary, float(max_step_size))

        click.echo("\n====== Scan finished ======\n", narrative=True)
        click.echo(format_elapsed("[time] Elapsed Time for Scan", time_start), narrative=True)

        if out_json:
            from mlmm.core.utils import write_result_json
            json_stages = []
            for srec in stages_summary:
                stage_entry: Dict[str, Any] = {
                    "index": srec["index"],
                    "n_steps": srec["num_steps"],
                    "converged": srec.get("converged"),
                    "bond_changes": srec.get("bond_change", {}),
                    "pairs_1based": srec.get("pairs_1based"),
                    "initial_distances_angstrom": srec.get("initial_distances_A"),
                    "target_distances_angstrom": srec.get("target_distances_A"),
                }
                # Final energy: use value captured directly from geometry object
                stage_entry["final_energy_hartree"] = srec.get("final_energy_hartree")
                # Per-step energy trajectory
                stage_entry["energies_hartree"] = srec.get("energies_hartree", [])
                json_stages.append(stage_entry)
            result_data: Dict[str, Any] = {
                "status": "completed",
                "charge": calc_cfg.get("model_charge"),
                "spin": calc_cfg.get("model_mult"),
                "backend": calc_cfg.get("backend", "uma"),
                "max_step_size_angstrom": float(max_step_size),
                "n_stages": len(stages_summary),
                "stages": json_stages,
                "files": {
                    "scan_trj_xyz": "scan_trj.xyz",
                },
            }
            for ext in (".pdb",):
                f = out_dir_path / f"scan{ext}"
                if f.exists():
                    result_data["files"][f"scan_{ext[1:]}"] = f.name
            write_result_json(
                out_dir_path, result_data,
                command="scan",
                elapsed_seconds=time.perf_counter() - time_start,
            )

    except KeyboardInterrupt:
        click.echo("\nInterrupted by user.", err=True)
        sys.exit(130)
    except Exception as e:
        render_cli_exception(e, label="scan", out_dir=out_dir, command="scan", time_start=time_start)
    finally:
        # Release GPU memory so subsequent pipeline stages don't OOM.
        # `= None` decref's the heavy refs; `del` then removes names from
        # the local frame so torch.nn.Module hooks / closures cannot retain.
        base_calc = biased = geom = optimizer = optimizer0 = end_optimizer = None
        del base_calc, biased, geom, optimizer, optimizer0, end_optimizer
        gc.collect()  # break cyclic refs inside torch.nn.Module
        if torch.cuda.is_available():
            torch.cuda.empty_cache()


if __name__ == "__main__":
    cli()
