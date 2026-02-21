# mlmm_toolkit/scan.py

"""
ML/MM staged bond-length scan with harmonic restraints.

Example:
    mlmm scan -i pocket.pdb --real-parm7 real.parm7 --model-pdb ml_region.pdb -q 0 --scan-lists "[(12,45,2.20)]"

For detailed documentation, see: docs/scan.md
"""

from __future__ import annotations

from copy import deepcopy
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional, Sequence, Tuple

import ast
import math
import sys
import textwrap
import traceback

import click
import numpy as np
import time

from pysisyphus.helpers import geom_loader
from pysisyphus.optimizers.LBFGS import LBFGS
from pysisyphus.optimizers.exceptions import OptimizationError, ZeroStepLength
from pysisyphus.constants import BOHR2ANG, ANG2BOHR

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
    convert_xyz_to_pdb,
    load_yaml_dict,
    apply_yaml_overrides,
    pretty_block,
    strip_inherited_keys,
    format_freeze_atoms_for_echo,
    format_elapsed,
    merge_freeze_atom_indices,
    prepare_input_structure,
    resolve_charge_spin_or_raise,
    collect_single_option_values,
    load_pdb_atom_metadata,
    parse_scan_list_triples,
    parse_scan_spec_stages,
    PDB_ATOM_META_HEADER,
    format_pdb_atom_metadata,
    parse_indices_string,
    build_model_pdb_from_bfactors,
    build_model_pdb_from_indices,
    snapshot_geometry,
)
from .bond_changes import compare_structures, summarize_changes


# --------------------------------------------------------------------------------------
# Defaults (merge order: defaults ← YAML ← CLI)
# --------------------------------------------------------------------------------------

# Geometry handling (Cartesian recommended for scans)
GEOM_KW: Dict[str, Any] = deepcopy(_OPT_GEOM_KW)

# ML/MM calculator defaults (shared with opt/path_*)
CALC_KW: Dict[str, Any] = deepcopy(_OPT_CALC_KW)

# Optimizer base (convergence, dumping, etc.)
OPT_BASE_KW: Dict[str, Any] = deepcopy(_OPT_BASE_KW)
OPT_BASE_KW.update({
    "max_cycles": 100,
    "out_dir": "./result_scan/",
})

# LBFGS specifics
LBFGS_KW: Dict[str, Any] = deepcopy(_OPT_LBFGS_KW)
LBFGS_KW.update({
    "out_dir": "./result_scan/",
})

# Bias (harmonic well) defaults; can be overridden via YAML: section "bias"
BIAS_KW: Dict[str, Any] = {
    "k": 100,  # float, harmonic bias strength in eV/Å^2
}

# Bond-change detection (as in path_search)
BOND_KW: Dict[str, Any] = {
    "device": "cuda",            # str, device used during UMA graph analysis in bond detection
    "bond_factor": 1.20,         # float, scaling of covalent radii for bond cutoff
    "margin_fraction": 0.05,     # float, fractional margin to tolerate small deviations
    "delta_fraction": 0.05,      # float, change threshold to flag bond formation/breaking
}


def _coords3d_to_xyz_string(geom, energy: Optional[float] = None) -> str:
    s = geom.as_xyz()
    lines = s.splitlines()
    if energy is not None and len(lines) >= 2 and lines[0].strip().isdigit():
        lines[1] = f"{energy:.12f}"
        s = "\n".join(lines)
    if not s.endswith("\n"):
        s += "\n"
    return s


def _parse_scan_lists(args: Sequence[str], one_based: bool) -> List[List[Tuple[int, int, float]]]:
    """
    Parse multiple Python-like list strings:
      ["[(0,1,1.5), (2,3,2.0)]", "[(5,7,1.2)]", ...]
    Returns: [[(i,j,t), ...], [(i,j,t), ...], ...] with 0-based indices.
    """
    if not args:
        raise click.BadParameter("--scan-lists must be provided at least once.")
    stages: List[List[Tuple[int, int, float]]] = []
    for idx, s in enumerate(args, start=1):
        try:
            obj = ast.literal_eval(s)
        except Exception as e:
            raise click.BadParameter(f"Invalid literal for --scan-lists #{idx}: {e}")
        if not isinstance(obj, (list, tuple)):
            raise click.BadParameter(f"--scan-lists #{idx} must be a list/tuple of (i,j,target).")
        tuples: List[Tuple[int, int, float]] = []
        for t in obj:
            if (
                isinstance(t, (list, tuple)) and len(t) == 3
                and isinstance(t[0], (int, np.integer))
                and isinstance(t[1], (int, np.integer))
                and isinstance(t[2], (int, float, np.floating))
            ):
                i, j, r = int(t[0]), int(t[1]), float(t[2])
                if one_based:
                    i -= 1
                    j -= 1
                if i < 0 or j < 0:
                    raise click.BadParameter(f"Negative atom index in --scan-lists #{idx}: {(i,j,r)} (0-based expected).")
                if r <= 0.0:
                    raise click.BadParameter(f"Non-positive target length in --scan-lists #{idx}: {(i,j,r)}.")
                tuples.append((i, j, r))
            else:
                raise click.BadParameter(f"--scan-lists #{idx} contains an invalid triple: {t}")
        stages.append(tuples)
    return stages


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


# --------------------------------------------------------------------------------------
# Bond‑change helpers
# --------------------------------------------------------------------------------------

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
    help="Bond-length driven scan with staged harmonic restraints and relaxation (UMA only).",
    context_settings={"help_option_names": ["-h", "--help"], "allow_extra_args": True},
)
@click.option(
    "-i", "--input",
    "input_path",
    type=click.Path(path_type=Path, exists=True, dir_okay=False),
    required=True,
    help="Full-enzyme PDB used by the ML/MM calculator and as reference for conversions.",
)
@click.option(
    "--real-parm7",
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
@click.option("-q", "--charge", type=int, default=None, show_default=False, help="ML region charge.")
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
    "--scan-lists",
    "--scan-list",
    "scan_lists_raw",
    type=str,
    multiple=True,
    required=False,
    help="Python-like list of (i,j,target) per stage. Pass a single --scan-list(s) followed by "
         "multiple literals to run sequential stages, e.g. --scan-lists '[(0,1,1.50)]' '[(5,7,1.20)]'.",
)
@click.option(
    "--spec",
    "spec_path",
    type=click.Path(path_type=Path, exists=True, dir_okay=False),
    required=False,
    help="YAML/JSON scan spec file (recommended). Use this instead of --scan-list(s).",
)
@click.option("--one-based/--zero-based", "one_based", default=True, show_default=True,
              help="Interpret (i,j) indices in --scan-list(s) as 1-based (default) or 0-based.")
@click.option(
    "--print-parsed/--no-print-parsed",
    "print_parsed",
    default=False,
    show_default=True,
    help="Print parsed scan targets after resolving --spec/--scan-list(s).",
)
@click.option("--max-step-size", type=float, default=0.20, show_default=True,
              help="Maximum change in any scanned bond length per step [Å].")
@click.option("--bias-k", type=float, default=100, show_default=True,
              help="Harmonic well strength k [eV/Å^2].")
@click.option(
    "--max-cycles",
    type=int,
    default=10000,
    show_default=True,
    help="Maximum LBFGS cycles per biased step and per (pre|end)opt stage.",
)
@click.option("--dump", type=click.BOOL, default=False, show_default=True,
              help="Write stage trajectory as scan.trj (and scan.pdb for PDB input).")
@click.option("--out-dir", type=str, default="./result_scan/", show_default=True,
              help="Base output directory.")
@click.option(
    "--thresh",
    type=str,
    default=None,
    help="Convergence preset for relaxations (gau_loose|gau|gau_tight|gau_vtight|baker|never).",
)
@click.option(
    "--args-yaml",
    type=click.Path(path_type=Path, exists=True, dir_okay=False),
    default=None,
    help="YAML file with extra args (sections: geom, calc/mlmm, opt, lbfgs, bias, bond).",
)
@click.option(
    "--ref-pdb",
    type=click.Path(path_type=Path, exists=True, dir_okay=False),
    default=None,
    help="Reference PDB topology to use when --input is XYZ (keeps XYZ coordinates).",
)
@click.option("--preopt", type=click.BOOL, default=True, show_default=True,
              help="Pre-optimize initial structure without bias before the scan.")
@click.option("--endopt", type=click.BOOL, default=True, show_default=True,
              help="After each stage, run an additional unbiased optimization of the stage result.")
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
    scan_lists_raw: Sequence[str],
    spec_path: Optional[Path],
    one_based: bool,
    print_parsed: bool,
    max_step_size: float,
    bias_k: Optional[float],
    max_cycles: int,
    dump: bool,
    out_dir: str,
    thresh: Optional[str],
    args_yaml: Optional[Path],
    ref_pdb: Optional[Path],
    preopt: bool,
    endopt: bool,
) -> None:
    time_start = time.perf_counter()

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

            yaml_cfg = load_yaml_dict(args_yaml)

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

            try:
                geom_freeze = _normalize_geom_freeze(geom_cfg.get("freeze_atoms"))
            except click.BadParameter as e:
                click.echo(f"ERROR: {e}", err=True)
                sys.exit(1)
            geom_cfg["freeze_atoms"] = geom_freeze
            if freeze_atoms_list:
                merge_freeze_atom_indices(geom_cfg, freeze_atoms_list)
            freeze_atoms_final = list(geom_cfg.get("freeze_atoms") or [])
            calc_cfg["freeze_atoms"] = freeze_atoms_final

            opt_cfg["out_dir"] = out_dir
            opt_cfg["dump"] = False
            opt_cfg["max_cycles"] = int(max_cycles)
            if thresh is not None:
                opt_cfg["thresh"] = str(thresh)
            lbfgs_cfg["max_cycles"] = int(max_cycles)

            if bias_k is not None:
                bias_cfg["k"] = float(bias_k)

            out_dir_path = Path(out_dir).resolve()

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
            echo_geom = format_freeze_atoms_for_echo(geom_cfg, key="freeze_atoms")
            echo_calc = strip_inherited_keys(calc_cfg, CALC_KW, mode="same")
            echo_calc = format_freeze_atoms_for_echo(echo_calc, key="freeze_atoms")
            echo_opt = strip_inherited_keys({**opt_cfg, "out_dir": str(out_dir_path)}, OPT_BASE_KW, mode="same")
            # Show only lbfgs-specific settings, not inherited from opt_cfg
            echo_lbfgs = strip_inherited_keys(lbfgs_cfg, opt_cfg)
            echo_bias = strip_inherited_keys(bias_cfg, BIAS_KW, mode="same")
            echo_bond = strip_inherited_keys(bond_cfg, BOND_KW, mode="same")
            click.echo(pretty_block("geom", echo_geom))
            click.echo(pretty_block("calc", echo_calc))
            click.echo(pretty_block("opt", echo_opt))
            click.echo(pretty_block("lbfgs", echo_lbfgs))
            click.echo(pretty_block("bias", echo_bias))
            click.echo(pretty_block("bond", echo_bond))

            pdb_atom_meta: List[Dict[str, Any]] = []
            if source_path.suffix.lower() == ".pdb":
                pdb_atom_meta = load_pdb_atom_metadata(source_path)

            cli_scan_values = collect_single_option_values(
                sys.argv[1:], ("--scan-lists", "--scan-list"), "--scan-list/--scan-lists"
            )
            if spec_path is not None and cli_scan_values:
                raise click.BadParameter("Use either --spec or --scan-list(s), not both.")

            stages: List[List[Tuple[int, int, float]]]
            scan_one_based = bool(one_based)
            scan_source = "--scan-list(s)"
            if spec_path is not None:
                stages, scan_one_based = parse_scan_spec_stages(
                    spec_path,
                    one_based_default=one_based,
                    atom_meta=pdb_atom_meta,
                    option_name="--spec",
                )
                scan_source = f"--spec ({spec_path})"
            else:
                if not cli_scan_values:
                    raise click.BadParameter("Provide either --spec or --scan-list(s).")
                stages = []
                for idx, raw in enumerate(cli_scan_values, start=1):
                    parsed, _ = parse_scan_list_triples(
                        raw,
                        one_based=scan_one_based,
                        atom_meta=pdb_atom_meta,
                        option_name=f"--scan-lists #{idx}",
                    )
                    for i, j, r in parsed:
                        if r <= 0.0:
                            raise click.BadParameter(
                                f"Non-positive target length in --scan-lists #{idx}: {(i, j, r)}."
                            )
                    stages.append(parsed)
            K = len(stages)
            click.echo(f"[scan] Received {K} stage(s).")
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

            if pdb_atom_meta:
                click.echo("[scan] PDB atom details for scanned pairs:")
                legend = PDB_ATOM_META_HEADER
                click.echo(f"        legend: {legend}")
                for stage_idx, tuples in enumerate(stages, start=1):
                    click.echo(f"  Stage {stage_idx}:")
                    for pair_idx, (i, j, _) in enumerate(tuples, start=1):
                        click.echo(
                            f"    pair {pair_idx} i: {format_pdb_atom_metadata(pdb_atom_meta, i)}"
                        )
                        click.echo(
                            f"           j: {format_pdb_atom_metadata(pdb_atom_meta, j)}"
                        )
            stages_summary: List[Dict[str, Any]] = []

            out_dir_path.mkdir(parents=True, exist_ok=True)
            coord_type = geom_cfg.get("coord_type", "cart")
            geom = geom_loader(geom_input_path, coord_type=coord_type)

            freeze = list(geom_cfg.get("freeze_atoms") or [])
            if freeze:
                try:
                    geom.freeze_atoms = np.array(freeze, dtype=int)
                except Exception:
                    pass

            base_calc = mlmm(**calc_cfg)

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
                optimizer0 = _make_lbfgs(pre_dir, "preopt_")
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

            for k, tuples in enumerate(stages, start=1):
                stage_dir = out_dir_path / f"stage_{k:02d}"
                stage_dir.mkdir(parents=True, exist_ok=True)
                click.echo(f"\n--- Stage {k}/{K} ---")
                click.echo(f"Targets (i,j,target Å): {tuples}")

                start_geom_for_stage = _snapshot_geometry(geom)

                R_bohr = np.array(geom.coords3d, dtype=float)
                R_ang = R_bohr * BOHR2ANG
                Nsteps, r0, rT, step_widths = _schedule_for_stage(R_ang, tuples, float(max_step_size))
                click.echo(f"[stage {k}] initial distances (Å) = {['{:.3f}'.format(x) for x in r0]}")
                click.echo(f"[stage {k}] target distances  (Å) = {['{:.3f}'.format(x) for x in rT]}")
                click.echo(f"[stage {k}] steps N = {Nsteps}")

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

                trj_blocks: List[str] = [] if dump else None
                pairs = [(i, j) for (i, j, _) in tuples]

                if Nsteps == 0:
                    if endopt:
                        geom.set_calculator(base_calc)
                        click.echo(f"[stage {k}] endopt (unbiased) ...")
                        try:
                            end_optimizer = _make_lbfgs(stage_dir, "endopt_")
                            end_optimizer.run()
                        except ZeroStepLength:
                            click.echo(f"[stage {k}] endopt ZeroStepLength — continuing.", err=True)
                        except OptimizationError as e:
                            click.echo(f"[stage {k}] endopt OptimizationError — {e}", err=True)
                        finally:
                            geom.set_calculator(biased)

                    try:
                        changed, summary = _has_bond_change(start_geom_for_stage, geom, bond_cfg)
                        click.echo(f"[stage {k}] Covalent-bond changes (start vs final): {'Yes' if changed else 'No'}")
                        if changed and summary and summary.strip():
                            click.echo(textwrap.indent(summary.strip(), prefix="  "))
                        if not changed:
                            click.echo("  (no covalent changes detected)")
                        try:
                            srec["bond_change"]["changed"] = bool(changed)
                            srec["bond_change"]["summary"] = (summary.strip() if (summary and summary.strip()) else "")
                        except Exception:
                            pass
                    except Exception as e:
                        click.echo(f"[stage {k}] WARNING: Failed to evaluate bond changes: {e}", err=True)

                    final_xyz = stage_dir / "result.xyz"
                    with open(final_xyz, "w") as f:
                        f.write(_coords3d_to_xyz_string(geom))
                    click.echo(f"[write] Wrote '{final_xyz}'.")
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

                    prefix = f"scan_s{s:04d}_"
                    optimizer = _make_lbfgs(stage_dir, prefix)
                    click.echo(f"[stage {k}] step {s}/{Nsteps}: relaxation (LBFGS) ...")
                    try:
                        optimizer.run()
                    except ZeroStepLength:
                        click.echo(f"[stage {k}] step {s}: ZeroStepLength — continuing to next step.", err=True)
                    except OptimizationError as e:
                        click.echo(f"[stage {k}] step {s}: OptimizationError — {e}", err=True)

                    if dump and trj_blocks is not None:
                        trj_blocks.append(_coords3d_to_xyz_string(geom))

                if endopt:
                    geom.set_calculator(base_calc)
                    click.echo(f"[stage {k}] endopt (unbiased) ...")
                    try:
                        end_optimizer = _make_lbfgs(stage_dir, "endopt_")
                        end_optimizer.run()
                    except ZeroStepLength:
                        click.echo(f"[stage {k}] endopt ZeroStepLength — continuing.", err=True)
                    except OptimizationError as e:
                        click.echo(f"[stage {k}] endopt OptimizationError — {e}", err=True)
                    finally:
                        geom.set_calculator(biased)

                try:
                    changed, summary = _has_bond_change(start_geom_for_stage, geom, bond_cfg)
                    click.echo(f"[stage {k}] Covalent-bond changes (start vs final): {'Yes' if changed else 'No'}")
                    if changed and summary and summary.strip():
                        click.echo(textwrap.indent(summary.strip(), prefix="  "))
                    if not changed:
                        click.echo("  (no covalent changes detected)")
                    try:
                        srec["bond_change"]["changed"] = bool(changed)
                        srec["bond_change"]["summary"] = (summary.strip() if (summary and summary.strip()) else "")
                    except Exception:
                        pass
                except Exception as e:
                    click.echo(f"[stage {k}] WARNING: Failed to evaluate bond changes: {e}", err=True)

                if dump and trj_blocks:
                    trj_path = stage_dir / "scan.trj"
                    with open(trj_path, "w") as f:
                        f.write("".join(trj_blocks))
                    click.echo(f"[write] Wrote '{trj_path}'.")
                    try:
                        convert_xyz_to_pdb(trj_path, source_path.resolve(), stage_dir / "scan.pdb")
                        click.echo(f"[convert] Wrote '{stage_dir / 'scan.pdb'}'.")
                    except Exception as e:
                        click.echo(f"[convert] WARNING: Failed to convert stage trajectory to PDB: {e}", err=True)

                final_xyz = stage_dir / "result.xyz"
                with open(final_xyz, "w") as f:
                    f.write(_coords3d_to_xyz_string(geom))
                click.echo(f"[write] Wrote '{final_xyz}'.")
                try:
                    convert_xyz_to_pdb(final_xyz, source_path.resolve(), stage_dir / "result.pdb")
                    click.echo(f"[convert] Wrote '{stage_dir / 'result.pdb'}'.")
                except Exception as e:
                    click.echo(f"[convert] WARNING: Failed to convert stage result to PDB: {e}", err=True)

        # ------------------------------------------------------------------
        # 5) Final summary echo (human‑friendly)
        # ------------------------------------------------------------------
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

                click.echo(f"[stage {idx}] Targets (i,j,target Å): { _targets_triplet_str(pairs_1b, rT) }")
                click.echo(f"[stage {idx}] initial distances (Å) = { _list_of_str_3f(r0) }")
                click.echo(f"[stage {idx}] target distances  (Å) = { _list_of_str_3f(rT) }")
                click.echo(f"[stage {idx}] per_pair_step     (Å) = { _list_of_str_3f(dA) }")
                click.echo(f"[stage {idx}] steps N = {N}")
                click.echo(f"[stage {idx}] Covalent-bond changes (start vs final): {'Yes' if changed else 'No'}")
                if changed and summary_txt:
                    click.echo(textwrap.indent(summary_txt, prefix="  "))
                if not changed:
                    click.echo("  (no covalent changes detected)")
                click.echo("")  # blank line between stages

        _echo_human_summary(stages_summary, float(max_step_size))
        # ------------------------------------------------------------------

        click.echo("\n=== Scan finished ===\n")

        click.echo(format_elapsed("[time] Elapsed Time for Scan", time_start))

    except KeyboardInterrupt:
        click.echo("\nInterrupted by user.", err=True)
        sys.exit(130)
    except Exception as e:
        tb = "".join(traceback.format_exception(type(e), e, e.__traceback__))
        click.echo("Unhandled error during scan:\n" + textwrap.indent(tb, "  "), err=True)
        sys.exit(1)


if __name__ == "__main__":
    cli()
