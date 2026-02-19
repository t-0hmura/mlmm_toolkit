# mlmm_toolkit/irc.py

"""
ML/MM IRC calculation using the EulerPC predictor-corrector integrator.

Example:
    mlmm irc -i ts.pdb --real-parm7 real.parm7 --model-pdb ml_region.pdb -q 0

For detailed documentation, see: docs/irc.md
"""

from __future__ import annotations

from pathlib import Path
from typing import Any, Dict, Optional, List

import sys
import textwrap

import click
import yaml
import time
import numpy as np

from pysisyphus.helpers import geom_loader
from pysisyphus.irc.EulerPC import EulerPC
from .mlmm_calc import mlmm
from .freq import _torch_device, _calc_full_hessian_torch
from .defaults import (
    GEOM_KW_DEFAULT,
    MLMM_CALC_KW as _UMA_CALC_KW,
    IRC_KW,
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
    parse_indices_string,
    build_model_pdb_from_bfactors,
    build_model_pdb_from_indices,
)


# --------------------------
# Default configuration
# --------------------------

CALC_KW_DEFAULT: Dict[str, Any] = dict(_UMA_CALC_KW)

IRC_KW_DEFAULT: Dict[str, Any] = {
    # Arguments for IRC.__init__ (forwarded to EulerPC via **kwargs)
    "step_length": 0.10,         # float, default step length in mass-weighted coordinates (overridden by CLI)
    "max_cycles": 125,           # int, maximum IRC steps (overridden by CLI)
    "downhill": False,           # bool, follow downhill potential (debug option)
    "forward": True,             # bool, integrate forward branch (CLI override --forward)
    "backward": True,            # bool, integrate backward branch (CLI override --backward)
    "root": 0,                   # int, imaginary mode index for initial displacement (CLI override --root)
    "hessian_init": "calc",      # str, initial Hessian source ("calc" = calculator-provided TS Hessian)
    "displ": "energy",          # str, displacement metric (energy|length)
    "displ_energy": 1.0e-3,      # float, energy step (Hartree) when displ == "energy"
    "displ_length": 0.10,        # float, length step in mass-weighted coordinates when displ == "length"
    "rms_grad_thresh": 1.0e-3,   # float, RMS gradient threshold for convergence (Hartree/bohr)
    "hard_rms_grad_thresh": None,# Optional[float], stricter RMS gradient cutoff
    "energy_thresh": 1.0e-6,     # float, energy-change threshold for convergence (Hartree)
    "imag_below": 0.0,           # float, treat imaginary frequency below this as zero
    "force_inflection": True,    # bool, stop when force inflection detected
    "check_bonds": False,        # bool, enable bond-change detection during IRC
    "out_dir": "./result_irc/",  # str, output directory
    "prefix": "",                # str, file name prefix
    "dump_fn": "irc_data.h5",    # str, HDF5 dump filename
    "dump_every": 5,             # int, write dump every N steps

    # EulerPC-specific options
    "hessian_update": "bofill",  # str, Hessian update algorithm
    "hessian_recalc": None,      # Optional[int], force Hessian recalculation every N steps
    "max_pred_steps": 500,       # int, predictor steps per segment
    "loose_cycles": 3,           # int, cycles using looser thresholds
    "corr_func": "mbs",          # str, correction function selection
}


def _echo_convert_trj_to_pdb_if_exists(trj_path: Path, ref_pdb: Path, out_path: Path) -> None:
    if trj_path.exists():
        try:
            convert_xyz_to_pdb(trj_path, ref_pdb, out_path)
            click.echo(f"[convert] Wrote '{out_path}'.")
        except Exception as e:
            click.echo(f"[convert] WARNING: Failed to convert '{trj_path.name}' to PDB: {e}", err=True)


# --------------------------
# CLI
# --------------------------

@click.command(
    help="Run an IRC calculation with EulerPC. Only the documented CLI options are accepted; all other settings come from YAML.",
    context_settings={"help_option_names": ["-h", "--help"]},
)
@click.option(
    "-i", "--input",
    "input_path",
    type=click.Path(path_type=Path, exists=True, dir_okay=False),
    required=True,
    help="Input structure file (.pdb, .xyz, .trj, etc.).",
)
@click.option(
    "--real-parm7",
    "real_parm7",
    type=click.Path(path_type=Path, exists=True, dir_okay=False),
    required=False,
    help="Amber parm7 topology for the whole enzyme (MM region). "
         "If omitted, must be provided in YAML as calc.real_parm7.",
)
@click.option(
    "--model-pdb",
    "model_pdb",
    type=click.Path(path_type=Path, exists=True, dir_okay=False),
    required=False,
    help="PDB defining atoms belonging to the ML region. Optional when --detect-layer is enabled.",
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
@click.option("-q", "--charge", type=int, default=None, show_default=False, help="Total charge; overrides calc.charge from YAML.")
@click.option(
    "-m",
    "--multiplicity",
    "spin",
    type=int,
    default=None,
    show_default=False,
    help="Spin multiplicity (2S+1); overrides calc.spin from YAML.",
)
@click.option("--max-cycles", type=int, default=None, help="Maximum number of IRC steps; overrides irc.max_cycles from YAML.")
@click.option("--step-size", type=float, default=None, help="Step length in mass-weighted coordinates; overrides irc.step_length from YAML.")
@click.option("--root", type=int, default=None, help="Imaginary mode index used for the initial displacement; overrides irc.root from YAML.")
@click.option("--forward", type=bool, default=None, help="Run the forward IRC; overrides irc.forward from YAML. Specify True/False explicitly.")
@click.option("--backward", type=bool, default=None, help="Run the backward IRC; overrides irc.backward from YAML. Specify True/False explicitly.")
@click.option("--out-dir", type=str, default="./result_irc/", show_default=True, help="Output directory; overrides irc.out_dir from YAML.")
@click.option(
    "--hessian-calc-mode",
    type=click.Choice(["Analytical", "FiniteDifference"], case_sensitive=False),
    default=None,
    help="How UMA builds the Hessian (Analytical or FiniteDifference); overrides calc.hessian_calc_mode from YAML.",
)
@click.option(
    "--args-yaml",
    type=click.Path(path_type=Path, exists=True, dir_okay=False),
    default=None,
    help="YAML file providing extra parameters (sections: geom, calc, irc).",
)
@click.option(
    "--ref-pdb",
    type=click.Path(path_type=Path, exists=True, dir_okay=False),
    default=None,
    help="Reference PDB topology to use when --input is XYZ (keeps XYZ coordinates).",
)
def cli(
    input_path: Path,
    real_parm7: Optional[Path],
    model_pdb: Optional[Path],
    model_indices_str: Optional[str],
    model_indices_one_based: bool,
    detect_layer: bool,
    charge: Optional[int],
    spin: Optional[int],
    max_cycles: Optional[int],
    step_size: Optional[float],
    root: Optional[int],
    forward: Optional[bool],
    backward: Optional[bool],
    out_dir: str,
    hessian_calc_mode: Optional[str],
    args_yaml: Optional[Path],
    ref_pdb: Optional[Path],
) -> None:
    prepared_input = prepare_input_structure(input_path)
    try:
        apply_ref_pdb_override(prepared_input, ref_pdb)
    except click.BadParameter as e:
        click.echo(f"ERROR: {e}", err=True)
        prepared_input.cleanup()
        sys.exit(1)
    geom_input_path = prepared_input.geom_path
    source_path = prepared_input.source_path
    charge, spin = resolve_charge_spin_or_raise(prepared_input, charge, spin)

    model_indices: Optional[List[int]] = None
    if model_indices_str:
        try:
            model_indices = parse_indices_string(model_indices_str, one_based=model_indices_one_based)
        except click.BadParameter as e:
            click.echo(f"ERROR: {e}", err=True)
            prepared_input.cleanup()
            sys.exit(1)
    try:
        time_start = time.perf_counter()

        # --------------------------
        # 1) Assemble configuration: defaults -> YAML overrides -> CLI overrides
        # --------------------------
        yaml_cfg = load_yaml_dict(args_yaml)

        geom_cfg: Dict[str, Any] = dict(GEOM_KW_DEFAULT)
        calc_cfg: Dict[str, Any] = dict(CALC_KW_DEFAULT)
        irc_cfg:  Dict[str, Any] = dict(IRC_KW_DEFAULT)

        apply_yaml_overrides(
            yaml_cfg,
            [
                (geom_cfg, (("geom",),)),
                (calc_cfg, (("calc",),)),
                (irc_cfg, (("irc",),)),
            ],
        )

        # CLI overrides
        calc_cfg["model_charge"] = int(charge)
        calc_cfg["model_mult"] = int(spin)
        if real_parm7 is not None:
            calc_cfg["real_parm7"] = str(real_parm7)
        # Default to partial Hessian to reduce VRAM unless explicitly overridden in YAML.
        yaml_calc = {}
        if isinstance(yaml_cfg, dict):
            yaml_calc = yaml_cfg.get("calc", {}) or yaml_cfg.get("mlmm", {}) or {}
        if "return_partial_hessian" not in (yaml_calc or {}):
            calc_cfg["return_partial_hessian"] = True

        if hessian_calc_mode is not None:
            # pass through exactly as chosen; uma_pysis normalizes internally
            calc_cfg["hessian_calc_mode"] = str(hessian_calc_mode)

        if max_cycles is not None:
            irc_cfg["max_cycles"] = int(max_cycles)
        if step_size is not None:
            irc_cfg["step_length"] = float(step_size)
        if root is not None:
            irc_cfg["root"] = int(root)
        if forward is not None:
            irc_cfg["forward"] = bool(forward)
        if backward is not None:
            irc_cfg["backward"] = bool(backward)
        if out_dir:
            irc_cfg["out_dir"] = str(out_dir)

        # Normalize any existing freeze list from YAML before wiring it to UMA
        merge_freeze_atom_indices(geom_cfg)

        # Ensure the calculator receives the freeze list used by geometry
        #      (so FD Hessian can skip frozen DOF, etc.)
        calc_cfg["freeze_atoms"] = list(geom_cfg.get("freeze_atoms", []))

        calc_cfg["input_pdb"] = str(source_path)
        if not calc_cfg.get("real_parm7"):
            raise click.BadParameter("Missing --real-parm7 (or calc.real_parm7 in YAML).")

        out_dir_path = Path(irc_cfg["out_dir"]).resolve()
        out_dir_path.mkdir(parents=True, exist_ok=True)

        layer_source_pdb = source_path
        if detect_layer and layer_source_pdb.suffix.lower() != ".pdb":
            raise click.BadParameter("--detect-layer requires a PDB input (or --ref-pdb).")

        model_pdb_cfg = model_pdb if model_pdb is not None else calc_cfg.get("model_pdb")
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
                if model_pdb_cfg is None and not model_indices:
                    raise click.BadParameter(str(e))
                click.echo(f"[layer] WARNING: {e} Falling back to explicit ML region.", err=True)
                detect_layer = False

        if not detect_layer:
            if model_pdb_cfg is None and not model_indices:
                raise click.BadParameter("Provide --model-pdb or --model-indices when --no-detect-layer.")
            if model_pdb_cfg is not None:
                model_pdb_path = Path(model_pdb_cfg)
            else:
                if layer_source_pdb.suffix.lower() != ".pdb":
                    raise click.BadParameter("--model-indices requires a PDB input (or --ref-pdb).")
                try:
                    model_pdb_path = build_model_pdb_from_indices(layer_source_pdb, out_dir_path, model_indices or [])
                except Exception as e:
                    raise click.BadParameter(str(e))
            calc_cfg["use_bfactor_layers"] = False

        if model_pdb_path is None:
            raise click.BadParameter("Failed to resolve model PDB for the ML region.")

        calc_cfg["model_pdb"] = str(model_pdb_path)
        _ = apply_layer_freeze_constraints(
            geom_cfg,
            calc_cfg,
            layer_info,
            echo_fn=click.echo,
        )

        # Pretty-print configuration (expand freeze_atoms for readability)
        click.echo(pretty_block("geom", format_freeze_atoms_for_echo(geom_cfg, key="freeze_atoms")))
        echo_calc = strip_inherited_keys(calc_cfg, CALC_KW_DEFAULT, mode="same")
        echo_calc = format_freeze_atoms_for_echo(echo_calc, key="freeze_atoms")
        click.echo(pretty_block("calc", echo_calc))
        echo_irc = strip_inherited_keys({**irc_cfg, "out_dir": str(out_dir_path)}, IRC_KW_DEFAULT, mode="same")
        click.echo(pretty_block("irc", echo_irc))

        # --------------------------
        # 2) Load geometry and configure UMA calculator
        # --------------------------
        coord_type = geom_cfg.get("coord_type", "cart")
        coord_kwargs = dict(geom_cfg)
        coord_kwargs.pop("coord_type", None)

        geometry = geom_loader(geom_input_path, coord_type=coord_type, **coord_kwargs)

        # Create mlmm calculator
        calc = mlmm(**calc_cfg)
        geometry.set_calculator(calc)
        # If using partial Hessian, freeze non-Hessian atoms so IRC evolves in the same subspace.
        if calc_cfg.get("return_partial_hessian"):
            calc_core = calc.core if hasattr(calc, "core") else calc
            hess_freeze = list(getattr(calc_core, "hess_freeze_atoms", []) or [])
            if hess_freeze:
                existing_freeze = list(getattr(geometry, "freeze_atoms", []))
                freeze_union = sorted(set(existing_freeze) | set(hess_freeze))
                geometry.freeze_atoms = np.array(freeze_union, dtype=int)
                for attr in ("_active_atom_indices", "_active_dof_indices"):
                    if hasattr(geometry, attr):
                        delattr(geometry, attr)
                try:
                    calc.freeze_atoms = freeze_union
                except Exception:
                    pass

        # Seed the initial Hessian via the shared freq backend so IRC reuses
        # the same Hessian path as frequency analysis.
        hess_device = _torch_device(calc_cfg.get("ml_device", "auto"))
        h_init, _ = _calc_full_hessian_torch(
            geometry,
            calc_cfg,
            hess_device,
            refresh_geom_meta=True,
        )
        geometry.cart_hessian = h_init

        # --------------------------
        # 3) Construct and run EulerPC
        # --------------------------
        # EulerPC.__init__ forwards **kwargs directly to IRC.__init__
        eulerpc = EulerPC(geometry, **irc_cfg)

        click.echo("\n=== IRC (EulerPC) started ===\n")
        eulerpc.run()
        click.echo("\n=== IRC (EulerPC) finished ===\n")

        # --------------------------
        # 4) Convert trajectories to PDB when the input was PDB (or --ref-pdb provided)
        # --------------------------
        if source_path.suffix.lower() == ".pdb":
            ref_pdb_path = source_path.resolve()

            # Whole IRC trajectory
            _echo_convert_trj_to_pdb_if_exists(
                out_dir_path / f"{irc_cfg.get('prefix','')}{'finished_irc.trj'}",
                ref_pdb_path,
                out_dir_path / f"{irc_cfg.get('prefix','')}{'finished_irc.pdb'}",
            )
            # Forward/backward trajectories
            _echo_convert_trj_to_pdb_if_exists(
                out_dir_path / f"{irc_cfg.get('prefix','')}{'forward_irc.trj'}",
                ref_pdb_path,
                out_dir_path / f"{irc_cfg.get('prefix','')}{'forward_irc.pdb'}",
            )
            _echo_convert_trj_to_pdb_if_exists(
                out_dir_path / f"{irc_cfg.get('prefix','')}{'backward_irc.trj'}",
                ref_pdb_path,
                out_dir_path / f"{irc_cfg.get('prefix','')}{'backward_irc.pdb'}",
            )

        click.echo(format_elapsed("[time] Elapsed Time for IRC", time_start))

    except KeyboardInterrupt:
        click.echo("\nInterrupted by user.", err=True)
        sys.exit(130)
    except click.BadParameter as e:
        click.echo(f"ERROR: {e}", err=True)
        sys.exit(1)
    except Exception as e:
        tb = textwrap.indent("".join(__import__("traceback").format_exception(type(e), e, e.__traceback__)), "  ")
        click.echo("Unhandled exception during IRC:\n" + tb, err=True)
        sys.exit(1)
    finally:
        prepared_input.cleanup()


# Script entry point
if __name__ == "__main__":
    cli()
