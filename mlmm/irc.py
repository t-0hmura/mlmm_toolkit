# mlmm/irc.py

"""
ML/MM IRC calculation using the EulerPC predictor-corrector integrator.

Example:
    mlmm irc -i ts.pdb --parm real.parm7 --model-pdb ml_region.pdb -q 0

For detailed documentation, see: docs/irc.md
"""

from __future__ import annotations

from pathlib import Path
from typing import Any, Dict, Optional, List

import gc
import logging
import sys
import textwrap

logger = logging.getLogger(__name__)

import click
import numpy as np
import time
import torch

from pysisyphus.helpers import geom_loader
from pysisyphus.irc.EulerPC import EulerPC
from .mlmm_calc import mlmm
from .freq import _torch_device, _calc_full_hessian_torch, _align_three_layer_hessian_targets
from .defaults import (
    GEOM_KW_DEFAULT,
    MLMM_CALC_KW as _UMA_CALC_KW,
    IRC_KW,
)
from .utils import (
    apply_ref_pdb_override,
    apply_layer_freeze_constraints,
    convert_xyz_to_pdb,
    set_convert_file_enabled,
    is_convert_file_enabled,
    convert_xyz_like_outputs,
    load_yaml_dict,
    deep_update,
    apply_yaml_overrides,
    pretty_block,
    strip_inherited_keys,
    filter_calc_for_echo,
    format_freeze_atoms_for_echo,
    format_elapsed,
    merge_freeze_atom_indices,
    prepare_input_structure,
    resolve_charge_spin_or_raise,
    parse_indices_string,
    build_model_pdb_from_bfactors,
    build_model_pdb_from_indices,
    yaml_section_has_key,
)
from .cli_utils import resolve_yaml_sources, load_merged_yaml_cfg, make_is_param_explicit


# --------------------------
# Default configuration
# --------------------------

CALC_KW_DEFAULT: Dict[str, Any] = dict(_UMA_CALC_KW)

IRC_KW_DEFAULT: Dict[str, Any] = {
    **IRC_KW,
    "dump_fn": "irc_data.h5",
    "dump_every": 5,
}


def _echo_convert_trj_to_pdb_if_exists(trj_path: Path, ref_pdb: Path, out_path: Path) -> None:
    if not is_convert_file_enabled():
        return
    if trj_path.exists():
        try:
            convert_xyz_to_pdb(trj_path, ref_pdb, out_path)
            click.echo(f"[convert] Wrote '{out_path}'.")
        except Exception as e:
            logger.debug("Failed to convert %s to PDB", trj_path.name, exc_info=True)
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
    help="Input structure file (.pdb, .xyz, _trj.xyz, etc.).",
)
@click.option(
    "--parm",
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
    help="Detect ML/MM layers from input PDB B-factors (ML=0, MovableMM=10, FrozenMM=20). "
         "If disabled, you must provide --model-pdb or --model-indices.",
)
@click.option("-q", "--charge", type=int, required=False,
              help="Total charge; overrides calc.charge from YAML. Required unless --ligand-charge is provided.")
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
    help="Spin multiplicity (2S+1); overrides calc.spin from YAML.",
)
@click.option(
    "--max-cycles", type=int, default=None, help="Maximum number of IRC steps; overrides irc.max_cycles from YAML."
)
@click.option("--step-size", type=float, default=None, help="Step length in Bohr (unweighted Cartesian coordinates). Default: 0.10 Bohr. Overrides irc.step_length from YAML.")
@click.option("--root", type=int, default=None, help="Imaginary mode index used for the initial displacement; overrides irc.root from YAML.")
@click.option(
    "--forward/--no-forward",
    "forward",
    default=None,
    help="Run the forward IRC; overrides irc.forward from YAML.",
)
@click.option(
    "--backward/--no-backward",
    "backward",
    default=None,
    help="Run the backward IRC; overrides irc.backward from YAML.",
)
@click.option("-o", "--out-dir", type=str, default="./result_irc/", show_default=True, help="Output directory; overrides irc.out_dir from YAML.")
@click.option(
    "--hessian-calc-mode",
    type=click.Choice(["Analytical", "FiniteDifference"], case_sensitive=False),
    default=None,
    help="How the ML backend builds the Hessian (Analytical or FiniteDifference); overrides calc.hessian_calc_mode from YAML. Default: 'FiniteDifference'. Use 'Analytical' when VRAM is sufficient.",
)
@click.option(
    "--config",
    "config_yaml",
    type=click.Path(path_type=Path, exists=True, dir_okay=False),
    default=None,
    help="Base YAML configuration file applied before explicit CLI options.",
)
@click.option(
    "--show-config/--no-show-config",
    "show_config",
    default=False,
    show_default=True,
    help="Print resolved configuration and continue execution.",
)
@click.option(
    "--dry-run/--no-dry-run",
    "dry_run",
    default=False,
    show_default=True,
    help="Validate options and print the execution plan without running IRC.",
)
@click.option(
    "--ref-pdb",
    type=click.Path(path_type=Path, exists=True, dir_okay=False),
    default=None,
    help="Reference PDB topology to use when --input is XYZ (keeps XYZ coordinates).",
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
    "--hess-device",
    "hess_device",
    type=click.Choice(["auto", "cuda", "cpu"], case_sensitive=False),
    default="auto",
    show_default=True,
    help="Device for initial Hessian storage and IRC operations (auto/cuda/cpu). "
         "Use 'cpu' for large unfrozen systems to avoid VRAM limits.",
)
@click.option(
    "--read-hess",
    "read_hess",
    type=click.Path(exists=True, dir_okay=False),
    default=None,
    show_default=False,
    help="Read initial Hessian from a .npz file (produced by 'mlmm freq --dump-hess'). "
         "Takes priority over the hessian_cache and fresh computation.",
)
@click.pass_context
def cli(
    ctx: click.Context,
    input_path: Path,
    real_parm7: Optional[Path],
    model_pdb: Optional[Path],
    model_indices_str: Optional[str],
    model_indices_one_based: bool,
    detect_layer: bool,
    charge: Optional[int],
    ligand_charge: Optional[str],
    spin: Optional[int],
    max_cycles: Optional[int],
    step_size: Optional[float],
    root: Optional[int],
    forward: Optional[bool],
    backward: Optional[bool],
    out_dir: str,
    hessian_calc_mode: Optional[str],
    config_yaml: Optional[Path],
    show_config: bool,
    dry_run: bool,
    ref_pdb: Optional[Path],
    convert_files: bool,
    backend: Optional[str],
    embedcharge: bool,
    embedcharge_cutoff: Optional[float],
    link_atom_method: Optional[str],
    mm_backend: Optional[str],
    use_cmap: Optional[bool],
    hess_device: str,
    read_hess: Optional[str],
) -> None:
    set_convert_file_enabled(convert_files)
    _is_param_explicit = make_is_param_explicit(ctx)

    config_yaml, override_yaml, used_legacy_yaml = resolve_yaml_sources(
        config_yaml=config_yaml,
        override_yaml=None,
        args_yaml_legacy=None,
    )
    merged_yaml_cfg, _, _ = load_merged_yaml_cfg(
        config_yaml=config_yaml,
        override_yaml=None,
    )

    prepared_input = prepare_input_structure(input_path)
    try:
        apply_ref_pdb_override(prepared_input, ref_pdb)
    except click.BadParameter as e:
        click.echo(f"ERROR: {e}", err=True)
        prepared_input.cleanup()
        sys.exit(1)
    geom_input_path = prepared_input.geom_path
    source_path = prepared_input.source_path
    charge, spin = resolve_charge_spin_or_raise(
        prepared_input, charge, spin,
        ligand_charge=ligand_charge, prefix="[irc]",
    )

    model_indices: Optional[List[int]] = None
    if model_indices_str:
        try:
            model_indices = parse_indices_string(model_indices_str, one_based=model_indices_one_based)
        except click.BadParameter as e:
            click.echo(f"ERROR: {e}", err=True)
            prepared_input.cleanup()
            sys.exit(1)
    calc = eulerpc = geometry = None
    try:
        time_start = time.perf_counter()

        # --------------------------
        # 1) Assemble configuration: defaults < config < CLI(explicit) < override
        # --------------------------
        config_layer_cfg = load_yaml_dict(config_yaml)
        override_layer_cfg = load_yaml_dict(override_yaml)

        geom_cfg: Dict[str, Any] = dict(GEOM_KW_DEFAULT)
        calc_cfg: Dict[str, Any] = dict(CALC_KW_DEFAULT)
        irc_cfg: Dict[str, Any] = dict(IRC_KW_DEFAULT)
        # Keep command-level default for detect-layer unless YAML/explicit CLI overrides it.
        calc_cfg["use_bfactor_layers"] = bool(detect_layer)

        apply_yaml_overrides(
            config_layer_cfg,
            [
                (geom_cfg, (("geom",),)),
                (calc_cfg, (("calc",), ("mlmm",))),
                (irc_cfg, (("irc",),)),
            ],
        )

        # CLI explicit overrides (after config YAML, before override YAML)
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

        if _is_param_explicit("hessian_calc_mode") and hessian_calc_mode is not None:
            calc_cfg["hessian_calc_mode"] = str(hessian_calc_mode)
        if _is_param_explicit("max_cycles") and max_cycles is not None:
            irc_cfg["max_cycles"] = int(max_cycles)
        if _is_param_explicit("step_size") and step_size is not None:
            irc_cfg["step_length"] = float(step_size)
        if _is_param_explicit("root") and root is not None:
            irc_cfg["root"] = int(root)
        if _is_param_explicit("forward") and forward is not None:
            irc_cfg["forward"] = bool(forward)
        if _is_param_explicit("backward") and backward is not None:
            irc_cfg["backward"] = bool(backward)
        if _is_param_explicit("out_dir"):
            irc_cfg["out_dir"] = str(out_dir)
        if _is_param_explicit("detect_layer"):
            calc_cfg["use_bfactor_layers"] = bool(detect_layer)

        model_charge_value = calc_cfg.get("model_charge", charge)
        if model_charge_value is None:
            model_charge_value = charge
        calc_cfg["model_charge"] = int(model_charge_value)
        if _is_param_explicit("charge"):
            calc_cfg["model_charge"] = int(charge)

        model_mult_value = calc_cfg.get("model_mult", spin)
        if model_mult_value is None:
            model_mult_value = spin
        calc_cfg["model_mult"] = int(model_mult_value)
        if _is_param_explicit("spin"):
            calc_cfg["model_mult"] = int(spin)

        calc_cfg["input_pdb"] = str(source_path)
        if real_parm7 is not None:
            calc_cfg["real_parm7"] = str(real_parm7)
        if model_pdb is not None:
            calc_cfg["model_pdb"] = str(model_pdb)

        apply_yaml_overrides(
            override_layer_cfg,
            [
                (geom_cfg, (("geom",),)),
                (calc_cfg, (("calc",), ("mlmm",))),
                (irc_cfg, (("irc",),)),
            ],
        )
        calc_paths = (("calc",), ("mlmm",))
        partial_explicit = (
            yaml_section_has_key(config_layer_cfg, calc_paths, "return_partial_hessian")
            or yaml_section_has_key(override_layer_cfg, calc_paths, "return_partial_hessian")
        )
        if not partial_explicit:
            calc_cfg["return_partial_hessian"] = True

        # Normalize any existing freeze list from YAML before wiring it to UMA
        merge_freeze_atom_indices(geom_cfg)
        calc_cfg["freeze_atoms"] = list(geom_cfg.get("freeze_atoms", []))
        if not calc_cfg.get("real_parm7"):
            raise click.BadParameter("Missing --parm (or calc.real_parm7 in YAML).")

        out_dir_path = Path(irc_cfg["out_dir"]).resolve()
        layer_source_pdb = source_path
        detect_layer_enabled = bool(calc_cfg.get("use_bfactor_layers", True))
        model_pdb_cfg = calc_cfg.get("model_pdb")

        if show_config:
            click.echo(
                pretty_block(
                    "yaml_layers",
                    {
                        "config": None if config_yaml is None else str(config_yaml),
                        "override_yaml": None if override_yaml is None else str(override_yaml),
                        "merged_keys": sorted(merged_yaml_cfg.keys()),
                    },
                )
            )

        if dry_run:
            model_region_source = "bfactor"
            if not detect_layer_enabled:
                if model_pdb_cfg is not None:
                    model_region_source = "model_pdb"
                elif model_indices:
                    model_region_source = "model_indices"
                else:
                    raise click.BadParameter("Provide --model-pdb or --model-indices when --no-detect-layer.")
            if detect_layer_enabled and layer_source_pdb.suffix.lower() != ".pdb":
                raise click.BadParameter("--detect-layer requires a PDB input (or --ref-pdb).")
            if (
                not detect_layer_enabled
                and model_pdb_cfg is None
                and model_indices
                and layer_source_pdb.suffix.lower() != ".pdb"
            ):
                raise click.BadParameter("--model-indices requires a PDB input (or --ref-pdb).")
            click.echo(
                pretty_block(
                    "dry_run_plan",
                    {
                        "input_geometry": str(geom_input_path),
                        "output_dir": str(out_dir_path),
                        "detect_layer": bool(detect_layer_enabled),
                        "model_region_source": model_region_source,
                        "model_indices_count": 0 if not model_indices else len(model_indices),
                        "will_run_irc": True,
                        "will_write_trajectories": True,
                        "backend": calc_cfg.get("backend", "uma"),
                        "embedcharge": bool(calc_cfg.get("embedcharge", False)),
                    },
                )
            )
            click.echo("[dry-run] Validation complete. IRC execution was skipped.")
            return

        out_dir_path.mkdir(parents=True, exist_ok=True)

        if detect_layer_enabled and layer_source_pdb.suffix.lower() != ".pdb":
            raise click.BadParameter("--detect-layer requires a PDB input (or --ref-pdb).")

        model_pdb_path: Optional[Path] = None
        layer_info: Optional[Dict[str, List[int]]] = None

        if detect_layer_enabled:
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
                detect_layer_enabled = False

        if not detect_layer_enabled:
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
        _align_three_layer_hessian_targets(calc_cfg, echo_fn=click.echo)

        # Pretty-print configuration (expand freeze_atoms for readability)
        click.echo(pretty_block("geom", format_freeze_atoms_for_echo(geom_cfg, key="freeze_atoms")))
        echo_calc = format_freeze_atoms_for_echo(filter_calc_for_echo(calc_cfg), key="freeze_atoms")
        click.echo(pretty_block("calc", echo_calc, defaults=_UMA_CALC_KW))
        echo_irc = strip_inherited_keys({**irc_cfg, "out_dir": str(out_dir_path)}, IRC_KW_DEFAULT, mode="same")
        click.echo(pretty_block("irc", echo_irc))

        # --------------------------
        # 2) Load geometry and configure UMA calculator
        # --------------------------
        geom_cfg["coord_type"] = "cart"  # IRC requires Cartesian coordinates
        coord_type = "cart"
        coord_kwargs = dict(geom_cfg)
        coord_kwargs.pop("coord_type", None)

        geometry = geom_loader(geom_input_path, coord_type=coord_type, **coord_kwargs)

        # Create mlmm calculator
        calc = mlmm(**calc_cfg)
        geometry.set_calculator(calc)

        # Seed the initial Hessian.
        # Priority: --read-hess file > hessian_cache > fresh computation.
        from .hessian_cache import load as _hess_load, store as _hess_store
        if hess_device.lower() == "auto":
            _hess_dev = _torch_device(calc_cfg.get("ml_device", "auto"))
        else:
            _hess_dev = _torch_device(hess_device.lower())
        if _hess_dev.type == "cpu":
            click.echo("[device] Hessian operations will run on CPU.")

        if read_hess:
            click.echo(f"[irc] Loading initial Hessian from {read_hess}")
            _data = np.load(read_hess)
            h_init = torch.as_tensor(_data["hessian"], dtype=torch.float64, device=_hess_dev)
            del _data
        elif (cached := _hess_load("ts")) is not None:
            click.echo("[irc] Reusing cached TS Hessian from tsopt.")
            active_dofs = cached.get("active_dofs")
            h_raw = cached["hessian"]
            if isinstance(h_raw, torch.Tensor):
                h_init = h_raw.to(device=_hess_dev)
            else:
                h_init = torch.as_tensor(h_raw, dtype=torch.float64, device=_hess_dev)
            if active_dofs is not None:
                geometry.within_partial_hessian = {
                    "active_n_dof": len(active_dofs),
                    "full_n_dof": geometry.cart_coords.size,
                    "active_dofs": active_dofs,
                    "active_atoms": sorted(set(d // 3 for d in active_dofs)),
                }
        else:
            click.echo("[irc] Seeding initial Hessian via shared freq backend.")
            h_init, _ = _calc_full_hessian_torch(
                geometry,
                calc_cfg,
                _hess_dev,
                refresh_geom_meta=True,
            )

        geometry.cart_hessian = h_init
        click.echo(f"[irc] Initial Hessian seeded (shape={h_init.shape[0]}x{h_init.shape[1]}).")
        del h_init

        # --------------------------
        # 3) Construct and run EulerPC
        # --------------------------
        # EulerPC.__init__ forwards **kwargs directly to IRC.__init__
        eulerpc = EulerPC(geometry, **irc_cfg)

        click.echo("\n=== IRC (EulerPC) started ===\n")
        eulerpc.run()
        click.echo("\n=== IRC (EulerPC) finished ===\n")

        # Cache IRC endpoint Hessians (Bofill-updated mw → Cartesian)
        def _unmw_and_store(mw_H, key):
            """Un-mass-weight active-DOF Hessian on device, store partial on CPU."""
            import numpy as np
            act = eulerpc._act_dofs
            m_sqrt = geometry.masses_rep ** 0.5
            ms_act = m_sqrt[act]
            if isinstance(mw_H, torch.Tensor):
                ms_t = torch.as_tensor(ms_act, dtype=mw_H.dtype, device=mw_H.device)
                H_cart_act = ms_t.unsqueeze(1) * mw_H * ms_t.unsqueeze(0)
                H_cart_act_np = H_cart_act.detach().cpu().numpy()
            else:
                H_cart_act_np = np.diag(ms_act) @ mw_H @ np.diag(ms_act)
            _hess_store(key, H_cart_act_np, active_dofs=list(act))

        if getattr(eulerpc, "forward_mw_hessian", None) is not None:
            _unmw_and_store(eulerpc.forward_mw_hessian, "irc_left")
            click.echo("[irc] Cached forward endpoint Hessian as 'irc_left'.")
        if getattr(eulerpc, "mw_hessian", None) is not None:
            _unmw_and_store(eulerpc.mw_hessian, "irc_right")
            click.echo("[irc] Cached backward endpoint Hessian as 'irc_right'.")

        # --------------------------
        # 4) Convert trajectories to PDB when the input was PDB (or --ref-pdb provided)
        # --------------------------
        if source_path.suffix.lower() == ".pdb":
            ref_pdb_path = source_path.resolve()

            # Whole IRC trajectory
            _echo_convert_trj_to_pdb_if_exists(
                out_dir_path / f"{irc_cfg.get('prefix','')}{'finished_irc_trj.xyz'}",
                ref_pdb_path,
                out_dir_path / f"{irc_cfg.get('prefix','')}{'finished_irc.pdb'}",
            )
            # Forward/backward trajectories
            _echo_convert_trj_to_pdb_if_exists(
                out_dir_path / f"{irc_cfg.get('prefix','')}{'forward_irc_trj.xyz'}",
                ref_pdb_path,
                out_dir_path / f"{irc_cfg.get('prefix','')}{'forward_irc.pdb'}",
            )
            _echo_convert_trj_to_pdb_if_exists(
                out_dir_path / f"{irc_cfg.get('prefix','')}{'backward_irc_trj.xyz'}",
                ref_pdb_path,
                out_dir_path / f"{irc_cfg.get('prefix','')}{'backward_irc.pdb'}",
            )
            # Single-frame endpoint PDBs (forward_last, backward_last)
            prefix = irc_cfg.get("prefix", "")
            for tag in ("forward_last", "backward_last"):
                endpoint_xyz = out_dir_path / f"{prefix}{tag}.xyz"
                endpoint_pdb = out_dir_path / f"{prefix}{tag}.pdb"
                if endpoint_xyz.exists() and not endpoint_pdb.exists():
                    try:
                        convert_xyz_to_pdb(endpoint_xyz, ref_pdb_path, endpoint_pdb)
                        click.echo(f"[convert] Wrote '{endpoint_pdb}'.")
                    except Exception as e:
                        logger.debug("Failed to convert %s to PDB", endpoint_xyz.name, exc_info=True)
                        click.echo(f"[convert] WARNING: Failed to convert '{tag}.xyz' to PDB: {e}", err=True)

        # summary.md and key_* outputs are disabled.
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
        # Release GPU memory (model + Hessian) so subsequent stages don't OOM
        calc = eulerpc = geometry = None
        gc.collect()  # break cyclic refs inside torch.nn.Module
        if torch.cuda.is_available():
            torch.cuda.empty_cache()


# Script entry point
if __name__ == "__main__":
    cli()
