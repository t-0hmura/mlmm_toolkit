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

logger = logging.getLogger(__name__)

import click
import numpy as np
import time
import torch

from pysisyphus.helpers import geom_loader
from pysisyphus.irc.EulerPC import EulerPC
from mlmm.backends.mlmm_calc import mlmm
from mlmm.workflows.freq import _torch_device, _calc_full_hessian_torch, _align_three_layer_hessian_targets
from mlmm.core.defaults import (
    GEOM_KW_DEFAULT,
    MLMM_CALC_KW as _UMA_CALC_KW,
    IRC_KW,
)
from mlmm.core.utils import (
    apply_ref_pdb_override,
    apply_layer_freeze_constraints,
    convert_xyz_to_pdb,
    set_convert_file_enabled,
    is_convert_file_enabled,
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
    parse_indices_string,
    build_model_pdb_from_bfactors,
    build_model_pdb_from_indices,
    yaml_section_has_key,
    echo_resolved_device,
)
from mlmm.cli.common_options import (
    add_ml_layer_detection_options,
    add_precision_option, add_backend_model_option,
    add_deterministic_option, add_allow_charge_mult_mismatch_option,
    add_irc_pos_def_option,
)
from mlmm.cli.decorators import resolve_yaml_sources, load_merged_yaml_cfg, make_is_param_explicit, _write_error_json, render_cli_exception



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
@click.option("-o", "--out-dir", type=str, default=IRC_KW["out_dir"], show_default=True, help="Output directory; overrides irc.out_dir from YAML.")
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
@click.option(
    "--freeze-atoms",
    "freeze_atoms_text",
    type=str,
    default=None,
    show_default=False,
    help="Comma-separated 1-based atom indices to freeze (e.g., '1,3,5').",
)
@click.option(
    "--out-json/--no-out-json",
    "out_json",
    default=False,
    show_default=True,
    help="Write machine-readable result.json to out_dir.",
)
@add_ml_layer_detection_options()
@add_precision_option()
@add_backend_model_option()
@add_deterministic_option()
@add_allow_charge_mult_mismatch_option()
@add_irc_pos_def_option()
@click.pass_context
def cli(
    ctx: click.Context,
    input_path: Path,
    real_parm7: Optional[Path],
    model_pdb: Optional[Path],
    model_indices_str: Optional[str],
    model_indices_one_based: bool,
    detect_layer: bool,
    freeze_atoms_text: Optional[str],
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
    out_json: bool,
    precision: Optional[str],
    backend_model: Optional[str],
    irc_pos_def: Optional[bool],
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

        config_layer_cfg = load_yaml_dict(config_yaml)
        override_layer_cfg = load_yaml_dict(override_yaml)

        geom_cfg: Dict[str, Any] = dict(GEOM_KW_DEFAULT)
        calc_cfg: Dict[str, Any] = dict(CALC_KW_DEFAULT)
        irc_cfg: Dict[str, Any] = dict(IRC_KW_DEFAULT)
        # DO NOT INLINE: (irc): detect-layer has different defaults across subcommands (freq/irc default True). Guard prevents YAML accidentally suppressing CLI-level default.
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
        if precision is not None:
            from mlmm.backends import apply_precision_to_calc_cfg, apply_backend_model_to_calc_cfg
            apply_precision_to_calc_cfg(calc_cfg, precision)
            apply_backend_model_to_calc_cfg(calc_cfg, backend_model)
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
        # CLI knobs → irc_cfg. require_pos_def_hessian = PSD-Hessian convergence guard.
        if _is_param_explicit("irc_pos_def") and irc_pos_def is not None:
            irc_cfg["require_pos_def_hessian"] = bool(irc_pos_def)
        if _is_param_explicit("detect_layer"):
            calc_cfg["use_bfactor_layers"] = bool(detect_layer)

        # CLI-resolved charge/spin (from -q / -l derivation, or -m / spin_default)
        # always wins over the CALC_KW default carried in calc_cfg.
        calc_cfg["model_charge"] = int(charge)
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
        if freeze_atoms_text:
            from mlmm.workflows.opt import _parse_freeze_atoms
            freeze_cli = _parse_freeze_atoms(freeze_atoms_text)
            merge_freeze_atom_indices(geom_cfg, freeze_cli)
        else:
            merge_freeze_atom_indices(geom_cfg)
        calc_cfg["freeze_atoms"] = list(geom_cfg.get("freeze_atoms", []))
        from mlmm.workflows.opt import _convert_yaml_layer_atoms_1to0
        _convert_yaml_layer_atoms_1to0(calc_cfg)
        if not calc_cfg.get("real_parm7"):
            raise click.BadParameter(
                "Missing --parm (or calc.real_parm7 in YAML).; "
                "recover: pass --parm /path/to/real.parm7 OR add 'calc:\\n  real_parm7: ...' to --config YAML."
            )

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

        # Default-verbosity entry summary (skipped in child mode).
        from mlmm.core.utils import echo_run_summary
        _backend = calc_cfg.get("backend") or "uma"
        _model = calc_cfg.get("uma_model") or calc_cfg.get("model")
        _precision = calc_cfg.get("uma_precision") or calc_cfg.get("precision", "fp32")
        echo_run_summary({
            "input": str(input_path),
            "backend": f"{_backend} ({_model}, {_precision})" if _model else _backend,
            "out": str(out_dir_path),
        })

        # Pretty-print configuration (expand freeze_atoms for readability)
        click.echo(pretty_block("geom", format_freeze_atoms_for_echo(geom_cfg, key="freeze_atoms")))
        echo_calc = format_freeze_atoms_for_echo(filter_calc_for_echo(calc_cfg), key="freeze_atoms")
        click.echo(pretty_block("calc", echo_calc))
        echo_irc = strip_inherited_keys({**irc_cfg, "out_dir": str(out_dir_path)}, IRC_KW_DEFAULT, mode="same")
        click.echo(pretty_block("irc", echo_irc))

        geom_cfg["coord_type"] = "cart"  # IRC requires Cartesian coordinates
        coord_type = "cart"
        coord_kwargs = dict(geom_cfg)
        coord_kwargs.pop("coord_type", None)

        geometry = geom_loader(geom_input_path, coord_type=coord_type, **coord_kwargs)

        # Create mlmm calculator
        calc = mlmm(**calc_cfg)
        geometry.set_calculator(calc)

        echo_resolved_device()

        # Seed the initial Hessian.
        # Priority: --read-hess file > hessian_cache > fresh computation.
        from mlmm.io.hessian_cache import load as _hess_load, store as _hess_store
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
            # Restore partial-Hessian metadata if freq --dump-hess saved it,
            # so a partial Hessian (active_n_dof != 3N) is consumed correctly
            # instead of tripping the Geometry cart_hessian shape assertion.
            if "wph_active_dofs" in getattr(_data, "files", []):
                _ad = [int(d) for d in _data["wph_active_dofs"].tolist()]
                if _ad:
                    geometry.within_partial_hessian = {
                        "active_n_dof": int(_data["wph_active_n_dof"]) if "wph_active_n_dof" in _data.files else len(_ad),
                        "full_n_dof": int(_data["wph_full_n_dof"]) if "wph_full_n_dof" in _data.files else int(geometry.cart_coords.size),
                        "active_dofs": _ad,
                        "active_atoms": sorted(set(d // 3 for d in _ad)),
                    }
                    click.echo(
                        f"[irc] Restored partial-Hessian metadata from npz "
                        f"(active_n_dof={len(_ad)})."
                    )
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

        # --- Fix A: reduce the seeded Hessian to the ML macro sub-block
        # (ML atoms + link-MM parents), matching the microiteration macro step
        # used by tsopt/opt. In 3-layer mode the TS path caches a
        # Hessian-target Hessian spanning ML+MovableMM (e.g. 14472 DOF);
        # running IRC's initial_displacement eigh on that dense matrix OOMs on
        # a 16 GB GPU. tsopt avoids this because its macro RFO step operates on
        # the ML sub-block (e.g. 2652 DOF) via opt.py's sub-block extraction.
        # Mirror that here so IRC handles the Hessian the same way as TS.
        # Robustness (opt.py-consistent): the cached-"ts" basis may omit a few
        # link-parent DOFs, so we reduce to the macro-DOF intersection of the
        # seeded-Hessian basis (still ~ML-sized → no OOM); and if a large
        # Hessian cannot be reduced from its current basis at all, recompute a
        # FRESH Hessian and reduce that — exactly how opt.py falls back —
        # instead of silently keeping the full (OOM-prone) matrix. Guarded:
        # non-microiter / 2-layer / already-small Hessians and any resolution
        # failure leave behavior unchanged (no regression). --read-hess is
        # never silently recomputed (the user supplied an explicit Hessian).
        try:
            from mlmm.workflows.freq import _collect_layer_atom_sets as _f_collect_layer_sets

            _layer_sets = _f_collect_layer_sets(calc_cfg)
            _ml_idx = set(int(i) for i in _layer_sets.get("ml", set()))
            _reduced = False
            if _ml_idx:
                _core = getattr(calc, "core", calc)
                _link_parents = set(
                    int(mm_1) - 1 for (_ml_1, mm_1) in getattr(_core, "mlmm_links", [])
                )
                _macro_atoms = sorted(_ml_idx | _link_parents)
                _macro_dofs = []
                for _a in _macro_atoms:
                    _macro_dofs.extend((3 * _a, 3 * _a + 1, 3 * _a + 2))
                _macro_set = set(_macro_dofs)

                def _try_reduce(_h):
                    """Reduce _h to the ML-macro sub-block of its own basis.

                    Returns (hessian, reduced?). Sets geometry.within_partial_hessian
                    to exactly the kept (global) DOFs so IRC's _act_dofs matches.
                    """
                    _w = getattr(geometry, "within_partial_hessian", None)
                    if _w is not None and _w.get("active_dofs") is not None:
                        _basis = [int(d) for d in _w["active_dofs"]]
                    elif _h.shape[0] == geometry.cart_coords.size:
                        _basis = list(range(int(geometry.cart_coords.size)))
                    else:
                        return _h, False
                    _sub = [i for i, d in enumerate(_basis) if d in _macro_set]
                    if not _sub or len(_sub) >= _h.shape[0]:
                        return _h, False
                    _kept = [_basis[i] for i in _sub]
                    _idx = torch.as_tensor(_sub, dtype=torch.long, device=_h.device)
                    _h = _h.index_select(0, _idx).index_select(1, _idx)
                    geometry.within_partial_hessian = {
                        "active_n_dof": len(_kept),
                        "full_n_dof": int(geometry.cart_coords.size),
                        "active_dofs": _kept,
                        "active_atoms": sorted(set(d // 3 for d in _kept)),
                    }
                    return _h, True

                h_init, _reduced = _try_reduce(h_init)
                if (not _reduced) and (not read_hess) and h_init.shape[0] > len(_macro_dofs):
                    click.echo(
                        "[irc] Cached/seeded Hessian could not be reduced to the "
                        "ML macro sub-block; recomputing a fresh Hessian "
                        "(opt.py-consistent fallback)."
                    )
                    geometry.within_partial_hessian = None
                    h_init, _ = _calc_full_hessian_torch(
                        geometry, calc_cfg, _hess_dev, refresh_geom_meta=True,
                    )
                    h_init, _reduced = _try_reduce(h_init)
                if _reduced:
                    click.echo(
                        f"[irc] Reduced seeded Hessian to ML macro sub-block "
                        f"(-> {h_init.shape[0]}x{h_init.shape[0]}; "
                        f"ML+link~{len(_macro_atoms)} atoms) to match the tsopt "
                        f"macro step and avoid IRC OOM."
                    )
        except Exception as _e:  # pragma: no cover - defensive guard
            click.echo(
                f"[irc] WARNING: ML-macro Hessian reduction skipped ({_e}); "
                f"using the originally seeded Hessian."
            )

        geometry.cart_hessian = h_init
        click.echo(f"[irc] Initial Hessian seeded (shape={h_init.shape[0]}x{h_init.shape[1]}).")
        del h_init

        eulerpc = EulerPC(geometry, **irc_cfg)

        eulerpc.run()

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
                # disk cache contract; free GPU copy after npy dump.
                del H_cart_act
                if torch.cuda.is_available():
                    torch.cuda.empty_cache()
            else:
                H_cart_act_np = np.diag(ms_act) @ mw_H @ np.diag(ms_act)
            _hess_store(key, H_cart_act_np, active_dofs=list(act))

        if getattr(eulerpc, "forward_mw_hessian", None) is not None:
            _unmw_and_store(eulerpc.forward_mw_hessian, "irc_left")
            click.echo("[irc] Cached forward endpoint Hessian as 'irc_left'.")
        if getattr(eulerpc, "mw_hessian", None) is not None:
            _unmw_and_store(eulerpc.mw_hessian, "irc_right")
            click.echo("[irc] Cached backward endpoint Hessian as 'irc_right'.")

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
        click.echo(format_elapsed("[time] Elapsed Time for IRC", time_start), narrative=True)

        if out_json:
            from mlmm.core.utils import write_result_json
            _all_e = eulerpc.all_energies
            _n_fwd = len(getattr(eulerpc, "forward_energies", [])) if hasattr(eulerpc, "forward_energies") else 0
            _n_bwd = len(getattr(eulerpc, "backward_energies", [])) if hasattr(eulerpc, "backward_energies") else 0
            _ts_e = float(eulerpc.ts_energy) if hasattr(eulerpc, "ts_energy") else None
            _e_reactant = float(_all_e[0]) if len(_all_e) > 0 else None
            _e_product = float(_all_e[-1]) if len(_all_e) > 0 else None
            _irc_files = {}
            prefix = irc_cfg.get("prefix", "")
            for _fn in ("finished_irc_trj.xyz", "forward_irc_trj.xyz", "backward_irc_trj.xyz"):
                _fp = out_dir_path / f"{prefix}{_fn}"
                if _fp.exists():
                    _irc_files[_fn.replace("_trj.xyz", "")] = _fp.name
            for _fn in ("finished_irc.pdb", "forward_irc.pdb", "backward_irc.pdb"):
                _fp = out_dir_path / f"{prefix}{_fn}"
                if _fp.exists():
                    _irc_files[_fn.replace(".pdb", "_pdb")] = _fp.name
            result_data = {
                "status": "completed",
                "n_frames_forward": _n_fwd,
                "n_frames_backward": _n_bwd,
                "n_frames_total": len(_all_e),
                "energy_reactant_hartree": _e_reactant,
                "energy_ts_hartree": _ts_e,
                "energy_product_hartree": _e_product,
                "forward_converged": getattr(eulerpc, 'forward_is_converged', None),
                "backward_converged": getattr(eulerpc, 'backward_is_converged', None),
                "backend": calc_cfg.get("backend", "uma"),
                "charge": calc_cfg.get("model_charge"),
                "spin": calc_cfg.get("model_mult"),
                "n_freeze_atoms": len(geom_cfg.get("freeze_atoms", [])),
                "step_length": irc_cfg.get("step_length"),
                "max_cycles": irc_cfg.get("max_cycles"),
                "input_file": str(source_path),
                "files": _irc_files,
            }

            # Bond changes between IRC endpoints
            try:
                from mlmm.domain.bond_changes import compare_structures
                _irc_first_xyz = out_dir_path / f"{prefix}forward_last.xyz"
                _irc_last_xyz = out_dir_path / f"{prefix}backward_last.xyz"
                if _irc_first_xyz.exists() and _irc_last_xyz.exists():
                    _g1 = geom_loader(str(_irc_first_xyz))
                    _g2 = geom_loader(str(_irc_last_xyz))
                    _bc = compare_structures(_g1, _g2, device="cpu")
                    _elems = [a.capitalize() for a in _g1.atoms]
                    result_data["bond_changes"] = {
                        "formed": [f"{_elems[i]}{i+1}-{_elems[j]}{j+1}" for i, j in sorted(_bc.formed_covalent)],
                        "broken": [f"{_elems[i]}{i+1}-{_elems[j]}{j+1}" for i, j in sorted(_bc.broken_covalent)],
                    }
            except Exception:
                logger.debug("irc: bond-changes enrichment skipped", exc_info=True)

            write_result_json(
                out_dir_path, result_data,
                command="irc",
                elapsed_seconds=time.perf_counter() - time_start,
            )

    except KeyboardInterrupt:
        click.echo("\nInterrupted by user.", err=True)
        sys.exit(130)
    except click.BadParameter as e:
        _write_error_json(Path(out_dir).resolve(), "irc", e, "BadParameter", time_start)
        click.echo(f"ERROR: {e}", err=True)
        sys.exit(1)
    except Exception as e:
        render_cli_exception(e, label="IRC", out_dir=out_dir, command="irc", time_start=time_start)
    finally:
        prepared_input.cleanup()
        # DO NOT INLINE: bare `= None` leaves the name bound; PyTorch hooks captured inside eulerpc closures could resurrect the model. Two-step pattern (=None then del) is mandatory; pure `= None` is insufficient for ONIOM stack.
        # Release GPU memory (model + Hessian) so subsequent stages don't OOM.
        # `= None` decref's the heavy refs; `del` then removes names from
        # the local frame so torch.nn.Module hooks / closures cannot retain.
        calc = eulerpc = geometry = None
        del calc, eulerpc, geometry
        gc.collect()  # break cyclic refs inside torch.nn.Module
        if torch.cuda.is_available():
            torch.cuda.empty_cache()


# Script entry point
if __name__ == "__main__":
    cli()
