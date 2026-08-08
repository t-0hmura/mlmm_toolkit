"""
ML/MM IRC calculation using the EulerPC predictor-corrector integrator.

Example:
    mlmm irc -i ts.pdb --parm real.parm7 --model-pdb ml_region.pdb -q 0

For detailed documentation, see: docs/irc.md
"""

from __future__ import annotations

from pathlib import Path
from typing import Any, Dict, Optional, List, Tuple

import gc
import logging
import sys

logger = logging.getLogger(__name__)

import click
from mlmm.core.output import emit
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
from mlmm.workflows.charge_prep import resolve_charge_spin_or_raise
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
    parse_indices_string,
    resolve_ml_layer_assignment,
    yaml_section_has_key,
    echo_resolved_device,
)
from mlmm.cli.common_options import (
    add_ml_layer_detection_options,
    add_precision_option, add_backend_model_option, add_calc_file_option,
    add_workers_options,
    add_deterministic_option, add_allow_charge_mult_mismatch_option,
    add_irc_pos_def_option,
)
from mlmm.cli.decorators import resolve_yaml_sources, load_merged_yaml_cfg, make_is_param_explicit, _write_error_json, render_cli_exception



CALC_KW_DEFAULT: Dict[str, Any] = dict(_UMA_CALC_KW)

IRC_KW_DEFAULT: Dict[str, Any] = dict(IRC_KW)


def _validate_irc_directions(irc_cfg: Dict[str, Any]) -> None:
    """Require at least one executable IRC direction."""

    if bool(irc_cfg.get("downhill", False)):
        raise click.BadParameter(
            "irc.downhill is not supported; select forward and/or backward."
        )
    if not bool(irc_cfg.get("forward", False)) and not bool(
        irc_cfg.get("backward", False)
    ):
        raise click.BadParameter(
            "Enable at least one IRC direction: forward or backward."
        )


def _directional_endpoint_energy_fields(
    all_energies: Any,
    ts_energy: Any,
) -> Dict[str, Any]:
    """Report standalone IRC endpoints without inventing reactant/product identity."""
    first = float(all_energies[0]) if len(all_energies) > 0 else None
    last = float(all_energies[-1]) if len(all_energies) > 0 else None
    ts = float(ts_energy) if ts_energy is not None else None
    return {
        "energy_first_hartree": first,
        "energy_ts_hartree": ts,
        "energy_last_hartree": last,
        "endpoint_energy_orientation": "finished_first_to_finished_last",
        # Retained for schema compatibility; their orientation is declared above.
        "energy_reactant_hartree": first,
        "energy_product_hartree": last,
    }


def _consume_mw_hessian_to_cartesian_active(
    mw_hessian: Any,
    mass_sqrt_active: Any,
) -> np.ndarray:
    """Mass-unweight a terminal active Hessian in place and return it on CPU.

    The input is consumed: callers must first detach it from the EulerPC owner
    and must not reuse it. This avoids a second dense device allocation at the
    IRC-to-cache stage boundary.
    """
    if isinstance(mw_hessian, torch.Tensor):
        ms_t = torch.as_tensor(
            mass_sqrt_active,
            dtype=mw_hessian.dtype,
            device=mw_hessian.device,
        )
        with torch.no_grad():
            mw_hessian.mul_(ms_t.unsqueeze(1))
            mw_hessian.mul_(ms_t.unsqueeze(0))
        return mw_hessian.detach().cpu().numpy()

    result = np.asarray(mw_hessian)
    masses = np.asarray(mass_sqrt_active, dtype=result.dtype)
    result *= masses[:, None]
    result *= masses[None, :]
    return result


def _irc_output_path(eulerpc: EulerPC, filename: str) -> Path:
    """Resolve an engine-authored IRC filename, including normalized prefix."""
    return Path(eulerpc.get_path_for_fn(filename))


_IRC_GENERATION_FILENAMES = tuple(
    f"{stem}{suffix}"
    for stem in (
        "finished_irc",
        "forward_irc",
        "backward_irc",
        "finished_first",
        "finished_last",
        "forward_first",
        "forward_last",
        "backward_first",
        "backward_last",
    )
    for suffix in (
        ("_trj.xyz", ".pdb", ".cif")
        if stem.endswith("_irc")
        else (".xyz", ".pdb", ".cif")
    )
)


class _IRCOutputCollisionError(click.UsageError):
    """An IRC output/input collision."""


def _prepare_irc_output_dir(
    path: Path,
    *,
    prefix: str = "",
    protected_inputs: Tuple[Optional[Path], ...] = (),
) -> Path:
    """Invalidate command-owned IRC artifacts before a real generation."""
    resolved = Path(path).resolve()
    resolved.mkdir(parents=True, exist_ok=True)
    normalized_prefix = f"{prefix}_" if prefix else ""
    owned = [
        *(resolved / f"{normalized_prefix}{name}" for name in _IRC_GENERATION_FILENAMES),
        resolved / "result.json",
        resolved / "summary.json",
    ]
    reserved = {candidate.resolve() for candidate in owned}
    for protected in protected_inputs:
        if protected is not None and Path(protected).resolve() in reserved:
            raise _IRCOutputCollisionError(
                f"Input {protected} collides with a reserved IRC output path "
                f"under {resolved}."
            )
    for candidate in owned:
        candidate.unlink(missing_ok=True)
    return resolved


def _collect_irc_output_files(eulerpc: EulerPC) -> Dict[str, str]:
    """Collect normalized-prefix XYZ/PDB/CIF trajectory and endpoint outputs."""
    specs = (
        ("finished_irc_trj.xyz", "finished_irc"),
        ("forward_irc_trj.xyz", "forward_irc"),
        ("backward_irc_trj.xyz", "backward_irc"),
        ("finished_irc.pdb", "finished_irc_pdb"),
        ("forward_irc.pdb", "forward_irc_pdb"),
        ("backward_irc.pdb", "backward_irc_pdb"),
        ("finished_irc.cif", "finished_irc_cif"),
        ("forward_irc.cif", "forward_irc_cif"),
        ("backward_irc.cif", "backward_irc_cif"),
        ("forward_last.xyz", "forward_last"),
        ("backward_last.xyz", "backward_last"),
        ("forward_last.pdb", "forward_last_pdb"),
        ("backward_last.pdb", "backward_last_pdb"),
        ("forward_last.cif", "forward_last_cif"),
        ("backward_last.cif", "backward_last_cif"),
        ("forward_first.xyz", "forward_endpoint"),
        ("backward_last.xyz", "backward_endpoint"),
        ("forward_first.pdb", "forward_endpoint_pdb"),
        ("backward_last.pdb", "backward_endpoint_pdb"),
        ("forward_first.cif", "forward_endpoint_cif"),
        ("backward_last.cif", "backward_endpoint_cif"),
    )
    files: Dict[str, str] = {}
    for filename, key in specs:
        path = _irc_output_path(eulerpc, filename)
        if path.exists():
            files[key] = path.name
    return files


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
    help="Input structure file (.pdb, .cif, .mmcif, .xyz, _trj.xyz, etc.).",
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
              help="Net charge of the ML region/model system; overrides calc.model_charge from YAML. "
                   "Required unless --ligand-charge is provided.")
@click.option("-l", "--ligand-charge", type=str, default=None, show_default=False,
              help="Total charge for unknown ligand residues or a per-resname mapping "
                   "(e.g., GPP:-3,SAM:1), used to derive the ML-region charge when -q "
                   "is omitted (requires PDB/mmCIF input or --ref-pdb).")
@click.option(
    "-m",
    "--multiplicity",
    "spin",
    type=int,
    default=None,
    show_default="1",
    help="Spin multiplicity (2S+1); overrides calc.model_mult from YAML.",
)
@click.option(
    "--max-cycles", type=int, default=None, show_default="125", help="Maximum number of IRC steps; overrides irc.max_cycles from YAML."
)
@click.option("--step-size", type=float, default=None, show_default="0.10", help="Step length in Bohr (unweighted Cartesian coordinates). Default: 0.10 Bohr. Overrides irc.step_length from YAML.")
@click.option("--root", type=int, default=None, show_default="0", help="Imaginary mode index used for the initial displacement; overrides irc.root from YAML.")
@click.option(
    "--forward/--no-forward",
    "forward",
    default=None,
    show_default="forward",
    help="Run the forward IRC; overrides irc.forward from YAML.",
)
@click.option(
    "--backward/--no-backward",
    "backward",
    default=None,
    show_default="backward",
    help="Run the backward IRC; overrides irc.backward from YAML.",
)
@click.option(
    "--never-stop/--no-never-stop",
    "never_stop",
    default=None,
    show_default="no-never-stop",
    help=(
        "Ignore RMS-gradient, hard-gradient, energy-rise, and energy-change "
        "stops and trace until max-cycles. Numerical/integration failures and "
        "external interruption still stop the run; default off."
    ),
)
@click.option("-o", "--out-dir", type=str, default=IRC_KW["out_dir"], show_default=True, help="Output directory; overrides irc.out_dir from YAML.")
@click.option(
    "--hessian-calc-mode",
    type=click.Choice(["Analytical", "FiniteDifference"], case_sensitive=False),
    default=None, show_default="FiniteDifference",
    help=("How the ML backend builds the Hessian (Analytical or "
          "FiniteDifference); overrides calc.hessian_calc_mode from YAML. "
          "Default: 'FiniteDifference'. Runtime and memory depend on the "
          "backend and system; compare both modes on a representative pilot."),
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
    help="Reference PDB/mmCIF topology to use when --input is XYZ (keeps XYZ coordinates).",
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
    show_default="uma",
    help="ML backend for the ONIOM high-level region (default: uma).",
)
@click.option(
    "--embedcharge/--no-embedcharge",
    "embedcharge",
    default=False,
    show_default=True,
    help="Unavailable in v0.3.3; retained so older commands fail with an actionable diagnostic.",
)
@click.option(
    "--embedcharge-cutoff",
    "embedcharge_cutoff",
    type=float,
    default=None,
    show_default="12.0",
    help="Unavailable in v0.3.3 together with the retired electronic-embedding path.",
)
@click.option(
    "--link-atom-method",
    "link_atom_method",
    type=click.Choice(["scaled", "fixed"], case_sensitive=False),
    default=None,
    show_default="scaled",
    help="Link-atom position mode: scaled (g-factor, default) or fixed (legacy 1.09/1.01 Å).",
)
@click.option(
    "--mm-backend",
    "mm_backend",
    type=click.Choice(["hessian_ff", "openmm"], case_sensitive=False),
    default=None,
    show_default="hessian_ff",
    help="MM backend (default: hessian_ff). MM Hessians use finite differences by default; set calc.mm_fd: false for the hessian_ff analytical path.",
)
@click.option(
    "--cmap/--no-cmap",
    "use_cmap",
    default=None,
    show_default="cmap",
    help="Preserve CMAP terms in both real and model MM layers. Default: enabled when present in parm7.",
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
    show_default="None",
    help="Read an identified initial Hessian from 'mlmm freq --dump-hess'. "
         "Geometry, atom order, active-DOF basis, charge, and multiplicity must "
         "match; the file takes priority over hessian_cache and fresh computation.",
)
@click.option(
    "--allow-unverified-hess-state/--no-allow-unverified-hess-state",
    "allow_unverified_hess_state",
    default=False,
    show_default=True,
    help="Allow a schema-1 Hessian file whose charge and multiplicity cannot be "
         "verified. Use only after independently checking the electronic state.",
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
@add_workers_options()
@add_backend_model_option()
@add_calc_file_option()
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
    never_stop: Optional[bool],
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
    allow_unverified_hess_state: bool,
    out_json: bool,
    precision: Optional[str],
    workers: Optional[int],
    workers_per_node: Optional[int],
    backend_model: Optional[str],
    calc_file: Optional[str],
    calc_factory: Optional[str],
    irc_pos_def: Optional[bool],
) -> None:
    set_convert_file_enabled(convert_files)
    _is_param_explicit = make_is_param_explicit(ctx)
    if allow_unverified_hess_state and not read_hess:
        raise click.UsageError(
            "--allow-unverified-hess-state requires --read-hess."
        )

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
        model_pdb=model_pdb,
        model_indices_spec=model_indices_str,
        detect_layer=detect_layer,
        yaml_cfg=merged_yaml_cfg,
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
    error_out_dir = Path(out_dir).resolve()
    try:
        time_start = time.perf_counter()

        config_layer_cfg = load_yaml_dict(config_yaml)
        override_layer_cfg = load_yaml_dict(override_yaml)

        geom_cfg: Dict[str, Any] = dict(GEOM_KW_DEFAULT)
        calc_cfg: Dict[str, Any] = dict(CALC_KW_DEFAULT)
        irc_cfg: Dict[str, Any] = dict(IRC_KW_DEFAULT)
        # Keep the command-level detect-layer default unless YAML or an explicit
        # CLI option overrides it.
        calc_cfg["use_bfactor_layers"] = bool(detect_layer)

        apply_yaml_overrides(
            config_layer_cfg,
            [
                (geom_cfg, (("geom",),)),
                (calc_cfg, (("calc",), ("mlmm",))),
                (irc_cfg, (("irc",),)),
            ],
        )
        error_out_dir = Path(irc_cfg["out_dir"]).resolve()

        # CLI explicit overrides (after config YAML, before override YAML)
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
        if link_atom_method is not None:
            calc_cfg["link_atom_method"] = str(link_atom_method).lower()
        if mm_backend is not None:
            calc_cfg["mm_backend"] = str(mm_backend).lower()
        if use_cmap is not None:
            calc_cfg["use_cmap"] = use_cmap

        if _is_param_explicit("hessian_calc_mode") and hessian_calc_mode is not None:
            calc_cfg["hessian_calc_mode"] = str(hessian_calc_mode)
        # Validate again after the explicit mode override; the routing pass
        # above primarily installs worker counts and validates YAML values.
        apply_workers_to_calc_cfg(calc_cfg, None, None)
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
        if _is_param_explicit("never_stop") and never_stop is not None:
            irc_cfg["never_stop"] = bool(never_stop)
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
        error_out_dir = Path(irc_cfg["out_dir"]).resolve()
        from pysisyphus.tr_projection import normalize_tr_projection_mode
        geom_cfg["tr_projection"] = normalize_tr_projection_mode(
            geom_cfg.get("tr_projection")
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
        _validate_irc_directions(irc_cfg)
        if not calc_cfg.get("real_parm7"):
            raise click.BadParameter(
                "Missing --parm (or calc.real_parm7 in YAML).; "
                "recover: pass --parm /path/to/real.parm7 OR add 'calc:\\n  real_parm7: ...' to --config YAML."
            )

        out_dir_path = Path(irc_cfg["out_dir"]).resolve()
        layer_source_pdb = source_path
        detect_layer_enabled = bool(calc_cfg.get("use_bfactor_layers", True))
        model_pdb_cfg = calc_cfg.get("model_pdb")

        from mlmm.core.embedcharge_policy import reject_retired_embedcharge_cli

        reject_retired_embedcharge_cli(
            calc_cfg,
            cutoff_requested=_is_param_explicit("embedcharge_cutoff"),
        )

        if show_config:
            click.echo(
                pretty_block(
                    "yaml_layers",
                    {
                        "config": None if config_yaml is None else str(config_yaml),
                        "override_yaml": None if override_yaml is None else str(override_yaml),
                        "merged_keys": sorted(merged_yaml_cfg.keys()),
                    },
                force=True)
            )

        if dry_run:
            if model_pdb_cfg is not None:
                model_region_source = "model_pdb"
            elif model_indices:
                model_region_source = "model_indices"
            elif detect_layer_enabled:
                model_region_source = "bfactor"
            else:
                raise click.BadParameter("Provide --model-pdb or --model-indices when B-factor layer detection is disabled in the configuration.")
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
                        "tr_projection": geom_cfg["tr_projection"],
                        "will_run_irc": True,
                        "will_write_trajectories": True,
                        "backend": calc_cfg.get("backend", "uma"),
                        "embedcharge": bool(calc_cfg.get("embedcharge", False)),
                    },
                )
            )
            click.echo("[dry-run] Validation complete. IRC execution was skipped.")
            emit(
                format_elapsed("[time] Elapsed Time for IRC", time_start),
                narrative=True,
            )
            return

        irc_protected_inputs = (
            input_path,
            prepared_input.source_path,
            geom_input_path,
            prepared_input.original_path,
            ref_pdb,
            (
                Path(calc_cfg["input_pdb"])
                if calc_cfg.get("input_pdb")
                else None
            ),
            config_yaml,
            override_yaml,
            Path(calc_cfg["real_parm7"]) if calc_cfg.get("real_parm7") else None,
            Path(model_pdb_cfg) if model_pdb_cfg else None,
            (
                Path(calc_cfg["calc_file"])
                if calc_cfg.get("calc_file")
                else None
            ),
            read_hess,
        )
        out_dir_path = _prepare_irc_output_dir(
            out_dir_path,
            prefix=str(irc_cfg.get("prefix") or ""),
            protected_inputs=irc_protected_inputs,
        )

        if detect_layer_enabled and layer_source_pdb.suffix.lower() != ".pdb":
            raise click.BadParameter("--detect-layer requires a PDB input (or --ref-pdb).")

        try:
            model_pdb_path, layer_info = resolve_ml_layer_assignment(
                source_path=layer_source_pdb,
                out_dir_path=out_dir_path,
                model_pdb=model_pdb_cfg,
                model_indices=model_indices,
                detect_layer=detect_layer_enabled,
                hess_cutoff=calc_cfg.get("hess_cutoff"),
                movable_cutoff=calc_cfg.get("movable_cutoff"),
                calc_cfg=calc_cfg,
                protected_inputs=irc_protected_inputs,
                echo_fn=click.echo,
            )
        except click.ClickException as exc:
            raise click.BadParameter(exc.message) from exc
        _ = apply_layer_freeze_constraints(
            geom_cfg,
            calc_cfg,
            layer_info,
            echo_fn=click.echo,
        )
        _align_three_layer_hessian_targets(calc_cfg, echo_fn=click.echo)

        # Default-verbosity entry summary (skipped in child mode).
        from mlmm.core.utils import calculator_run_label, echo_run_summary
        echo_run_summary({
            "input": str(input_path),
            "backend": calculator_run_label(calc_cfg),
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

        def _current_hessian_active_dofs() -> np.ndarray:
            """Return the exact Cartesian basis produced by the current calculator."""
            full_n_dof = int(geometry.cart_coords.size)
            core = getattr(calc, "core", None)
            if core is None or not bool(getattr(core, "return_partial_hessian", False)):
                return np.arange(full_n_dof, dtype=np.int64)
            active_atoms = np.asarray(
                getattr(core, "hess_active_atoms", []), dtype=np.int64
            ).reshape(-1)
            if active_atoms.size == 0:
                raise click.ClickException(
                    "Current calculator resolved an empty Hessian active-atom basis."
                )
            return np.concatenate(
                [3 * active_atoms + axis for axis in range(3)]
            ).reshape(3, -1).T.reshape(-1)

        _expected_hessian_dofs = _current_hessian_active_dofs()

        echo_resolved_device()

        # Seed the initial Hessian.
        # Priority: --read-hess file > hessian_cache > fresh computation.
        from mlmm.io.hessian_cache import (
            discard as _hess_discard,
            load_matching as _hess_load_matching,
            store as _hess_store,
            identity_from_context as _hess_identity,
        )
        # calibrated GPU-first IRC Hessian/integration device policy.
        # ``auto`` keeps the resolved backend device (GPU-first when a CUDA
        # device is present), so the integration Hessian and large tensors stay
        # GPU-resident by default. An explicit ``cuda`` request stays on CUDA or
        # errors -- it is never silently moved to CPU. ``cpu`` is an explicit,
        # calibrated offload for large unfrozen systems, not a silent fallback.
        from mlmm.workflows._microiteration import resolve_hessian_device
        _requested_hess_device = (hess_device or "auto").strip().lower()
        if _requested_hess_device == "auto":
            _hess_dev = _torch_device(calc_cfg.get("ml_device", "auto"))
            _hess_dev_reason = "auto_gpu_first" if _hess_dev.type == "cuda" else "auto_cpu"
        else:
            try:
                _eff_hess_device, _hess_dev_reason = resolve_hessian_device(
                    _requested_hess_device, torch.cuda.is_available()
                )
            except ValueError as exc:
                raise click.ClickException(str(exc)) from exc
            _hess_dev = _torch_device(_eff_hess_device)
        click.echo(
            f"[device] IRC Hessian device: requested={_requested_hess_device}, "
            f"effective={_hess_dev.type} ({_hess_dev_reason})."
        )
        if _hess_dev.type == "cpu":
            click.echo("[device] Hessian operations will run on CPU.")

        if read_hess:
            _initial_hessian_source = "file"
            click.echo(f"[irc] Loading initial Hessian from {read_hess}")
            from mlmm.io.hessian_file import load_hessian_file
            from mlmm.io.hessian_cache import persistent_identity_from_context

            try:
                _loaded_hessian = load_hessian_file(
                    read_hess,
                    cart_coords_bohr=geometry.cart_coords,
                    atomic_numbers=geometry.atomic_numbers,
                    expected_model_charge=int(calc_cfg["model_charge"]),
                    expected_model_mult=int(calc_cfg["model_mult"]),
                    expected_potential_identity=persistent_identity_from_context(
                        geometry,
                        calc_cfg,
                    ),
                    expected_active_dofs=_expected_hessian_dofs,
                    allow_unverified_state=allow_unverified_hess_state,
                    allow_unverified_pes=allow_unverified_hess_state,
                )
            except ValueError as exc:
                raise click.ClickException(str(exc)) from exc
            h_init = torch.as_tensor(
                _loaded_hessian["hessian"], dtype=torch.float64, device=_hess_dev
            )
            # Restore partial-Hessian metadata if freq --dump-hess saved it,
            # so a partial Hessian (active_n_dof != 3N) is consumed correctly
            # instead of tripping the Geometry cart_hessian shape assertion.
            _partial_metadata = _loaded_hessian["partial_metadata"]
            _hessian_state_verified = bool(
                _loaded_hessian["electronic_state_verified"]
            )
            _hessian_pes_verified = bool(
                _loaded_hessian["potential_identity_verified"]
            )
            _hessian_file_schema = int(_loaded_hessian["schema_version"])
            import hashlib

            _hessian_hasher = hashlib.sha256()
            with Path(read_hess).open("rb") as _hessian_stream:
                for _hessian_block in iter(
                    lambda: _hessian_stream.read(1 << 20),
                    b"",
                ):
                    _hessian_hasher.update(_hessian_block)
            _hessian_file_sha256 = _hessian_hasher.hexdigest()
            if not _hessian_state_verified:
                click.echo(
                    "[irc] WARNING: the schema-1 Hessian does not identify "
                    "charge or multiplicity; proceeding by explicit opt-in.",
                    err=True,
                )
            if not _hessian_pes_verified:
                click.echo(
                    "[irc] WARNING: the legacy Hessian does not identify its "
                    "generating PES; proceeding by explicit opt-in.",
                    err=True,
                )
            if _partial_metadata is not None:
                geometry.within_partial_hessian = dict(_partial_metadata)
                click.echo(
                    f"[irc] Restored partial-Hessian metadata from npz "
                    f"(active_n_dof={_partial_metadata['active_n_dof']})."
                )
            del _loaded_hessian
        else:
            _hessian_state_verified = True
            _hessian_pes_verified = True
            _hessian_file_schema = None
            _hessian_file_sha256 = None
            # reuse the tsopt TS Hessian only on a full evaluation-identity
            # match; the all workflow may round-trip the TS through a
            # three-decimal PDB, so the coordinate field keeps the wider bohr
            # tolerance.  The layer-specific active-DOF basis check below is an
            # additional guard that identity matching does not replace.
            cached = _hess_load_matching(
                "ts",
                _hess_identity(geometry, calc_cfg, role="ts"),
                atol=1.1e-3,
            )
            if cached is not None:
                _cached_dofs = cached.get("active_dofs")
                if _cached_dofs is None:
                    _cached_dofs = np.arange(geometry.cart_coords.size, dtype=np.int64)
                if not np.array_equal(
                    np.asarray(_cached_dofs, dtype=np.int64).reshape(-1),
                    _expected_hessian_dofs,
                ):
                    click.echo(
                        "[irc] Cached TS Hessian active-DOF basis does not match "
                        "the current layer selection; calculating a fresh Hessian.",
                        err=True,
                    )
                    cached = None
            if cached is not None:
                _initial_hessian_source = "cache"
                emit("[irc] Reusing cached TS Hessian from tsopt.", narrative=True)
                active_dofs = cached.get("active_dofs")
                h_raw = cached["hessian"]
                if isinstance(h_raw, torch.Tensor):
                    h_init = h_raw.to(device=_hess_dev)
                else:
                    h_init = torch.as_tensor(
                        h_raw, dtype=torch.float64, device=_hess_dev
                    )
                if active_dofs is not None:
                    geometry.within_partial_hessian = {
                        "active_n_dof": len(active_dofs),
                        "full_n_dof": geometry.cart_coords.size,
                        "active_dofs": active_dofs,
                        "active_atoms": sorted(set(d // 3 for d in active_dofs)),
                    }
            else:
                _initial_hessian_source = "fresh"
                click.echo("[irc] Seeding initial Hessian via shared freq backend.")
                h_init, _ = _calc_full_hessian_torch(
                    geometry,
                    calc_cfg,
                    _hess_dev,
                    refresh_geom_meta=True,
                )

        # preserve the declared ordered active basis (ML + MovableMM).
        # A promised ML+MovableMM IRC active space must NOT be silently cropped
        # to ML + link-parent DOFs: that drops every MovableMM coordinate from
        # the physical IRC path, so the bundled integrator writes zero
        # displacement into those atoms. VRAM pressure on a large dense Hessian
        # is managed by the explicit, logged --hess-device policy (GPU-first by
        # default; `cpu` is a calibrated user offload), never by a silent basis
        # change. We only VALIDATE the seeded matrix against the already-declared
        # ordered basis here; we never reduce it or rebuild its active map.
        _full_n_dof = int(geometry.cart_coords.size)
        _seeded_n = int(h_init.shape[0])
        _within = getattr(geometry, "within_partial_hessian", None)
        if _within is not None and _within.get("active_dofs") is not None:
            _declared = np.asarray(_within["active_dofs"], dtype=np.int64).reshape(-1)
            if set(_declared.tolist()) != set(int(d) for d in _expected_hessian_dofs.tolist()):
                raise click.ClickException(
                    "Seeded Hessian active-DOF basis does not match the declared "
                    "layer selection; refusing to run IRC on an inconsistent basis."
                )
            _expected_n = int(_declared.size)
        else:
            _expected_n = int(_expected_hessian_dofs.size)
        if _seeded_n not in (_expected_n, _full_n_dof):
            raise click.ClickException(
                f"Seeded Hessian dimension {_seeded_n} matches neither the "
                f"declared active basis ({_expected_n}) nor the full Cartesian "
                f"space ({_full_n_dof}); refusing to run IRC on a cropped basis."
            )
        del _expected_hessian_dofs

        geometry.cart_hessian = h_init
        click.echo(f"[irc] Initial Hessian seeded (shape={h_init.shape[0]}x{h_init.shape[1]}).")
        del h_init

        eulerpc = EulerPC(geometry, **irc_cfg)
        from pysisyphus.tr_projection import active_tr_basis
        _basis, _rigid_info = active_tr_basis(
            torch.as_tensor(geometry.coords3d, dtype=torch.float64),
            torch.as_tensor(geometry.masses, dtype=torch.float64),
            eulerpc._act_atoms,
            mode=geometry.tr_projection,
        )
        del _basis
        eulerpc.rigid_projection_info = _rigid_info
        click.echo(
            "[irc] Rigid projection: "
            f"treatment={_rigid_info.treatment}, "
            f"rank={_rigid_info.effective_rank}, "
            f"full_rigid_rank={_rigid_info.full_rigid_rank}."
        )

        # A failed or one-sided IRC must not expose an endpoint Hessian left by
        # an earlier segment in this process.
        _hess_discard("irc_left")
        _hess_discard("irc_right")
        eulerpc.run()

        quick_directions = []
        for direction in ("forward", "backward"):
            if not getattr(eulerpc, direction, False):
                continue
            n_frames = len(getattr(eulerpc, f"{direction}_energies", []))
            last_cycle = getattr(eulerpc, f"{direction}_cycle", None)
            reached_cycle_cap = (
                last_cycle is not None
                and int(last_cycle) + 1 >= int(eulerpc.max_cycles)
            )
            if 0 < n_frames <= 3 and not reached_cycle_cap:
                quick_directions.append(direction)
        if quick_directions:
            warning = (
                "[irc] IRC stopped after only a few frames in "
                + ", ".join(quick_directions)
                + ". Retry with a smaller maximum step, for example "
                "--step-size 0.05."
            )
            if not eulerpc.never_stop:
                warning += (
                    " If a small uphill/flat section is intentional, also "
                    "consider --never-stop; it is opt-in."
                )
            click.echo(
                warning,
                err=True,
            )

        # Cache IRC endpoint Hessians (Bofill-updated mw → Cartesian)
        def _unmw_and_store(mw_H, key, endpoint_cart_coords, direction):
            """Un-mass-weight active-DOF Hessian on device, store partial on CPU."""
            act = eulerpc._act_dofs
            m_sqrt = geometry.masses_rep ** 0.5
            ms_act = m_sqrt[act]
            H_cart_act_np = _consume_mw_hessian_to_cartesian_active(
                mw_H, ms_act
            )
            _hess_store(
                key,
                H_cart_act_np,
                active_dofs=list(act),
                meta={
                    "cart_coords": endpoint_cart_coords,
                    "irc_direction": direction,
                },
                identity=_hess_identity(
                    geometry,
                    calc_cfg,
                    role=key,
                    cart_coords=endpoint_cart_coords,
                ),
            )

        # cache an endpoint Hessian ONLY for a requested direction that
        # explicitly CONVERGED. A nonconverged (max-cycle) direction may still
        # carry a Bofill-updated Hessian, but promoting it would let a
        # nonconverged endpoint seed a downstream RFO as if it were a real
        # minimum. The keys were discarded before eulerpc.run; keep them
        # discarded when the direction did not converge so no stale/never-
        # converged Hessian is reused.
        from mlmm.workflows._outcomes import (
            irc_hessian_cache_eligible as _irc_hess_eligible,
        )
        _fwd_conv = _irc_hess_eligible(eulerpc, "forward_is_converged")
        _bwd_conv = _irc_hess_eligible(eulerpc, "backward_is_converged")
        if (
            eulerpc.forward
            and _fwd_conv
            and getattr(eulerpc, "forward_mw_hessian", None) is not None
        ):
            forward_mw_hessian = eulerpc.forward_mw_hessian
            eulerpc.forward_mw_hessian = None
            forward_endpoint = (
                np.asarray(eulerpc.forward_mw_coords[0], dtype=float)
                / np.asarray(eulerpc.m_sqrt, dtype=float)
            )
            try:
                _unmw_and_store(
                    forward_mw_hessian,
                    "irc_left",
                    forward_endpoint,
                    "forward",
                )
            finally:
                del forward_mw_hessian
                if torch.cuda.is_available():
                    torch.cuda.empty_cache()
            click.echo("[irc] Cached forward endpoint Hessian as 'irc_left'.")
        else:
            eulerpc.forward_mw_hessian = None
            _hess_discard("irc_left")
        if (
            eulerpc.backward
            and _bwd_conv
            and getattr(eulerpc, "mw_hessian", None) is not None
        ):
            backward_mw_hessian = eulerpc.mw_hessian
            eulerpc.mw_hessian = None
            backward_endpoint = (
                np.asarray(eulerpc.backward_mw_coords[-1], dtype=float)
                / np.asarray(eulerpc.m_sqrt, dtype=float)
            )
            try:
                _unmw_and_store(
                    backward_mw_hessian,
                    "irc_right",
                    backward_endpoint,
                    "backward",
                )
            finally:
                del backward_mw_hessian
                if torch.cuda.is_available():
                    torch.cuda.empty_cache()
            click.echo("[irc] Cached backward endpoint Hessian as 'irc_right'.")
        else:
            eulerpc.mw_hessian = None
            _hess_discard("irc_right")
        eulerpc.mw_hessian = None
        eulerpc.forward_mw_hessian = None
        if torch.cuda.is_available():
            torch.cuda.empty_cache()

        if source_path.suffix.lower() == ".pdb":
            ref_pdb_path = source_path.resolve()

            # Whole IRC trajectory
            for stem in ("finished", "forward", "backward"):
                _echo_convert_trj_to_pdb_if_exists(
                    _irc_output_path(eulerpc, f"{stem}_irc_trj.xyz"),
                    ref_pdb_path,
                    _irc_output_path(eulerpc, f"{stem}_irc.pdb"),
                )
            # Forward arrays are reversed for the stitched IRC, so the
            # direction-semantic endpoints are forward_first/backward_last.
            for tag in ("forward_first", "backward_last"):
                endpoint_xyz = _irc_output_path(eulerpc, f"{tag}.xyz")
                endpoint_pdb = _irc_output_path(eulerpc, f"{tag}.pdb")
                if (
                    is_convert_file_enabled()
                    and endpoint_xyz.exists()
                    and not endpoint_pdb.exists()
                ):
                    try:
                        convert_xyz_to_pdb(endpoint_xyz, ref_pdb_path, endpoint_pdb)
                        click.echo(f"[convert] Wrote '{endpoint_pdb}'.")
                    except Exception as e:
                        logger.debug("Failed to convert %s to PDB", endpoint_xyz.name, exc_info=True)
                        click.echo(f"[convert] WARNING: Failed to convert '{tag}.xyz' to PDB: {e}", err=True)

        if out_json:
            from mlmm.core.utils import calculator_provenance, write_result_json
            _all_e = eulerpc.all_energies
            _n_fwd = len(getattr(eulerpc, "forward_energies", [])) if hasattr(eulerpc, "forward_energies") else 0
            _n_bwd = len(getattr(eulerpc, "backward_energies", [])) if hasattr(eulerpc, "backward_energies") else 0
            _ts_e = float(eulerpc.ts_energy) if hasattr(eulerpc, "ts_energy") else None
            _irc_files = _collect_irc_output_files(eulerpc)
            result_data = {
                "status": "completed",
                "n_frames_forward": _n_fwd,
                "n_frames_backward": _n_bwd,
                "n_frames_total": len(_all_e),
                "forward_converged": getattr(eulerpc, 'forward_is_converged', None),
                "backward_converged": getattr(eulerpc, 'backward_is_converged', None),
                **calculator_provenance(calc_cfg),
                "charge": calc_cfg.get("model_charge"),
                "spin": calc_cfg.get("model_mult"),
                "n_freeze_atoms": len(geom_cfg.get("freeze_atoms", [])),
                "step_length": irc_cfg.get("step_length"),
                "max_cycles": irc_cfg.get("max_cycles"),
                "never_stop": bool(irc_cfg.get("never_stop", False)),
                "never_stop_energy_bypasses": int(
                    getattr(eulerpc, "never_stop_energy_increase_bypasses", 0)
                    + getattr(eulerpc, "never_stop_energy_convergence_bypasses", 0)
                ),
                "rigid_projection": {
                    **_rigid_info.as_dict(),
                    "hessian_space": (
                        "active" if len(eulerpc._act_atoms) < len(geometry.atoms) else "full"
                    ),
                    "hessian_shape": list(eulerpc.init_hessian_shape),
                    "hessian_source": _initial_hessian_source,
                    "electronic_state_verified": _hessian_state_verified,
                    "pes_identity_verified": _hessian_pes_verified,
                    "hessian_file_schema": _hessian_file_schema,
                    "hessian_file_sha256": _hessian_file_sha256,
                    "hessian_representation": "cartesian-unweighted-unprojected",
                },
                "input_file": str(source_path),
                "files": _irc_files,
            }
            result_data.update(_directional_endpoint_energy_fields(_all_e, _ts_e))

            # one LeafOutcome per requested IRC direction. A
            # requested direction is usable only when it explicitly converged; a
            # disabled direction is optional (not a failure). Legacy ``status``
            # stays "completed" (the IRC process ran).
            from mlmm.workflows._outcomes import (
                aggregate_workflow_truth as _agg_truth,
                attach_outcomes as _attach,
                irc_direction_leaves as _irc_dir_leaves,
            )
            _dir_leaves, _dir_expected = _irc_dir_leaves(
                (
                    (
                        "forward",
                        bool(getattr(eulerpc, "forward", False)),
                        getattr(eulerpc, "forward_is_converged", None),
                        _n_fwd,
                        [_irc_files["forward_irc"]] if "forward_irc" in _irc_files else [],
                    ),
                    (
                        "backward",
                        bool(getattr(eulerpc, "backward", False)),
                        getattr(eulerpc, "backward_is_converged", None),
                        _n_bwd,
                        [_irc_files["backward_irc"]] if "backward_irc" in _irc_files else [],
                    ),
                )
            )
            _attach(
                result_data,
                truth=_agg_truth(_dir_leaves, _dir_expected),
                stage_outcomes=_dir_leaves,
            )

            # Bond changes between IRC endpoints
            try:
                from mlmm.domain.bond_changes import compare_structures
                _irc_first_xyz = _irc_output_path(
                    eulerpc, "finished_first.xyz"
                )
                _irc_last_xyz = _irc_output_path(
                    eulerpc, "finished_last.xyz"
                )
                if _irc_first_xyz.exists() and _irc_last_xyz.exists():
                    _g1 = geom_loader(str(_irc_first_xyz))
                    _g2 = geom_loader(str(_irc_last_xyz))
                    _bc = compare_structures(_g1, _g2, device="cpu")
                    _elems = [a.capitalize() for a in _g1.atoms]
                    result_data["bond_changes"] = {
                        "formed": [f"{_elems[i]}{i+1}-{_elems[j]}{j+1}" for i, j in sorted(_bc.formed_covalent)],
                        "broken": [f"{_elems[i]}{i+1}-{_elems[j]}{j+1}" for i, j in sorted(_bc.broken_covalent)],
                    }
                    result_data["bond_changes_direction"] = (
                        "finished_first_to_finished_last"
                    )
            except Exception:
                logger.debug("irc: bond-changes enrichment skipped", exc_info=True)

            write_result_json(
                out_dir_path, result_data,
                command="irc",
                elapsed_seconds=time.perf_counter() - time_start,
            )

        # summary.md and key_* outputs are disabled.
        emit(
            format_elapsed("[time] Elapsed Time for IRC", time_start),
            narrative=True,
        )

    except KeyboardInterrupt:
        click.echo("\nInterrupted by user.", err=True)
        sys.exit(130)
    except _IRCOutputCollisionError:
        raise
    except click.BadParameter as e:
        _write_error_json(
            error_out_dir, "irc", e, "BadParameter", time_start
        )
        raise
    except Exception as e:
        render_cli_exception(
            e,
            label="IRC",
            out_dir=error_out_dir,
            command="irc",
            time_start=time_start,
        )
    finally:
        prepared_input.cleanup()
        # Drop local references before collecting cycles held by calculator and
        # torch module objects.
        calc = eulerpc = geometry = None
        del calc, eulerpc, geometry
        gc.collect()  # break cyclic refs inside torch.nn.Module
        if torch.cuda.is_available():
            torch.cuda.empty_cache()


# Script entry point
if __name__ == "__main__":
    cli()
