"""
Single-point ML/MM ONIOM energy / forces (and optional Hessian) calculation.

Example:
    mlmm sp -i structure.pdb --real-parm7 real.parm7 -q 0 -m 1
    mlmm sp -i structure.pdb --real-parm7 real.parm7 -q 0 --hess

For detailed documentation, see: docs/sp.md
"""
# DOMAIN_PURE

from __future__ import annotations

import gc
import logging
import time
from pathlib import Path
from typing import Optional

import click
import numpy as np
import torch
import yaml

from pysisyphus.helpers import geom_loader
from pysisyphus.constants import AU2EV

from mlmm.backends.mlmm_calc import mlmm
from mlmm.core.defaults import GEOM_KW_DEFAULT, MLMM_CALC_KW, OUT_DIR_SP
from mlmm.core.utils import (
    apply_yaml_overrides,
    build_model_pdb_from_bfactors,
    build_model_pdb_from_indices,
    format_elapsed,
    merge_freeze_atom_indices,
    parse_indices_string,
    prepare_input_structure,
    resolve_charge_spin_or_raise,
    set_convert_file_enabled,
)
from mlmm.cli.common_options import (
    add_ml_layer_detection_options,
    add_precision_option, add_backend_model_option,
    add_deterministic_option, add_allow_charge_mult_mismatch_option,
    add_print_every_option,
)
from mlmm.cli.decorators import (
    load_merged_yaml_cfg,
    make_is_param_explicit,
    render_cli_exception,
    resolve_yaml_sources,
)
logger = logging.getLogger(__name__)

EV2AU = 1.0 / AU2EV


@click.command(
    name="sp",
    short_help="Single-point ML/MM ONIOM energy / forces (and optional Hessian).",
)
@click.option(
    "-i", "--input", "input_path",
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
    required=True,
    help="Layered PDB (or XYZ) defining the ML/MM/Frozen system.",
)
@click.option(
    "--parm", "--real-parm7", "real_parm7",
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
    required=True,
    help="Amber parm7 of the full enzyme (canonical flag is --parm; --real-parm7 retained as alias).",
)
# ML-layer selection options (mirror opt.py / dft.py surface).
# These five options must remain on sp.py's decorator stack to match the
# cli() signature; without them every invocation fails with
# `TypeError: cli() missing 5 required positional arguments`.
@click.option(
    "--model-pdb", "model_pdb",
    type=click.Path(path_type=Path, exists=True, dir_okay=False),
    required=False, default=None,
    help="PDB defining atoms that belong to the ML (high-level) region. "
         "Optional when --detect-layer is enabled.",
)
@click.option(
    "--model-indices", "model_indices_str",
    type=str, default=None, show_default=False,
    help="Comma-separated atom indices for the ML region (ranges allowed like 1-5). "
         "Used when --model-pdb is omitted.",
)
@click.option(
    "--freeze-atoms", "freeze_atoms_cli",
    type=str, default=None, show_default=False,
    help="Comma-separated 1-based atom indices to freeze (e.g., '1,3,5').",
)
@click.option(
    "--radius-partial-hessian", "--hess-cutoff", "hess_cutoff",
    type=float, default=None, show_default=False,
    help="Distance cutoff (Å) from ML region for MM atoms to include in Hessian "
         "calculation. Applied to movable MM atoms; combinable with --detect-layer.",
)
@click.option(
    "--radius-freeze", "--movable-cutoff", "movable_cutoff",
    type=float, default=None, show_default=False,
    help="Distance cutoff (Å) from ML region for movable MM atoms. "
         "MM atoms beyond this are frozen.",
)
@click.option(
    "-q", "--charge", type=int, default=None,
    help="ML region total charge.",
)
@click.option(
    "-l", "--ligand-charge", "ligand_charge",
    type=str, default=None,
    help="Per-ligand charge mapping, e.g. 'SAM:1,GPP:-3'.",
)
@click.option(
    "-m", "--multiplicity", "spin", type=int, default=None,
    help="ML region spin multiplicity (2S+1).",
)
@click.option(
    "-o", "--out-dir", type=str, default=OUT_DIR_SP,
    show_default=True, help="Output directory.",
)
@click.option(
    "--hess/--no-hess", "do_hess", default=False, show_default=True,
    help="Also compute the full ONIOM Hessian and save to hessian.npy.",
)
@click.option(
    "--hessian-calc-mode", "hessian_calc_mode",
    type=click.Choice(["Analytical", "FiniteDifference"], case_sensitive=False),
    default=None, show_default=False,
    help="Hessian backend when --hess is set. Analytical only works for UMA; other backends fall back to FiniteDifference.",
)
@click.option(
    "--convert-files/--no-convert-files", "convert_files",
    default=True, show_default=True,
    help="Auto-convert output XYZ-like files into matching PDB beside them.",
)
@click.option(
    "--config", "config_yaml",
    type=click.Path(path_type=Path, exists=True, dir_okay=False),
    default=None, help="YAML config file with sections (calc:, geom:, …).",
)
@click.option(
    "--show-config/--no-show-config", "show_config",
    default=False, help="Print effective merged config and exit.",
)
@click.option(
    "--dry-run/--no-dry-run", "dry_run",
    default=False, help="Validate options and print the plan without running.",
)
@click.option(
    "--out-json/--no-out-json", "out_json",
    default=False, show_default=True,
    help="Write machine-readable result.json to out_dir.",
)
@click.option(
    "-b", "--backend",
    type=click.Choice(["uma", "orb", "mace", "aimnet2"], case_sensitive=False),
    default=None, show_default=False, help="ML backend for the ONIOM high-level region (default: uma).",
)
@click.option(
    "--embedcharge/--no-embedcharge", "embedcharge",
    default=False, show_default=True,
    help="Enable xTB point-charge embedding correction for MM->ML environmental effects.",
)
@click.option(
    "--embedcharge-cutoff", "embedcharge_cutoff",
    type=float, default=None, show_default=False,
    help="Embed-charge cutoff radius (Å) around the ML region.",
)
@click.option(
    "--link-atom-method", "link_atom_method",
    type=click.Choice(["scaled", "fixed"], case_sensitive=False),
    default=None, show_default=False,
    help="Link-atom positioning: scaled (g-factor) or fixed (1.09/1.01 Å).",
)
@click.option(
    "--mm-backend", "mm_backend",
    type=click.Choice(["hessian_ff", "openmm"], case_sensitive=False),
    default=None, show_default=False,
    help="MM backend (default: hessian_ff).",
)
@click.option(
    "--cmap/--no-cmap", "use_cmap",
    default=None, show_default=False,
    help="Enable CMAP (backbone cross-map) terms in model parm7.",
)
@click.option(
    "--use-cmap/--no-use-cmap", "use_cmap_legacy",
    default=None, show_default=False, hidden=True,
    help="Legacy alias for --cmap/--no-cmap. Prefer --cmap.",
)
@add_ml_layer_detection_options()
@add_precision_option()
@add_backend_model_option()
@add_deterministic_option()
@add_allow_charge_mult_mismatch_option()
@add_print_every_option()
@click.pass_context
def cli(
    ctx: click.Context,
    input_path: Path,
    real_parm7: Path,
    model_pdb: Optional[Path],
    model_indices_str: Optional[str],
    model_indices_one_based: bool,
    detect_layer: bool,
    freeze_atoms_cli: Optional[str],
    hess_cutoff: Optional[float],
    movable_cutoff: Optional[float],
    charge: Optional[int],
    ligand_charge: Optional[str],
    spin: Optional[int],
    out_dir: str,
    do_hess: bool,
    hessian_calc_mode: Optional[str],
    convert_files: bool,
    config_yaml: Optional[Path],
    show_config: bool,
    dry_run: bool,
    out_json: bool,
    backend: Optional[str],
    embedcharge: bool,
    embedcharge_cutoff: Optional[float],
    link_atom_method: Optional[str],
    mm_backend: Optional[str],
    use_cmap: Optional[bool],
    use_cmap_legacy: Optional[bool],
    precision: Optional[str],
    backend_model: Optional[str],
    print_every: Optional[int],
) -> None:
    """Compute a single-point ML/MM ONIOM energy + forces (and optionally Hessian)."""
    set_convert_file_enabled(convert_files)
    _is_param_explicit = make_is_param_explicit(ctx)
    # Legacy alias resolution: --use-cmap/--no-use-cmap → --cmap/--no-cmap.
    # If the canonical flag was not given, fall back to the legacy alias value.
    if use_cmap is None and use_cmap_legacy is not None:
        use_cmap = use_cmap_legacy

    config_yaml, _override, _legacy = resolve_yaml_sources(
        config_yaml=config_yaml, override_yaml=None, args_yaml_legacy=None,
    )
    merged_yaml_cfg, config_layer_cfg, _override_layer_cfg = load_merged_yaml_cfg(
        config_yaml=config_yaml, override_yaml=None,
    )

    prepared = prepare_input_structure(input_path)
    out_dir_path: Optional[Path] = None
    time_start: float = time.perf_counter()

    try:
        geom_cfg: dict = dict(GEOM_KW_DEFAULT)
        calc_cfg: dict = dict(MLMM_CALC_KW)
        sp_cfg: dict = {"out_dir": out_dir, "hess": False, "hessian_calc_mode": None}

        apply_yaml_overrides(
            config_layer_cfg,
            [
                (geom_cfg, (("geom",),)),
                (calc_cfg, (("calc",),)),
                (sp_cfg, (("sp",),)),
            ],
        )

        # Required parm7 / model selection
        calc_cfg["input_pdb"] = str(prepared.source_path)
        calc_cfg["real_parm7"] = str(real_parm7)
        if model_pdb is not None:
            calc_cfg["model_pdb"] = str(model_pdb)
        # Note: model_pdb auto-derive (from --model-indices or B-factor layers)
        # happens later after out_dir_path is mkdir'd — see comment "Auto-derive
        # model_pdb" before the mlmm(**calc_cfg) call.
        if _is_param_explicit("detect_layer"):
            calc_cfg["use_bfactor_layers"] = bool(detect_layer)
        if hess_cutoff is not None:
            calc_cfg["hess_cutoff"] = float(hess_cutoff)
        if movable_cutoff is not None:
            calc_cfg["movable_cutoff"] = float(movable_cutoff)
        if backend is not None:
            calc_cfg["backend"] = str(backend).lower()
        if _is_param_explicit("embedcharge"):
            calc_cfg["embedcharge"] = bool(embedcharge)
        if embedcharge_cutoff is not None:
            calc_cfg["embedcharge_cutoff"] = float(embedcharge_cutoff)
        if link_atom_method is not None:
            calc_cfg["link_atom_method"] = str(link_atom_method).lower()
        if mm_backend is not None:
            calc_cfg["mm_backend"] = str(mm_backend).lower()
        if use_cmap is not None:
            calc_cfg["use_cmap"] = bool(use_cmap)
        if _is_param_explicit("precision") and precision is not None:
            from mlmm.backends import apply_precision_to_calc_cfg, apply_backend_model_to_calc_cfg
            apply_precision_to_calc_cfg(calc_cfg, str(precision))
            apply_backend_model_to_calc_cfg(calc_cfg, backend_model)
        if _is_param_explicit("print_every") and print_every is not None:
            calc_cfg["print_every"] = int(print_every)

        # SP-specific CLI overrides
        if _is_param_explicit("out_dir"):
            sp_cfg["out_dir"] = out_dir
        if _is_param_explicit("do_hess"):
            sp_cfg["hess"] = bool(do_hess)
        if _is_param_explicit("hessian_calc_mode") and hessian_calc_mode is not None:
            sp_cfg["hessian_calc_mode"] = str(hessian_calc_mode)

        # Charge/spin resolution
        resolved_charge, resolved_spin = resolve_charge_spin_or_raise(
            prepared, charge, spin,
            ligand_charge=ligand_charge,
            prefix="[sp]",
        )
        calc_cfg["charge"] = int(resolved_charge)
        calc_cfg["spin"] = int(resolved_spin)

        out_dir_path = Path(sp_cfg["out_dir"]).resolve()

        if show_config:
            click.echo(yaml.safe_dump(
                {"calc": calc_cfg, "geom": geom_cfg, "sp": sp_cfg},
                sort_keys=False, allow_unicode=True,
            ).rstrip())
            # Help text already says "Print effective merged config and exit."
            # Honor that contract by returning before any SCF/Hessian work.
            return

        if dry_run:
            click.echo(f"[sp] dry-run: would compute ONIOM SP on {input_path} -> {out_dir_path}")
            return

        out_dir_path.mkdir(parents=True, exist_ok=True)

        # Apply optional CLI freeze_atoms extension
        if freeze_atoms_cli is not None:
            from mlmm.core.utils import _parse_freeze_atoms
            extra = _parse_freeze_atoms(freeze_atoms_cli)
            merge_freeze_atom_indices(geom_cfg, extra)

        coord_type = geom_cfg.get("coord_type", "cart")
        coord_kwargs = dict(geom_cfg)
        coord_kwargs.pop("coord_type", None)
        geom = geom_loader(prepared.geom_path, coord_type=coord_type, **coord_kwargs)

        click.echo(f"[sp] {input_path} -> {out_dir_path} (backend={calc_cfg.get('backend','uma')}, charge={calc_cfg['charge']}, spin={calc_cfg['spin']})")
        calc_cfg.setdefault("freeze_atoms", list(geom_cfg.get("freeze_atoms", [])))
        calc_cfg["return_partial_hessian"] = bool(sp_cfg["hess"])

        # Auto-derive model_pdb when not user-supplied. MLMMCore unconditionally
        # copies model_pdb to its tmpdir (mlmm_calc.py:1137), so it must be a
        # valid path. Mirror opt.py / tsopt.py / irc.py: B-factor layers first,
        # then --model-indices, fall back to input PDB itself for layered files.
        # NOTE: MLMM_CALC_KW seeds calc_cfg with "model_pdb": None, so use
        # explicit None check rather than `not in`.
        if calc_cfg.get("model_pdb") is None:
            _src_pdb = Path(prepared.source_path)
            if model_indices_str is not None:
                calc_cfg["model_indices_str"] = str(model_indices_str)
                calc_cfg["model_indices_one_based"] = bool(model_indices_one_based)
                # Parse --model-indices into an explicit atom list before
                # building the model PDB. Mirror opt.py / tsopt.py: without
                # this the ML region defaults to the entire structure because
                # write_model_pdb_from_indices([]) raises and we fall back to
                # the full input PDB.
                _model_indices = parse_indices_string(
                    str(model_indices_str), one_based=bool(model_indices_one_based)
                )
                try:
                    _model_pdb_path = build_model_pdb_from_indices(
                        _src_pdb, out_dir_path, _model_indices or []
                    )
                    calc_cfg["model_pdb"] = str(_model_pdb_path)
                except Exception:
                    calc_cfg["model_pdb"] = str(_src_pdb)
            else:
                try:
                    _model_pdb_path, _ = build_model_pdb_from_bfactors(
                        _src_pdb, out_dir_path
                    )
                    calc_cfg["model_pdb"] = str(_model_pdb_path)
                except Exception:
                    # Input PDB already layered serves as model_pdb fallback.
                    calc_cfg["model_pdb"] = str(_src_pdb)

        # Rename CLI-style keys to mlmm constructor kwargs to avoid duplicate-
        # value TypeError at super().__init__(charge=model_charge, ...).
        calc_cfg["model_charge"] = calc_cfg.pop("charge")
        calc_cfg["model_mult"] = calc_cfg.pop("spin")
        calc = mlmm(**calc_cfg)
        geom.set_calculator(calc)

        # Energy + forces
        t0 = time.perf_counter()
        energy_au = float(geom.energy)
        gradient_au = np.asarray(geom.gradient, dtype=float)
        forces_au = -gradient_au.reshape(-1, 3)
        elapsed_ef = time.perf_counter() - t0
        click.echo(f"[sp] energy = {energy_au:.10f} a.u.  |force|_max = {np.max(np.abs(forces_au)):.4e} a.u./bohr  ({elapsed_ef:.2f} s)")

        forces_path = out_dir_path / "forces.npy"
        np.save(forces_path, forces_au)

        # Optional Hessian
        hessian_path: Optional[Path] = None
        if sp_cfg["hess"]:
            mode = sp_cfg.get("hessian_calc_mode") or ("Analytical" if calc_cfg.get("backend", "uma") == "uma" else "FiniteDifference")
            click.echo(f"[sp] computing full ONIOM Hessian (mode={mode}) ...")
            t0 = time.perf_counter()
            # geom.hessian may be a CUDA torch.Tensor (UMA analytical path);
            # detach + cpu first so np.asarray doesn't trip the "can't convert
            # cuda:0 device type tensor to numpy" TypeError.
            _h_obj = geom.hessian
            if hasattr(_h_obj, "detach"):
                _h_obj = _h_obj.detach().cpu()
            H_np = np.asarray(_h_obj, dtype=float)
            elapsed_h = time.perf_counter() - t0
            hessian_path = out_dir_path / "hessian.npy"
            np.save(hessian_path, H_np)
            click.echo(f"[sp] Hessian {H_np.shape} written to {hessian_path}  ({elapsed_h:.2f} s)")

        elapsed_total = format_elapsed("[sp] Elapsed Time", time_start)
        summary = {
            "stage": "sp",
            "status": "ok",
            "input": str(prepared.source_path),
            "real_parm7": str(real_parm7),
            "backend": calc_cfg.get("backend", "uma"),
            # charge/spin were popped + renamed to model_charge/model_mult
            # before mlmm() construction (see "Rename CLI-style keys" block).
            "charge": calc_cfg.get("model_charge"),
            "spin": calc_cfg.get("model_mult"),
            "energy_au": energy_au,
            "forces_path": str(forces_path),
            "hessian_path": str(hessian_path) if hessian_path else None,
            "elapsed": elapsed_total,
        }
        if out_json:
            from mlmm.core.utils import write_result_json
            write_result_json(
                out_dir_path, summary,
                command="sp",
                elapsed_seconds=time.perf_counter() - time_start,
            )

        click.echo(f"[sp] {elapsed_total}.")

    except Exception as exc:
        render_cli_exception(
            exc, label="single-point",
            out_dir=out_dir_path or Path(sp_cfg.get("out_dir", OUT_DIR_SP)).resolve(),
            command="sp", time_start=time_start,
        )
    finally:
        try:
            del geom  # type: ignore[name-defined]
        except NameError:
            pass
        try:
            del calc  # type: ignore[name-defined]
        except NameError:
            pass
        gc.collect()
        if torch.cuda.is_available():
            torch.cuda.empty_cache()
