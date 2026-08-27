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
import tempfile
import time
from pathlib import Path
from typing import Any, Dict, Optional, Sequence

import click
import numpy as np
import torch
import yaml

from pysisyphus.helpers import geom_loader
from pysisyphus.constants import AU2EV

from mlmm.backends.mlmm_calc import mlmm
from mlmm.core.defaults import GEOM_KW_DEFAULT, MLMM_CALC_KW, OUT_DIR_SP
from mlmm.workflows.charge_prep import resolve_charge_spin_or_raise
from mlmm.workflows._opt_freq_common import (
    _convert_yaml_layer_atoms_1to0,
    _normalize_geom_freeze,
)
from mlmm.core.utils import (
    apply_yaml_overrides,
    apply_ref_pdb_override,
    calculator_provenance,
    emit_dry_run_complete,
    format_elapsed,
    merge_freeze_atom_indices,
    parse_indices_string,
    prepare_input_structure,
    read_bfactors_from_pdb,
    resolve_ml_layer_assignment,
)
from mlmm.cli.common_options import (
    add_ml_layer_detection_options,
    add_precision_option, add_backend_model_option, add_calc_file_option,
    add_workers_options,
    add_deterministic_option, add_allow_charge_mult_mismatch_option,
)
from mlmm.cli.decorators import (
    load_merged_yaml_cfg,
    make_is_param_explicit,
    render_cli_exception,
    resolve_yaml_sources,
)
logger = logging.getLogger(__name__)

EV2AU = 1.0 / AU2EV


class _SPOutputCollisionError(click.UsageError):
    """An SP output/input collision that bypasses envelope writing."""


def _reject_sp_output_collisions(
    out_dir: Path, protected_inputs: Sequence[Optional[Path]]
) -> None:
    destinations = [out_dir / name for name in ("result.json", "summary.json")]
    for protected in protected_inputs:
        if protected is None:
            continue
        source = Path(protected).expanduser().resolve()
        for destination in destinations:
            target = destination.resolve(strict=False)
            aliases = source == target
            if not aliases and source.exists() and destination.exists():
                aliases = source.samefile(destination)
            if aliases:
                raise _SPOutputCollisionError(
                    f"Input {protected} physically aliases reserved SP output {destination}."
                )


def _resolve_sp_ml_region(
    *,
    source_path: Path,
    out_dir_path: Path,
    calc_cfg: Dict[str, Any],
    model_indices_str: Optional[str],
    model_indices_one_based: bool,
    protected_inputs: Sequence[Optional[Path]],
) -> Dict[str, Any]:
    """Resolve one explicit/layered SP ML region and return its provenance."""

    atom_count = len(read_bfactors_from_pdb(source_path))
    if atom_count <= 0:
        raise click.ClickException(f"No atoms found in input PDB: {source_path}")

    configured_model = calc_cfg.get("model_pdb")
    model_indices = None
    if configured_model is None and model_indices_str is not None:
        model_indices = parse_indices_string(
            str(model_indices_str), one_based=bool(model_indices_one_based)
        )
        if not model_indices:
            raise click.BadParameter("--model-indices must select at least one atom.")
        invalid = [
            index for index in model_indices if index < 0 or index >= atom_count
        ]
        if invalid:
            raise click.BadParameter(
                "--model-indices contains indices outside the input atom bounds "
                f"[0, {atom_count}): {invalid}"
            )

    layer_detection_requested = bool(calc_cfg.get("use_bfactor_layers", True))
    if configured_model is not None:
        source = "model_pdb"
    elif model_indices is not None:
        source = "model_indices"
    else:
        source = "bfactor"

    model_pdb_path, layer_info = resolve_ml_layer_assignment(
        source_path=source_path,
        out_dir_path=out_dir_path,
        model_pdb=configured_model,
        model_indices=model_indices,
        detect_layer=layer_detection_requested,
        hess_cutoff=calc_cfg.get("hess_cutoff"),
        movable_cutoff=calc_cfg.get("movable_cutoff"),
        calc_cfg=calc_cfg,
        protected_inputs=protected_inputs,
        echo_fn=click.echo,
    )
    if source == "model_indices":
        region_count = len(model_indices or [])
    elif source == "bfactor" and layer_info is not None:
        region_count = len(layer_info.get("ml_indices", []))
    else:
        region_count = len(read_bfactors_from_pdb(model_pdb_path))
    if region_count <= 0:
        raise click.ClickException(
            f"The resolved ML region contains no atoms: {model_pdb_path}"
        )

    return {
        "ml_region_source": source,
        "ml_region_atom_count": int(region_count),
        "full_system_ml": bool(region_count == atom_count),
        "ml_region_model_pdb": str(model_pdb_path),
        "ml_region_indices": (
            [int(index) for index in model_indices]
            if source == "model_indices" and model_indices is not None
            else None
        ),
    }


@click.command(
    name="sp",
    short_help="Single-point ML/MM ONIOM energy / forces (and optional Hessian).",
)
@click.option(
    "-i", "--input", "input_path",
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
    required=True,
    help="Layered PDB/mmCIF, or XYZ with --ref-pdb, defining the ML/MM/Frozen system.",
)
@click.option(
    "--ref-pdb",
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
    default=None,
    help="Full-system PDB/mmCIF topology required when --input is XYZ.",
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
    help="ML-only, link-H-free PDB subset; atom identity/order must match the "
         "full PDB/parm7. When provided, it defines ML membership; "
         "--detect-layer still reads valid movable/frozen MM B-factors.",
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
    "--hess-cutoff", "hess_cutoff",
    type=float, default=None, show_default="all movable MM atoms",
    help="Distance cutoff (Å) from ML region for MM atoms to include in Hessian "
         "calculation. Applied to movable MM atoms; combinable with --detect-layer.",
)
@click.option(
    "--movable-cutoff", "movable_cutoff",
    type=float, default=None, show_default="use freeze_atoms",
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
    show_default="1",
    help="ML region spin multiplicity (2S+1).",
)
@click.option(
    "-o", "--out-dir", type=str, default=OUT_DIR_SP,
    show_default=True, help="Output directory.",
)
@click.option(
    "--hess/--no-hess", "do_hess", default=False, show_default=True,
    help="Also compute the active-coordinate ONIOM Hessian and save to hessian.npy.",
)
@click.option(
    "--hessian-calc-mode", "hessian_calc_mode",
    type=click.Choice(["Analytical", "FiniteDifference"], case_sensitive=False),
    default=None, show_default="FiniteDifference",
    help=(
        "Hessian backend when --hess is set. Analytical is supported by UMA, "
        "ORB, MACE, and AIMNet2; custom calculators use FiniteDifference. "
        "Analytical cannot be combined with --workers > 1."
    ),
)
@click.option(
    "--config", "config_yaml",
    type=click.Path(path_type=Path, exists=True, dir_okay=False),
    default=None, help="YAML config file with sections (calc:, geom:, …).",
)
@click.option(
    "--show-config/--no-show-config", "show_config",
    default=False, show_default=True, help="Print effective merged config and exit.",
)
@click.option(
    "--dry-run/--no-dry-run", "dry_run",
    default=False, show_default=True, help="Validate options and print the plan without running.",
)
@click.option(
    "--out-json/--no-out-json", "out_json",
    default=False, show_default=True,
    help="Write machine-readable result.json to out_dir.",
)
@click.option(
    "-b", "--backend",
    type=click.Choice(["uma", "orb", "mace", "aimnet2"], case_sensitive=False),
    default=None, show_default="uma", help="ML backend for the ONIOM high-level region.",
)
@click.option(
    "--embedcharge/--no-embedcharge", "embedcharge",
    default=False, show_default=True,
    help="Enable the experimental, computationally expensive xTB point-charge delta correction for MLIP/MM.",
)
@click.option(
    "--embedcharge-cutoff", "embedcharge_cutoff",
    type=float, default=None, show_default="12.0",
    help="Distance cutoff (Å) from the ML region for MM point charges used by the xTB delta correction.",
)
@click.option(
    "--link-atom-method", "link_atom_method",
    type=click.Choice(["scaled", "fixed"], case_sensitive=False),
    default=None, show_default="scaled",
    help="Link-atom positioning: scaled (g-factor) or fixed (1.09/1.01 Å).",
)
@click.option(
    "--mm-backend", "mm_backend",
    type=click.Choice(["hessian_ff", "openmm"], case_sensitive=False),
    default=None, show_default="hessian_ff",
    help="MM backend.",
)
@click.option(
    "--cmap/--no-cmap", "use_cmap",
    default=None, show_default="cmap",
    help="Preserve CMAP terms in both real and model MM layers when present in parm7.",
)
@add_ml_layer_detection_options()
@add_precision_option()
@add_workers_options()
@add_backend_model_option()
@add_calc_file_option()
@add_deterministic_option()
@add_allow_charge_mult_mismatch_option()
@click.pass_context
def cli(
    ctx: click.Context,
    input_path: Path,
    ref_pdb: Optional[Path],
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
    precision: Optional[str],
    workers: Optional[int],
    workers_per_node: Optional[int],
    backend_model: Optional[str],
    calc_file: Optional[str],
    calc_factory: Optional[str],
) -> None:
    """Compute a single-point ML/MM ONIOM energy + forces (and optionally Hessian)."""
    _is_param_explicit = make_is_param_explicit(ctx)

    config_yaml, _override, _legacy = resolve_yaml_sources(
        config_yaml=config_yaml, override_yaml=None, args_yaml_legacy=None,
    )
    merged_yaml_cfg, config_layer_cfg, _override_layer_cfg = load_merged_yaml_cfg(
        config_yaml=config_yaml, override_yaml=None,
    )

    if input_path.suffix.lower() == ".xyz" and ref_pdb is None:
        raise click.BadParameter(
            "XYZ input requires --ref-pdb topology.",
            param_hint="--ref-pdb",
        )
    prepared = prepare_input_structure(input_path)
    apply_ref_pdb_override(prepared, ref_pdb)
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
                (calc_cfg, (("calc",), ("mlmm",))),
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
        from mlmm.backends import apply_precision_to_calc_cfg
        from mlmm.backends import apply_workers_to_calc_cfg
        from mlmm.backends import apply_backend_model_to_calc_cfg
        apply_backend_model_to_calc_cfg(calc_cfg, backend_model)
        # --calc-file overrides --backend with a user ASE Calculator (custom backend).
        from mlmm.backends import apply_calc_file_to_calc_cfg
        apply_calc_file_to_calc_cfg(calc_cfg, calc_file, calc_factory)
        apply_precision_to_calc_cfg(
            calc_cfg,
            str(precision)
            if _is_param_explicit("precision") and precision is not None
            else None,
        )
        # SP-specific CLI overrides
        if _is_param_explicit("out_dir"):
            sp_cfg["out_dir"] = out_dir
        if _is_param_explicit("do_hess"):
            sp_cfg["hess"] = bool(do_hess)
        if _is_param_explicit("hessian_calc_mode") and hessian_calc_mode is not None:
            sp_cfg["hessian_calc_mode"] = str(hessian_calc_mode)
        if sp_cfg.get("hessian_calc_mode"):
            # ``geom.hessian`` reads the mode from the calculator, not from
            # the reporting-only SP config.
            calc_cfg["hessian_calc_mode"] = str(sp_cfg["hessian_calc_mode"])
        # Validate only after all YAML and CLI overrides have been resolved.
        try:
            apply_workers_to_calc_cfg(calc_cfg, workers, workers_per_node)
        except ValueError as exc:
            raise click.ClickException(str(exc)) from exc

        # Charge/spin resolution
        resolved_charge, resolved_spin = resolve_charge_spin_or_raise(
            prepared, charge, spin,
            ligand_charge=ligand_charge,
            prefix="[sp]",
            model_pdb=model_pdb,
            model_indices_spec=model_indices_str,
            detect_layer=detect_layer,
            yaml_cfg=merged_yaml_cfg,
        )
        calc_cfg["charge"] = int(resolved_charge)
        calc_cfg["spin"] = int(resolved_spin)

        out_dir_path = Path(sp_cfg["out_dir"]).resolve()
        geom_cfg["freeze_atoms"] = _normalize_geom_freeze(
            geom_cfg.get("freeze_atoms")
        )
        _convert_yaml_layer_atoms_1to0(calc_cfg)
        protected_inputs = (
            input_path,
            prepared.original_path,
            prepared.source_path,
            prepared.geom_path,
            ref_pdb,
            real_parm7,
            model_pdb,
            config_yaml,
            Path(calc_cfg["calc_file"]) if calc_cfg.get("calc_file") else None,
        )
        _reject_sp_output_collisions(out_dir_path, protected_inputs)

        if show_config:
            click.echo(yaml.safe_dump(
                {"calc": calc_cfg, "geom": geom_cfg, "sp": sp_cfg},
                sort_keys=False, allow_unicode=True,
            ).rstrip())
            if not dry_run:
                # Help text says show-config exits before any calculation.
                click.echo(format_elapsed("[time] Elapsed Time for SP", time_start))
                return

        # Validate the ML-region contract without creating the requested output
        # directory or constructing a calculator. The real run repeats this
        # resolution below and publishes its generated model PDB there.
        with tempfile.TemporaryDirectory(prefix="mlmm_sp_validate_") as tmp_dir:
            validation_cfg = dict(calc_cfg)
            _resolve_sp_ml_region(
                source_path=Path(prepared.source_path),
                out_dir_path=Path(tmp_dir),
                calc_cfg=validation_cfg,
                model_indices_str=model_indices_str,
                model_indices_one_based=model_indices_one_based,
                protected_inputs=(),
            )

        if dry_run:
            click.echo(f"[sp] dry-run: would compute ONIOM SP on {input_path} -> {out_dir_path}")
            emit_dry_run_complete()
            return

        out_dir_path.mkdir(parents=True, exist_ok=True)
        if not sp_cfg["hess"]:
            (out_dir_path / "hessian.npy").unlink(missing_ok=True)
        if not out_json:
            for name in ("result.json", "summary.json"):
                (out_dir_path / name).unlink(missing_ok=True)

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

        ml_region_provenance = _resolve_sp_ml_region(
            source_path=Path(prepared.source_path),
            out_dir_path=out_dir_path,
            calc_cfg=calc_cfg,
            model_indices_str=model_indices_str,
            model_indices_one_based=model_indices_one_based,
            protected_inputs=protected_inputs,
        )

        # Rename CLI-style keys to mlmm constructor kwargs to avoid duplicate-
        # value TypeError at super().__init__(charge=model_charge, ...).
        calc_cfg["model_charge"] = calc_cfg.pop("charge")
        calc_cfg["model_mult"] = calc_cfg.pop("spin")
        calc = mlmm(**calc_cfg)
        geom.set_calculator(calc)
        # pysisyphus' set_calculator REPLACES calc.freeze_atoms with the geometry's list, which
        # drops the FrozenMM layer the calculator itself unioned in. Restore the union so `sp`
        # force-masks the frozen layer like opt/freq do — otherwise |force|_max spans atoms the
        # user marked frozen and is not comparable with the value opt converged on.
        _frozen_layer = getattr(getattr(calc, "core", None), "frozen_layer_indices", None)
        if _frozen_layer:
            calc.freeze_atoms = sorted(set(calc.freeze_atoms) | set(_frozen_layer))

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
            mode = sp_cfg.get("hessian_calc_mode") or calc_cfg.get(
                "hessian_calc_mode", "FiniteDifference"
            )
            click.echo(f"[sp] computing active-coordinate ONIOM Hessian (mode={mode}) ...")
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

        elapsed_total = format_elapsed("[time] Elapsed Time for SP", time_start)
        summary = {
            "stage": "sp",
            "status": "ok",
            "input": str(prepared.source_path),
            "real_parm7": str(real_parm7),
            **calculator_provenance(calc_cfg),
            # charge/spin were popped + renamed to model_charge/model_mult
            # before mlmm() construction (see "Rename CLI-style keys" block).
            "charge": calc_cfg.get("model_charge"),
            "spin": calc_cfg.get("model_mult"),
            **ml_region_provenance,
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

        click.echo(format_elapsed("[time] Elapsed Time for SP", time_start))

    except _SPOutputCollisionError:
        raise
    except Exception as exc:
        render_cli_exception(
            exc, label="single-point",
            out_dir=out_dir_path or Path(sp_cfg.get("out_dir", OUT_DIR_SP)).resolve(),
            command="sp", time_start=time_start,
        )
    finally:
        prepared.cleanup()
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
