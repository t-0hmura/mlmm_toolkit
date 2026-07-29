"""End-to-end ML/MM extraction, path, TS, IRC, frequency, and DFT workflow."""

from __future__ import annotations

from collections.abc import Mapping
from copy import deepcopy
from dataclasses import dataclass
from pathlib import Path
from types import MappingProxyType
from typing import List, Sequence, Optional, Tuple, Dict, Any
import shutil
import tempfile

import gc
import logging
import sys
import math
import click
from mlmm.core.output import emit
from mlmm.cli.common_options import add_coord_type_option, add_precision_option, add_workers_options, add_backend_model_option, add_calc_file_option, add_deterministic_option, add_allow_charge_mult_mismatch_option
from mlmm.cli.decorators import canonicalize_calculator_section, make_is_param_explicit
# presentation dependency (workflow -> cli). One advanced-help callback +
# one visibility loop, shared with the lazily-loaded subcommands.
from mlmm.cli.help_pages import _show_advanced_subcommand_help, _hide_advanced_options
import time
import json
import yaml
import numpy as np
import torch

logger = logging.getLogger(__name__)


def _validate_postprocessing_dependencies(
    *, do_tsopt: bool, do_thermo: bool, do_dft: bool
) -> None:
    if (do_thermo or do_dft) and not do_tsopt:
        raise click.UsageError(
            "`all --thermo` and `all --dft` require `--tsopt`; an unoptimized "
            "MEP highest-energy image is not a validated transition state."
        )

# Biopython for PDB parsing (post-processing helpers)
from Bio import PDB

# pysisyphus helpers/constants
from pysisyphus.helpers import geom_loader
from pysisyphus.constants import BOHR2ANG, AU2KCALPERMOL

# Local imports from the package
from mlmm.workflows.extract import extract_api, compute_charge_summary, log_charge_summary
from mlmm.workflows.charge_prep import (
    configured_model_charge_spin,
    infer_present_terminal_cap_ids,
)
from mlmm.workflows import path_search as _path_search
from mlmm.workflows import path_opt as _path_opt
from mlmm.workflows import opt as _opt_cli
from mlmm.workflows import tsopt as _ts_opt
from mlmm.workflows import freq as _freq_cli
from mlmm.workflows import irc as _irc_cli

from mlmm.io.trj2fig import run_trj2fig
from mlmm.io.summary import (
    emit_method_citations,
    method_references,
    write_summary_log,
)
from mlmm.io.structure_formats import (
    coordinate_template_for,
    register_output_template_and_write_cif,
)
from mlmm.workflows.align_freeze import (
    align_and_refine_sequence_inplace,
    alignment_failed_pair_indices,
)
from mlmm.core.defaults import (
    GEOM_KW_DEFAULT,
    MLMM_CALC_KW,
    OUT_DIR_ALL,
    SEGMENTS_DIRNAME,
    THRESH_CHOICES,
    WORK_DIRNAME,
    fresh_dmf_config,
)
from mlmm.core.utils import (
    apply_ref_pdb_override,
    build_energy_diagram,
    close_matplotlib_figures,
    convert_xyz_to_pdb,
    ensure_dir,
    format_elapsed,
    prepare_input_structure,
    PreparedInputStructure,
    load_yaml_dict,
    load_pdb_atom_metadata,
    parse_scan_list_triples,
    read_bfactors_from_pdb,
    read_xyz_as_blocks,
    read_xyz_first_last,
    verbose_level,
    xyz_blocks_first_last,
)
from mlmm.core.result_commit import (
    commit_exact_bytes,
    commit_json_exact,
    with_current_run_id,
)
from mlmm.workflows._run_session import (
    CalculatorLease,
    InvocationManifest,
    RunSession,
    claim_public_output as _claim_public_output,
    current_key_output_files as _current_key_output_files,
    current_output_paths as _current_output_paths,
    declare_public_output as _declare_public_output,
    public_output_key as _public_output_key,
    refresh_current_public_outputs as _refresh_current_public_outputs,
)
from mlmm.cli.decorators import resolve_yaml_sources, load_merged_yaml_cfg
from mlmm.cli.preflight import validate_existing_files
from mlmm.workflows import scan as _scan_cli
from mlmm.domain.add_elem_info import assign_elements as _assign_elem_info
from mlmm.workflows.define_layer import define_layers as _define_layers
from mlmm.backends.mlmm_calc import mlmm as _mlmm_calc
from mlmm.workflows.mm_parm import (
    Args as _AutoMMArgs,
    ambertools_command_paths as _ambertools_command_paths,
    missing_ambertools_commands as _missing_ambertools_commands,
    parse_ligand_charge as _mm_parse_ligand_charge,
    parse_ligand_mult as _mm_parse_ligand_mult,
    run_pipeline as _mm_run,
)

AtomKey = Tuple[str, str, str, str, str, str]


_CALC_CONFIG_ONLY_KEYS = frozenset(
    {
        "backend_model",
        "charge",
        "model_indices_one_based",
        "model_indices_str",
        "precision",
        "spin",
    }
)


def _freeze_calc_value(value: Any) -> Any:
    """Recursively freeze a resolved calculator value for request-local reuse."""

    if isinstance(value, Mapping):
        return MappingProxyType(
            {str(key): _freeze_calc_value(item) for key, item in value.items()}
        )
    if isinstance(value, (list, tuple)):
        return tuple(_freeze_calc_value(item) for item in value)
    if isinstance(value, set):
        return frozenset(_freeze_calc_value(item) for item in value)
    return deepcopy(value)


def _thaw_calc_value(value: Any) -> Any:
    """Return an independent mutable value from :func:`_freeze_calc_value`."""

    if isinstance(value, Mapping):
        return {str(key): _thaw_calc_value(item) for key, item in value.items()}
    if isinstance(value, tuple):
        return [_thaw_calc_value(item) for item in value]
    if isinstance(value, frozenset):
        return {_thaw_calc_value(item) for item in value}
    return deepcopy(value)


@dataclass(frozen=True)
class _ResolvedCalculatorTemplate:
    """Immutable calculator configuration shared by one ``all`` invocation."""

    values: Mapping[str, Any]

    @classmethod
    def from_mapping(cls, values: Mapping[str, Any]) -> "_ResolvedCalculatorTemplate":
        cleaned = {
            str(key): value
            for key, value in values.items()
            if key not in _CALC_CONFIG_ONLY_KEYS
        }
        return cls(_freeze_calc_value(cleaned))

    def materialize(self) -> Dict[str, Any]:
        return _thaw_calc_value(self.values)


def _resolve_calculator_template(
    args_yaml: Optional[Path],
    *,
    backend: Optional[str],
    embedcharge: bool,
    embedcharge_explicit: bool,
    embedcharge_cutoff: Optional[float],
    link_atom_method: Optional[str],
    mm_backend: Optional[str],
    use_cmap: Optional[bool],
) -> _ResolvedCalculatorTemplate:
    """Materialize the effective calculator mapping exactly once.

    ``args_yaml`` is already the canonical C1 output and includes translated
    precision, worker, model, and custom-calculator CLI values.  The remaining
    ``all``-level calculator flags are overlaid once with Click's explicitness
    semantics; child/post-stage evaluators only derive state identity from the
    returned immutable template.
    """

    resolved: Dict[str, Any] = deepcopy(MLMM_CALC_KW)
    if args_yaml is not None:
        effective = canonicalize_calculator_section(load_yaml_dict(args_yaml))
        calc_section = effective.get("calc")
        if calc_section is not None and not isinstance(calc_section, Mapping):
            raise click.ClickException("The effective 'calc' section must be a mapping.")
        if isinstance(calc_section, Mapping):
            resolved.update(deepcopy(dict(calc_section)))

    # ``--calc-file`` is the documented final calculator selector and switches
    # the effective YAML to ``backend: custom``.  Do not let an accompanying
    # ``--backend`` overlay split the run so child processes use the custom
    # calculator while in-process pre-alignment/endpoint stages use an MLIP.
    if backend is not None and not resolved.get("calc_file"):
        resolved["backend"] = str(backend).lower()
    if embedcharge_explicit:
        resolved["embedcharge"] = bool(embedcharge)
    if embedcharge_cutoff is not None:
        resolved["embedcharge_cutoff"] = float(embedcharge_cutoff)
    if link_atom_method is not None:
        resolved["link_atom_method"] = str(link_atom_method).lower()
    if mm_backend is not None:
        resolved["mm_backend"] = str(mm_backend).lower()
    if use_cmap is not None:
        resolved["use_cmap"] = bool(use_cmap)

    return _ResolvedCalculatorTemplate.from_mapping(resolved)


def _stage_calc_kwargs(
    template: _ResolvedCalculatorTemplate,
    *,
    input_pdb: Path | str,
    real_parm7: Path | str,
    model_pdb: Path | str,
    charge: int,
    spin: int,
    use_bfactor_layers: bool,
) -> Dict[str, Any]:
    """Derive an evaluator config by changing state identity only."""

    kwargs = template.materialize()
    kwargs.update(
        {
            "input_pdb": str(input_pdb),
            "real_parm7": str(real_parm7),
            "model_pdb": str(model_pdb),
            "model_charge": int(charge),
            "model_mult": int(spin),
            "use_bfactor_layers": bool(use_bfactor_layers),
        }
    )
    return kwargs


def _resolve_mlip_provenance(
    *,
    backend: Optional[str],
    backend_model: Optional[str],
    calc_file: Optional[str],
    calc_factory: Optional[str],
    precision: Optional[str] = None,
    merged_yaml_cfg: Optional[Dict[str, Any]] = None,
) -> Tuple[str, str, Optional[str]]:
    """Return effective backend, model, and canonical precision metadata."""
    merged = merged_yaml_cfg or {}
    calc_yaml = merged.get("calc") or merged.get("mlmm") or {}
    if not isinstance(calc_yaml, dict):
        calc_yaml = {}
    effective_calc_file = calc_file or calc_yaml.get("calc_file")
    if effective_calc_file:
        effective_factory = (
            calc_factory or calc_yaml.get("calc_factory") or "get_calculator"
        )
        return (
            "custom",
            f"{Path(str(effective_calc_file)).name}:{effective_factory}",
            None,
        )

    backend_name = str(
        backend or calc_yaml.get("backend") or MLMM_CALC_KW["backend"]
    ).lower()
    model_keys = {
        "uma": "uma_model",
        "orb": "orb_model",
        "mace": "mace_model",
        "aimnet2": "aimnet2_model",
    }
    model_key = model_keys.get(backend_name)
    model = backend_model or calc_yaml.get("backend_model")
    if model is None and model_key is not None:
        model = calc_yaml.get(model_key, MLMM_CALC_KW.get(model_key))
    from mlmm.core.utils import calculator_provenance

    provenance_cfg: Dict[str, Any] = {"backend": backend_name}
    if model_key is not None and model is not None:
        provenance_cfg[model_key] = model
    precision_keys = {
        "uma": "uma_precision",
        "orb": "orb_precision",
        "mace": "mace_dtype",
    }
    precision_key = precision_keys.get(backend_name)
    if precision_key is not None:
        raw_precision = precision
        if raw_precision is None:
            raw_precision = calc_yaml.get(precision_key, calc_yaml.get("precision"))
        if raw_precision is not None and str(raw_precision).lower() != "auto":
            provenance_cfg[precision_key] = raw_precision
    resolved = calculator_provenance(provenance_cfg)
    return backend_name, str(model or "-"), resolved["mlip_precision"]


class _EchoState:
    """Encapsulate CLI output state for section-spacing logic."""

    def __init__(self) -> None:
        self._started = False

    def reset(self) -> None:
        self._started = False

    def echo(self, *args, **kwargs) -> None:
        kwargs.setdefault("narrative", False)
        emit(*args, **kwargs)
        self._started = True

    def section(self, message: str, **kwargs) -> None:
        # Section banners form the narrative backbone of the pipeline log, so
        # they default to narrative (visible at default verbosity). The leading
        # blank carries the same flag to preserve spacing around a shown banner.
        narrative = kwargs.setdefault("narrative", True)
        if self._started:
            emit(narrative=narrative)
        emit(message, **kwargs)
        self._started = True


_echo_state = _EchoState()


def _echo(*args, **kwargs) -> None:
    """Echo a line with local output tracking for section spacing.

    Untagged by default (visible at ``-v 3`` inside ``all``). Use
    ``_echo_detail`` for default ``-v 2`` stage details and ``narrative=True``
    for milestone lines.
    """
    _echo_state.echo(*args, **kwargs)


def _echo_detail(*args, **kwargs) -> None:
    """Echo a level-2 detail line with local output tracking."""
    kwargs.setdefault("detail", True)
    _echo_state.echo(*args, **kwargs)


def _echo_section(message: str, **kwargs) -> None:
    """Echo a section header (narrative) with a leading blank line unless first."""
    _echo_state.section(message, **kwargs)


def _emit_final_summary(
    out_dir: Path | None,
    time_start: float,
    manifest: Optional[InvocationManifest] = None,
    citation_payload: Optional[Dict[str, Any]] = None,
) -> None:
    """Print a visual `====== Pipeline summary ======` block + Elapsed line.

    Reads ``summary.json`` if present and lifts the most-asked-for numbers
    (status, highest local barrier, reaction energy, reactive-segment count,
    output dir) so the user sees them at the bottom of the log without
    scrolling back through `[diagram] Wrote ...` / `[time] Elapsed Time
    for X:` clutter. Falls back to just the Elapsed line when summary.json
    is absent (dry-run, early failure, TSOPT-only without aggregation).

    A stale ``summary.json`` from an earlier invocation reusing the same
    out_dir is never surfaced: when a ``manifest`` is supplied the summary is
    read only if the current run declared+claimed it as a public output.
    """
    summary: Dict[str, Any] = {}
    if out_dir is not None:
        summary_path = (Path(out_dir) / "summary.json").resolve(strict=False)
        current_summary = manifest is None or any(
            path == summary_path for path in manifest.paths("output.public.")
        )
        if current_summary and summary_path.exists():
            try:
                _loaded = json.loads(summary_path.read_text(encoding="utf-8"))
                if isinstance(_loaded, dict):
                    summary = _loaded
            except (OSError, json.JSONDecodeError):
                summary = {}
    if summary:
        _echo_section("====== Pipeline summary ======")
        status = summary.get("status")
        if status:
            _echo(f"Status: {status}", narrative=True)
        rls = summary.get("rate_limiting_step")
        if isinstance(rls, dict):
            barrier = rls.get("barrier_kcal")
            seg_idx = rls.get("segment")
            method = rls.get("method", "?")
            if barrier is not None:
                _echo(
                    f"Highest local barrier: {float(barrier):.2f} kcal/mol (segment {seg_idx}, method {method})",
                    narrative=True,
                )
        rxn_e = summary.get("overall_reaction_energy_kcal")
        if rxn_e is not None:
            _echo(f"Reaction energy: {float(rxn_e):.2f} kcal/mol", narrative=True)
        n_reactive = summary.get("n_segments_reactive")
        if n_reactive is not None:
            _echo(f"Reactive segments: {n_reactive}", narrative=True)
        # Report the pipeline out-dir the user passed (-o), not the stage
        # sub-dir that summary.json happens to record (e.g. <out>/path_search).
        out_dir_show = str(out_dir) if out_dir is not None else summary.get("out_dir")
        if out_dir_show:
            _echo(f"Output dir: {out_dir_show}", narrative=True)
        _echo(narrative=True)
    if citation_payload:
        emit_method_citations(citation_payload)
    _echo(format_elapsed("[all] Elapsed for Whole Pipeline", time_start), narrative=True)


def _run_cli_main(
    cmd_name: str,
    cli_obj,
    args: Sequence[str],
    *,
    on_nonzero: str = "warn",
    on_exception: str = "raise",
    prefix: Optional[str] = None,
) -> Optional[int]:
    """Run a Click command with temporary argv and consistent error handling.

    Returns the child's exit code (``0`` on success). A caller that must gate a
    downstream artifact on the child's success (for example, FREQ thermochemistry
    parsing) reads this instead of inferring success from a written file.
    """
    saved = list(sys.argv)
    label = prefix or cmd_name
    # In-proc subcommand dispatch — flag the child's banner / device echo to
    # stay silent so a 4-stage `all` pipeline doesn't reprint the same lines
    # `mlmm-toolkit ver. X` / `[calc] Resolved device: cuda` once per stage.
    from mlmm.core.utils import set_child_mode
    set_child_mode(True)
    rc: Optional[int] = 0
    try:
        sys.argv = ["mlmm", cmd_name] + list(args)
        _echo("")
        cli_obj.main(args=list(args), standalone_mode=False)
    except SystemExit as e:
        code = getattr(e, "code", 1)
        if code not in (None, 0):
            rc = code
            if code == 130:
                raise
            if on_nonzero == "raise":
                raise click.ClickException(f"[{label}] {cmd_name} exit code {code}.")
            _echo(f"[{label}] WARNING: {cmd_name} exited with code {code}")
    except Exception as e:
        rc = 1
        if on_exception == "raise":
            raise click.ClickException(f"[{label}] {cmd_name} failed: {e}")
        _echo(f"[{label}] WARNING: {cmd_name} failed: {e}")
    finally:
        sys.argv = saved
        set_child_mode(False)
        # Release GPU memory between pipeline stages to prevent OOM.
        # Subcommand finally blocks unbind their heavy locals (= None).
        # gc.collect is needed to break cyclic refs inside torch.nn.Module,
        # then empty_cache reclaims the CUDA allocator cache.
        gc.collect()
        if torch.cuda.is_available():
            torch.cuda.empty_cache()
        _echo("")
    return rc



def _append_cli_arg(args: List[str], flag: str, value: Any | None) -> None:
    """Append ``flag`` and ``value`` (converted to string) to ``args`` when ``value`` is not ``None``."""
    if value is None:
        return
    args.extend([flag, str(value)])


def _append_toggle_arg(args: List[str], flag: str, value: Any | None) -> None:
    """Append Click bool-toggle option as ``--flag`` / ``--no-flag`` when value is not ``None``."""
    if value is None:
        return
    if not isinstance(value, bool):
        raise TypeError(f"Toggle flag '{flag}' requires bool value, got {type(value).__name__}.")
    base = flag if not flag.startswith("--no-") else f"--{flag[5:]}"
    neg = f"--no-{base[2:]}"
    args.append(base if value else neg)


def _resolve_override_dir(default: Path, override: Path | None) -> Path:
    """Return ``override`` when provided (respecting absolute paths); otherwise ``default``."""
    if override is None:
        return default
    if override.is_absolute():
        return override
    return default.parent / override


def _build_effective_args_yaml(
    config_yaml: Optional[Path],
    override_yaml: Optional[Path],
    *,
    tmp_prefix: str,
) -> Tuple[Optional[Path], Dict[str, Any]]:
    """
    Build an effective args-yaml file path.

    Precedence for file layering:
      config_yaml < override_yaml
    """
    merged_raw, base_cfg, override_cfg = load_merged_yaml_cfg(config_yaml, override_yaml)
    merged = canonicalize_calculator_section(merged_raw)

    if config_yaml is None and override_yaml is None:
        return None, {}
    if config_yaml is None and merged == override_cfg:
        return override_yaml, merged
    if override_yaml is None and merged == base_cfg:
        return config_yaml, merged

    with tempfile.NamedTemporaryFile(
        mode="w",
        encoding="utf-8",
        suffix=".yaml",
        prefix=tmp_prefix,
        delete=False,
    ) as tf:
        yaml.safe_dump(merged, tf, sort_keys=False, allow_unicode=True)
        effective = Path(tf.name).resolve()

    # Register cleanup so the temp file is removed when the process exits.
    import atexit
    atexit.register(lambda p=effective: p.unlink(missing_ok=True))

    return effective, merged


def _inject_coord_type_into_args_yaml(
    args_yaml: Optional[Path],
    coord_type: Optional[str],
    tr_projection: Optional[str] = None,
    backend: Optional[str] = None,
    precision: Optional[str] = None,
    workers: Optional[int] = None,
    workers_per_node: Optional[int] = None,
    backend_model: Optional[str] = None,
    calc_file: Optional[str] = None,
    calc_factory: Optional[str] = None,
) -> Optional[Path]:
    """Inject geometry and backend-native calculator overrides into args YAML.

    Used by ``mlmm all --coord-type cart|dlc`` and the backend/model/precision
    options to propagate the choice through the all-pipeline args YAML. Only the opt/tsopt
    stages honour ``coord_type`` (DLC is meaningful there via microiteration);
    freq/scan/path stages are fixed to cartesian and ignore it. Returns the
    original ``args_yaml`` unchanged when there are no injected values.
    """
    cfg = {} if args_yaml is None else load_yaml_dict(args_yaml)
    if not isinstance(cfg, dict):
        cfg = {}
    cfg = canonicalize_calculator_section(cfg)
    existing_calc = cfg.get("calc")
    has_generic_calc_alias = isinstance(existing_calc, dict) and any(
        key in existing_calc
        for key in ("precision", "backend_model", "calc_file", "calc_factory")
    )
    if (
        coord_type is None
        and tr_projection is None
        and backend is None
        and precision is None
        and workers is None
        and workers_per_node is None
        and backend_model is None
        and calc_file is None
        and not has_generic_calc_alias
    ):
        return args_yaml
    if coord_type is not None or tr_projection is not None:
        geom_cfg = cfg.get("geom")
        if not isinstance(geom_cfg, dict):
            geom_cfg = {}
        geom_cfg = dict(geom_cfg)
        if coord_type is not None:
            geom_cfg["coord_type"] = coord_type
        if tr_projection is not None:
            geom_cfg["tr_projection"] = tr_projection
        cfg["geom"] = geom_cfg
    if (
        backend is not None
        or precision is not None
        or workers is not None
        or workers_per_node is not None
        or backend_model is not None
        or calc_file is not None
        or has_generic_calc_alias
    ):
        calc_cfg = cfg.get("calc")
        if not isinstance(calc_cfg, dict):
            calc_cfg = {}
        calc_cfg = dict(calc_cfg)
        # Set the explicit CLI backend before translating generic model and
        # precision tokens.  The translators dispatch from calc.backend; without
        # this, a config-less ``all --backend orb|mace|aimnet2`` is mistaken for
        # UMA and the child silently runs its backend defaults.
        if backend is not None:
            calc_cfg["backend"] = str(backend).strip().lower()
        # Translate --precision into the per-backend NATIVE kwarg
        # (uma_precision / orb_precision / mace_dtype) HERE and write THAT into
        # the args YAML. Writing the raw ``precision`` token instead leaks an
        # unknown kwarg into a sub-stage's Calculator(**calc_cfg) whenever that
        # stage's own --precision is unset (the stage does not re-translate the
        # YAML token), raising "Calculator.__init__ got an unexpected keyword
        # argument 'precision'". apply_precision_to_calc_cfg also pops any stray
        # raw ``precision`` key, keeping calc_cfg Calculator-clean.
        from mlmm.backends import apply_precision_to_calc_cfg, apply_backend_model_to_calc_cfg, apply_calc_file_to_calc_cfg, apply_workers_to_calc_cfg
        apply_workers_to_calc_cfg(calc_cfg, workers, workers_per_node)
        apply_backend_model_to_calc_cfg(calc_cfg, backend_model)
        apply_calc_file_to_calc_cfg(calc_cfg, calc_file, calc_factory)
        apply_precision_to_calc_cfg(calc_cfg, precision)
        cfg["calc"] = calc_cfg

    with tempfile.NamedTemporaryFile(
        mode="w",
        encoding="utf-8",
        suffix=".yaml",
        prefix="mlmm_all_coord_type_",
        delete=False,
    ) as tf:
        yaml.safe_dump(cfg, tf, sort_keys=False, allow_unicode=True)
        new_path = Path(tf.name).resolve()
    import atexit
    atexit.register(lambda p=new_path: p.unlink(missing_ok=True))
    return new_path


def _write_ml_region_definition(pocket_pdb: Path, dest: Path) -> Path:
    """
    Write a link-free atom-selection PDB for downstream ML/MM commands.

    Extractor-only ``HL/LKH`` atoms are removed because the calculator generates
    link H from the parm7 bonds crossing this selection.
    """
    dest.parent.mkdir(parents=True, exist_ok=True)
    try:
        lines = pocket_pdb.read_text(encoding="utf-8", errors="replace").splitlines(
            keepends=True
        )
    except FileNotFoundError:
        raise click.ClickException(f"[all] Pocket PDB not found while building ML region: {pocket_pdb}")
    kept: List[str] = []
    for line in lines:
        if line.startswith(("ATOM  ", "HETATM")):
            if line[12:16].strip() == "HL" and line[17:20].strip() == "LKH":
                continue
        kept.append(line)
    dest.write_text("".join(kept), encoding="utf-8")
    return dest.resolve()


def _write_bfactor_ml_subset(src_pdb: Path, dest: Path) -> Optional[Path]:
    """Write only the ML-layer (B-factor ≈ 0) ATOM/HETATM records of *src_pdb* to *dest*.

    Used as the ``--model-pdb`` ML-region definition when extraction is skipped but
    the input carries B-factor layers (detect-layer). Without it ml_region.pdb is the
    FULL input (skip_extract copies the whole system), so any downstream stage whose
    detect-layer can't read B-factors (e.g. an XYZ geometry) falls back to --model-pdb
    and treats the ENTIRE system as the ML/QM region (sum_Z huge → electron-count
    error). Returns ``None`` if the input has no B≈0 atoms (caller falls back).
    """
    try:
        ml_lines: List[str] = []
        for ln in open(src_pdb, "r", encoding="utf-8", errors="ignore"):
            if ln.startswith(("ATOM", "HETATM")):
                try:
                    bf = float(ln[60:66])
                except ValueError:
                    continue
                if abs(bf) < 0.5:
                    ml_lines.append(ln)
        if not ml_lines:
            _echo(
                f"[all] NOTE: {src_pdb} carries no B-factor ML layer (no B≈0 atoms); "
                f"the ML region falls back to the full input.",
                err=True,
            )
            return None
        dest.parent.mkdir(parents=True, exist_ok=True)
        with open(dest, "w", encoding="utf-8") as fh:
            fh.writelines(ml_lines)
            fh.write("END\n")
        return dest.resolve()
    except Exception as exc:
        # Never let a read/write failure look like "no ML atoms": the caller falls back to the
        # FULL input as the ML region, which is exactly the electron-count defect this helper
        # exists to prevent. Say so instead of returning a silent None.
        _echo(
            f"[all] WARNING: could not build the B≈0 ML-region subset from {src_pdb} ({exc}); "
            f"falling back to the FULL input as the ML region — downstream ML-region checks "
            f"will count the entire system.",
            err=True,
        )
        return None


def _summarize_existing_bfactor_layers(pdb_path: Path) -> Dict[str, int]:
    """Count atoms per B-factor layer (ML=0 / MovableMM=10 / FrozenMM=20).

    Atoms whose B-factor is none of these landmark values are reported under
    ``"other"`` so users can spot non-layered inputs quickly.
    """
    counts = {"ml": 0, "movable": 0, "frozen": 0, "other": 0}
    try:
        with open(pdb_path, "r") as fh:
            for ln in fh:
                if not (ln.startswith("ATOM") or ln.startswith("HETATM")):
                    continue
                try:
                    bf = float(ln[60:66])
                except ValueError:
                    counts["other"] += 1
                    continue
                if abs(bf) < 0.5:
                    counts["ml"] += 1
                elif abs(bf - 10.0) < 0.5:
                    counts["movable"] += 1
                elif abs(bf - 20.0) < 0.5:
                    counts["frozen"] += 1
                else:
                    counts["other"] += 1
    except FileNotFoundError:
        pass
    return counts


def _ml_region_atom_summary(pdb_path: Path) -> Optional[Tuple[int, int]]:
    """Return ``(atom_count, sum_Z)`` for the ATOM/HETATM records of *pdb_path*.

    Mirrors ``validate_charge_spin`` (pysisyphus ``ATOMIC_NUMBERS``); element is taken
    from PDB columns 77-78 and falls back to ``guess_element`` when blank. Used only for
    a human-readable one-line diagnostic on the ML-region definition, so any failure
    returns ``None`` — a summary-compute error must never break the run.
    """
    try:
        from pysisyphus.elem_data import ATOMIC_NUMBERS
        from mlmm.domain.add_elem_info import guess_element

        n_atoms = 0
        sum_z = 0
        with open(pdb_path, "r", encoding="utf-8", errors="ignore") as fh:
            for ln in fh:
                if not ln.startswith(("ATOM", "HETATM")):
                    continue
                n_atoms += 1
                elem = ln[76:78].strip()
                if not elem:
                    atname = ln[12:16].strip()
                    resn = ln[17:20].strip()
                    elem = guess_element(atname, resn, ln.startswith("HETATM"))
                z = ATOMIC_NUMBERS.get(str(elem).lower()) if elem else None
                if z is not None:
                    sum_z += int(z)
        return n_atoms, sum_z
    except Exception:
        return None


def _ml_region_summary_suffix(pdb_path: Path) -> str:
    """Return ``" (atoms=N, sumZ=Z)"`` for *pdb_path*, or ``""`` if it cannot be computed."""
    summary = _ml_region_atom_summary(pdb_path)
    if summary is None:
        return ""
    return f" (atoms={summary[0]}, sumZ={summary[1]})"


def _element_fix_path(root: Path, source: Path, ordinal: int) -> Path:
    """Allocate a collision-free private element-repair path."""

    return (Path(root) / f"{int(ordinal):03d}_{Path(source).name}").resolve()


def _materialize_all_coordinate_inputs(
    prepared_inputs: Sequence[PreparedInputStructure],
    work_dir: Path,
) -> Tuple[Path, ...]:
    """Overlay XYZ coordinates on private PDB topology inputs."""

    coordinate_input_dir = work_dir / "coordinate_inputs"
    coordinate_destinations = {
        input_ordinal: (
            coordinate_input_dir / f"endpoint_{input_ordinal:02d}.pdb"
        ).resolve()
        for input_ordinal, prepared in enumerate(prepared_inputs, start=1)
        if prepared.geom_path.suffix.lower() == ".xyz"
    }
    protected_paths = {
        path.resolve()
        for prepared in prepared_inputs
        for path in (prepared.source_path, prepared.geom_path)
    }
    collisions = set(coordinate_destinations.values()) & protected_paths
    if collisions:
        collision = sorted(collisions, key=str)[0]
        raise click.BadParameter(
            "Input and --ref-pdb files must be outside the managed "
            f"all-workflow path {collision}."
        )

    coordinate_inputs: List[Path] = []
    for input_ordinal, prepared in enumerate(prepared_inputs, start=1):
        if prepared.geom_path.suffix.lower() != ".xyz":
            coordinate_inputs.append(prepared.source_path.resolve())
            continue
        ensure_dir(coordinate_input_dir)
        coordinate_pdb = coordinate_destinations[input_ordinal]
        convert_xyz_to_pdb(
            prepared.geom_path,
            prepared.source_path,
            coordinate_pdb,
        )
        coordinate_inputs.append(coordinate_pdb.resolve())
    return tuple(coordinate_inputs)


def _mm_charge_mapping(expr: Optional[str]) -> Dict[str, int]:
    """Return a ligand-charge mapping for mm_parm when ``expr`` uses RES=Q or RES:Q syntax."""
    if not expr:
        return {}
    if ("=" not in expr) and (":" not in expr):
        return {}
    try:
        return _mm_parse_ligand_charge(expr)
    except Exception as exc:  # pragma: no cover - defensive
        raise click.ClickException(f"[all] Invalid --ligand-charge mapping for mm_parm: {exc}")


def _mm_mult_mapping(expr: Optional[str]) -> Dict[str, int]:
    """Return a ligand-multiplicity mapping for mm_parm when ``expr`` uses RES=M or RES:M syntax."""
    if not expr:
        return {}
    try:
        return _mm_parse_ligand_mult(expr)
    except Exception as exc:  # pragma: no cover - defensive
        raise click.ClickException(f"[all] Invalid --auto-mm-ligand-mult mapping for mm_parm: {exc}")


def _build_mm_parm7(
    pdb: Path,
    ligand_charge_expr: Optional[str],
    ligand_mult_expr: Optional[str],
    out_dir: Path,
    ff_set: str,
    add_ter: bool,
    keep_temp: bool,
    auto_disulfide: bool = True,
) -> Tuple[Path, Path]:
    """Run mm_parm on ``pdb`` and return (parm7, rst7)."""
    out_dir.mkdir(parents=True, exist_ok=True)
    out_prefix = (out_dir / pdb.stem).resolve()
    out_prefix.parent.mkdir(parents=True, exist_ok=True)
    args = _AutoMMArgs(
        pdb=pdb.resolve(),
        out_prefix=str(out_prefix),
        ligand_charge=_mm_charge_mapping(ligand_charge_expr),
        ligand_mult=_mm_mult_mapping(ligand_mult_expr),
        keep_temp=bool(keep_temp),
        add_ter=bool(add_ter),
        auto_disulfide=bool(auto_disulfide),
        add_h=False,
        ph=7.0,
        ff_set=str(ff_set),
        out_prefix_given=True,
    )
    try:
        _mm_run(args)
    except SystemExit as exc:  # pragma: no cover - click exit translation
        code = getattr(exc, "code", 1)
        raise click.ClickException(f"[all] mm_parm exited with code {code}.")
    except Exception as exc:
        raise click.ClickException(f"[all] mm_parm failed: {exc}")

    parm7 = Path(f"{args.out_prefix}.parm7").resolve()
    rst7 = Path(f"{args.out_prefix}.rst7").resolve()
    if not parm7.exists():
        raise click.ClickException(f"[all] mm_parm did not produce parm7 at {parm7}")
    if not rst7.exists():
        raise click.ClickException(f"[all] mm_parm did not produce rst7 at {rst7}")
    return parm7, rst7


def _ts_imag_record(n_imag, imag_freqs_cm=None) -> dict:
    """Build the `ts_imag` record for summary.json / summary.log.

    The formatter has always been able to render the frequency and to warn on a
    soft mode, but nothing populated it, so `nu_imag (max)` printed `-` and the
    warning never fired: a -15 cm^-1 soft mode and a -450 cm^-1 reaction
    coordinate were indistinguishable downstream. Carry the values through.
    """
    record: dict = {"n_imag": int(n_imag)}
    freqs = [float(x) for x in (imag_freqs_cm or [])]
    if freqs:
        record["imag_freqs_cm"] = freqs
        record["nu_imag_max_cm"] = min(freqs)          # most negative = the certifying mode
        record["min_abs_imag_cm"] = min(abs(f) for f in freqs)
    return record


def _parse_atom_key_from_line(line: str) -> Optional[AtomKey]:
    """Extract a structural identity key from a PDB ATOM/HETATM record."""
    if not (line.startswith("ATOM") or line.startswith("HETATM")):
        return None
    atomname = line[12:16].strip()
    altloc = (line[16] if len(line) > 16 else " ").strip()
    resname = line[17:20].strip()
    chain = (line[21] if len(line) > 21 else " ").strip()
    resseq = line[22:26].strip()
    icode = (line[26] if len(line) > 26 else " ").strip()
    return (chain, resname, resseq, icode, atomname, altloc)


def _key_variants(key: AtomKey) -> List[AtomKey]:
    """Return key variants with progressively relaxed identity fields (deduplicated)."""
    chain, resn, resseq, icode, atom, alt = key
    raw_variants = [
        (chain, resn, resseq, icode, atom, alt),
        (chain, resn, resseq, icode, atom, ""),
        (chain, resn, resseq, "", atom, alt),
        (chain, resn, resseq, "", atom, ""),
    ]
    seen: set[AtomKey] = set()
    variants: List[AtomKey] = []
    for variant in raw_variants:
        if variant in seen:
            continue
        seen.add(variant)
        variants.append(variant)
    return variants


def _parse_scan_lists_literals(
    scan_lists_raw: Sequence[str],
    atom_meta: Optional[Sequence[Dict[str, Any]]] = None,
    one_based: bool = True,
) -> List[List[Tuple[int, int, float]]]:
    """Parse ``--scan-lists`` literals without re-basing atom indices.

    Parameters
    ----------
    one_based : bool, default True
        Honour the CLI ``--scan-one-based`` toggle so users can pass 0-based
        indices via ``all`` and have them forwarded unchanged to ``scan``.
    """
    stages: List[List[Tuple[int, int, float]]] = []
    for idx_stage, literal in enumerate(scan_lists_raw, start=1):
        tuples, _ = parse_scan_list_triples(
            literal,
            one_based=one_based,
            atom_meta=atom_meta,
            option_name=f"--scan-lists #{idx_stage}",
            return_one_based=one_based,
        )
        if not tuples:
            raise click.BadParameter(
                f"--scan-lists #{idx_stage} must contain at least one (i,j,target) triple."
            )
        stages.append(tuples)
    return stages


def _format_scan_stage(stage: List[Tuple[int, int, float]]) -> str:
    """Serialize a scan stage back into a Python-like literal string."""
    return "[" + ", ".join(f"({i},{j},{target})" for (i, j, target) in stage) + "]"


def _round_charge_with_note(q: float) -> int:
    """
    Cast the extractor's ML-region charge (float) to an integer suitable for the path search.
    If it is not already an integer within 1e-6, round to the nearest integer with a console note.
    """
    q_rounded = int(round(float(q)))
    if not math.isfinite(q):
        raise click.BadParameter(f"Computed total charge is non-finite: {q!r}")
    if abs(float(q) - q_rounded) > 1e-6:
        click.echo(f"[all] NOTE: extractor ML-region charge = {q:g} → rounded to integer {q_rounded} for the path search.")
    return q_rounded


def _derive_charge_from_ligand_charge_when_extract_skipped(
    pdb_path: Path,
    ligand_charge: Optional[str],
) -> Optional[int]:
    """Derive the ML-region charge from a PDB using extract-style charge summary.

    *pdb_path* may be a full-complex PDB or a --model-pdb pocket.
    """
    if ligand_charge is None:
        return None
    try:
        parser = PDB.PDBParser(QUIET=True)
        complex_struct = parser.get_structure("complex", str(pdb_path))
        selected_ids = {res.get_full_id() for res in complex_struct.get_residues()}
        keep_ncap_ids, keep_ccap_ids = infer_present_terminal_cap_ids(
            complex_struct,
            selected_ids,
        )
        summary = compute_charge_summary(
            complex_struct,
            selected_ids,
            set(),
            ligand_charge,
            keep_ncap_ids=keep_ncap_ids,
            keep_ccap_ids=keep_ccap_ids,
        )
        log_charge_summary("[all]", summary)
        q_total = float(summary.get("total_charge", 0.0))
        click.echo(f"[all] Charge summary from {pdb_path.name} (--ligand-charge without extraction):")
        click.echo(
            f"  Protein: {summary.get('protein_charge', 0.0):+g},  "
            f"Ligand: {summary.get('ligand_total_charge', 0.0):+g},  "
            f"Ions: {summary.get('ion_total_charge', 0.0):+g},  "
            f"Total: {q_total:+g}"
        )
        return _round_charge_with_note(q_total)
    except Exception as e:
        click.echo(
            f"[all] NOTE: failed to derive ML-region charge from --ligand-charge: {e}",
            err=True,
        )
        return None


def _derive_ml_charge_from_layered_pdb(
    pdb_path: Path,
    ligand_charge: Optional[str],
) -> Optional[int]:
    """Derive the ML-region (B≈0) charge from a B-factor-layered PDB when extraction
    is skipped (ts-only / no -c), reusing extract's ``compute_charge_summary`` WITH
    terminal-cap correction at the ML/MM boundary.

    Needed because ML ⊊ system in ONIOM: summing the whole input gives the total
    system charge (mis-applied as the ML model charge — the ts-only charge bug),
    while summing only the B≈0 atoms misses the backbone-cut terminal caps (off by
    the number of cut termini). The cut residues (peptide neighbor not in the ML
    set) are flagged as N-/C-caps, exactly as the extract path does. Validated:
    Returns ``None`` on any failure so the caller can
    fall back to the full-input derivation.
    """
    if ligand_charge is None:
        return None
    try:
        from mlmm.workflows.extract import compute_charge_summary

        parser = PDB.PDBParser(QUIET=True)
        st = parser.get_structure("complex", str(pdb_path))
        ml_ids = {
            r.get_full_id()
            for r in st.get_residues()
            if any(abs(a.get_bfactor()) < 0.5 for a in r.get_atoms())
        }
        if not ml_ids:
            return None
        keep_ncap, keep_ccap = infer_present_terminal_cap_ids(st, ml_ids)
        summary = compute_charge_summary(
            st, ml_ids, set(), ligand_charge,
            keep_ncap_ids=keep_ncap, keep_ccap_ids=keep_ccap,
        )
        q_total = float(summary.get("total_charge", 0.0))
        click.echo(
            f"[all] ML-region charge from {pdb_path.name} "
            f"(detect-layer, extraction skipped; cap-corrected): "
            f"Protein: {summary.get('protein_charge', 0.0):+g},  "
            f"Ligand: {summary.get('ligand_total_charge', 0.0):+g},  "
            f"Total: {q_total:+g}"
        )
        return _round_charge_with_note(q_total)
    except Exception as e:
        click.echo(
            f"[all] NOTE: cap-aware ML-region charge derivation failed: {e}",
            err=True,
        )
        return None


def _pdb_needs_elem_fix(p: Path) -> bool:
    """
    Return True if the PDB has at least one ATOM/HETATM record whose element field (cols 77–78) is empty.
    This is a light-weight check to decide whether to run add_elem_info.
    """
    try:
        with p.open("r", encoding="utf-8", errors="ignore") as fh:
            for line in fh:
                if line.startswith("ATOM") or line.startswith("HETATM"):
                    if len(line) < 78 or not line[76:78].strip():
                        return True
        return False
    except Exception:
        # On I/O errors, skip fixing (use original)
        return False


# ---------- Post-processing helpers (minimal, reuse internals) ----------


# Severity order for composing the legacy completeness axis with the
# convergence-gated aggregate: scientific_status is never LESS severe than the
# legacy ``status`` (so a nonconverged leaf can only demote, never promote).
_STATUS_SEVERITY = {"success": 0, "partial": 1, "failed": 2}


def _is_reactive_segment(item: Any) -> bool:
    """Return whether a segment legitimately requires TS post-processing."""
    if not isinstance(item, dict):
        return False
    kind = item.get("kind", "seg")
    if kind == "tsopt":
        return True
    if kind != "seg":
        return False
    # Legacy/directly constructed segment records predate bond-change
    # serialization and remain reactive.  Only an explicit no-change result
    # suppresses post-processing.
    if "bond_changes" not in item:
        return True
    changes = str(item.get("bond_changes", "")).strip()
    return bool(changes and changes != "(no covalent changes detected)")


def _read_irc_outcome(irc_dir: Path) -> Dict[str, Any]:
    """Read the IRC child's ``result.json`` into a fail-closed usability record.

    The IRC leaf is *usable* only when the child reports ``scientific_status ==
    "success"`` — i.e. every requested direction explicitly converged. A
    missing / unreadable result, or any nonconverged requested direction, yields
    ``usable=False`` while the endpoint trajectory remains a reportable artifact.
    """

    outcome: Dict[str, Any] = {
        "usable": False,
        "reason": "irc_result_missing",
        "scientific_status": None,
        "forward_converged": None,
        "backward_converged": None,
        "n_frames_forward": None,
        "n_frames_backward": None,
        "traj": None,
    }
    result_path = irc_dir / "result.json"
    if not result_path.exists():
        return outcome
    try:
        data = json.loads(result_path.read_text(encoding="utf-8"))
    except (OSError, ValueError):
        outcome["reason"] = "irc_result_unreadable"
        return outcome
    if not isinstance(data, dict):
        outcome["reason"] = "irc_result_unreadable"
        return outcome

    sci = data.get("scientific_status")
    outcome["scientific_status"] = sci
    outcome["forward_converged"] = data.get("forward_converged")
    outcome["backward_converged"] = data.get("backward_converged")
    outcome["n_frames_forward"] = data.get("n_frames_forward")
    outcome["n_frames_backward"] = data.get("n_frames_backward")
    _files = data.get("files") if isinstance(data.get("files"), dict) else {}
    outcome["traj"] = _files.get("finished_irc")

    if sci == "success":
        outcome["usable"] = True
        outcome["reason"] = "ok"
    elif isinstance(sci, str):
        outcome["usable"] = False
        reasons = data.get("scientific_status_reasons")
        outcome["reason"] = (
            ";".join(str(r) for r in reasons)
            if isinstance(reasons, list) and reasons
            else f"irc_{sci}"
        )
    else:
        # No explicit status field: fail closed rather than trust file existence.
        outcome["usable"] = False
        outcome["reason"] = "irc_status_unknown"
    return outcome


def _read_path_opt_segment_converged(seg_dir: Path) -> Optional[bool]:
    """Read a path-opt segment child's reported MEP convergence (tri-state).

    Reads the additive ``stage_outcomes`` leaf ``converged`` bit from the child's
    ``result.json`` — sourced from the optimizer for both GSM and DMF (IPOPT
    bit). Returns ``None`` when no readable signal exists (fail-closed: a
    missing or unreadable child result never promotes the segment to converged).
    """
    try:
        rj = seg_dir / "result.json"
        if not rj.exists():
            return None
        data = json.loads(rj.read_text(encoding="utf-8")) or {}
        for leaf in data.get("stage_outcomes") or []:
            if isinstance(leaf, dict) and isinstance(leaf.get("converged"), bool):
                return bool(leaf["converged"])
        return None
    except Exception as exc:
        logger.debug("Failed to read path-opt segment convergence %s: %s", seg_dir, exc)
        return None


def _read_opt_endpoint_converged(opt_dir: Path) -> Optional[bool]:
    """Read an endpoint-opt child's reported convergence (tri-state).

    The ``opt`` subcommand writes the final optimizer's ``is_converged`` bit to
    its ``result.json`` as ``status`` = ``"converged"`` / ``"not_converged"``
    (mirroring how :func:`_read_irc_outcome` reads the IRC child). Returns
    ``None`` when no readable signal exists (fail-closed: a missing / unreadable
    child result never promotes the endpoint to converged, so a nonconverged or
    unrun endpoint optimization cannot silently promote its segment to success).
    """
    try:
        rj = opt_dir / "result.json"
        if not rj.exists():
            return None
        data = json.loads(rj.read_text(encoding="utf-8")) or {}
        if not isinstance(data, dict):
            return None
        status = data.get("status")
        if status == "converged":
            return True
        if status == "not_converged":
            return False
        return None
    except Exception as exc:
        logger.debug("Failed to read endpoint-opt convergence %s: %s", opt_dir, exc)
        return None


def _pipeline_aggregate_truth(
    summary: dict,
    *,
    post_segments: Optional[list],
    config: Optional[dict],
    legacy_status: str,
    legacy_reasons: Optional[Sequence[str]] = None,
):
    """Compose the ``all``-pipeline aggregate from per-segment leaves.

    One required :class:`LeafOutcome` is built per reactive segment. A path
    segment is usable only when its MEP and every post-processing convergence
    signal are explicitly ``True``. A direct TSOPT segment has no MEP stage, so
    it is gated by its IRC and, when present, both endpoint optimizations. A
    dict-present / trajectory-present but nonconverged leaf never counts toward
    fail-closed completeness — a never_stop / max-cycle IRC therefore
    cannot yield ``scientific_status == "success"``.

    The convergence-gated aggregate is then composed with the legacy completeness
    axis (``legacy_status`` from :func:`_derive_pipeline_status`, which already
    covers DFT / thermo / n_imag): ``scientific_status`` is the MORE severe of the
    two so the new field carries at least as much information as the legacy
    ``status``. The
    legacy ``status`` string itself is untouched (byte-compatible).
    """

    from mlmm.workflows._outcomes import (
        AggregateTruth,
        aggregate_workflow_truth,
        make_leaf,
    )

    segments = summary.get("segments") or []
    reactive = [s for s in segments if _is_reactive_segment(s)]
    cfg = config or {}
    post_requested = any(
        bool(cfg.get(name)) for name in ("tsopt", "thermo", "dft")
    )
    tsopt_requested = bool(cfg.get("tsopt"))
    legacy_reasons = list(legacy_reasons or [])

    post_by_idx: Dict[Any, dict] = {}
    if post_segments is not None:
        for ps in post_segments:
            if isinstance(ps, dict) and ps.get("index") is not None:
                post_by_idx[ps.get("index")] = ps

    def _and3(a: Optional[bool], b: Optional[bool]) -> Optional[bool]:
        if a is False or b is False:
            return False
        if a is None or b is None:
            return None
        return True

    leaves: List[Any] = []
    expected: List[str] = []
    for s in reactive:
        idx = s.get("index")
        if idx is None:
            continue
        seg_id = f"segment_{idx}"
        expected.append(seg_id)
        post = post_by_idx.get(idx)
        reason = ""
        artifacts: List[str] = []
        # The segment's own reported convergence, threaded from path_search's
        # SegmentReport. A missing field is None (fail-closed), never a silent
        # True.
        _seg_conv = s.get("converged")
        seg_converged: Optional[bool] = _seg_conv if isinstance(_seg_conv, bool) else None
        # ``kind=tsopt`` is the direct-TS branch: no MEP child runs and therefore
        # no MEP convergence field exists. Only that explicit kind bypasses the
        # MEP gate; unknown/future kinds remain fail-closed like path segments.
        mep_converged: Optional[bool] = (
            True if s.get("kind") == "tsopt" else seg_converged
        )
        if post is not None:
            # Post-processing ran: compose its explicit IRC / endpoint records
            # with the MEP engine's own convergence.  Successful downstream
            # work must never promote a nonconverged/unknown path segment.
            converged: Optional[bool] = mep_converged
            if converged is not True:
                reason = (
                    "mep_not_converged"
                    if converged is False
                    else "mep_convergence_unknown"
                )
            irc = post.get("irc")
            if isinstance(irc, dict):
                _u = irc.get("usable")
                _irc_conv = True if _u is True else (False if _u is False else None)
                converged = _and3(converged, _irc_conv)
                if _irc_conv is not True and not reason:
                    reason = f"irc:{irc.get('reason') or 'not_usable'}"
                _traj = irc.get("traj")
                if _traj:
                    artifacts.append(str(_traj))
            elif tsopt_requested:
                # IRC requested but no directional record: fail closed
                # rather than trust the trajectory file's existence.
                converged = _and3(converged, None)
                if not reason:
                    reason = "irc_missing"
            if tsopt_requested and s.get("kind") != "tsopt":
                endpoint_assignment = post.get("endpoint_assignment")
                connectivity = (
                    endpoint_assignment.get("connectivity_validated")
                    if isinstance(endpoint_assignment, dict)
                    else None
                )
                connectivity_truth = (
                    connectivity if isinstance(connectivity, bool) else None
                )
                converged = _and3(converged, connectivity_truth)
                if connectivity_truth is not True and not reason:
                    reason = "irc_endpoint_connectivity_unvalidated"
            eo = post.get("endpoint_opt")
            if isinstance(eo, dict):
                for _k in ("reactant_converged", "product_converged"):
                    if _k in eo:
                        _v = eo.get(_k)
                        converged = _and3(converged, _v if isinstance(_v, bool) else None)
                        if not (isinstance(_v, bool) and _v) and not reason:
                            reason = f"endpoint_opt:{_k}"
        elif post_requested:
            # tsopt was requested but this segment's IRC/endpoint post-processing
            # has not run yet (the intermediate MEP summary is written before
            # post-processing). Fail closed rather than promote a reactive leaf on
            # the MEP trajectory's existence alone.
            converged = None
            reason = "post_missing"
        else:
            # Path-only final summary (no tsopt): the segment's own reported
            # convergence is the whole truth. A missing/unknown field fails closed
            # (None) — never default to True.
            converged = mep_converged
            if converged is not True and not reason:
                reason = "not_converged" if converged is False else "convergence_unknown"
        leaves.append(
            make_leaf(
                "all",
                seg_id,
                required=True,
                executed=(post is not None) if post_requested else True,
                converged=converged,
                reason=reason,
                artifacts=artifacts,
            )
        )

    if leaves:
        agg = aggregate_workflow_truth(leaves, expected)
        agg_sci = agg.scientific_status
        agg_exec = agg.execution_status
        agg_reasons = list(agg.status_reasons)
        observed = (
            [
                item_id
                for item_id in expected
                if item_id.removeprefix("segment_") in {
                    str(index) for index in post_by_idx
                }
            ]
            if post_requested
            else list(agg.observed_item_ids)
        )
    else:
        # No reactive-segment leaves to gate on (degenerate/endpoint-only
        # summary): mirror the legacy completeness axis rather than manufacture a
        # spurious failure.
        agg_sci = legacy_status
        agg_exec = "failed" if legacy_status == "failed" else "completed"
        agg_reasons = []
        observed = list(expected)

    # Compose with the legacy completeness axis: keep the MORE severe verdict.
    if _STATUS_SEVERITY.get(legacy_status, 0) >= _STATUS_SEVERITY.get(agg_sci, 0):
        scientific = legacy_status
    else:
        scientific = agg_sci
    execution = "failed" if (legacy_status == "failed" or agg_exec == "failed") else "completed"
    reasons = legacy_reasons + [r for r in agg_reasons if r not in legacy_reasons]

    return AggregateTruth(
        execution_status=execution,
        scientific_status=scientific,
        status_reasons=tuple(reasons),
        expected_item_ids=tuple(expected),
        observed_item_ids=tuple(observed),
    )


def _apply_pipeline_truth(
    summary: dict,
    *,
    post_segments: Optional[list],
    config: Optional[dict],
    legacy_status: str,
    legacy_reasons: Optional[Sequence[str]] = None,
) -> None:
    """Write the outcome axes onto ``summary`` in place.

    Never touches the legacy overloaded ``status`` field; only adds
    ``execution_status`` / ``scientific_status`` / expected+observed IDs and the
    distinct ``scientific_status_reasons`` key.
    """

    truth = _pipeline_aggregate_truth(
        summary,
        post_segments=post_segments,
        config=config,
        legacy_status=legacy_status,
        legacy_reasons=legacy_reasons,
    )
    summary["execution_status"] = truth.execution_status
    summary["scientific_status"] = truth.scientific_status
    summary["expected_item_ids"] = list(truth.expected_item_ids)
    summary["observed_item_ids"] = list(truth.observed_item_ids)
    if truth.scientific_status != "success" and truth.status_reasons:
        summary["scientific_status_reasons"] = list(truth.status_reasons)
    else:
        summary.pop("scientific_status_reasons", None)


def _derive_pipeline_status(
    summary: dict,
    *,
    post_segments: Optional[list],
    config: Optional[dict],
) -> Tuple[str, List[str]]:
    """Return aggregate pipeline status and machine-readable reasons."""
    segments = summary.get("segments") or []
    has_diagrams = bool(summary.get("energy_diagrams"))
    if not segments and not has_diagrams:
        return "failed", ["no usable path segments or energy diagrams were produced"]

    reasons: List[str] = []
    if not has_diagrams:
        reasons.append("no usable energy diagram was produced")

    cfg = config or {}
    requested = any(bool(cfg.get(name)) for name in ("tsopt", "thermo", "dft"))
    if post_segments is not None and requested:
        logs = [item for item in post_segments if isinstance(item, dict)]
        reactive_ids = {
            item.get("index")
            for item in segments
            if _is_reactive_segment(item) and item.get("index") is not None
        }
        if not logs and reactive_ids:
            reasons.append("requested post-processing produced no segment records")
        observed_ids = {
            item.get("index")
            for item in logs
            if item.get("index") is not None
        }
        for missing_idx in sorted(reactive_ids - observed_ids, key=str):
            reasons.append(
                f"segment {missing_idx}: requested post-processing record is missing"
            )
        for ordinal, item in enumerate(logs, start=1):
            prefix = f"segment {item.get('index', ordinal)}"
            if cfg.get("tsopt"):
                if not isinstance(item.get("mlip"), dict):
                    reasons.append(f"{prefix}: TSOPT/IRC refined MLIP energies are missing")
                if not item.get("irc_traj"):
                    reasons.append(f"{prefix}: IRC trajectory is missing")
                ts_imag = item.get("ts_imag")
                if not isinstance(ts_imag, dict) or ts_imag.get("n_imag") is None:
                    reasons.append(f"{prefix}: TS imaginary-mode validation is missing")
                elif int(ts_imag["n_imag"]) != 1:
                    reasons.append(
                        f"{prefix}: TS imaginary-mode validation found "
                        f"n_imag={int(ts_imag['n_imag'])}, expected 1"
                    )
            if cfg.get("thermo"):
                if not isinstance(item.get("gibbs_mlip"), dict):
                    reasons.append(f"{prefix}: MLIP thermochemistry result is missing")
            if cfg.get("dft"):
                if not isinstance(item.get("dft"), dict):
                    reasons.append(f"{prefix}: DFT result is missing")
                if cfg.get("thermo") and not isinstance(
                    item.get("gibbs_dft_mlip"), dict
                ):
                    reasons.append(f"{prefix}: DFT//MLIP/MM thermochemistry result is missing")

    if cfg.get("dft") and cfg.get("dft_status") == "failed":
        reasons.append("DFT failed for one or more TS-only states")
    reasons = list(dict.fromkeys(reasons))
    return ("partial" if reasons else "success"), reasons


def _enrich_summary(
    summary: dict,
    *,
    version: str,
    pipeline_mode: str,
    mlip_backend: str,
    mlip_precision: Optional[str] = None,
    charge: int,
    spin: int,
    command: str = "",
    post_segments: Optional[list] = None,
    config: Optional[dict] = None,
    freeze_atoms: Optional[str] = None,
    out_dir: Optional[Path] = None,
    mlip_model: Optional[str] = None,
    manifest: Optional[InvocationManifest] = None,
) -> dict:
    """Add machine-readable metadata to summary dict for AI agent consumption.

    The resulting dict is written as summary.json and is the machine-readable
    pipeline output consumed by MCP tools and other clients. Formatted tables
    and the filesystem tree remain specific to summary.log.
    """
    from mlmm import __version__
    from mlmm.core.utils import RESULT_JSON_SCHEMA_VERSION

    segments = summary.get("segments", [])
    reactive = [s for s in segments if _is_reactive_segment(s)]
    n_reactive = len(reactive)
    ts_only = pipeline_mode == "tsopt-only"

    status, status_reasons = _derive_pipeline_status(
        summary,
        post_segments=post_segments,
        config=config,
    )

    post_by_idx = {
        item.get("index"): item
        for item in (post_segments or [])
        if isinstance(item, dict)
    }
    best_method = None
    method_key = None
    rls = None
    if reactive:
        # ``rate_limiting_step`` is a legacy schema key for the highest
        # independently referenced local barrier, not a kinetic RLS assignment.
        def _has_finite_barrier(block: Any) -> bool:
            if not isinstance(block, dict) or block.get("barrier_kcal") is None:
                return False
            try:
                return bool(np.isfinite(float(block["barrier_kcal"])))
            except (TypeError, ValueError):
                return False

        for candidate_method, candidate_key in (
            ("DFT//MLIP/MM_Gibbs", "gibbs_dft_mlip"),
            ("DFT", "dft"),
            ("MLIP_Gibbs", "gibbs_mlip"),
            ("MLIP", "mlip"),
        ):
            if all(
                _has_finite_barrier(
                    (post_by_idx.get(seg.get("index")) or {}).get(candidate_key)
                )
                for seg in reactive
            ):
                best_method = candidate_method
                method_key = candidate_key
                break
        if best_method is None:
            best_method = "MLIP" if ts_only else "MEP"

        max_barrier = -1e9
        for s in reactive:
            refined = (
                (post_by_idx.get(s.get("index")) or {}).get(method_key)
                if method_key
                else None
            )
            if isinstance(refined, dict) and refined.get("barrier_kcal") is not None:
                b = float(refined["barrier_kcal"])
                method = best_method
            else:
                b = float(s.get("barrier_kcal", 0) or 0)
                method = (
                    "MLIP"
                    if ts_only and s.get("kind") == "tsopt"
                    else "MEP"
                )
            if b > max_barrier:
                max_barrier = b
                rls = {
                    "segment": s.get("index"),
                    "barrier_kcal": round(b, 2),
                    "method": method,
                }
                if not (ts_only and s.get("kind") == "tsopt"):
                    rls["mep_barrier_kcal"] = round(
                        float(s.get("barrier_kcal", 0) or 0), 2
                    )

    overall_rxn_e = None
    overall_rxn_method = None
    diagrams_by_name = {
        str(diag.get("name", "")): diag
        for diag in summary.get("energy_diagrams", [])
        if isinstance(diag, dict)
    }
    ranks = {"MEP": 0, "MLIP": 1, "MLIP_Gibbs": 2, "DFT": 3, "DFT//MLIP/MM_Gibbs": 4}
    max_rank = ranks.get(best_method or "MEP", 0)
    for diagram_name, method in (
        ("energy_diagram_G_DFT_plus_MLIP_all", "DFT//MLIP/MM_Gibbs"),
        ("energy_diagram_DFT_all", "DFT"),
        ("energy_diagram_G_MLIP_all", "MLIP_Gibbs"),
        ("energy_diagram_MLIP_all", "MLIP"),
        ("energy_diagram_MEP", "MEP"),
        ("MEP", "MEP"),
    ):
        if ranks[method] > max_rank:
            continue
        energies = (diagrams_by_name.get(diagram_name) or {}).get("energies_kcal", [])
        if len(energies) >= 2:
            try:
                first, last = float(energies[0]), float(energies[-1])
            except (TypeError, ValueError):
                continue
            if np.isfinite(first) and np.isfinite(last):
                overall_rxn_e = round(last - first, 2)
                overall_rxn_method = method
                break

    summary["mlmm_toolkit_version"] = __version__
    summary["schema_version"] = RESULT_JSON_SCHEMA_VERSION
    summary["pipeline_mode"] = pipeline_mode
    summary["status"] = status
    if status_reasons:
        summary["status_reasons"] = status_reasons
    else:
        summary.pop("status_reasons", None)
    # Keep the legacy overloaded ``status`` intact and
    # expose the execution/scientific split plus expected/observed segment IDs so
    # a forward-compatible consumer can tell "the pipeline ran" from "the science
    # is complete and usable". ``scientific_status`` is computed from explicit
    # per-segment LeafOutcomes (IRC directional + endpoint-opt convergence)
    # composed with the legacy completeness axis, so a nonconverged IRC/endpoint
    # leaf whose trajectory still exists cannot make the pipeline a success.
    _apply_pipeline_truth(
        summary,
        post_segments=post_segments,
        config=config,
        legacy_status=status,
        legacy_reasons=status_reasons,
    )
    summary["mlip_backend"] = mlip_backend
    summary["mlip_precision"] = mlip_precision
    # Record the resolved MLIP model used by the high-level ML/MM calculator.
    if mlip_model is not None:
        summary["mlip_model"] = mlip_model
    summary["charge"] = charge
    summary["spin"] = spin
    if manifest is not None:
        summary["run_id"] = manifest.run_id
    summary["n_segments_reactive"] = n_reactive
    if rls:
        summary["rate_limiting_step"] = rls
    else:
        summary.pop("rate_limiting_step", None)
    if overall_rxn_e is not None:
        summary["overall_reaction_energy_kcal"] = overall_rxn_e
        summary["overall_reaction_energy_method"] = overall_rxn_method
    else:
        summary.pop("overall_reaction_energy_kcal", None)
        summary.pop("overall_reaction_energy_method", None)
    if command:
        summary["command"] = command
    if config:
        summary["config"] = config
    citation_config = config or {}
    summary["references"] = method_references(
        {
            "pipeline_mode": pipeline_mode,
            "opt_mode": citation_config.get("opt_mode"),
            "opt_mode_post": citation_config.get("opt_mode_post"),
            "path_opt_mode": citation_config.get("path_opt_mode"),
            "post_opt_mode": citation_config.get("post_opt_mode"),
            "ts_opt_mode": citation_config.get("ts_opt_mode"),
            "endpoint_opt_mode": citation_config.get("endpoint_opt_mode"),
            "mep_mode": citation_config.get("mep_mode"),
            "dmf_correlated": citation_config.get("dmf_correlated"),
            "post_segments": post_segments or [],
        }
    )
    if freeze_atoms:
        summary["freeze_atoms"] = freeze_atoms
    if post_segments:
        summary["post_segments"] = _json_safe(post_segments)

    # Key output file paths for AI agent consumption
    if "out_dir" in summary:
        # Real pipeline root. Fall back to the legacy module_dir.parent for
        # any caller that does not pass out_dir explicitly.
        root = Path(out_dir) if out_dir is not None else Path(summary["out_dir"]).parent
        if manifest is not None:
            # Producer-declared, current-run outputs only (no discovery): a
            # stale file from an earlier invocation is never surfaced.
            current_paths = _current_output_paths(manifest, root)
            key_files: Dict[str, Any] = _current_key_output_files(manifest, root)
        else:
            current_paths = []
            key_files = {}
            # Root-level deliverables (MEP products + authored/mirrored summaries live at root)
            for name, desc in [
                ("summary.log", "Human-readable results summary"),
                ("summary.json", "Machine-readable results summary"),
                ("mep_trj.xyz", "Full MEP trajectory"),
                ("mep.pdb", "Full MEP as PDB"),
                ("mep.cif", "Full MEP with original mmCIF identifiers"),
                ("energy_diagram_MEP.png", "MEP energy plot"),
                ("mep_plot.png", "MEP energy plot (trj2fig)"),
                ("irc_plot_all.png", "Aggregated IRC plot"),
            ]:
                if (root / name).exists():
                    key_files.setdefault(name, desc)
            # Per-segment deliverables under segments/seg_NN/
            seg_parent = root / SEGMENTS_DIRNAME
            if seg_parent.exists():
                for child in sorted(seg_parent.iterdir()):
                    if child.is_dir() and child.name.startswith("seg_"):
                        seg_files = [f.name for f in sorted(child.iterdir()) if f.is_file()]
                        key_files[child.name] = {
                            "description": f"Per-segment results for {child.name}",
                            "files": seg_files,
                        }
        if current_paths:
            summary["current_output_paths"] = current_paths
        else:
            summary.pop("current_output_paths", None)
        if key_files:
            summary["key_output_files"] = key_files
        else:
            summary.pop("key_output_files", None)

    try:
        from mlmm.core.utils import _collect_environment_info
        summary.setdefault("environment", _collect_environment_info())
    except Exception:
        pass

    return summary


def _json_safe(obj):
    """Recursively convert Path objects to strings for JSON serialization."""
    if isinstance(obj, Path):
        return str(obj)
    if isinstance(obj, dict):
        return {k: _json_safe(v) for k, v in obj.items()}
    if isinstance(obj, (list, tuple)):
        return [_json_safe(item) for item in obj]
    return obj


def _required_xyz_block_energies(
    blocks: Sequence[Sequence[str]],
    *,
    path: Path,
    context: str,
) -> List[float]:
    """Return one finite Hartree energy per XYZ block or fail the segment."""
    from mlmm.io.xyz_trajectory import parse_xyz_energy_comment

    energies: List[float] = []
    for frame_index, block in enumerate(blocks, start=1):
        comment = block[1] if len(block) >= 2 else ""
        energy, provenance = parse_xyz_energy_comment(comment)
        if energy is None or not np.isfinite(energy):
            raise click.ClickException(
                f"[all] {context} trajectory {path} has no unambiguous finite "
                f"energy in frame {frame_index} ({provenance}); write "
                "E=<value> with an optional unit."
            )
        energies.append(float(energy))
    return energies


def _commit_summary_json(
    primary: Path,
    summary: Dict[str, Any],
    *,
    mirrors: Sequence[Path] = (),
) -> Path:
    """Publish one aggregate summary generation without relabeling stale data."""

    identified = with_current_run_id(summary)
    summary.clear()
    summary.update(identified)
    return commit_json_exact(primary, summary, mirrors=mirrors)


def _persist_run_manifest(manifest: InvocationManifest, out_dir: Path) -> Path:
    """Persist current-run ownership under the private ``_work`` tree."""

    return manifest.write_internal(
        Path(out_dir) / WORK_DIRNAME / "_run_manifest.json"
    )


def _publish_manifest_summary(
    primary: Path,
    summary: Dict[str, Any],
    *,
    manifest: InvocationManifest,
    out_dir: Path,
    mirrors: Sequence[Path] = (),
) -> Path:
    """Atomically publish, claim, and record one aggregate summary generation.

    The primary destination is the run's public ``summary.json``; the mirrors
    are internal ``_work`` copies and are never tracked as public outputs.
    """

    destination = Path(primary).resolve(strict=False)
    key = _public_output_key(out_dir, destination)
    if key not in manifest.expected:
        manifest.declare(key, [destination])
    summary["run_id"] = manifest.run_id
    published = _commit_summary_json(destination, summary, mirrors=mirrors)
    manifest.claim_one(key)
    _persist_run_manifest(manifest, out_dir)
    return published


def _finalize_current_summary(
    primary: Path,
    summary: Dict[str, Any],
    *,
    manifest: InvocationManifest,
    out_dir: Path,
    mirrors: Sequence[Path] = (),
) -> Path:
    """Republish a summary after every late public producer has finished.

    ``all`` writes and promotes a few public artifacts after its first
    ``summary.json`` generation (notably ``summary.log`` and aggregate plots).
    Refresh ownership after those producers, then publish the final path list.
    A second generation is needed only when the first publication itself adds
    ``summary.json`` to a previously empty manifest.
    """

    def refresh_metadata() -> tuple[list[str], Dict[str, Any]]:
        current_paths = _current_output_paths(manifest, out_dir)
        key_files = _current_key_output_files(manifest, out_dir)
        if current_paths:
            summary["current_output_paths"] = current_paths
        else:
            summary.pop("current_output_paths", None)
        if key_files:
            summary["key_output_files"] = key_files
        else:
            summary.pop("key_output_files", None)
        return current_paths, key_files

    before = refresh_metadata()
    published = _publish_manifest_summary(
        primary,
        summary,
        manifest=manifest,
        out_dir=out_dir,
        mirrors=mirrors,
    )
    after = refresh_metadata()
    if after != before:
        published = _publish_manifest_summary(
            primary,
            summary,
            manifest=manifest,
            out_dir=out_dir,
            mirrors=mirrors,
        )
    _refresh_current_public_outputs(manifest, out_dir)
    _persist_run_manifest(manifest, out_dir)
    return published


def _copy_structures_to_seg_dir(
    state_structs, out_dir, seg_idx, input_suffix,
    prepared_input=None, ref_pdb_path=None, manifest=None,
):
    """Copy R/TS/P structures to out_dir/segments/seg_XX/ in the input format."""
    seg_dir = out_dir / SEGMENTS_DIRNAME / f"seg_{seg_idx:02d}"
    seg_dir.mkdir(parents=True, exist_ok=True)
    name_map = {"R": "reactant", "TS": "ts", "P": "product"}

    def declare_destination(path: Path) -> None:
        if manifest is not None:
            _declare_public_output(manifest, out_dir, path)

    def claim_destination(path: Path) -> None:
        if manifest is not None:
            _claim_public_output(manifest, out_dir, path)

    def copy_destination(src_path: Path, dst_path: Path) -> None:
        declare_destination(dst_path)
        shutil.copy2(src_path, dst_path)
        claim_destination(dst_path)

    for key, src_xyz in state_structs.items():
        src = Path(src_xyz)
        if not src.exists():
            continue
        dst_name = name_map.get(key, key.lower())
        if input_suffix in {".pdb", ".cif", ".mmcif"}:
            src_pdb = src.with_suffix(".pdb")
            if src_pdb.exists():
                dst_pdb = seg_dir / f"{dst_name}.pdb"
                copy_destination(src_pdb, dst_pdb)
                template = coordinate_template_for(src_pdb)
                if template is not None:
                    dst_cif = dst_pdb.with_suffix(".cif")
                    declare_destination(dst_cif)
                    register_output_template_and_write_cif(dst_pdb, template)
                    claim_destination(dst_cif)
                elif src_pdb.with_suffix(".cif").exists():
                    copy_destination(
                        src_pdb.with_suffix(".cif"),
                        dst_pdb.with_suffix(".cif"),
                    )
            else:
                copy_destination(src, seg_dir / f"{dst_name}.xyz")
        elif input_suffix == ".gjf":
            if (prepared_input is not None and getattr(prepared_input, "gjf_template", None) is not None):
                dst_gjf = seg_dir / f"{dst_name}.gjf"
                declare_destination(dst_gjf)
                try:
                    from mlmm.core.utils import convert_xyz_to_gjf
                    convert_xyz_to_gjf(src, prepared_input.gjf_template, dst_gjf)
                    claim_destination(dst_gjf)
                except Exception:
                    copy_destination(src, seg_dir / f"{dst_name}.xyz")
            else:
                copy_destination(src, seg_dir / f"{dst_name}.xyz")
        else:
            copy_destination(src, seg_dir / f"{dst_name}.xyz")
    return seg_dir


def _read_summary(summary_json: Path) -> List[Dict[str, Any]]:
    """
    Read path_search/summary.json and return segments list (empty if not found).
    """
    try:
        if not summary_json.exists():
            return []
        data = json.loads(summary_json.read_text(encoding="utf-8")) or {}
        segs = data.get("segments", []) or []
        if not isinstance(segs, list):
            return []
        return segs
    except Exception:
        return []


def _pdb_models_to_coords_and_elems(pdb_path: Path) -> Tuple[List[np.ndarray], List[str]]:
    """
    Return ([coords_model1, coords_model2, ...] in Å), [elements] from a multi-model PDB.
    """
    parser = PDB.PDBParser(QUIET=True)
    st = parser.get_structure("seg", str(pdb_path))
    models = list(st.get_models())
    if not models:
        raise click.ClickException(f"[post] No MODEL found in PDB: {pdb_path}")
    # atom order taken from first model
    atoms0 = [a for a in models[0].get_atoms()]
    elems: List[str] = []
    for a in atoms0:
        el = (a.element or "").strip()
        if not el:
            # fall back: derive from atom name
            nm = a.get_name().strip()
            el = "".join([c for c in nm if c.isalpha()])[:2].title() or "C"
        elems.append(el)
    coords_list: List[np.ndarray] = []
    for m in models:
        atoms = [a for a in m.get_atoms()]
        if len(atoms) != len(atoms0):
            raise click.ClickException(f"[post] Atom count mismatch across models in {pdb_path}")
        coords = np.array([a.get_coord() for a in atoms], dtype=float)
        coords_list.append(coords)
    return coords_list, elems


def _geom_from_angstrom(elems: Sequence[str],
                        coords_ang: np.ndarray,
                        freeze_atoms: Sequence[int]) -> Any:
    """
    Create a Geometry from Å coordinates using _path_search._new_geom_from_coords (expects Bohr).
    """
    coords_bohr = np.asarray(coords_ang, dtype=float) / BOHR2ANG
    return _path_search._new_geom_from_coords(elems, coords_bohr, coord_type="cart", freeze_atoms=freeze_atoms)


def _load_segment_end_geoms(seg_pdb: Path, freeze_atoms: Sequence[int]) -> Tuple[Any, Any]:
    """
    Load first/last model as Geometries from a per-segment pocket PDB.
    """
    coords_list, elems = _pdb_models_to_coords_and_elems(seg_pdb)
    gL = _geom_from_angstrom(elems, coords_list[0], freeze_atoms)
    gR = _geom_from_angstrom(elems, coords_list[-1], freeze_atoms)
    return gL, gR


def _orient_irc_endpoint_geometries(
    g_left: Any,
    g_right: Any,
    g_mep_left: Any,
    g_mep_right: Any,
) -> Tuple[Any, Any, str, str, bool, Dict[str, Any]]:
    """Orient IRC endpoints; topology decides only an unambiguous XOR."""

    bond_cfg = dict(_path_search.BOND_KW)

    def _matches(x: Any, y: Any) -> bool:
        try:
            changed, _ = _path_search._has_bond_change(x, y, bond_cfg)
            return not changed
        except Exception:
            return False

    def _rmsd_cart(g1: Any, g2: Any) -> float:
        c1 = np.asarray(g1.coords).reshape(-1, 3)
        c2 = np.asarray(g2.coords).reshape(-1, 3)
        n = min(len(c1), len(c2))
        return float(np.sqrt(np.mean((c1[:n] - c2[:n]) ** 2)))

    match_LL = _matches(g_left, g_mep_left)
    match_LR = _matches(g_left, g_mep_right)
    match_RL = _matches(g_right, g_mep_left)
    match_RR = _matches(g_right, g_mep_right)
    direct_match = bool(match_LL and match_RR)
    swapped_match = bool(match_LR and match_RL)
    reverse_irc = False
    assignment: Dict[str, Any] = {
        "match_matrix": {
            "left_to_mep_left": bool(match_LL),
            "left_to_mep_right": bool(match_LR),
            "right_to_mep_left": bool(match_RL),
            "right_to_mep_right": bool(match_RR),
        },
        "reversed": False,
    }

    if direct_match ^ swapped_match:
        assignment["method"] = "bond_topology"
        assignment["connectivity_validated"] = True
        reverse_irc = swapped_match
    else:
        assignment["method"] = (
            "rmsd_topology_tie"
            if direct_match and swapped_match
            else "rmsd_topology_unmatched"
        )
        assignment["connectivity_validated"] = bool(
            direct_match and swapped_match
        )
        try:
            direct_score = float(
                _rmsd_cart(g_left, g_mep_left)
                + _rmsd_cart(g_right, g_mep_right)
            )
            swapped_score = float(
                _rmsd_cart(g_left, g_mep_right)
                + _rmsd_cart(g_right, g_mep_left)
            )
            assignment["rmsd_direct"] = direct_score
            assignment["rmsd_swapped"] = swapped_score
            reverse_irc = swapped_score < direct_score
        except Exception as exc:
            assignment["method"] = "unresolved"
            assignment["connectivity_validated"] = False
            assignment["reason"] = f"rmsd_failed:{exc}"
        if not assignment["connectivity_validated"]:
            assignment.setdefault(
                "reason",
                "neither endpoint pairing matched the MEP bond topology",
            )

    if reverse_irc:
        g_left, g_right = g_right, g_left
        left_tag, right_tag = "backward", "forward"
    else:
        left_tag, right_tag = "forward", "backward"
    assignment["reversed"] = bool(reverse_irc)
    return (
        g_left,
        g_right,
        left_tag,
        right_tag,
        reverse_irc,
        assignment,
    )


def _irc_and_match(seg_idx: int,
                   seg_dir: Path,
                   ref_pdb_for_seg: Path,
                   seg_pocket_pdb: Path,
                   g_ts: Any,
                   q_int: int,
                   spin: int,
                   *,
                   resolved_calc_template: _ResolvedCalculatorTemplate,
                   mep_dir: Optional[Path] = None,
                   real_parm7: Optional[Path] = None,
                   model_pdb: Optional[Path] = None,
                   detect_layer: bool = False,
                   backend: Optional[str] = None,
                   embedcharge: bool = False,
                   embedcharge_cutoff: Optional[float] = None,
                   embedcharge_explicit: bool = False,
                   link_atom_method: Optional[str] = None,
                   mm_backend: Optional[str] = None,
                   use_cmap: Optional[bool] = None,
                   irc_step_size: Optional[float] = None,
                   irc_never_stop: Optional[bool] = None,
                   session: Optional[RunSession] = None,
                   args_yaml: Optional[Path] = None) -> Dict[str, Any]:
    """
    Run EulerPC IRC from a TS geometry, then map the IRC endpoints to (left, right)
    by comparing bond states with the MEP segment endpoints (when available).
    Falls back to raw IRC orientation in TSOPT-only mode.

    Endpoint matching logic (when MEP endpoints exist):
      - Compute bond change sets at IRC's two endpoints (`bond_changes.compare_structures`).
      - Score each pairing (IRC.fwd, IRC.bwd) ↔ (MEP.left, MEP.right) by symmetric-diff
        bond change count, pick the orientation with minimum total diff.
      - On tie, prefer the orientation whose forward endpoint shares more atoms
        with the MEP reactant side (= side selected by `seg_idx`-based ordering convention).

    TSOPT-only fallback: when no MEP endpoints (= TS-only pipeline), IRC's raw
    forward/backward orientation is preserved as (left, right) without remapping;
    the caller can post-hoc swap if needed.

    GPU memory handling: caller-supplied `g_ts` pins TS-stage allocator pages.
    The fix-C `gc.collect()` + `torch.cuda.empty_cache()` at function entry frees
    the previous stage's residency so IRC's `initial_displacement eigh` (large
    contiguous ~9 GiB block) succeeds.
    """
    # Fix C: free GPU memory carried over from the preceding TS-opt stage
    # before IRC. IRC's initial_displacement eigh needs a large contiguous
    # block; ~9 GiB of TS-stage allocator residency otherwise leaves too little
    # free even after the Fix-A ML-macro Hessian reduction. The per-stage
    # _run_cli_main finally also does this, but orchestrator locals (g_ts, etc.)
    # still pin memory here at the TS->IRC boundary.
    gc.collect()
    if torch.cuda.is_available():
        torch.cuda.empty_cache()

    irc_dir = seg_dir / "irc"
    ensure_dir(irc_dir)

    # Build irc CLI arguments
    irc_args: List[str] = [
        "-i", str(ref_pdb_for_seg),
        "--parm", str(real_parm7),
        "--model-pdb", str(model_pdb),
        "-q", str(int(q_int)),
        "-m", str(int(spin)),
        "--out-dir", str(irc_dir),
    ]
    irc_args.append("--detect-layer" if detect_layer else "--no-detect-layer")
    if irc_step_size is not None:
        irc_args.extend(["--step-size", str(float(irc_step_size))])
    if irc_never_stop is not None:
        irc_args.append(
            "--never-stop" if irc_never_stop else "--no-never-stop"
        )
    from mlmm.workflows._all_helpers import append_backend_forwarding_args
    append_backend_forwarding_args(
        irc_args,
        backend=backend,
        embedcharge=embedcharge,
        embedcharge_cutoff=embedcharge_cutoff,
        embedcharge_explicit=embedcharge_explicit,
        link_atom_method=link_atom_method,
        mm_backend=mm_backend,
        use_cmap=use_cmap,
        args_yaml=args_yaml,
    )
    # request the child's machine-readable result.json so the aggregate
    # can gate on reported per-direction IRC convergence instead of trajectory-
    # file existence. A never_stop / max-cycle direction still writes its
    # trajectory, but the child reports it as nonconverged and we must not promote
    # it.
    irc_args.append("--out-json")

    _echo_detail(f"[irc] Running EulerPC IRC → out={irc_dir}")
    try:
        _run_cli_main("irc", _irc_cli.cli, irc_args, on_nonzero="raise", prefix="irc")
    except BaseException:
        # Make the consequence explicit instead of dying mid-recovery with a
        # bare stack: IRC is a hard prerequisite for this segment's IRC
        # endpoint Hessians and the subsequent freq/thermo/DFT, so those are
        # not produced for this segment when IRC fails.
        _echo(
            f"[all] IRC failed for segment {seg_idx}; freq/thermochemistry/DFT "
            f"post-processing for this segment is skipped (the pipeline will "
            f"now abort with the IRC error above).",
            err=True,
        )
        raise

    # Read IRC endpoints
    finished_trj = irc_dir / "finished_irc_trj.xyz"
    finished_pdb = irc_dir / "finished_irc.pdb"
    irc_plot = irc_dir / "irc_plot.png"

    if not finished_trj.exists():
        raise click.ClickException(f"[irc] IRC trajectory not found: {finished_trj}")

    # Convert to PDB if not already done
    if not finished_pdb.exists():
        _path_search._maybe_convert_to_pdb(finished_trj, ref_pdb_path=seg_pocket_pdb, out_path=finished_pdb)

    elems, c_first, c_last = read_xyz_first_last(finished_trj)

    # Create geometries from IRC endpoints
    _irc_calc_kwargs = _stage_calc_kwargs(
        resolved_calc_template,
        input_pdb=ref_pdb_for_seg,
        real_parm7=real_parm7,
        model_pdb=model_pdb,
        charge=q_int,
        spin=spin,
        use_bfactor_layers=detect_layer,
    )
    # One heavy ML/MM core serves this segment's IRC endpoints; the parent owns
    # its release at the phase
    # boundary via the returned lease.
    calc = _mlmm_calc(**_irc_calc_kwargs)
    lease = CalculatorLease(calc)
    if session is not None:
        session.resources.add(lease.release)
    try:

        g_left = _path_search._new_geom_from_coords(
            elems, c_first / BOHR2ANG, coord_type="cart", freeze_atoms=[])
        g_right = _path_search._new_geom_from_coords(
            elems, c_last / BOHR2ANG, coord_type="cart", freeze_atoms=[])
        lease.attach(g_left)
        lease.attach(g_right)
        _ = float(g_left.energy)
        _ = float(g_right.energy)

        # Reload TS geometry with energy
        if g_ts.calculator is None:
            lease.attach(g_ts)
        _ = float(g_ts.energy)

        left_tag = "forward"
        right_tag = "backward"
        reverse_irc = False
        expects_mep_endpoints = mep_dir is not None
        endpoint_assignment: Dict[str, Any] = {
            "method": "raw",
            "reversed": False,
            "connectivity_validated": (
                False if expects_mep_endpoints else None
            ),
        }

        # Try to load segment endpoints for mapping.
        # mep_seg_NN.pdb is written by the MEP engine under path_dir (now _work/path_*);
        # seg_dir moved to segments/, so read from mep_dir when provided.
        gL_end = None
        gR_end = None
        mep_root = mep_dir if mep_dir is not None else seg_dir.parent
        seg_pocket_path = mep_root / f"mep_seg_{seg_idx:02d}.pdb"
        if seg_pocket_path.exists():
            try:
                gL_end, gR_end = _load_segment_end_geoms(seg_pocket_path, [])
            except Exception as e:
                endpoint_assignment["method"] = "unresolved"
                endpoint_assignment["reason"] = f"mep_endpoint_load_failed:{e}"
                click.echo(f"[post] WARNING: failed to load segment endpoints: {e}", err=True)
        elif expects_mep_endpoints:
            endpoint_assignment["reason"] = "mep_endpoint_trajectory_missing"

        # Map IRC endpoints to left/right using bond-change analysis
        if gL_end is not None and gR_end is not None:
            (
                g_left,
                g_right,
                left_tag,
                right_tag,
                reverse_irc,
                endpoint_assignment,
            ) = _orient_irc_endpoint_geometries(
                g_left,
                g_right,
                gL_end,
                gR_end,
            )

        return {
            "left_min_geom": g_left,
            "right_min_geom": g_right,
            "ts_geom": g_ts,
            "left_tag": left_tag,
            "right_tag": right_tag,
            "irc_trj": str(finished_trj) if finished_trj.exists() else None,
            "irc_plot": str(irc_plot) if irc_plot.exists() else None,
            "reverse_irc": reverse_irc,
            "endpoint_assignment": endpoint_assignment,
            "calculator_lease": lease,
            # the child's per-direction convergence. The IRC leaf
            # is usable only when EVERY requested direction explicitly converged; a
            # trajectory can exist for a nonconverged (never_stop / max-cycle)
            # direction, so aggregate promotion must gate on this, not file
            # existence.
            "irc_outcome": _read_irc_outcome(irc_dir),
        }
    except BaseException:
        lease.release()
        raise


def _save_single_geom_for_tools(g: Any, ref_pdb: Path, out_dir: Path, name: str) -> Tuple[Path, Path]:
    """
    Write XYZ (primary, full precision) + PDB (companion) for a single geometry.
    Returns (xyz_path, pdb_path).
    """
    out_dir.mkdir(parents=True, exist_ok=True)
    # XYZ — full precision
    xyz_out = out_dir / f"{name}.xyz"
    with open(xyz_out, "w") as f:
        f.write(g.as_xyz() + "\n")
    # TRJ with energy (for PDB conversion and trajectory viewers)
    xyz_trj = out_dir / f"{name}_trj.xyz"
    _path_search._write_xyz_trj_with_energy([g], [float(g.energy)], xyz_trj)
    # PDB companion
    pdb_out = out_dir / f"{name}.pdb"
    _path_search._maybe_convert_to_pdb(xyz_trj, ref_pdb_path=ref_pdb, out_path=pdb_out)
    return xyz_out, pdb_out


def _validate_tsopt_result_payload(
    payload: Dict[str, Any], *, skip_final_freq: bool
) -> None:
    """Reject a TS result that is not a verified first-order saddle."""
    status = str(payload.get("status") or "unknown")
    n_imag = payload.get("n_imaginary_modes")
    if status == "unverified" and skip_final_freq:
        return
    if status != "converged" or n_imag != 1:
        raise click.ClickException(
            "[tsopt] TS optimization did not produce a validated first-order "
            f"saddle (status={status!r}, n_imag={n_imag!r}); IRC was not started."
        )


def _run_tsopt_on_hei(hei_pdb: Path,
                      charge: int,
                      spin: int,
                      real_parm7: Path,
                      model_pdb: Path,
                      detect_layer: bool,
                      args_yaml: Optional[Path],
                      out_dir: Path,
                      opt_mode_default: str,
                      *,
                      resolved_calc_template: _ResolvedCalculatorTemplate,
                      overrides: Optional[Dict[str, Any]] = None,
                      backend: Optional[str] = None,
                      embedcharge: bool = False,
                      embedcharge_cutoff: Optional[float] = None,
                      embedcharge_explicit: bool = False,
                      link_atom_method: Optional[str] = None,
                      mm_backend: Optional[str] = None,
                      use_cmap: Optional[bool] = None,
                      ref_pdb: Optional[Path] = None) -> Tuple[Path, Any]:
    """
    Run tsopt CLI on a HEI structure; return (final_ts_pdb_path, ts_geom).

    When *ref_pdb* (layered PDB with B-factor layer info) is given, the HEI XYZ
    is used as input and *ref_pdb* is passed via ``--ref-pdb`` so that the
    calculator correctly detects ML/MM layers from B-factors.
    """
    overrides = overrides or {}
    # Prefer HEI XYZ (full precision) + layered ref-pdb (B-factor layer info)
    hei_xyz = hei_pdb.with_suffix(".xyz")
    if ref_pdb is not None and hei_xyz.exists():
        input_file = hei_xyz
        topology_pdb = ref_pdb
    else:
        input_file = hei_pdb
        topology_pdb = hei_pdb
    prepared_input = prepare_input_structure(input_file)
    if input_file.suffix.lower() == ".xyz" and ref_pdb is not None:
        apply_ref_pdb_override(prepared_input, ref_pdb)
    try:
        ts_dir = _resolve_override_dir(out_dir / "ts", overrides.get("out_dir"))
        ensure_dir(ts_dir)

        opt_mode = overrides.get("opt_mode", opt_mode_default)

        ts_args: List[str] = [
            "-i", str(prepared_input.geom_path),
        ]
        if input_file.suffix.lower() == ".xyz" and ref_pdb is not None:
            ts_args.extend(["--ref-pdb", str(ref_pdb)])
        ts_args.extend([
            "--parm", str(real_parm7),
            "--model-pdb", str(model_pdb),
            "-q", str(int(charge)),
            "-m", str(int(spin)),
            "--out-dir", str(ts_dir),
        ])
        ts_args.append("--detect-layer" if detect_layer else "--no-detect-layer")

        if opt_mode is not None:
            ts_args.extend(["--opt-mode", str(opt_mode)])

        reference_mode = overrides.get("reference_mode")
        if reference_mode is not None:
            ts_args.extend(["--ref-mode", str(reference_mode)])

        _append_cli_arg(ts_args, "--max-cycles", overrides.get("max_cycles"))
        _append_toggle_arg(ts_args, "--dump", overrides.get("dump"))
        _append_toggle_arg(ts_args, "--convert-files", overrides.get("convert_files"))
        _append_cli_arg(ts_args, "--thresh", overrides.get("thresh"))
        _append_toggle_arg(ts_args, "--flatten", overrides.get("flatten"))

        hess_mode = overrides.get("hessian_calc_mode")
        if hess_mode:
            ts_args.extend(["--hessian-calc-mode", str(hess_mode)])

        if args_yaml is not None:
            ts_args.extend(["--config", str(args_yaml)])

        from mlmm.workflows._all_helpers import append_backend_forwarding_args
        append_backend_forwarding_args(
            ts_args,
            backend=backend,
            embedcharge=embedcharge,
            embedcharge_cutoff=embedcharge_cutoff,
            embedcharge_explicit=embedcharge_explicit,
            link_atom_method=link_atom_method,
            mm_backend=mm_backend,
            use_cmap=use_cmap,
        )
        _append_toggle_arg(
            ts_args, "--skip-final-freq", overrides.get("skip_final_freq")
        )
        ts_args.append("--out-json")

        _echo_detail(f"[tsopt] Running tsopt on HEI → out={ts_dir}")
        _run_cli_main("tsopt", _ts_opt.cli, ts_args, on_nonzero="raise", prefix="tsopt")

        result_path = ts_dir / "result.json"
        if not result_path.exists():
            raise click.ClickException(
                f"[tsopt] Missing machine-readable TS validation result: {result_path}"
            )
        try:
            tsopt_result = json.loads(result_path.read_text(encoding="utf-8"))
        except (OSError, json.JSONDecodeError) as exc:
            raise click.ClickException(
                f"[tsopt] Could not read TS validation result '{result_path}': {exc}"
            ) from exc
        ts_status = str(tsopt_result.get("status") or "unknown")
        if ts_status == "unverified" and overrides.get("skip_final_freq"):
            _echo(
                "[tsopt] WARNING: saddle order is unverified because final frequency "
                "analysis was explicitly skipped.",
                err=True,
            )
        _validate_tsopt_result_payload(
            tsopt_result,
            skip_final_freq=bool(overrides.get("skip_final_freq")),
        )

        # Prefer XYZ (full precision) for geometry loading; PDB for topology
        final_xyz = ts_dir / "final_geometry.xyz"
        ts_pdb = ts_dir / "final_geometry.pdb"
        if not ts_pdb.exists() and final_xyz.exists():
            _path_search._maybe_convert_to_pdb(final_xyz, topology_pdb, ts_pdb)
        if not final_xyz.exists() and not ts_pdb.exists():
            raise click.ClickException("[tsopt] TS outputs not found.")
        geom_src = final_xyz if final_xyz.exists() else ts_pdb
        g_ts = geom_loader(geom_src, coord_type="cart")
        g_ts._tsopt_result = tsopt_result

        # Ensure calculator to have energy on g_ts
        _ts_calc_kwargs = _stage_calc_kwargs(
            resolved_calc_template,
            input_pdb=topology_pdb,
            real_parm7=real_parm7,
            model_pdb=model_pdb,
            charge=charge,
            spin=spin,
            use_bfactor_layers=detect_layer,
        )
        calc = _mlmm_calc(**_ts_calc_kwargs)
        g_ts.set_calculator(calc)
        _ = float(g_ts.energy)

        # release this probe core before returning.  The caller's next
        # phase builds its own leased core and attaches ``g_ts`` to it
        # (``_irc_and_match`` -> ``if g_ts.calculator is None: lease.attach``);
        # leaving this one attached both defeats that lease and keeps two heavy
        # ML/MM cores resident across the TSOPT->IRC handoff.  Detach by direct
        # assignment, not ``set_calculator(None)``, which would ``clear`` the
        # energy just computed.  Mirrors ``CalculatorLease.release``.
        g_ts.calculator = None
        _close = getattr(calc, "close", None)
        if callable(_close):
            try:
                _close()
            except Exception as exc:  # pragma: no cover - defensive
                logger.debug("TS probe calculator close failed: %s", exc)
        del calc
        gc.collect()
        try:
            if torch.cuda.is_available():
                torch.cuda.empty_cache()
        except Exception as exc:  # pragma: no cover - defensive
            logger.debug("CUDA cache release unavailable: %s", exc)

        return ts_pdb, g_ts
    finally:
        prepared_input.cleanup()


def _ensure_hei_path_tangent(
    mep_trj: Path,
    hei_path: Path,
    mode_path: Path,
) -> Optional[Path]:
    """Write the MEP tangent at the trajectory image matching the HEI."""
    if not mep_trj.exists() or not hei_path.exists():
        return None
    try:
        from ase.io import read as ase_read

        images = list(ase_read(str(mep_trj), index=":"))
        hei = ase_read(str(hei_path), index=0)
        if len(images) < 2:
            return None
        hei_numbers = np.asarray(hei.numbers)
        if any(
            len(image) != len(hei)
            or not np.array_equal(np.asarray(image.numbers), hei_numbers)
            for image in images
        ):
            return None
        hei_positions = np.asarray(hei.positions, dtype=float)
        image_index = int(
            np.argmin(
                [
                    float(np.linalg.norm(np.asarray(image.positions) - hei_positions))
                    for image in images
                ]
            )
        )
        energies = None
        blocks = read_xyz_as_blocks(mep_trj)
        if len(blocks) == len(images):
            parsed = []
            for block in blocks:
                from mlmm.io.xyz_trajectory import parse_xyz_energy_comment

                comment = block[1] if len(block) >= 2 else ""
                energy, _provenance = parse_xyz_energy_comment(comment)
                parsed.append(energy)
            if all(
                energy is not None and np.isfinite(energy)
                for energy in parsed
            ):
                energies = parsed
        tangent = _path_search._normalized_path_tangent(
            [np.asarray(image.positions, dtype=float) for image in images],
            image_index,
            energies=energies,
        )
        if tangent is None:
            return None
        mode_path.parent.mkdir(parents=True, exist_ok=True)
        np.savetxt(mode_path, tangent, fmt="%.17e")
        _echo(f"[tsopt] Derived HEI path-tangent reference mode → {mode_path}")
        return mode_path
    except Exception as exc:
        _echo(f"[tsopt] WARNING: Could not derive HEI path tangent: {exc}", err=True)
        return None


def _write_segment_energy_diagram(
    prefix: Path,
    labels: List[str],
    energies_eh: List[float],
    title_note: str,
    ylabel: str = "ΔE (kcal/mol)",
    write_html: bool = False,
) -> Optional[Dict[str, Any]]:
    """
    Write energy diagram (PNG only) using utils.build_energy_diagram.
    """
    if not energies_eh:
        return None
    e0 = energies_eh[0]
    energies_kcal = [(e - e0) * AU2KCALPERMOL for e in energies_eh]
    fig = build_energy_diagram(
        energies=energies_kcal,
        labels=labels,
        ylabel=ylabel,
        baseline=True,
        showgrid=False,
    )
    if title_note:
        fig.update_layout(title=title_note)
    png = prefix.with_suffix(".png")
    try:
        fig.write_image(str(png), scale=2)
    except Exception as e:
        click.echo(f"[diagram] NOTE: PNG export skipped (install 'kaleido' to enable): {e}", err=True)
    else:
        emit(f"[diagram] Wrote energy diagram → {png.name}", detail=True)

    payload: Dict[str, Any] = {
        "name": prefix.stem,
        "labels": labels,
        "energies_kcal": energies_kcal,
        "ylabel": ylabel,
        "energies_au": list(energies_eh),
        "image": str(png),
    }
    if title_note:
        payload["title"] = title_note
    return payload


def _build_global_segment_labels(n_segments: int) -> List[str]:
    """
    Build R/TS/P labels for an aggregated multi-segment MEP diagram.

    Pattern:
      - n = 1: ["R", "TS1", "P"]
      - n >= 2: R, TS1, IM1_1, IM1_2, TS2, IM2_1, IM2_2, ..., TSN, P
    """
    if n_segments <= 0:
        return []
    if n_segments == 1:
        return ["R", "TS1", "P"]

    labels: List[str] = []
    for seg_idx in range(1, n_segments + 1):
        if seg_idx == 1:
            labels.extend(["R", "TS1", "IM1_1"])
        elif seg_idx == n_segments:
            labels.extend([f"IM{seg_idx - 1}_2", f"TS{seg_idx}", "P"])
        else:
            labels.extend(
                [f"IM{seg_idx - 1}_2", f"TS{seg_idx}", f"IM{seg_idx}_1"]
            )
    return labels


def _merge_irc_trajectories_to_single_plot(
    trj_and_flags: Sequence[Tuple[Path, bool]],
    out_png: Path,
) -> None:
    """
    Build a single IRC plot over all reactive segments using trj2fig.
    """
    all_blocks: List[str] = []
    for trj_path, reverse in trj_and_flags:
        if not isinstance(trj_path, Path) or not trj_path.exists():
            continue
        try:
            blocks = read_xyz_as_blocks(trj_path)
        except click.ClickException as e:
            click.echo(str(e), err=True)
            continue
        if not blocks:
            continue
        if reverse:
            blocks = list(reversed(blocks))
        all_blocks.extend("\n".join(b) for b in blocks)

    if not all_blocks:
        return

    tmp_trj = out_png.with_name(f"{out_png.stem}_trj.xyz")
    ensure_dir(tmp_trj.parent)
    try:
        tmp_trj.write_text("\n".join(all_blocks) + "\n", encoding="utf-8")
    except Exception as e:
        click.echo(f"[irc_all] WARNING: Failed to write concatenated IRC trajectory: {e}", err=True)
        return

    try:
        run_trj2fig(tmp_trj, [out_png], unit="kcal", reference="init", reverse_x=False)
        click.echo(f"[irc_all] Wrote aggregated IRC plot → {out_png}")
    except Exception as e:
        click.echo(f"[irc_all] WARNING: failed to plot concatenated IRC trajectory: {e}", err=True)
    finally:
        try:
            tmp_trj.unlink()
        except Exception:
            logger.debug("Failed to unlink temp trajectory file", exc_info=True)


def _run_freq_for_state(pdb_path: Path,
                        q_int: int,
                        spin: int,
                        real_parm7: Path,
                        model_pdb: Path,
                        detect_layer: bool,
                        out_dir: Path,
                        args_yaml: Optional[Path],
                        overrides: Optional[Dict[str, Any]] = None,
                        backend: Optional[str] = None,
                        embedcharge: bool = False,
                        embedcharge_cutoff: Optional[float] = None,
                        embedcharge_explicit: bool = False,
                        link_atom_method: Optional[str] = None,
                        mm_backend: Optional[str] = None,
                        use_cmap: Optional[bool] = None,
                        xyz_path: Optional[Path] = None) -> Dict[str, Any]:
    """
    Run freq CLI; return parsed thermo dict (may be empty).
    When *xyz_path* is given, use it for full-precision coordinates with
    *pdb_path* as topology reference (--ref-pdb).
    """
    fdir = out_dir
    ensure_dir(fdir)
    overrides = overrides or {}

    dump_use = overrides.get("dump")
    # `all --thermo` assembles its Gibbs diagram from the child freq stage's
    # thermoanalysis.yaml, which freq writes only under --dump. Default the
    # child to dump so requesting thermochemistry actually yields it; an
    # explicit --no-dump still reaches the child through *overrides*.
    if dump_use is None:
        dump_use = True

    # Prefer XYZ (full precision) with --ref-pdb for topology
    if xyz_path is not None and xyz_path.exists():
        args = ["-i", str(xyz_path), "--ref-pdb", str(pdb_path)]
    else:
        args = ["-i", str(pdb_path)]
    args.extend([
        "--parm", str(real_parm7),
        "--model-pdb", str(model_pdb),
        "-q", str(int(q_int)),
        "-m", str(int(spin)),
        "--out-dir", str(fdir),
    ])
    args.append("--detect-layer" if detect_layer else "--no-detect-layer")

    _append_cli_arg(args, "--max-write", overrides.get("max_write"))
    _append_cli_arg(args, "--amplitude-ang", overrides.get("amplitude_ang"))
    _append_cli_arg(args, "--n-frames", overrides.get("n_frames"))
    if overrides.get("sort") is not None:
        args.extend(["--sort", str(overrides.get("sort"))])
    _append_cli_arg(args, "--temperature", overrides.get("temperature"))
    _append_cli_arg(args, "--pressure", overrides.get("pressure"))
    _append_cli_arg(args, "--symmetry-number", overrides.get("symmetry_number"))
    _append_toggle_arg(args, "--dump", dump_use)
    _append_toggle_arg(args, "--convert-files", overrides.get("convert_files"))

    hess_mode = overrides.get("hessian_calc_mode")
    if hess_mode:
        args.extend(["--hessian-calc-mode", str(hess_mode)])

    from mlmm.workflows._all_helpers import append_backend_forwarding_args
    append_backend_forwarding_args(
        args,
        backend=backend,
        embedcharge=embedcharge,
        embedcharge_cutoff=embedcharge_cutoff,
        embedcharge_explicit=embedcharge_explicit,
        link_atom_method=link_atom_method,
        mm_backend=mm_backend,
        use_cmap=use_cmap,
        args_yaml=args_yaml,
    )
    _freq_rc = _run_cli_main("freq", _freq_cli.cli, args, on_nonzero="warn", on_exception="raise", prefix="freq")
    # a nonzero freq exit means the thermochemistry is NOT usable, even if
    # a thermoanalysis.yaml (from a prior run or a partial write) exists with
    # finite fields. Never infer FREQ success from the filename or a finite
    # number — return {} so no Gibbs diagram/dict can be built from it.
    if _freq_rc not in (None, 0):
        _echo(
            f"[freq] WARNING: freq exited with code {_freq_rc}; thermochemistry is "
            "unusable and will not enter any Gibbs diagram.",
            err=True,
        )
        return {}
    if not bool(dump_use):
        return {}
    # parse thermoanalysis.yaml if any
    y = fdir / "thermoanalysis.yaml"
    if y.exists():
        try:
            return yaml.safe_load(y.read_text(encoding="utf-8")) or {}
        except Exception:
            return {}
    return {}


def _thermo_gibbs_ha(payload: Any) -> Optional[float]:
    """Return ``sum_EE_and_thermal_free_energy_ha`` only when finite; else None.

    A missing/nonfinite FREQ free-energy field must never be replaced by
    a MLIP electronic energy in a Gibbs-named result. The caller builds a Gibbs
    diagram/dict only when EVERY requested state returns a finite value here.
    """

    if not isinstance(payload, dict):
        return None
    val = payload.get("sum_EE_and_thermal_free_energy_ha")
    try:
        f = float(val)
    except (TypeError, ValueError):
        return None
    return f if math.isfinite(f) else None


def _thermo_correction_ha(payload: Any) -> Optional[float]:
    """Return ``thermal_correction_free_energy_ha`` only when finite; else None.

    A missing/nonfinite thermal correction must never be replaced by 0.0
    in a DFT//MLIP/MM Gibbs result (that would silently report the electronic DFT
    energy as a Gibbs free energy).
    """

    if not isinstance(payload, dict):
        return None
    val = payload.get("thermal_correction_free_energy_ha")
    try:
        f = float(val)
    except (TypeError, ValueError):
        return None
    return f if math.isfinite(f) else None


def _run_opt_for_state(
    pdb_path: Path,
    q_int: int,
    spin: int,
    real_parm7: Path,
    model_pdb: Path,
    detect_layer: bool,
    out_dir: Path,
    args_yaml: Optional[Path],
    opt_mode_default: str,
    *,
    resolved_calc_template: _ResolvedCalculatorTemplate,
    convert_files: Optional[bool] = None,
    backend: Optional[str] = None,
    embedcharge: bool = False,
    embedcharge_cutoff: Optional[float] = None,
    embedcharge_explicit: bool = False,
    link_atom_method: Optional[str] = None,
    mm_backend: Optional[str] = None,
    use_cmap: Optional[bool] = None,
    thresh: Optional[str] = None,
    reject_uphill: Optional[bool] = None,
    xyz_path: Optional[Path] = None,
) -> Tuple[Any, Path, Optional[bool]]:
    """
    Run opt CLI for a single endpoint and return
    ``(optimized Geometry, final geometry path, converged)``.

    ``converged`` is the fail-closed tri-state convergence bit of the endpoint
    opt child, read from its ``result.json``: an endpoint whose
    optimization did not explicitly converge is still retained as a geometry /
    artifact but must not promote its segment to a usable success. ``--out-json``
    is forced on so that the convergence bit is always emitted.

    When *xyz_path* is given, pass it as ``-i`` with ``--ref-pdb pdb_path`` to
    preserve full coordinate precision.
    """
    opt_dir = out_dir
    ensure_dir(opt_dir)

    # Use XYZ (full precision) when available; fall back to PDB
    if xyz_path is not None and xyz_path.exists():
        prepared_input = prepare_input_structure(xyz_path)
        apply_ref_pdb_override(prepared_input, pdb_path)
        input_label = xyz_path.name
    else:
        prepared_input = prepare_input_structure(pdb_path)
        input_label = pdb_path.name
    try:
        opt_mode = str(opt_mode_default or "heavy").lower()
        args = [
            "-i", str(prepared_input.geom_path),
        ]
        # Add --ref-pdb when input is XYZ
        if prepared_input.geom_path.suffix.lower() == ".xyz":
            args.extend(["--ref-pdb", str(prepared_input.source_path)])
        args.extend([
            "--parm", str(real_parm7),
            "--model-pdb", str(model_pdb),
            "-q", str(int(q_int)),
            "-m", str(int(spin)),
            "--out-dir", str(opt_dir),
            "--opt-mode", opt_mode,
            # Emit result.json so the endpoint opt child's explicit
            # convergence bit can gate the segment (never inferred from files).
            "--out-json",
        ])
        args.append("--detect-layer" if detect_layer else "--no-detect-layer")
        _append_toggle_arg(args, "--convert-files", convert_files)
        _append_cli_arg(args, "--thresh", thresh)
        _append_toggle_arg(args, "--reject-uphill", reject_uphill)

        if args_yaml is not None:
            args.extend(["--config", str(args_yaml)])

        from mlmm.workflows._all_helpers import append_backend_forwarding_args
        append_backend_forwarding_args(
            args,
            backend=backend,
            embedcharge=embedcharge,
            embedcharge_cutoff=embedcharge_cutoff,
            embedcharge_explicit=embedcharge_explicit,
            link_atom_method=link_atom_method,
            mm_backend=mm_backend,
            use_cmap=use_cmap,
        )

        _echo_detail(f"[endpoint-opt] Running opt on {input_label} (mode={opt_mode}) → out={opt_dir}")
        _run_cli_main("opt", _opt_cli.cli, args, on_nonzero="raise", on_exception="raise", prefix="endpoint-opt")

        # Read the endpoint opt child's explicit convergence bit from its
        # result.json (fail-closed tri-state) so a nonconverged endpoint cannot
        # silently promote its segment to a usable success.
        endpoint_converged = _read_opt_endpoint_converged(opt_dir)

        final_pdb = opt_dir / "final_geometry.pdb"
        final_xyz = opt_dir / "final_geometry.xyz"
        # Prefer XYZ (full precision) for geometry loading
        if final_xyz.exists():
            final_geom_path = final_xyz
        elif final_pdb.exists():
            final_geom_path = final_pdb
        else:
            raise click.ClickException(f"[endpoint-opt] opt outputs not found under {opt_dir}")

        g_opt = geom_loader(final_geom_path, coord_type="cart")
        calc_input_pdb = final_pdb if final_pdb.exists() else pdb_path
        _opt_calc_kwargs = _stage_calc_kwargs(
            resolved_calc_template,
            input_pdb=calc_input_pdb,
            real_parm7=real_parm7,
            model_pdb=model_pdb,
            charge=q_int,
            spin=spin,
            use_bfactor_layers=detect_layer,
        )
        calc = _mlmm_calc(**_opt_calc_kwargs)
        g_opt.set_calculator(calc)
        _ = float(g_opt.energy)

        return g_opt, final_geom_path, endpoint_converged
    finally:
        prepared_input.cleanup()


def _dft_succeeded(result: Dict[str, Any]) -> bool:
    """Return True only if DFT converged and produced a valid energy."""
    return bool(result) and not result.get("_dft_failed", True)


def _dft_energy_ha(result: Dict[str, Any]) -> Optional[float]:
    """Extract DFT energy in hartree, or None if DFT failed or the value is not finite.

    Non-finite is reported as None at this single chokepoint so that every consumer's
    ``is not None`` check is sufficient; a NaN/inf must never reach a diagram, summary.json
    or logged state energies.
    """
    if not _dft_succeeded(result):
        return None
    return _finite_float((result.get("energy") or {}).get("hartree"))


def _finite_float(value: Any) -> Optional[float]:
    try:
        fval = float(value)
    except (TypeError, ValueError):
        return None
    if not np.isfinite(fval):
        return None
    return fval


def _format_state_values(values: Dict[str, Optional[float]], *, precision: int) -> str:
    parts: List[str] = []
    for label in ("R", "TS", "P"):
        val = values.get(label)
        if val is None:
            parts.append(f"{label}=n/a")
        else:
            parts.append(f"{label}={val:.{precision}f}")
    return " ".join(parts)


def _scale_energy_values(
    values_ha: Dict[str, Optional[float]],
    scale: float,
) -> Dict[str, Optional[float]]:
    return {label: (val * scale if val is not None else None) for label, val in values_ha.items()}


def _relative_energy_values_kcal(values_ha: Dict[str, Optional[float]]) -> Optional[Dict[str, Optional[float]]]:
    ref = values_ha.get("R")
    if ref is None:
        return None
    return {
        label: ((val - ref) * AU2KCALPERMOL if val is not None else None)
        for label, val in values_ha.items()
    }


def _echo_state_energies(
    tag: str,
    seg_idx: int,
    label: str,
    values: Dict[str, Optional[float]],
    *,
    unit: str,
    precision: int,
) -> None:
    if not any(val is not None for val in values.values()):
        return
    _echo_detail(
        f"[{tag}] Segment {seg_idx:02d} {label} ({unit}): "
        f"{_format_state_values(values, precision=precision)}"
    )


def _thermo_correction_values(
    payloads: Dict[str, Dict[str, Any]],
    key: str,
) -> Dict[str, Optional[float]]:
    return {
        label: _finite_float((payloads.get(label) or {}).get(key))
        for label in ("R", "TS", "P")
    }


def _dft_total_mlmm_energy_ha(result: Dict[str, Any]) -> Optional[float]:
    if not _dft_succeeded(result):
        return None
    return _finite_float((result.get("mlmm_energy") or {}).get("E_total_ml_dft_mm_hartree"))


def _run_dft_for_state(pdb_path: Path,
                       q_int: int,
                       spin: int,
                       real_parm7: Path,
                       model_pdb: Path,
                       detect_layer: bool,
                       out_dir: Path,
                       args_yaml: Optional[Path],
                       func_basis: Optional[str] = None,
                       overrides: Optional[Dict[str, Any]] = None,
                       backend: Optional[str] = None,
                       embedcharge: bool = False,
                       embedcharge_cutoff: Optional[float] = None,
                       embedcharge_explicit: bool = False,
                       link_atom_method: Optional[str] = None,
                       mm_backend: Optional[str] = None,
                       use_cmap: Optional[bool] = None,
                       xyz_path: Optional[Path] = None) -> Dict[str, Any]:
    """
    Run dft CLI; return parsed result.yaml dict (may be empty).
    When *xyz_path* is given, use it for full-precision coordinates with
    *pdb_path* as topology reference (--ref-pdb).
    """
    ddir = out_dir
    ensure_dir(ddir)
    overrides = overrides or {}

    func_basis_use = overrides.get("func_basis", func_basis)

    # Prefer XYZ (full precision) with --ref-pdb for topology
    if xyz_path is not None and xyz_path.exists():
        args = ["-i", str(xyz_path), "--ref-pdb", str(pdb_path)]
    else:
        args = ["-i", str(pdb_path)]
    args.extend([
        "--parm", str(real_parm7),
        "--model-pdb", str(model_pdb),
        "-q", str(int(q_int)),
        "-m", str(int(spin)),
        "--out-dir", str(ddir),
    ])
    _append_cli_arg(args, "--func-basis", func_basis_use)
    args.append("--detect-layer" if detect_layer else "--no-detect-layer")

    _append_cli_arg(args, "--max-cycle", overrides.get("max_cycle"))
    _append_cli_arg(args, "--conv-tol", overrides.get("conv_tol"))
    _append_cli_arg(args, "--grid-level", overrides.get("grid_level"))
    _append_cli_arg(args, "--engine", overrides.get("engine"))
    _append_toggle_arg(args, "--convert-files", overrides.get("convert_files"))

    from mlmm.workflows._all_helpers import append_backend_forwarding_args
    append_backend_forwarding_args(
        args,
        backend=backend,
        embedcharge=embedcharge,
        embedcharge_cutoff=embedcharge_cutoff,
        embedcharge_explicit=embedcharge_explicit,
        link_atom_method=link_atom_method,
        mm_backend=mm_backend,
        use_cmap=use_cmap,
        args_yaml=args_yaml,
    )
    # Run DFT as a real subprocess to avoid libcusolver conflict with torch.
    # The MLIP stack (UMA / ORB / MACE via torch) and gpu4pyscf both link
    # against libcusolver but pin different versions; running DFT in the same
    # Python process triggers a dynamic-loader clash. A fresh interpreter is
    # the only reliable isolation.
    # Free GPU memory before spawning the DFT subprocess so it can claim VRAM.
    gc.collect()
    if torch.cuda.is_available():
        torch.cuda.empty_cache()
    import subprocess as _sp
    cmd = [sys.executable, "-m", "mlmm", "dft"] + list(args)
    _echo(f"\n[dft] subprocess: {' '.join(cmd)}")
    proc = _sp.run(cmd, capture_output=True, text=True)
    if proc.stdout:
        _echo(proc.stdout.rstrip())
    if proc.returncode != 0:
        _echo(f"[dft] WARNING: dft exited with code {proc.returncode}", err=True)
        if proc.stderr:
            _echo(proc.stderr.rstrip(), err=True)
    y = out_dir / "result.yaml"
    if y.exists():
        try:
            data = yaml.safe_load(y.read_text(encoding="utf-8")) or {}
        except Exception as exc:
            logger.debug("Failed to parse DFT result YAML %s: %s", y, exc)
            data = {}
    else:
        data = {}
    converged = (data.get("energy") or {}).get("converged", False)
    data["_dft_converged"] = bool(converged)
    data["_dft_failed"] = not bool(converged) or proc.returncode != 0
    return data



_ALL_PRIMARY_HELP_OPTIONS = frozenset(
    {
        "-i",
        "--input",
        "-c",
        "--center",
        "-l",
        "--ligand-charge",
        "-q",
        "--charge",
        "--out-dir",
        "--tsopt",
        "--thermo",
        "--dft",
        "--dft-func-basis",
        "--config",
        "--dry-run",
        "--embedcharge",
        "-s",
        "--scan-lists",
        "-b",
        "--backend",
        "--refine-path",
        "-o",
        "--help-advanced",
    }
)


def _configure_all_help_visibility(command: click.Command) -> None:
    """Hide advanced options from default --help while keeping them functional.

    Routes through the single ``help_pages`` visibility implementation so
    ``all`` and the lazily-loaded subcommands share one callback + one loop.
    """
    _hide_advanced_options(command, _ALL_PRIMARY_HELP_OPTIONS)


@click.command(
    help="Run pocket extraction → (optional single-structure staged scan) → MEP search in one shot.\n"
         "If exactly one input is provided: (a) with --scan-lists, stage results feed into path-opt (or path_search with --refine-path); "
         "(b) with --tsopt and no --scan-lists, run TSOPT-only mode.",
    context_settings={
        "help_option_names": ["-h", "--help"],
        "ignore_unknown_options": True,
        "allow_extra_args": True,
    },
)
@click.option(
    "--help-advanced",
    is_flag=True,
    is_eager=True,
    expose_value=False,
    callback=_show_advanced_subcommand_help,
    help="Show all options (including advanced settings) and exit.",
)
# ===== Inputs =====
@click.option(
    "-i", "--input", "input_paths",
    type=click.Path(path_type=Path, exists=True, dir_okay=False),
    multiple=True, required=True,
    help=("Two or more full PDB/mmCIF structures in reaction order "
          "(reactant [intermediates ...] product), or one full structure with "
          "--scan-lists or --tsopt. A single '-i' may be followed by multiple "
          "space-separated files (for example, '-i A.pdb B.pdb C.pdb').")
)
@click.option(
    "-c", "--center", "center_spec",
    type=str, required=False, default=None,
    help=("Substrate specification for the extractor: "
          "a PDB path, a residue-ID list like '123,124' or 'A:123,B:456' "
          "(insertion codes OK: '123A' / 'A:123A'), "
          "or a residue-name list like 'GPP,MMT'. "
          "When omitted, extraction is skipped and full structures are used directly.")
)
@click.option(
    "-o", "--out-dir", "out_dir",
    type=click.Path(path_type=Path, file_okay=False),
    default=Path(OUT_DIR_ALL), show_default=True,
    help="Top-level output directory for the pipeline."
)
# ===== Extractor knobs (subset of extract.parse_args) =====
@click.option("-r", "--radius", type=float, default=2.6, show_default=True,
              help="Inclusion cutoff (Å) around substrate atoms.")
@click.option("--radius-het2het", type=float, default=0.0, show_default=True,
              help="Independent hetero–hetero cutoff (Å) for non‑C/H pairs.")
@click.option("--include-h2o/--no-include-h2o", "include_h2o", default=True, show_default=True,
              help="Include waters (HOH/WAT/H2O/DOD/TIP/TIP3/SOL) in the pocket.")
@click.option("--exclude-backbone/--no-exclude-backbone", "exclude_backbone", default=False, show_default=True,
              help="Remove backbone atoms on non‑substrate amino acids (with PRO/HYP safeguards).")
@click.option("--add-linkh/--no-add-linkh", "add_linkh", default=False, show_default=True,
              help=("Add extractor-only link H to scratch pocket PDBs. The ML/MM "
                    "model selection remains link-free; runtime link H are generated "
                    "from parm7 boundary bonds."))
@click.option("--selected-resn", type=str, default="", show_default=True,
              help="Force-include residues (comma/space separated; chain/insertion codes allowed).")
@click.option("--modified-residue", type=str, default="", show_default=True,
              help=("Comma-separated residue names (with optional charge) to treat as amino acids "
                    "for backbone truncation and charge assignment. "
                    "Examples: 'HD1,HD2,HD3' or 'HD1:0,SEP:-2'."))
@click.option("-l", "--ligand-charge", type=str, default=None,
              help=("Either a total charge (number) to distribute across unknown residues "
                    "or a mapping like 'GPP:-3,MMT:-1'."))
@click.option(
    "-q",
    "--charge",
    "charge_override",
    type=int,
    default=None,
    help=(
        "Override the net charge of the ML region/model atoms. "
        "Highest priority over the charge derived by the workflow."
    ),
)
@click.option(
    "--parm",
    "parm7_override",
    type=click.Path(path_type=Path, exists=True, dir_okay=False),
    default=None,
    help="Pre-built AMBER parm7 topology file. When provided, mm_parm generation is skipped.",
)
@click.option(
    "--model-pdb",
    "model_pdb_override",
    type=click.Path(path_type=Path, exists=True, dir_okay=False),
    default=None,
    help=("ML-only atom-selection PDB. It must be an unchanged, link-H-free subset "
          "of the full PDB/parm7 in the same atom order. It takes precedence "
          "over ML membership from -c/--center or input B-factors."),
)
@click.option("--auto-mm-ff-set", "mm_ff_set",
              type=click.Choice(["ff19SB", "ff14SB"], case_sensitive=False),
              default="ff19SB", show_default=True,
              help="Force-field set forwarded to mm_parm (ff19SB uses OPC3; ff14SB uses TIP3P).")
@click.option("--auto-mm-add-ter/--auto-mm-no-add-ter", "mm_add_ter",
              default=True, show_default=True,
              help="Control mm_parm TER insertion around ligand/water/ion blocks and disconnected peptide blocks.")
@click.option("--auto-mm-keep-temp", "mm_keep_temp", is_flag=True, default=False, show_default=True,
              help="Keep the mm_parm temporary working directory (for debugging).")
@click.option(
    "--auto-mm-ligand-mult",
    "mm_ligand_mult",
    type=str,
    default=None,
    help=("Spin multiplicity mapping forwarded to mm_parm (e.g., 'GPP:2,SAM:1'). "
          "If omitted, mm_parm defaults to 1 for all ligands.")
)
@click.option("--auto-mm-disulfide/--auto-mm-no-disulfide", "mm_auto_disulfide",
              default=True, show_default=True,
              help="Forwarded to mm_parm: detect disulfides from SG-SG geometry across "
                   "CYS/CYM/CYX and bond them (renaming a bonded CYS to CYX). With "
                   "--auto-mm-no-disulfide only residues already named CYX are bonded "
                   "and CYS is left untouched.")
# ===== Path search knobs (subset of path_search.cli) =====
@click.option("-m", "--multiplicity", "spin", type=int, default=1, show_default=True, help="Multiplicity (2S+1).")
@click.option(
    "--tr-projection",
    type=click.Choice(["constrained", "legacy-active"], case_sensitive=False),
    default=GEOM_KW_DEFAULT["tr_projection"],
    show_default=True,
    help=(
        "Rigid translation/rotation treatment forwarded to TSopt, IRC, freq, "
        "and flatten PHVA. The default respects frozen anchors; 'legacy-active' "
        "is deprecated and must not be used for pass/HOSP transition-state certification."
    ),
)
@click.option(
    "--mep-mode",
    type=click.Choice(["gsm", "dmf"], case_sensitive=False),
    default="gsm",
    show_default=True,
    help="MEP optimizer: Growing String Method (gsm) or Direct Max Flux (dmf).",
)
@click.option(
    "--dmf-backend",
    type=click.Choice(["cpu", "gpu"], case_sensitive=False),
    default="gpu",
    show_default=True,
    help=(
        "DMF compute backend (--mep-mode dmf only): gpu (dmf.torch / CUDA) "
        "or cpu (dmf / NumPy). On a GPU out-of-memory error, retry with cpu."
    ),
)
@click.option("--max-nodes", type=int, default=_path_opt.GS_KW["max_nodes"], show_default=True,
              help="Max internal nodes per GSM/DMF segment (max_nodes+2 images including endpoints).")
@click.option("--max-cycles", type=int, default=300, show_default=True, help="Maximum MEP optimization cycles.")
@click.option("--climb/--no-climb", default=True, show_default=True,
              help="Enable transition-state climbing after growth for the *first* segment in each pair.")
@click.option(
    "--opt-mode",
    type=click.Choice(["grad", "hess"], case_sensitive=False),
    default="grad",
    show_default=True,
    help=(
        "Optimizer mode forwarded to scan/path-search and used for single optimizations: "
        "grad (=L-BFGS/Dimer) or hess (=RFO/RSIRFO)."
    ),
)
@click.option(
    "--opt-mode-post",
    type=click.Choice(["grad", "hess"], case_sensitive=False),
    default="hess",
    show_default=True,
    help=(
        "Optimizer mode for TSOPT and post-IRC endpoint optimizations. "
        "Takes precedence over --opt-mode for these stages."
    ),
)
@click.option("--dump/--no-dump", default=False, show_default=True,
              help="Dump MEP / single-structure trajectories during the run, forwarding the same flag to scan/tsopt/freq.")
@click.option(
    "--refine-path/--no-refine-path",
    "refine_path",
    default=False,
    show_default=True,
    help=(
        "If False (default), run single-pass path-opt with the selected MEP optimizer between each adjacent pair and concatenate the "
        "segments (no path_search); if True, run recursive path_search on the full ordered series for "
        "automatic multistep discovery."
    ),
)
@click.option(
    "--thresh",
    type=click.Choice(THRESH_CHOICES, case_sensitive=False),
    default=None,
    show_default=False,
    help=(
        "Convergence preset (gau_loose|gau|gau_tight|gau_vtight|baker|never). "
        "Defaults to 'gau_loose' for path-opt, 'gau' for scan."
    ),
)
@click.option(
    "--thresh-post",
    type=click.Choice(THRESH_CHOICES, case_sensitive=False),
    default="baker",
    show_default=True,
    help=(
        "Convergence preset for post-IRC endpoint optimizations "
        "(gau_loose|gau|gau_tight|gau_vtight|baker|never)."
    ),
)
@click.option("--config", "config_yaml", type=click.Path(path_type=Path, exists=True, dir_okay=False),
              default=None, help="Base YAML configuration file applied before explicit CLI options.")
@click.option("--show-config/--no-show-config", "show_config", default=False, show_default=True,
              help="Print resolved configuration and continue execution.")
@click.option("--dry-run/--no-dry-run", "dry_run", default=False, show_default=True,
              help="Run input preparation and preflight checks in a temporary directory, "
                   "print the execution plan, and skip calculation stages.")
@click.option("--preopt/--no-preopt", "pre_opt", default=True, show_default=True,
              help="Run initial single-structure optimizations of the pocket inputs.")
@click.option("--hessian-calc-mode",
              type=click.Choice(["Analytical", "FiniteDifference"], case_sensitive=False),
              default=None,
              help="Common MLIP Hessian calculation mode forwarded to tsopt and freq. Default: 'FiniteDifference'. Use 'Analytical' when VRAM is sufficient.")
@click.option(
    "--detect-layer/--no-detect-layer",
    "detect_layer",
    default=True,
    show_default=True,
    help="Detect ML/MM layers from input PDB B-factors (ML=0, MovableMM=10, FrozenMM=20) in downstream tools. "
         "If disabled, mlmm all requires --model-pdb.",
)
# ===== Post-processing toggles =====
@click.option("--tsopt/--no-tsopt", "do_tsopt", default=False, show_default=True,
              help="TS optimization + EulerPC IRC per reactive segment (or TSOPT-only mode for single-structure), and build energy diagrams.")
@click.option("--thermo/--no-thermo", "do_thermo", default=False, show_default=True,
              help="Run freq on (R,TS,P) per reactive segment (or TSOPT-only mode) and build Gibbs free-energy diagram (MLIP).")
@click.option("--dft/--no-dft", "do_dft", default=False, show_default=True,
              help="Run DFT single-point on (R,TS,P) and build a DFT energy diagram. With --thermo, also generate a DFT//MLIP/MM Gibbs diagram.")
@click.option("--tsopt-max-cycles", type=int, default=None,
              help="Override tsopt --max-cycles value.")
@click.option(
    "--flatten/--no-flatten",
    "flatten",
    default=False,
    show_default=True,
    help="Enable the extra-imaginary-mode flattening loop in tsopt (grad: dimer loop, hess: post-RSIRFO); --no-flatten forces flatten_max_iter=0.",
)
@click.option(
    "--reject-uphill/--no-reject-uphill",
    "reject_uphill",
    default=True,
    show_default=True,
    help=(
        "Reject uphill RFO trials during post-IRC endpoint re-optimization only "
        "and final-check the retained endpoint at the emergency floor. Does not "
        "affect TS optimization or path search."
    ),
)
@click.option(
    "--irc-step-size",
    type=float,
    default=None,
    help=(
        "Override IRC --step-size (Bohr). If an IRC stops after only a few "
        "frames, retry with a smaller value such as 0.05."
    ),
)
@click.option(
    "--irc-never-stop/--no-irc-never-stop",
    "irc_never_stop",
    default=None,
    help=(
        "Forward IRC never-stop mode to every post-TS IRC. It ignores "
        "energy-rise/plateau stops but retains physical/integrator stops; "
        "default follows irc.never_stop (off)."
    ),
)
@click.option(
    "--skip-final-freq/--no-skip-final-freq",
    "skip_final_freq",
    default=False,
    show_default=True,
    help="Skip post-convergence frequency analysis in tsopt. Useful for large unfrozen systems.",
)
@click.option("--tsopt-out-dir", type=click.Path(path_type=Path, file_okay=False), default=None,
              help="Override tsopt output subdirectory (relative paths are resolved against the default).")
@click.option("--freq-out-dir", type=click.Path(path_type=Path, file_okay=False), default=None,
              help="Override freq output base directory (relative paths resolved against the default).")
@click.option("--freq-max-write", type=int, default=None,
              help="Override freq --max-write value.")
@click.option("--freq-amplitude-ang", type=float, default=None,
              help="Override freq --amplitude-ang (Å).")
@click.option("--freq-n-frames", type=int, default=None,
              help="Override freq --n-frames value.")
@click.option("--freq-sort", type=click.Choice(["value", "abs"], case_sensitive=False), default=None,
              help="Override freq mode sorting.")
@click.option("--freq-temperature", type=float, default=None,
              help="Override freq thermochemistry temperature (K).")
@click.option("--freq-pressure", type=float, default=None,
              help="Override freq thermochemistry pressure (atm).")
@click.option(
    "--freq-symmetry-number",
    type=click.IntRange(min=1),
    default=None,
    help=(
        "Use one rotational symmetry number for every R/TS/P frequency job. "
        "When omitted, each child follows its YAML/default setting."
    ),
)
@click.option("--dft-out-dir", type=click.Path(path_type=Path, file_okay=False), default=None,
              help="Override dft output base directory (relative paths resolved against the default).")
@click.option("--dft-func-basis", type=str, default=None,
              help="Override dft --func-basis value.")
@click.option("--dft-max-cycle", type=int, default=None,
              help="Override dft --max-cycle value.")
@click.option("--dft-conv-tol", type=float, default=None,
              help="Override dft --conv-tol value.")
@click.option("--dft-grid-level", type=int, default=None,
              help="Override dft --grid-level value.")
@click.option("--dft-engine", type=click.Choice(["gpu", "cpu"]), default=None,
              help="Override dft --engine value.")
# ===== Staged scan specification for single-structure route =====
@click.option(
    "-s", "--scan-lists",
    "scan_lists_raw",
    type=str, multiple=True, required=False,
    help='Scan targets: inline Python literal or a YAML/JSON spec file path. '
         'Multiple inline literals define sequential stages, e.g. '
         '"[(12,45,1.35)]" "[(10,55,2.20),(23,34,1.80)]". '
         'Indices refer to the original full PDB (1-based) or PDB atom selectors like "TYR,285,CA"; '
         'they are auto-mapped to the pocket after extraction.',
)
@click.option("--scan-out-dir", type=click.Path(path_type=Path, file_okay=False), default=None,
              help="Override the scan output directory (default: <out-dir>/scan/). Relative paths are resolved against the default parent.")
@click.option("--scan-one-based/--scan-zero-based", default=None,
              help="Override scan indexing interpretation (one-based or zero-based).")
@click.option("--scan-max-step-size", type=float, default=None,
              help="Override scan --max-step-size (Å).")
@click.option("--scan-bias-k", type=float, default=None,
              help="Override scan harmonic bias strength k (eV/Å^2).")
@click.option("--scan-relax-max-cycles", type=int, default=None,
              help="Override scan relaxation max cycles per step.")
@click.option("--scan-preopt/--no-scan-preopt", "scan_preopt_override", default=None,
              help="Override scan --preopt flag.")
@click.option("--scan-endopt/--no-scan-endopt", "scan_endopt_override", default=None,
              help="Override scan --endopt flag.")
@click.option("--convert-files/--no-convert-files", "convert_files", default=True, show_default=True,
              help="Convert XYZ/TRJ outputs to PDB format using reference topology; forwarded to all subcommands.")
@click.option(
    "--ref-pdb",
    "ref_pdb_cli",
    type=click.Path(path_type=Path, exists=True, dir_okay=False),
    default=None,
    help=(
        "Reference PDB for topology/B-factor layer information when -i provides XYZ inputs. "
        "Used for define-layer, mm_parm, ml_region, and forwarded to downstream tools "
        "(tsopt, irc, freq, path_search) as --ref-pdb."
    ),
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
    help="Unavailable in v0.3.3; retained so older commands fail with an actionable diagnostic.",
)
@click.option(
    "--embedcharge-cutoff",
    "embedcharge_cutoff",
    type=float,
    default=None,
    show_default=False,
    help="Unavailable in v0.3.3 together with the retired electronic-embedding path.",
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
    help="MM backend (default: hessian_ff). MM Hessians use finite differences by default; set calc.mm_fd: false for the hessian_ff analytical path.",
)
@click.option(
    "--cmap/--no-cmap",
    "use_cmap",
    default=None,
    show_default=False,
    help="Preserve CMAP terms in both real and model MM layers. Default: enabled when present in parm7.",
)
@add_coord_type_option(choices=("cart", "dlc"))
@add_precision_option()
@add_workers_options()
@add_backend_model_option()
@add_calc_file_option()
@add_deterministic_option()
@add_allow_charge_mult_mismatch_option()
@click.pass_context
def cli(
    ctx: click.Context,
    input_paths: Sequence[Path],
    center_spec: Optional[str],
    out_dir: Path,
    radius: float,
    radius_het2het: float,
    include_h2o: bool,
    exclude_backbone: bool,
    add_linkh: bool,
    selected_resn: str,
    modified_residue: str,
    ligand_charge: Optional[str],
    charge_override: Optional[int],
    parm7_override: Optional[Path],
    model_pdb_override: Optional[Path],
    mm_ff_set: str,
    mm_add_ter: bool,
    mm_auto_disulfide: bool,
    mm_keep_temp: bool,
    mm_ligand_mult: Optional[str],
    spin: int,
    tr_projection: str,
    mep_mode: str,
    dmf_backend: str,
    max_nodes: int,
    max_cycles: int,
    climb: bool,
    opt_mode: str,
    opt_mode_post: Optional[str],
    dump: bool,
    refine_path: bool,
    thresh: Optional[str],
    thresh_post: str,
    config_yaml: Optional[Path],
    show_config: bool,
    dry_run: bool,
    pre_opt: bool,
    hessian_calc_mode: Optional[str],
    detect_layer: bool,
    do_tsopt: bool,
    do_thermo: bool,
    do_dft: bool,
    scan_lists_raw: Sequence[str],
    scan_out_dir: Optional[Path],
    scan_one_based: Optional[bool],
    scan_max_step_size: Optional[float],
    scan_bias_k: Optional[float],
    scan_relax_max_cycles: Optional[int],
    scan_preopt_override: Optional[bool],
    scan_endopt_override: Optional[bool],
    convert_files: bool,
    ref_pdb_cli: Optional[Path],
    backend: Optional[str],
    embedcharge: bool,
    embedcharge_cutoff: Optional[float],
    link_atom_method: Optional[str],
    mm_backend: Optional[str],
    use_cmap: Optional[bool],
    tsopt_max_cycles: Optional[int],
    flatten: bool,
    reject_uphill: bool,
    irc_step_size: Optional[float],
    irc_never_stop: Optional[bool],
    skip_final_freq: bool,
    tsopt_out_dir: Optional[Path],
    freq_out_dir: Optional[Path],
    freq_max_write: Optional[int],
    freq_amplitude_ang: Optional[float],
    freq_n_frames: Optional[int],
    freq_sort: Optional[str],
    freq_temperature: Optional[float],
    freq_pressure: Optional[float],
    freq_symmetry_number: Optional[int],
    dft_out_dir: Optional[Path],
    dft_func_basis: Optional[str],
    dft_max_cycle: Optional[int],
    dft_conv_tol: Optional[float],
    dft_grid_level: Optional[int],
    dft_engine: Optional[str],
    cli_coord_type: Optional[str],
    precision: Optional[str],
    workers: Optional[int],
    workers_per_node: Optional[int],
    backend_model: Optional[str],
    calc_file: Optional[str],
    calc_factory: Optional[str],
) -> None:
    """
    The **all** command composes `extract` → (optional `scan` on pocket) → MEP search (single-pass `path-opt` by default,
    or recursive `path_search` with ``--refine-path``) and hides ref-template bookkeeping.
    It also accepts the sloppy `-i A B C` style like `path_search` does. With single input:
      - with --scan-lists: run staged scan on the pocket and use stage results as inputs for path-opt (or path_search),
      - with --tsopt and no --scan-lists: run TSOPT-only mode (no MEP search).
    """
    from mlmm.core.utils import (
        collect_option_values,
        current_cli_args,
        reject_option_like_extra_args,
    )

    _argv = current_cli_args(ctx)

    reject_option_like_extra_args(
        ctx.args,
        allowed_values=collect_option_values(
            _argv, ("-i", "--input", "-s", "--scan-lists")
        ),
        consumed_values=[*input_paths, *scan_lists_raw],
    )
    # Turn on pipeline-scoped default-verbosity suppression for this `all` run
    # (reset per-invocation in DefaultGroup.parse_args). Standalone leaf/report
    # commands are unaffected and keep full output at default verbosity.
    from mlmm.core import utils as _mlmm_utils
    from mlmm.core.utils import set_pipeline_mode

    # Register one invocation owner before the first process-global mutation so
    # every exit path (success, failure, Ctrl-C) restores the exact prior state
    # and releases owned resources in LIFO order via ctx.call_on_close.
    session = RunSession()
    ctx.call_on_close(session.close)
    session.own_run_id_environment()
    ctx.meta["mlmm_run_session"] = session
    manifest = session.manifest
    prior_pipeline_mode = bool(_mlmm_utils._PIPELINE_MODE)
    prior_echo_started = bool(_echo_state._started)

    def _restore_invocation_state() -> None:
        _echo_state._started = prior_echo_started
        set_pipeline_mode(prior_pipeline_mode)

    session.resources.add(_restore_invocation_state)

    set_pipeline_mode(True)
    _echo_state.reset()

    time_start = time.perf_counter()
    command_str = "mlmm " + " ".join(_argv)

    _is_param_explicit = make_is_param_explicit(ctx)
    # Post-IRC endpoint re-optimization uphill-rejection toggle, forwarded to the
    # opt child. ``None`` unless the flag was explicitly passed, so the default
    # path keeps the opt child's own RFO_KW reject_uphill (unchanged behavior).
    _reject_uphill_eff = (
        bool(reject_uphill) if _is_param_explicit("reject_uphill") else None
    )
    explicit_params = frozenset(
        parameter.name
        for parameter in ctx.command.params
        if parameter.name and _is_param_explicit(parameter.name)
    )
    dump_override_requested = _is_param_explicit("dump")
    opt_mode_set = _is_param_explicit("opt_mode")
    opt_mode_post_set = _is_param_explicit("opt_mode_post")
    # Gate --no-embedcharge forwarding: when CLI default False is not user-supplied,
    # do not emit --no-embedcharge to downstream subprocesses (otherwise a `calc.embedcharge: true`
    # in --config YAML is silently overridden by the CLI default).
    embedcharge_explicit = _is_param_explicit("embedcharge")
    from mlmm.core.embedcharge_policy import reject_retired_embedcharge_cli

    # Reject explicit activation before resolving inputs or calculator state.
    # A later check covers activation inherited from YAML.
    if (embedcharge_explicit and embedcharge) or _is_param_explicit(
        "embedcharge_cutoff"
    ):
        reject_retired_embedcharge_cli(
            {"embedcharge": embedcharge},
            cutoff_requested=_is_param_explicit("embedcharge_cutoff"),
        )

    config_yaml, override_yaml, _ = resolve_yaml_sources(config_yaml, None, None)
    args_yaml, merged_yaml_cfg = _build_effective_args_yaml(
        config_yaml=config_yaml,
        override_yaml=None,
        tmp_prefix="mlmm_all_merged_",
    )
    yaml_model_charge, yaml_model_spin = configured_model_charge_spin(
        merged_yaml_cfg
    )
    if not _is_param_explicit("spin") and yaml_model_spin is not None:
        spin = int(yaml_model_spin)
    _dmf_yaml_cfg = (
        merged_yaml_cfg.get("dmf", {})
        if isinstance(merged_yaml_cfg, dict)
        else {}
    )
    dmf_backend_effective = str(
        dmf_backend
        if "dmf_backend" in explicit_params
        else (
            _dmf_yaml_cfg.get("backend", dmf_backend)
            if isinstance(_dmf_yaml_cfg, dict)
            else dmf_backend
        )
    ).lower()
    dmf_correlated_effective = bool(
        fresh_dmf_config(_dmf_yaml_cfg).get("correlated", False)
    )
    _injected_coord = (
        str(cli_coord_type).lower()
        if _is_param_explicit("cli_coord_type") and cli_coord_type is not None
        else None
    )
    _injected_tr_projection = (
        str(tr_projection).lower()
        if _is_param_explicit("tr_projection")
        else None
    )
    if (
        args_yaml is not None
        or _injected_coord is not None
        or _injected_tr_projection is not None
        or precision is not None
        or workers is not None
        or workers_per_node is not None
        or backend_model is not None
        or calc_file is not None
    ):
        args_yaml = _inject_coord_type_into_args_yaml(
            args_yaml, _injected_coord, tr_projection=_injected_tr_projection,
            backend=backend,
            precision=precision, workers=workers, workers_per_node=workers_per_node, backend_model=backend_model,
            calc_file=(str(Path(calc_file).resolve()) if calc_file else None), calc_factory=calc_factory,
        )
    (
        mlip_backend_resolved,
        mlip_model_resolved,
        mlip_precision_resolved,
    ) = _resolve_mlip_provenance(
        backend=backend,
        backend_model=backend_model,
        calc_file=calc_file,
        calc_factory=calc_factory,
        precision=precision,
        merged_yaml_cfg=merged_yaml_cfg,
    )
    resolved_calc_template = _resolve_calculator_template(
        args_yaml,
        backend=backend,
        embedcharge=embedcharge,
        embedcharge_explicit=embedcharge_explicit,
        embedcharge_cutoff=embedcharge_cutoff,
        link_atom_method=link_atom_method,
        mm_backend=mm_backend,
        use_cmap=use_cmap,
    )
    reject_retired_embedcharge_cli(
        resolved_calc_template.materialize(),
        cutoff_requested=_is_param_explicit("embedcharge_cutoff"),
    )

    mm_ff_set = "ff14SB" if str(mm_ff_set).lower().startswith("ff14") else "ff19SB"

    # --- Robustly accept a single "-i" followed by multiple paths (like path_search.cli) ---
    argv_all = _argv
    i_vals = collect_option_values(argv_all, ("-i", "--input"))
    if i_vals:
        i_parsed = validate_existing_files(
            i_vals,
            option_name="-i/--input",
            hint="When using '-i', list only existing file paths (multiple paths may follow a single '-i').",
        )
        input_paths = tuple(i_parsed)

    scan_vals = collect_option_values(argv_all, ("-s", "--scan-lists"))
    if scan_vals:
        scan_lists_raw = tuple(scan_vals)

    is_single = (len(input_paths) == 1)
    has_scan = bool(scan_lists_raw)
    if has_scan and not is_single:
        raise click.BadParameter(
            "--scan-lists requires exactly one input structure; provide one "
            "reactant for a staged scan or omit --scan-lists for an endpoint path."
        )
    single_tsopt_mode = (is_single and (not has_scan) and do_tsopt)

    if (len(input_paths) < 2) and (not (is_single and (has_scan or do_tsopt))):
        raise click.BadParameter(
            "Provide at least two structures with -i/--input in reaction order, "
            "or use one structure with --scan-lists or --tsopt."
        )

    # Normalize mmCIF and oversized PDB inputs once and keep the bridge alive
    # for the complete in-process pipeline. Computational stages continue to
    # consume safe internal PDB files; public outputs regain original IDs via
    # the registered coordinate template.
    _prepared_all_inputs: List[Any] = []
    _original_input_paths = tuple(Path(path) for path in input_paths)
    for path in _original_input_paths:
        suffix = path.suffix.lower()
        if suffix not in {".pdb", ".cif", ".mmcif", ".xyz"}:
            raise click.BadParameter(
                f"Unsupported input format '{suffix}' for {path.name}; use "
                ".pdb/.cif/.mmcif or .xyz with --ref-pdb."
            )
        prepared = prepare_input_structure(path)
        if suffix == ".xyz":
            if ref_pdb_cli is None:
                prepared.cleanup()
                raise click.BadParameter("XYZ input requires --ref-pdb topology.")
            apply_ref_pdb_override(prepared, ref_pdb_cli)
        _prepared_all_inputs.append(prepared)
        # Own the temp-tree cleanup on the session so it runs on the real
        # success/exception path too (not only the dry-run/validation branches).
        session.resources.own_cleanup(prepared)
    input_paths = tuple(
        prepared.source_path for prepared in _prepared_all_inputs
    )

    _prepared_ref_pdb = None
    if ref_pdb_cli is not None:
        _prepared_ref_pdb = prepare_input_structure(ref_pdb_cli)
        session.resources.own_cleanup(_prepared_ref_pdb)
        ref_pdb_cli = _prepared_ref_pdb.source_path

    _prepared_model_pdb = None
    if model_pdb_override is not None:
        _prepared_model_pdb = prepare_input_structure(model_pdb_override)
        session.resources.own_cleanup(_prepared_model_pdb)
        model_pdb_override = _prepared_model_pdb.source_path

    if single_tsopt_mode:
        all_mode = "tsopt-only"
    elif has_scan:
        all_mode = "scan-to-path-search" if refine_path else "scan-to-path-opt"
    else:
        all_mode = "path-search" if refine_path else "path-opt"
    all_mode_label = "ts-only" if single_tsopt_mode else ("scan-lists" if has_scan else "mep")
    if verbose_level() >= 2:
        _echo(
            f"[mode] all ({all_mode_label}) inputs={len(input_paths)} "
            f"extract={'yes' if center_spec is not None and str(center_spec).strip() else 'no'} "
            f"scan={'yes' if has_scan else 'no'} tsopt={'yes' if do_tsopt else 'no'} "
            f"thermo={'yes' if do_thermo else 'no'} dft={'yes' if do_dft else 'no'} "
            f"dry_run={'yes' if dry_run else 'no'} internal={all_mode}",
            narrative=True,
        )

    _validate_postprocessing_dependencies(
        do_tsopt=do_tsopt,
        do_thermo=do_thermo,
        do_dft=do_dft,
    )

    _mode_alias = {
        "grad": "grad",
        "hess": "hess",
        "light": "grad",
        "heavy": "hess",
    }
    opt_mode_norm = _mode_alias.get(str(opt_mode).strip().lower(), "grad")
    mep_mode_kind = str(mep_mode).strip().lower()
    path_search_opt_mode = opt_mode_norm
    opt_mode_post_norm = (
        None
        if opt_mode_post is None
        else _mode_alias.get(str(opt_mode_post).strip().lower(), "hess")
    )
    endpoint_opt_mode_default = (
        opt_mode_post_norm if (opt_mode_post_set and opt_mode_post_norm is not None)
        else (opt_mode_norm if opt_mode_set else "hess")
    )
    if opt_mode_post_norm in {"grad", "hess"}:
        tsopt_opt_mode_default = opt_mode_post_norm
    elif opt_mode_set:
        tsopt_opt_mode_default = opt_mode_norm
    else:
        tsopt_opt_mode_default = "hess"

    citation_post_segments: List[Dict[str, Any]] = []

    def _all_method_citation_payload() -> Dict[str, Any]:
        return {
            "pipeline_mode": all_mode,
            "path_opt_mode": path_search_opt_mode,
            "post_opt_mode": tsopt_opt_mode_default,
            "ts_opt_mode": tsopt_opt_mode_default,
            "endpoint_opt_mode": endpoint_opt_mode_default,
            "mep_mode": mep_mode_kind,
            "dmf_correlated": dmf_correlated_effective,
            "post_segments": citation_post_segments,
        }

    from mlmm.workflows._all_helpers import (
        build_path_child_argv as _build_path_child_argv,
        build_scan_child_argv as _build_scan_child_argv,
        build_tsopt_overrides as _build_tsopt_overrides,
        build_freq_overrides as _build_freq_overrides,
        build_dft_overrides as _build_dft_overrides,
        resolve_dft_func_basis_forwarding as _resolve_dft_func_basis_forwarding,
        resolve_post_thresh_forwarding as _resolve_post_thresh_forwarding,
    )
    post_thresh_forward = _resolve_post_thresh_forwarding(
        explicit_params,
        thresh_post=thresh_post,
        yaml_cfg=merged_yaml_cfg,
    )
    (
        dft_func_basis_use,
        dft_method_fallback,
    ) = _resolve_dft_func_basis_forwarding(
        explicit_params,
        dft_func_basis=dft_func_basis,
        yaml_cfg=merged_yaml_cfg,
    )
    tsopt_overrides = _build_tsopt_overrides(
        tsopt_max_cycles=tsopt_max_cycles,
        dump=dump,
        dump_override_requested=dump_override_requested,
        tsopt_out_dir=tsopt_out_dir,
        hessian_calc_mode=hessian_calc_mode,
        opt_mode_post_norm=opt_mode_post_norm,
        opt_mode_post_set=opt_mode_post_set,
        opt_mode_set=opt_mode_set,
        tsopt_opt_mode_default=tsopt_opt_mode_default,
        convert_files=convert_files,
        convert_files_explicit=("convert_files" in explicit_params),
        thresh_post_forward=post_thresh_forward,
        flatten_explicit=_is_param_explicit("flatten"),
        flatten=flatten,
        skip_final_freq=skip_final_freq,
        skip_final_freq_explicit=("skip_final_freq" in explicit_params),
    )
    from mlmm.workflows.freq import _validated_thermo_condition

    if freq_temperature is not None:
        freq_temperature = _validated_thermo_condition(
            freq_temperature, name="temperature"
        )
    if freq_pressure is not None:
        freq_pressure = _validated_thermo_condition(
            freq_pressure, name="pressure_atm"
        )
    freq_overrides = _build_freq_overrides(
        freq_max_write=freq_max_write,
        freq_amplitude_ang=freq_amplitude_ang,
        freq_n_frames=freq_n_frames,
        freq_sort=freq_sort,
        freq_temperature=freq_temperature,
        freq_pressure=freq_pressure,
        freq_symmetry_number=freq_symmetry_number,
        dump_override_requested=dump_override_requested,
        dump=dump,
        require_thermo_artifact=do_thermo,
        hessian_calc_mode=hessian_calc_mode,
        convert_files=convert_files,
        convert_files_explicit=("convert_files" in explicit_params),
    )
    dft_overrides = _build_dft_overrides(
        dft_max_cycle=dft_max_cycle,
        dft_conv_tol=dft_conv_tol,
        dft_grid_level=dft_grid_level,
        dft_engine=dft_engine,
        dft_func_basis_forward=dft_func_basis_use,
        convert_files=convert_files,
        convert_files_explicit=("convert_files" in explicit_params),
    )

    post_convert_files_forward = (
        convert_files if "convert_files" in explicit_params else None
    )

    if show_config or (dry_run and verbose_level() >= 3):
        config_payload: Dict[str, Any] = {
            "yaml": {
                "config": str(config_yaml) if config_yaml else None,
                "override_yaml": str(override_yaml) if override_yaml else None,
                "effective_args_yaml": str(args_yaml) if args_yaml else None,
            },
            "all": {
                "inputs": [str(p) for p in input_paths],
                "center": center_spec,
                "charge_override": charge_override,
                "skip_extract": bool(center_spec is None or str(center_spec).strip() == ""),
                "out_dir": str(out_dir),
                "spin": int(spin),
                "mep_mode": mep_mode_kind,
                "dmf_backend": dmf_backend_effective,
                "max_nodes": int(max_nodes),
                "max_cycles": int(max_cycles),
                "climb": bool(climb),
                "opt_mode": str(opt_mode),
                "opt_mode_post": (None if opt_mode_post is None else str(opt_mode_post)),
                "path_search_opt_mode": str(path_search_opt_mode),
                "endpoint_opt_mode": str(endpoint_opt_mode_default),
                "dump": bool(dump),
                "refine_path": bool(refine_path),
                "thresh": thresh,
                "thresh_post": thresh_post,
                "flatten": bool(flatten),
                "pre_opt": bool(pre_opt),
                "detect_layer": bool(detect_layer),
                "tsopt": bool(do_tsopt),
                "thermo": bool(do_thermo),
                "dft": bool(do_dft),
            },
            "overrides": {
                "tsopt": tsopt_overrides,
                "freq": freq_overrides,
                "dft": dft_overrides,
            },
        }
        if merged_yaml_cfg:
            config_payload["effective_yaml"] = merged_yaml_cfg
        _echo_section("====== [all] Effective configuration ======")
        # `--show-config` is an explicit output request; dry-run's automatic
        # config dump is level-3 debug context so -v 1/2 stay compact.
        emit(
            yaml.safe_dump(config_payload, sort_keys=False, allow_unicode=True).rstrip(),
            narrative=show_config,
        )

    if dry_run:
        # Dry-run performs cheap structure/topology/layer validation, but never
        # invokes AmberTools or loads an ML model.
        if parm7_override is not None:
            import parmed as pmd
            from mlmm.backends.mlmm_calc import validate_parmed_atom_order

            top = pmd.load_file(str(parm7_override))
            for prepared in _prepared_all_inputs:
                structure_path = prepared.source_path
                structure = pmd.load_file(str(structure_path))
                validate_parmed_atom_order(
                    structure,
                    top,
                    input_label=str(prepared.display_path),
                    topology_label=str(parm7_override),
                )
        elif _missing_ambertools_commands(_ambertools_command_paths()):
            raise click.ClickException(
                "[all] AmberTools commands tleap, antechamber, and parmchk2 are "
                "required when --parm is not supplied."
            )
        if (
            center_spec is None
            and detect_layer
            and model_pdb_override is None
        ):
            counts = _summarize_existing_bfactor_layers(
                _prepared_all_inputs[0].source_path
            )
            if counts.get("movable", 0) == 0 and counts.get("frozen", 0) == 0:
                raise click.ClickException(
                    "[all] --detect-layer requires 0/10/20 B-factor layers when "
                    "extraction is skipped and --model-pdb is absent."
                )
        _echo(
            "[all] Dry-run validation passed: structure normalization, layer "
            "metadata, topology atom count/order, and required tools were checked. "
            "No calculation stage was executed.",
            narrative=True,
        )
        _echo(
            "[all] Planned stages: extract -> mm_parm -> optional scan -> path_opt/path_search -> optional tsopt/freq/dft.",
            narrative=True,
        )
        _emit_final_summary(out_dir, time_start, manifest)
        # Prepared-input temp trees are session-owned (own_cleanup, above); the
        # run's ctx.call_on_close(session.close) frees them on this dry-run
        # return as well, so no branch-local cleanup is needed here.
        return

    out_dir = out_dir.resolve()
    work_dir = out_dir / WORK_DIRNAME  # pipeline-wide scratch (safe to rm -rf)
    session.resources.own_exclusive_lock(work_dir / ".run.lock")
    input_paths = _materialize_all_coordinate_inputs(
        _prepared_all_inputs,
        work_dir,
    )
    # Declare the run's public root deliverables up front so their pre-run
    # baseline is captured before any producer writes.  A stale file from an
    # earlier invocation that this run does not rewrite stays unclaimed and is
    # excluded from the current-run key_output_files.
    for _public_name in (
        "summary.log",
        "summary.json",
        "mep_trj.xyz",
        "mep.xyz",
        "mep.pdb",
        "mep.cif",
        "ml_region.pdb",
        "ml_region_without_linkH.xyz",
        "ml_region_with_linkH.xyz",
        "ml_region_without_linkH.pdb",
        "ml_region_with_linkH.pdb",
        "energy_diagram_MEP.png",
        "mep_plot.png",
        "irc_plot_all.png",
        "energy_diagram_MLIP_all.png",
        "energy_diagram_G_MLIP_all.png",
        "energy_diagram_DFT_all.png",
        "energy_diagram_G_DFT_plus_MLIP_all.png",
    ):
        _declare_public_output(manifest, out_dir, out_dir / _public_name)

    def _write_public_segment_diagram(
        prefix: Path,
        **diagram_kwargs: Any,
    ) -> Optional[Dict[str, Any]]:
        """Write a per-segment energy diagram as a current-run public output.

        Mirrors p2r's ``_write_public_energy_diagram``: the exact ``.png``
        destination is producer-declared (and claimed) before it is written so
        the current run's segment diagrams reappear in summary.json's
        ``key_output_files`` (which now surfaces only declared current-run
        outputs, no directory discovery). Every prefix routed here is a
        descendant of ``<out_dir>/<SEGMENTS_DIRNAME>/seg_NN/`` (verified at each
        call site), so it satisfies the public-output layout requirement; the
        root aggregate/MEP diagrams are written outside this layout and are not
        routed through this helper.
        """
        destination = Path(prefix).with_suffix(".png")
        if manifest is not None:
            _declare_public_output(manifest, out_dir, destination)
        payload = _write_segment_energy_diagram(prefix, **diagram_kwargs)
        if manifest is not None:
            _claim_public_output(manifest, out_dir, destination)
        return payload

    pockets_dir = work_dir / "pockets"
    # MEP-engine raw output is scratch under _work/; only its moved products reach root.
    path_dir = work_dir / ("path_search" if refine_path else "path_opt")
    scan_dir = _resolve_override_dir(work_dir / "scan", scan_out_dir)  # for single-structure scan mode
    # One monotonic Stage numbering for the whole pipeline: Stage 1 extraction
    # (+ lettered preparation sub-stages 1b/1c/1d), 2 MEP search, 3 merge,
    # 4 post-processing. Banners read "Stage N/{stage_total}". Stage 4 only
    # runs when at least one of tsopt/thermo/dft is requested, so the
    # denominator drops to 3 on a default run to avoid an "N/4" that never
    # reaches 4.
    stage_total = 4 if (do_tsopt or do_thermo or do_dft) else 3
    ensure_dir(out_dir)
    if not single_tsopt_mode:
        ensure_dir(path_dir)  # path_search might be skipped only in tsopt-only mode

    # Preflight: add_elem_info only for inputs lacking element fields
    # → Create fixed copies under a temporary folder inside out_dir (used ONLY for extraction)
    elem_tmp_dir = work_dir / "add_elem_info"
    inputs_for_extract: List[Path] = []
    elem_fix_echo=False
    for input_ordinal, p in enumerate(input_paths, start=1):
        if _pdb_needs_elem_fix(p):
            if elem_fix_echo==False:
                _echo_section("====== [all] Preflight — add_elem_info (only when element fields are missing) ======")
                elem_fix_echo=True
            ensure_dir(elem_tmp_dir)
            out_p = _element_fix_path(elem_tmp_dir, p, input_ordinal)
            try:
                _assign_elem_info(str(p), str(out_p), overwrite=False)
                _echo(f"[all] add_elem_info: fixed elements → {out_p}")
                inputs_for_extract.append(out_p)
            except SystemExit as e:
                code = getattr(e, "code", 1)
                _echo(f"[all] WARNING: add_elem_info exited with code {code} for {p}; using original.", err=True)
                inputs_for_extract.append(p.resolve())
            except Exception as e:
                _echo(f"[all] WARNING: add_elem_info failed for {p}: {e} — using original file.", err=True)
                inputs_for_extract.append(p.resolve())
        else:
            inputs_for_extract.append(p.resolve())

    extract_inputs = tuple(inputs_for_extract)
    skip_extract = center_spec is None or str(center_spec).strip() == ""

    # OOM hazard guard: skip_extract + --no-detect-layer + no --model-pdb collapses the
    # ML region to the entire input PDB, causing downstream ML/MM ONIOM to scale ML over
    # all atoms (OOM on enzyme-sized systems). Hard-fail with an actionable message
    # rather than silently running the doomed configuration.
    if skip_extract and (not detect_layer) and model_pdb_override is None:
        raise click.ClickException(
            "[all] Skipping extraction (no -c/--center) with --no-detect-layer requires "
            "--model-pdb. Otherwise the ML region collapses to the entire input PDB and "
            "downstream ML/MM ONIOM will treat every atom as ML (OOM hazard on enzyme-sized "
            "systems). Provide --model-pdb <ml_region.pdb>, or enable --detect-layer to use "
            "B-factor layer information from the input PDB."
        )

    # When inputs are XYZ and --ref-pdb is provided, use it for topology-requiring steps
    ref_pdb_for_topology: Optional[Path] = None
    if ref_pdb_cli is not None:
        ref_pdb_for_topology = ref_pdb_cli.resolve()
        _echo(f"[all] --ref-pdb provided: {ref_pdb_for_topology}")

    resolved_charge: Optional[int] = None
    pocket_outputs: List[Path] = []

    if skip_extract:
        _echo_section(
            f"====== [all] Stage 1/{stage_total} — Extraction skipped (no -c/--center); using full structures as pockets ======"
        )
        pocket_outputs = [p.resolve() for p in extract_inputs]
        _echo("[all] Pocket inputs (full structures):")
        for op in pocket_outputs:
            _echo(f"  - {op}")
        # Charge derivation when extraction is skipped:
        #  - --model-pdb provided → derive over that ML pocket PDB.
        #  - detect-layer with a layered input that actually has MM atoms (ML ⊊
        #    system) → derive the ML-region (B≈0) charge WITH cap correction.
        #    Deriving over the full input would misapply the whole-system charge
        #    to the ML region, and an electron-parity check cannot detect every
        #    such mismatch. Reuse extract's cap-corrected charge summary.
        #  - otherwise (the whole input is the model) → full input PDB.
        resolved_charge = None
        if model_pdb_override is not None:
            resolved_charge = _derive_charge_from_ligand_charge_when_extract_skipped(
                model_pdb_override, ligand_charge
            )
        elif detect_layer:
            _layer_counts = _summarize_existing_bfactor_layers(extract_inputs[0])
            if _layer_counts.get("movable", 0) > 0 or _layer_counts.get("frozen", 0) > 0:
                resolved_charge = _derive_ml_charge_from_layered_pdb(
                    extract_inputs[0], ligand_charge
                )
        if resolved_charge is None:
            resolved_charge = _derive_charge_from_ligand_charge_when_extract_skipped(
                extract_inputs[0], ligand_charge
            )
    else:
        _echo_section(
            f"====== [all] Stage 1/{stage_total} — Active-site pocket extraction ======"
        )
        ensure_dir(pockets_dir)
        pocket_stems = [p.stem for p in extract_inputs]
        for idx, p in enumerate(extract_inputs, start=1):
            suffix = f"_{idx:02d}" if pocket_stems.count(p.stem) > 1 else ""
            pocket_outputs.append(
                (pockets_dir / f"pocket_{p.stem}{suffix}.pdb").resolve()
            )

        try:
            ex_res = extract_api(
                complex_pdb=[str(p) for p in extract_inputs],
                center=center_spec,
                output=[str(p) for p in pocket_outputs],
                radius=float(radius),
                radius_het2het=float(radius_het2het),
                include_h2o=bool(include_h2o),
                exclude_backbone=bool(exclude_backbone),
                add_linkh=bool(add_linkh),
                selected_resn=selected_resn or "",
                modified_residue=modified_residue or "",
                ligand_charge=ligand_charge,
                verbose=True,  # extractor INFO now gated by the unified global -v level
            )
        except Exception as e:
            raise click.ClickException(f"[all] Extractor failed: {e}")

        _echo("[all] Pocket files:")
        for op in pocket_outputs:
            _echo(f"  - {op}")

        try:
            cs = ex_res.get("charge_summary", {})
            q_total = float(cs.get("total_charge", 0.0))
            q_prot = float(cs.get("protein_charge", 0.0))
            q_lig = float(cs.get("ligand_total_charge", 0.0))
            q_ion = float(cs.get("ion_total_charge", 0.0))
            _echo("")
            _echo("[all] Charge summary from extractor (model #1):")
            _echo(
                f"  Protein: {q_prot:+g},  Ligand: {q_lig:+g},  Ions: {q_ion:+g},  Total: {q_total:+g}"
            )
            resolved_charge = _round_charge_with_note(q_total)
        except Exception as e:
            raise click.ClickException(f"[all] Could not obtain ML-region charge from extractor: {e}")

    if charge_override is not None:
        q_int = int(charge_override)
        override_msg = (
            "[all] WARNING: -q/--charge override supplied; "
            f"forcing ML-region charge to {q_int:+d}"
        )
        if resolved_charge is not None:
            override_msg += f" (would otherwise use {int(resolved_charge):+d} from workflow)"
        _echo(override_msg)
    else:
        if resolved_charge is None:
            if yaml_model_charge is None:
                raise click.ClickException(
                    "[all] ML-region charge could not be resolved. Provide "
                    "-q/--charge, --ligand-charge, or calc.model_charge in YAML."
                )
            q_int = int(yaml_model_charge)
            _echo(
                "[all] ML-region charge from YAML "
                f"calc.model_charge: {q_int:+d}"
            )
        else:
            q_int = int(resolved_charge)

    # Stage 1b: ML-region definition (copy first pocket) and mm_parm on the first full input
    _echo_section("====== [all] Stage 1b — ML/MM preparation — ML region + parm7 ======")
    first_pocket = pocket_outputs[0]
    first_full_input = extract_inputs[0]
    pocket_for_ml_region = first_pocket
    pdb_for_mm_parm = first_full_input

    # ML region definition: use --model-pdb if provided, otherwise generate from pocket.
    # When extraction was skipped + detect-layer, the "pocket" is the whole input, so
    # write only the B≈0 ML atoms — otherwise ml_region.pdb is the full system and a
    # downstream stage that can't read B-factors collapses ML to the entire input
    # (sum_Z=45875 electron-count error at freq/dft).
    if model_pdb_override is not None:
        model_pdb_source = model_pdb_override.resolve()
        # Safeguard (a): when detect-layer is active, an override --model-pdb should be an
        # ML-only region (B≈0 atoms only). If it is instead a FULL layered system (B≈0 ML
        # atoms alongside MovableMM=10 / FrozenMM=20 atoms), downstream ML-region checks
        # count the ENTIRE system (huge sum_Z → electron-count error). Warn, then
        # materialize the link-free public selection copy below.
        if detect_layer:
            _layers = _summarize_existing_bfactor_layers(model_pdb_source)
            if _layers.get("ml", 0) > 0 and (_layers.get("movable", 0) + _layers.get("frozen", 0)) > 0:
                _echo(
                    f"[all] WARNING: --model-pdb {model_pdb_source} looks like a full layered system "
                    f"(ML={_layers['ml']}, MovableMM={_layers['movable']}, FrozenMM={_layers['frozen']}); "
                    f"downstream ML-region checks will count the ENTIRE system. Pass an ML-only PDB "
                    f"(B≈0 atoms only) as --model-pdb, or omit --model-pdb to auto-generate one."
                )
        ml_region_pdb = _write_ml_region_definition(
            model_pdb_source,
            out_dir / "ml_region.pdb",
        )
        _echo_detail(
            f"[all] ML region selection (--model-pdb {model_pdb_source}) → "
            f"{ml_region_pdb}{_ml_region_summary_suffix(ml_region_pdb)}"
        )
    else:
        ml_region_pdb = None
        if skip_extract and detect_layer:
            ml_region_pdb = _write_bfactor_ml_subset(pocket_for_ml_region, out_dir / "ml_region.pdb")
            if ml_region_pdb is not None:
                _echo_detail(f"[all] ML region definition (B≈0 subset from layered input) → {ml_region_pdb}{_ml_region_summary_suffix(ml_region_pdb)}")
        if ml_region_pdb is None:
            ml_region_pdb = _write_ml_region_definition(pocket_for_ml_region, out_dir / "ml_region.pdb")  # reusable deliverable (--model-pdb input for follow-up runs)
            _echo_detail(f"[all] ML region definition → {ml_region_pdb}{_ml_region_summary_suffix(ml_region_pdb)}")

    # mm_parm: use --parm if provided, otherwise run tleap
    if parm7_override is not None:
        real_parm7_path = parm7_override.resolve()
        _echo_detail(f"[all] parm7 (--parm override) → {real_parm7_path}")
    else:
        _echo_detail(f"[all] mm_parm source PDB → {pdb_for_mm_parm}")
        mm_dir = out_dir / "mm_parm"  # reusable deliverable (parm7/rst7 = --parm input for follow-up runs)
        missing_ambertools = _missing_ambertools_commands(
            _ambertools_command_paths()
        )
        if missing_ambertools:
            raise click.ClickException(
                "[preflight] Missing required command(s) for mm_parm "
                f"(AmberTools): {', '.join(missing_ambertools)}. "
                "Install them in the active environment or make them available "
                "on PATH."
            )
        real_parm7_path, real_rst7_path = _build_mm_parm7(
            pdb=pdb_for_mm_parm,
            ligand_charge_expr=ligand_charge,
            ligand_mult_expr=mm_ligand_mult,
            out_dir=mm_dir,
            ff_set=mm_ff_set,
            add_ter=mm_add_ter,
            keep_temp=mm_keep_temp,
            auto_disulfide=mm_auto_disulfide,
        )
        _echo_detail(f"[all] mm_parm outputs → parm7: {real_parm7_path.name}, rst7: {real_rst7_path.name}")

    # Write both model systems before any ML backend is loaded. The ML selection
    # is link-free; link H atoms are generated exclusively from parm7 bonds that
    # cross the ML/MM boundary.
    from mlmm.workflows.dft import (
        _prepare_ml_region_workspace,
        write_ml_region_pdb_pair,
        write_ml_region_xyz_pair,
    )

    structure_calc_cfg = resolved_calc_template.materialize()
    structure_calc_cfg.update(
        {
            "model_charge": int(q_int),
            "model_mult": int(spin),
        }
    )
    region_workspace = _prepare_ml_region_workspace(
        input_pdb=pdb_for_mm_parm,
        coordinate_path=(
            _prepared_all_inputs[0].geom_path
            if _prepared_all_inputs[0].geom_path.suffix.lower() == ".xyz"
            else None
        ),
        real_parm7=real_parm7_path,
        model_pdb=ml_region_pdb,
        link_mlmm=structure_calc_cfg.get("link_mlmm"),
        link_atom_method=str(
            structure_calc_cfg.get("link_atom_method") or "scaled"
        ).lower(),
        use_cmap=bool(structure_calc_cfg.get("use_cmap", True)),
        calc_kwargs=structure_calc_cfg,
    )
    try:
        ml_without_link, ml_with_link = write_ml_region_xyz_pair(
            region_workspace,
            out_dir,
        )
        n_model_atoms = len(region_workspace.atoms_model)
        n_link_atoms = len(region_workspace.link_pairs)
        link_source = (
            "manual link_mlmm"
            if structure_calc_cfg.get("link_mlmm") is not None
            else "parm7 boundary bonds"
        )
        _echo_detail(
            f"[all] ML structure without link H ({n_model_atoms} atoms) → "
            f"{ml_without_link}"
        )
        _echo_detail(
            f"[all] ML structure with link H ({n_model_atoms} + {n_link_atoms}; "
            f"links from {link_source}) → {ml_with_link}"
        )
        if _original_input_paths[0].suffix.lower() == ".pdb":
            ml_without_link_pdb, ml_with_link_pdb = write_ml_region_pdb_pair(
                region_workspace,
                out_dir,
                xyz_paths=(ml_without_link, ml_with_link),
            )
            _echo_detail(
                "[all] PDB companions for the two ML structures → "
                f"{ml_without_link_pdb}, {ml_with_link_pdb}"
            )
        if region_workspace.link_pairs:
            _echo_detail(
                "[all] Link-H boundary pairs (1-based full-system parm7 ML-MM): "
                + ", ".join(
                    f"{ml_idx}-{mm_idx}"
                    for ml_idx, mm_idx in region_workspace.link_pairs
                )
            )
    finally:
        region_workspace.cleanup()

    # define-layer: assign 3-layer B-factors to each full-system PDB
    _echo_section("====== [all] Stage 1c — define-layer — assign 3-layer B-factors to full-system PDBs ======")
    layered_dir = out_dir / "layered"  # deliverable (B-factor-layered PDBs for inspection / reuse)
    ensure_dir(layered_dir)
    layer_source_dir = work_dir / "layer_sources"
    layered_inputs: List[Path] = []

    def _coordinate_layer_source(full_path: Path, index: int) -> Path:
        """Return a PDB with this endpoint's coordinates and shared metadata."""
        if full_path.suffix.lower() == ".pdb":
            return full_path
        if ref_pdb_for_topology is None:
            raise click.ClickException(
                f"[all] {full_path.name} requires --ref-pdb for layer assignment."
            )
        ensure_dir(layer_source_dir)
        source_pdb = layer_source_dir / f"endpoint_{index + 1:02d}.pdb"
        convert_xyz_to_pdb(full_path, ref_pdb_for_topology, source_pdb)
        return source_pdb

    # With extraction skipped, --detect-layer always reads the input layers.
    # An explicit --model-pdb owns ML membership while the input B-factors
    # retain their movable/frozen MM assignments. Recomputing either case here
    # would replace the user's layer boundary with define-layer's default
    # radius.
    honor_input_bfactors = bool(skip_extract and detect_layer)
    if honor_input_bfactors:
        _echo_detail(
            "[all] Extraction skipped and --detect-layer is on; honoring input PDB "
            "B-factor layer encoding (ML=0/MovableMM=10/FrozenMM=20)."
        )
        for idx, full_pdb in enumerate(extract_inputs):
            pdb_for_layer = _coordinate_layer_source(full_pdb, idx)
            counts = _summarize_existing_bfactor_layers(pdb_for_layer)
            _echo_detail(
                f"[all] define-layer [{idx}]: {full_pdb.name} (input B-factor layers honored)  "
                f"(ML={counts['ml']}, MovableMM={counts['movable']}, FrozenMM={counts['frozen']})"
            )
            layered_inputs.append(pdb_for_layer)
    else:
        layer_sources = [
            _coordinate_layer_source(full_pdb, idx)
            for idx, full_pdb in enumerate(extract_inputs)
        ]
        layer_stems = [source.stem for source in layer_sources]
        for idx, (full_pdb, pdb_for_layer) in enumerate(
            zip(extract_inputs, layer_sources)
        ):
            # When --ref-pdb is given and input is not PDB, use ref_pdb for define-layer
            suffix = (
                f"_{idx + 1:02d}"
                if layer_stems.count(pdb_for_layer.stem) > 1
                else ""
            )
            out_layered = (
                layered_dir / f"{pdb_for_layer.stem}{suffix}_layered.pdb"
            )
            try:
                layer_info = _define_layers(
                    input_pdb=pdb_for_layer,
                    output_pdb=out_layered,
                    model_pdb=ml_region_pdb,
                )
                _echo_detail(f"[all] define-layer [{idx}]: {full_pdb.name} → {out_layered.name}  "
                             f"(ML={len(layer_info.get('ml_indices', []))}, "
                             f"MovableMM={len(layer_info.get('movable_mm_indices', []))}, "
                             f"FrozenMM={len(layer_info.get('frozen_indices', []))})")
                layered_inputs.append(out_layered)
            except Exception as e:
                _echo(f"[all] WARNING: define-layer failed for {full_pdb.name}: {e}", err=True)
                _echo("[all] Falling back to original PDB (no B-factor layers).", err=True)
                layered_inputs.append(full_pdb)

    # Other path: single-structure + --tsopt True (and NO scan-lists) → TSOPT-only mode
    if single_tsopt_mode:
        _echo_section("====== [all] TSOPT-only single-structure mode ======")
        irc_trj_for_all: List[Tuple[Path, bool]] = []
        tsroot = out_dir / SEGMENTS_DIRNAME / "seg_01"
        ensure_dir(tsroot)

        # Use the layered full-system PDB as TS initial guess
        layered_pdb = layered_inputs[0]
        # When --ref-pdb is given and input is XYZ, copy the XYZ next to the layered PDB
        # so that _run_tsopt_on_hei can use XYZ (full precision) + layered PDB (topology)
        if (
            ref_pdb_for_topology is not None
            and _original_input_paths[0].suffix.lower() == ".xyz"
        ):
            xyz_companion = layered_pdb.with_suffix(".xyz")
            if not xyz_companion.exists():
                shutil.copy2(_prepared_all_inputs[0].geom_path, xyz_companion)
                _echo(f"[all] Copied XYZ input → {xyz_companion} (full precision for tsopt)")
        # TS optimization
        ts_pdb, g_ts = _run_tsopt_on_hei(
            layered_pdb,
            q_int,
            spin,
            real_parm7_path,
            ml_region_pdb,
            detect_layer,
            args_yaml,
            tsroot,
            tsopt_opt_mode_default,
            resolved_calc_template=resolved_calc_template,
            overrides=tsopt_overrides,
            backend=backend,
            embedcharge=embedcharge,
            embedcharge_cutoff=embedcharge_cutoff,
            embedcharge_explicit=embedcharge_explicit,
            link_atom_method=link_atom_method,
            mm_backend=mm_backend,
            use_cmap=use_cmap,
            ref_pdb=layered_pdb,
        )

        # EulerPC IRC & map endpoints (no segment endpoints exist → fallback mapping)
        irc_pocket_ref = ref_pdb_for_topology if ref_pdb_for_topology is not None else first_pocket
        irc_res = _irc_and_match(seg_idx=1,
                                 seg_dir=tsroot,
                                 ref_pdb_for_seg=ts_pdb,
                                 seg_pocket_pdb=irc_pocket_ref,
                                 g_ts=g_ts,
                                 q_int=q_int,
                                 spin=spin,
                                 resolved_calc_template=resolved_calc_template,
                                 real_parm7=real_parm7_path,
                                 model_pdb=ml_region_pdb,
                                 detect_layer=detect_layer,
                                 backend=backend,
                                 embedcharge=embedcharge,
                                 embedcharge_cutoff=embedcharge_cutoff,
                                 embedcharge_explicit=embedcharge_explicit,
                                 link_atom_method=link_atom_method,
                                 mm_backend=mm_backend,
                                 use_cmap=use_cmap,
                                 irc_step_size=irc_step_size,
                                 irc_never_stop=irc_never_stop,
                                 session=session,
                                 args_yaml=args_yaml)
        gL = irc_res["left_min_geom"]
        gR = irc_res["right_min_geom"]
        gT = irc_res["ts_geom"]
        irc_plot_path = irc_res.get("irc_plot")
        irc_trj_path = irc_res.get("irc_trj")
        if irc_trj_path:
            try:
                irc_trj_for_all.append((Path(irc_trj_path), bool(irc_res.get("reverse_irc", False))))
            except Exception:
                logger.debug("Failed to append IRC trajectory path", exc_info=True)

        # Ensure MLIP energies
        eL = float(gL.energy)
        eT = float(gT.energy)
        eR = float(gR.energy)

        # In this mode ONLY: assign Reactant/Product so that higher-energy end is the Reactant
        if eL >= eR:
            g_react, e_react = gL, eL
            g_prod,  e_prod  = gR, eR
        else:
            g_react, e_react = gR, eR
            g_prod,  e_prod  = gL, eL

        endpoint_assignment = {
            "policy": "higher_energy_endpoint_as_reactant",
            "chemical_direction_known": False,
            "left_role": "reactant" if eL >= eR else "product",
            "right_role": "product" if eL >= eR else "reactant",
            "left_energy_hartree": float(eL),
            "right_energy_hartree": float(eR),
            "tie_rule": (
                "on equal energy (eL == eR) the left endpoint is the reactant"
            ),
        }

        # Save XYZ (full precision) + PDB (companion) and run endpoint-opt
        struct_dir = tsroot / "structures"
        ensure_dir(struct_dir)
        pocket_ref = ref_pdb_for_topology if ref_pdb_for_topology is not None else first_pocket
        xR_irc, pR_irc = _save_single_geom_for_tools(g_react, pocket_ref, struct_dir, "reactant_irc")
        xT, pT         = _save_single_geom_for_tools(gT,       pocket_ref, struct_dir, "ts")
        xP_irc, pP_irc = _save_single_geom_for_tools(g_prod,   pocket_ref, struct_dir, "product_irc")

        endpoint_opt_dir = tsroot / "endpoint_opt"
        ensure_dir(endpoint_opt_dir)

        # Map IRC left/right Hessians → R/P endpoint (left=forward, right=backward)
        from mlmm.io.hessian_cache import (
            clear as _clear_hess_cache,
            discard as _hess_discard,
            load as _hess_load,
            store as _hess_store,
        )
        _react_hk = "irc_left" if eL >= eR else "irc_right"
        _prod_hk  = "irc_right" if eL >= eR else "irc_left"

        _hess_discard("irc_endpoint")
        _c = _hess_load(_react_hk)
        if _c:
            _hess_store("irc_endpoint", _c["hessian"], active_dofs=_c.get("active_dofs"), meta=_c.get("meta"), identity=_c.get("identity"))
        # Fail-closed endpoint-opt convergence (None if the opt could not run).
        _react_opt_conv: Optional[bool] = None
        try:
            g_react, _, _react_opt_conv = _run_opt_for_state(
                pR_irc, q_int, spin, real_parm7_path, ml_region_pdb, detect_layer,
                endpoint_opt_dir / "R", args_yaml, endpoint_opt_mode_default,
                resolved_calc_template=resolved_calc_template,
                convert_files=post_convert_files_forward,
                backend=backend,
                embedcharge=embedcharge,
                embedcharge_cutoff=embedcharge_cutoff,
                embedcharge_explicit=embedcharge_explicit,
                link_atom_method=link_atom_method,
                mm_backend=mm_backend,
                use_cmap=use_cmap,
                thresh=post_thresh_forward,
                reject_uphill=_reject_uphill_eff,
                xyz_path=xR_irc,
            )
        except Exception as e:
            _echo(
                f"[post] WARNING: Reactant endpoint optimization failed in TSOPT-only mode: {e}",
                err=True,
            )
            _react_opt_conv = None

        _hess_discard("irc_endpoint")
        _c = _hess_load(_prod_hk)
        if _c:
            _hess_store("irc_endpoint", _c["hessian"], active_dofs=_c.get("active_dofs"), meta=_c.get("meta"), identity=_c.get("identity"))
        _prod_opt_conv: Optional[bool] = None
        try:
            g_prod, _, _prod_opt_conv = _run_opt_for_state(
                pP_irc, q_int, spin, real_parm7_path, ml_region_pdb, detect_layer,
                endpoint_opt_dir / "P", args_yaml, endpoint_opt_mode_default,
                resolved_calc_template=resolved_calc_template,
                convert_files=post_convert_files_forward,
                backend=backend,
                embedcharge=embedcharge,
                embedcharge_cutoff=embedcharge_cutoff,
                embedcharge_explicit=embedcharge_explicit,
                link_atom_method=link_atom_method,
                mm_backend=mm_backend,
                use_cmap=use_cmap,
                thresh=post_thresh_forward,
                reject_uphill=_reject_uphill_eff,
                xyz_path=xP_irc,
            )
        except Exception as e:
            _echo(
                f"[post] WARNING: Product endpoint optimization failed in TSOPT-only mode: {e}",
                err=True,
            )
            _prod_opt_conv = None
        shutil.rmtree(endpoint_opt_dir, ignore_errors=True)
        _echo_detail("[endpoint-opt] Clean endpoint-opt working dir.")

        xR, pR = _save_single_geom_for_tools(g_react, pocket_ref, struct_dir, "reactant")
        xP, pP = _save_single_geom_for_tools(g_prod,   pocket_ref, struct_dir, "product")
        e_react = float(g_react.energy)
        e_prod = float(g_prod.energy)

        # ML/MM energy diagram (R, TS, P)
        mlip_prefix = tsroot / "energy_diagram_MLIP"
        mlip_diag = _write_public_segment_diagram(
            mlip_prefix,
            labels=["R", "TS", "P"],
            energies_eh=[e_react, eT, e_prod],
            title_note="(MLIP, TSOPT/IRC)",
        )
        g_mlip_diag = None
        dft_diag = None
        g_dft_mlip_diag = None

        # ── Release GPU memory before freq/thermo/DFT ──
        _irc_lease = irc_res.get("calculator_lease")
        if _irc_lease is not None:
            _irc_lease.release()
        for _g in (gL, gR, gT, g_react, g_prod):
            if _g is not None and hasattr(_g, "calculator"):
                _g.calculator = None
        gc.collect()
        if torch.cuda.is_available():
            torch.cuda.empty_cache()

        # Thermochemistry (ML/MM) Gibbs
        thermo_payloads: Dict[str, Dict[str, Any]] = {}
        GR = GT = GP = None
        eR_dft = eT_dft = eP_dft = None
        _dft_all_ok = not do_dft
        GR_dftMLIP = GT_dftMLIP = GP_dftMLIP = None
        freq_root = _resolve_override_dir(tsroot / "freq", freq_out_dir)
        dft_root = _resolve_override_dir(tsroot / "dft", dft_out_dir)

        if do_thermo:
            _echo_detail("[thermo] Single TSOPT: freq on TS/R/P")
            tT = _run_freq_for_state(pT, q_int, spin, real_parm7_path, ml_region_pdb, detect_layer,
                                     freq_root / "TS", args_yaml, overrides=freq_overrides,
                                     backend=backend, embedcharge=embedcharge, embedcharge_cutoff=embedcharge_cutoff,
                                     embedcharge_explicit=embedcharge_explicit,
                                     link_atom_method=link_atom_method, mm_backend=mm_backend, use_cmap=use_cmap, xyz_path=xT)
            _clear_hess_cache()  # TS Hessian consumed; R/P need exact computation
            tR = _run_freq_for_state(pR, q_int, spin, real_parm7_path, ml_region_pdb, detect_layer,
                                     freq_root / "R", args_yaml, overrides=freq_overrides,
                                     backend=backend, embedcharge=embedcharge, embedcharge_cutoff=embedcharge_cutoff,
                                     embedcharge_explicit=embedcharge_explicit,
                                     link_atom_method=link_atom_method, mm_backend=mm_backend, use_cmap=use_cmap, xyz_path=xR)
            tP = _run_freq_for_state(pP, q_int, spin, real_parm7_path, ml_region_pdb, detect_layer,
                                     freq_root / "P", args_yaml, overrides=freq_overrides,
                                     backend=backend, embedcharge=embedcharge, embedcharge_cutoff=embedcharge_cutoff,
                                     embedcharge_explicit=embedcharge_explicit,
                                     link_atom_method=link_atom_method, mm_backend=mm_backend, use_cmap=use_cmap, xyz_path=xP)
            thermo_payloads = {"R": tR, "TS": tT, "P": tP}
            GR = _thermo_gibbs_ha(tR)
            GT = _thermo_gibbs_ha(tT)
            GP = _thermo_gibbs_ha(tP)
            # build the MLIP Gibbs diagram ONLY when every requested state
            # returned a finite FREQ free energy. A failed/partial FREQ (empty or
            # missing field) must NOT be substituted by the MLIP electronic energy
            # (e_react/eT/e_prod) into a Gibbs-named result.
            if all(value is not None for value in (GR, GT, GP)):
                try:
                    g_mlip_diag = _write_public_segment_diagram(
                        tsroot / "energy_diagram_G_MLIP",
                        labels=["R", "TS", "P"],
                        energies_eh=[GR, GT, GP],
                        title_note="(Gibbs, MLIP)",
                        ylabel="ΔG (kcal/mol)",
                    )
                except Exception as e:
                    _echo(f"[thermo] WARNING: failed to build Gibbs diagram: {e}", err=True)
            else:
                _echo(
                    "[thermo] WARNING: one or more R/TS/P FREQ free energies are "
                    "unavailable; MLIP Gibbs diagram skipped (no MLIP-energy "
                    "substitution).",
                    err=True,
                )

        # DFT & DFT//MLIP/MM
        if do_dft:
            # Frequency subprocess parsing may
            # have re-bound heavy refs onto cli-frame Geometry locals.
            # Two layers (null calculator + rebind local to None) prevent
            # closure/hook capture from resurrecting the model before the
            # DFT subprocess fork. `del locals[name]` is a CPython no-op
            # (locals returns a *copy* of the frame namespace), so the
            # rebind to None below is the only mechanism that actually
            # decrements the heavy refs before gc.collect + empty_cache.
            # NOTE: thermo_payloads is deliberately NOT nulled here — the
            # DFT//MLIP/MM Gibbs-diagram block below and the
            # segment_log num_imag write both read from it.
            for _g in (gL, gR, gT, g_react, g_prod):
                if _g is not None and hasattr(_g, "calculator"):
                    _g.calculator = None
            gL = gR = gT = g_react = g_prod = None
            gc.collect()
            if torch.cuda.is_available():
                torch.cuda.empty_cache()

            _echo_detail("[dft] Single TSOPT: DFT on R/TS/P")
            dR = _run_dft_for_state(pR, q_int, spin, real_parm7_path, ml_region_pdb, detect_layer,
                                     dft_root / "R", args_yaml, func_basis=dft_func_basis_use, overrides=dft_overrides,
                                     backend=backend, embedcharge=embedcharge, embedcharge_cutoff=embedcharge_cutoff,
                                     embedcharge_explicit=embedcharge_explicit,
                                     link_atom_method=link_atom_method, mm_backend=mm_backend, use_cmap=use_cmap, xyz_path=xR)
            dT = _run_dft_for_state(pT, q_int, spin, real_parm7_path, ml_region_pdb, detect_layer,
                                     dft_root / "TS", args_yaml, func_basis=dft_func_basis_use, overrides=dft_overrides,
                                     backend=backend, embedcharge=embedcharge, embedcharge_cutoff=embedcharge_cutoff,
                                     embedcharge_explicit=embedcharge_explicit,
                                     link_atom_method=link_atom_method, mm_backend=mm_backend, use_cmap=use_cmap, xyz_path=xT)
            dP = _run_dft_for_state(pP, q_int, spin, real_parm7_path, ml_region_pdb, detect_layer,
                                     dft_root / "P", args_yaml, func_basis=dft_func_basis_use, overrides=dft_overrides,
                                     backend=backend, embedcharge=embedcharge, embedcharge_cutoff=embedcharge_cutoff,
                                     embedcharge_explicit=embedcharge_explicit,
                                     link_atom_method=link_atom_method, mm_backend=mm_backend, use_cmap=use_cmap, xyz_path=xP)
            eR_dft = _dft_energy_ha(dR)
            eT_dft = _dft_energy_ha(dT)
            eP_dft = _dft_energy_ha(dP)
            # Match the per-segment gate: a non-finite DFT energy must not reach a diagram.
            _dft_all_ok = all(
                e is not None and np.isfinite(e) for e in (eR_dft, eT_dft, eP_dft)
            )
            if not _dft_all_ok:
                _failed_states = [s for s, e in zip(["R", "TS", "P"], [eR_dft, eT_dft, eP_dft]) if e is None]
                _echo(f"[dft] WARNING: DFT failed for state(s): {', '.join(_failed_states)}. Skipping DFT diagrams.", err=True)
            if _dft_all_ok:
                try:
                    dft_diag = _write_public_segment_diagram(
                        tsroot / "energy_diagram_DFT",
                        labels=["R", "TS", "P"],
                        energies_eh=[eR_dft, eT_dft, eP_dft],
                        title_note=f"({dft_method_fallback})",
                    )
                except Exception as e:
                    _echo(f"[dft] WARNING: failed to build DFT diagram: {e}", err=True)

            # Build the DFT//MLIP/MM Gibbs diagram only from the subtractive
            # DFT//MLIP/MM electronic total plus a finite FREQ thermal correction.
            # Raw model-region DFT energy and 0.0 are never substitutes.
            eR_dft_mlmm = _dft_total_mlmm_energy_ha(dR)
            eT_dft_mlmm = _dft_total_mlmm_energy_ha(dT)
            eP_dft_mlmm = _dft_total_mlmm_energy_ha(dP)
            dG_R = _thermo_correction_ha(thermo_payloads.get("R"))
            dG_T = _thermo_correction_ha(thermo_payloads.get("TS"))
            dG_P = _thermo_correction_ha(thermo_payloads.get("P"))
            if do_thermo and all(
                value is not None
                for value in (
                    eR_dft_mlmm,
                    eT_dft_mlmm,
                    eP_dft_mlmm,
                    dG_R,
                    dG_T,
                    dG_P,
                )
            ):
                try:
                    GR_dftMLIP = eR_dft_mlmm + dG_R
                    GT_dftMLIP = eT_dft_mlmm + dG_T
                    GP_dftMLIP = eP_dft_mlmm + dG_P
                    g_dft_mlip_diag = _write_public_segment_diagram(
                        tsroot / "energy_diagram_G_DFT_plus_MLIP",
                        labels=["R", "TS", "P"],
                        energies_eh=[GR_dftMLIP, GT_dftMLIP, GP_dftMLIP],
                        title_note="(Gibbs, DFT//MLIP/MM)",
                        ylabel="ΔG (kcal/mol)",
                    )
                except Exception as e:
                    _echo(f"[dft//mlip] WARNING: failed to build DFT//MLIP/MM Gibbs diagram: {e}", err=True)
            elif do_thermo:
                _echo(
                    "[dft//mlip] WARNING: a subtractive DFT//MLIP/MM electronic "
                    "total or FREQ thermal correction is unusable; Gibbs diagram "
                    "skipped (no raw-model DFT or 0.0 substitution).",
                    err=True,
                )

        # Summary.yaml / summary.log for TSOPT-only mode
        bond_cfg = dict(_path_search.BOND_KW)
        bond_summary = ""
        try:
            changed, bond_summary = _path_search._has_bond_change(g_react, g_prod, bond_cfg)
            if not changed:
                bond_summary = "(no covalent changes detected)"
        except Exception as exc:
            # Without this the same string is published for "checked, nothing changed" and
            # "the check itself failed", so summary.json/summary.log assert no reaction
            # occurred for a run whose bond analysis never ran.
            _echo(
                f"[all] WARNING: bond-change detection failed ({exc}); reporting "
                "'(no covalent changes detected)' without having compared the endpoints.",
                err=True,
            )
            bond_summary = "(no covalent changes detected)"

        barrier = (eT - e_react) * AU2KCALPERMOL
        delta = (e_prod - e_react) * AU2KCALPERMOL

        from mlmm.workflows._all_helpers import promote_diag_for_root
        energy_diagrams: List[Dict[str, Any]] = []
        for diag, stem in (
            (mlip_diag, "energy_diagram_MLIP"),
            (g_mlip_diag, "energy_diagram_G_MLIP"),
            (dft_diag, "energy_diagram_DFT"),
            (g_dft_mlip_diag, "energy_diagram_G_DFT_plus_MLIP"),
        ):
            promoted = promote_diag_for_root(diag, stem, out_dir)
            if promoted is not None:
                energy_diagrams.append(promoted)

        n_irc_frames: Optional[int] = None
        if irc_trj_path:
            try:
                n_irc_frames = len(
                    read_xyz_as_blocks(Path(irc_trj_path), strict=True)
                )
            except Exception as exc:
                _echo(
                    "[post] WARNING: could not count complete IRC trajectory "
                    f"frames: {exc}",
                    err=True,
                )

        summary = {
            "out_dir": str(tsroot),
            "n_images": n_irc_frames,
            "n_segments": 1,
            "endpoint_assignment": endpoint_assignment,
            "segments": [
                {
                    "index": 1,
                    "tag": "seg_01",
                    "kind": "tsopt",
                    "barrier_kcal": float(barrier),
                    "delta_kcal": float(delta),
                    "bond_changes": bond_summary,
                    "endpoint_assignment": endpoint_assignment,
                }
            ],
            "energy_diagrams": list(energy_diagrams),
        }
        _enrich_summary(
            summary,
            version="",
            pipeline_mode="tsopt-only",
            out_dir=out_dir,
            manifest=manifest,
            mlip_backend=mlip_backend_resolved,
            mlip_model=mlip_model_resolved,
            mlip_precision=mlip_precision_resolved,
            charge=q_int,
            spin=spin,
            command=command_str,
            config={
                "refine_path": bool(refine_path),
                "tsopt": do_tsopt,
                "thermo": do_thermo,
                "dft": do_dft,
                "opt_mode": tsopt_opt_mode_default,
                "path_opt_mode": path_search_opt_mode,
                "post_opt_mode": tsopt_opt_mode_default,
                "ts_opt_mode": tsopt_opt_mode_default,
                "endpoint_opt_mode": endpoint_opt_mode_default,
                "mep_mode": mep_mode_kind,
                "dmf_correlated": dmf_correlated_effective,
            },
        )
        try:
            _publish_manifest_summary(
                out_dir / "summary.json",
                summary,
                manifest=manifest,
                out_dir=out_dir,
                mirrors=(tsroot / "summary.json",),
            )
        except Exception as e:
            _echo(f"[write] WARNING: failed to write summary.json: {e}", err=True)

        # Copy R/TS/P structures to out_dir/seg_01/
        try:
            _state_structs = {"R": pR, "TS": pT, "P": pP}
            _input_suffix = (
                _original_input_paths[0].suffix.lower()
                if _original_input_paths
                else ".xyz"
            )
            _seg_out = _copy_structures_to_seg_dir(
                _state_structs, out_dir, 1, _input_suffix,
                manifest=manifest,
            )
            _echo(f"[all] Wrote R/TS/P for segment 01 → {_seg_out}", narrative=True)
        except Exception as e:
            _echo(f"[all] WARNING: Failed to copy R/TS/P structures: {e}", err=True)

        segment_log: Dict[str, Any] = {
            "index": 1,
            "tag": "seg_01",
            "kind": "tsopt",
            "bond_changes": bond_summary,
            "post_dir": str(tsroot),
            "endpoint_assignment": endpoint_assignment,
        }
        if irc_plot_path:
            segment_log["irc_plot"] = str(irc_plot_path)
        if irc_trj_path:
            segment_log["irc_traj"] = str(irc_trj_path)
        # thread the per-direction IRC outcome so the aggregate
        # gates on convergence, not trajectory-file existence.
        _irc_outcome_seg = irc_res.get("irc_outcome")
        if isinstance(_irc_outcome_seg, dict):
            segment_log["irc"] = _irc_outcome_seg
        segment_log["endpoint_assignment"] = irc_res.get(
            "endpoint_assignment"
        )
        # record endpoint-opt convergence so a nonconverged endpoint
        # (whose geometry is still used for the diagram) does not silently
        # promote its segment to a usable success.
        segment_log["endpoint_opt"] = {
            "reactant_converged": _react_opt_conv,
            "product_converged": _prod_opt_conv,
        }
        if do_tsopt:
            tsopt_n_imag = (getattr(gT, "_tsopt_result", {}) or {}).get(
                "n_imaginary_modes"
            )
            if tsopt_n_imag is not None:
                segment_log["ts_imag"] = _ts_imag_record(
                    tsopt_n_imag,
                    (getattr(gT, "_tsopt_result", {}) or {}).get(
                        "imaginary_frequencies_cm"
                    ),
                )
        if do_thermo:
            n_imag = None
            try:
                n_imag = int(thermo_payloads.get("TS", {}).get("num_imag_freq"))
            except Exception:
                n_imag = None
            if n_imag is not None:
                # thermoanalysis.yaml carries `num_imag_freq` but no frequency
                # list, so this branch must not clobber the frequencies the
                # tsopt branch already published above.
                _prior_freqs = (segment_log.get("ts_imag") or {}).get(
                    "imag_freqs_cm"
                )
                segment_log["ts_imag"] = _ts_imag_record(
                    n_imag,
                    (thermo_payloads.get("TS") or {}).get(
                        "imaginary_frequencies_cm"
                    )
                    or _prior_freqs,
                )

        from mlmm.workflows._all_helpers import (
            build_energy_level_dict,
            build_thermo_symmetry_provenance,
        )
        _thermo_symmetry = build_thermo_symmetry_provenance(thermo_payloads)
        if _thermo_symmetry:
            segment_log["thermo_symmetry"] = _thermo_symmetry
        _structs = {"R": pR, "TS": pT, "P": pP}
        segment_log["mlip"] = build_energy_level_dict(
            labels=["R", "TS", "P"],
            energies_au=[e_react, eT, e_prod],
            ref_energy=e_react,
            au_to_kcal=AU2KCALPERMOL,
            diagram_path=str((tsroot / "energy_diagram_MLIP").with_suffix(".png")),
            structures=_structs,
        )
        if GR is not None and GT is not None and GP is not None:
            segment_log["gibbs_mlip"] = build_energy_level_dict(
                labels=["R", "TS", "P"],
                energies_au=[GR, GT, GP],
                ref_energy=GR,
                au_to_kcal=AU2KCALPERMOL,
                diagram_path=str((tsroot / "energy_diagram_G_MLIP").with_suffix(".png")),
                structures=_structs,
            )
        if eR_dft is not None and eT_dft is not None and eP_dft is not None:
            segment_log["dft"] = build_energy_level_dict(
                labels=["R", "TS", "P"],
                energies_au=[eR_dft, eT_dft, eP_dft],
                ref_energy=eR_dft,
                au_to_kcal=AU2KCALPERMOL,
                diagram_path=str((tsroot / "energy_diagram_DFT").with_suffix(".png")),
                structures=_structs,
            )
        if GR_dftMLIP is not None and GT_dftMLIP is not None and GP_dftMLIP is not None:
            segment_log["gibbs_dft_mlip"] = build_energy_level_dict(
                labels=["R", "TS", "P"],
                energies_au=[GR_dftMLIP, GT_dftMLIP, GP_dftMLIP],
                ref_energy=GR_dftMLIP,
                au_to_kcal=AU2KCALPERMOL,
                diagram_path=str((tsroot / "energy_diagram_G_DFT_plus_MLIP").with_suffix(".png")),
                structures=_structs,
            )

        summary_payload = {
            "root_out_dir": str(out_dir),
            "path_dir": str(tsroot),
            "path_module_dir": "tsopt_single",
            "pipeline_mode": "tsopt-only",
            "n_images": n_irc_frames,
            "n_segments": 1,
            "refine_path": bool(refine_path),
            "thresh": thresh,
            "thresh_post": thresh_post,
            "flatten": bool(flatten),
            "tsopt": do_tsopt,
            "thermo": do_thermo,
            "dft": do_dft,
            "dft_status": (
                "failed"
                if do_dft and not _dft_all_ok
                else ("converged" if do_dft else None)
            ),
            "opt_mode": tsopt_opt_mode_default,
            "post_opt_mode": tsopt_opt_mode_default,
            "ts_opt_mode": tsopt_opt_mode_default,
            "endpoint_opt_mode": endpoint_opt_mode_default,
            "mep_mode": "tsopt-only",
            "dmf_correlated": dmf_correlated_effective,
            "mlip_backend": mlip_backend_resolved,
            "mlip_model": mlip_model_resolved,
            "mlip_precision": mlip_precision_resolved,
            "status": summary.get("status"),
            "status_reasons": summary.get("status_reasons", []),
            "execution_status": summary.get("execution_status"),
            "scientific_status": summary.get("scientific_status"),
            "scientific_status_reasons": summary.get(
                "scientific_status_reasons", []
            ),
            "command": command_str,
            "charge": q_int,
            "spin": spin,
            "mep": {"n_images": 0, "n_segments": 1},
            "segments": summary.get("segments", []),
            "energy_diagrams": list(energy_diagrams),
            "post_segments": [segment_log],
            "key_files": {},
        }
        # Refresh summary.json with post_segments and key_output_files
        try:
            _enrich_summary(
                summary,
                version="",
                pipeline_mode="tsopt-only",
                out_dir=out_dir,
                manifest=manifest,
                mlip_backend=mlip_backend_resolved,
                mlip_model=mlip_model_resolved,
                mlip_precision=mlip_precision_resolved,
                charge=q_int,
                spin=spin,
                command=command_str,
                post_segments=[segment_log],
                config={
                    "refine_path": bool(refine_path),
                    "tsopt": do_tsopt,
                    "thermo": do_thermo,
                    "dft": do_dft,
                    "dft_status": summary_payload["dft_status"],
                    "opt_mode": tsopt_opt_mode_default,
                    "path_opt_mode": path_search_opt_mode,
                    "post_opt_mode": tsopt_opt_mode_default,
                    "ts_opt_mode": tsopt_opt_mode_default,
                    "endpoint_opt_mode": endpoint_opt_mode_default,
                    "dmf_correlated": dmf_correlated_effective,
                },
            )
            for key in (
                "status",
                "status_reasons",
                "execution_status",
                "scientific_status",
                "scientific_status_reasons",
            ):
                summary_payload[key] = summary.get(
                    key,
                    [] if key.endswith("_reasons") else None,
                )
            summary["post_segments"] = _json_safe([segment_log])
            # key_output_files is rebuilt from the current-run manifest inside
            # _enrich_summary (producer-declared claims only); no filesystem
            # discovery rebuild is needed here.
            _publish_manifest_summary(
                out_dir / "summary.json",
                summary,
                manifest=manifest,
                out_dir=out_dir,
                mirrors=(tsroot / "summary.json",),
            )
        except Exception as e:
            _echo(f"[write] WARNING: failed to refresh summary.json: {e}", err=True)

        try:
            write_summary_log(tsroot / "summary.log", summary_payload)
            commit_exact_bytes(
                out_dir / "summary.log",
                (tsroot / "summary.log").read_bytes(),
            )
        except Exception as e:
            _echo(f"[write] WARNING: failed to write summary.log: {e}", err=True)

        try:
            for stem in (
                "energy_diagram_MLIP",
                "energy_diagram_G_MLIP",
                "energy_diagram_DFT",
                "energy_diagram_G_DFT_plus_MLIP",
            ):
                src = tsroot / f"{stem}.png"
                if src.exists():
                    shutil.copy2(src, out_dir / f"{stem}_all.png")
        except Exception as e:
            _echo(f"[all] WARNING: failed to copy *_all diagrams: {e}", err=True)

        try:
            if irc_plot_path:
                irc_plot_src = Path(irc_plot_path)
                if irc_plot_src.exists():
                    shutil.copy2(irc_plot_src, out_dir / "irc_plot_all.png")
        except Exception as e:
            _echo(f"[all] WARNING: failed to copy irc_plot_all.png: {e}", err=True)

        _finalize_current_summary(
            out_dir / "summary.json",
            summary,
            manifest=manifest,
            out_dir=out_dir,
            mirrors=(tsroot / "summary.json",),
        )

        _echo_section("====== [all] TSOPT-only pipeline finished successfully ======")
        citation_post_segments = [segment_log]
        _emit_final_summary(
            out_dir,
            time_start,
            manifest,
            citation_payload=_all_method_citation_payload(),
        )
        return

    # Stage 1d: Optional scan (single-structure only) to build ordered pocket inputs
    pockets_for_path: List[Path]
    if is_single and has_scan:
        _echo_section("====== [all] Stage 1d — Staged scan on layered full-system PDB (single-structure mode) ======")
        ensure_dir(scan_dir)
        layered_pdb = Path(layered_inputs[0]).resolve()
        full_input_pdb = Path(input_paths[0]).resolve()
        # Use the layered full-system PDB for scan (no pocket index remapping needed)
        full_atom_meta = load_pdb_atom_metadata(full_input_pdb)
        # Honour --scan-one-based CLI toggle (None defaults to 1-based for backward compat).
        scan_one_based_use = True if scan_one_based is None else bool(scan_one_based)
        converted_scan_stages = _parse_scan_lists_literals(
            scan_lists_raw, atom_meta=full_atom_meta, one_based=scan_one_based_use,
        )
        scan_stage_literals: List[str] = []
        for stage in converted_scan_stages:
            scan_stage_literals.append(_format_scan_stage(stage))
        _echo("[all] Remapped --scan-lists indices from the full PDB to the pocket ordering.")
        scan_preopt_use = pre_opt if scan_preopt_override is None else bool(scan_preopt_override)
        scan_endopt_use = False if scan_endopt_override is None else bool(scan_endopt_override)
        scan_opt_mode_use = path_search_opt_mode

        scan_args: List[str] = [
            "-i", str(layered_pdb),
            "--parm", str(real_parm7_path),
            "-q", str(int(q_int)),
            "-m", str(int(spin)),
            "--out-dir", str(scan_dir),
            "--preopt" if scan_preopt_use else "--no-preopt",
            "--endopt" if scan_endopt_use else "--no-endopt",
            "--opt-mode", str(scan_opt_mode_use),
        ]
        scan_args.append("--detect-layer" if detect_layer else "--no-detect-layer")

        if dump_override_requested:
            scan_args.append("--dump" if dump else "--no-dump")

        # Forward the scan-indexing convention selected by --scan-one-based
        scan_args.append("--one-based" if scan_one_based_use else "--zero-based")

        _append_cli_arg(scan_args, "--max-step-size", scan_max_step_size)
        _append_cli_arg(scan_args, "--bias-k", scan_bias_k)
        _append_cli_arg(scan_args, "--relax-max-cycles", scan_relax_max_cycles)
        scan_args.extend(
            _build_scan_child_argv(
                explicit_params,
                convert_files=convert_files,
                thresh=thresh,
            )
        )
        if args_yaml is not None:
            scan_args.extend(["--config", str(args_yaml)])
        # Forward all converted --scan-lists (aligned to the pocket atom order)
        if scan_stage_literals:
            scan_args.append("--scan-lists")
            scan_args.extend(scan_stage_literals)

        from mlmm.workflows._all_helpers import append_backend_forwarding_args
        append_backend_forwarding_args(
            scan_args,
            backend=backend,
            embedcharge=embedcharge,
            embedcharge_cutoff=embedcharge_cutoff,
            embedcharge_explicit=embedcharge_explicit,
            link_atom_method=link_atom_method,
            mm_backend=mm_backend,
            use_cmap=use_cmap,
        )

        _echo_detail(
            f"[all] dispatch scan: input={layered_pdb.name}, "
            f"stages={len(scan_stage_literals)}, preopt={'yes' if scan_preopt_use else 'no'}, "
            f"endopt={'yes' if scan_endopt_use else 'no'}, out={scan_dir}"
        )
        _echo("[all] mlmm scan " + " ".join(scan_args))

        _run_cli_main("scan", _scan_cli.cli, scan_args, on_nonzero="raise", on_exception="raise", prefix="all")

        # Collect stage results — prefer XYZ (full precision), keep PDB as ref for topology
        stage_results: List[Path] = []
        stage_refs: List[Path] = []
        for st in sorted(scan_dir.glob("stage_*")):
            if not st.is_dir():
                continue
            xyz = st / "result.xyz"
            pdb = st / "result.pdb"
            if xyz.exists():
                stage_results.append(xyz.resolve())
                stage_refs.append(pdb.resolve() if pdb.exists() else layered_pdb)
            elif pdb.exists():
                stage_results.append(pdb.resolve())
                stage_refs.append(pdb.resolve())
        if not stage_results:
            raise click.ClickException("[all] No stage result files found under scan/.")
        _echo_detail("[all] Collected scan stage files:")
        for p in stage_results:
            _echo_detail(f"  - {p}")

        # Input series to path_search: [preopt result (if available), scan stage results ...]
        # When scan ran with --preopt, its optimized reactant geometry lives in
        # scan/preopt/result.xyz (full precision) or result.pdb.  Using this
        # avoids a redundant ~2000-cycle re-optimization inside path_search.
        preopt_xyz = scan_dir / "preopt" / "result.xyz"
        preopt_pdb = scan_dir / "preopt" / "result.pdb"
        if preopt_xyz.exists():
            init0_geom = preopt_xyz.resolve()
            init0_ref = layered_pdb          # layered PDB has authoritative B-factor layers
        elif preopt_pdb.exists():
            init0_geom = preopt_pdb.resolve()
            init0_ref = layered_pdb
        else:
            # No preopt output — fall back to original layered PDB
            init0_geom = layered_pdb
            init0_ref = layered_pdb
        pockets_for_path = [init0_geom] + stage_results
        refs_for_path = [init0_ref] + stage_refs
        _echo_detail(f"[all] Using scan initial endpoint: {init0_geom}")
    else:
        # Multi-structure standard route: use layered full-system PDBs
        pockets_for_path = list(layered_inputs)

    # --- Global pre-alignment for coordinate continuity across segments ---
    if not refine_path and len(pockets_for_path) >= 2:
        try:
            _echo("[all] Pre-aligning all input structures to first frame...")
            _align_dir = path_dir / "pre_align"
            ensure_dir(_align_dir)
            _bfs = read_bfactors_from_pdb(pockets_for_path[0])
            _fa = [i for i, bf in enumerate(_bfs) if bf >= 15.0]
            if _fa:
                _geoms = [geom_loader(str(p), coord_type="cart") for p in pockets_for_path]
                for _g in _geoms:
                    _g.freeze_atoms = np.array(_fa, dtype=int)
                _calc_kw = _stage_calc_kwargs(
                    resolved_calc_template,
                    input_pdb=pockets_for_path[0],
                    real_parm7=real_parm7_path,
                    model_pdb=ml_region_pdb,
                    charge=q_int,
                    spin=spin,
                    use_bfactor_layers=True,
                )
                _align_calc = _mlmm_calc(**_calc_kw)
                alignment_results = align_and_refine_sequence_inplace(
                    _geoms, shared_calc=_align_calc,
                    out_dir=_align_dir / "refine", verbose=True,
                )
                failed_pairs = alignment_failed_pair_indices(alignment_results)
                if failed_pairs:
                    raise click.ClickException(
                        "Input alignment did not converge for pair(s): "
                        + ", ".join(str(index) for index in failed_pairs)
                    )
                del _align_calc
                _new_pockets: List[Path] = []
                for _i, (_g, _orig) in enumerate(zip(_geoms, pockets_for_path)):
                    _xyz = _align_dir / f"{_i:03d}.xyz"
                    _xyz.write_text(_g.as_xyz() + "\n")
                    _pdb = _align_dir / f"{_i:03d}.pdb"
                    convert_xyz_to_pdb(_xyz, _orig, _pdb)
                    _new_pockets.append(_pdb)
                pockets_for_path = _new_pockets
                _echo("[all] Pre-alignment completed.")
        except Exception as e:
            raise click.ClickException(
                f"Path pre-alignment failed: {e}"
            ) from e

    # Stage 2: Path search on full-system layered PDBs
    mep_mode_label = mep_mode_kind.upper()
    if refine_path:
        _echo_section(
            f"====== [all] Stage 2/{stage_total} — MEP search on full-system "
            f"layered PDBs (recursive {mep_mode_label}) ======"
        )

        # Build path_search CLI args using *repeated* options (robust for Click)
        ps_args: List[str] = []

        # Inputs: single -i followed by all layered full-system PDBs
        ps_args.append("-i")
        for p in pockets_for_path:
            ps_args.append(str(p))

        # Charge & spin
        ps_args.extend(["-q", str(q_int)])
        ps_args.extend(["-m", str(int(spin))])
        ps_args.extend(["--parm", str(real_parm7_path)])
        # Layered PDBs have B-factors → --detect-layer is True by default in
        # the option declaration. Honor the user's explicit `--no-detect-layer`
        # by forwarding the chosen toggle instead of hardcoding `--detect-layer`.
        ps_args.append("--detect-layer" if detect_layer else "--no-detect-layer")

        # User-tunable parent defaults stay absent so child YAML remains the
        # effective middle layer. Pipeline-owned output paths are always set.
        ps_args.extend(
            _build_path_child_argv(
                explicit_params,
                include_opt_mode=True,
                mep_mode=mep_mode_kind,
                dmf_backend=dmf_backend,
                max_nodes=max_nodes,
                max_cycles=max_cycles,
                climb=climb,
                opt_mode=path_search_opt_mode,
                dump=dump,
                pre_opt=pre_opt,
                convert_files=convert_files,
                thresh=thresh,
            )
        )
        ps_args.extend(["--out-dir", str(path_dir)])
        if args_yaml is not None:
            ps_args.extend(["--config", str(args_yaml)])

        # Provide --ref-pdb for topology/B-factor info (one per input)
        # MUST use layered PDBs (with B-factor layer info) so that downstream
        # PDB conversion preserves ML/MovableMM/FrozenMM layer encoding.
        ps_args.append("--ref-pdb")
        if is_single and has_scan:
            # single+scan: use refs_for_path which maps to each pocket (XYZ→PDB ref)
            for ref in refs_for_path:
                ps_args.append(str(ref))
        else:
            for lp in layered_inputs:
                ps_args.append(str(lp))

        from mlmm.workflows._all_helpers import append_backend_forwarding_args
        append_backend_forwarding_args(
            ps_args,
            backend=backend,
            embedcharge=embedcharge,
            embedcharge_cutoff=embedcharge_cutoff,
            embedcharge_explicit=embedcharge_explicit,
            link_atom_method=link_atom_method,
            mm_backend=mm_backend,
            use_cmap=use_cmap,
        )

        _echo_detail(
            f"[all] dispatch path-search: inputs={len(pockets_for_path)}, "
            f"mode=recursive-{mep_mode_kind}, preopt={'yes' if pre_opt else 'no'}, "
            f"detect_layer={'yes' if detect_layer else 'no'}, out={path_dir}"
        )
        _echo("[all] mlmm path-search " + " ".join(ps_args))

        _run_cli_main("path_search", _path_search.cli, ps_args, on_nonzero="raise", on_exception="raise", prefix="all")
    else:
        # --no-refine-path: run path-opt between each adjacent pair and concatenate.
        _echo_section(
            f"====== [all] Stage 2/{stage_total} — MEP path-opt on full-system "
            f"layered PDBs (single-pass {mep_mode_label} per pair) ======"
        )

        if len(pockets_for_path) < 2:
            raise click.ClickException("[all] Need at least two structures for path-opt MEP concatenation.")

        ensure_dir(path_dir)
        combined_blocks: List[str] = []
        path_opt_segments: List[Dict[str, Any]] = []

        for pair_pos in range(len(pockets_for_path) - 1):
            # Array access remains zero-based; every public segment identifier
            # follows the documented one-based seg_01, seg_02, ... contract.
            seg_idx = pair_pos + 1
            p_left = pockets_for_path[pair_pos]
            p_right = pockets_for_path[pair_pos + 1]
            seg_tag = f"seg_{seg_idx:02d}"
            seg_out = path_dir / f"{seg_tag}_mep"
            ensure_dir(seg_out)

            po_args: List[str] = [
                "-i", str(p_left), str(p_right),
                "-q", str(q_int),
                "-m", str(int(spin)),
                "--parm", str(real_parm7_path),
            ]
            # When the single+scan route handed over XYZ pockets, forward the
            # matching layered-template ref PDBs so path-opt can overlay the XYZ
            # coordinates onto full ML/MM topology (path-search receives these too).
            if is_single and has_scan:
                po_args.extend([
                    "--ref-pdb", str(refs_for_path[pair_pos]),
                    "--ref-pdb", str(refs_for_path[pair_pos + 1]),
                ])
            # Forward the chosen --detect-layer/--no-detect-layer toggle
            # (default True). Hardcoded "--detect-layer" silently overrode
            # user's `--no-detect-layer` request.
            po_args.append("--detect-layer" if detect_layer else "--no-detect-layer")
            po_args.extend(
                _build_path_child_argv(
                    explicit_params,
                    include_opt_mode=False,
                    mep_mode=mep_mode_kind,
                    dmf_backend=dmf_backend,
                    max_nodes=max_nodes,
                    max_cycles=max_cycles,
                    climb=climb,
                    opt_mode=None,
                    dump=dump,
                    pre_opt=pre_opt,
                    convert_files=convert_files,
                    thresh=thresh,
                )
            )
            po_args.extend(["--out-dir", str(seg_out)])
            # Pipeline-owned machine contract: the aggregate reads this child's
            # real MEP convergence from result.json.
            po_args.append("--out-json")
            from mlmm.workflows._all_helpers import append_backend_forwarding_args
            append_backend_forwarding_args(
                po_args,
                backend=backend,
                embedcharge=embedcharge,
                embedcharge_cutoff=embedcharge_cutoff,
                embedcharge_explicit=embedcharge_explicit,
                link_atom_method=link_atom_method,
                mm_backend=mm_backend,
                use_cmap=use_cmap,
                args_yaml=args_yaml,
            )

            _echo_detail(
                f"[all] dispatch path-opt pair {seg_idx}/{len(pockets_for_path) - 1}: "
                f"mode={mep_mode_kind}, preopt={'yes' if pre_opt else 'no'}, "
                f"climb={'yes' if climb else 'no'}, out={seg_out}"
            )
            _echo("[all] mlmm path-opt " + " ".join(po_args))
            _run_cli_main("path_opt", _path_opt.cli, po_args, on_nonzero="raise", on_exception="raise", prefix="all")

            seg_converged = _read_path_opt_segment_converged(seg_out)

            # --- Post-processing per segment ---
            seg_trj = seg_out / "final_geometries_trj.xyz"
            if not seg_trj.exists():
                raise click.ClickException(
                    f"[all] path-opt segment {seg_idx} did not produce final_geometries_trj.xyz"
                )

            # Copy per-segment trajectory to path_dir
            try:
                seg_mep_trj = path_dir / f"mep_seg_{seg_idx:02d}_trj.xyz"
                shutil.copy2(seg_trj, seg_mep_trj)
                if pockets_for_path[0].suffix.lower() == ".pdb":
                    _path_search._maybe_convert_to_pdb(
                        seg_mep_trj,
                        ref_pdb_path=pockets_for_path[0],
                        out_path=path_dir / f"mep_seg_{seg_idx:02d}.pdb",
                    )
            except Exception as e:
                _echo(
                    f"[all] WARNING: failed to emit per-segment trajectory copies for segment {seg_idx:02d}: {e}",
                    err=True,
                )

            # Mirror HEI artifacts
            hei_src = seg_out / "hei.xyz"
            if hei_src.exists():
                try:
                    shutil.copy2(hei_src, path_dir / f"hei_seg_{seg_idx:02d}.xyz")
                    hei_pdb_src = seg_out / "hei.pdb"
                    if hei_pdb_src.exists():
                        shutil.copy2(hei_pdb_src, path_dir / f"hei_seg_{seg_idx:02d}.pdb")
                except Exception as e:
                    _echo(
                        f"[all] WARNING: failed to prepare HEI artifacts for segment {seg_idx:02d}: {e}",
                        err=True,
                    )

            # Parse trajectory blocks for concatenation and energy extraction
            raw_blocks = read_xyz_as_blocks(seg_trj, strict=True)
            blocks = ["\n".join(b) + "\n" for b in raw_blocks]
            if not blocks:
                raise click.ClickException(
                    f"[all] No frames read from path-opt segment {seg_idx} trajectory: {seg_trj}"
                )
            # Skip duplicate first frame for subsequent segments
            if pair_pos > 0:
                blocks = blocks[1:]
            combined_blocks.extend(blocks)

            # Segment energetics require one finite energy for every frame.
            energies_seg = _required_xyz_block_energies(
                raw_blocks,
                path=seg_trj,
                context=f"path-opt segment {seg_idx}",
            )

            # Parse first/last frame coordinates for bond-change detection
            first_last = None
            try:
                first_last = xyz_blocks_first_last(raw_blocks, path=seg_trj)
            except Exception as e:
                _echo(
                    f"[all] WARNING: failed to parse first/last frames for segment {seg_idx:02d}: {e}",
                    err=True,
                )

            path_opt_segments.append(
                {
                    "tag": seg_tag,
                    "energies": energies_seg,
                    "traj": seg_trj,
                    "inputs": (p_left, p_right),
                    "first_last": first_last,
                    # Child MEP signal; False/unknown remains
                    # fail-closed even if later IRC/endpoint work succeeds.
                    "converged": seg_converged,
                }
            )

        # --- Concatenated MEP trajectory ---
        final_trj = path_dir / "mep_trj.xyz"
        try:
            final_trj.write_text("".join(combined_blocks), encoding="utf-8")
            _echo(f"[all] Wrote concatenated MEP trajectory: {final_trj}", narrative=True)
        except Exception as e:
            raise click.ClickException(f"[all] Failed to write concatenated MEP: {e}")

        # Energy plot for concatenated trajectory
        try:
            run_trj2fig(final_trj, [path_dir / "mep_plot.png"], unit="kcal", reference="init", reverse_x=False)
            close_matplotlib_figures()
            _echo_detail(f"[plot] Saved energy plot → '{path_dir / 'mep_plot.png'}'")
        except Exception as e:
            _echo(f"[plot] WARNING: Failed to plot concatenated MEP: {e}", err=True)

        # PDB conversion of concatenated trajectory
        try:
            if pockets_for_path[0].suffix.lower() == ".pdb":
                mep_pdb_path = path_dir / "mep.pdb"
                _path_search._maybe_convert_to_pdb(
                    final_trj, ref_pdb_path=pockets_for_path[0], out_path=mep_pdb_path
                )
                if mep_pdb_path.exists():
                    shutil.copy2(mep_pdb_path, out_dir / mep_pdb_path.name)
                    _echo_detail(f"[all] Copied concatenated MEP PDB → {out_dir / mep_pdb_path.name}")
        except Exception as e:
            _echo(
                f"[all] WARNING: Failed to convert/copy concatenated MEP to PDB: {e}",
                err=True,
            )

        # --- Energy diagram ---
        energy_diagrams_po: List[Dict[str, Any]] = []
        try:
            labels = _build_global_segment_labels(len(path_opt_segments))
            energies_chain: List[float] = []
            for si, seg_info in enumerate(path_opt_segments):
                Es = [float(x) for x in seg_info.get("energies", [])]
                if not Es:
                    continue
                if si == 0:
                    energies_chain.append(Es[0])
                energies_chain.append(max(Es))
                energies_chain.append(Es[-1])
            if labels and energies_chain and len(labels) == len(energies_chain):
                title_note = (
                    f"({mep_mode_label}; all segments)"
                    if len(path_opt_segments) > 1
                    else f"({mep_mode_label})"
                )
                diag_payload = _write_segment_energy_diagram(
                    path_dir / "energy_diagram_MEP",
                    labels=labels,
                    energies_eh=energies_chain,
                    title_note=title_note,
                )
                if diag_payload:
                    energy_diagrams_po.append(diag_payload)
        except Exception as e:
            _echo(
                f"[diagram] WARNING: Failed to build {mep_mode_label} diagram "
                f"for path-opt branch: {e}",
                err=True,
            )

        # --- Bond change detection and summary.json ---
        segments_summary: List[Dict[str, Any]] = []
        bond_cfg = dict(_path_search.BOND_KW)
        for seg_idx, info in enumerate(path_opt_segments, start=1):
            Es = [float(x) for x in info.get("energies", [])]
            if not Es:
                continue
            barrier = (max(Es) - Es[0]) * AU2KCALPERMOL
            delta = (Es[-1] - Es[0]) * AU2KCALPERMOL
            bond_summary = ""
            try:
                first_last = info.get("first_last")
                if first_last:
                    elems, c_first, c_last = first_last
                else:
                    elems, c_first, c_last = read_xyz_first_last(Path(info["traj"]))
                gL = _geom_from_angstrom(elems, c_first, [])
                gR = _geom_from_angstrom(elems, c_last, [])
                changed, bond_summary = _path_search._has_bond_change(gL, gR, bond_cfg)
                if not changed:
                    bond_summary = "(no covalent changes detected)"
            except Exception as e:
                _echo(
                    f"[all] WARNING: Failed to detect bond changes for segment {seg_idx:02d}: {e}",
                    err=True,
                )
                bond_summary = "(no covalent changes detected)"

            segments_summary.append(
                {
                    "index": seg_idx,
                    "tag": info.get("tag", f"seg_{seg_idx:02d}"),
                    "kind": "seg",
                    "converged": info.get("converged"),
                    "barrier_kcal": float(barrier),
                    "delta_kcal": float(delta),
                    "bond_changes": bond_summary,
                }
            )

        po_summary: Dict[str, Any] = {
            "out_dir": str(path_dir),
            "n_images": len(read_xyz_as_blocks(final_trj)),
            "n_segments": len(segments_summary),
            "segments": segments_summary,
        }
        if energy_diagrams_po:
            po_summary["energy_diagrams"] = list(energy_diagrams_po)
        _enrich_summary(
            po_summary,
            version="",
            pipeline_mode="path-search" if refine_path else "path-opt",
            out_dir=out_dir,
            manifest=manifest,
            mlip_backend=mlip_backend_resolved,
            mlip_model=mlip_model_resolved,
            mlip_precision=mlip_precision_resolved,
            charge=q_int,
            spin=spin,
            command=command_str,
            config={
                "refine_path": bool(refine_path),
                "tsopt": do_tsopt,
                "thermo": do_thermo,
                "dft": do_dft,
                "opt_mode": tsopt_opt_mode_default,
                "path_opt_mode": path_search_opt_mode,
                "post_opt_mode": tsopt_opt_mode_default,
                "ts_opt_mode": tsopt_opt_mode_default,
                "endpoint_opt_mode": endpoint_opt_mode_default,
                "mep_mode": mep_mode_kind,
                "dmf_correlated": dmf_correlated_effective,
                "dmf_backend": dmf_backend_effective,
            },
        )
        try:
            _publish_manifest_summary(
                out_dir / "summary.json",
                po_summary,
                manifest=manifest,
                out_dir=out_dir,
                mirrors=(path_dir / "summary.json",),
            )
            _echo_detail(f"[write] Wrote '{path_dir / 'summary.json'}'.")
        except Exception as e:
            _echo(f"[write] WARNING: Failed to write summary.json for path-opt branch: {e}", err=True)

        # Copy key outputs to out_dir root
        try:
            for name in ("mep_plot.png", "energy_diagram_MEP.png"):
                src = path_dir / name
                if src.exists():
                    shutil.copy2(src, out_dir / name)
            for ext in ("_trj.xyz", ".xyz"):
                src = path_dir / f"mep{ext}"
                if src.exists():
                    shutil.copy2(src, out_dir / src.name)
        except Exception as e:
            _echo(f"[all] WARNING: Failed to relocate path-opt summary files: {e}", err=True)

    # Stage 3: Merge (performed by path_search when --ref-pdb was supplied)
    _echo_section(f"====== [all] Stage 3/{stage_total} — Core MEP outputs ======")
    _echo_detail(f"[all] Final products can be found under: {out_dir}")
    _echo_detail("  - mep_trj.xyz              (concatenated MEP trajectory)")
    _echo_detail("  - mep.pdb / mep.cif        (coordinate companions when topology is available)")
    _echo_detail("  - summary.json             (segment barriers, ΔE, bond changes)")
    _echo_detail("  - mep_plot.png / energy_diagram_MEP.png / summary.log")
    _echo_detail(f"[all] Raw per-segment MEP-engine files stay under: {path_dir}")
    _echo_detail("  - mep_seg_XX_trj.xyz       (per-segment trajectories)")
    _echo_detail("  - hei_seg_XX.xyz/.pdb      (HEI per segment)")
    _echo_section("====== [all] Core MEP pipeline finished successfully ======")

    summary_json_path = path_dir / "summary.json"
    summary_loaded = {}
    if summary_json_path.exists():
        try:
            summary_loaded = json.loads(summary_json_path.read_text(encoding="utf-8")) or {}
        except Exception:
            summary_loaded = {}
    summary: Dict[str, Any] = summary_loaded if isinstance(summary_loaded, dict) else {}
    segments = _read_summary(summary_json_path)
    energy_diagrams: List[Dict[str, Any]] = []
    existing_diagrams = summary.get("energy_diagrams", [])
    if isinstance(existing_diagrams, list):
        energy_diagrams.extend(existing_diagrams)

    def _copy_path_outputs_to_root() -> None:
        # Thin wrapper preserving the closure-captured helper signature.
        # Body extracted to mlmm.workflows._all_helpers so it is unit-
        # testable and the cli() body shrinks one slot.
        from mlmm.workflows._all_helpers import copy_path_outputs_to_root
        copy_path_outputs_to_root(
            path_dir,
            out_dir,
            warn_fn=lambda msg: _echo(msg, err=True),
        )

    def _write_pipeline_summary_log(post_segment_logs: Sequence[Dict[str, Any]]) -> None:
        # Payload assembly extracted to mlmm.workflows._all_helpers; the
        # I/O wrapper here keeps the original closure capture + error
        # routing semantics so callers do not change.
        from mlmm.workflows._all_helpers import build_pipeline_summary_payload
        nonlocal citation_post_segments
        citation_post_segments = list(post_segment_logs)
        try:
            summary_payload = build_pipeline_summary_payload(
                out_dir=out_dir,
                path_dir=path_dir,
                summary=summary,
                refine_path=refine_path,
                thresh=thresh,
                thresh_post=thresh_post,
                flatten=flatten,
                do_tsopt=do_tsopt,
                do_thermo=do_thermo,
                do_dft=do_dft,
                opt_mode_norm=opt_mode_norm,
                opt_mode_post=opt_mode_post,
                path_opt_mode=path_search_opt_mode,
                post_opt_mode=tsopt_opt_mode_default,
                ts_opt_mode=tsopt_opt_mode_default,
                endpoint_opt_mode=endpoint_opt_mode_default,
                mep_mode=mep_mode_kind,
                dmf_backend=dmf_backend_effective,
                dmf_correlated=dmf_correlated_effective,
                command_str=command_str,
                q_int=q_int,
                spin=spin,
                post_segment_logs=post_segment_logs,
                mlip_backend=mlip_backend_resolved,
                mlip_model=mlip_model_resolved,
                mlip_precision=mlip_precision_resolved,
            )
            write_summary_log(path_dir / "summary.log", summary_payload)
            _copy_path_outputs_to_root()
        except (OSError, KeyError, ValueError, TypeError) as e:
            _echo(f"[write] WARNING: Failed to write summary.log: {e}", err=True)

    # Optional Stage 4: TSOPT / THERMO / DFT (per reactive segment)
    if not (do_tsopt or do_thermo or do_dft):
        _write_pipeline_summary_log([])
        _finalize_current_summary(
            out_dir / "summary.json",
            summary,
            manifest=manifest,
            out_dir=out_dir,
            mirrors=(path_dir / "summary.json",),
        )
        # Elapsed time
        _emit_final_summary(
            out_dir,
            time_start,
            manifest,
            citation_payload=_all_method_citation_payload(),
        )
        return

    _echo_section(f"====== [all] Stage 4/{stage_total} — Post-processing per reactive segment ======")

    # Use segment summary from path_search / path-opt
    if not segments:
        _echo("[post] No segments found in summary; nothing to do.", narrative=True)
        _write_pipeline_summary_log([])
        _finalize_current_summary(
            out_dir / "summary.json",
            summary,
            manifest=manifest,
            out_dir=out_dir,
            mirrors=(path_dir / "summary.json",),
        )
        _emit_final_summary(
            out_dir,
            time_start,
            manifest,
            citation_payload=_all_method_citation_payload(),
        )
        return

    # Iterate only bond-change segments (kind='seg' and bond_changes not empty and not '(no covalent...)')
    reactive = [s for s in segments if _is_reactive_segment(s)]
    if not reactive:
        _echo("[post] No bond-change segments. Skipping TS/thermo/DFT.", narrative=True)
        _write_pipeline_summary_log([])
        _finalize_current_summary(
            out_dir / "summary.json",
            summary,
            manifest=manifest,
            out_dir=out_dir,
            mirrors=(path_dir / "summary.json",),
        )
        _emit_final_summary(
            out_dir,
            time_start,
            manifest,
            citation_payload=_all_method_citation_payload(),
        )
        return

    post_segment_logs: List[Dict[str, Any]] = []
    tsopt_seg_energies: List[Tuple[float, float, float]] = []
    g_mlip_seg_energies: List[Tuple[float, float, float]] = []
    dft_seg_energies: List[Tuple[float, float, float]] = []
    g_dft_mlip_seg_energies: List[Tuple[float, float, float]] = []
    irc_trj_for_all: List[Tuple[Path, bool]] = []

    # For each reactive segment
    for s in reactive:
        seg_idx = int(s.get("index", 0) or 0)
        seg_tag = s.get("tag", f"seg_{seg_idx:02d}")
        _echo_section(f"--- [post] seg_{seg_idx:02d} ({seg_tag}) ---")

        seg_root = path_dir  # MEP-engine scratch root (hei_seg_/mep_seg_ live here, under _work/)
        seg_dir = out_dir / SEGMENTS_DIRNAME / f"seg_{seg_idx:02d}"  # per-segment deliverables
        ensure_dir(seg_dir)

        # HEI pocket file prepared by path_search (only for bond-change segments)
        hei_pocket_pdb = seg_root / f"hei_seg_{seg_idx:02d}.pdb"
        if not hei_pocket_pdb.exists():
            _echo(f"[post] WARNING: HEI pocket PDB not found for segment {seg_idx:02d}; skipping TSOPT.", err=True)
            continue

        # 4.1 TS optimization (optional; still needed to drive IRC & diagrams)
        if do_tsopt:
            segment_tsopt_overrides = dict(tsopt_overrides)
            reference_mode_path = seg_root / f"hei_mode_seg_{seg_idx:02d}.txt"
            reference_mode_path = _ensure_hei_path_tangent(
                seg_root / f"mep_seg_{seg_idx:02d}_trj.xyz",
                hei_pocket_pdb,
                reference_mode_path,
            )
            if reference_mode_path is not None:
                segment_tsopt_overrides["reference_mode"] = reference_mode_path
            ts_pdb, g_ts = _run_tsopt_on_hei(
                hei_pocket_pdb,
                q_int,
                spin,
                real_parm7_path,
                ml_region_pdb,
                detect_layer,
                args_yaml,
                seg_dir,
                tsopt_opt_mode_default,
                resolved_calc_template=resolved_calc_template,
                overrides=segment_tsopt_overrides,
                backend=backend,
                embedcharge=embedcharge,
                embedcharge_cutoff=embedcharge_cutoff,
                embedcharge_explicit=embedcharge_explicit,
                link_atom_method=link_atom_method,
                mm_backend=mm_backend,
                use_cmap=use_cmap,
                ref_pdb=layered_inputs[0] if layered_inputs else None,
            )
        else:
            # If TSOPT is off, use the MEP highest-energy image as TS geometry.
            ts_pdb = hei_pocket_pdb
            g_ts = geom_loader(ts_pdb, coord_type="cart")
            _hei_calc_kwargs = _stage_calc_kwargs(
                resolved_calc_template,
                input_pdb=ts_pdb,
                real_parm7=real_parm7_path,
                model_pdb=ml_region_pdb,
                charge=q_int,
                spin=spin,
                use_bfactor_layers=detect_layer,
            )
            calc = _mlmm_calc(**_hei_calc_kwargs)
            g_ts.set_calculator(calc); _ = float(g_ts.energy)

        # 4.2 EulerPC IRC & mapping to (left,right)
        irc_plot_path = None
        irc_trj_path = None
        irc_res = _irc_and_match(seg_idx=seg_idx,
                                 seg_dir=seg_dir,
                                 mep_dir=path_dir,
                                 ref_pdb_for_seg=ts_pdb,
                                 seg_pocket_pdb=hei_pocket_pdb,
                                 g_ts=g_ts,
                                 q_int=q_int,
                                 spin=spin,
                                 resolved_calc_template=resolved_calc_template,
                                 real_parm7=real_parm7_path,
                                 model_pdb=ml_region_pdb,
                                 detect_layer=detect_layer,
                                 backend=backend,
                                 embedcharge=embedcharge,
                                 embedcharge_cutoff=embedcharge_cutoff,
                                 embedcharge_explicit=embedcharge_explicit,
                                 link_atom_method=link_atom_method,
                                 mm_backend=mm_backend,
                                 use_cmap=use_cmap,
                                 irc_step_size=irc_step_size,
                                 irc_never_stop=irc_never_stop,
                                 session=session,
                                 args_yaml=args_yaml)
        irc_plot_path = irc_res.get("irc_plot")
        irc_trj_path = irc_res.get("irc_trj")
        if irc_trj_path:
            try:
                irc_trj_for_all.append((Path(irc_trj_path), bool(irc_res.get("reverse_irc", False))))
            except Exception:
                logger.debug("Failed to append IRC trajectory path", exc_info=True)

        gL = irc_res["left_min_geom"]
        gR = irc_res["right_min_geom"]
        gT = irc_res["ts_geom"]
        # Save IRC endpoints (XYZ primary), run endpoint-opt, then save optimized structures
        struct_dir = seg_dir / "structures"
        ensure_dir(struct_dir)
        xL_irc, pL_irc = _save_single_geom_for_tools(gL, hei_pocket_pdb, struct_dir, "reactant_irc")
        xT, pT         = _save_single_geom_for_tools(gT, hei_pocket_pdb, struct_dir, "ts")
        xR_irc, pR_irc = _save_single_geom_for_tools(gR, hei_pocket_pdb, struct_dir, "product_irc")

        endpoint_opt_dir = seg_dir / "endpoint_opt"
        ensure_dir(endpoint_opt_dir)

        # Map IRC left/right Hessians → R/P endpoint
        # When reverse_irc is True, _irc_and_match swapped left/right to match MEP endpoints,
        # so "irc_left" (=forward) now corresponds to gR and "irc_right" (=backward) to gL.
        from mlmm.io.hessian_cache import (
            clear as _clear_hess_cache,
            discard as _hess_discard,
            load as _hess_load,
            store as _hess_store,
        )
        _reversed = bool(irc_res.get("reverse_irc", False))
        _left_hk  = "irc_right" if _reversed else "irc_left"
        _right_hk = "irc_left"  if _reversed else "irc_right"

        _hess_discard("irc_endpoint")
        _c = _hess_load(_left_hk)
        if _c:
            _hess_store("irc_endpoint", _c["hessian"], active_dofs=_c.get("active_dofs"), meta=_c.get("meta"), identity=_c.get("identity"))
        # Fail-closed endpoint-opt convergence (None if the opt could not run).
        _react_opt_conv: Optional[bool] = None
        try:
            gL, _, _react_opt_conv = _run_opt_for_state(
                pL_irc, q_int, spin, real_parm7_path, ml_region_pdb, detect_layer,
                endpoint_opt_dir / "R", args_yaml, endpoint_opt_mode_default,
                resolved_calc_template=resolved_calc_template,
                convert_files=post_convert_files_forward,
                backend=backend,
                embedcharge=embedcharge,
                embedcharge_cutoff=embedcharge_cutoff,
                embedcharge_explicit=embedcharge_explicit,
                link_atom_method=link_atom_method,
                mm_backend=mm_backend,
                use_cmap=use_cmap,
                thresh=post_thresh_forward,
                reject_uphill=_reject_uphill_eff,
                xyz_path=xL_irc,
            )
        except Exception as e:
            _echo(
                f"[post] WARNING: Reactant endpoint optimization failed for segment {seg_idx:02d}: {e}",
                err=True,
            )
            _react_opt_conv = None

        _hess_discard("irc_endpoint")
        _c = _hess_load(_right_hk)
        if _c:
            _hess_store("irc_endpoint", _c["hessian"], active_dofs=_c.get("active_dofs"), meta=_c.get("meta"), identity=_c.get("identity"))
        _prod_opt_conv: Optional[bool] = None
        try:
            gR, _, _prod_opt_conv = _run_opt_for_state(
                pR_irc, q_int, spin, real_parm7_path, ml_region_pdb, detect_layer,
                endpoint_opt_dir / "P", args_yaml, endpoint_opt_mode_default,
                resolved_calc_template=resolved_calc_template,
                convert_files=post_convert_files_forward,
                backend=backend,
                embedcharge=embedcharge,
                embedcharge_cutoff=embedcharge_cutoff,
                embedcharge_explicit=embedcharge_explicit,
                link_atom_method=link_atom_method,
                mm_backend=mm_backend,
                use_cmap=use_cmap,
                thresh=post_thresh_forward,
                reject_uphill=_reject_uphill_eff,
                xyz_path=xR_irc,
            )
        except Exception as e:
            _echo(
                f"[post] WARNING: Product endpoint optimization failed for segment {seg_idx:02d}: {e}",
                err=True,
            )
            _prod_opt_conv = None
        shutil.rmtree(endpoint_opt_dir, ignore_errors=True)
        _echo_detail("[endpoint-opt] Clean endpoint-opt working dir.")

        xL, pL = _save_single_geom_for_tools(gL, hei_pocket_pdb, struct_dir, "reactant")
        xR, pR = _save_single_geom_for_tools(gR, hei_pocket_pdb, struct_dir, "product")

        # Copy R/TS/P structures to out_dir/seg_XX/
        try:
            _state_structs = {"R": pL, "TS": pT, "P": pR}
            _input_suffix = (
                _original_input_paths[0].suffix.lower()
                if _original_input_paths
                else ".xyz"
            )
            _seg_out = _copy_structures_to_seg_dir(
                _state_structs, out_dir, seg_idx, _input_suffix,
                manifest=manifest,
            )
            _echo(f"[all] Wrote R/TS/P for segment {seg_idx:02d} → {_seg_out}", narrative=True)
        except Exception as e:
            _echo(f"[all] WARNING: Failed to copy R/TS/P structures for segment {seg_idx:02d}: {e}", err=True)

        # 4.3 Segment-level ML/MM energy diagram (R, TS, P)
        eR = float(gL.energy)
        eT = float(gT.energy)
        eP = float(gR.energy)
        tsopt_seg_energies.append((eR, eT, eP))
        mlip_prefix = seg_dir / "energy_diagram_MLIP"
        _write_public_segment_diagram(
            mlip_prefix,
            labels=["R", f"TS{seg_idx}", "P"],
            energies_eh=[eR, eT, eP],
            title_note="(MLIP, TSOPT/IRC)",
        )

        # ── Release GPU memory before freq/thermo/DFT ──
        _irc_lease = irc_res.get("calculator_lease")
        if _irc_lease is not None:
            _irc_lease.release()
        for _g in (gL, gR, gT):
            if _g is not None and hasattr(_g, "calculator"):
                _g.calculator = None
        gc.collect()
        if torch.cuda.is_available():
            torch.cuda.empty_cache()

        # 4.4 Thermochemistry (ML/MM frequencies) and Gibbs diagram
        thermo_payloads: Dict[str, Dict[str, Any]] = {}
        GR = GT = GP = None
        freq_seg_root = _resolve_override_dir(seg_dir / "freq", freq_out_dir)
        dft_seg_root = _resolve_override_dir(seg_dir / "dft", dft_out_dir)

        if do_thermo:
            _echo_detail(f"[thermo] Segment {seg_idx:02d}: freq on TS/R/P")
            tT = _run_freq_for_state(
                pT, q_int, spin, real_parm7_path, ml_region_pdb, detect_layer,
                freq_seg_root / "TS", args_yaml, overrides=freq_overrides,
                backend=backend, embedcharge=embedcharge, embedcharge_cutoff=embedcharge_cutoff,
                embedcharge_explicit=embedcharge_explicit,
                link_atom_method=link_atom_method, mm_backend=mm_backend, use_cmap=use_cmap, xyz_path=xT,
            )
            _clear_hess_cache()  # TS Hessian consumed; R/P need exact computation
            tR = _run_freq_for_state(
                pL, q_int, spin, real_parm7_path, ml_region_pdb, detect_layer,
                freq_seg_root / "R", args_yaml, overrides=freq_overrides,
                backend=backend, embedcharge=embedcharge, embedcharge_cutoff=embedcharge_cutoff,
                embedcharge_explicit=embedcharge_explicit,
                link_atom_method=link_atom_method, mm_backend=mm_backend, use_cmap=use_cmap, xyz_path=xL,
            )
            tP = _run_freq_for_state(
                pR, q_int, spin, real_parm7_path, ml_region_pdb, detect_layer,
                freq_seg_root / "P", args_yaml, overrides=freq_overrides,
                backend=backend, embedcharge=embedcharge, embedcharge_cutoff=embedcharge_cutoff,
                embedcharge_explicit=embedcharge_explicit,
                link_atom_method=link_atom_method, mm_backend=mm_backend, use_cmap=use_cmap, xyz_path=xR,
            )
            thermo_payloads = {"R": tR, "TS": tT, "P": tP}
            GR = _thermo_gibbs_ha(tR)
            GT = _thermo_gibbs_ha(tT)
            GP = _thermo_gibbs_ha(tP)
            _echo_state_energies(
                "thermo",
                seg_idx,
                "ZPE correction",
                _scale_energy_values(
                    _thermo_correction_values(thermo_payloads, "zpe_correction_ha"),
                    AU2KCALPERMOL,
                ),
                unit="kcal/mol",
                precision=2,
            )
            _echo_state_energies(
                "thermo",
                seg_idx,
                "thermal energy correction",
                _scale_energy_values(
                    _thermo_correction_values(thermo_payloads, "thermal_correction_energy_ha"),
                    AU2KCALPERMOL,
                ),
                unit="kcal/mol",
                precision=2,
            )
            _echo_state_energies(
                "thermo",
                seg_idx,
                "thermal free-energy correction",
                _scale_energy_values(
                    _thermo_correction_values(thermo_payloads, "thermal_correction_free_energy_ha"),
                    AU2KCALPERMOL,
                ),
                unit="kcal/mol",
                precision=2,
            )
            # build the per-segment MLIP Gibbs diagram ONLY when every
            # requested state returned a finite FREQ free energy. A failed/partial
            # FREQ must NOT be substituted by the MLIP electronic energy
            # (eR/eT/eP) into a Gibbs-named result.
            if all(value is not None for value in (GR, GT, GP)):
                try:
                    g_mlip_seg_energies.append((GR, GT, GP))
                    _g_rel = _relative_energy_values_kcal({"R": GR, "TS": GT, "P": GP})
                    if _g_rel is not None:
                        _echo_state_energies(
                            "thermo",
                            seg_idx,
                            "G_MLIP relative",
                            _g_rel,
                            unit="kcal/mol",
                            precision=2,
                        )
                    _write_public_segment_diagram(
                        seg_dir / "energy_diagram_G_MLIP",
                        labels=["R", f"TS{seg_idx}", "P"],
                        energies_eh=[GR, GT, GP],
                        title_note="(Gibbs, MLIP)",
                        ylabel="ΔG (kcal/mol)",
                    )
                except Exception as e:
                    _echo(f"[thermo] WARNING: failed to build Gibbs diagram: {e}", err=True)
            else:
                _echo(
                    f"[thermo] WARNING: seg {seg_idx}: one or more R/TS/P FREQ "
                    "free energies are unavailable; MLIP Gibbs diagram skipped "
                    "(no MLIP-energy "
                    "substitution).",
                    err=True,
                )

        # 4.5 DFT single-point and (optionally) DFT//MLIP/MM Gibbs
        eR_dft = eT_dft = eP_dft = None
        GR_dftMLIP = GT_dftMLIP = GP_dftMLIP = None
        if do_dft:
            _echo_detail(f"[dft] Segment {seg_idx:02d}: DFT on R/TS/P")
            dR = _run_dft_for_state(
                pL, q_int, spin, real_parm7_path, ml_region_pdb, detect_layer,
                dft_seg_root / "R", args_yaml, func_basis=dft_func_basis_use, overrides=dft_overrides,
                backend=backend, embedcharge=embedcharge, embedcharge_cutoff=embedcharge_cutoff,
                embedcharge_explicit=embedcharge_explicit,
                link_atom_method=link_atom_method, mm_backend=mm_backend, use_cmap=use_cmap, xyz_path=xL,
            )
            dT = _run_dft_for_state(
                pT, q_int, spin, real_parm7_path, ml_region_pdb, detect_layer,
                dft_seg_root / "TS", args_yaml, func_basis=dft_func_basis_use, overrides=dft_overrides,
                backend=backend, embedcharge=embedcharge, embedcharge_cutoff=embedcharge_cutoff,
                embedcharge_explicit=embedcharge_explicit,
                link_atom_method=link_atom_method, mm_backend=mm_backend, use_cmap=use_cmap, xyz_path=xT,
            )
            dP = _run_dft_for_state(
                pR, q_int, spin, real_parm7_path, ml_region_pdb, detect_layer,
                dft_seg_root / "P", args_yaml, func_basis=dft_func_basis_use, overrides=dft_overrides,
                backend=backend, embedcharge=embedcharge, embedcharge_cutoff=embedcharge_cutoff,
                embedcharge_explicit=embedcharge_explicit,
                link_atom_method=link_atom_method, mm_backend=mm_backend, use_cmap=use_cmap, xyz_path=xR,
            )
            try:
                # read the DFT energy through _dft_energy_ha, which returns
                # None when the DFT child failed (_dft_failed) — a finite hartree in
                # A failed payload must not enter a diagram.
                eR_dft = _dft_energy_ha(dR)
                eT_dft = _dft_energy_ha(dT)
                eP_dft = _dft_energy_ha(dP)
                if all(e is not None and np.isfinite(e) for e in (eR_dft, eT_dft, eP_dft)):
                    _dft_values = {"R": eR_dft, "TS": eT_dft, "P": eP_dft}
                    _echo_state_energies(
                        "dft",
                        seg_idx,
                        "E_DFT",
                        _dft_values,
                        unit="Hartree",
                        precision=8,
                    )
                    _dft_rel = _relative_energy_values_kcal(_dft_values)
                    if _dft_rel is not None:
                        _echo_state_energies(
                            "dft",
                            seg_idx,
                            "E_DFT relative",
                            _dft_rel,
                            unit="kcal/mol",
                            precision=2,
                        )
                    _dft_mlmm_values = {
                        "R": _dft_total_mlmm_energy_ha(dR),
                        "TS": _dft_total_mlmm_energy_ha(dT),
                        "P": _dft_total_mlmm_energy_ha(dP),
                    }
                    _echo_state_energies(
                        "dft",
                        seg_idx,
                        "E_total ML(dft)/MM",
                        _dft_mlmm_values,
                        unit="Hartree",
                        precision=8,
                    )
                    _dft_mlmm_rel = _relative_energy_values_kcal(_dft_mlmm_values)
                    if _dft_mlmm_rel is not None:
                        _echo_state_energies(
                            "dft",
                            seg_idx,
                            "E_total ML(dft)/MM relative",
                            _dft_mlmm_rel,
                            unit="kcal/mol",
                            precision=2,
                        )
                    dft_seg_energies.append((eR_dft, eT_dft, eP_dft))
                    _write_public_segment_diagram(
                        seg_dir / "energy_diagram_DFT",
                        labels=["R", f"TS{seg_idx}", "P"],
                        energies_eh=[eR_dft, eT_dft, eP_dft],
                        title_note=f"({dft_method_fallback})",
                    )
                else:
                    _echo("[dft] WARNING: some DFT energies missing; diagram skipped.", err=True)
            except Exception as e:
                _echo(f"[dft] WARNING: failed to build DFT diagram: {e}", err=True)

            # DFT//MLIP/MM thermal Gibbs uses the subtractive electronic total
            # plus the ML/MM thermal correction. Raw model-region DFT, MLIP, and
            # 0.0 are never substitutes for a missing component.
            if do_thermo:
                eR_dft_mlmm = _dft_total_mlmm_energy_ha(dR)
                eT_dft_mlmm = _dft_total_mlmm_energy_ha(dT)
                eP_dft_mlmm = _dft_total_mlmm_energy_ha(dP)
                dG_R = _thermo_correction_ha(tR)
                dG_T = _thermo_correction_ha(tT)
                dG_P = _thermo_correction_ha(tP)
                if all(
                    value is not None
                    for value in (
                        eR_dft_mlmm,
                        eT_dft_mlmm,
                        eP_dft_mlmm,
                        dG_R,
                        dG_T,
                        dG_P,
                    )
                ):
                    try:
                        GR_dftMLIP = eR_dft_mlmm + dG_R
                        GT_dftMLIP = eT_dft_mlmm + dG_T
                        GP_dftMLIP = eP_dft_mlmm + dG_P
                        g_dft_mlip_seg_energies.append((GR_dftMLIP, GT_dftMLIP, GP_dftMLIP))
                        _g_dft_mlip_rel = _relative_energy_values_kcal(
                            {"R": GR_dftMLIP, "TS": GT_dftMLIP, "P": GP_dftMLIP}
                        )
                        if _g_dft_mlip_rel is not None:
                            _echo_state_energies(
                                "dft//mlip",
                                seg_idx,
                                "G_DFT+thermo relative",
                                _g_dft_mlip_rel,
                                unit="kcal/mol",
                                precision=2,
                            )
                        _write_public_segment_diagram(
                            seg_dir / "energy_diagram_G_DFT_plus_MLIP",
                            labels=["R", f"TS{seg_idx}", "P"],
                            energies_eh=[GR_dftMLIP, GT_dftMLIP, GP_dftMLIP],
                            title_note="(Gibbs, DFT//MLIP/MM)",
                            ylabel="ΔG (kcal/mol)",
                        )
                    except Exception as e:
                        _echo(f"[dft//mlip] WARNING: failed to build DFT//MLIP/MM Gibbs diagram: {e}", err=True)
                else:
                    _echo(
                        f"[dft//mlip] WARNING: seg {seg_idx}: a subtractive "
                        "DFT//MLIP/MM electronic total or FREQ thermal correction "
                        "is unusable; Gibbs diagram skipped (no raw-model DFT, "
                        "MLIP, or 0.0 substitution).",
                        err=True,
                    )

        segment_log: Dict[str, Any] = {
            "index": seg_idx,
            "tag": seg_tag,
            "kind": s.get("kind", "seg"),
            "bond_changes": s.get("bond_changes", ""),
            "mep_barrier_kcal": s.get("barrier_kcal"),
            "mep_delta_kcal": s.get("delta_kcal"),
            "post_dir": str(seg_dir),
        }
        if irc_plot_path:
            segment_log["irc_plot"] = str(irc_plot_path)
        if irc_trj_path:
            segment_log["irc_traj"] = str(irc_trj_path)
        # thread the per-direction IRC outcome so the aggregate
        # gates on convergence, not trajectory-file existence.
        _irc_outcome_seg = irc_res.get("irc_outcome")
        if isinstance(_irc_outcome_seg, dict):
            segment_log["irc"] = _irc_outcome_seg
        segment_log["endpoint_assignment"] = irc_res.get(
            "endpoint_assignment"
        )
        # record endpoint-opt convergence so a nonconverged endpoint
        # (whose geometry is still used for the diagram) does not silently
        # promote its segment to a usable success.
        segment_log["endpoint_opt"] = {
            "reactant_converged": _react_opt_conv,
            "product_converged": _prod_opt_conv,
        }
        if do_tsopt:
            tsopt_n_imag = (getattr(gT, "_tsopt_result", {}) or {}).get(
                "n_imaginary_modes"
            )
            if tsopt_n_imag is not None:
                segment_log["ts_imag"] = _ts_imag_record(
                    tsopt_n_imag,
                    (getattr(gT, "_tsopt_result", {}) or {}).get(
                        "imaginary_frequencies_cm"
                    ),
                )
        if do_thermo:
            n_imag = None
            try:
                n_imag = int(thermo_payloads.get("TS", {}).get("num_imag_freq"))
            except Exception:
                n_imag = None
            if n_imag is not None:
                # thermoanalysis.yaml carries `num_imag_freq` but no frequency
                # list, so this branch must not clobber the frequencies the
                # tsopt branch already published above.
                _prior_freqs = (segment_log.get("ts_imag") or {}).get(
                    "imag_freqs_cm"
                )
                segment_log["ts_imag"] = _ts_imag_record(
                    n_imag,
                    (thermo_payloads.get("TS") or {}).get(
                        "imaginary_frequencies_cm"
                    )
                    or _prior_freqs,
                )
        from mlmm.workflows._all_helpers import (
            build_energy_level_dict,
            build_thermo_symmetry_provenance,
        )
        _thermo_symmetry = build_thermo_symmetry_provenance(thermo_payloads)
        if _thermo_symmetry:
            segment_log["thermo_symmetry"] = _thermo_symmetry
        _structs_seg = {"R": pL, "TS": pT, "P": pR}
        segment_log["mlip"] = build_energy_level_dict(
            labels=["R", "TS", "P"],
            energies_au=[eR, eT, eP],
            ref_energy=eR,
            au_to_kcal=AU2KCALPERMOL,
            diagram_path=str((seg_dir / "energy_diagram_MLIP").with_suffix(".png")),
            structures=_structs_seg,
        )
        if GR is not None and GT is not None and GP is not None:
            segment_log["gibbs_mlip"] = build_energy_level_dict(
                labels=["R", "TS", "P"],
                energies_au=[GR, GT, GP],
                ref_energy=GR,
                au_to_kcal=AU2KCALPERMOL,
                diagram_path=str((seg_dir / "energy_diagram_G_MLIP").with_suffix(".png")),
                structures=_structs_seg,
            )
        if eR_dft is not None and eT_dft is not None and eP_dft is not None and all(
            map(np.isfinite, [eR_dft, eT_dft, eP_dft])
        ):
            segment_log["dft"] = build_energy_level_dict(
                labels=["R", "TS", "P"],
                energies_au=[eR_dft, eT_dft, eP_dft],
                ref_energy=eR_dft,
                au_to_kcal=AU2KCALPERMOL,
                diagram_path=str((seg_dir / "energy_diagram_DFT").with_suffix(".png")),
                structures=_structs_seg,
            )
        if GR_dftMLIP is not None and GT_dftMLIP is not None and GP_dftMLIP is not None:
            segment_log["gibbs_dft_mlip"] = build_energy_level_dict(
                labels=["R", "TS", "P"],
                energies_au=[GR_dftMLIP, GT_dftMLIP, GP_dftMLIP],
                ref_energy=GR_dftMLIP,
                au_to_kcal=AU2KCALPERMOL,
                diagram_path=str((seg_dir / "energy_diagram_G_DFT_plus_MLIP").with_suffix(".png")),
                structures=_structs_seg,
            )

        post_segment_logs.append(segment_log)

    _all_diagram_specs = [
        (True, tsopt_seg_energies, "energy_diagram_MLIP_all",
         "(MLIP, TSOPT + IRC; all segments)", None),
        (do_thermo, g_mlip_seg_energies, "energy_diagram_G_MLIP_all",
         "(MLIP + Thermal Correction; all segments)", "ΔG (kcal/mol)"),
        (do_dft, dft_seg_energies, "energy_diagram_DFT_all",
         f"({dft_method_fallback}; all segments)", None),
        (do_dft and do_thermo, g_dft_mlip_seg_energies, "energy_diagram_G_DFT_plus_MLIP_all",
         f"({dft_method_fallback} // MLIP + Thermal Correction; all segments)", "ΔG (kcal/mol)"),
    ]
    from mlmm.workflows._all_helpers import has_complete_segment_energy_series

    for cond, seg_energies, fname_stem, title_note, ylabel in _all_diagram_specs:
        if not cond or not has_complete_segment_energy_series(
            seg_energies, expected_segments=len(reactive)
        ):
            continue
        all_energies = [e for triple in seg_energies for e in triple]
        all_labels = _build_global_segment_labels(len(seg_energies))
        if not (all_labels and len(all_labels) == len(all_energies)):
            continue
        extra_kwargs = {"ylabel": ylabel} if ylabel is not None else {}
        diag_payload = _write_segment_energy_diagram(
            out_dir / fname_stem,
            labels=all_labels,
            energies_eh=all_energies,
            title_note=title_note,
            write_html=False,
            **extra_kwargs,
        )
        if diag_payload:
            energy_diagrams.append(diag_payload)

    if irc_trj_for_all:
        _merge_irc_trajectories_to_single_plot(
            irc_trj_for_all, out_dir / "irc_plot_all.png"
        )

    # Refresh summary.json with final energy diagram metadata
    try:
        summary["energy_diagrams"] = list(energy_diagrams)
        _enrich_summary(
            summary,
            version="",
            pipeline_mode="path-search" if refine_path else "path-opt",
            out_dir=out_dir,
            manifest=manifest,
            mlip_backend=mlip_backend_resolved,
            mlip_model=mlip_model_resolved,
            mlip_precision=mlip_precision_resolved,
            charge=q_int,
            spin=spin,
            command=command_str,
            post_segments=post_segment_logs,
            config={
                "refine_path": bool(refine_path),
                "tsopt": do_tsopt,
                "thermo": do_thermo,
                "dft": do_dft,
                "opt_mode": tsopt_opt_mode_default,
                "path_opt_mode": path_search_opt_mode,
                "post_opt_mode": tsopt_opt_mode_default,
                "ts_opt_mode": tsopt_opt_mode_default,
                "endpoint_opt_mode": endpoint_opt_mode_default,
                "mep_mode": mep_mode_kind,
                "dmf_correlated": dmf_correlated_effective,
                "dmf_backend": dmf_backend_effective,
            },
        )
        _publish_manifest_summary(
            out_dir / "summary.json",
            summary,
            manifest=manifest,
            out_dir=out_dir,
            mirrors=(path_dir / "summary.json",),
        )
    except Exception as e:
        _echo(f"[write] WARNING: Failed to refresh summary.json with energy diagram metadata: {e}", err=True)

    _write_pipeline_summary_log(post_segment_logs)
    _finalize_current_summary(
        out_dir / "summary.json",
        summary,
        manifest=manifest,
        out_dir=out_dir,
        mirrors=(path_dir / "summary.json",),
    )
    _emit_final_summary(
        out_dir,
        time_start,
        manifest,
        citation_payload=_all_method_citation_payload(),
    )


_configure_all_help_visibility(cli)


if __name__ == "__main__":
    cli()
