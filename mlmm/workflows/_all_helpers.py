"""Helpers for the `mlmm all` subcommand.

Holds side-effect-free helpers for `mlmm/workflows/all.py:cli()`.

Anything imported here must be safe to use from ``cli()`` body callers
without changing observable behavior.
"""

from __future__ import annotations

import shutil
from pathlib import Path
from typing import Any, Callable, Collection, Dict, Mapping, Optional, Sequence, Tuple


def copy_path_outputs_to_root(
    path_dir: Path,
    out_dir: Path,
    *,
    warn_fn: Optional[Callable[[str], None]] = None,
) -> None:
    """Copy MEP path-search outputs from ``path_dir`` up to ``out_dir``.

    Extracted from the nested ``_copy_path_outputs_to_root`` in
    ``workflows/all.py:cli()``. Behavior-equivalent: best-effort copy
    of the canonical MEP artifacts (plot / pdb / xyz / summary), with
    any copy failure routed to ``warn_fn`` instead of propagating.

    When ``warn_fn`` is None (default), failures are swallowed silently
    — the original nested helper printed a "[all] WARNING: ..." line,
    so callers that want that behavior pass an `_echo`-equivalent.
    """
    try:
        from mlmm.core.result_commit import commit_exact_bytes

        # MEP deliverables: MOVE to root (no _work/ duplicate). The summary
        # payload reads these from path_dir before this runs, so the move is
        # safe; summary.{json,log} stay COPY because they are re-written under
        # path_dir downstream.
        for name in (
            "mep_plot.png",
            "energy_diagram_MEP.png",
            "mep.pdb",
            "mep.cif",
        ):
            src = path_dir / name
            if src.exists():
                shutil.move(str(src), str(out_dir / name))
        for stem in ("mep",):
            for ext in ("_trj.xyz", ".xyz"):
                src = path_dir / f"{stem}{ext}"
                if src.exists():
                    shutil.move(str(src), str(out_dir / src.name))
        for name in ("summary.json", "summary.log"):
            src = path_dir / name
            if src.exists():
                commit_exact_bytes(out_dir / name, src.read_bytes())
    except (OSError, shutil.Error) as exc:
        if warn_fn is not None:
            warn_fn(f"[all] WARNING: Failed to copy path_search outputs: {exc}")


def build_energy_level_dict(
    *,
    labels: Sequence[str],
    energies_au: Sequence[float],
    ref_energy: float,
    au_to_kcal: float,
    diagram_path: str,
    structures: Dict[str, Any],
) -> Dict[str, Any]:
    """Assemble one ML/MM, Gibbs, or DFT R/TS/P energy payload.

    The shared layout contains absolute and relative energies, barrier and
    reaction energies, the diagram path, and R/TS/P structure mappings.

    Parameters
    ----------
    labels : sequence of str
        Per-image labels (typically ``["R", "TS", "P"]``).
    energies_au : sequence of float
        Energies (Hartree) for each label, same order as ``labels``.
    ref_energy : float
        Reference energy (Hartree) — typically the first entry of
        ``energies_au``. Used as the kcal/mol zero point.
    au_to_kcal : float
        ``pysisyphus.constants.AU2KCALPERMOL`` multiplier; passed in
        so this helper doesn't reach back into the pysisyphus import.
    diagram_path : str
        Path to the rendered energy diagram PNG.
    structures : dict
        ``{label: pdb_path}`` mapping for downstream consumers.
    """
    kcal = [(e - ref_energy) * au_to_kcal for e in energies_au]
    payload = {
        "labels": list(labels),
        "energies_au": list(energies_au),
        "energies_kcal": kcal,
        "diagram": diagram_path,
        "structures": dict(structures),
    }
    if list(labels) == ["R", "TS", "P"]:
        payload["barrier_kcal"] = kcal[1]
        payload["delta_kcal"] = kcal[-1]
    elif list(labels) == ["E1", "TS", "E2"]:
        payload["barrier_from_endpoint_1_kcal"] = kcal[1] - kcal[0]
        payload["barrier_from_endpoint_2_kcal"] = kcal[1] - kcal[2]
    return payload


def promote_diag_for_root(
    diag: Optional[Dict[str, Any]],
    stem: str,
    out_dir: Path,
) -> Optional[Dict[str, Any]]:
    """Re-tag an energy-diagram dict so it points at the ``_all`` rendering.

    Extracted from the nested ``_promote_all`` in
    ``workflows/all.py:cli()``. Pure function: returns a NEW dict with
    `name` rewritten to ``f"{stem}_all"`` and `image` pointing at
    ``out_dir / f"{stem}_all.png"``, leaving the caller to append it
    into the energy_diagrams list. Returns ``None`` when ``diag`` is
    empty / None so the caller can `if promoted := ...:` skip cleanly.
    """
    if not diag:
        return None
    promoted = dict(diag)
    promoted["name"] = f"{stem}_all"
    promoted["image"] = str(out_dir / f"{stem}_all.png")
    return promoted


def build_pipeline_summary_payload(
    *,
    out_dir: Path,
    path_dir: Path,
    summary: Dict[str, Any],
    refine_path: bool,
    thresh: Optional[str],
    thresh_post: str,
    flatten: bool,
    do_tsopt: bool,
    do_thermo: bool,
    do_dft: bool,
    opt_mode_norm: str,
    opt_mode_post: Optional[str],
    path_opt_mode: Optional[str],
    post_opt_mode: Optional[str],
    ts_opt_mode: Optional[str],
    endpoint_opt_mode: Optional[str],
    mep_mode: str,
    dmf_backend: str,
    dmf_correlated: bool,
    command_str: str,
    q_int: int,
    spin: int,
    post_segment_logs: Sequence[Dict[str, Any]],
    mlip_backend: str = "uma",
    mlip_model: Optional[str] = None,
    mlip_precision: Optional[str] = None,
) -> Dict[str, Any]:
    """Assemble the summary_log payload for the `all` pipeline.

    Extracted from the inner body of the nested
    ``_write_pipeline_summary_log`` in ``workflows/all.py:cli()``;
    splitting the dict construction out makes it easy to unit-test
    without exercising the surrounding I/O.

    The returned dict matches the shape consumed by
    ``mlmm.io.summary.write_summary_log``.
    """
    diag_for_log: Dict[str, Any] = {}
    for diag in summary.get("energy_diagrams", []) or []:
        if isinstance(diag, dict) and str(diag.get("name", "")).lower().endswith("mep"):
            diag_for_log = diag
            break
    mep_info = {
        "n_images": summary.get("n_images"),
        "n_segments": summary.get("n_segments"),
        "traj_pdb": str(out_dir / "mep.pdb") if (path_dir / "mep.pdb").exists() else None,
        "mep_plot": str(out_dir / "mep_plot.png") if (path_dir / "mep_plot.png").exists() else None,
        "diagram": diag_for_log,
    }
    return {
        "root_out_dir": str(out_dir),
        "path_dir": str(path_dir),
        "path_module_dir": path_dir.name,
        "pipeline_mode": "path-search" if refine_path else "path-opt",
        "refine_path": bool(refine_path),
        "thresh": thresh,
        "thresh_post": thresh_post,
        "flatten": bool(flatten),
        "tsopt": do_tsopt,
        "thermo": do_thermo,
        "dft": do_dft,
        "opt_mode": opt_mode_norm,
        "opt_mode_post": opt_mode_post.lower() if opt_mode_post else None,
        "path_opt_mode": (
            path_opt_mode.lower() if path_opt_mode else None
        ),
        "post_opt_mode": (
            post_opt_mode.lower() if post_opt_mode else None
        ),
        "ts_opt_mode": ts_opt_mode.lower() if ts_opt_mode else None,
        "endpoint_opt_mode": (
            endpoint_opt_mode.lower() if endpoint_opt_mode else None
        ),
        "mep_mode": mep_mode,
        "dmf_backend": dmf_backend,
        "dmf_correlated": bool(dmf_correlated),
        "mlip_backend": mlip_backend,
        "mlip_model": mlip_model,
        "mlip_precision": mlip_precision,
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
        "mep": mep_info,
        "segments": summary.get("segments", []),
        "energy_diagrams": summary.get("energy_diagrams", []),
        "post_segments": list(post_segment_logs),
        "key_files": {},
    }


def build_tsopt_overrides(
    *,
    tsopt_max_cycles: Optional[int],
    dump: bool,
    dump_override_requested: bool,
    tsopt_out_dir: Optional[Path],
    hessian_calc_mode: Optional[str],
    opt_mode_post_norm: Optional[str],
    opt_mode_post_set: bool,
    opt_mode_set: bool,
    tsopt_opt_mode_default: Optional[str],
    convert_files: bool,
    convert_files_explicit: bool,
    thresh_post_forward: Optional[str],
    flatten_explicit: bool,
    flatten: Optional[bool],
    skip_final_freq: bool,
    skip_final_freq_explicit: bool,
    stop_plateau: Optional[bool] = None,
    stop_plateau_thresh: Optional[float] = None,
    stop_plateau_window: Optional[int] = None,
) -> Dict[str, Any]:
    """Assemble the `tsopt_overrides` dict consumed by the post-MEP TSOPT call.

    Extracted from the inline `if x is not None: tsopt_overrides[k] = ...`
    ladder in ``workflows/all.py:cli()`` (~20 LOC) so the dict-construction
    is independently unit-testable.  Default-valued parent options are
    omitted unless their parameter source is explicit or the pipeline has
    deliberately resolved a value for the child stage.
    """
    overrides: Dict[str, Any] = {}
    if tsopt_max_cycles is not None:
        overrides["max_cycles"] = int(tsopt_max_cycles)
    if dump_override_requested:
        overrides["dump"] = bool(dump)
    if tsopt_out_dir is not None:
        overrides["out_dir"] = tsopt_out_dir
    if hessian_calc_mode is not None:
        overrides["hessian_calc_mode"] = hessian_calc_mode
    if opt_mode_post_set and opt_mode_post_norm in {"grad", "hess"}:
        overrides["opt_mode"] = opt_mode_post_norm
    elif opt_mode_set:
        overrides["opt_mode"] = tsopt_opt_mode_default
    if convert_files_explicit:
        overrides["convert_files"] = bool(convert_files)
    if thresh_post_forward is not None:
        overrides["thresh"] = str(thresh_post_forward)
    if flatten_explicit:
        overrides["flatten"] = bool(flatten)
    if skip_final_freq_explicit:
        overrides["skip_final_freq"] = bool(skip_final_freq)
    if stop_plateau is not None:
        overrides["stop_plateau"] = bool(stop_plateau)
    if stop_plateau_thresh is not None:
        overrides["stop_plateau_thresh"] = float(stop_plateau_thresh)
    if stop_plateau_window is not None:
        overrides["stop_plateau_window"] = int(stop_plateau_window)
    return overrides


def build_freq_overrides(
    *,
    freq_max_write: Optional[int],
    freq_amplitude_ang: Optional[float],
    freq_n_frames: Optional[int],
    freq_sort: Optional[str],
    freq_temperature: Optional[float],
    freq_pressure: Optional[float],
    dump_override_requested: bool,
    dump: bool,
    require_thermo_artifact: bool,
    hessian_calc_mode: Optional[str],
    convert_files: bool,
    convert_files_explicit: bool,
) -> Dict[str, Any]:
    """Assemble the `freq_overrides` dict for the post-TSOPT FREQ call.

    Mirror of :func:`build_tsopt_overrides` for the freq stage.
    """
    overrides: Dict[str, Any] = {}
    if freq_max_write is not None:
        overrides["max_write"] = int(freq_max_write)
    if freq_amplitude_ang is not None:
        overrides["amplitude_ang"] = float(freq_amplitude_ang)
    if freq_n_frames is not None:
        overrides["n_frames"] = int(freq_n_frames)
    if freq_sort is not None:
        overrides["sort"] = freq_sort.lower()
    if freq_temperature is not None:
        overrides["temperature"] = float(freq_temperature)
    if freq_pressure is not None:
        overrides["pressure"] = float(freq_pressure)
    if require_thermo_artifact:
        # ``all --thermo`` consumes this child artifact as its stage hand-off.
        # Dump is forced to True to ensure the artifact is produced.
        overrides["dump"] = True
    elif dump_override_requested:
        overrides["dump"] = bool(dump)
    if hessian_calc_mode is not None:
        overrides["hessian_calc_mode"] = hessian_calc_mode
    if convert_files_explicit:
        overrides["convert_files"] = bool(convert_files)
    return overrides


def build_thermo_symmetry_provenance(
    thermo_payloads: Mapping[str, Mapping[str, Any]],
) -> Dict[str, Dict[str, Any]]:
    """Copy complete, child-reported state symmetry provenance."""

    provenance: Dict[str, Dict[str, Any]] = {}
    for label in thermo_payloads:
        payload = thermo_payloads.get(label)
        if not isinstance(payload, Mapping):
            continue
        value = payload.get("symmetry_number")
        source = payload.get("symmetry_number_source")
        if (
            isinstance(value, int)
            and not isinstance(value, bool)
            and value > 0
            and isinstance(source, str)
            and source.strip()
        ):
            state_provenance = {
                "symmetry_number": value,
                "symmetry_number_source": source,
            }
            point_group = payload.get("point_group")
            point_group_source = payload.get("point_group_source")
            if (
                isinstance(point_group, str)
                and point_group.strip()
                and isinstance(point_group_source, str)
                and point_group_source.strip()
            ):
                state_provenance.update(
                    point_group=point_group,
                    point_group_source=point_group_source,
                )
            provenance[label] = state_provenance
    return provenance


def has_complete_segment_energy_series(
    segment_energies: Sequence[Sequence[Any]],
    *,
    expected_segments: int,
) -> bool:
    """Return whether every segment has exactly one R/TS/P energy series."""

    return (
        expected_segments > 0
        and len(segment_energies) == expected_segments
        and all(len(values) == 3 for values in segment_energies)
    )


ChildArgSpec = Tuple[str, str, Any, bool]


def build_explicit_child_argv(
    explicit_params: Collection[str],
    specs: Sequence[ChildArgSpec],
) -> list[str]:
    """Return canonical child tokens for explicitly supplied parent options.

    Each spec is ``(parameter_name, option, value, is_toggle)``.  Omitted
    parent defaults produce no token, so the child's effective YAML remains
    authoritative.  Pipeline-owned paths, charge, topology, and other
    deliberately resolved values are assembled separately by the caller.
    """

    argv: list[str] = []
    for parameter_name, option, value, is_toggle in specs:
        if parameter_name not in explicit_params or value is None:
            continue
        if is_toggle:
            if not isinstance(value, bool):
                raise TypeError(
                    f"Toggle option {option!r} requires bool, got {type(value).__name__}."
                )
            positive = option if not option.startswith("--no-") else f"--{option[5:]}"
            negative = f"--no-{positive[2:]}"
            argv.append(positive if value else negative)
        else:
            argv.extend([option, str(value)])
    return argv


def build_path_child_argv(
    explicit_params: Collection[str],
    *,
    mep_mode: str,
    dmf_backend: str,
    max_nodes: int,
    max_cycles_gsm: Optional[int],
    max_cycles_dmf: Optional[int],
    climb: bool,
    dump: bool,
    pre_opt: bool,
    convert_files: bool,
    thresh: Optional[str],
    thresh_gsm: Optional[str] = None,
    thresh_dmf: Optional[str] = None,
) -> list[str]:
    """Build parent-controlled argv shared by path-search and path-opt children.

    The selected MEP algorithm is always forwarded because it is an ``all``
    workflow selector rather than a YAML setting.  Other parent defaults stay
    absent so child YAML remains authoritative; in particular, the DMF backend
    is forwarded only when explicitly supplied. Pipeline-owned input, charge,
    topology, output, and config tokens remain at the dispatch call site.
    """

    mode = str(mep_mode).strip().lower()
    cycle_parameter = "max_cycles_dmf" if mode == "dmf" else "max_cycles_gsm"
    cycle_value = max_cycles_dmf if mode == "dmf" else max_cycles_gsm

    specs: list[ChildArgSpec] = [
        ("dmf_backend", "--dmf-backend", dmf_backend, False),
        ("max_nodes", "--max-nodes", max_nodes, False),
        ("climb", "--climb", climb, True),
    ]
    specs.extend(
        [
            ("dump", "--dump", dump, True),
            ("pre_opt", "--preopt", pre_opt, True),
            ("convert_files", "--convert-files", convert_files, True),
            ("thresh", "--thresh", thresh, False),
            ("thresh_gsm", "--thresh-gsm", thresh_gsm, False),
            ("thresh_dmf", "--thresh-dmf", thresh_dmf, False),
        ]
    )
    argv = ["--mep-mode", mode]
    if cycle_value is not None and cycle_parameter in explicit_params:
        # Forward only the budget for the selected MEP algorithm; an explicit
        # GSM-only budget must never cap a DMF child (and vice versa).
        argv.extend([f"--max-cycles-{mode}", str(cycle_value)])
    argv.extend(build_explicit_child_argv(explicit_params, specs))
    return argv


def build_scan_child_argv(
    explicit_params: Collection[str],
    *,
    convert_files: bool,
    thresh: Optional[str],
) -> list[str]:
    """Build the explicit-only parent options forwarded to staged scan."""

    return build_explicit_child_argv(
        explicit_params,
        (
            ("convert_files", "--convert-files", convert_files, True),
            ("thresh", "--thresh", thresh, False),
        ),
    )


def resolve_post_thresh_forwarding(
    explicit_params: Collection[str],
    *,
    thresh_post: Optional[str],
    yaml_cfg: Mapping[str, Any],
) -> Optional[str]:
    """Return a post-stage threshold only when it should become a CLI token."""

    opt_cfg = yaml_cfg.get("opt")
    if "thresh_post" in explicit_params:
        return None if thresh_post is None else str(thresh_post)
    if isinstance(opt_cfg, Mapping) and "thresh" in opt_cfg:
        return None
    return None if thresh_post is None else str(thresh_post)


def resolve_dft_func_basis_forwarding(
    explicit_params: Collection[str],
    *,
    dft_func_basis: Optional[str],
    yaml_cfg: Mapping[str, Any],
    default: str = "wb97m-v/def2-tzvpd",
) -> Tuple[Optional[str], str]:
    """Return the child CLI value and effective DFT method label.

    A YAML-only method is used for reporting but remains absent from child
    argv, allowing the DFT command to resolve its own effective config.
    """

    forwarded = (
        str(dft_func_basis)
        if "dft_func_basis" in explicit_params and dft_func_basis is not None
        else None
    )
    dft_cfg = yaml_cfg.get("dft")
    yaml_value = dft_cfg.get("func_basis") if isinstance(dft_cfg, Mapping) else None
    effective = str(forwarded or yaml_value or default)
    return forwarded, effective


def append_backend_forwarding_args(
    args: list,
    *,
    backend: Optional[str],
    embedcharge: bool,
    embedcharge_cutoff: Optional[float],
    embedcharge_explicit: bool,
    link_atom_method: Optional[str],
    mm_backend: Optional[str],
    use_cmap: Optional[bool],
    args_yaml: Optional[Any] = None,
) -> None:
    """Append shared calculator options to a child CLI argv list.

    Mutates ``args`` in place; matches the existing ``_append_cli_arg``
    / ``_append_toggle_arg`` style.
    """
    if backend is not None:
        args.extend(["--backend", str(backend)])
    if embedcharge:
        args.append("--embedcharge")
        if embedcharge_cutoff is not None:
            args.extend(["--embedcharge-cutoff", str(embedcharge_cutoff)])
    elif embedcharge_explicit:
        args.append("--no-embedcharge")
    if link_atom_method is not None:
        args.extend(["--link-atom-method", str(link_atom_method)])
    if mm_backend is not None:
        args.extend(["--mm-backend", str(mm_backend)])
    if use_cmap is not None:
        args.extend(["--cmap" if use_cmap else "--no-cmap"])
    if args_yaml is not None:
        args.extend(["--config", str(args_yaml)])


def build_dft_overrides(
    *,
    dft_max_cycle: Optional[int],
    dft_conv_tol: Optional[float],
    dft_grid_level: Optional[int],
    dft_engine: Optional[str],
    dft_func_basis_forward: Optional[str],
    convert_files: bool,
    convert_files_explicit: bool,
) -> Dict[str, Any]:
    """Assemble the `dft_overrides` dict for the post-FREQ DFT call.

    Mirror of :func:`build_tsopt_overrides` for the dft stage.
    """
    overrides: Dict[str, Any] = {}
    if dft_max_cycle is not None:
        overrides["max_cycle"] = int(dft_max_cycle)
    if dft_conv_tol is not None:
        overrides["conv_tol"] = float(dft_conv_tol)
    if dft_grid_level is not None:
        overrides["grid_level"] = int(dft_grid_level)
    if dft_engine is not None:
        overrides["engine"] = str(dft_engine)
    if dft_func_basis_forward is not None:
        overrides["func_basis"] = str(dft_func_basis_forward)
    if convert_files_explicit:
        overrides["convert_files"] = bool(convert_files)
    return overrides


__all__ = [
    "copy_path_outputs_to_root",
    "promote_diag_for_root",
    "build_energy_level_dict",
    "build_pipeline_summary_payload",
    "build_tsopt_overrides",
    "build_freq_overrides",
    "has_complete_segment_energy_series",
    "build_thermo_symmetry_provenance",
    "build_dft_overrides",
    "build_explicit_child_argv",
    "build_path_child_argv",
    "build_scan_child_argv",
    "resolve_post_thresh_forwarding",
    "resolve_dft_func_basis_forwarding",
    "append_backend_forwarding_args",
]
