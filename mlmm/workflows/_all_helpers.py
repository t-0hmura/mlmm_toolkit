"""Helpers for the `mlmm all` subcommand.

Holds side-effect-free helpers plus a frozen parameter context
(``AllContext``) for `mlmm/workflows/all.py:cli()`.

Anything imported here must be safe to use from ``cli()`` body callers
without changing observable behavior.
"""

from __future__ import annotations

import shutil
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Callable, Dict, Optional, Sequence


@dataclass(frozen=True)
class AllContext:
    """Frozen bundle of the `mlmm all` CLI parameters.

    Mirrors the ``cli()`` signature in `mlmm/workflows/all.py` so that
    helper functions can accept a single argument instead of re-listing
    74 keyword parameters.

    A drift guard in ``tests/test_all_helpers.py`` keeps the field names
    in lockstep with the ``cli.callback`` signature.

    Field order matches the CLI option declaration order in
    `workflows/all.py` so reading the two side-by-side is easy.
    """

    # Inputs / output
    input_paths: Sequence[Path]
    center_spec: Optional[str]
    out_dir: Path
    # Extract knobs
    radius: float
    radius_het2het: float
    include_h2o: bool
    exclude_backbone: bool
    add_linkh: bool
    selected_resn: str
    modified_residue: str
    ligand_charge: Optional[str]
    charge_override: Optional[int]
    parm7_override: Optional[Path]
    model_pdb_override: Optional[Path]
    # MM topology
    mm_ff_set: str
    mm_add_ter: bool
    mm_keep_temp: bool
    mm_ligand_mult: Optional[str]
    # Charges
    spin: int
    # MEP search
    max_nodes: int
    max_cycles: int
    climb: bool
    opt_mode: str
    opt_mode_post: Optional[str]
    dump: bool
    refine_path: bool
    thresh: Optional[str]
    thresh_post: str
    config_yaml: Optional[Path]
    show_config: bool
    dry_run: bool
    pre_opt: bool
    hessian_calc_mode: Optional[str]
    detect_layer: bool
    # Stage toggles
    do_tsopt: bool
    do_thermo: bool
    do_dft: bool
    # Scan
    scan_lists_raw: Sequence[str]
    scan_out_dir: Optional[Path]
    scan_one_based: Optional[bool]
    scan_max_step_size: Optional[float]
    scan_bias_k: Optional[float]
    scan_relax_max_cycles: Optional[int]
    scan_preopt_override: Optional[bool]
    scan_endopt_override: Optional[bool]
    convert_files: bool
    ref_pdb_cli: Optional[Path]
    # MLIP backend
    backend: Optional[str]
    embedcharge: bool
    embedcharge_cutoff: Optional[float]
    link_atom_method: Optional[str]
    mm_backend: Optional[str]
    use_cmap: Optional[bool]
    # TSOPT / FREQ / DFT
    tsopt_max_cycles: Optional[int]
    flatten: bool
    skip_final_freq: bool
    tsopt_out_dir: Optional[Path]
    freq_out_dir: Optional[Path]
    freq_max_write: Optional[int]
    freq_amplitude_ang: Optional[float]
    freq_n_frames: Optional[int]
    freq_sort: Optional[str]
    freq_temperature: Optional[float]
    freq_pressure: Optional[float]
    dft_out_dir: Optional[Path]
    dft_func_basis: Optional[str]
    dft_max_cycle: Optional[int]
    dft_conv_tol: Optional[float]
    dft_grid_level: Optional[int]
    dft_engine: Optional[str]
    cli_coord_type: Optional[str]
    precision: Optional[str]
    backend_model: Optional[str]


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
        # MEP deliverables: MOVE to root (no _work/ duplicate). The summary
        # payload reads these from path_dir before this runs, so the move is
        # safe; summary.{json,log} stay COPY because they are re-written under
        # path_dir downstream.
        for name in ("mep_plot.png", "energy_diagram_MEP.png", "mep.pdb"):
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
                shutil.copy2(src, out_dir / name)
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
    """Assemble a segment_log sub-dict for one energy level (UMA / Gibbs / DFT).

    Extracted from the repeated R/TS/P payload pattern in
    ``workflows/all.py:cli()`` (TSOPT-only mode segment_log). The 3
    inline copies that built ``segment_log["uma"]`` /
    ``segment_log["gibbs_uma"]`` / ``segment_log["dft"]`` all shared
    the same layout: energies in Hartree, energies in kcal/mol relative
    to the reactant, barrier_kcal, delta_kcal, plus the diagram + the
    R/TS/P structure mapping. Pulling them through this helper drops
    the in-place duplication.

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
    return {
        "labels": list(labels),
        "energies_au": list(energies_au),
        "energies_kcal": kcal,
        "diagram": diagram_path,
        "structures": dict(structures),
        "barrier_kcal": kcal[1] if len(kcal) > 1 else 0.0,
        "delta_kcal": kcal[-1] if kcal else 0.0,
    }


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
    command_str: str,
    q_int: int,
    spin: int,
    post_segment_logs: Sequence[Dict[str, Any]],
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
        "mep_mode": "path-search" if refine_path else "path-opt",
        "uma_model": None,
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
    opt_mode_set: bool,
    tsopt_opt_mode_default: Optional[str],
    convert_files: bool,
    thresh_post: Optional[str],
    flatten_explicit: bool,
    flatten: Optional[bool],
    skip_final_freq: bool,
) -> Dict[str, Any]:
    """Assemble the `tsopt_overrides` dict consumed by the post-MEP TSOPT call.

    Extracted from the inline `if x is not None: tsopt_overrides[k] = ...`
    ladder in ``workflows/all.py:cli()`` (~20 LOC) so the dict-construction
    is independently unit-testable. Behavior preserved: a key only lands
    in the returned dict when the corresponding CLI flag was either
    explicitly set or has a non-None value to forward.
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
    if opt_mode_post_norm in {"grad", "hess"}:
        overrides["opt_mode"] = opt_mode_post_norm
    elif opt_mode_set:
        overrides["opt_mode"] = tsopt_opt_mode_default
    overrides["convert_files"] = bool(convert_files)
    if thresh_post is not None:
        overrides["thresh"] = str(thresh_post)
    if flatten_explicit:
        overrides["flatten"] = bool(flatten)
    if skip_final_freq:
        overrides["skip_final_freq"] = True
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
    hessian_calc_mode: Optional[str],
    convert_files: bool,
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
    if dump_override_requested:
        overrides["dump"] = bool(dump)
    if hessian_calc_mode is not None:
        overrides["hessian_calc_mode"] = hessian_calc_mode
    overrides["convert_files"] = bool(convert_files)
    return overrides


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
    """Append the backend / charge-embedding / link / mm flags to a CLI argv list.

    Extracted from the 7 near-identical inline blocks that each cli()
    invocation (_run_tsopt_on_hei / _run_freq_for_state / _run_opt_for_state
    / _run_dft_for_state / path_search / path_opt / build_irc_args) used
    to build before forwarding to the matching subcommand. All seven
    copies shared the same truthiness/None-skip semantics (most
    importantly: ``--no-embedcharge`` is only emitted when the user
    explicitly typed ``--no-embedcharge``, not when the CLI default
    False is in effect, so a YAML ``calc.embedcharge: true`` is not
    silently overridden).

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
    convert_files: bool,
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
    overrides["convert_files"] = bool(convert_files)
    return overrides


__all__ = [
    "AllContext",
    "copy_path_outputs_to_root",
    "promote_diag_for_root",
    "build_energy_level_dict",
    "build_pipeline_summary_payload",
    "build_tsopt_overrides",
    "build_freq_overrides",
    "build_dft_overrides",
    "append_backend_forwarding_args",
]
