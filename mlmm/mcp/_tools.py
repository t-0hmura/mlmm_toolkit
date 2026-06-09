"""MCP tool definitions — 1 tool per mlmm CLI subcommand.

Each tool wraps `mlmm <subcmd>` via subprocess + `summary.json` parsing.
Tool coverage includes the AMBER MM topology, ONIOM layer assignment, and
ONIOM input I/O subcommands.
"""
from __future__ import annotations

import tempfile
from pathlib import Path
from typing import Any, Optional

from mlmm.mcp._runner import run_subcmd


def _resolve_out_dir(out_dir: Optional[str], prefix: str) -> Path:
    """Pick a usable output directory: explicit > unique temp."""
    if out_dir:
        path = Path(out_dir).expanduser().resolve()
        path.mkdir(parents=True, exist_ok=True)
        return path
    return Path(tempfile.mkdtemp(prefix=f"mlmm_mcp_{prefix}_"))


def _shared_calc_flags(
    *,
    backend: Optional[str],
    precision: Optional[str],
    embedcharge: Optional[bool],
    embedcharge_cutoff: Optional[float],
    link_atom_method: Optional[str],
    mm_backend: Optional[str],
    use_cmap: Optional[bool] = None,
) -> list[str]:
    """Common backend / precision / embedcharge / mm_backend / cmap flags shared by stage runners."""
    args: list[str] = []
    if backend:
        args.extend(["-b", str(backend)])
    if precision:
        args.extend(["--precision", str(precision)])
    if embedcharge is not None:
        args.append("--embedcharge" if embedcharge else "--no-embedcharge")
    if embedcharge_cutoff is not None:
        args.extend(["--embedcharge-cutoff", str(embedcharge_cutoff)])
    if link_atom_method:
        args.extend(["--link-atom-method", str(link_atom_method)])
    if mm_backend:
        args.extend(["--mm-backend", str(mm_backend)])
    if use_cmap is not None:
        args.append("--cmap" if use_cmap else "--no-cmap")
    return args


def register_all(mcp) -> None:
    """Register every mlmm subcommand as an MCP tool on `mcp`."""

    # Topology / layer-prep helpers (mlmm-specific)

    @mcp.tool()
    def prepare_amber_topology(
        input_pdb: str,
        output_prefix: str,
        *,
        ligand_charge: Optional[str] = None,
        ligand_mult: Optional[str] = None,
        add_h: bool = False,
        ph: float = 7.0,
        add_ter: bool = True,
        ff_set: str = "ff19SB",
        keep_temp: bool = False,
        extra_args: Optional[list[str]] = None,
        timeout_seconds: Optional[float] = None,
    ) -> dict[str, Any]:
        """Generate AMBER parm7/rst7 topology files (CLI: `mlmm mm-parm`).

        Calls AmberTools (antechamber/parmchk2/tleap) under the hood.
        """
        argv: list[str] = ["mlmm", "mm-parm", "-i", input_pdb, "-o", output_prefix]
        if ligand_charge:
            argv.extend(["--ligand-charge", ligand_charge])
        if ligand_mult:
            argv.extend(["--ligand-mult", ligand_mult])
        argv.append("--add-h" if add_h else "--no-add-h")
        argv.extend(["--ph", str(ph)])
        argv.append("--add-ter" if add_ter else "--no-add-ter")
        argv.extend(["--ff-set", str(ff_set)])
        if keep_temp:
            argv.append("--keep-temp")
        if extra_args:
            argv.extend(extra_args)
        return run_subcmd(argv, out_dir=None, timeout=timeout_seconds).to_dict()

    @mcp.tool()
    def define_layer(
        input_pdb: str,
        output_pdb: str,
        *,
        model_pdb: Optional[str] = None,
        model_indices: Optional[str] = None,
        radius_freeze: Optional[float] = None,
        radius_partial_hessian: Optional[float] = None,
        one_based: Optional[bool] = None,
        extra_args: Optional[list[str]] = None,
        timeout_seconds: Optional[float] = None,
    ) -> dict[str, Any]:
        """Assign ML / MM-movable / MM-frozen B-factor layers (CLI: `mlmm define-layer`).

        Requires either `model_pdb` (PDB containing only the ML-region atoms)
        or `model_indices` (comma-separated atom indices, ranges allowed,
        e.g. ``"1,5,8-12"``). To carve the ML region from a ligand residue
        first, call `extract_pocket` (CLI: `mlmm extract`) and feed its output
        as `model_pdb` here.
        """
        argv: list[str] = ["mlmm", "define-layer", "-i", input_pdb, "-o", output_pdb]
        if model_pdb:
            argv.extend(["--model-pdb", model_pdb])
        if model_indices:
            argv.extend(["--model-indices", model_indices])
        if radius_freeze is not None:
            argv.extend(["--radius-freeze", str(radius_freeze)])
        if radius_partial_hessian is not None:
            argv.extend(["--radius-partial-hessian", str(radius_partial_hessian)])
        if one_based is not None:
            argv.append("--one-based" if one_based else "--zero-based")
        if extra_args:
            argv.extend(extra_args)
        return run_subcmd(argv, out_dir=None, timeout=timeout_seconds).to_dict()

    @mcp.tool()
    def extract_pocket(
        complex_pdb: str,
        ligand_id: str,
        radius_angstrom: float,
        output_pdb: str,
        *,
        ligand_charge: Optional[str] = None,
        exclude_backbone: Optional[bool] = None,
        extra_args: Optional[list[str]] = None,
        timeout_seconds: Optional[float] = None,
    ) -> dict[str, Any]:
        """Extract a binding pocket model (CLI: `mlmm extract`).

        Cuts a sphere of radius `radius_angstrom` around the ligand. Used to
        prepare a focused active-site model for subsequent ONIOM analysis.
        `exclude_backbone=None` (default) defers to the CLI's own default
        (`--no-exclude-backbone` semantics for the bool flag in shipped
        defaults); pass True / False to override explicitly.
        """
        argv: list[str] = ["mlmm", "extract",
                           "-i", complex_pdb, "-c", ligand_id,
                           "-r", str(radius_angstrom), "-o", output_pdb]
        if ligand_charge:
            argv.extend(["--ligand-charge", ligand_charge])
        if exclude_backbone is True:
            argv.append("--exclude-backbone")
        elif exclude_backbone is False:
            argv.append("--no-exclude-backbone")
        if extra_args:
            argv.extend(extra_args)
        return run_subcmd(argv, out_dir=None, timeout=timeout_seconds).to_dict()

    # Core stage runners (opt / tsopt / irc / freq) — ONIOM-aware

    @mcp.tool()
    def optimize_geometry(
        input_pdb: str,
        parm7: str,
        charge: int,
        multiplicity: int,
        *,
        ligand_charge: Optional[str] = None,
        opt_mode: str = "grad",
        max_cycles: int = 10000,
        thresh: Optional[str] = None,
        coord_type: Optional[str] = None,
        microiter: Optional[bool] = None,
        backend: Optional[str] = None,
        precision: Optional[str] = None,
        embedcharge: Optional[bool] = None,
        embedcharge_cutoff: Optional[float] = None,
        link_atom_method: Optional[str] = None,
        mm_backend: Optional[str] = None,
        print_every: Optional[int] = None,
        dump_trajectory: bool = False,
        out_dir: Optional[str] = None,
        extra_args: Optional[list[str]] = None,
        timeout_seconds: Optional[float] = None,
    ) -> dict[str, Any]:
        """Optimize an ONIOM-layered geometry (CLI: `mlmm opt`).

        `input_pdb` must already carry ML/MM B-factor assignments
        (e.g. produced by `define_layer` or `extract_pocket` followed by
        manual / scripted layer setup). Default `opt_mode="grad"` matches
        the CLI default; pass `"hess"` for RFO Hessian-based optimization.
        """
        od = _resolve_out_dir(out_dir, "opt")
        argv: list[str] = ["mlmm", "opt", "-i", input_pdb,
                           "--parm", parm7, "-q", str(charge), "-m", str(multiplicity)]
        if ligand_charge:
            argv.extend(["-l", ligand_charge])
        argv.extend(["--opt-mode", opt_mode, "--max-cycles", str(max_cycles)])
        if thresh:
            argv.extend(["--thresh", thresh])
        if coord_type:
            argv.extend(["--coord-type", coord_type])
        if microiter is not None:
            argv.append("--microiter" if microiter else "--no-microiter")
        if print_every is not None:
            argv.extend(["--print-every", str(print_every)])
        if dump_trajectory:
            argv.append("--dump")
        argv.extend(_shared_calc_flags(
            backend=backend, precision=precision,
            embedcharge=embedcharge, embedcharge_cutoff=embedcharge_cutoff,
            link_atom_method=link_atom_method, mm_backend=mm_backend,
        ))
        argv.extend(["--out-json"])
        argv.extend(["--out-dir", str(od)])
        if extra_args:
            argv.extend(extra_args)
        return run_subcmd(argv, out_dir=od, timeout=timeout_seconds).to_dict()

    @mcp.tool()
    def find_transition_state(
        input_pdb: str,
        parm7: str,
        charge: int,
        multiplicity: int,
        *,
        ligand_charge: Optional[str] = None,
        opt_mode: str = "hess",
        max_cycles: int = 10000,
        thresh: Optional[str] = None,
        coord_type: Optional[str] = None,
        microiter: Optional[bool] = None,
        flatten: Optional[bool] = None,
        hessian_calc_mode: Optional[str] = None,
        backend: Optional[str] = None,
        precision: Optional[str] = None,
        embedcharge: Optional[bool] = None,
        embedcharge_cutoff: Optional[float] = None,
        link_atom_method: Optional[str] = None,
        mm_backend: Optional[str] = None,
        print_every: Optional[int] = None,
        skip_final_freq: bool = False,
        out_dir: Optional[str] = None,
        extra_args: Optional[list[str]] = None,
        timeout_seconds: Optional[float] = None,
    ) -> dict[str, Any]:
        """ONIOM TS optimization (CLI: `mlmm tsopt`).

        Notes:
        - opt_mode 'trim' (Helgaker 1991) / 'rsprfo' (Banerjee 1985) are
          non-microiter TS optimizers; the server passes --no-microiter
          automatically when those modes are set (the mlmm CLI emits a
          warning otherwise).
        """
        od = _resolve_out_dir(out_dir, "tsopt")
        argv: list[str] = ["mlmm", "tsopt", "-i", input_pdb,
                           "--parm", parm7, "-q", str(charge), "-m", str(multiplicity)]
        if ligand_charge:
            argv.extend(["-l", ligand_charge])
        argv.extend(["--opt-mode", opt_mode, "--max-cycles", str(max_cycles)])
        if thresh:
            argv.extend(["--thresh", thresh])
        if coord_type:
            argv.extend(["--coord-type", coord_type])
        if microiter is not None:
            argv.append("--microiter" if microiter else "--no-microiter")
        elif opt_mode in ("trim", "rsprfo"):
            argv.append("--no-microiter")  # mlmm CLI ignores microiter for these anyway
        if flatten is not None:
            argv.append("--flatten" if flatten else "--no-flatten")
        if hessian_calc_mode:
            argv.extend(["--hessian-calc-mode", hessian_calc_mode])
        if print_every is not None:
            argv.extend(["--print-every", str(print_every)])
        if skip_final_freq:
            argv.append("--skip-final-freq")
        argv.extend(_shared_calc_flags(
            backend=backend, precision=precision,
            embedcharge=embedcharge, embedcharge_cutoff=embedcharge_cutoff,
            link_atom_method=link_atom_method, mm_backend=mm_backend,
        ))
        argv.extend(["--out-json"])
        argv.extend(["--out-dir", str(od)])
        if extra_args:
            argv.extend(extra_args)
        return run_subcmd(argv, out_dir=od, timeout=timeout_seconds).to_dict()

    @mcp.tool()
    def run_irc(
        input_pdb: str,
        parm7: str,
        charge: int,
        multiplicity: int,
        *,
        ligand_charge: Optional[str] = None,
        max_cycles: Optional[int] = None,
        step_size: Optional[float] = None,
        root: Optional[int] = None,
        forward: Optional[bool] = None,
        backward: Optional[bool] = None,
        backend: Optional[str] = None,
        precision: Optional[str] = None,
        embedcharge: Optional[bool] = None,
        embedcharge_cutoff: Optional[float] = None,
        link_atom_method: Optional[str] = None,
        mm_backend: Optional[str] = None,
        hessian_calc_mode: Optional[str] = None,
        irc_pos_def: Optional[bool] = None,
        out_dir: Optional[str] = None,
        extra_args: Optional[list[str]] = None,
        timeout_seconds: Optional[float] = None,
    ) -> dict[str, Any]:
        """ONIOM IRC integration from a TS geometry (CLI: `mlmm irc`)."""
        od = _resolve_out_dir(out_dir, "irc")
        argv: list[str] = ["mlmm", "irc", "-i", input_pdb,
                           "--parm", parm7, "-q", str(charge), "-m", str(multiplicity)]
        if ligand_charge:
            argv.extend(["-l", ligand_charge])
        if max_cycles is not None:
            argv.extend(["--max-cycles", str(max_cycles)])
        if step_size is not None:
            argv.extend(["--step-size", str(step_size)])
        if root is not None:
            argv.extend(["--root", str(root)])
        if forward is not None:
            argv.append("--forward" if forward else "--no-forward")
        if backward is not None:
            argv.append("--backward" if backward else "--no-backward")
        if hessian_calc_mode:
            argv.extend(["--hessian-calc-mode", hessian_calc_mode])
        if irc_pos_def is not None:
            argv.append("--irc-pos-def" if irc_pos_def else "--no-irc-pos-def")
        argv.extend(_shared_calc_flags(
            backend=backend, precision=precision,
            embedcharge=embedcharge, embedcharge_cutoff=embedcharge_cutoff,
            link_atom_method=link_atom_method, mm_backend=mm_backend,
        ))
        argv.extend(["--out-json"])
        argv.extend(["--out-dir", str(od)])
        if extra_args:
            argv.extend(extra_args)
        return run_subcmd(argv, out_dir=od, timeout=timeout_seconds).to_dict()

    @mcp.tool()
    def compute_frequencies(
        input_pdb: str,
        parm7: str,
        charge: int,
        multiplicity: int,
        *,
        ligand_charge: Optional[str] = None,
        temperature: Optional[float] = None,
        backend: Optional[str] = None,
        precision: Optional[str] = None,
        hessian_calc_mode: Optional[str] = None,
        out_dir: Optional[str] = None,
        extra_args: Optional[list[str]] = None,
        timeout_seconds: Optional[float] = None,
    ) -> dict[str, Any]:
        """ONIOM vibrational analysis + thermochemistry (CLI: `mlmm freq`)."""
        od = _resolve_out_dir(out_dir, "freq")
        argv: list[str] = ["mlmm", "freq", "-i", input_pdb,
                           "--parm", parm7, "-q", str(charge), "-m", str(multiplicity)]
        if ligand_charge:
            argv.extend(["-l", ligand_charge])
        if temperature is not None:
            argv.extend(["--temperature", str(temperature)])
        if hessian_calc_mode:
            argv.extend(["--hessian-calc-mode", hessian_calc_mode])
        if backend:
            argv.extend(["-b", backend])
        if precision:
            argv.extend(["--precision", precision])
        argv.extend(["--out-json"])
        argv.extend(["--out-dir", str(od)])
        if extra_args:
            argv.extend(extra_args)
        return run_subcmd(argv, out_dir=od, timeout=timeout_seconds).to_dict()

    @mcp.tool()
    def run_single_point_oniom(
        input_pdb: str,
        parm7: str,
        *,
        charge: Optional[int] = None,
        multiplicity: Optional[int] = None,
        ligand_charge: Optional[str] = None,
        model_pdb: Optional[str] = None,
        model_indices: Optional[str] = None,
        detect_layer: Optional[bool] = None,
        freeze_atoms: Optional[str] = None,
        hess_cutoff: Optional[float] = None,
        movable_cutoff: Optional[float] = None,
        do_hess: bool = False,
        hessian_calc_mode: Optional[str] = None,
        backend: Optional[str] = None,
        precision: Optional[str] = None,
        embedcharge: Optional[bool] = None,
        embedcharge_cutoff: Optional[float] = None,
        link_atom_method: Optional[str] = None,
        mm_backend: Optional[str] = None,
        use_cmap: Optional[bool] = None,
        convert_files: Optional[bool] = None,
        print_every: Optional[int] = None,
        out_dir: Optional[str] = None,
        extra_args: Optional[list[str]] = None,
        timeout_seconds: Optional[float] = None,
    ) -> dict[str, Any]:
        """Single-point ONIOM energy / forces (and optional Hessian) (CLI: `mlmm sp`).

        `input_pdb` must already carry ML/MM B-factor assignments
        (e.g. produced by `define_layer` or `extract_pocket`). Either
        `--model-pdb` / `--model-indices` or `--detect-layer` must resolve the
        ML region. Set `do_hess=True` to also compute and save the full ONIOM
        Hessian to `hessian.npy`.
        """
        od = _resolve_out_dir(out_dir, "sp")
        argv: list[str] = ["mlmm", "sp", "-i", input_pdb, "--parm", parm7]
        if charge is not None:
            argv.extend(["-q", str(charge)])
        if multiplicity is not None:
            argv.extend(["-m", str(multiplicity)])
        if ligand_charge:
            argv.extend(["-l", ligand_charge])
        if model_pdb:
            argv.extend(["--model-pdb", model_pdb])
        if model_indices:
            argv.extend(["--model-indices", model_indices])
        if detect_layer is not None:
            argv.append("--detect-layer" if detect_layer else "--no-detect-layer")
        if freeze_atoms:
            argv.extend(["--freeze-atoms", freeze_atoms])
        if hess_cutoff is not None:
            argv.extend(["--radius-partial-hessian", str(hess_cutoff)])
        if movable_cutoff is not None:
            argv.extend(["--radius-freeze", str(movable_cutoff)])
        if do_hess:
            argv.append("--hess")
        if hessian_calc_mode:
            argv.extend(["--hessian-calc-mode", hessian_calc_mode])
        if convert_files is not None:
            argv.append("--convert-files" if convert_files else "--no-convert-files")
        if print_every is not None:
            argv.extend(["--print-every", str(print_every)])
        argv.extend(_shared_calc_flags(
            backend=backend, precision=precision,
            embedcharge=embedcharge, embedcharge_cutoff=embedcharge_cutoff,
            link_atom_method=link_atom_method, mm_backend=mm_backend,
            use_cmap=use_cmap,
        ))
        argv.extend(["--out-json"])
        argv.extend(["--out-dir", str(od)])
        if extra_args:
            argv.extend(extra_args)
        return run_subcmd(argv, out_dir=od, timeout=timeout_seconds).to_dict()


    @mcp.tool()
    def scan_1d(
        input_pdb: str,
        parm7: str,
        charge: int,
        multiplicity: int,
        scan_lists: str,
        *,
        ligand_charge: Optional[str] = None,
        max_step_size: Optional[float] = None,
        relax_max_cycles: Optional[int] = None,
        opt_mode: Optional[str] = None,
        thresh: Optional[str] = None,
        backend: Optional[str] = None,
        precision: Optional[str] = None,
        embedcharge: Optional[bool] = None,
        embedcharge_cutoff: Optional[float] = None,
        out_dir: Optional[str] = None,
        extra_args: Optional[list[str]] = None,
        timeout_seconds: Optional[float] = None,
    ) -> dict[str, Any]:
        """1D ONIOM scan with harmonic restraints (CLI: `mlmm scan`)."""
        od = _resolve_out_dir(out_dir, "scan")
        argv: list[str] = ["mlmm", "scan", "-i", input_pdb,
                           "--parm", parm7, "-q", str(charge), "-m", str(multiplicity),
                           "--scan-lists", scan_lists]
        if ligand_charge:
            argv.extend(["-l", ligand_charge])
        if max_step_size is not None:
            argv.extend(["--max-step-size", str(max_step_size)])
        if relax_max_cycles is not None:
            argv.extend(["--relax-max-cycles", str(relax_max_cycles)])
        if opt_mode:
            argv.extend(["--opt-mode", opt_mode])
        if thresh:
            argv.extend(["--thresh", thresh])
        argv.extend(_shared_calc_flags(
            backend=backend, precision=precision,
            embedcharge=embedcharge, embedcharge_cutoff=embedcharge_cutoff,
            link_atom_method=None, mm_backend=None,
        ))
        argv.extend(["--out-json"])
        argv.extend(["--out-dir", str(od)])
        if extra_args:
            argv.extend(extra_args)
        return run_subcmd(argv, out_dir=od, timeout=timeout_seconds).to_dict()

    @mcp.tool()
    def scan_2d(
        input_pdb: str,
        parm7: str,
        charge: int,
        multiplicity: int,
        scan_lists: str,
        *,
        ligand_charge: Optional[str] = None,
        max_step_size: Optional[float] = None,
        relax_max_cycles: Optional[int] = None,
        thresh: Optional[str] = None,
        backend: Optional[str] = None,
        precision: Optional[str] = None,
        out_dir: Optional[str] = None,
        extra_args: Optional[list[str]] = None,
        timeout_seconds: Optional[float] = None,
    ) -> dict[str, Any]:
        """2D ONIOM scan (CLI: `mlmm scan2d`)."""
        od = _resolve_out_dir(out_dir, "scan2d")
        argv: list[str] = ["mlmm", "scan2d", "-i", input_pdb,
                           "--parm", parm7, "-q", str(charge), "-m", str(multiplicity),
                           "--scan-lists", scan_lists]
        if ligand_charge:
            argv.extend(["-l", ligand_charge])
        if max_step_size is not None:
            argv.extend(["--max-step-size", str(max_step_size)])
        if relax_max_cycles is not None:
            argv.extend(["--relax-max-cycles", str(relax_max_cycles)])
        if thresh:
            argv.extend(["--thresh", thresh])
        if backend:
            argv.extend(["-b", backend])
        if precision:
            argv.extend(["--precision", precision])
        argv.extend(["--out-json"])
        argv.extend(["--out-dir", str(od)])
        if extra_args:
            argv.extend(extra_args)
        return run_subcmd(argv, out_dir=od, timeout=timeout_seconds).to_dict()

    @mcp.tool()
    def scan_3d(
        input_pdb: str,
        parm7: str,
        charge: int,
        multiplicity: int,
        scan_lists: str,
        *,
        ligand_charge: Optional[str] = None,
        max_step_size: Optional[float] = None,
        relax_max_cycles: Optional[int] = None,
        thresh: Optional[str] = None,
        backend: Optional[str] = None,
        precision: Optional[str] = None,
        out_dir: Optional[str] = None,
        extra_args: Optional[list[str]] = None,
        timeout_seconds: Optional[float] = None,
    ) -> dict[str, Any]:
        """3D ONIOM scan (CLI: `mlmm scan3d`)."""
        od = _resolve_out_dir(out_dir, "scan3d")
        argv: list[str] = ["mlmm", "scan3d", "-i", input_pdb,
                           "--parm", parm7, "-q", str(charge), "-m", str(multiplicity),
                           "--scan-lists", scan_lists]
        if ligand_charge:
            argv.extend(["-l", ligand_charge])
        if max_step_size is not None:
            argv.extend(["--max-step-size", str(max_step_size)])
        if relax_max_cycles is not None:
            argv.extend(["--relax-max-cycles", str(relax_max_cycles)])
        if thresh:
            argv.extend(["--thresh", thresh])
        if backend:
            argv.extend(["-b", backend])
        if precision:
            argv.extend(["--precision", precision])
        argv.extend(["--out-json"])
        argv.extend(["--out-dir", str(od)])
        if extra_args:
            argv.extend(extra_args)
        return run_subcmd(argv, out_dir=od, timeout=timeout_seconds).to_dict()

    @mcp.tool()
    def optimize_path(
        reactant_pdb: str,
        product_pdb: str,
        parm7: str,
        charge: int,
        multiplicity: int,
        *,
        ligand_charge: Optional[str] = None,
        max_nodes: Optional[int] = None,
        max_cycles: Optional[int] = None,
        mep_mode: Optional[str] = None,
        backend: Optional[str] = None,
        precision: Optional[str] = None,
        out_dir: Optional[str] = None,
        extra_args: Optional[list[str]] = None,
        timeout_seconds: Optional[float] = None,
    ) -> dict[str, Any]:
        """ONIOM MEP optimization (CLI: `mlmm path-opt`)."""
        od = _resolve_out_dir(out_dir, "path_opt")
        argv: list[str] = ["mlmm", "path-opt", "-i", reactant_pdb, product_pdb,
                           "--parm", parm7, "-q", str(charge), "-m", str(multiplicity)]
        if ligand_charge:
            argv.extend(["-l", ligand_charge])
        if max_nodes is not None:
            argv.extend(["--max-nodes", str(max_nodes)])
        if max_cycles is not None:
            argv.extend(["--max-cycles", str(max_cycles)])
        if mep_mode:
            argv.extend(["--mep-mode", mep_mode])
        if backend:
            argv.extend(["-b", backend])
        if precision:
            argv.extend(["--precision", precision])
        argv.append("--out-json")
        argv.extend(["--out-dir", str(od)])
        if extra_args:
            argv.extend(extra_args)
        return run_subcmd(argv, out_dir=od, timeout=timeout_seconds).to_dict()

    @mcp.tool()
    def search_paths(
        input_pdb: str,
        parm7: str,
        charge: int,
        multiplicity: int,
        *,
        ligand_charge: Optional[str] = None,
        max_nodes: Optional[int] = None,
        mep_mode: Optional[str] = None,
        refine_mode: Optional[str] = None,
        backend: Optional[str] = None,
        precision: Optional[str] = None,
        out_dir: Optional[str] = None,
        extra_args: Optional[list[str]] = None,
        timeout_seconds: Optional[float] = None,
    ) -> dict[str, Any]:
        """Recursive ONIOM pathway search (CLI: `mlmm path-search`)."""
        od = _resolve_out_dir(out_dir, "path_search")
        argv: list[str] = ["mlmm", "path-search", "-i", input_pdb,
                           "--parm", parm7, "-q", str(charge), "-m", str(multiplicity)]
        if ligand_charge:
            argv.extend(["-l", ligand_charge])
        if max_nodes is not None:
            argv.extend(["--max-nodes", str(max_nodes)])
        if mep_mode:
            argv.extend(["--mep-mode", mep_mode])
        if refine_mode:
            argv.extend(["--refine-mode", refine_mode])
        if backend:
            argv.extend(["-b", backend])
        if precision:
            argv.extend(["--precision", precision])
        argv.extend(["--out-dir", str(od)])
        if extra_args:
            argv.extend(extra_args)
        return run_subcmd(argv, out_dir=od, timeout=timeout_seconds).to_dict()

    @mcp.tool()
    def run_full_pipeline(
        reactant_complex_pdb: str,
        product_complex_pdb: Optional[str] = None,
        *,
        charge: Optional[int] = None,
        ligand_charge: Optional[str] = None,
        multiplicity: Optional[int] = None,
        extract_ligand: Optional[str] = None,
        extract_radius: Optional[float] = None,
        max_cycles: Optional[int] = None,
        thresh: Optional[str] = None,
        thresh_post: Optional[str] = None,
        backend: Optional[str] = None,
        precision: Optional[str] = None,
        refine_path: Optional[bool] = None,
        do_tsopt: Optional[bool] = None,
        do_dft: Optional[bool] = None,
        do_thermo: Optional[bool] = None,
        out_dir: Optional[str] = None,
        extra_args: Optional[list[str]] = None,
        timeout_seconds: Optional[float] = None,
    ) -> dict[str, Any]:
        """End-to-end ONIOM workflow extract → MEP → TS → IRC → freq → DFT (CLI: `mlmm all`).

        Inputs may be raw complex PDBs (with `extract_ligand` + `extract_radius`
        to auto-extract the active site) or pre-extracted ONIOM-layered PDBs.
        """
        od = _resolve_out_dir(out_dir, "all")
        argv: list[str] = ["mlmm", "all", "-i", reactant_complex_pdb]
        if product_complex_pdb:
            argv.append(product_complex_pdb)
        if charge is not None:
            argv.extend(["-q", str(charge)])
        if ligand_charge:
            argv.extend(["-l", ligand_charge])
        if multiplicity is not None:
            argv.extend(["-m", str(multiplicity)])
        if extract_ligand:
            argv.extend(["-c", extract_ligand])
        if extract_radius is not None:
            argv.extend(["-r", str(extract_radius)])
        if max_cycles is not None:
            argv.extend(["--max-cycles", str(max_cycles)])
        if thresh:
            argv.extend(["--thresh", thresh])
        if thresh_post:
            argv.extend(["--thresh-post", thresh_post])
        if backend:
            argv.extend(["-b", backend])
        if precision:
            argv.extend(["--precision", precision])
        if refine_path is not None:
            argv.append("--refine-path" if refine_path else "--no-refine-path")
        if do_tsopt is not None:
            argv.extend(["--tsopt", "true" if do_tsopt else "false"])
        if do_dft is not None:
            argv.extend(["--dft", "true" if do_dft else "false"])
        if do_thermo is not None:
            argv.extend(["--thermo", "true" if do_thermo else "false"])
        argv.extend(["--out-dir", str(od)])
        if extra_args:
            argv.extend(extra_args)
        return run_subcmd(argv, out_dir=od, timeout=timeout_seconds).to_dict()

    @mcp.tool()
    def run_single_point_dft(
        input_pdb: str,
        parm7: str,
        charge: int,
        multiplicity: int,
        *,
        ligand_charge: Optional[str] = None,
        func_basis: Optional[str] = None,
        out_dir: Optional[str] = None,
        extra_args: Optional[list[str]] = None,
        timeout_seconds: Optional[float] = None,
    ) -> dict[str, Any]:
        """ONIOM-embedded single-point DFT (CLI: `mlmm dft`).

        ``func_basis`` is the single combined CLI argument (``FUNC/BASIS``,
        e.g. ``"wb97m-v/def2-tzvpd"``).
        """
        od = _resolve_out_dir(out_dir, "dft")
        argv: list[str] = ["mlmm", "dft", "-i", input_pdb,
                           "--parm", parm7, "-q", str(charge), "-m", str(multiplicity)]
        if ligand_charge:
            argv.extend(["-l", ligand_charge])
        if func_basis:
            argv.extend(["--func-basis", func_basis])
        argv.extend(["--out-json"])
        argv.extend(["--out-dir", str(od)])
        if extra_args:
            argv.extend(extra_args)
        return run_subcmd(argv, out_dir=od, timeout=timeout_seconds).to_dict()

    # ONIOM I/O (Gaussian / ORCA)

    @mcp.tool()
    def export_oniom_input(
        input_layered_pdb: str,
        parm7: str,
        charge: int,
        multiplicity: int,
        output_file: str,
        *,
        model_pdb: Optional[str] = None,
        format_engine: str = "g16",
        extra_args: Optional[list[str]] = None,
        timeout_seconds: Optional[float] = None,
    ) -> dict[str, Any]:
        """Export an ONIOM input deck for Gaussian or ORCA (CLI: `mlmm oniom-export`)."""
        argv: list[str] = ["mlmm", "oniom-export",
                           "--parm", parm7, "-i", input_layered_pdb,
                           "-q", str(charge), "-m", str(multiplicity),
                           "-o", output_file]
        if model_pdb:
            argv.extend(["--model-pdb", model_pdb])
        if format_engine:
            argv.extend(["--mode", format_engine])
        if extra_args:
            argv.extend(extra_args)
        return run_subcmd(argv, out_dir=None, timeout=timeout_seconds).to_dict()

    @mcp.tool()
    def import_oniom_input(
        input_file: str,
        output_prefix: str,
        *,
        extra_args: Optional[list[str]] = None,
        timeout_seconds: Optional[float] = None,
    ) -> dict[str, Any]:
        """Import an ONIOM input deck and reconstruct XYZ / layered PDB (CLI: `mlmm oniom-import`)."""
        argv: list[str] = ["mlmm", "oniom-import", "-i", input_file, "-o", output_prefix]
        if extra_args:
            argv.extend(extra_args)
        return run_subcmd(argv, out_dir=None, timeout=timeout_seconds).to_dict()

    # Structure / I/O helpers (no out_dir, no summary.json)

    @mcp.tool()
    def add_element_info(
        input_pdb: str,
        output_pdb: str,
        *,
        overwrite: bool = False,
        extra_args: Optional[list[str]] = None,
        timeout_seconds: Optional[float] = None,
    ) -> dict[str, Any]:
        """Repair / add PDB element columns (CLI: `mlmm add-elem-info`)."""
        argv: list[str] = ["mlmm", "add-elem-info", "-i", input_pdb, "-o", output_pdb]
        if overwrite:
            argv.append("--overwrite")
        if extra_args:
            argv.extend(extra_args)
        return run_subcmd(argv, out_dir=None, timeout=timeout_seconds).to_dict()

    @mcp.tool()
    def fix_altloc(
        input_pdb: str,
        output_pdb: str,
        *,
        extra_args: Optional[list[str]] = None,
        timeout_seconds: Optional[float] = None,
    ) -> dict[str, Any]:
        """Resolve PDB alternate locations (CLI: `mlmm fix-altloc`)."""
        argv: list[str] = ["mlmm", "fix-altloc", "-i", input_pdb, "-o", output_pdb]
        if extra_args:
            argv.extend(extra_args)
        return run_subcmd(argv, out_dir=None, timeout=timeout_seconds).to_dict()

    @mcp.tool()
    def plot_trajectory(
        input_trj_xyz: str,
        output_png: str,
        *,
        extra_args: Optional[list[str]] = None,
        timeout_seconds: Optional[float] = None,
    ) -> dict[str, Any]:
        """Plot an energy profile from a trajectory (CLI: `mlmm trj2fig`)."""
        argv: list[str] = ["mlmm", "trj2fig", "-i", input_trj_xyz, "-o", output_png]
        if extra_args:
            argv.extend(extra_args)
        return run_subcmd(argv, out_dir=None, timeout=timeout_seconds).to_dict()

    @mcp.tool()
    def plot_energy_diagram(
        energies: str,
        output_png: str,
        *,
        extra_args: Optional[list[str]] = None,
        timeout_seconds: Optional[float] = None,
    ) -> dict[str, Any]:
        """Draw an energy diagram (CLI: `mlmm energy-diagram`).

        `energies`: Python-literal list of floats, e.g. "[0, 12.5, 4.3, 18.7, -1.2]".
        """
        argv: list[str] = ["mlmm", "energy-diagram", "-i", energies, "-o", output_png]
        if extra_args:
            argv.extend(extra_args)
        return run_subcmd(argv, out_dir=None, timeout=timeout_seconds).to_dict()

    @mcp.tool()
    def detect_bond_changes(
        reactant_pdb: str,
        product_pdb: str,
        *,
        extra_args: Optional[list[str]] = None,
        timeout_seconds: Optional[float] = None,
    ) -> dict[str, Any]:
        """Detect bond changes between two PDB structures (CLI: `mlmm bond-summary`)."""
        argv: list[str] = ["mlmm", "bond-summary", "-i", reactant_pdb, product_pdb]
        if extra_args:
            argv.extend(extra_args)
        return run_subcmd(argv, out_dir=None, timeout=timeout_seconds).to_dict()
