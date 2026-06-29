"""
ML/MM minimum-energy path optimization via the Growing String Method (GSM) or Direct Max Flux (DMF).

Example:
    mlmm path-opt -i reac.pdb prod.pdb --parm real.parm7 --model-pdb ml_region.pdb -q 0
    mlmm path-opt -i reac.pdb prod.pdb --parm real.parm7 -q 0 --mep-mode dmf

For detailed documentation, see: docs/path_opt.md
"""

from __future__ import annotations

from copy import deepcopy
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence, Set, Tuple

import gc
import logging
import sys
import traceback
import textwrap

logger = logging.getLogger(__name__)

import click
import numpy as np
import time
import torch

from pysisyphus.helpers import geom_loader
from pysisyphus.cos.GrowingString import GrowingString
from pysisyphus.optimizers.StringOptimizer import StringOptimizer
from pysisyphus.optimizers.exceptions import OptimizationError
from pysisyphus.optimizers.LBFGS import LBFGS

from mlmm.backends.mlmm_calc import mlmm, MLMMASECalculator
from mlmm.workflows.opt import (
    GEOM_KW as OPT_GEOM_KW,
    CALC_KW as OPT_CALC_KW,
    LBFGS_KW as OPT_LBFGS_KW,
    _parse_freeze_atoms as _parse_freeze_atoms_opt,
    _normalize_geom_freeze as _normalize_geom_freeze_opt,
)
from mlmm.workflows.opt import _convert_yaml_layer_atoms_1to0
from mlmm.core.utils import (
    apply_layer_freeze_constraints,
    convert_xyz_to_pdb,
    set_convert_file_enabled,
    is_convert_file_enabled,
    deep_update,
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
    PreparedInputStructure,
    parse_indices_string,
    build_model_pdb_from_bfactors,
    build_model_pdb_from_indices,
    echo_resolved_device,
)
from mlmm.cli.common_options import add_coord_type_option, add_ml_layer_detection_options, add_precision_option, add_backend_model_option, add_deterministic_option, add_allow_charge_mult_mismatch_option
from mlmm.cli.decorators import resolve_yaml_sources, load_merged_yaml_cfg, make_is_param_explicit, _write_error_json, render_cli_exception
from mlmm.workflows.align_freeze import align_and_refine_sequence_inplace
from mlmm.core.defaults import (
    BFACTOR_FROZEN,
    BFACTOR_ML,
    BFACTOR_MOVABLE_MM,
    DMF_KW as _DMF_KW_DEFAULT,
    GS_KW as _GS_KW_DEFAULT,
    OUT_DIR_PATH_OPT,
    STOPT_KW as _STOPT_KW_DEFAULT,
    THRESH_CHOICES,
)


# Defaults (overridden by YAML/CLI)

# Geometry (input handling) — share defaults with opt.py
GEOM_KW: Dict[str, Any] = deepcopy(OPT_GEOM_KW)

# ML/MM calculator settings — share defaults with opt.py
CALC_KW: Dict[str, Any] = deepcopy(OPT_CALC_KW)

# LBFGS (used for optional endpoint pre-optimization)
LBFGS_KW: Dict[str, Any] = deepcopy(OPT_LBFGS_KW)

# DMF (Direct Max Flux) defaults
DMF_KW: Dict[str, Any] = deepcopy(_DMF_KW_DEFAULT)

# GrowingString (path representation)
GS_KW: Dict[str, Any] = deepcopy(_GS_KW_DEFAULT)

# StringOptimizer (optimization control)
STOPT_KW: Dict[str, Any] = deepcopy(_STOPT_KW_DEFAULT)

def _load_two_endpoints(
    inputs: Sequence[PreparedInputStructure],
    coord_type: str,
    base_freeze: Sequence[int],
) -> Sequence:
    """
    Load the two endpoint structures and set `freeze_atoms` as needed.
    """
    geoms = []
    for prepared in inputs:
        geom_path = prepared.geom_path
        cfg: Dict[str, Any] = {"freeze_atoms": list(base_freeze)}
        freeze = merge_freeze_atom_indices(cfg)
        g = geom_loader(geom_path, coord_type=coord_type, freeze_atoms=freeze)
        g.freeze_atoms = np.array(freeze, dtype=int)
        geoms.append(g)
    return geoms


# Helpers shared with opt.py (imported for consistency)
_parse_freeze_atoms = _parse_freeze_atoms_opt
_normalize_geom_freeze = _normalize_geom_freeze_opt


# B-factor annotation helpers

def _parse_pdb_atoms_for_indexing(pdb_path: Path) -> List[Dict[str, Any]]:
    """Parse PDB ATOM/HETATM records for indexing/matching."""
    atoms: List[Dict[str, Any]] = []
    with open(pdb_path, "r") as f:
        for line in f:
            if line.startswith(("ATOM  ", "HETATM")):
                # Ensure line is long enough
                s = line.rstrip("\n")
                s = s + (" " * (80 - len(s))) if len(s) < 80 else s
                serial_str = s[6:11].strip()
                resseq_str = s[22:26].strip()
                try:
                    serial = int(serial_str) if serial_str else None
                except ValueError:
                    serial = None
                try:
                    resseq = int(resseq_str) if resseq_str else None
                except ValueError:
                    resseq = None
                atom = {
                    "line": s,
                    "serial": serial,
                    "name": s[12:16].strip(),
                    "altloc": s[16].strip() if len(s) > 16 else "",
                    "resname": s[17:20].strip(),
                    "chain": s[21].strip() if len(s) > 21 else "",
                    "resseq": resseq,
                    "icode": s[26].strip() if len(s) > 26 else "",
                }
                atoms.append(atom)
    return atoms


def _compute_ml_indices_from_model_and_ref(ref_pdb: Path, model_pdb: Path) -> Set[int]:
    """
    Compute 0-based atom indices (in the order of ATOM/HETATM records of ref_pdb)
    that belong to the ML region defined by model_pdb.

    Matching strategy (robust, fallback-based):
      1) Exact key match on (name, altloc, resname, chain, resseq, icode)
      2) Key match ignoring altloc when model altloc is blank
      3) Serial number match
    """
    ref_atoms = _parse_pdb_atoms_for_indexing(ref_pdb)
    model_atoms = _parse_pdb_atoms_for_indexing(model_pdb)

    # Build maps for the reference structure
    key_to_indices: Dict[Tuple[str, str, str, str, Optional[int], str], List[int]] = {}
    key_wo_alt_to_indices: Dict[Tuple[str, str, str, Optional[int], str], List[int]] = {}
    serial_to_index: Dict[int, int] = {}

    for idx, a in enumerate(ref_atoms):
        key = (a["name"], a["altloc"], a["resname"], a["chain"], a["resseq"], a["icode"])
        key_wo = (a["name"], a["resname"], a["chain"], a["resseq"], a["icode"])
        key_to_indices.setdefault(key, []).append(idx)
        key_wo_alt_to_indices.setdefault(key_wo, []).append(idx)
        if a["serial"] is not None:
            # If duplicated serials exist, keep the first occurrence
            serial_to_index.setdefault(a["serial"], idx)

    ml_indices: Set[int] = set()
    misses = 0

    for ma in model_atoms:
        key = (ma["name"], ma["altloc"], ma["resname"], ma["chain"], ma["resseq"], ma["icode"])
        key_wo = (ma["name"], ma["resname"], ma["chain"], ma["resseq"], ma["icode"])

        idx: Optional[int] = None

        candidates = key_to_indices.get(key)
        if candidates and len(candidates) == 1:
            idx = candidates[0]
        elif candidates and len(candidates) > 1:
            # Multiple; try to disambiguate via serial number if possible
            if ma["serial"] is not None:
                si = serial_to_index.get(ma["serial"])
                if si in candidates:
                    idx = si

        if idx is None and (ma["altloc"] == "" or ma["altloc"] == " "):
            candidates2 = key_wo_alt_to_indices.get(key_wo)
            if candidates2:
                idx = candidates2[0]  # pick first

        if idx is None and ma["serial"] is not None:
            idx = serial_to_index.get(ma["serial"])

        if idx is not None:
            ml_indices.add(idx)
        else:
            misses += 1

    if misses:
        click.echo(f"[annotate] WARNING: {misses} ML atoms from '{model_pdb.name}' could not be mapped to '{ref_pdb.name}'.", err=True)

    return ml_indices


def _apply_bfactor_annotations_inplace(
    pdb_path: Path,
    ml_indices: Set[int],
    freeze_indices: Sequence[int],
    beta_ml: float = BFACTOR_ML,
    beta_freeze: float = BFACTOR_FROZEN,
    beta_both: float = BFACTOR_ML,
) -> None:
    """
    In-place update of B-factors for PDB ATOM/HETATM records.

    Rules (new 3-layer encoding):
      - ML only:       0  (BFACTOR_ML)
      - Freeze only:  20  (BFACTOR_FROZEN)
      - ML ∩ Freeze:   0  (ML takes precedence)
      - Others:       10  (BFACTOR_MOVABLE_MM)

    The index for lookups is the 0-based position among ATOM/HETATM
    records and resets at each MODEL record for multi-model PDBs.
    """
    freeze_set: Set[int] = set(int(i) for i in (freeze_indices or []))
    ml_set: Set[int] = set(int(i) for i in (ml_indices or set()))

    def _format_b(b: float) -> str:
        # PDB tempFactor field is 6.2 width
        return f"{b:6.2f}"

    # Read and process lines
    lines_out: List[str] = []
    atom_idx = 0  # 0-based within a MODEL (or entire file if no MODEL)

    with open(pdb_path, "r") as f:
        lines = f.readlines()

    for line in lines:
        if line.startswith("MODEL"):
            # Reset index at each model
            atom_idx = 0
            lines_out.append(line)
            continue

        if line.startswith(("ATOM  ", "HETATM")):
            s = line.rstrip("\n")
            # Pad to at least 66 chars so we can safely replace tempFactor (cols 61-66)
            if len(s) < 66:
                s = s + (" " * (66 - len(s)))

            # Decide B-factor for this atom index
            if (atom_idx in ml_set) and (atom_idx in freeze_set):
                b = beta_both
            elif atom_idx in ml_set:
                b = beta_ml
            elif atom_idx in freeze_set:
                b = beta_freeze
            else:
                b = BFACTOR_MOVABLE_MM

            s = s[:60] + _format_b(b) + s[66:]
            # Ensure trailing newline
            s = s if s.endswith("\n") else s + "\n"
            lines_out.append(s)

            atom_idx += 1
        else:
            lines_out.append(line)

    with open(pdb_path, "w") as f:
        f.writelines(lines_out)

    click.echo(
        f"[annotate] Updated B-factors in '{pdb_path}' "
        f"(ML={BFACTOR_ML:.0f}, MovableMM={BFACTOR_MOVABLE_MM:.0f}, "
        f"FrozenMM={BFACTOR_FROZEN:.0f}; {len(ml_set)} ML, {len(freeze_set)} frozen)."
    )



def _select_hei_index(energies: Sequence[float]) -> int:
    """Pick an HEI index preferring internal local maxima."""
    E = np.array(energies, dtype=float)
    nE = int(len(E))
    hei_idx = None
    if nE >= 3:
        candidates = [i for i in range(1, nE - 1)
                      if E[i] > E[i - 1] and E[i] > E[i + 1]]
        if candidates:
            hei_idx = int(max(candidates, key=lambda i: E[i]))
        else:
            hei_idx = 1 + int(np.argmax(E[1:-1]))
    if hei_idx is None:
        hei_idx = int(np.argmax(E))
    return hei_idx


# DMF (Direct Max Flux) MEP optimization

def _is_cuda_oom(exc: BaseException) -> bool:
    """True if `exc` looks like a CUDA out-of-memory (torch.cuda.OutOfMemoryError or a
    RuntimeError carrying 'out of memory'), so the DMF gpu backend can advise --dmf-backend cpu."""
    if type(exc).__name__ == "OutOfMemoryError":
        return True
    msg = str(exc).lower()
    return "out of memory" in msg or "cuda oom" in msg


def _run_dmf_mep(
    geoms: Sequence,
    calc_cfg: Dict[str, Any],
    out_dir_path: Path,
    input_paths: Sequence[Path],
    max_nodes: int,
    fix_atoms: Sequence[int],
    dmf_cfg: Optional[Dict[str, Any]] = None,
    ml_indices_set: Optional[Set[int]] = None,
    freeze_atoms_final: Optional[Sequence[int]] = None,
) -> None:
    """Run Direct Max Flux (DMF) MEP optimization between two endpoints.

    Uses pydmf with harmonic constraints for frozen atoms; the ML/MM ONIOM calculator is
    wrapped as an ASE calculator. The backend is selected by ``dmf_cfg["backend"]``: ``"gpu"``
    (default) imports the PyTorch backend ``dmf.torch`` (CUDA), ``"cpu"`` imports ``dmf`` (NumPy).

    References:
    [1] S.-i. Koda and S. Saito, JCTC, 20, 2798-2811 (2024). doi: 10.1021/acs.jctc.3c01246
    [2] S.-i. Koda and S. Saito, JCTC, 20, 7176-7187 (2024). doi: 10.1021/acs.jctc.4c00792
    [3] S.-i. Koda and S. Saito, JCTC, 21, 3513-3522 (2025). doi: 10.1021/acs.jctc.4c01549
    """
    dmf_backend = str((dmf_cfg or {}).get("backend", "gpu")).strip().lower()
    try:
        from ase.io import read as ase_read, write as ase_write
        from ase.calculators.mixing import SumCalculator
        if dmf_backend == "cpu":
            from dmf import DirectMaxFlux, interpolate_fbenm
        else:
            from dmf.torch import DirectMaxFlux, interpolate_fbenm
    except Exception as e:
        raise RuntimeError(
            "DMF mode (--mep-mode dmf) requires ase, cyipopt, and pydmf>=1.2 "
            "(`pip install 'pydmf[torch]'` for the default GPU backend; the `cpu` "
            "backend needs only `pip install pydmf`). "
            f"Import error: {e}"
        ) from e

    from mlmm.workflows.restraints import HarmonicFixAtoms

    def _geom_to_ase(g):
        from io import StringIO
        return ase_read(StringIO(g.as_xyz()), format="xyz")

    fix_atoms = list(sorted(set(map(int, fix_atoms))))

    ref_images = [_geom_to_ase(g) for g in geoms]
    charge = int(calc_cfg.get("model_charge", 0))
    spin = int(calc_cfg.get("model_mult", 1))
    for img in ref_images:
        img.info["charge"] = charge
        img.info["spin"] = spin

    # Build the ONIOM ASE calculator
    shared_pysis_calc = mlmm(**calc_cfg)
    ase_calc = MLMMASECalculator(core=shared_pysis_calc.core)

    dmf_cfg = deep_update(dict(DMF_KW), dmf_cfg)
    fbenm_opts: Dict[str, Any] = dict(dmf_cfg.get("fbenm_options", {}))
    cfbenm_opts: Dict[str, Any] = dict(dmf_cfg.get("cfbenm_options", {}))
    dmf_opts: Dict[str, Any] = dict(dmf_cfg.get("dmf_options", {}))
    update_teval = bool(dmf_opts.pop("update_teval", False))
    k_fix = float(dmf_cfg.get("k_fix", DMF_KW["k_fix"]))

    # Default-mode IPOPT options: print_level=0 silences the per-iteration
    # IPOPT banner + MUMPS license header that would otherwise print 16+
    # times per pipeline (one per align+scan refinement). Under `-v` we
    # let dmf use its own default (print_level=5 = full table). User-
    # supplied `dmf_cfg["ipopt_options"]` always wins.
    from mlmm.core.utils import is_verbose
    ipopt_opts: Dict[str, Any] = dict(dmf_cfg.get("ipopt_options", {}))
    if "print_level" not in ipopt_opts and not is_verbose():
        ipopt_opts["print_level"] = 0

    # Run FB-ENM interpolation
    click.echo("\n====== DMF: FB-ENM interpolation ======\n", narrative=True)
    mxflx_fbenm = interpolate_fbenm(
        ref_images,
        nmove=max(1, int(max_nodes)),
        fbenm_only_endpoints=bool(dmf_cfg.get("fbenm_only_endpoints", False)),
        correlated=bool(dmf_cfg.get("correlated", False)),
        sequential=bool(dmf_cfg.get("sequential", False)),
        output_file=str(out_dir_path / "dmf_fbenm_ipopt.out"),
        fbenm_options=fbenm_opts,
        cfbenm_options=cfbenm_opts,
        dmf_options=dmf_opts,
        ipopt_options=ipopt_opts,
    )

    initial_trj = out_dir_path / "dmf_initial_trj.xyz"
    ase_write(initial_trj, mxflx_fbenm.images, format="xyz")
    click.echo(f"[write] Wrote '{initial_trj}' ({len(mxflx_fbenm.images)} images).")

    # Convert initial trajectory to PDB if possible
    if input_paths[0].suffix.lower() == ".pdb" and is_convert_file_enabled():
        try:
            initial_pdb = initial_trj.with_suffix(".pdb")
            convert_xyz_to_pdb(initial_trj, input_paths[0].resolve(), initial_pdb)
            click.echo(f"[convert] Wrote '{initial_pdb}'.")
        except Exception as e:
            click.echo(f"[convert] WARNING: {e}", err=True)

    coefs = mxflx_fbenm.coefs.copy()

    # Create DirectMaxFlux object
    click.echo("\n====== DMF: Direct Max Flux optimization ======\n", narrative=True)
    mxflx = DirectMaxFlux(
        ref_images,
        coefs=coefs,
        nmove=max(1, int(max_nodes)),
        update_teval=update_teval,
        remove_rotation_and_translation=bool(
            dmf_opts.get("remove_rotation_and_translation", False)
        ),
        mass_weighted=bool(dmf_opts.get("mass_weighted", False)),
        parallel=bool(dmf_opts.get("parallel", False)),
        eps_vel=float(dmf_opts.get("eps_vel", DMF_KW["dmf_options"]["eps_vel"])),
        eps_rot=float(dmf_opts.get("eps_rot", DMF_KW["dmf_options"]["eps_rot"])),
        beta=float(dmf_opts.get("beta", DMF_KW["dmf_options"]["beta"])),
    )

    # Endpoint anchor: take ref_positions ONCE from the reactant endpoint so the
    # HarmonicFixAtoms restraint pulls every image back to a shared reference
    # (per-image self-reference made the restraint a no-op — atoms drifted freely
    # along the path).
    fix_ref_positions = (
        ref_images[0].get_positions()[fix_atoms] if fix_atoms else None
    )

    # Assign calculators to images
    for image in mxflx.images:
        if "charge" not in image.info:
            image.info["charge"] = charge
        if "spin" not in image.info:
            image.info["spin"] = spin

        if fix_atoms:
            harmonic_calc = HarmonicFixAtoms(
                indices=fix_atoms,
                ref_positions=fix_ref_positions,
                k_fix=k_fix,
            )
            image.calc = SumCalculator([ase_calc, harmonic_calc])
        else:
            image.calc = ase_calc

    mxflx.add_ipopt_options({"output_file": str(out_dir_path / "dmf_ipopt.out")})
    # Same default-mode IPOPT quiet as the FB-ENM interpolation above:
    # the banner + MUMPS license header would otherwise reprint per solve.
    if not is_verbose() and "print_level" not in ipopt_opts:
        mxflx.add_ipopt_options({"print_level": 0})
    max_cycles = dmf_cfg.get("max_cycles") if isinstance(dmf_cfg, dict) else None
    if max_cycles is not None:
        try:
            max_iter = int(max_cycles)
            if max_iter > 0:
                mxflx.add_ipopt_options({"max_iter": max_iter})
        except Exception:
            logger.debug("Failed to set ipopt max_iter option", exc_info=True)
    mxflx.solve(tol="tight")
    click.echo("\n====== DMF: optimization finished ======\n", narrative=True)

    # Evaluate final energies using the PySisyphus calculator for consistency
    from pysisyphus.constants import ANG2BOHR
    energies = []
    for image in mxflx.images:
        elems = image.get_chemical_symbols()
        coords_bohr = np.asarray(image.get_positions(), dtype=float).reshape(-1, 3) * ANG2BOHR
        energies.append(float(shared_pysis_calc.get_energy(elems, coords_bohr)["energy"]))
    hei_idx = _select_hei_index(energies)

    # Write final trajectory
    final_trj = out_dir_path / "final_geometries_trj.xyz"
    blocks = []
    for _idx, (image, E) in enumerate(zip(mxflx.images, energies)):
        from io import StringIO
        buf = StringIO()
        ase_write(buf, image, format="xyz")
        s = buf.getvalue()
        lines = s.splitlines()
        if len(lines) >= 2 and lines[0].strip().isdigit():
            lines[1] = f"{E:.12f}"
        blocks.append("\n".join(lines) + "\n")
    with open(final_trj, "w") as f:
        f.write("".join(blocks))
    click.echo(f"[write] Wrote '{final_trj}' with energy.")

    # Convert to PDB
    if input_paths[0].suffix.lower() == ".pdb" and is_convert_file_enabled():
        ref_pdb = input_paths[0].resolve()
        try:
            final_pdb = out_dir_path / "final_geometries.pdb"
            convert_xyz_to_pdb(final_trj, ref_pdb, final_pdb)
            click.echo(f"[convert] Wrote '{final_pdb}'.")
            _apply_bfactor_annotations_inplace(
                final_pdb,
                ml_indices=ml_indices_set or set(),
                freeze_indices=freeze_atoms_final or [],
            )
        except Exception as e:
            click.echo(f"[convert] WARNING: {e}", err=True)

    # Write HEI
    hei_geom = mxflx.images[hei_idx]
    hei_E = energies[hei_idx]
    hei_xyz = out_dir_path / "hei.xyz"
    from io import StringIO
    buf = StringIO()
    ase_write(buf, hei_geom, format="xyz")
    s = buf.getvalue()
    lines = s.splitlines()
    if len(lines) >= 2 and lines[0].strip().isdigit():
        lines[1] = f"{hei_E:.12f}"
        s = "\n".join(lines) + "\n"
    with open(hei_xyz, "w") as f:
        f.write(s)
    click.echo(f"[write] Wrote '{hei_xyz}' (HEI index={hei_idx}).")

    if input_paths[0].suffix.lower() == ".pdb" and is_convert_file_enabled():
        ref_pdb = input_paths[0].resolve()
        hei_pdb = out_dir_path / "hei.pdb"
        try:
            convert_xyz_to_pdb(hei_xyz, ref_pdb, hei_pdb)
            click.echo(f"[convert] Wrote '{hei_pdb}'.")
            _apply_bfactor_annotations_inplace(
                hei_pdb,
                ml_indices=ml_indices_set or set(),
                freeze_indices=freeze_atoms_final or [],
            )
        except Exception as e:
            click.echo(f"[convert] WARNING: {e}", err=True)



@click.command(
    help="MEP optimization via the Growing String method or Direct Max Flux.",
    context_settings={"help_option_names": ["-h", "--help"]},
)
@click.option(
    "-i", "--input",
    "input_paths",
    type=click.Path(path_type=Path, exists=True, dir_okay=False),
    nargs=2,
    required=True,
    help="Two endpoint structures (reactant/product); both must be full-enzyme PDBs.",
)
@click.option(
    "-q",
    "--charge",
    type=int,
    required=False,
    help="ML region charge. Required unless --ligand-charge is provided.",
)
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
    help="Spin multiplicity (2S+1). Defaults to 1 when omitted.",
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
    help="DMF compute backend (--mep-mode dmf only): gpu (dmf.torch / CUDA) or cpu (dmf / NumPy). "
    "On a GPU out-of-memory error, retry with cpu.",
)
@click.option("--max-nodes", type=int, default=GS_KW["max_nodes"], show_default=True,
              help="Number of internal nodes (for GSM: string has max_nodes+2 images including endpoints; for DMF: number of path waypoints).")
@click.option("--max-cycles", type=int, default=300, show_default=True, help="Maximum optimization cycles.")
@click.option(
    "--climb/--no-climb",
    default=True,
    show_default=True,
    help="Search for a transition state (climbing image) after path growth.",
)
@click.option(
    "--preopt/--no-preopt",
    # Default True matches GS_KW.fix_first/fix_last semantics and the
    # mlmm-all.py forwarding: endpoints are typically pre-relaxed before
    # string growth to avoid GSM step inflation.
    default=True,
    show_default=True,
    help="Pre-optimize the two endpoint structures with LBFGS before string growth.",
)
@click.option("--preopt-max-cycles", "preopt_max_cycles", type=int, default=10000, show_default=True,
              help="Maximum LBFGS cycles for endpoint pre-optimization when --preopt is enabled.")
@click.option(
    "--fix-ends/--no-fix-ends",
    default=True,
    show_default=True,
    help="Fix endpoint structures during path growth.",
)
@click.option(
    "--dump/--no-dump",
    default=False,
    show_default=True,
    help="Dump optimizer trajectory/restarts during the run.",
)
@click.option("-o", "--out-dir", "out_dir", type=str, default=OUT_DIR_PATH_OPT, show_default=True,
              help="Output directory.")
@click.option(
    "--thresh",
    type=click.Choice(THRESH_CHOICES, case_sensitive=False),
    default=None,
    help="Convergence preset for the string optimizer.",
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
    help="Validate options and print the execution plan without running path optimization.",
)
@click.option(
    "--parm",
    "real_parm7",
    type=click.Path(path_type=Path, exists=True, dir_okay=False),
    required=True,
    help="Amber parm7 topology for the enzyme complex (MM layers).",
)
@click.option(
    "--model-pdb",
    type=click.Path(path_type=Path, exists=True, dir_okay=False),
    required=False,
    help="PDB defining the ML region (atom IDs used by the ML/MM calculator). "
         "Optional when --detect-layer is enabled.",
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
    "--freeze-atoms",
    "freeze_atoms_cli",
    type=str,
    default=None,
    help="Comma-separated 1-based indices to freeze (applied to every image).",
)
@click.option(
    "--hess-cutoff",
    "hess_cutoff",
    type=float,
    default=None,
    show_default=False,
    help="Distance cutoff (Å) from ML region for MM atoms to include in Hessian calculation. "
         "Applied to movable MM atoms and can be combined with --detect-layer.",
)
@click.option(
    "--movable-cutoff",
    "movable_cutoff",
    type=float,
    default=None,
    show_default=False,
    help="Distance cutoff (Å) from ML region for movable MM atoms. MM atoms beyond this are frozen. "
         "Providing --movable-cutoff disables --detect-layer.",
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
    help="Enable xTB point-charge embedding correction for MM→ML environmental effects (experimental).",
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
    "--out-json/--no-out-json",
    "out_json",
    default=False,
    show_default=True,
    help="Write machine-readable result.json to out_dir.",
)
# Full template PDBs for XYZ→PDB conversion and topology reference
@click.option(
    "--ref-pdb",
    "ref_pdb_paths",
    type=click.Path(path_type=Path, exists=True, dir_okay=False),
    multiple=True,
    default=None,
    help=("Full-size template PDBs in the same order as --input. "
          "Required when using XYZ inputs to provide topology and B-factor information."),
)
@add_ml_layer_detection_options()
@add_precision_option()
@add_backend_model_option()
@add_deterministic_option()
@add_allow_charge_mult_mismatch_option()
@click.pass_context
def cli(
    ctx: click.Context,
    input_paths: Sequence[Path],
    ref_pdb_paths: Sequence[Path],
    charge: Optional[int],
    ligand_charge: Optional[str],
    spin: Optional[int],
    mep_mode: str,
    dmf_backend: str,
    max_nodes: int,
    max_cycles: int,
    climb: bool,
    preopt: bool,
    preopt_max_cycles: int,
    fix_ends: bool,
    dump: bool,
    out_dir: str,
    thresh: Optional[str],
    config_yaml: Optional[Path],
    show_config: bool,
    dry_run: bool,
    real_parm7: Path,
    model_pdb: Optional[Path],
    model_indices_str: Optional[str],
    model_indices_one_based: bool,
    detect_layer: bool,
    freeze_atoms_cli: Optional[str],
    hess_cutoff: Optional[float],
    movable_cutoff: Optional[float],
    convert_files: bool,
    backend: Optional[str],
    embedcharge: bool,
    embedcharge_cutoff: Optional[float],
    link_atom_method: Optional[str],
    mm_backend: Optional[str],
    use_cmap: Optional[bool],
    out_json: bool,
    precision: Optional[str],
    backend_model: Optional[str],
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

    input_paths = tuple(Path(p) for p in input_paths)
    prepared_inputs = [prepare_input_structure(p) for p in input_paths]
    try:
        time_start = time.perf_counter()

        if len(prepared_inputs) != 2:
            click.echo("ERROR: Provide exactly two endpoint structures (-i reactant product).", err=True)
            sys.exit(1)

        # Accept XYZ inputs when a matching --ref-pdb is supplied: overlay the XYZ
        # coordinates onto the reference PDB topology so path-opt runs on full
        # ML/MM PDBs (mirrors the --ref-pdb support already in path-search/scan).
        ref_list = list(ref_pdb_paths) if ref_pdb_paths else []
        _resolved_inputs: List[Path] = []
        for i, src in enumerate(input_paths):
            suffix = src.suffix.lower()
            if suffix == ".pdb":
                _resolved_inputs.append(src)
            elif suffix == ".xyz":
                if i >= len(ref_list):
                    click.echo(
                        f"ERROR: XYZ input '{src.name}' requires a corresponding --ref-pdb "
                        "for topology/B-factor info.",
                        err=True,
                    )
                    sys.exit(1)
                Path(out_dir).mkdir(parents=True, exist_ok=True)
                materialized = (Path(out_dir) / f"path_opt_input_{i:02d}.pdb").resolve()
                convert_xyz_to_pdb(src.resolve(), Path(ref_list[i]).resolve(), materialized)
                _resolved_inputs.append(materialized)
            else:
                click.echo(
                    f"ERROR: '{src.name}': unsupported format. Use .pdb or .xyz (with --ref-pdb).",
                    err=True,
                )
                sys.exit(1)
        input_paths = tuple(_resolved_inputs)
        prepared_inputs = [prepare_input_structure(p) for p in input_paths]

        config_layer_cfg = load_yaml_dict(config_yaml)
        override_layer_cfg = load_yaml_dict(override_yaml)

        mep_mode_kind = mep_mode.strip().lower()

        geom_cfg = dict(GEOM_KW)
        calc_cfg = dict(CALC_KW)
        gs_cfg = dict(GS_KW)
        stopt_cfg = dict(STOPT_KW)
        lbfgs_cfg = dict(LBFGS_KW)
        dmf_cfg = dict(DMF_KW)

        apply_yaml_overrides(
            config_layer_cfg,
            [
                (geom_cfg, (("geom",),)),
                (calc_cfg, (("calc",), ("mlmm",))),
                (gs_cfg, (("gs",),)),
                (stopt_cfg, (("stopt",), ("opt",))),
                (lbfgs_cfg, (("opt", "lbfgs"), ("lbfgs",), ("stopt", "lbfgs"))),
                (dmf_cfg, (("dmf",),)),
            ],
        )

        # CLI explicit overrides (after config YAML, before override YAML)
        if backend is not None:
            calc_cfg["backend"] = str(backend).lower()
        if precision is not None:
            from mlmm.backends import apply_precision_to_calc_cfg
            apply_precision_to_calc_cfg(calc_cfg, precision)
        if backend_model is not None:
            from mlmm.backends import apply_backend_model_to_calc_cfg
            apply_backend_model_to_calc_cfg(calc_cfg, backend_model)
        if _is_param_explicit("embedcharge"):
            calc_cfg["embedcharge"] = bool(embedcharge)
        if _is_param_explicit("embedcharge_cutoff"):
            calc_cfg["embedcharge_cutoff"] = embedcharge_cutoff
        geom_cfg["coord_type"] = "cart"  # no microiteration: DLC over MM atoms is meaningless; fixed to cart
        if link_atom_method is not None:
            calc_cfg["link_atom_method"] = str(link_atom_method).lower()
        if mm_backend is not None:
            calc_cfg["mm_backend"] = str(mm_backend).lower()
        if use_cmap is not None:
            calc_cfg["use_cmap"] = use_cmap

        if _is_param_explicit("max_nodes"):
            gs_cfg["max_nodes"] = int(max_nodes)
        if _is_param_explicit("max_cycles"):
            stopt_cfg["max_cycles"] = int(max_cycles)
            stopt_cfg["stop_in_when_full"] = int(max_cycles)
            dmf_cfg["max_cycles"] = int(max_cycles)
        if _is_param_explicit("dmf_backend"):
            dmf_cfg["backend"] = str(dmf_backend).lower()
        if _is_param_explicit("climb"):
            gs_cfg["climb"] = bool(climb)
            gs_cfg["climb_lanczos"] = bool(climb)
        if _is_param_explicit("fix_ends"):
            gs_cfg["fix_first"] = bool(fix_ends)
            gs_cfg["fix_last"] = bool(fix_ends)
        if _is_param_explicit("dump"):
            stopt_cfg["dump"] = bool(dump)
            lbfgs_cfg["dump"] = bool(dump)
        if _is_param_explicit("out_dir"):
            stopt_cfg["out_dir"] = out_dir
            lbfgs_cfg["out_dir"] = out_dir
        if _is_param_explicit("thresh") and thresh is not None:
            stopt_cfg["thresh"] = str(thresh)
            lbfgs_cfg["thresh"] = str(thresh)
        if _is_param_explicit("detect_layer"):
            calc_cfg["use_bfactor_layers"] = bool(detect_layer)
        if _is_param_explicit("hess_cutoff") and hess_cutoff is not None:
            calc_cfg["hess_cutoff"] = float(hess_cutoff)
        if _is_param_explicit("movable_cutoff") and movable_cutoff is not None:
            calc_cfg["movable_cutoff"] = float(movable_cutoff)
            calc_cfg["use_bfactor_layers"] = False
        if _is_param_explicit("preopt_max_cycles"):
            lbfgs_cfg["max_cycles"] = int(preopt_max_cycles)

        resolved_charge = charge
        resolved_spin = spin
        for prepared in prepared_inputs:
            resolved_charge, resolved_spin = resolve_charge_spin_or_raise(
                prepared,
                resolved_charge,
                resolved_spin,
                ligand_charge=ligand_charge,
                prefix="[path-opt]",
            )
        # CLI-resolved charge/spin (from -q / -l derivation, or -m / spin_default)
        # always wins over the CALC_KW default carried in calc_cfg.
        # Same pattern as opt.py fix (commit f9a2799).
        calc_cfg["model_charge"] = int(resolved_charge)
        calc_cfg["model_mult"] = int(resolved_spin)

        if model_pdb is not None:
            calc_cfg["model_pdb"] = str(model_pdb)
        calc_cfg["input_pdb"] = str(input_paths[0])
        calc_cfg["real_parm7"] = str(real_parm7)

        apply_yaml_overrides(
            override_layer_cfg,
            [
                (geom_cfg, (("geom",),)),
                (calc_cfg, (("calc",), ("mlmm",))),
                (gs_cfg, (("gs",),)),
                (stopt_cfg, (("stopt",), ("opt",))),
                (lbfgs_cfg, (("opt", "lbfgs"), ("lbfgs",), ("stopt", "lbfgs"))),
                (dmf_cfg, (("dmf",),)),
            ],
        )

        try:
            geom_freeze = _normalize_geom_freeze(geom_cfg.get("freeze_atoms"))
        except click.BadParameter as e:
            click.echo(f"ERROR: {e}", err=True)
            sys.exit(1)
        geom_cfg["freeze_atoms"] = geom_freeze
        _convert_yaml_layer_atoms_1to0(calc_cfg)

        try:
            cli_freeze = _parse_freeze_atoms(freeze_atoms_cli)
        except click.BadParameter as e:
            click.echo(f"ERROR: {e}", err=True)
            sys.exit(1)

        model_indices: Optional[List[int]] = None
        if model_indices_str:
            try:
                model_indices = parse_indices_string(model_indices_str, one_based=model_indices_one_based)
            except click.BadParameter as e:
                click.echo(f"ERROR: {e}", err=True)
                sys.exit(1)
        if cli_freeze:
            merge_freeze_atom_indices(geom_cfg, cli_freeze)

        freeze_atoms_final = list(geom_cfg.get("freeze_atoms") or [])
        calc_cfg["freeze_atoms"] = freeze_atoms_final

        # Keep optimizer alignment policy deterministic.
        stopt_cfg["align"] = False
        stopt_cfg["stop_in_when_full"] = int(stopt_cfg.get("max_cycles", STOPT_KW["max_cycles"]))

        out_dir_path = Path(stopt_cfg["out_dir"]).resolve()
        preopt_max_cycles_effective = int(lbfgs_cfg.get("max_cycles", preopt_max_cycles))

        # movable_cutoff implies full distance-based layer assignment.
        # hess_cutoff alone can be combined with --detect-layer.
        detect_layer_enabled = bool(calc_cfg.get("use_bfactor_layers", True))
        model_pdb_cfg = calc_cfg.get("model_pdb")
        if calc_cfg.get("movable_cutoff") is not None:
            if detect_layer_enabled:
                click.echo("[layer] movable_cutoff is set; disabling --detect-layer.", err=True)
            detect_layer_enabled = False
            calc_cfg["use_bfactor_layers"] = False

        layer_source_pdb = input_paths[0]
        if detect_layer_enabled and layer_source_pdb.suffix.lower() != ".pdb":
            click.echo("ERROR: --detect-layer requires a PDB input.", err=True)
            sys.exit(1)

        # Preflight: with --detect-layer, a PDB that carries no B-factor
        # layering (every atom B≈ML) makes the WHOLE system the ML region,
        # which silently OOMs the MLIP on large inputs. Warn with the remedy
        # before any GPU work instead of dying later with a CUDA OOM.
        if detect_layer_enabled and layer_source_pdb.suffix.lower() == ".pdb":
            try:
                _n_atoms = 0
                _n_nonml = 0
                with open(layer_source_pdb, encoding="utf-8", errors="ignore") as _fh:
                    for _ln in _fh:
                        if _ln.startswith(("ATOM", "HETATM")):
                            _n_atoms += 1
                            try:
                                _b = float(_ln[60:66])
                            except ValueError:
                                continue
                            if abs(_b - BFACTOR_ML) > 1.0:
                                _n_nonml += 1
                if _n_atoms > 0 and _n_nonml == 0:
                    click.echo(
                        f"[layer] WARNING: --detect-layer is on but "
                        f"'{layer_source_pdb.name}' has no B-factor layering "
                        f"(all {_n_atoms} atoms are ML, B≈{BFACTOR_ML:.0f}). The "
                        f"ENTIRE system will be treated as the ML region, which "
                        f"typically exhausts GPU memory. Provide a layered PDB "
                        f"(B={BFACTOR_ML:.0f}/{BFACTOR_MOVABLE_MM:.0f}/"
                        f"{BFACTOR_FROZEN:.0f} via `mlmm define-layer`), or "
                        f"`--model-pdb`/`--model-indices`, or `--movable-cutoff`.",
                        err=True,
                    )
            except OSError:
                pass

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
                    click.echo("ERROR: Provide --model-pdb or --model-indices when --no-detect-layer.", err=True)
                    sys.exit(1)
            if (
                not detect_layer_enabled
                and model_pdb_cfg is None
                and model_indices
                and layer_source_pdb.suffix.lower() != ".pdb"
            ):
                click.echo("ERROR: --model-indices requires a PDB input.", err=True)
                sys.exit(1)
            click.echo(
                pretty_block(
                    "dry_run_plan",
                    {
                        "input_endpoints": [str(p) for p in input_paths],
                        "output_dir": str(out_dir_path),
                        "mep_mode": mep_mode_kind,
                        "fix_ends": bool(gs_cfg.get("fix_first", False) and gs_cfg.get("fix_last", False)),
                        "detect_layer": bool(detect_layer_enabled),
                        "model_region_source": model_region_source,
                        "model_indices_count": 0 if not model_indices else len(model_indices),
                        "preopt": bool(preopt),
                        "preopt_max_cycles": int(preopt_max_cycles_effective),
                        "will_run_path_opt": True,
                        "will_write_summary": True,
                        "backend": calc_cfg.get("backend", "uma"),
                        "embedcharge": bool(calc_cfg.get("embedcharge", False)),
                    },
                )
            )
            click.echo("[dry-run] Validation complete. Path optimization execution was skipped.")
            return

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
                    click.echo(f"ERROR: {e}", err=True)
                    sys.exit(1)
                click.echo(f"[layer] WARNING: {e} Falling back to explicit ML region.", err=True)
                detect_layer_enabled = False

        if not detect_layer_enabled:
            if model_pdb_cfg is None and not model_indices:
                click.echo("ERROR: Provide --model-pdb or --model-indices when --no-detect-layer.", err=True)
                sys.exit(1)
            if model_pdb_cfg is not None:
                model_pdb_path = Path(model_pdb_cfg)
            else:
                if layer_source_pdb.suffix.lower() != ".pdb":
                    click.echo("ERROR: --model-indices requires a PDB input.", err=True)
                    sys.exit(1)
                try:
                    model_pdb_path = build_model_pdb_from_indices(layer_source_pdb, out_dir_path, model_indices or [])
                except Exception as e:
                    click.echo(f"ERROR: {e}", err=True)
                    sys.exit(1)
            calc_cfg["use_bfactor_layers"] = False

        if model_pdb_path is None:
            click.echo("ERROR: Failed to resolve model PDB for the ML region.", err=True)
            sys.exit(1)

        calc_cfg["model_pdb"] = str(model_pdb_path)
        freeze_atoms_final = apply_layer_freeze_constraints(
            geom_cfg,
            calc_cfg,
            layer_info,
            echo_fn=click.echo,
        )

        for key in ("input_pdb", "real_parm7", "model_pdb", "mm_fd_dir"):
            val = calc_cfg.get(key)
            if val:
                calc_cfg[key] = str(Path(val).expanduser().resolve())

        # For display: resolved configuration (show only non-default values)
        echo_geom = format_freeze_atoms_for_echo(geom_cfg, key="freeze_atoms")
        echo_calc = format_freeze_atoms_for_echo(filter_calc_for_echo(calc_cfg), key="freeze_atoms")
        echo_gs = strip_inherited_keys(gs_cfg, GS_KW, mode="same")
        echo_stopt = strip_inherited_keys({**stopt_cfg, "out_dir": str(out_dir_path)}, STOPT_KW, mode="same")
        echo_lbfgs = strip_inherited_keys({**lbfgs_cfg, "out_dir": stopt_cfg.get("out_dir")}, LBFGS_KW, mode="same")

        click.echo(pretty_block("geom", echo_geom))
        click.echo(pretty_block("calc", echo_calc))
        if mep_mode_kind == "gsm":
            click.echo(pretty_block("gs", echo_gs))
            click.echo(pretty_block("stopt", echo_stopt))
            click.echo(pretty_block("lbfgs", echo_lbfgs))
        elif mep_mode_kind == "dmf":
            click.echo(pretty_block("dmf", dmf_cfg))
        click.echo(
            pretty_block(
                "run_flags",
                {
                    "mep_mode": mep_mode_kind,
                    "preopt": bool(preopt),
                    "preopt_max_cycles": int(preopt_max_cycles_effective),
                    "fix_ends": bool(gs_cfg.get("fix_first", False) and gs_cfg.get("fix_last", False)),
                },
            )
        )

        if int(stopt_cfg.get("max_cycles", 0)) <= 0:
            click.echo("[INFO] max_cycles <= 0: skipping path optimization.")
            return

        out_dir_path.mkdir(parents=True, exist_ok=True)

        source_paths = [prep.source_path for prep in prepared_inputs]

        # Pre-compute ML-region indices (0-based in ref PDB atom order) for later PDB annotation
        ml_indices_set: Set[int] = set()
        try:
            ref_pdb_for_map = source_paths[0]
            if ref_pdb_for_map.suffix.lower() == ".pdb":
                ml_indices_set = _compute_ml_indices_from_model_and_ref(
                    ref_pdb_for_map.resolve(),
                    Path(calc_cfg["model_pdb"]).resolve(),
                )
                click.echo(f"[annotate] ML-region atoms mapped: {len(ml_indices_set)}")
        except Exception as e:
            click.echo(f"[annotate] WARNING: Failed to pre-compute ML-region indices: {e}", err=True)

        # Load endpoints (if PDB, merge in link-parent freezing)
        geoms = _load_two_endpoints(
            inputs=prepared_inputs,
            coord_type=geom_cfg.get("coord_type", "cart"),
            base_freeze=geom_cfg.get("freeze_atoms", []),
        )

        # Shared ML/MM calculator (reuse the same instance for all images)
        shared_calc = mlmm(**calc_cfg)

        echo_resolved_device()

        # optional endpoint pre-optimization
        if preopt:
            try:
                click.echo("\n====== Pre-optimizing endpoints (LBFGS) ======\n", narrative=True)
                pre_dir_base = out_dir_path / "preopt"
                for i, g in enumerate(geoms):
                    try:
                        g.set_calculator(shared_calc)
                    except Exception:
                        logger.debug("Failed to set calculator on geometry", exc_info=True)
                    subdir = pre_dir_base / f"end{i:02d}"
                    subdir.mkdir(parents=True, exist_ok=True)
                    lbfgs_args = dict(lbfgs_cfg)
                    lbfgs_args.update({
                        "out_dir": str(subdir),
                        "max_cycles": int(preopt_max_cycles_effective),
                    })
                    optimizer = LBFGS(g, **lbfgs_args)
                    optimizer.run()
                    try:
                        final_xyz_path = optimizer.final_fn if isinstance(optimizer.final_fn, Path) else Path(optimizer.final_fn)
                        g_new = geom_loader(
                            final_xyz_path,
                            coord_type=geom_cfg.get("coord_type", "cart"),
                            freeze_atoms=getattr(g, "freeze_atoms", []),
                        )
                        try:
                            g_new.freeze_atoms = np.array(getattr(g, "freeze_atoms", []), dtype=int)
                        except Exception:
                            logger.debug("Failed to set freeze_atoms on new geometry", exc_info=True)
                        geoms[i] = g_new
                    except Exception as e:
                        click.echo(f"[preopt] WARNING: Failed to reload optimized endpoint #{i}: {e}", err=True)
                click.echo("[preopt] Completed endpoint pre-optimization.")
            except Exception as e:
                click.echo(f"[preopt] WARNING: Pre-optimization skipped due to error: {e}", err=True)

        # By default, apply external Kabsch alignment (if freeze_atoms exist, use only them)
        align_thresh = str(stopt_cfg.get("thresh", "gau"))
        try:
            click.echo("\n====== Aligning all inputs to the first structure (freeze-guided scan + relaxation) ======\n", narrative=True)
            _ = align_and_refine_sequence_inplace(
                geoms,
                thresh=align_thresh,
                shared_calc=shared_calc,
                out_dir=out_dir_path / "align_refine",
                verbose=True,
            )
            click.echo("[align] Completed input alignment.")
        except Exception as e:
            click.echo(f"[align] WARNING: alignment skipped: {e}", err=True)

        # Collect freeze_atoms for DMF
        fix_atoms: List[int] = []
        try:
            fix_atoms = sorted(
                {int(i) for g in geoms for i in getattr(g, "freeze_atoms", [])}
            )
        except Exception:
            logger.debug("Failed to extract freeze_atoms from geometries", exc_info=True)

        if mep_mode_kind == "dmf":
            try:
                _run_dmf_mep(
                    geoms,
                    calc_cfg,
                    out_dir_path,
                    input_paths,
                    max_nodes,
                    fix_atoms,
                    dmf_cfg=dmf_cfg,
                    ml_indices_set=ml_indices_set,
                    freeze_atoms_final=freeze_atoms_final,
                )
            except Exception as e:
                if str(dmf_cfg.get("backend", "gpu")).lower() != "cpu" and _is_cuda_oom(e):
                    click.echo(
                        "[dmf] GPU out of memory. Retry with `--dmf-backend cpu` "
                        "(NumPy backend; slower but not limited by GPU memory).",
                        err=True,
                    )
                else:
                    tb = "".join(traceback.format_exception(type(e), e, e.__traceback__))
                    click.echo(f"[dmf] ERROR: DMF optimization failed:\n{textwrap.indent(tb, '  ')}", err=True)
                sys.exit(3)
            click.echo(format_elapsed("[time] Elapsed Time for Path Opt (DMF)", time_start), narrative=True)

            if out_json:
                from mlmm.core.utils import write_result_json
                from pysisyphus.constants import AU2KCALPERMOL as _AU2KCAL
                # _run_dmf_mep writes hei.xyz; re-read energies from the trajectory
                _dmf_trj = out_dir_path / "final_geometries_trj.xyz"
                _dmf_energies: list = []
                if _dmf_trj.exists():
                    try:
                        # Comment line is a bare float; parse it directly
                        with open(_dmf_trj) as _f:
                            _lines = _f.readlines()
                        _i = 0
                        while _i < len(_lines):
                            _natoms_line = _lines[_i].strip()
                            if _natoms_line.isdigit():
                                _comment = _lines[_i + 1].strip()
                                try:
                                    _dmf_energies.append(float(_comment))
                                except ValueError:
                                    pass
                                _i += int(_natoms_line) + 2
                            else:
                                _i += 1
                    except Exception:
                        pass
                _dmf_hei_idx = 0
                _dmf_barrier = None
                _dmf_delta = None
                if _dmf_energies:
                    _dmf_hei_idx = int(np.argmax(_dmf_energies))
                    _dmf_e0 = _dmf_energies[0]
                    _dmf_barrier = (_dmf_energies[_dmf_hei_idx] - _dmf_e0) * _AU2KCAL
                    _dmf_delta = (_dmf_energies[-1] - _dmf_e0) * _AU2KCAL
                # DMF convergence is not directly returned by _run_dmf_mep;
                # infer from trajectory: if energies were parsed, mark as completed.
                _dmf_converged = None
                result_data_dmf: Dict[str, Any] = {
                    "status": "converged" if _dmf_converged else ("not_converged" if _dmf_converged is False else "completed"),
                    "converged": _dmf_converged,
                    "mep_mode": "dmf",
                    "backend": calc_cfg.get("backend", "uma"),
                    "charge": calc_cfg.get("model_charge"),
                    "spin": calc_cfg.get("model_mult"),
                    "reactant_energy_hartree": float(_dmf_energies[0]) if _dmf_energies else None,
                    "product_energy_hartree": float(_dmf_energies[-1]) if _dmf_energies else None,
                    "image_energies_hartree": [float(e) for e in _dmf_energies] if _dmf_energies else [],
                    "n_images": len(_dmf_energies) if _dmf_energies else None,
                    "hei_index": _dmf_hei_idx,
                    "hei_energy_hartree": float(_dmf_energies[_dmf_hei_idx]) if _dmf_energies else None,
                    "barrier_kcal": round(_dmf_barrier, 6) if _dmf_barrier is not None else None,
                    "delta_kcal": round(_dmf_delta, 6) if _dmf_delta is not None else None,
                    "files": {
                        "final_geometries_trj_xyz": "final_geometries_trj.xyz",
                        "hei_xyz": "hei.xyz",
                    },
                }
                for ext in (".pdb", ".gjf"):
                    f = out_dir_path / f"hei{ext}"
                    if f.exists():
                        result_data_dmf["files"][f"hei_{ext[1:]}"] = f.name
                write_result_json(
                    out_dir_path, result_data_dmf,
                    command="path-opt",
                    elapsed_seconds=time.perf_counter() - time_start,
                )

            return

        for g in geoms:
            g.set_calculator(shared_calc)

        def calc_getter():
            # Used when GrowingString generates new nodes
            return shared_calc

        gs = GrowingString(
            images=geoms,
            calc_getter=calc_getter,
            **gs_cfg,
        )

        # StringOptimizer expects 'out_dir' under the key "out_dir"
        opt_args = dict(stopt_cfg)
        opt_args["out_dir"] = str(out_dir_path)

        optimizer = StringOptimizer(
            geometry=gs,
            **{k: v for k, v in opt_args.items() if k != "type"}  # 'type' is just a tag
        )

        optimizer.run()

        final_trj = out_dir_path / "final_geometries_trj.xyz"
        try:
            try:
                energies = np.array(gs.energy, dtype=float)
                blocks = []
                for _idx, (geom, E) in enumerate(zip(gs.images, energies)):
                    s = geom.as_xyz()
                    lines = s.splitlines()
                    if len(lines) >= 2 and lines[0].strip().isdigit():
                        lines[1] = f"{E:.12f}"
                    s_mod = "\n".join(lines)
                    if not s_mod.endswith("\n"):
                        s_mod += "\n"
                    blocks.append(s_mod)
                annotated = "".join(blocks)
                with open(final_trj, "w") as f:
                    f.write(annotated)
                click.echo(f"[write] Wrote '{final_trj}' with energy.")
            except Exception:
                with open(final_trj, "w") as f:
                    f.write(gs.as_xyz())
                click.echo(f"[write] Wrote '{final_trj}'.")

            if input_paths[0].suffix.lower() == ".pdb" and is_convert_file_enabled():
                ref_pdb = input_paths[0].resolve()

                try:
                    out_pdb = out_dir_path / "final_geometries.pdb"
                    convert_xyz_to_pdb(final_trj, ref_pdb, out_pdb)
                    click.echo(f"[convert] Wrote '{out_pdb}'.")
                    # Annotate B-factors for ML & freeze atoms
                    _apply_bfactor_annotations_inplace(
                        out_pdb,
                        ml_indices=ml_indices_set,
                        freeze_indices=freeze_atoms_final,
                    )
                except Exception as e:
                    click.echo(f"[convert] WARNING: Failed to convert MEP path trajectory to PDB: {e}", err=True)

        except Exception as e:
            click.echo(f"[write] ERROR: Failed to write final trajectory: {e}", err=True)
            sys.exit(4)

        try:
            energies = np.array(gs.energy, dtype=float)
            # --- HEI identification logic ---
            # Choose the internal local maximum (exclude endpoints) with the highest energy,
            # i.e., nodes whose immediate neighbors have lower energy.
            # Fallback 1: if none exist, pick the maximum among internal nodes (exclude endpoints).
            # Fallback 2: if internal nodes are unavailable, pick the global maximum.
            nE = int(len(energies))
            hei_idx = None
            if nE >= 3:
                # Strict internal local maxima (both neighbors lower)
                candidates = [i for i in range(1, nE - 1)
                              if energies[i] > energies[i - 1] and energies[i] > energies[i + 1]]
                if candidates:
                    cand_es = energies[candidates]
                    rel = int(np.argmax(cand_es))
                    hei_idx = int(candidates[rel])
                else:
                    # Fallback 1: maximum over internal nodes (exclude endpoints)
                    if nE > 2:
                        rel = int(np.argmax(energies[1:-1]))
                        hei_idx = 1 + rel
            if hei_idx is None:
                # Fallback 2: global maximum
                hei_idx = int(np.argmax(energies))

            hei_geom = gs.images[hei_idx]
            hei_E = float(energies[hei_idx])

            hei_xyz = out_dir_path / "hei.xyz"
            s = hei_geom.as_xyz()
            lines = s.splitlines()
            if len(lines) >= 2 and lines[0].strip().isdigit():
                lines[1] = f"{hei_E:.12f}"
                s = "\n".join(lines) + ("\n" if not s.endswith("\n") else "")
            with open(hei_xyz, "w") as f:
                f.write(s)
            click.echo(f"[write] Wrote '{hei_xyz}'.")

            ref_pdb = None
            if source_paths[0].suffix.lower() == ".pdb":
                ref_pdb = source_paths[0].resolve()
            if ref_pdb is not None and is_convert_file_enabled():
                hei_pdb = out_dir_path / "hei.pdb"
                convert_xyz_to_pdb(hei_xyz, ref_pdb, hei_pdb)
                click.echo(f"[convert] Wrote '{hei_pdb}'.")
                # Annotate B-factors for ML & freeze atoms
                _apply_bfactor_annotations_inplace(
                    hei_pdb,
                    ml_indices=ml_indices_set,
                    freeze_indices=freeze_atoms_final,
                )
            else:
                click.echo("[convert] Skipped 'hei.pdb' (no PDB reference among inputs).")

        except Exception as e:
            click.echo(f"[HEI] ERROR: Failed to dump HEI: {e}", err=True)
            sys.exit(5)

        # summary.md and key_* outputs are disabled.
        click.echo(format_elapsed("[time] Elapsed Time for Path Opt", time_start), narrative=True)

        if out_json:
            from mlmm.core.utils import write_result_json
            from pysisyphus.constants import AU2KCALPERMOL as _AU2KCAL
            _gsm_energies = list(map(float, energies))
            _gsm_hei = int(hei_idx)
            _gsm_hei_E = float(_gsm_energies[_gsm_hei])
            _gsm_e0 = float(_gsm_energies[0])
            _gsm_eN = float(_gsm_energies[-1])
            _barrier = (_gsm_hei_E - _gsm_e0) * _AU2KCAL
            _delta = (_gsm_eN - _gsm_e0) * _AU2KCAL
            _converged = getattr(optimizer, 'is_converged', None) if 'optimizer' in dir() else None
            result_data_gsm: Dict[str, Any] = {
                "status": "converged" if _converged else ("not_converged" if _converged is False else "completed"),
                "converged": _converged,
                "mep_mode": "gsm",
                "backend": calc_cfg.get("backend", "uma"),
                "charge": calc_cfg.get("model_charge"),
                "spin": calc_cfg.get("model_mult"),
                "reactant_energy_hartree": float(_gsm_e0),
                "product_energy_hartree": float(_gsm_eN),
                "image_energies_hartree": [float(e) for e in _gsm_energies],
                "n_images": len(_gsm_energies),
                "hei_index": _gsm_hei,
                "hei_energy_hartree": _gsm_hei_E,
                "barrier_kcal": round(_barrier, 6),
                "delta_kcal": round(_delta, 6),
                "files": {
                    "final_geometries_trj_xyz": "final_geometries_trj.xyz",
                    "hei_xyz": "hei.xyz",
                },
            }
            for ext in (".pdb", ".gjf"):
                f = out_dir_path / f"hei{ext}"
                if f.exists():
                    result_data_gsm["files"][f"hei_{ext[1:]}"] = f.name
            write_result_json(
                out_dir_path, result_data_gsm,
                command="path-opt",
                elapsed_seconds=time.perf_counter() - time_start,
            )

    except OptimizationError as e:
        _write_error_json(Path(out_dir).resolve(), "path-opt", e, "OptimizationError", time_start)
        click.echo(f"ERROR: Path optimization failed — {e}", err=True)
        sys.exit(3)
    except KeyboardInterrupt:
        click.echo("\nInterrupted by user.", err=True)
        sys.exit(130)
    except Exception as e:
        render_cli_exception(e, label="path optimization", out_dir=out_dir, command="path-opt", time_start=time_start)
    finally:
        for prepared in prepared_inputs:
            prepared.cleanup()
        # Release GPU memory so subsequent pipeline stages don't OOM.
        # `= None` decref's the heavy refs; `del` then removes names from
        # the local frame so torch.nn.Module hooks / closures cannot retain.
        shared_calc = gs = geoms = None
        del shared_calc, gs, geoms
        gc.collect()  # break cyclic refs inside torch.nn.Module
        if torch.cuda.is_available():
            torch.cuda.empty_cache()


if __name__ == "__main__":
    cli()
