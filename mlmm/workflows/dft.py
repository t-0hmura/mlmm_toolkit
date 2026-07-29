"""ML/MM-aware single-point DFT for the ML region with energy recombination.

Example:
    mlmm dft -i enzyme.pdb --parm real.parm7 --model-pdb ml_region.pdb -q 0

For detailed documentation, see: docs/dft.md
"""

from __future__ import annotations

from copy import deepcopy
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence, Tuple

import logging
import sys
import tempfile
import time

logger = logging.getLogger(__name__)

import click
from mlmm.core.output import emit
import numpy as np
import yaml

from ase import Atoms

from pysisyphus.constants import AU2EV, AU2KCALPERMOL

from mlmm.backends.mlmm_calc import MLMMCore, mlmm as MLMMCalculator
from mlmm.workflows.opt import (
    GEOM_KW as OPT_GEOM_KW,
    CALC_KW as OPT_CALC_KW,
    _parse_freeze_atoms as _parse_freeze_atoms_opt,
    _normalize_geom_freeze as _normalize_geom_freeze_opt,
)
from mlmm.workflows.charge_prep import resolve_charge_spin_or_raise
from mlmm.core.utils import (
    apply_layer_freeze_constraints,
    apply_ref_pdb_override,
    apply_yaml_overrides,
    convert_xyz_to_pdb,
    pretty_block,
    format_freeze_atoms_for_echo,
    format_elapsed,
    merge_freeze_atom_indices,
    prepare_input_structure,
    parse_indices_string,
    resolve_ml_layer_assignment,
    set_convert_file_enabled,
)
from mlmm.cli.common_options import (
    add_allow_charge_mult_mismatch_option,
    add_ml_layer_detection_options,
)
from mlmm.cli.decorators import resolve_yaml_sources, load_merged_yaml_cfg, make_is_param_explicit, render_cli_exception
from mlmm.core.defaults import DFT_KW as _DFT_KW_DEFAULT
from mlmm.io.pdb_indexing import (
    PDBOrdinalAtom,
    parse_pdb_ordinal_atoms,
)

from functools import reduce

EV2AU = 1.0 / AU2EV

# Module-level alias (deepcopy at use-site for mutable safety)
DFT_KW = _DFT_KW_DEFAULT


@dataclass
class MLRegionWorkspace:
    calculator: MLMMCalculator
    high_level_backend: "_DFTEnergyOnlyBackend"
    input_pdb: Path
    real_parm7: Path
    real_rst7: Path
    model_pdb: Path
    model_parm7: Path
    model_rst7: Path
    selection_indices: List[int]
    link_pairs: List[Tuple[int, int]]  # 1-based REAL indices (ml_idx, mm_idx)
    atoms_real: Atoms
    atoms_model: Atoms
    atoms_model_lh: Atoms

    def cleanup(self) -> None:
        self.calculator.core.cleanup()

    @property
    def core(self) -> MLMMCore:
        return self.calculator.core


class _DFTEnergyOnlyBackend:
    """High-level slot used by MLMMCore for a precomputed DFT energy."""

    name = "dft"

    def __init__(self) -> None:
        self._energy_ev: Optional[float] = None

    def set_energy_hartree(self, energy_hartree: float) -> None:
        self._energy_ev = float(energy_hartree) * AU2EV

    def energy(self, atoms: Atoms) -> float:
        if self._energy_ev is None:
            raise RuntimeError("DFT energy has not been evaluated.")
        return self._energy_ev

    def eval(self, atoms: Atoms, need_grad: bool = True):
        raise RuntimeError(
            "The DFT high-level backend is energy-only; forces are not available."
        )


def _parse_func_basis(s: str) -> Tuple[str, str]:
    if not s or "/" not in s:
        raise click.BadParameter("Expected 'FUNC/BASIS' (e.g., 'wb97m-v/def2-tzvpd').")
    func, basis = s.split("/", 1)
    func = func.strip()
    basis = basis.strip()
    if not func or not basis:
        raise click.BadParameter("Functional or basis is empty. Example: --func-basis 'wb97m-v/6-31g**'")
    return func, basis


def _atoms_to_xyz_string(atoms: Atoms, comment: str) -> str:
    lines = [str(len(atoms)), comment]
    for sym, (x, y, z) in zip(atoms.get_chemical_symbols(), atoms.get_positions()):
        lines.append(f"{sym:<2s} {x:15.8f} {y:15.8f} {z:15.8f}")
    return "\n".join(lines) + "\n"


def write_ml_region_xyz_pair(
    workspace: "MLRegionWorkspace",
    out_dir: Path,
) -> Tuple[Path, Path]:
    """Write directly inspectable ML structures without and with generated link H."""

    destination = Path(out_dir)
    destination.mkdir(parents=True, exist_ok=True)
    without_link = destination / "ml_region_without_linkH.xyz"
    with_link = destination / "ml_region_with_linkH.xyz"
    without_link.write_text(
        _atoms_to_xyz_string(workspace.atoms_model, "ML region without link H"),
        encoding="utf-8",
    )
    with_link.write_text(
        _atoms_to_xyz_string(
            workspace.atoms_model_lh,
            "ML region with "
            f"{len(workspace.link_pairs)} link H; full-system ML-MM pairs: "
            + (
                ", ".join(
                    f"{ml_idx}-{mm_idx}"
                    for ml_idx, mm_idx in workspace.link_pairs
                )
                if workspace.link_pairs
                else "none"
            ),
        ),
        encoding="utf-8",
    )
    return without_link, with_link


def write_ml_region_pdb_pair(
    workspace: "MLRegionWorkspace",
    out_dir: Path,
    *,
    xyz_paths: Optional[Tuple[Path, Path]] = None,
) -> Tuple[Path, Path]:
    """Write topology-bearing PDB companions for the two ML snapshots."""

    without_xyz, with_xyz = (
        write_ml_region_xyz_pair(workspace, out_dir)
        if xyz_paths is None
        else xyz_paths
    )
    destination = Path(out_dir)
    without_pdb = destination / "ml_region_without_linkH.pdb"
    with_pdb = destination / "ml_region_with_linkH.pdb"

    convert_xyz_to_pdb(without_xyz, workspace.model_pdb, without_pdb)

    from mlmm.workflows.extract import _format_linkH_block, _max_serial_from_pdb_text

    link_coordinates = workspace.atoms_model_lh.get_positions()[
        len(workspace.atoms_model) :
    ]
    template_text = without_pdb.read_text(encoding="utf-8")
    link_template = workspace.model_pdb.parent / "model_with_linkH_template.pdb"
    link_template.write_text(
        template_text.rstrip()
        + "\n"
        + _format_linkH_block(
            [tuple(map(float, xyz)) for xyz in link_coordinates],
            _max_serial_from_pdb_text(template_text),
        ),
        encoding="utf-8",
    )
    convert_xyz_to_pdb(with_xyz, link_template, with_pdb)
    return without_pdb, with_pdb


def _atoms_to_pyscf_atoms(atoms: Atoms) -> List[Tuple[str, Tuple[float, float, float]]]:
    entries: List[Tuple[str, Tuple[float, float, float]]] = []
    for sym, coord in zip(atoms.get_chemical_symbols(), atoms.get_positions()):
        entries.append((sym, (float(coord[0]), float(coord[1]), float(coord[2]))))
    return entries


def _load_input_atoms(input_pdb: Path) -> List[PDBOrdinalAtom]:
    """Compatibility wrapper around the product-local ordinal parser."""

    return parse_pdb_ordinal_atoms(input_pdb)


def _prepare_ml_region_workspace(
    *,
    input_pdb: Path,
    coordinate_path: Optional[Path] = None,
    real_parm7: Path,
    model_pdb: Path,
    link_mlmm: Optional[Sequence[Sequence[str]]],
    link_atom_method: str = "scaled",
    use_cmap: bool = True,
    calc_kwargs: Optional[Dict[str, Any]] = None,
) -> MLRegionWorkspace:
    high_level_backend = _DFTEnergyOnlyBackend()
    calculator_kwargs = dict(calc_kwargs or {})
    calculator_kwargs.update(
        {
            "input_pdb": str(Path(input_pdb)),
            "coordinate_path": (
                None if coordinate_path is None else str(Path(coordinate_path))
            ),
            "real_parm7": str(Path(real_parm7)),
            "model_pdb": str(Path(model_pdb)),
            "link_mlmm": link_mlmm,
            "link_atom_method": str(link_atom_method or "scaled").lower(),
            "use_cmap": bool(use_cmap),
            "_high_level_backend": high_level_backend,
        }
    )
    calculator = MLMMCalculator(**calculator_kwargs)
    core = calculator.core
    try:
        atoms_real, atoms_model, atoms_model_lh = core.prepare_atoms()
    except Exception:
        core.cleanup()
        raise

    return MLRegionWorkspace(
        calculator=calculator,
        high_level_backend=high_level_backend,
        input_pdb=Path(core.input_pdb),
        real_parm7=Path(core.real_parm7),
        real_rst7=Path(core.real_rst7),
        model_pdb=Path(core.model_pdb),
        model_parm7=Path(core.model_parm7),
        model_rst7=Path(core.model_rst7),
        selection_indices=list(core.selection_indices),
        link_pairs=list(core.mlmm_links),
        atoms_real=atoms_real,
        atoms_model=atoms_model,
        atoms_model_lh=atoms_model_lh,
    )


def _hartree_to_kcalmol(Eh: float) -> float:
    return float(Eh * AU2KCALPERMOL)


def _build_dft_result_payload(
    *,
    converged: bool,
    energy_hartree: float,
    energy_kcal_per_mol: float,
    xc: str,
    basis: str,
    engine_label: str,
    using_gpu: bool,
    using_lowmem: bool,
    dft_kw: Dict[str, Any],
    calc_kw: Dict[str, Any],
    n_atoms: int,
    input_path: Path,
    charges: Dict[str, Any],
    spin_densities: Dict[str, Any],
) -> Dict[str, Any]:
    """Build the DFT JSON payload from the effective runtime request."""

    from mlmm.core.utils import calculator_provenance

    did_converge = bool(converged)
    return {
        "status": "converged" if did_converge else "not_converged",
        "converged": did_converge,
        "energy_hartree": energy_hartree,
        "energy_kcal_per_mol": energy_kcal_per_mol,
        "xc_functional": xc,
        "basis_set": basis,
        "engine": engine_label,
        "used_gpu": bool(using_gpu),
        "used_lowmem": bool(using_lowmem),
        **calculator_provenance(calc_kw),
        "charge": calc_kw.get("model_charge"),
        "spin": calc_kw.get("model_mult"),
        "n_atoms": int(n_atoms),
        "grid_level": dft_kw["grid_level"],
        "conv_tol": dft_kw["conv_tol"],
        "max_cycle": dft_kw["max_cycle"],
        "input_file": str(input_path),
        "charges": dict(charges),
        "spin_densities": dict(spin_densities),
        "files": {"result_yaml": "result.yaml"},
    }


def _apply_explicit_dft_overrides(
    dft_kw: Dict[str, Any],
    *,
    is_param_explicit,
    conv_tol: float,
    max_cycle: int,
    grid_level: int,
    out_dir: Path,
    lowmem: bool,
) -> Dict[str, Any]:
    """Apply only explicit CLI values over an already YAML-resolved DFT map."""

    resolved = dict(dft_kw)
    if is_param_explicit("conv_tol"):
        resolved["conv_tol"] = float(conv_tol)
    if is_param_explicit("max_cycle"):
        resolved["max_cycle"] = int(max_cycle)
    if is_param_explicit("grid_level"):
        resolved["grid_level"] = int(grid_level)
    if is_param_explicit("out_dir"):
        resolved["out_dir"] = str(out_dir)
    if is_param_explicit("lowmem"):
        resolved["lowmem"] = bool(lowmem)
    return resolved


def _finalize_dft_result(
    *,
    out_json: bool,
    out_dir: Path,
    payload: Dict[str, Any],
    elapsed_seconds: float,
) -> None:
    """Commit the result payload before signalling SCF nonconvergence."""

    if out_json:
        from mlmm.core.utils import write_result_json

        write_result_json(
            out_dir,
            payload,
            command="dft",
            elapsed_seconds=elapsed_seconds,
        )
    if not bool(payload["converged"]):
        raise SystemExit(3)


def _prepare_dft_output_dir(path: Path) -> Path:
    """Create the output directory and invalidate prior public DFT results."""

    resolved = Path(path).resolve()
    resolved.mkdir(parents=True, exist_ok=True)
    for name in (
        "result.yaml",
        "result.json",
        "summary.json",
        "ml_region_without_linkH.xyz",
        "ml_region_with_linkH.xyz",
        "ml_region_without_linkH.pdb",
        "ml_region_with_linkH.pdb",
    ):
        (resolved / name).unlink(missing_ok=True)
    return resolved


class FlowList(list):
    pass


def _flow_seq_representer(dumper, data):
    return dumper.represent_sequence('tag:yaml.org,2002:seq', data, flow_style=True)


yaml.SafeDumper.add_representer(FlowList, _flow_seq_representer)


def _format_row_for_echo(row: List[Any]) -> str:
    def _fmt(x: Any) -> str:
        if x is None:
            return "null"
        if isinstance(x, float):
            return f"{x:.10g}"
        return str(x)

    return "[" + ", ".join(_fmt(v) for v in row) + "]"


# PySCF population-analysis helpers
def fast_iao_mullikan_spin_pop(mol, dm, iaos, verbose=None):
    import numpy
    from pyscf.lib import logger as pyscf_logger
    from pyscf.lo.iao import reference_mol
    from pyscf.scf import uhf as scf_uhf

    if verbose is None:
        verbose = pyscf_logger.DEBUG

    pmol = reference_mol(mol)
    if getattr(mol, 'pbc_intor', None):
        ovlpS = mol.pbc_intor('int1e_ovlp')
    else:
        ovlpS = mol.intor_symmetric('int1e_ovlp')

    cs = numpy.dot(iaos.T.conj(), ovlpS)
    s_iao = numpy.dot(cs, iaos)
    iao_inv = numpy.linalg.solve(s_iao, cs)

    if isinstance(dm, numpy.ndarray) and dm.ndim == 2:
        spin_pop_ao = numpy.zeros(s_iao.shape[0], dtype=float)
        Ms = numpy.zeros(pmol.natm, dtype=float)
        return spin_pop_ao, Ms

    dm_a = reduce(numpy.dot, (iao_inv, dm[0], iao_inv.conj().T))
    dm_b = reduce(numpy.dot, (iao_inv, dm[1], iao_inv.conj().T))
    return scf_uhf.mulliken_spin_pop(pmol, [dm_a, dm_b], s_iao, verbose)


def _build_cpu_surrogate_for_analysis(mol, mf, xc: str, spin2s: int):
    """Construct a CPU `pyscf.dft` SCF object carrying converged orbitals from `mf`.

    Required when the GPU SCF object's `to_cpu()` is not implemented
    (e.g. `gpu4pyscf.dft.rks_lowmem.RKS`). MOs are transferred from CuPy if
    needed so that `pyscf.scf.hf` / `pyscf.scf.uhf` population analysis
    routines can operate on plain numpy arrays.
    """
    from pyscf import dft as pdft

    def _to_np(x):
        try:
            import cupy as _cp
            if isinstance(x, _cp.ndarray):
                return _cp.asnumpy(x)
        except Exception:
            pass
        return np.asarray(x)

    surrogate = pdft.RKS(mol) if spin2s == 0 else pdft.UKS(mol)
    surrogate.xc = xc
    surrogate.chkfile = None
    mo_coeff = mf.mo_coeff
    mo_occ = mf.mo_occ
    mo_energy = getattr(mf, "mo_energy", None)
    if spin2s == 0:
        surrogate.mo_coeff = _to_np(mo_coeff)
        surrogate.mo_occ = _to_np(mo_occ)
        if mo_energy is not None:
            surrogate.mo_energy = _to_np(mo_energy)
    else:
        surrogate.mo_coeff = (_to_np(mo_coeff[0]), _to_np(mo_coeff[1]))
        surrogate.mo_occ = (_to_np(mo_occ[0]), _to_np(mo_occ[1]))
        if mo_energy is not None:
            surrogate.mo_energy = (_to_np(mo_energy[0]), _to_np(mo_energy[1]))
    surrogate.converged = bool(getattr(mf, "converged", False))
    surrogate.e_tot = float(getattr(mf, "e_tot", float("nan")))
    return surrogate


def _compute_atomic_charges(mol, mf) -> Dict[str, Optional[List[float]]]:
    from pyscf.scf import hf as scf_hf
    from pyscf.lo import iao as lo_iao

    dm = mf.make_rdm1()
    S = mf.get_ovlp()
    dm_tot = dm[0] + dm[1] if (isinstance(dm, np.ndarray) and dm.ndim == 3) else dm

    try:
        _, mull_chg = scf_hf.mulliken_pop(mol, dm_tot, s=S, verbose=0)
        mull_q = np.asarray(mull_chg, dtype=float).tolist()
    except Exception as e:
        click.echo(f"[Mulliken] WARNING: Failed to compute Mulliken charges: {e}", err=True)
        mull_q = None

    try:
        _, low_chg = scf_hf.mulliken_pop_meta_lowdin_ao(mol, dm_tot, verbose=0, s=S)
        low_q = np.asarray(low_chg, dtype=float).tolist()
    except Exception as e:
        click.echo(f"[Löwdin] WARNING: Failed to compute meta-Löwdin charges: {e}", err=True)
        low_q = None

    iao_q: Optional[List[float]] = None
    try:
        mo = mf.mo_coeff
        mo_occ = mf.mo_occ
        if isinstance(mo, np.ndarray) and mo.ndim == 2:
            occ_idx = np.asarray(mo_occ) > 0
            orbocc = mo[:, occ_idx]
        else:
            occ_idx = np.asarray(mo_occ[0]) > 0
            orbocc = mo[0][:, occ_idx]
        iaos = lo_iao.iao(mol, orbocc, minao="minao")
        _, iao_chg = lo_iao.fast_iao_mullikan_pop(mol, dm, iaos, verbose=0)
        iao_q = np.asarray(iao_chg, dtype=float).tolist()
    except Exception as e:
        click.echo(f"[IAO] WARNING: Failed to compute IAO charges: {e}", err=True)
        iao_q = None

    return {
        "mulliken": mull_q,
        "lowdin": low_q,
        "iao": iao_q,
    }


def _compute_atomic_spin_densities(mol, mf) -> Dict[str, Optional[List[float]]]:
    from pyscf.scf import uhf as scf_uhf
    from pyscf.lo import iao as lo_iao

    dm = mf.make_rdm1()
    S = mf.get_ovlp()
    nat = mol.natm

    if not (isinstance(dm, np.ndarray) and dm.ndim == 3):
        zeros = [0.0] * nat
        return {"mulliken": zeros, "lowdin": zeros, "iao": zeros}

    try:
        _, Ms_mull = scf_uhf.mulliken_spin_pop(mol, dm, s=S, verbose=0)
        mull = np.asarray(Ms_mull, dtype=float).tolist()
    except Exception as e:
        click.echo(f"[Spin Mulliken] WARNING: Failed to compute Mulliken spin densities: {e}", err=True)
        mull = None

    try:
        _, Ms_low = scf_uhf.mulliken_spin_pop_meta_lowdin_ao(mol, dm, verbose=0, s=S)
        low = np.asarray(Ms_low, dtype=float).tolist()
    except Exception as e:
        click.echo(f"[Spin Löwdin] WARNING: Failed to compute meta-Löwdin spin densities: {e}", err=True)
        low = None

    iao_ms: Optional[List[float]] = None
    try:
        mo = mf.mo_coeff
        mo_occ = mf.mo_occ
        if isinstance(mo, np.ndarray) and mo.ndim == 2:
            occ_idx = np.asarray(mo_occ) > 0
            orbocc = mo[:, occ_idx]
        else:
            occ_idx = np.asarray(mo_occ[0]) > 0
            orbocc = mo[0][:, occ_idx]
        iaos = lo_iao.iao(mol, orbocc, minao="minao")
        _, Ms_iao = fast_iao_mullikan_spin_pop(mol, dm, iaos, verbose=0)
        iao_ms = np.asarray(Ms_iao, dtype=float).tolist()
    except Exception as e:
        click.echo(f"[Spin IAO] WARNING: Failed to compute IAO spin densities: {e}", err=True)
        iao_ms = None

    return {"mulliken": mull, "lowdin": low, "iao": iao_ms}


@click.command(
    help="Single-point ML-region DFT with ML(dft)/MM energy recombination.",
    context_settings={"help_option_names": ["-h", "--help"]},
)
@click.option(
    "-i",
    "--input",
    "input_path",
    type=click.Path(path_type=Path, exists=True, dir_okay=False),
    required=True,
    help="Full enzyme structure (PDB or XYZ). If XYZ, use --ref-pdb for topology.",
)
@click.option(
    "--ref-pdb",
    "ref_pdb",
    type=click.Path(path_type=Path, exists=True, dir_okay=False),
    default=None,
    show_default=False,
    help="Reference PDB topology when input is XYZ. XYZ coordinates are used (higher precision) "
         "while PDB provides atom ordering and residue information.",
)
@click.option(
    "--parm",
    "real_parm7",
    type=click.Path(path_type=Path, exists=True, dir_okay=False),
    required=True,
    help="Amber parm7 topology for the full system.",
)
@click.option(
    "--model-pdb",
    "model_pdb",
    type=click.Path(path_type=Path, exists=True, dir_okay=False),
    required=False,
    help="ML-only, link-H-free PDB subset; atom identity/order must match the "
         "full PDB/parm7. Optional when --detect-layer is enabled.",
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
    "-q",
    "--charge",
    type=int,
    required=False,
    help=(
        "ML region charge. Required unless --ligand-charge is provided "
        "(PDB inputs or XYZ with --ref-pdb)."
    ),
)
@click.option(
    "-l",
    "--ligand-charge",
    type=str,
    default=None,
    show_default=False,
    help=(
        "Total charge or per-resname mapping (e.g., 'SAM:1,GPP:-3') used to derive "
        "ML region charge when -q is omitted (requires PDB input or --ref-pdb)."
    ),
)
@click.option(
    "-m",
    "--multiplicity",
    "spin",
    type=int,
    default=None,
    show_default=False,
    help="Spin multiplicity (2S+1) for the ML region; defaults to YAML or 1.",
)
@click.option(
    "--freeze-atoms",
    "freeze_atoms_text",
    type=str,
    default=None,
    help="Comma-separated 1-based indices to freeze (e.g., '1,3,5').",
)
@click.option(
    "--func-basis",
    "func_basis",
    type=str,
    default="wb97m-v/def2-tzvpd",
    show_default=True,
    help='Exchange-correlation functional and basis set as "FUNC/BASIS".',
)
@click.option("--max-cycle", type=int, default=DFT_KW["max_cycle"], show_default=True, help="Maximum SCF iterations.")
@click.option("--conv-tol", type=float, default=DFT_KW["conv_tol"], show_default=True, help="SCF energy convergence threshold (ΔE in Hartree between SCF cycles).")
@click.option("--grid-level", type=int, default=DFT_KW["grid_level"], show_default=True, help="DFT integration grid level (0=coarse, 3=default, 5=fine, 9=very fine).")
@click.option(
    "--engine",
    type=click.Choice(["gpu", "cpu"], case_sensitive=False),
    default="gpu",
    show_default=True,
    help="SCF backend: gpu (GPU4PySCF, raises error if unavailable) or cpu (PySCF).",
)
@click.option(
    "--lowmem/--no-lowmem",
    "lowmem",
    default=DFT_KW["lowmem"],
    show_default=True,
    help="Use gpu4pyscf rks_lowmem.RKS for closed-shell GPU runs "
         "(memory-efficient direct JK; mlmm dft does not call density_fit() on either path). "
         "Open-shell or CPU engines fall back to standard RKS/UKS automatically.",
)
@click.option(
    "-o", "--out-dir",
    type=click.Path(path_type=Path, dir_okay=True, file_okay=False),
    default=Path(DFT_KW["out_dir"]),
    show_default=True,
    help="Output directory.",
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
    help="Validate options and print the execution plan without running DFT.",
)
@click.option(
    "--convert-files/--no-convert-files",
    "convert_files",
    default=True,
    show_default=True,
    help="Toggle XYZ/TRJ to PDB companions when a PDB template is available.",
)
@click.option(
    "-b", "--backend",
    type=click.Choice(["uma", "orb", "mace", "aimnet2"], case_sensitive=False),
    default=None,
    show_default=False,
    help="Backend label recorded in output metadata; the ML region in dft is computed with DFT (PySCF), so this does not select a calculator.",
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
    help="Link-atom placement: 'scaled' (g-factor, Gaussian ONIOM standard, default) or "
         "'fixed' (legacy 1.09 Å for C, 1.01 Å for N).",
)
@click.option(
    "--mm-backend",
    "mm_backend",
    type=click.Choice(["hessian_ff", "openmm"], case_sensitive=False),
    default=None,
    show_default=False,
    help="MM backend for the low-level ONIOM evaluation: 'hessian_ff' (default) or 'openmm'.",
)
@click.option(
    "--cmap/--no-cmap",
    "use_cmap",
    default=None,
    show_default=False,
    help="Preserve CMAP terms in both real and model MM layers. Default: enabled when present in parm7.",
)
@click.option(
    "--out-json/--no-out-json",
    "out_json",
    default=False,
    show_default=True,
    help="Write machine-readable result.json to out_dir.",
)
@add_allow_charge_mult_mismatch_option()
@add_ml_layer_detection_options()
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
    charge: Optional[int],
    ligand_charge: Optional[str],
    spin: Optional[int],
    freeze_atoms_text: Optional[str],
    func_basis: str,
    max_cycle: int,
    conv_tol: float,
    grid_level: int,
    engine: str,
    lowmem: bool,
    out_dir: Path,
    config_yaml: Optional[Path],
    show_config: bool,
    dry_run: bool,
    convert_files: bool,
    backend: Optional[str],
    embedcharge: bool,
    embedcharge_cutoff: Optional[float],
    link_atom_method: Optional[str],
    mm_backend: Optional[str],
    use_cmap: Optional[bool],
    out_json: bool,
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

    prepared_input = None
    workspace: Optional[MLRegionWorkspace] = None
    # time_start must be bound BEFORE the try-block (same pattern as sp.py).
    # Otherwise any exception raised between the `try:` and the original
    # assignment site reaches the `except Exception as exc:` branch with
    # time_start unbound, masking the real error with UnboundLocalError.
    time_start: float = time.perf_counter()
    out_dir_path = Path(out_dir).resolve()

    model_indices: Optional[List[int]] = None
    if model_indices_str:
        try:
            model_indices = parse_indices_string(model_indices_str, one_based=model_indices_one_based)
        except click.BadParameter as e:
            raise click.ClickException(str(e))

    try:
        # One prepared input owns the complete DFT lifetime.  Reference
        # topology is overlaid onto this same object, so layer selection and
        # ligand-charge derivation cannot observe different structures.
        prepared_input = prepare_input_structure(input_path)
        if input_path.suffix.lower() == ".xyz":
            if ref_pdb is None:
                raise click.BadParameter(
                    "Provide --ref-pdb for topology when using XYZ input."
                )
            apply_ref_pdb_override(prepared_input, ref_pdb)
        source_pdb = prepared_input.source_path

        geom_kw = deepcopy(OPT_GEOM_KW)
        calc_kw = deepcopy(OPT_CALC_KW)
        dft_kw = dict(DFT_KW)

        apply_yaml_overrides(
            merged_yaml_cfg,
            [
                (geom_kw, (("geom",),)),
                (calc_kw, (("calc",), ("mlmm",))),
                (dft_kw, (("dft",),)),
            ],
        )

        # CLI explicit overrides (after config YAML)
        if backend is not None:
            calc_kw["backend"] = str(backend).lower()
        if _is_param_explicit("embedcharge"):
            calc_kw["embedcharge"] = bool(embedcharge)
        if _is_param_explicit("embedcharge_cutoff"):
            calc_kw["embedcharge_cutoff"] = embedcharge_cutoff
        if link_atom_method is not None:
            calc_kw["link_atom_method"] = str(link_atom_method).lower()
        if mm_backend is not None:
            calc_kw["mm_backend"] = str(mm_backend).lower()
        if use_cmap is not None:
            calc_kw["use_cmap"] = use_cmap

        dft_kw = _apply_explicit_dft_overrides(
            dft_kw,
            is_param_explicit=_is_param_explicit,
            conv_tol=conv_tol,
            max_cycle=max_cycle,
            grid_level=grid_level,
            out_dir=out_dir,
            lowmem=lowmem,
        )
        out_dir_path = Path(dft_kw["out_dir"]).resolve()

        func_basis_value = str(dft_kw.get("func_basis", func_basis))
        if _is_param_explicit("func_basis"):
            func_basis_value = func_basis
        xc, basis = _parse_func_basis(func_basis_value)

        geom_kw["coord_type"] = "cart"
        geom_kw["freeze_atoms"] = _normalize_geom_freeze_opt(geom_kw.get("freeze_atoms"))
        freeze_atoms_cli = _parse_freeze_atoms_opt(freeze_atoms_text)
        calc_kw["freeze_atoms"] = merge_freeze_atom_indices(geom_kw, freeze_atoms_cli)

        detect_layer_enabled = bool(calc_kw.get("use_bfactor_layers", True))
        if _is_param_explicit("detect_layer"):
            detect_layer_enabled = bool(detect_layer)

        layer_source_pdb = source_pdb.resolve()
        model_pdb_cfg = calc_kw.get("model_pdb")
        if model_pdb is not None:
            model_pdb_cfg = model_pdb
        model_pdb_path: Optional[Path] = None
        layer_info: Optional[Dict[str, List[int]]] = None

        # Hoist charge/spin resolution above any `int(charge)` use so that
        # `--ligand-charge`-only invocations (-q omitted) do not crash with
        # TypeError(None→int). resolve_charge_spin_or_raise handles both
        # charge=None + ligand_charge=... (derives charge from PDB residues)
        # and charge=None + ligand_charge=None (raises a clean ClickException
        # "ML-region charge is unresolved" instead of a TypeError).
        charge, spin = resolve_charge_spin_or_raise(
            prepared_input, charge, spin,
            ligand_charge=ligand_charge, prefix="[dft]",
            model_pdb=model_pdb,
            model_indices_spec=model_indices_str,
            detect_layer=detect_layer,
            yaml_cfg=merged_yaml_cfg,
        )
        calc_kw["model_mult"] = int(spin)
        calc_kw["model_charge"] = int(charge)
        engine_name = (
            str(engine).strip().lower()
            if _is_param_explicit("engine")
            else str(dft_kw.get("engine", "gpu")).strip().lower()
        )
        if engine_name not in {"cpu", "gpu"}:
            raise click.BadParameter(
                "dft.engine must be either 'cpu' or 'gpu', "
                f"got {engine_name!r}."
            )
        dft_kw["engine"] = engine_name

        dft_block = {
            "charge": int(charge),
            "multiplicity": int(calc_kw["model_mult"]),
            "xc": xc,
            "basis": basis,
            "conv_tol": dft_kw["conv_tol"],
            "max_cycle": dft_kw["max_cycle"],
            "grid_level": dft_kw["grid_level"],
            "out_dir": str(Path(dft_kw["out_dir"]).resolve()),
            "engine": engine_name,
            "lowmem": bool(dft_kw.get("lowmem", True)),
        }
        from mlmm.core.embedcharge_policy import reject_retired_embedcharge_cli

        reject_retired_embedcharge_cli(
            calc_kw,
            cutoff_requested=_is_param_explicit("embedcharge_cutoff"),
        )

        click.echo(pretty_block("geom", format_freeze_atoms_for_echo(geom_kw, key="freeze_atoms")))
        click.echo(pretty_block("calc", {k: calc_kw[k] for k in sorted(calc_kw.keys()) if k not in {"freeze_atoms"}}))
        click.echo(pretty_block("dft", dft_block))

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
            with tempfile.TemporaryDirectory(prefix="mlmm_dft_validate_") as tmp:
                validation_cfg = dict(calc_kw)
                validation_model, _ = resolve_ml_layer_assignment(
                    source_path=layer_source_pdb,
                    out_dir_path=Path(tmp),
                    model_pdb=model_pdb_cfg,
                    model_indices=model_indices,
                    detect_layer=detect_layer_enabled,
                    hess_cutoff=validation_cfg.get("hess_cutoff"),
                    movable_cutoff=validation_cfg.get("movable_cutoff"),
                    calc_cfg=validation_cfg,
                    echo_fn=click.echo,
                )
                validation_workspace = _prepare_ml_region_workspace(
                    input_pdb=source_pdb,
                    coordinate_path=prepared_input.geom_path,
                    real_parm7=real_parm7,
                    model_pdb=validation_model,
                    link_mlmm=validation_cfg.get("link_mlmm"),
                    link_atom_method=str(
                        validation_cfg.get("link_atom_method") or "scaled"
                    ).lower(),
                    use_cmap=bool(calc_kw.get("use_cmap", True)),
                    calc_kwargs=validation_cfg,
                )
                validation_workspace.cleanup()
                click.echo(
                    pretty_block(
                        "dry_run_plan",
                        {
                            "will_prepare_input": True,
                            "detect_layer": bool(detect_layer_enabled),
                            "model_region_source": (
                                "model_pdb"
                                if model_pdb_cfg is not None
                                else (
                                    "model_indices"
                                    if model_indices
                                    else "bfactor"
                                )
                            ),
                            "model_indices_count": (
                                0 if not model_indices else len(model_indices)
                            ),
                            "output_dir": str(Path(dft_kw["out_dir"]).resolve()),
                            "backend": calc_kw.get("backend", "uma"),
                            "embedcharge": bool(
                                calc_kw.get("embedcharge", False)
                            ),
                        },
                    )
                )
                click.echo(
                    "[dry-run] Validation complete. DFT execution was skipped."
                )
                return

        # `prepared_input` and `charge/spin` already resolved above
        # (hoisted to fix `int(charge)` TypeError with --ligand-charge only).
        out_dir_path = _prepare_dft_output_dir(out_dir_path)

        model_pdb_path, layer_info = resolve_ml_layer_assignment(
            source_path=layer_source_pdb,
            out_dir_path=out_dir_path,
            model_pdb=model_pdb_cfg,
            model_indices=model_indices,
            detect_layer=detect_layer_enabled,
            hess_cutoff=calc_kw.get("hess_cutoff"),
            movable_cutoff=calc_kw.get("movable_cutoff"),
            calc_cfg=calc_kw,
            echo_fn=click.echo,
        )

        calc_kw["input_pdb"] = str(source_pdb.resolve())
        calc_kw["real_parm7"] = str(real_parm7.resolve())
        calc_kw["model_pdb"] = str(model_pdb_path.resolve())
        apply_layer_freeze_constraints(
            geom_kw,
            calc_kw,
            layer_info,
            echo_fn=click.echo,
        )
        # Re-baseline so the wall-clock figure excludes YAML resolution
        # and ML-region preparation overhead.
        time_start = time.perf_counter()

        workspace = _prepare_ml_region_workspace(
            input_pdb=Path(calc_kw["input_pdb"]),
            coordinate_path=prepared_input.geom_path,
            real_parm7=Path(calc_kw["real_parm7"]),
            model_pdb=Path(calc_kw["model_pdb"]),
            link_mlmm=calc_kw.get("link_mlmm"),
            link_atom_method=str(calc_kw.get("link_atom_method") or "scaled").lower(),
            use_cmap=bool(calc_kw.get("use_cmap", True)),
            calc_kwargs=calc_kw,
        )
        model_charge = int(calc_kw["model_charge"])
        model_mult = int(calc_kw["model_mult"])
        model_spin2s = model_mult - 1

        without_link_path, with_link_path = write_ml_region_xyz_pair(
            workspace,
            out_dir_path,
        )
        click.echo(f"[write] Wrote '{without_link_path}'.")
        click.echo(f"[write] Wrote '{with_link_path}'.")
        if Path(input_path).suffix.lower() == ".pdb":
            without_link_pdb, with_link_pdb = write_ml_region_pdb_pair(
                workspace,
                out_dir_path,
                xyz_paths=(without_link_path, with_link_path),
            )
            click.echo(f"[write] Wrote '{without_link_pdb}'.")
            click.echo(f"[write] Wrote '{with_link_pdb}'.")

        try:
            from pyscf import gto
        except Exception as exc:
            raise click.ClickException(f"PySCF import failed: {exc}") from exc

        from mlmm.core.utils import is_verbose
        mol = gto.Mole()
        # PySCF verbose level: 0 silent / 2 warnings / 3 info / 4 prints the
        # full `[INPUT]` per-atom coordinate dump + SCF iteration banner.
        # Default runs stay quiet (DFT_KW seeds a low value); `-v` raises the
        # level to >=4 to restore the dump for debugging (~600 lines/pipeline).
        mol.verbose = int(dft_kw.pop("verbose", 0))
        if is_verbose():
            mol.verbose = max(mol.verbose, 4)
        # CHEMISTRY-RULE:5 def2 family auto-ECP injection.
        # def2 family includes Stuttgart ECPs for heavy elements (Z>=21);
        # must set mol.ecp explicitly or PySCF uses all-electron treatment.
        _ecp = dft_kw.get("ecp", None)
        if _ecp is None and basis.lower().startswith("def2"):
            _ecp = basis
        _build_kw: Dict[str, Any] = dict(
            atom=_atoms_to_pyscf_atoms(workspace.atoms_model_lh),
            unit="Angstrom",
            charge=model_charge,
            spin=model_spin2s,
            basis=basis,
        )
        if _ecp:
            _build_kw["ecp"] = _ecp
            click.echo(f"[dft] Using ECP: {_ecp}")
        mol.build(**_build_kw)

        using_gpu = False
        using_lowmem = False
        engine_label = "pyscf(cpu)"
        _engine = engine_name
        lowmem_requested = bool(dft_kw.get("lowmem", True))
        if _engine != "cpu":
            try:
                from gpu4pyscf import dft as gdf

                # CHEMISTRY-RULE:4 gpu4pyscf rks_lowmem triple-guard
                # (lowmem flag + closed-shell + module-import success).
                rks_lowmem_mod = None
                if lowmem_requested and model_spin2s == 0:
                    try:
                        from gpu4pyscf.dft import rks_lowmem as rks_lowmem_mod  # type: ignore
                    except ImportError:
                        click.echo(
                            "[lowmem] WARNING: gpu4pyscf.dft.rks_lowmem is not available "
                            "in this gpu4pyscf install; falling back to standard RKS.",
                            err=True,
                        )
                        rks_lowmem_mod = None

                if rks_lowmem_mod is not None:
                    mf = rks_lowmem_mod.RKS(mol, xc=xc)
                    using_lowmem = True
                    engine_label = "gpu4pyscf(rks_lowmem)"
                else:
                    if lowmem_requested and model_spin2s != 0:
                        click.echo(
                            "[lowmem] NOTE: rks_lowmem only supports closed-shell; "
                            "open-shell run uses standard UKS.",
                            err=True,
                        )
                    mf = gdf.RKS(mol) if model_spin2s == 0 else gdf.UKS(mol)
                    engine_label = "gpu4pyscf"
                using_gpu = True
            except Exception as e:
                raise click.ClickException(
                    f"[gpu] GPU backend failed: {e}. "
                    "Set dft.engine: cpu in YAML config to explicitly run on CPU."
                )
        if _engine == "cpu":
            from pyscf import dft as pdft

            mf = pdft.RKS(mol) if model_spin2s == 0 else pdft.UKS(mol)

        mf.xc = xc
        mf.max_cycle = int(dft_kw["max_cycle"])
        mf.conv_tol = float(dft_kw["conv_tol"])
        try:
            mf.grids.level = int(dft_kw["grid_level"])
        except Exception as exc:
            click.echo(f"[grids] WARNING: Could not set grids.level={dft_kw['grid_level']}: {exc}", err=True)
        try:
            mf.chkfile = None
        except Exception:
            logger.debug("Failed to disable chkfile", exc_info=True)
        if xc.lower().endswith("-v") or "vv10" in xc.lower():
            mf.nlc = "vv10"

        # --- Electrostatic embedding (--embedcharge) ---
        n_mm_charges = 0
        if calc_kw.get("embedcharge", False):
            import parmed as pmd
            from pyscf import qmmm as pyscf_qmmm

            real_top = pmd.load_file(str(workspace.real_parm7))
            ml_set = set(workspace.selection_indices)
            mm_indices = [i for i in range(len(workspace.atoms_real)) if i not in ml_set]

            _cutoff = calc_kw.get("embedcharge_cutoff", None)
            if mm_indices and _cutoff is not None:
                from scipy.spatial.distance import cdist
                _ml_ref = workspace.atoms_real.get_positions()[sorted(ml_set)]
                _mm_all = workspace.atoms_real.get_positions()[mm_indices]
                _dists = cdist(_mm_all, _ml_ref).min(axis=1)
                mm_indices = [mm_indices[j] for j in range(len(mm_indices)) if _dists[j] <= _cutoff]

            if mm_indices:
                mm_coords = workspace.atoms_real.get_positions()[mm_indices]
                mm_charges = np.array([real_top.atoms[i].charge for i in mm_indices])
                mf = pyscf_qmmm.mm_charge(mf, mm_coords, mm_charges, unit="Angstrom")
                n_mm_charges = len(mm_indices)
                click.echo(f"[embedcharge] {n_mm_charges} MM point charges embedded into QM Hamiltonian.")
            else:
                click.echo("[embedcharge] No MM atoms found; skipping embedding.")

        tic_scf = time.time()
        e_tot = mf.kernel()
        toc_scf = time.time()

        converged = bool(getattr(mf, "converged", False))
        if e_tot is None:
            e_tot = float(getattr(mf, "e_tot", np.nan))
        e_h = float(e_tot)
        e_kcal = _hartree_to_kcalmol(e_h)

        if using_lowmem:
            mf_for_analysis = _build_cpu_surrogate_for_analysis(mol, mf, xc, model_spin2s)
        else:
            try:
                mf_for_analysis = mf.to_cpu()
            except Exception:
                mf_for_analysis = mf
        charges = _compute_atomic_charges(mol, mf_for_analysis)
        spins = _compute_atomic_spin_densities(mol, mf_for_analysis)

        def _round(xs: Optional[List[float]]) -> Optional[List[float]]:
            if xs is None:
                return None
            return [0.0 if (x == x) and abs(x) < 1e-10 else float(x) for x in xs]

        charges = {k: _round(v) for k, v in charges.items()}
        spins = {k: _round(v) for k, v in spins.items()}

        charges_table: List[List[Any]] = []
        spins_table: List[List[Any]] = []
        for i in range(mol.natm):
            elem = mol.atom_symbol(i)
            charges_table.append([
                i,
                elem,
                None if charges["mulliken"] is None else charges["mulliken"][i],
                None if charges["lowdin"] is None else charges["lowdin"][i],
                None if charges["iao"] is None else charges["iao"][i],
            ])
            spins_table.append([
                i,
                elem,
                None if spins["mulliken"] is None else spins["mulliken"][i],
                None if spins["lowdin"] is None else spins["lowdin"][i],
                None if spins["iao"] is None else spins["iao"][i],
            ])

        click.echo("\ncharges [index, element, mulliken, lowdin, iao]:")
        for row in charges_table:
            click.echo(f"- {_format_row_for_echo(row)}")

        # Skip the per-atom spin_densities echo for closed-shell systems —
        # spin=0 forces every row to all-zero and the resulting 77 zero rows
        # contribute nothing but noise. The data still lands in result.yaml.
        if model_spin2s != 0:
            click.echo("\nspin_densities [index, element, mulliken, lowdin, iao]:")
            for row in spins_table:
                click.echo(f"- {_format_row_for_echo(row)}")

        workspace.high_level_backend.set_energy_hartree(e_h)
        mlmm_result = workspace.core.compute(
            workspace.atoms_real.get_positions(),
            return_forces=False,
            return_hessian=False,
        )
        energy_components = mlmm_result["energy_components"]
        e_real_low = float(energy_components["real_low"])
        e_model_low = float(energy_components["model_low"])

        e_real_low_au = e_real_low * EV2AU
        e_model_low_au = e_model_low * EV2AU
        e_total_au = float(mlmm_result["energy"]) * EV2AU
        e_total_kcal = _hartree_to_kcalmol(e_total_au)

        result_yaml = {
            "input": {
                "charge": model_charge,
                "multiplicity": model_mult,
                "xc": xc,
                "basis": basis,
                "conv_tol": dft_kw["conv_tol"],
                "max_cycle": dft_kw["max_cycle"],
                "grid_level": dft_kw["grid_level"],
                "out_dir": str(out_dir_path),
                "embedcharge": bool(calc_kw.get("embedcharge", False)),
                "n_mm_charges": n_mm_charges,
            },
            "energy": {
                "hartree": e_h,
                "kcal_per_mol": e_kcal,
                "converged": converged,
                "scf_time_sec": round(toc_scf - tic_scf, 3),
                "engine": engine_label,
                "used_gpu": bool(using_gpu),
                "used_lowmem": bool(using_lowmem),
            },
            "mlmm_energy": {
                "E_real_low_eV": e_real_low,
                "E_model_low_eV": e_model_low,
                "E_real_low_hartree": e_real_low_au,
                "E_model_low_hartree": e_model_low_au,
                "E_total_ml_dft_mm_hartree": e_total_au,
                "E_total_ml_dft_mm_kcal_per_mol": e_total_kcal,
            },
            "charges [index, element, mulliken, lowdin, iao]": [FlowList(r) for r in charges_table],
            "spin_densities [index, element, mulliken, lowdin, iao]": [FlowList(r) for r in spins_table],
        }

        result_file = out_dir_path / "result.yaml"
        result_file.write_text(yaml.safe_dump(result_yaml, sort_keys=False, allow_unicode=True))
        click.echo(f"[write] Wrote '{result_file}'.")
        # summary.md and key_* outputs are disabled.
        click.echo(f"\nE_DFT (Hartree): {e_h:.12f}")
        click.echo(f"E_DFT (kcal/mol): {e_kcal:.6f}")
        click.echo(f"E_total ML(dft)/MM (Hartree): {e_total_au:.12f}")
        click.echo(f"E_total ML(dft)/MM (kcal/mol): {e_total_kcal:.6f}")

        if not converged:
            click.echo("WARNING: SCF did not converge.", err=True)

        emit(format_elapsed("[time] Elapsed Time for DFT", time_start), narrative=True)

        result_data = _build_dft_result_payload(
            converged=converged,
            energy_hartree=e_h,
            energy_kcal_per_mol=e_kcal,
            xc=xc,
            basis=basis,
            engine_label=engine_label,
            using_gpu=using_gpu,
            using_lowmem=using_lowmem,
            dft_kw=dft_kw,
            calc_kw=calc_kw,
            n_atoms=mol.natm,
            input_path=input_path,
            charges=charges,
            spin_densities=spins,
        )
        _finalize_dft_result(
            out_json=out_json,
            out_dir=out_dir_path,
            payload=result_data,
            elapsed_seconds=time.perf_counter() - time_start,
        )

    except KeyboardInterrupt:
        click.echo("\nInterrupted by user.", err=True)
        sys.exit(130)
    except click.ClickException:
        raise
    except Exception as exc:
        render_cli_exception(
            exc,
            label="ML/MM DFT",
            out_dir=out_dir_path,
            command="dft",
            time_start=time_start,
        )
    finally:
        if prepared_input is not None:
            prepared_input.cleanup()
        if workspace is not None:
            workspace.cleanup()


if __name__ == "__main__":
    cli()
