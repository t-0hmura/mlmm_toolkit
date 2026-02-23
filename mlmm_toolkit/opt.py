# mlmm_toolkit/opt.py

"""
ML/MM geometry optimization (LBFGS or RFO) with UMA + hessian_ff calculator.

Example:
    mlmm opt -i pocket.pdb --real-parm7 real.parm7 --model-pdb ml_region.pdb -q 0

For detailed documentation, see: docs/opt.md
"""

from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence, Tuple, Set

import ast
import os
import shutil
import sys
import textwrap
import traceback

import click
import numpy as np
import yaml
import time
from click.core import ParameterSource

from pysisyphus.helpers import geom_loader
from pysisyphus.optimizers.LBFGS import LBFGS
from pysisyphus.optimizers.RFOptimizer import RFOptimizer
from pysisyphus.optimizers.exceptions import OptimizationError, ZeroStepLength
from pysisyphus.constants import ANG2BOHR, BOHR2ANG, AU2EV

from .mlmm_calc import mlmm
from .defaults import (
    GEOM_KW_DEFAULT,
    MLMM_CALC_KW,
    OPT_BASE_KW,
    LBFGS_KW,
    RFO_KW,
    LAYEROPT_KW,
    OPT_MODE_ALIASES,
    BFACTOR_ML,
    BFACTOR_MOVABLE_MM,
    BFACTOR_FROZEN,
)
from .utils import (
    append_xyz_trajectory as _append_xyz_trajectory,
    convert_xyz_to_pdb,
    deep_update,
    load_yaml_dict,
    apply_yaml_overrides,
    pretty_block,
    strip_inherited_keys,
    format_freeze_atoms_for_echo,
    format_elapsed,
    merge_freeze_atom_indices,
    prepare_input_structure,
    apply_ref_pdb_override,
    resolve_charge_spin_or_raise,
    parse_indices_string,
    build_model_pdb_from_bfactors,
    build_model_pdb_from_indices,
    update_pdb_bfactors_from_layers,
    normalize_choice,
)

EV2AU = 1.0 / AU2EV                 # eV → Hartree
H_EVAA_2_AU = EV2AU / (ANG2BOHR * ANG2BOHR)  # (eV/Å^2) → (Hartree/Bohr^2)


def _first_existing_artifact(out_dir: Path, patterns: Sequence[str]) -> Optional[Path]:
    """Resolve the first existing artifact for a list of relative patterns."""
    for pattern in patterns:
        if any(ch in pattern for ch in "*?[]"):
            for candidate in sorted(out_dir.glob(pattern)):
                if candidate.is_file():
                    return candidate.resolve()
            continue
        candidate = out_dir / pattern
        if candidate.is_file():
            return candidate.resolve()
    return None


def _link_or_copy_file(src: Path, dst: Path) -> bool:
    """Create a symlink when possible; fall back to copy."""
    try:
        if dst.exists() or dst.is_symlink():
            if dst.is_dir():
                return False
            dst.unlink()
        rel = os.path.relpath(src, start=dst.parent)
        dst.symlink_to(rel)
        return True
    except Exception:
        try:
            shutil.copy2(src, dst)
            return True
        except Exception:
            return False


def _write_output_summary_md(out_dir: Path) -> None:
    """Write summary.md and expose key outputs at out_dir root."""
    try:
        out_dir = out_dir.resolve()
        if not out_dir.exists():
            return

        root_specs: List[Tuple[str, Sequence[str]]] = [
            ("Optimized geometry (XYZ)", ["final_geometry.xyz"]),
            ("Optimized geometry (PDB)", ["final_geometry.pdb"]),
            ("Optimization trajectory", ["optimization_all_trj.xyz", "optimization_trj.xyz"]),
            ("Optimization trajectory (PDB)", ["optimization_all.pdb", "optimization.pdb"]),
            ("Restart snapshot", ["restart*.yml"]),
        ]
        root_lines: List[str] = []
        for label, patterns in root_specs:
            src = _first_existing_artifact(out_dir, patterns)
            if src is None:
                continue
            rel = os.path.relpath(src, start=out_dir)
            root_lines.append(f"- {label}: [`{rel}`]({rel})")

        shortcut_specs: List[Tuple[str, str, Sequence[str]]] = [
            ("key_opt.xyz", "Optimized geometry (XYZ)", ["final_geometry.xyz"]),
            ("key_opt.pdb", "Optimized geometry (PDB)", ["final_geometry.pdb"]),
            ("key_opt_trj.xyz", "Optimization trajectory", ["optimization_all_trj.xyz", "optimization_trj.xyz"]),
            ("key_opt_traj.pdb", "Optimization trajectory (PDB)", ["optimization_all.pdb", "optimization.pdb"]),
            ("key_restart.yml", "Restart snapshot", ["restart*.yml"]),
        ]
        shortcut_lines: List[str] = []
        for name, label, patterns in shortcut_specs:
            src = _first_existing_artifact(out_dir, patterns)
            if src is None:
                continue
            dst = out_dir / name
            try:
                same = src.resolve() == dst.resolve()
            except Exception:
                same = False
            if not same and not _link_or_copy_file(src, dst):
                continue
            src_rel = os.path.relpath(src, start=out_dir)
            shortcut_lines.append(
                f"- {label}: [`{name}`]({name}) (source: `{src_rel}`)"
            )

        lines: List[str] = [
            "# Opt Summary",
            "",
            f"- Generated: `{time.strftime('%Y-%m-%d %H:%M:%S %Z')}`",
            f"- Output directory: `{out_dir}`",
            "",
            "## Primary Artifacts",
        ]
        if root_lines:
            lines.extend(root_lines)
        else:
            lines.append("- No primary artifacts detected yet.")
        lines.extend(
            [
                "",
                "## Root Shortcuts",
                "- `key_*` files are symlinks when possible, otherwise copied files.",
            ]
        )
        if shortcut_lines:
            lines.extend(shortcut_lines)
        else:
            lines.append("- No shortcuts were generated.")
        lines.extend(
            [
                "",
                "## Notes",
                "- Start from `key_opt.xyz` (or `key_opt.pdb`) and inspect `key_opt_trj.xyz` when available.",
            ]
        )

        summary_md = out_dir / "summary.md"
        summary_md.write_text("\n".join(lines).rstrip() + "\n", encoding="utf-8")
        click.echo(f"[write] Wrote '{summary_md}'.")
    except Exception as e:
        click.echo(f"[write] WARNING: Failed to write summary.md: {e}", err=True)


def _resolve_yaml_sources(
    config_yaml: Optional[Path],
    override_yaml: Optional[Path],
    args_yaml_legacy: Optional[Path],
) -> Tuple[Optional[Path], Optional[Path], bool]:
    if override_yaml is not None and args_yaml_legacy is not None:
        raise click.BadParameter(
            "Use a single YAML source option."
        )
    if args_yaml_legacy is not None:
        return config_yaml, args_yaml_legacy, True
    return config_yaml, override_yaml, False


def _load_merged_yaml_cfg(
    config_yaml: Optional[Path],
    override_yaml: Optional[Path],
) -> Dict[str, Any]:
    merged: Dict[str, Any] = {}
    deep_update(merged, load_yaml_dict(config_yaml))
    deep_update(merged, load_yaml_dict(override_yaml))
    return merged


# -----------------------------------------------
# Default settings (imported from defaults.py, aliased for compatibility)
# -----------------------------------------------

GEOM_KW: Dict[str, Any] = dict(GEOM_KW_DEFAULT)
CALC_KW: Dict[str, Any] = dict(MLMM_CALC_KW)

# Note: OPT_BASE_KW, LBFGS_KW, RFO_KW are imported from defaults.py


# -----------------------------------------------
# MicroIteration optimizer for LayerOpt-style optimization
# -----------------------------------------------

class MicroIterationOptimizer:
    """
    Microiteration optimizer implementing a LayerOpt-style pattern:
    1. Fix ML region, optimize outer region with LBFGS
    2. Fix outer region, optimize ML region with RFO
    3. Repeat until convergence

    This provides a simplified LayerOpt implementation suitable for ML/MM
    calculations where the ML region indices are known.
    """

    def __init__(
        self,
        geometry,
        ml_indices: List[int],
        outer_lbfgs_kwargs: Dict[str, Any],
        inner_rfo_kwargs: Dict[str, Any],
        max_macro_cycles: int = 100,
        outer_max_cycles: int = 1500,
        outer_thresh: str = "gau_loose",
        echo_fn=None,
        dump: bool = False,
    ):
        """
        Parameters
        ----------
        geometry : pysisyphus Geometry
            Geometry with calculator attached
        ml_indices : List[int]
            0-based indices of ML region atoms
        outer_lbfgs_kwargs : dict
            LBFGS kwargs for outer optimization
        inner_rfo_kwargs : dict
            RFO kwargs for inner optimization
        max_macro_cycles : int
            Maximum macrocycles (outer+inner pairs)
        outer_max_cycles : int
            Max LBFGS cycles per outer optimization
        outer_thresh : str
            Convergence threshold for outer optimization
        echo_fn : callable
            Function to echo messages (default: click.echo)
        """
        self.geometry = geometry
        self.ml_indices = sorted(set(int(i) for i in ml_indices))
        self.n_atoms = len(geometry.atoms)
        self.outer_indices = [i for i in range(self.n_atoms) if i not in self.ml_indices]

        self.outer_lbfgs_kwargs = dict(outer_lbfgs_kwargs)
        self.inner_rfo_kwargs = dict(inner_rfo_kwargs)
        self.max_macro_cycles = max_macro_cycles
        self.outer_max_cycles = outer_max_cycles
        self.outer_thresh = outer_thresh
        self.echo = echo_fn or click.echo
        self.dump = bool(dump)
        self.out_dir = Path(self.inner_rfo_kwargs.get("out_dir", "./result_opt/"))
        self.optim_all_path = self.out_dir / "optimization_all_trj.xyz"

        self.is_converged = False
        self.cur_cycle = 0
        self.micro_cycles = []

    def run(self) -> None:
        """Run microiteration optimization."""
        from pysisyphus.Geometry import Geometry

        # Store original freeze_atoms
        calc = self.geometry.calculator
        base_calc = calc.base if hasattr(calc, 'base') else calc
        original_freeze = list(getattr(base_calc.core, 'freeze_atoms', []) or [])

        self.echo("\n=== MicroIteration Optimization ===\n")
        self.echo(f"ML region: {len(self.ml_indices)} atoms")
        self.echo(f"Outer region: {len(self.outer_indices)} atoms")
        if self.dump and self.optim_all_path.exists():
            self.optim_all_path.unlink()

        for macro in range(self.max_macro_cycles):
            self.cur_cycle = macro
            self.echo(f"\n--- Macrocycle {macro + 1}/{self.max_macro_cycles} ---\n")

            # Step 1: Optimize outer region (freeze ML region)
            self.echo("[MicroIter] Optimizing outer region (ML frozen)...")

            # Set freeze atoms to ML region + original freeze
            outer_freeze = sorted(set(self.ml_indices + original_freeze))
            base_calc.core.freeze_atoms = outer_freeze

            # Create fresh geometry for outer optimization
            outer_geom = Geometry(
                atoms=self.geometry.atoms,
                coords=self.geometry.coords.copy(),
                coord_type="cart",
                freeze_atoms=outer_freeze,
            )
            outer_geom.set_calculator(calc)

            outer_kwargs = dict(self.outer_lbfgs_kwargs)
            outer_kwargs["max_cycles"] = self.outer_max_cycles
            outer_kwargs["thresh"] = self.outer_thresh
            outer_kwargs["prefix"] = f"mc{macro:03d}_outer"
            outer_kwargs["dump"] = self.dump

            outer_opt = LBFGS(outer_geom, **outer_kwargs)
            outer_opt.run()
            outer_cycles = outer_opt.cur_cycle + 1
            if self.dump:
                _append_xyz_trajectory(
                    self.optim_all_path,
                    outer_opt.get_path_for_fn("optimization_trj.xyz"),
                )

            # Update main geometry with outer-optimized coords
            self.geometry.coords = outer_geom.coords.copy()

            # Step 2: Optimize ML region (freeze outer region)
            self.echo("[MicroIter] Optimizing ML region (outer frozen)...")

            # Set freeze atoms to outer region + original freeze
            inner_freeze = sorted(set(self.outer_indices + original_freeze))
            base_calc.core.freeze_atoms = inner_freeze

            # Create fresh geometry for inner optimization
            inner_geom = Geometry(
                atoms=self.geometry.atoms,
                coords=self.geometry.coords.copy(),
                coord_type="cart",
                freeze_atoms=inner_freeze,
            )
            inner_geom.set_calculator(calc)

            inner_kwargs = dict(self.inner_rfo_kwargs)
            inner_kwargs["max_cycles"] = 1  # Single RFO step per macrocycle
            inner_kwargs["prefix"] = f"mc{macro:03d}_inner"
            inner_kwargs["dump"] = self.dump

            inner_opt = RFOptimizer(inner_geom, **inner_kwargs)
            inner_opt.run()
            inner_cycles = inner_opt.cur_cycle + 1
            if self.dump:
                _append_xyz_trajectory(
                    self.optim_all_path,
                    inner_opt.get_path_for_fn("optimization_trj.xyz"),
                )

            # Update main geometry with inner-optimized coords
            self.geometry.coords = inner_geom.coords.copy()

            self.micro_cycles.append((outer_cycles, inner_cycles))
            self.echo(f"[MicroIter] Macrocycle {macro + 1}: outer={outer_cycles}, inner={inner_cycles}")

            # Check convergence using inner optimizer's convergence
            if inner_opt.is_converged:
                self.is_converged = True
                self.echo("\n[MicroIter] Converged!")
                break

        # Restore original freeze_atoms
        base_calc.core.freeze_atoms = original_freeze

        total_outer = sum(o for o, _ in self.micro_cycles)
        total_inner = sum(i for _, i in self.micro_cycles)
        self.echo(f"\n[MicroIter] Total cycles: outer={total_outer}, inner={total_inner}")
        self.echo(f"[MicroIter] Macrocycles: {self.cur_cycle + 1}, converged={self.is_converged}")

    @property
    def final_fn(self) -> Path:
        """Return final geometry path (compatible with standard optimizer interface)."""
        return self.out_dir / "final_geometry.xyz"

    def get_path_for_fn(self, fn: str) -> Path:
        """Get path for output files (compatible interface)."""
        return self.out_dir / fn


class PartialHessianMicroIterationOptimizer:
    """
    Microiteration optimizer for partial-Hessian workflows:
    1) Freeze Hessian-target atoms, optimize non-Hessian movable atoms (LBFGS).
    2) Freeze non-Hessian movable atoms, optimize Hessian-target atoms (RFO).
    3) Repeat outer->inner cycles until a fixed-point criterion is met or macrocycle limit.
    """

    def __init__(
        self,
        geometry,
        hess_indices: List[int],
        outer_indices: List[int],
        outer_lbfgs_kwargs: Dict[str, Any],
        inner_rfo_kwargs: Dict[str, Any],
        max_macro_cycles: int = 100,
        outer_max_cycles: int = 1500,
        inner_max_cycles: Optional[int] = None,
        macro_disp_thresh_bohr: float = 5e-4,
        outer_thresh: str = "gau_loose",
        echo_fn=None,
        inner_opt_cls=RFOptimizer,
        dump: bool = False,
    ):
        self.geometry = geometry
        self.hess_indices = sorted({int(i) for i in hess_indices})
        self.outer_indices = sorted({int(i) for i in outer_indices if i not in set(self.hess_indices)})
        self.n_atoms = len(geometry.atoms)

        self.outer_lbfgs_kwargs = dict(outer_lbfgs_kwargs)
        self.inner_rfo_kwargs = dict(inner_rfo_kwargs)
        self.max_macro_cycles = max_macro_cycles
        self.outer_max_cycles = outer_max_cycles
        self.inner_max_cycles = int(inner_max_cycles) if inner_max_cycles is not None else None
        self.macro_disp_thresh_bohr = float(macro_disp_thresh_bohr)
        self.outer_thresh = outer_thresh
        self.echo = echo_fn or click.echo
        self.inner_opt_cls = inner_opt_cls
        self.dump = bool(dump)
        self.out_dir = Path(self.inner_rfo_kwargs.get("out_dir", "./result_opt/"))
        self.optim_all_path = self.out_dir / "optimization_all_trj.xyz"

        self.is_converged = False
        self.cur_cycle = 0
        self.micro_cycles = []

    def run(self) -> None:
        from pysisyphus.Geometry import Geometry

        geom_calc = self.geometry.calculator
        base_calc = geom_calc.base if hasattr(geom_calc, "base") else geom_calc
        if hasattr(base_calc, "freeze_atoms"):
            original_freeze = list(base_calc.freeze_atoms or [])
        elif hasattr(base_calc, "core"):
            original_freeze = list(getattr(base_calc.core, "freeze_atoms", []) or [])
        else:
            original_freeze = []

        self.echo("\n=== Partial-Hessian MicroIteration Optimization ===\n")
        self.echo(f"Hessian-target atoms: {len(self.hess_indices)}")
        self.echo(f"Non-Hessian movable atoms: {len(self.outer_indices)}")
        resolved_inner_max_cycles = self.inner_max_cycles
        if resolved_inner_max_cycles is None:
            try:
                resolved_inner_max_cycles = int(self.inner_rfo_kwargs.get("max_cycles", 1))
            except Exception:
                resolved_inner_max_cycles = 1
        resolved_inner_max_cycles = max(1, int(resolved_inner_max_cycles))
        self.echo(f"Inner RFO max_cycles per macro: {resolved_inner_max_cycles}")
        self.echo(f"Macro fixed-point disp threshold: {self.macro_disp_thresh_bohr:.2e} bohr")
        if self.dump and self.optim_all_path.exists():
            self.optim_all_path.unlink()

        for macro in range(self.max_macro_cycles):
            self.cur_cycle = macro
            self.echo(f"\n--- Macrocycle {macro + 1}/{self.max_macro_cycles} ---\n")
            macro_start_coords = np.array(self.geometry.coords, copy=True)

            outer_cycles = 0
            inner_cycles = 0
            outer_converged = not self.outer_indices

            # Step 1: Optimize non-Hessian movable atoms (freeze Hessian-target atoms)
            if self.outer_indices:
                self.echo("[PartialMicroIter] Optimizing non-Hessian movable atoms (Hessian atoms frozen)...")
                outer_freeze = sorted(set(self.hess_indices + original_freeze))
                if hasattr(base_calc, "freeze_atoms"):
                    base_calc.freeze_atoms = list(outer_freeze)
                elif hasattr(base_calc, "core"):
                    base_calc.core.freeze_atoms = list(outer_freeze)
                    if hasattr(base_calc.core, "_update_active_dof_mappings"):
                        base_calc.core._update_active_dof_mappings()

                outer_geom = Geometry(
                    atoms=self.geometry.atoms,
                    coords=self.geometry.coords.copy(),
                    coord_type="cart",
                    freeze_atoms=outer_freeze,
                )
                outer_geom.set_calculator(self.geometry.calculator)

                outer_kwargs = dict(self.outer_lbfgs_kwargs)
                outer_kwargs["max_cycles"] = self.outer_max_cycles
                outer_kwargs["thresh"] = self.outer_thresh
                outer_kwargs["prefix"] = f"mc{macro:03d}_outer"
                outer_kwargs["dump"] = self.dump

                outer_opt = LBFGS(outer_geom, **outer_kwargs)
                outer_opt.run()
                outer_cycles = outer_opt.cur_cycle + 1
                outer_converged = bool(getattr(outer_opt, "is_converged", False))
                if self.dump:
                    _append_xyz_trajectory(
                        self.optim_all_path,
                        outer_opt.get_path_for_fn("optimization_trj.xyz"),
                    )

                # Update main geometry
                self.geometry.coords = outer_geom.coords.copy()

            # Step 2: Optimize Hessian-target atoms (freeze non-Hessian movable atoms)
            self.echo("[PartialMicroIter] Optimizing Hessian-target atoms (non-Hessian movable frozen)...")
            inner_freeze = sorted(set(self.outer_indices + original_freeze))
            if hasattr(base_calc, "freeze_atoms"):
                base_calc.freeze_atoms = list(inner_freeze)
            elif hasattr(base_calc, "core"):
                base_calc.core.freeze_atoms = list(inner_freeze)
                if hasattr(base_calc.core, "_update_active_dof_mappings"):
                    base_calc.core._update_active_dof_mappings()

            inner_geom = Geometry(
                atoms=self.geometry.atoms,
                coords=self.geometry.coords.copy(),
                coord_type="cart",
                freeze_atoms=inner_freeze,
            )
            inner_geom.set_calculator(self.geometry.calculator)

            inner_kwargs = dict(self.inner_rfo_kwargs)
            inner_kwargs["max_cycles"] = resolved_inner_max_cycles
            inner_kwargs["prefix"] = f"mc{macro:03d}_inner"
            inner_kwargs["dump"] = self.dump

            inner_opt = self.inner_opt_cls(inner_geom, **inner_kwargs)
            inner_opt.run()
            inner_cycles = inner_opt.cur_cycle + 1
            inner_converged = bool(getattr(inner_opt, "is_converged", False))
            if self.dump:
                _append_xyz_trajectory(
                    self.optim_all_path,
                    inner_opt.get_path_for_fn("optimization_trj.xyz"),
                )

            # Update main geometry
            self.geometry.coords = inner_geom.coords.copy()
            macro_disp = float(
                np.max(np.abs(np.asarray(self.geometry.coords) - np.asarray(macro_start_coords)))
            )
            macro_fixed_point = macro_disp <= self.macro_disp_thresh_bohr

            self.micro_cycles.append((outer_cycles, inner_cycles))
            self.echo(
                f"[PartialMicroIter] Macrocycle {macro + 1}: "
                f"outer={outer_cycles} (conv={outer_converged}), "
                f"inner={inner_cycles} (conv={inner_converged}), "
                f"max|Δx|={macro_disp:.3e} bohr (fixed_point={macro_fixed_point})"
            )

            if outer_converged and inner_converged and macro_fixed_point:
                self.is_converged = True
                self.echo("\n[PartialMicroIter] Converged!")
                break

        # Restore original freeze_atoms
        if hasattr(base_calc, "freeze_atoms"):
            base_calc.freeze_atoms = list(original_freeze)
        elif hasattr(base_calc, "core"):
            base_calc.core.freeze_atoms = list(original_freeze)
            if hasattr(base_calc.core, "_update_active_dof_mappings"):
                base_calc.core._update_active_dof_mappings()

        total_outer = sum(o for o, _ in self.micro_cycles)
        total_inner = sum(i for _, i in self.micro_cycles)
        self.echo(f"\n[PartialMicroIter] Total cycles: outer={total_outer}, inner={total_inner}")
        self.echo(f"[PartialMicroIter] Macrocycles: {self.cur_cycle + 1}, converged={self.is_converged}")

    @property
    def final_fn(self) -> Path:
        return self.out_dir / "final_geometry.xyz"

    def get_path_for_fn(self, fn: str) -> Path:
        return self.out_dir / fn


class HarmonicBiasCalculator:
    """Wrap a base UMA calculator with harmonic distance restraints."""

    def __init__(self, base_calc, k: float = 10.0, pairs: Optional[List[Tuple[int, int, float]]] = None):
        self.base = base_calc
        self.k_evAA = float(k)
        self.k_au_bohr2 = self.k_evAA * H_EVAA_2_AU
        self._pairs: List[Tuple[int, int, float]] = list(pairs or [])

    def set_pairs(self, pairs: List[Tuple[int, int, float]]) -> None:
        self._pairs = [(int(i), int(j), float(t)) for (i, j, t) in pairs]

    def _bias_energy_forces_bohr(self, coords_bohr: np.ndarray) -> Tuple[float, np.ndarray]:
        coords = np.array(coords_bohr, dtype=float).reshape(-1, 3)
        n = coords.shape[0]
        E_bias = 0.0
        F_bias = np.zeros((n, 3), dtype=float)
        k = self.k_au_bohr2
        for (i, j, target_ang) in self._pairs:
            if not (0 <= i < n and 0 <= j < n):
                continue
            rij_vec = coords[i] - coords[j]
            rij = float(np.linalg.norm(rij_vec))
            if rij < 1e-14:
                continue
            target_bohr = float(target_ang) * ANG2BOHR
            diff_bohr = rij - target_bohr
            E_bias += 0.5 * k * diff_bohr * diff_bohr
            u = rij_vec / max(rij, 1e-14)
            Fi = -k * diff_bohr * u
            F_bias[i] += Fi
            F_bias[j] -= Fi
        return E_bias, F_bias.reshape(-1)

    def get_forces(self, elem, coords):
        coords_bohr = np.asarray(coords, dtype=float).reshape(-1, 3)
        base = self.base.get_forces(elem, coords_bohr)
        E0 = float(base["energy"])
        F0 = np.asarray(base["forces"], dtype=float).reshape(-1)
        Ebias, Fbias = self._bias_energy_forces_bohr(coords_bohr)
        return {"energy": E0 + Ebias, "forces": F0 + Fbias}

    def get_energy(self, elem, coords):
        coords_bohr = np.asarray(coords, dtype=float).reshape(-1, 3)
        E0 = float(self.base.get_energy(elem, coords_bohr)["energy"])
        Ebias, _ = self._bias_energy_forces_bohr(coords_bohr)
        return {"energy": E0 + Ebias}

    def get_energy_and_forces(self, elem, coords):
        res = self.get_forces(elem, coords)
        return res["energy"], res["forces"]

    def get_energy_and_gradient(self, elem, coords):
        res = self.get_forces(elem, coords)
        return res["energy"], -np.asarray(res["forces"], dtype=float).reshape(-1)

    def __getattr__(self, name: str):
        return getattr(self.base, name)


def _parse_freeze_atoms(arg: Optional[str]) -> List[int]:
    """Parse comma-separated 1-based indices (e.g., "1,3,5") into 0-based ints."""
    if arg is None:
        return []

    items = [chunk.strip() for chunk in str(arg).split(",")]
    indices: List[int] = []
    for idx, chunk in enumerate(items, start=1):
        if not chunk:
            continue
        try:
            value = int(chunk)
        except ValueError as exc:
            raise click.BadParameter(
                f"Invalid integer in --freeze-atoms entry #{idx}: '{chunk}'"
            ) from exc
        if value <= 0:
            raise click.BadParameter(
                f"--freeze-atoms expects 1-based positive indices; got {value}"
            )
        indices.append(value - 1)
    return sorted(set(indices))


def _normalize_geom_freeze(value: Any) -> List[int]:
    """Normalize YAML-provided geom.freeze_atoms to a sorted 0-based list."""
    if value is None:
        return []
    if isinstance(value, str):
        tokens = [tok.strip() for tok in value.split(",") if tok.strip()]
        try:
            return sorted({int(tok) for tok in tokens})
        except ValueError as exc:
            raise click.BadParameter(
                "geom.freeze_atoms must contain integers (string form)."
            ) from exc
    try:
        return sorted({int(idx) for idx in value})
    except TypeError as exc:
        raise click.BadParameter("geom.freeze_atoms must be iterable of integers.") from exc


def _parse_dist_freeze(
    args: Sequence[str],
    one_based: bool,
) -> List[Tuple[int, int, Optional[float]]]:
    """Parse --dist-freeze arguments into 0-based pairs with optional targets."""
    parsed: List[Tuple[int, int, Optional[float]]] = []
    for idx, raw in enumerate(args, start=1):
        try:
            obj = ast.literal_eval(raw)
        except Exception as e:
            raise click.BadParameter(f"Invalid literal for --dist-freeze #{idx}: {e}")
        if isinstance(obj, (list, tuple)) and not obj:
            iterable = []
        elif isinstance(obj, (list, tuple)) and isinstance(obj[0], (list, tuple)):
            iterable = obj
        else:
            iterable = [obj]
        for entry in iterable:
            if not (isinstance(entry, (list, tuple)) and len(entry) in (2, 3)):
                raise click.BadParameter(
                    f"--dist-freeze #{idx} entries must be (i,j) or (i,j,target_A): {entry}"
                )
            if not (
                isinstance(entry[0], (int, np.integer))
                and isinstance(entry[1], (int, np.integer))
            ):
                raise click.BadParameter(f"Atom indices in --dist-freeze #{idx} must be integers: {entry}")
            i = int(entry[0])
            j = int(entry[1])
            target = None
            if len(entry) == 3:
                if not isinstance(entry[2], (int, float, np.floating)):
                    raise click.BadParameter(f"Target distance must be numeric in --dist-freeze #{idx}: {entry}")
                target = float(entry[2])
                if target <= 0.0:
                    raise click.BadParameter(f"Target distance must be > 0 in --dist-freeze #{idx}: {entry}")
            if one_based:
                i -= 1
                j -= 1
            if i < 0 or j < 0:
                raise click.BadParameter(f"--dist-freeze #{idx} produced negative index after conversion: {entry}")
            parsed.append((i, j, target))
    return parsed


def _resolve_dist_freeze_targets(
    geometry,
    tuples: List[Tuple[int, int, Optional[float]]],
) -> List[Tuple[int, int, float]]:
    coords_bohr = np.array(geometry.coords3d, dtype=float).reshape(-1, 3)
    coords_ang = coords_bohr * BOHR2ANG
    n = coords_ang.shape[0]
    resolved: List[Tuple[int, int, float]] = []
    for (i, j, target) in tuples:
        if not (0 <= i < n and 0 <= j < n):
            raise click.BadParameter(
                f"--dist-freeze indices {(i, j)} are out of bounds for the loaded geometry (N={n})."
            )
        if target is None:
            vec = coords_ang[i] - coords_ang[j]
            dist = float(np.linalg.norm(vec))
        else:
            dist = float(target)
        resolved.append((i, j, dist))
    return resolved


# -------------------------------
# PDB helpers for B-factor patch
# -------------------------------

def _pdb_keys_from_line(line: str) -> Tuple[Tuple, Tuple]:
    """
    Extract robust keys from a PDB ATOM/HETATM record.

    Returns:
        key_full: (chain, resseq, icode, resname, atomname, altloc)
        key_simple: (chain, resseq, icode, atomname)
    """
    atom_name = line[12:16].strip()
    altloc = line[16:17].strip()
    resname = line[17:20].strip()
    chain = line[21:22].strip()
    resseq_str = line[22:26].strip()
    try:
        resseq = int(resseq_str)
    except ValueError:
        resseq = -10**9  # unlikely sentinel when missing
    icode = line[26:27].strip()
    key_full = (chain, resseq, icode, resname, atom_name, altloc)
    key_simple = (chain, resseq, icode, atom_name)
    return key_full, key_simple


def _collect_ml_atom_keys(model_pdb: Path) -> Tuple[Set[Tuple], Set[Tuple]]:
    """Collect ML-region atom keys from model_pdb."""
    keys_full: Set[Tuple] = set()
    keys_simple: Set[Tuple] = set()
    try:
        with model_pdb.open("r") as fh:
            for line in fh:
                if line.startswith("ATOM") or line.startswith("HETATM"):
                    kf, ks = _pdb_keys_from_line(line)
                    keys_full.add(kf)
                    keys_simple.add(ks)
    except Exception:
        # If anything goes wrong, leave sets empty; caller will handle gracefully.
        pass
    return keys_full, keys_simple


def _format_with_bfactor(line: str, b: float) -> str:
    """Return PDB line with B-factor field (cols 61-66) set to b (6.2f)."""
    if len(line) < 66:
        line = line.rstrip("\n")
        line = line + " " * max(0, 66 - len(line))
        line = line + "\n"
    bf_str = f"{b:6.2f}"
    # Preserve occupancy (cols 55-60), overwrite tempFactor (61-66).
    new_line = line[:60] + bf_str + line[66:]
    return new_line


def _annotate_b_factors_inplace(
    pdb_path: Path,
    model_pdb: Path,
    freeze_indices_0based: Sequence[int],
    beta_ml: float = 100.0,
    beta_frz: float = 50.0,
    beta_both: float = 150.0,
) -> None:
    """
    Overwrite B-factors in-place:
      - ML-region atoms: 100.00
      - frozen atoms: 50.00
      - ML ∩ frozen: 150.00
    Indexing for 'frozen' is 0-based and resets at each MODEL.
    """
    ml_full, ml_simple = _collect_ml_atom_keys(model_pdb)
    frozen_set = set(int(i) for i in (freeze_indices_0based or []))

    try:
        lines = pdb_path.read_text().splitlines(keepends=True)
    except Exception:
        return

    out_lines: List[str] = []
    atom_idx = 0  # resets per MODEL

    for line in lines:
        rec = line[:6]
        if rec.startswith("MODEL"):
            # reset atom counter for each model
            atom_idx = 0
            out_lines.append(line)
            continue
        if rec.startswith("ATOM  ") or rec.startswith("HETATM"):
            kf, ks = _pdb_keys_from_line(line)
            is_ml = (kf in ml_full) or (ks in ml_simple)
            is_frz = (atom_idx in frozen_set)
            if is_ml and is_frz:
                out_lines.append(_format_with_bfactor(line, beta_both))
            elif is_ml:
                out_lines.append(_format_with_bfactor(line, beta_ml))
            elif is_frz:
                out_lines.append(_format_with_bfactor(line, beta_frz))
            else:
                out_lines.append(line)
            atom_idx += 1
        else:
            out_lines.append(line)

    try:
        pdb_path.write_text("".join(out_lines))
    except Exception:
        # Silently ignore if we cannot write; conversion outputs are still present.
        pass


def _maybe_convert_outputs_to_pdb(
    input_path: Path,
    out_dir: Path,
    dump: bool,
    get_trj_fn,
    final_xyz_path: Path,
    model_pdb: Path,
    freeze_indices_0based: Sequence[int],
    ml_indices: Optional[List[int]] = None,
    hess_mm_indices: Optional[List[int]] = None,
    movable_mm_indices: Optional[List[int]] = None,
    frozen_layer_indices: Optional[List[int]] = None,
) -> None:
    """
    If the input is a PDB, convert outputs (final_geometry.xyz and, if dump, optimization_all_trj.xyz /
    optimization_trj.xyz) to PDB,
    and annotate B-factors for the 3-layer ML/MM system.

    B-factor encoding (3-layer system):
        ML atoms: 0.0
        Movable MM atoms: 10.0
        Frozen MM atoms: 20.0

    If layer indices are not provided, falls back to legacy encoding:
        ML atoms: 100.0
        Frozen atoms: 50.0
        ML ∩ frozen: 150.0
    """
    if input_path.suffix.lower() != ".pdb":
        return

    # Determine if we should use the layer-based B-factor encoding
    use_layer_bfactors = ml_indices is not None

    ref_pdb = input_path.resolve()
    # final_geometry.xyz → final_geometry.pdb
    final_pdb = out_dir / "final_geometry.pdb"
    try:
        convert_xyz_to_pdb(final_xyz_path, ref_pdb, final_pdb)
        click.echo(f"[convert] Wrote '{final_pdb}'.")

        if use_layer_bfactors:
            update_pdb_bfactors_from_layers(
                final_pdb,
                ml_indices=ml_indices or [],
                hess_mm_indices=hess_mm_indices,
                movable_mm_indices=movable_mm_indices,
                frozen_indices=frozen_layer_indices,
            )
            click.echo(
                f"[annot]   B-factors set in '{final_pdb}' "
                f"(ML={BFACTOR_ML:.0f}, MovableMM={BFACTOR_MOVABLE_MM:.0f}, "
                f"FrozenMM={BFACTOR_FROZEN:.0f})."
            )
        else:
            # Fall back to legacy encoding
            _annotate_b_factors_inplace(
                final_pdb,
                model_pdb=model_pdb,
                freeze_indices_0based=freeze_indices_0based,
            )
            click.echo(f"[annot]   B-factors set in '{final_pdb}' (ML=100, frozen=50, both=150).")
    except Exception as e:
        click.echo(f"[convert] WARNING: Failed to convert final geometry to PDB: {e}", err=True)

    # optimization_all_trj.xyz / optimization_trj.xyz → PDB (if dump)
    if dump:
        try:
            wrote_any = False
            all_trj_path = get_trj_fn("optimization_all_trj.xyz")
            if all_trj_path.exists():
                all_opt_pdb = out_dir / "optimization_all.pdb"
                convert_xyz_to_pdb(all_trj_path, ref_pdb, all_opt_pdb)
                click.echo(f"[convert] Wrote '{all_opt_pdb}'.")
                wrote_any = True

                if use_layer_bfactors:
                    update_pdb_bfactors_from_layers(
                        all_opt_pdb,
                        ml_indices=ml_indices or [],
                        hess_mm_indices=hess_mm_indices,
                        movable_mm_indices=movable_mm_indices,
                        frozen_indices=frozen_layer_indices,
                    )
                    click.echo(
                        f"[annot]   B-factors set in '{all_opt_pdb}' "
                        f"(ML={BFACTOR_ML:.0f}, MovableMM={BFACTOR_MOVABLE_MM:.0f}, "
                        f"FrozenMM={BFACTOR_FROZEN:.0f})."
                    )
                else:
                    _annotate_b_factors_inplace(
                        all_opt_pdb,
                        model_pdb=model_pdb,
                        freeze_indices_0based=freeze_indices_0based,
                    )
                    click.echo(f"[annot]   B-factors set in '{all_opt_pdb}' (ML=100, frozen=50, both=150).")

            trj_path = get_trj_fn("optimization_trj.xyz")
            if trj_path.exists():
                opt_pdb = out_dir / "optimization.pdb"
                convert_xyz_to_pdb(trj_path, ref_pdb, opt_pdb)
                click.echo(f"[convert] Wrote '{opt_pdb}'.")
                wrote_any = True

                if use_layer_bfactors:
                    update_pdb_bfactors_from_layers(
                        opt_pdb,
                        ml_indices=ml_indices or [],
                        hess_mm_indices=hess_mm_indices,
                        movable_mm_indices=movable_mm_indices,
                        frozen_indices=frozen_layer_indices,
                    )
                    click.echo(
                        f"[annot]   B-factors set in '{opt_pdb}' "
                        f"(ML={BFACTOR_ML:.0f}, MovableMM={BFACTOR_MOVABLE_MM:.0f}, "
                        f"FrozenMM={BFACTOR_FROZEN:.0f})."
                    )
                else:
                    _annotate_b_factors_inplace(
                        opt_pdb,
                        model_pdb=model_pdb,
                        freeze_indices_0based=freeze_indices_0based,
                    )
                    click.echo(f"[annot]   B-factors set in '{opt_pdb}' (ML=100, frozen=50, both=150).")

            if not wrote_any:
                click.echo(
                    "[convert] WARNING: neither 'optimization_all_trj.xyz' nor 'optimization_trj.xyz' was found; "
                    "skipping trajectory PDB conversion.",
                    err=True,
                )
        except Exception as e:
            click.echo(f"[convert] WARNING: Failed to convert optimization trajectory to PDB: {e}", err=True)


# -----------------------------------------------
# CLI
# -----------------------------------------------

@click.command(
    help="ML/MM geometry optimization with LBFGS (light) or RFO (heavy).",
    context_settings={"help_option_names": ["-h", "--help"]},
)
@click.option(
    "-i", "--input",
    "input_path",
    type=click.Path(path_type=Path, exists=True, dir_okay=False),
    required=True,
    help="Input structure file (PDB, XYZ). XYZ provides higher coordinate precision. "
         "If XYZ, use --ref-pdb to specify PDB topology for atom ordering and output conversion.",
)
@click.option(
    "--ref-pdb",
    "ref_pdb",
    type=click.Path(path_type=Path, exists=True, dir_okay=False),
    default=None,
    show_default=False,
    help="Reference PDB topology when input is XYZ. XYZ coordinates are used (higher precision) "
         "while PDB provides atom ordering and residue information for output conversion.",
)
@click.option(
    "--real-parm7",
    "real_parm7",
    type=click.Path(path_type=Path, exists=True, dir_okay=False),
    required=True,
    help="Amber parm7 topology covering the whole enzyme complex.",
)
@click.option(
    "--model-pdb",
    "model_pdb",
    type=click.Path(path_type=Path, exists=True, dir_okay=False),
    required=False,
    help="PDB defining atoms that belong to the ML (high-level) region. "
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
    "--model-indices-one-based/--model-indices-zero-based",
    "model_indices_one_based",
    default=True,
    show_default=True,
    help="Interpret --model-indices as 1-based (default) or 0-based.",
)
@click.option(
    "--detect-layer/--no-detect-layer",
    "detect_layer",
    default=True,
    show_default=True,
    help="Detect ML/MM layers from input PDB B-factors (B=0/10/20). "
         "If disabled, you must provide --model-pdb or --model-indices.",
)
@click.option("-q", "--charge", type=int, default=None, show_default=False, help="ML region charge.")
@click.option(
    "-m",
    "--multiplicity",
    "spin",
    type=int,
    default=None,
    show_default=False,
    help="Spin multiplicity (2S+1) for the ML region. Defaults to 1 when omitted.",
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
    "--radius-partial-hessian",
    "--hess-cutoff",
    "radius_partial_hessian",
    type=float,
    default=None,
    show_default=False,
    help="Distance cutoff (Å) from ML region for MM atoms to include in Hessian calculation. "
         "Applied to movable MM atoms and can be combined with --detect-layer. "
         "`--hess-cutoff` is a compatibility alias.",
)
@click.option(
    "--radius-freeze",
    "--movable-cutoff",
    "radius_freeze",
    type=float,
    default=None,
    show_default=False,
    help="Distance cutoff (Å) from ML region for movable MM atoms. "
     "MM atoms beyond this are frozen. "
         "Providing --radius-freeze disables --detect-layer and uses distance-based layer assignment. "
         "`--movable-cutoff` is a compatibility alias.",
)
@click.option(
    "--dist-freeze",
    "dist_freeze_raw",
    type=str,
    multiple=True,
    default=(),
    show_default=False,
    help="Python-like list(s) of (i,j,target_A) to restrain distances (target optional).",
)
@click.option(
    "--one-based/--zero-based",
    "one_based",
    default=True,
    show_default=True,
    help="Interpret --dist-freeze indices as 1-based (default) or 0-based.",
)
@click.option(
    "--bias-k",
    type=float,
    default=10.0,
    show_default=True,
    help="Harmonic restraint strength k [eV/Å^2] for --dist-freeze.",
)
@click.option("--max-cycles", type=int, default=10000, show_default=True, help="Maximum number of optimization cycles.")
@click.option(
    "--dump/--no-dump",
    default=False,
    show_default=True,
    help="Write optimization trajectories ('optimization_trj.xyz' and 'optimization_all_trj.xyz').",
)
@click.option("--out-dir", type=str, default="./result_opt/", show_default=True, help="Output directory.")
@click.option(
    "--thresh",
    type=str,
    default=None,
    help="Convergence preset (gau_loose|gau|gau_tight|gau_vtight|baker|never).",
)
@click.option(
    "--opt-mode",
    type=click.Choice(["light", "heavy", "lbfgs", "rfo"], case_sensitive=False),
    default="light",
    show_default=True,
    help="Optimizer mode: light (LBFGS) or heavy (RFO with Hessian).",
)
@click.option(
    "--layer-opt/--no-layer-opt",
    "layer_opt",
    default=False,
    show_default=True,
    help="Enable LayerOpt-style microiteration (heavy mode only). "
         "Alternates between outer (LBFGS) and inner (RFO) optimization.",
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
    help="Validate options and print the execution plan without running optimization.",
)
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
    spin: Optional[int],
    freeze_atoms_text: Optional[str],
    radius_partial_hessian: Optional[float],
    radius_freeze: Optional[float],
    dist_freeze_raw: Sequence[str],
    one_based: bool,
    bias_k: float,
    max_cycles: int,
    dump: bool,
    out_dir: str,
    thresh: Optional[str],
    opt_mode: str,
    layer_opt: bool,
    config_yaml: Optional[Path],
    show_config: bool,
    dry_run: bool,
) -> None:
    time_start = time.perf_counter()
    prepared_input = None

    def _is_param_explicit(name: str) -> bool:
        try:
            source = ctx.get_parameter_source(name)
            return source not in (None, ParameterSource.DEFAULT)
        except Exception:
            return False

    config_yaml, override_yaml, used_legacy_yaml = _resolve_yaml_sources(
        config_yaml=config_yaml,
        override_yaml=None,
        args_yaml_legacy=None,
    )
    merged_yaml_cfg = _load_merged_yaml_cfg(
        config_yaml=config_yaml,
        override_yaml=None,
    )

    # Handle input: PDB directly, or XYZ with --ref-pdb for topology
    suffix = input_path.suffix.lower()
    if suffix == ".pdb":
        # PDB input: use directly
        prepared_input = prepare_input_structure(input_path)
    elif suffix == ".xyz":
        # XYZ input: require --ref-pdb for topology
        if ref_pdb is None:
            click.echo("ERROR: XYZ/TRJ input requires --ref-pdb to specify PDB topology.", err=True)
            sys.exit(1)
        prepared_input = prepare_input_structure(input_path)
        apply_ref_pdb_override(prepared_input, ref_pdb)
        click.echo(f"[input] Using XYZ coordinates from {input_path.name}, PDB topology from {ref_pdb.name}")
    else:
        click.echo(f"ERROR: Unsupported input format: {suffix}. Use .pdb or .xyz (with --ref-pdb).", err=True)
        sys.exit(1)

    geom_input_path = prepared_input.geom_path
    charge, spin = resolve_charge_spin_or_raise(prepared_input, charge, spin)

    try:
        freeze_atoms_cli = _parse_freeze_atoms(freeze_atoms_text)
    except click.BadParameter as e:
        click.echo(f"ERROR: {e}", err=True)
        prepared_input.cleanup()
        sys.exit(1)

    model_indices: Optional[List[int]] = None
    if model_indices_str:
        try:
            model_indices = parse_indices_string(model_indices_str, one_based=model_indices_one_based)
        except click.BadParameter as e:
            click.echo(f"ERROR: {e}", err=True)
            prepared_input.cleanup()
            sys.exit(1)

    try:
        dist_freeze = _parse_dist_freeze(dist_freeze_raw, one_based=bool(one_based))
    except click.BadParameter as e:
        click.echo(f"ERROR: {e}", err=True)
        prepared_input.cleanup()
        sys.exit(1)

    # Resolve optimizer mode
    mode_resolved = normalize_choice(
        opt_mode,
        param="--opt-mode",
        alias_groups=OPT_MODE_ALIASES,
        allowed_hint="light, heavy",
    )
    use_rfo = (mode_resolved == "rfo")

    try:
        config_layer_cfg = load_yaml_dict(config_yaml)
        override_layer_cfg = load_yaml_dict(override_yaml)
        geom_cfg = dict(GEOM_KW)
        calc_cfg = dict(CALC_KW)
        opt_cfg = dict(OPT_BASE_KW)
        lbfgs_cfg = dict(LBFGS_KW)
        rfo_cfg = dict(RFO_KW)

        apply_yaml_overrides(
            config_layer_cfg,
            [
                (geom_cfg, (("geom",),)),
                (calc_cfg, (("calc",), ("mlmm",))),
                (opt_cfg, (("opt",),)),
                (lbfgs_cfg, (("lbfgs",), ("opt", "lbfgs"))),
                (rfo_cfg, (("rfo",), ("opt", "rfo"))),
            ],
        )

        if _is_param_explicit("max_cycles"):
            opt_cfg["max_cycles"] = int(max_cycles)
        if _is_param_explicit("dump"):
            opt_cfg["dump"] = bool(dump)
        if _is_param_explicit("out_dir"):
            opt_cfg["out_dir"] = out_dir
        if _is_param_explicit("thresh") and thresh is not None:
            opt_cfg["thresh"] = str(thresh)

        if _is_param_explicit("detect_layer"):
            calc_cfg["use_bfactor_layers"] = bool(detect_layer)
        if _is_param_explicit("radius_partial_hessian") and radius_partial_hessian is not None:
            calc_cfg["hess_cutoff"] = float(radius_partial_hessian)
        if _is_param_explicit("radius_freeze") and radius_freeze is not None:
            calc_cfg["movable_cutoff"] = float(radius_freeze)
            calc_cfg["use_bfactor_layers"] = False

        model_charge_value = calc_cfg.get("model_charge", charge)
        if model_charge_value is None:
            model_charge_value = charge
        calc_cfg["model_charge"] = int(model_charge_value)
        if _is_param_explicit("charge"):
            calc_cfg["model_charge"] = int(charge)
        model_mult_value = calc_cfg.get("model_mult", spin)
        if model_mult_value is None:
            model_mult_value = spin
        calc_cfg["model_mult"] = int(model_mult_value)
        if _is_param_explicit("spin"):
            calc_cfg["model_mult"] = int(spin)
        if model_pdb is not None:
            calc_cfg["model_pdb"] = str(model_pdb)
        calc_cfg["input_pdb"] = str(prepared_input.source_path)
        calc_cfg["real_parm7"] = str(real_parm7)

        apply_yaml_overrides(
            override_layer_cfg,
            [
                (geom_cfg, (("geom",),)),
                (calc_cfg, (("calc",), ("mlmm",))),
                (opt_cfg, (("opt",),)),
                (lbfgs_cfg, (("lbfgs",), ("opt", "lbfgs"))),
                (rfo_cfg, (("rfo",), ("opt", "rfo"))),
            ],
        )

        try:
            geom_freeze = _normalize_geom_freeze(geom_cfg.get("freeze_atoms"))
        except click.BadParameter as e:
            click.echo(f"ERROR: {e}", err=True)
            prepared_input.cleanup()
            sys.exit(1)
        geom_cfg["freeze_atoms"] = geom_freeze
        if freeze_atoms_cli:
            merge_freeze_atom_indices(geom_cfg, freeze_atoms_cli)
        freeze_atoms_final = list(geom_cfg.get("freeze_atoms") or [])
        calc_cfg["freeze_atoms"] = freeze_atoms_final

        out_dir_path = Path(opt_cfg["out_dir"]).resolve()

        # radius_freeze implies full distance-based layer assignment.
        # radius_partial_hessian alone can be combined with --detect-layer.
        detect_layer_enabled = bool(calc_cfg.get("use_bfactor_layers", True))
        model_pdb_cfg = calc_cfg.get("model_pdb")
        if radius_freeze is not None:
            if detect_layer_enabled:
                click.echo("[layer] --radius-freeze provided; disabling --detect-layer.", err=True)
            detect_layer_enabled = False
            calc_cfg["use_bfactor_layers"] = False

        layer_source_pdb = prepared_input.source_path
        if detect_layer_enabled and layer_source_pdb.suffix.lower() != ".pdb":
            click.echo("ERROR: --detect-layer requires a PDB input (or --ref-pdb).", err=True)
            prepared_input.cleanup()
            sys.exit(1)

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
                    prepared_input.cleanup()
                    sys.exit(1)
            if (
                not detect_layer_enabled
                and model_pdb_cfg is None
                and model_indices
                and layer_source_pdb.suffix.lower() != ".pdb"
            ):
                click.echo("ERROR: --model-indices requires a PDB input (or --ref-pdb).", err=True)
                prepared_input.cleanup()
                sys.exit(1)
            click.echo(
                pretty_block(
                    "dry_run_plan",
                    {
                        "input_geometry": str(geom_input_path),
                        "output_dir": str(out_dir_path),
                        "optimizer_mode": "rfo" if use_rfo else "lbfgs",
                        "detect_layer": bool(detect_layer_enabled),
                        "model_region_source": model_region_source,
                        "model_indices_count": 0 if not model_indices else len(model_indices),
                        "will_run_optimization": True,
                        "will_convert_outputs": True,
                    },
                )
            )
            click.echo("[dry-run] Validation complete. Optimization execution was skipped.")
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
                    prepared_input.cleanup()
                    sys.exit(1)
                click.echo(f"[layer] WARNING: {e} Falling back to explicit ML region.", err=True)
                detect_layer_enabled = False

        if not detect_layer_enabled:
            if model_pdb_cfg is None and not model_indices:
                click.echo("ERROR: Provide --model-pdb or --model-indices when --no-detect-layer.", err=True)
                prepared_input.cleanup()
                sys.exit(1)
            if model_pdb_cfg is not None:
                model_pdb_path = Path(model_pdb_cfg)
            else:
                if layer_source_pdb.suffix.lower() != ".pdb":
                    click.echo("ERROR: --model-indices requires a PDB input (or --ref-pdb).", err=True)
                    prepared_input.cleanup()
                    sys.exit(1)
                try:
                    model_pdb_path = build_model_pdb_from_indices(layer_source_pdb, out_dir_path, model_indices or [])
                except Exception as e:
                    click.echo(f"ERROR: {e}", err=True)
                    prepared_input.cleanup()
                    sys.exit(1)
            calc_cfg["use_bfactor_layers"] = False

        if model_pdb_path is None:
            click.echo("ERROR: Failed to resolve model PDB for the ML region.", err=True)
            prepared_input.cleanup()
            sys.exit(1)

        calc_cfg["model_pdb"] = str(model_pdb_path)

        # When layer detection is enabled, also freeze frozen-layer atoms at the
        # optimizer geometry level (not only inside the calculator).
        # Otherwise LBFGS may still move those coordinates through coupled
        # inverse-Hessian updates, even if raw forces are zeroed there.
        if layer_info is not None:
            frozen_from_layer = [int(i) for i in layer_info.get("frozen_indices", [])]
            if frozen_from_layer:
                before = set(freeze_atoms_final)
                merged = sorted(before | set(frozen_from_layer))
                added = len(set(merged) - before)
                freeze_atoms_final = merged
                geom_cfg["freeze_atoms"] = freeze_atoms_final
                calc_cfg["freeze_atoms"] = freeze_atoms_final
                click.echo(
                    f"[layer] Applied optimizer freeze constraints: "
                    f"total={len(freeze_atoms_final)} (added_from_layer={added})"
                )

        # Distance-based overrides for Hessian-target and movable MM selection.
        hess_cutoff_final = calc_cfg.get("hess_cutoff")
        movable_cutoff_final = calc_cfg.get("movable_cutoff")
        if hess_cutoff_final is not None or movable_cutoff_final is not None:
            click.echo(
                f"[layer] Applied distance cutoffs: "
                f"hess={hess_cutoff_final} Å, freeze={movable_cutoff_final} Å"
            )

        for key in ("input_pdb", "real_parm7", "model_pdb", "mm_fd_dir"):
            val = calc_cfg.get(key)
            if val:
                calc_cfg[key] = str(Path(val).expanduser().resolve())

        mode_str = "LBFGS (light)"
        if use_rfo:
            mode_str = "RFO (heavy)" if not layer_opt else "RFO + LayerOpt (heavy)"
        click.echo(f"\n[mode] Optimizer: {mode_str}\n")
        click.echo(pretty_block("geom", format_freeze_atoms_for_echo(geom_cfg, key="freeze_atoms")))
        # Show only non-default calc settings for concise logging
        echo_calc = strip_inherited_keys(calc_cfg, CALC_KW, mode="same")
        echo_calc = format_freeze_atoms_for_echo(echo_calc, key="freeze_atoms")
        click.echo(pretty_block("calc", echo_calc))
        # Show only non-default opt settings
        echo_opt = strip_inherited_keys({**opt_cfg, "out_dir": str(out_dir_path)}, OPT_BASE_KW, mode="same")
        click.echo(pretty_block("opt", echo_opt))
        # Show only optimizer-specific settings, not inherited from opt_cfg
        if use_rfo:
            echo_rfo = strip_inherited_keys(rfo_cfg, opt_cfg)
            click.echo(pretty_block("rfo", echo_rfo))
        else:
            echo_lbfgs = strip_inherited_keys(lbfgs_cfg, opt_cfg)
            click.echo(pretty_block("lbfgs", echo_lbfgs))
        if dist_freeze:
            display_pairs = []
            for (i, j, target) in dist_freeze:
                label = (f"{target:.4f}" if target is not None else "<current>")
                display_pairs.append((int(i) + 1, int(j) + 1, label))
            click.echo(
                pretty_block(
                    "dist_freeze (input)",
                    {
                        "k (eV/Å^2)": float(bias_k),
                        "pairs_1based": display_pairs,
                    },
                )
            )

        out_dir_path.mkdir(parents=True, exist_ok=True)
        coord_type = geom_cfg.get("coord_type", "cart")
        coord_kwargs = dict(geom_cfg)
        coord_kwargs.pop("coord_type", None)
        geometry = geom_loader(
            geom_input_path,
            coord_type=coord_type,
            **coord_kwargs,
        )

        base_calc = mlmm(**calc_cfg)
        geometry.set_calculator(base_calc)

        resolved_dist_freeze: List[Tuple[int, int, float]] = []
        if dist_freeze:
            try:
                resolved_dist_freeze = _resolve_dist_freeze_targets(geometry, dist_freeze)
            except click.BadParameter as e:
                click.echo(f"ERROR: {e}", err=True)
                sys.exit(1)
            click.echo(
                pretty_block(
                    "dist_freeze (active)",
                    {
                        "k (eV/Å^2)": float(bias_k),
                        "pairs_1based": [
                            (int(i) + 1, int(j) + 1, float(f"{t:.4f}"))
                            for (i, j, t) in resolved_dist_freeze
                        ],
                    },
                )
            )
            bias_calc = HarmonicBiasCalculator(base_calc, k=float(bias_k))
            bias_calc.set_pairs(resolved_dist_freeze)
            geometry.set_calculator(bias_calc)

        common_kwargs = dict(opt_cfg)
        common_kwargs["out_dir"] = str(out_dir_path)
        used_partial_micro = False

        # Check layer_opt compatibility
        if layer_opt and not use_rfo:
            click.echo("WARNING: --layer-opt requires heavy mode. Ignoring --layer-opt.", err=True)
            layer_opt = False

        if layer_opt:
            # LayerOpt mode: microiteration optimization
            # Get ML region indices from calculator
            ml_indices = list(base_calc.core.selection_indices)
            if not ml_indices:
                click.echo("ERROR: No ML region indices found for LayerOpt.", err=True)
                sys.exit(1)

            # Load LayerOpt defaults
            layeropt_cfg = dict(LAYEROPT_KW)

            # Prepare outer LBFGS kwargs
            outer_lbfgs_kwargs = {**lbfgs_cfg, **common_kwargs}
            outer_lbfgs_kwargs["thresh"] = layeropt_cfg.get("outer_thresh", "gau_loose")

            # Prepare inner RFO kwargs
            inner_rfo_kwargs = {**rfo_cfg, **common_kwargs}
            inner_rfo_kwargs["thresh"] = opt_cfg.get("thresh", "gau")

            click.echo(pretty_block("layeropt", layeropt_cfg))

            optimizer = MicroIterationOptimizer(
                geometry=geometry,
                ml_indices=ml_indices,
                outer_lbfgs_kwargs=outer_lbfgs_kwargs,
                inner_rfo_kwargs=inner_rfo_kwargs,
                max_macro_cycles=int(opt_cfg.get("max_cycles", 10000)),
                outer_max_cycles=layeropt_cfg.get("outer_max_cycles", 1500),
                outer_thresh=layeropt_cfg.get("outer_thresh", "gau_loose"),
                echo_fn=click.echo,
                dump=bool(opt_cfg["dump"]),
            )

            click.echo("\n=== LayerOpt Optimization started ===\n")
            optimizer.run()
            click.echo("\n=== LayerOpt Optimization finished ===\n")

            # Write final geometry
            final_xyz_path = out_dir_path / "final_geometry.xyz"
            final_xyz_path.write_text(geometry.as_xyz(), encoding="utf-8")
        elif use_rfo:
            # Use partial-Hessian microiterations when non-Hessian movable atoms exist.
            calc_core = base_calc.core if hasattr(base_calc, "core") else base_calc
            movable_mm_indices = list(getattr(calc_core, "movable_mm_indices", []) or [])
            hess_active_atoms = list(getattr(calc_core, "hess_active_atoms", []) or [])
            use_partial_micro = bool(movable_mm_indices) and bool(hess_active_atoms)

            if use_partial_micro:
                used_partial_micro = True
                layeropt_cfg = dict(LAYEROPT_KW)
                outer_lbfgs_kwargs = {**lbfgs_cfg, **common_kwargs}
                inner_rfo_kwargs = {**rfo_cfg, **common_kwargs}
                outer_max_cycles = int(layeropt_cfg.get("outer_max_cycles", 1500))
                if "max_cycles" in opt_cfg:
                    outer_max_cycles = min(outer_max_cycles, int(opt_cfg["max_cycles"]))

                # Inner RFO operates on the active DOF block (Hessian-target atoms).
                # Ensure the calculator returns partial Hessians to keep Hessian/step
                # dimensions consistent inside RFOptimizer updates.
                if hasattr(calc_core, "return_partial_hessian"):
                    calc_core.return_partial_hessian = True

                optimizer = PartialHessianMicroIterationOptimizer(
                    geometry=geometry,
                    hess_indices=hess_active_atoms,
                    outer_indices=movable_mm_indices,
                    outer_lbfgs_kwargs=outer_lbfgs_kwargs,
                    inner_rfo_kwargs=inner_rfo_kwargs,
                    max_macro_cycles=int(opt_cfg.get("max_cycles", 10000)),
                    outer_max_cycles=outer_max_cycles,
                    outer_thresh=layeropt_cfg.get("outer_thresh", "gau_loose"),
                    echo_fn=click.echo,
                    dump=bool(opt_cfg["dump"]),
                )

                click.echo("\n=== Partial-Hessian MicroIteration (RFO) started ===\n")
                optimizer.run()
                click.echo("\n=== Partial-Hessian MicroIteration (RFO) finished ===\n")
            else:
                rfo_args = {**rfo_cfg, **common_kwargs}
                optimizer = RFOptimizer(geometry, **rfo_args)

                click.echo("\n=== Optimization started ===\n")
                optimizer.run()
                click.echo("\n=== Optimization finished ===\n")
        else:
            lbfgs_args = {**lbfgs_cfg, **common_kwargs}
            optimizer = LBFGS(geometry, **lbfgs_args)

            click.echo("\n=== Optimization started ===\n")
            optimizer.run()
            click.echo("\n=== Optimization finished ===\n")

        # Get final geometry path
        if layer_opt or used_partial_micro:
            final_xyz_path = out_dir_path / "final_geometry.xyz"
            final_xyz_path.write_text(geometry.as_xyz(), encoding="utf-8")
        else:
            final_xyz_path = optimizer.final_fn if isinstance(optimizer.final_fn, Path) else Path(optimizer.final_fn)

        if bool(opt_cfg["dump"]):
            optim_all_path = out_dir_path / "optimization_all_trj.xyz"
            if not optim_all_path.exists():
                trj_path = optimizer.get_path_for_fn("optimization_trj.xyz")
                _append_xyz_trajectory(optim_all_path, trj_path, reset=True)

        # Extract layer indices from calculator for layer-based B-factor encoding
        calc_core = base_calc.core if hasattr(base_calc, 'core') else base_calc
        ml_indices = getattr(calc_core, 'ml_indices', None)
        hess_mm_indices = getattr(calc_core, 'hess_mm_indices', None)
        movable_mm_indices = getattr(calc_core, 'movable_mm_indices', None)
        frozen_layer_indices = getattr(calc_core, 'frozen_layer_indices', None)

        _maybe_convert_outputs_to_pdb(
            input_path=prepared_input.source_path,  # Use PDB topology for conversion
            out_dir=out_dir_path,
            dump=bool(opt_cfg["dump"]),
            get_trj_fn=optimizer.get_path_for_fn,
            final_xyz_path=final_xyz_path,
            model_pdb=Path(calc_cfg["model_pdb"]),
            freeze_indices_0based=freeze_atoms_final,
            ml_indices=ml_indices,
            hess_mm_indices=hess_mm_indices,
            movable_mm_indices=movable_mm_indices,
            frozen_layer_indices=frozen_layer_indices,
        )

        _write_output_summary_md(out_dir_path)
        click.echo(format_elapsed("[time] Elapsed Time for Opt", time_start))

    except ZeroStepLength:
        click.echo("ERROR: Step length fell below the minimum allowed (ZeroStepLength).", err=True)
        sys.exit(2)
    except OptimizationError as e:
        click.echo(f"ERROR: Optimization failed - {e}", err=True)
        sys.exit(3)
    except KeyboardInterrupt:
        click.echo("\nInterrupted by user.", err=True)
        sys.exit(130)
    except Exception as e:
        tb = "".join(traceback.format_exception(type(e), e, e.__traceback__))
        click.echo("Unhandled exception during optimization:\n" + textwrap.indent(tb, "  "), err=True)
        sys.exit(1)
    finally:
        if prepared_input is not None:
            prepared_input.cleanup()


# Allow `python -m mlmm_toolkit.opt` direct execution
if __name__ == "__main__":
    cli()
