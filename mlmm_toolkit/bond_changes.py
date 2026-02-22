# bond_changes.py

"""
bond_changes — Bond-change detection and reporting utilities
====================================================================

Usage (API)
-----
    from <package>.bond_changes import compare_structures, summarize_changes

Examples::
    >>> result = compare_structures(geom_reactant, geom_product, device="cpu")
    >>> print(summarize_changes(geom_product, result))
    Bond formed (1)
      C1-O2 1.50 Å --> 1.36 Å

Description
-----
This module compares two molecular geometries with identical atom types and ordering and reports covalent bonds that are formed
or broken between the structures. Bond perception uses element-specific covalent radii and configurable tolerances; distances are computed with PyTorch on CPU or CUDA.

Algorithm (core logic):
- Inputs: `geom1`, `geom2` with attributes `atoms: Iterable[str]` and `coords3d: (N, 3) float` (in Bohr in pysisyphus). Atoms must match exactly (`assert geom1.atoms == geom2.atoms`).
- Per-element radii: from `pysisyphus.elem_data.COVALENT_RADII`.
- Threshold per pair (i, j): `T_cov = bond_factor * (r_i + r_j)`; conservative margin `eps = margin_fraction * T_cov`.
- Bond adjacency in a geometry: `A = (D <= T_cov - eps)` evaluated only for `i < j` (upper triangle).
- Change gating: only pairs with `|D2 - D1| >= delta_fraction * T_cov` are considered.
- Classification:
    - Formed: `(~A1) & A2 & need_change`
    - Broken:  `A1 & (~A2) & need_change`
- Distances: pairwise matrices `D1`, `D2` via `torch.cdist`.

Public API:
- `compare_structures(geom1, geom2, device='cuda', bond_factor=1.20, margin_fraction=0.05, delta_fraction=0.05) -> BondChangeResult`
  Detects formed and broken covalent bonds. Returns sets of zero-based index pairs and (by default) the full distance matrices (`numpy.ndarray`) in the same units as the inputs (Bohr in pysisyphus).
- `summarize_changes(geom, result, one_based: bool = True) -> str`
  Builds a human-readable report:
  - Sections: “Bond formed” and “Bond broken” (with counts).
  - Lines formatted as `ElemI-ElemJ` with atom indices (1-based by default, e.g., `C1-O2`).
  - If `result.distances_1/2` are present, prints bond lengths as `D1 Å --> D2 Å`, converting from Bohr using `pysisyphus.constants.BOHR2ANG`.
- Helper utilities (internal):
  - `_resolve_device(device: str) -> torch.device`: chooses the requested device; falls back to CPU with a warning if unavailable.
  - `_element_arrays(atoms) -> (elems, cov_radii)`: normalizes element symbols and looks up covalent radii.
  - `_upper_pairs_from_mask(mask) -> Set[Tuple[int, int]]`: returns index pairs where `mask` is True (assumes an upper-triangular mask).
  - `_bond_str(i, j, elems, one_based=True) -> str`: formats `ElemI-ElemJ` labels.
- Data container:
  - `BondChangeResult`:
    - `formed_covalent: Set[Tuple[int, int]]` — zero-based pairs for bonds formed.
    - `broken_covalent: Set[Tuple[int, int]]` — zero-based pairs for bonds broken.
    - `distances_1: Optional[np.ndarray]`, `distances_2: Optional[np.ndarray]` — square distance matrices (shape N×N).

Outputs (& Directory Layout)
-----
- No files or directories are created.
- `compare_structures` returns a `BondChangeResult` as described above.
- `summarize_changes` returns a multi-line string; typical headings:
  - `Bond formed (k):` followed by `ElemI-ElemJ : D1 Å --> D2 Å` (if distances available).
  - `Bond broken (m):` followed by lines in the same format.
  - If a set is empty, the section reads `None`.

Notes:
-----
- Units: In pysisyphus, `coords3d` are Bohr; the summary converts to Å with `BOHR2ANG`. If your inputs use different units, adjust accordingly.
- Atom identity & ordering must be identical between structures; otherwise the comparison is invalid.
- Only unique pairs with `i < j` are considered (upper triangle); indices in results are zero-based.
- The three tolerances control sensitivity:
  - `bond_factor` (default 1.20): global scaling of the covalent radii sum.
  - `margin_fraction` (default 0.05): conservative shrinkage of the bond cutoff to avoid borderline matches.
  - `delta_fraction` (default 0.05): minimum relative distance change required to count a bond event.
- Device selection: pass `'cpu'`, `'cuda'`, or `'cuda:0'` etc. If the requested device is not available or invalid, the code falls back to CPU and issues a `RuntimeWarning`. Computations use `float64` for stability and run under `torch.no_grad()`.
- The method detects binary bond formation/breakage; it does not estimate bond orders, angles, or multi-center bonding, and it ignores periodic boundary conditions.
- Numerical caveats: near-threshold pairs may toggle with small geometry noise; tune tolerances to your system size and sampling noise.
"""

from __future__ import annotations
from dataclasses import dataclass
from typing import Iterable, Tuple, Set, List, Optional, Dict, Any

import warnings
import torch
import numpy as np

from pysisyphus.elem_data import (
    COVALENT_RADII as CR,
)
from pysisyphus.constants import BOHR2ANG  # Convert Bohr distances to Ångström for reporting

Pair = Tuple[int, int]


@dataclass
class BondChangeResult:
    formed_covalent: Set[Pair]
    broken_covalent: Set[Pair]
    distances_1: Optional[np.ndarray] = None
    distances_2: Optional[np.ndarray] = None


def _upper_pairs_from_mask(mask: torch.Tensor) -> Set[Pair]:
    idx = torch.nonzero(mask, as_tuple=False).detach().cpu().numpy()
    return set(map(tuple, idx))


def _element_arrays(atoms: Iterable[str]) -> Tuple[List[str], np.ndarray]:
    elems = [a.capitalize() for a in atoms]
    cov = np.array([CR[a.lower()] for a in elems], dtype=float)
    return elems, cov


def _resolve_device(device: str) -> torch.device:
    dev_str = (device or "cpu").lower()
    if dev_str.startswith("cuda"):
        if torch.cuda.is_available():
            try:
                _ = torch.device(dev_str)
                return torch.device(dev_str)
            except Exception:
                warnings.warn(
                    f"Requested device '{device}' is not available. Falling back to CPU.",
                    RuntimeWarning,
                )
                return torch.device("cpu")
        else:
            warnings.warn(
                "CUDA is not available. Falling back to CPU.",
                RuntimeWarning,
            )
            return torch.device("cpu")
    try:
        return torch.device(dev_str)
    except Exception:
        warnings.warn(
            f"Requested device '{device}' is not recognized. Falling back to CPU.",
            RuntimeWarning,
        )
        return torch.device("cpu")


@torch.no_grad()
def compare_structures(
    geom1,
    geom2,
    device: str = "cuda",
    bond_factor: float = 1.20,
    margin_fraction: float = 0.05,
    delta_fraction: float = 0.05,
) -> BondChangeResult:

    assert geom1.atoms == geom2.atoms, "Atom types and ordering must be identical."
    N = len(geom1.atoms)

    elems, cov_np = _element_arrays(geom1.atoms)
    dev = _resolve_device(device)

    dtype = torch.float64
    R1 = torch.as_tensor(geom1.coords3d, dtype=dtype, device=dev)
    R2 = torch.as_tensor(geom2.coords3d, dtype=dtype, device=dev)
    cov = torch.as_tensor(cov_np, dtype=dtype, device=dev)

    T_cov = bond_factor * (cov[:, None] + cov[None, :])
    eps_cov = margin_fraction * T_cov

    D1 = torch.cdist(R1, R1)
    D2 = torch.cdist(R2, R2)

    up = torch.triu(torch.ones((N, N), dtype=torch.bool, device=dev), diagonal=1)

    A1 = (D1 <= (T_cov - eps_cov)) & up
    A2 = (D2 <= (T_cov - eps_cov)) & up

    dD = D2 - D1
    need_change = (dD.abs() >= (delta_fraction * T_cov)) & up

    formed_cov_mask = (~A1) & A2 & need_change
    broken_cov_mask = A1 & (~A2) & need_change

    formed_covalent = _upper_pairs_from_mask(formed_cov_mask)
    broken_covalent = _upper_pairs_from_mask(broken_cov_mask)

    return BondChangeResult(
        formed_covalent=formed_covalent,
        broken_covalent=broken_covalent,
        distances_1=D1.detach().cpu().numpy(),
        distances_2=D2.detach().cpu().numpy(),
    )


def _bond_str(i: int, j: int, elems: List[str], one_based: bool = True) -> str:
    ii = i + 1 if one_based else i
    jj = j + 1 if one_based else j
    return f"{elems[i]}{ii}-{elems[j]}{jj}"


def summarize_changes(geom, result: BondChangeResult, one_based: bool = True) -> str:
    """
    List bond formations and dissociations and report bond-length changes in Å.
    """
    elems = [a.capitalize() for a in geom.atoms]
    lines: List[str] = []

    # Use distance matrices (Bohr) converted to Å when available
    D1 = result.distances_1
    D2 = result.distances_2
    have_lengths = (
        isinstance(D1, np.ndarray)
        and isinstance(D2, np.ndarray)
        and D1.shape == D2.shape
    )

    def _len_str(i: int, j: int) -> str:
        if not have_lengths:
            return ""
        # ``coords3d`` is given in Bohr; convert to Å
        d1 = float(D1[i, j]) * BOHR2ANG
        d2 = float(D2[i, j]) * BOHR2ANG
        return f" : {d1:.3f} Å --> {d2:.3f} Å"

    def pairs_to_lines(title: str, pairs: Set[Pair]):
        if not pairs:
            lines.append(f"{title}: None")
            return
        lines.append(f"{title} ({len(pairs)}):")
        for i, j in sorted(pairs):
            lines.append(f"  - {_bond_str(i, j, elems, one_based)}{_len_str(i, j)}")

    pairs_to_lines("Bond formed", result.formed_covalent)
    pairs_to_lines("Bond broken", result.broken_covalent)

    return "\n".join(lines)


def has_bond_change(
    geom_start,
    geom_end,
    bond_cfg: Dict[str, Any],
    *,
    one_based: bool = True,
) -> Tuple[bool, str]:
    """
    Convenience wrapper to check if bond changes occur between two geometries.

    Parameters
    ----------
    geom_start : Geometry
        Starting geometry.
    geom_end : Geometry
        Ending geometry.
    bond_cfg : Dict[str, Any]
        Configuration dict with keys: 'device', 'bond_factor', 'margin_fraction', 'delta_fraction'.
    one_based : bool
        Use 1-based atom indices in summary (default True).

    Returns
    -------
    Tuple[bool, str]
        (has_changes, summary) where has_changes is True if any bonds formed or broken.
    """
    device = bond_cfg.get("device", "cuda")
    bond_factor = bond_cfg.get("bond_factor", 1.20)
    margin_fraction = bond_cfg.get("margin_fraction", 0.05)
    delta_fraction = bond_cfg.get("delta_fraction", 0.05)

    result = compare_structures(
        geom_start,
        geom_end,
        device=device,
        bond_factor=bond_factor,
        margin_fraction=margin_fraction,
        delta_fraction=delta_fraction,
    )

    has_changes = bool(result.formed_covalent) or bool(result.broken_covalent)
    summary = summarize_changes(geom_start, result, one_based=one_based)

    return has_changes, summary
