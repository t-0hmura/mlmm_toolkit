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
- Distances: row-chunked `torch.cdist` (rows processed in blocks) to bound peak memory to O(N·chunk); the full N×N matrices are not materialized.

Public API:
- `compare_structures(geom1, geom2, device='cuda', bond_factor=1.20, margin_fraction=0.05, delta_fraction=0.05) -> BondChangeResult`
  Detects formed and broken covalent bonds. Returns sets of zero-based index pairs plus `changed_lengths` (per-pair lengths in Bohr for the changed pairs only). The full distance matrices (`distances_1/2`) are no longer materialized or returned (memory-bounded path).
- `summarize_changes(geom, result, one_based: bool = True) -> str`
  Builds a human-readable report:
  - Sections: “Bond formed” and “Bond broken” (with counts).
  - Lines formatted as `ElemI-ElemJ` with atom indices (1-based by default, e.g., `C1-O2`).
  - Bond lengths come from `result.changed_lengths` (sparse per-pair map, preferred) or, when present, the full `result.distances_1/2` matrices; printed as `D1 Å --> D2 Å`, converting from Bohr using `pysisyphus.constants.BOHR2ANG`.
- Helper utilities (internal):
  - `_resolve_device(device: str) -> torch.device`: chooses the requested device; falls back to CPU with a warning if unavailable.
  - `_element_arrays(atoms) -> (elems, cov_radii)`: normalizes element symbols and looks up covalent radii.
  - `_bond_str(i, j, elems, one_based=True) -> str`: formats `ElemI-ElemJ` labels.
- Data container:
  - `BondChangeResult`:
    - `formed_covalent: Set[Tuple[int, int]]` — zero-based pairs for bonds formed.
    - `broken_covalent: Set[Tuple[int, int]]` — zero-based pairs for bonds broken.
    - `distances_1: Optional[np.ndarray]`, `distances_2: Optional[np.ndarray]` — full N×N distance matrices; `None` by default (memory-bounded path).
    - `changed_lengths: Optional[Dict[Tuple[int, int], Tuple[float, float]]]` — per-pair (d1, d2) lengths in Bohr for the changed pairs only.

Outputs (& Directory Layout)
-----
- No files or directories are created.
- `compare_structures` returns a `BondChangeResult` as described above.
- `summarize_changes` returns a multi-line string; typical headings:
  - `Bond formed (k):` followed by `ElemI-ElemJ : D1 Å --> D2 Å` (when lengths available).
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
from typing import Iterable, Tuple, Set, List, Optional, Dict

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
    # Sparse per-pair lengths (Bohr) for the changed pairs only; populated when
    # the full distance matrices are not materialized (memory-bounded path).
    changed_lengths: Optional[Dict[Pair, Tuple[float, float]]] = None


def _element_arrays(atoms: Iterable[str]) -> Tuple[List[str], np.ndarray]:
    elems = [a.capitalize() for a in atoms]
    cov = np.array([CR[a.lower()] for a in elems], dtype=float)
    return elems, cov


def _resolve_device(device: str) -> torch.device:
    dev_str = (device or "cpu").lower()
    if dev_str == "auto":
        dev_str = "cuda" if torch.cuda.is_available() else "cpu"
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

    if geom1.atoms != geom2.atoms:
        raise ValueError("Atom types and ordering must be identical.")
    N = len(geom1.atoms)

    elems, cov_np = _element_arrays(geom1.atoms)
    dev = _resolve_device(device)

    dtype = torch.float64
    R1 = torch.as_tensor(geom1.coords3d, dtype=dtype, device=dev)
    R2 = torch.as_tensor(geom2.coords3d, dtype=dtype, device=dev)
    cov = torch.as_tensor(cov_np, dtype=dtype, device=dev)

    # Row-chunked covalent bond-change detection. A dense formulation needs ~10
    # simultaneous N×N float64 tensors (D1, D2, masks); that is O(N^2) memory and
    # OOMs on large systems (full solvated clusters of ~20k+ atoms exceed 16-24 GB
    # GPUs via torch.cdist). Processing rows in blocks bounds peak memory to
    # O(N*chunk). Only formed/broken covalent pairs are consumed downstream, so the
    # full N×N distance matrices are no longer materialized or returned.
    formed_idx: List[torch.Tensor] = []
    broken_idx: List[torch.Tensor] = []
    chunk = max(1, min(N, 1024))   # ~chunk*N float64 peak; small enough to coexist with a resident ML model on a 16 GB GPU
    cols_g = torch.arange(N, device=dev)[None, :]
    for i0 in range(0, N, chunk):
        i1 = min(i0 + chunk, N)
        T = bond_factor * (cov[i0:i1, None] + cov[None, :])      # (b, N)
        eps = margin_fraction * T
        d1 = torch.cdist(R1[i0:i1], R1)                          # (b, N)
        d2 = torch.cdist(R2[i0:i1], R2)
        up = cols_g > torch.arange(i0, i1, device=dev)[:, None]  # strict upper triangle
        a1 = (d1 <= (T - eps)) & up
        a2 = (d2 <= (T - eps)) & up
        need = ((d2 - d1).abs() >= (delta_fraction * T)) & up
        fi = torch.nonzero((~a1) & a2 & need, as_tuple=False)
        bi = torch.nonzero(a1 & (~a2) & need, as_tuple=False)
        if fi.numel():
            fi[:, 0] += i0
            formed_idx.append(fi.detach().cpu())
        if bi.numel():
            bi[:, 0] += i0
            broken_idx.append(bi.detach().cpu())

    def _to_pairs(chunks: List[torch.Tensor]) -> Set[Pair]:
        if not chunks:
            return set()
        return set(map(tuple, torch.cat(chunks).numpy()))

    formed_covalent = _to_pairs(formed_idx)
    broken_covalent = _to_pairs(broken_idx)

    # Per-pair bond lengths (Bohr) for the changed pairs only — keeps the
    # "d1 Å --> d2 Å" reporting in summarize_changes without an O(N^2) matrix.
    changed_lengths: Dict[Pair, Tuple[float, float]] = {}
    for (i, j) in formed_covalent | broken_covalent:
        d1 = float(torch.linalg.norm(R1[i] - R1[j]))
        d2 = float(torch.linalg.norm(R2[i] - R2[j]))
        changed_lengths[(i, j)] = (d1, d2)

    return BondChangeResult(
        formed_covalent=formed_covalent,
        broken_covalent=broken_covalent,
        distances_1=None,
        distances_2=None,
        changed_lengths=changed_lengths,
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

    # Bond lengths (Bohr) → Å. Prefer the sparse per-pair map (memory-bounded
    # path); fall back to the full distance matrices when present.
    D1 = result.distances_1
    D2 = result.distances_2
    lengths = result.changed_lengths
    have_matrices = (
        isinstance(D1, np.ndarray)
        and isinstance(D2, np.ndarray)
        and D1.shape == D2.shape
    )

    def _len_str(i: int, j: int) -> str:
        if lengths is not None and (i, j) in lengths:
            d1, d2 = lengths[(i, j)]
        elif have_matrices:
            d1 = float(D1[i, j])
            d2 = float(D2[i, j])
        else:
            return ""
        # ``coords3d`` is given in Bohr; convert to Å
        return f" : {d1 * BOHR2ANG:.3f} Å --> {d2 * BOHR2ANG:.3f} Å"

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

