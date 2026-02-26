# mlmm_toolkit/align_freeze_atoms.py

"""
align_freeze_atoms — Rigid alignment and staged “scan + relaxation” utilities for pysisyphus Geometry objects
=====================================================================

Usage (API)
-----
    from mlmm_toolkit.align_freeze_atoms import (
        align_and_refine_pair_inplace,
        align_and_refine_sequence_inplace,
        align_second_to_first_kabsch_inplace,
        kabsch_R_t,
        scan_freeze_atoms_toward_target_inplace,
    )

Examples::
    >>> from mlmm_toolkit.align_freeze_atoms import align_and_refine_pair_inplace
    >>> result = align_and_refine_pair_inplace(g_ref, g_mob, shared_calc=mlmm_calc)
    >>> result["align"]["mode"]
    'kabsch'

Description
-----
API-only utilities to co-align and refine pre-optimized `pysisyphus.Geometry` objects—typically adjacent
images along a reaction path—using `freeze_atoms`. A rigid alignment is performed first (with special
handling when 1 or 2 atoms are frozen), followed by a staged “scan toward reference + local relaxation”
that moves the frozen atoms toward the reference in small steps while relaxing the surroundings.
All updates are applied in place on the mobile geometry.

The module now targets **ML/MM-only** workflows: a shared ML/MM calculator instance (created via
`mlmm_toolkit.mlmm_calc.mlmm`) must already be attached to the geometries or supplied through the
`shared_calc` keyword.

Provided functionality (concise):
- kabsch_R_t(P, Q)
    Row-vector Kabsch solver for two (N, 3) point sets. Returns (R, t) minimizing ||(Q @ R + t) − P||.
- align_second_to_first_kabsch_inplace(g_ref, g_mob, *, verbose=True)
    Rigidly aligns `g_mob` to `g_ref` in place. Special cases:
      - freeze_atoms length = 1: translate to match the anchor; restrict rotations about that anchor;
        minimizes all-atom RMSD; RMSD reported on all atoms.
      - freeze_atoms length = 2: align the axis defined by the two anchors; optimize rotation around
        that axis; RMSD reported on all atoms.
      - otherwise: Kabsch on the union of `freeze_atoms` from both geometries (or all atoms if empty);
        apply the transform to all atoms; RMSD reported on the Kabsch selection.
    Returns: dict(before_A, after_A, n_used, mode) with RMSD in Å and mode ∈ {"one_anchor","two_anchor","kabsch"}.
- scan_freeze_atoms_toward_target_inplace(
      g_ref, g_mob, *, step_A=0.1, per_step_cycles=50, final_cycles=200, max_steps=1000,
      shared_calc=None, out_dir=Path("./result_align_refine/"), thresh="gau", verbose=True)
    Moves `g_mob.freeze_atoms` toward `g_ref` by `step_A` Å per iteration. At each step the frozen atoms
    are held fixed and the surroundings are relaxed with LBFGS. When the maximum remaining distance is
    below one step, enforce exact coincidence of the frozen atoms and run a finishing relaxation.
    Returns: dict(max_remaining_A, n_steps, converged).
- align_and_refine_pair_inplace(
      g_ref, g_mob, *, shared_calc=None, out_dir=Path("./result_align_refine/"),
      step_A=0.1, per_step_cycles=50, final_cycles=200, max_steps=1000,
      thresh="gau", verbose=True)
    High-level pair API: (1) rigid alignment, then (2) scan + relaxation toward the reference.
    Returns: {"align": {...}, "scan": {...}}.
- align_and_refine_sequence_inplace(
      geoms, *, shared_calc=None, out_dir=Path("./result_align_refine/"),
      step_A=0.1, per_step_cycles=1000, final_cycles=1000, max_steps=10000,
      thresh="gau", verbose=True)
    Applies the pair procedure along [g0, g1, g2, ...] as (g0←g1), (g1←g2), ... and returns a list of per-pair results.

Outputs (& Directory Layout)
-----
- Default base directory: `./result_align_refine/`
  - Pair API: uses `out_dir/` for stepwise and final LBFGS runs.
  - Sequence API: creates subdirectories `pair_00/`, `pair_01/`, ... under `out_dir/`, one per adjacent pair.
- Directories are created even when optimizer dumping is disabled; LBFGS is invoked with `dump=False`.

Notes:
-----
- Units & conventions:
  - User-facing distances/thresholds are in Å; internal coordinates are in bohr (conversion via `BOHR2ANG`).
  - Row-vector convention for rigid transforms: points `Q` are mapped as `Q @ R + t`.
  - Indices in `freeze_atoms` are 0-based.
- Calculator handling:
  - A shared ML/MM calculator must already be attached to each geometry (e.g., via `mlmm(...)`) or supplied through
    `shared_calc`. UMA/link-atom fallbacks are not available.
- Rigid alignment priority and RMSD reporting:
  - freeze=1 → rotations about the single anchor only; RMSD minimized and reported on all atoms.
  - freeze=2 → align anchor-defined axis; optimize rotation around that axis; RMSD reported on all atoms.
  - else → Kabsch on the union of `freeze_atoms` (or all atoms if empty); transform applied to all atoms;
    RMSD reported on the same selection used by Kabsch.
- Scan + relaxation details:
  - At each step: move frozen atoms by `step_A` Å toward reference, hold them fixed, run a short LBFGS on the surroundings.
  - Finalization: enforce exact coincidence of frozen atoms, then run a finishing LBFGS.
  - The original `freeze_atoms` of `g_mob` is restored after the scan, even on exceptions.
- Error handling:
  - `LBFGS` exceptions (`ZeroStepLength`, `OptimizationError`) are caught and logged; the procedure continues when reasonable.
- Internals:
  - Uses Rodrigues-based rotations, vector-alignment, and planar projection helpers to implement axis-constrained motions.
"""


from __future__ import annotations

from pathlib import Path
from typing import List, Optional, Sequence, Tuple, Dict, Any

import click
import numpy as np

# pysisyphus
from pysisyphus.optimizers.LBFGS import LBFGS
from pysisyphus.optimizers.exceptions import OptimizationError, ZeroStepLength
from pysisyphus.constants import BOHR2ANG


# =============================================================================
# Math utilities (row-vector convention: Q @ R + t)
# =============================================================================

def kabsch_R_t(P: np.ndarray, Q: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    """
    Row-vector Kabsch: find R, t minimizing ||(Q @ R + t) - P|| for (N, 3).
    """
    P = np.asarray(P, float)
    Q = np.asarray(Q, float)
    if P.shape != Q.shape or P.ndim != 2 or P.shape[1] != 3:
        raise ValueError("Kabsch expects P, Q with shape (N, 3).")
    mu_P, mu_Q = P.mean(0), Q.mean(0)
    Pc, Qc = P - mu_P, Q - mu_Q
    H = Pc.T @ Qc
    U, _, Vt = np.linalg.svd(H)
    R = Vt.T @ U.T
    if np.linalg.det(R) < 0.0:
        Vt[-1] *= -1.0
        R = Vt.T @ U.T
    t = mu_P - mu_Q @ R
    return R, t


def _rodrigues(axis_unit: np.ndarray, theta: float) -> np.ndarray:
    """
    3×3 rotation matrix for a rotation of angle `theta` about `axis_unit`.
    """
    u = np.asarray(axis_unit, float)
    n = np.linalg.norm(u)
    if n < 1e-16:
        return np.eye(3)
    u /= n
    ux, uy, uz = u
    K = np.array([[0.0, -uz,  uy],
                  [uz,  0.0, -ux],
                  [-uy, ux,  0.0]], float)
    I = np.eye(3)
    return I + np.sin(theta) * K + (1.0 - np.cos(theta)) * (K @ K)


def _rotation_align_vectors(a: np.ndarray, b: np.ndarray) -> np.ndarray:
    """
    3×3 rotation (column-vector convention) that maps vector `a` to `b`.
    When applying to row-vectors, use the transpose of the returned matrix.
    """
    a = np.asarray(a, float)
    b = np.asarray(b, float)
    an = np.linalg.norm(a)
    bn = np.linalg.norm(b)
    if an < 1e-16 or bn < 1e-16:
        return np.eye(3)
    a /= an
    b /= bn
    v = np.cross(a, b)
    c = float(np.clip(a.dot(b), -1.0, 1.0))
    s = np.linalg.norm(v)
    if s < 1e-12:
        if c > 0.0:
            return np.eye(3)
        # 180°: choose any axis orthogonal to `a`
        tmp = np.array([1.0, 0.0, 0.0]) if abs(a[0]) <= 0.9 else np.array([0.0, 1.0, 0.0])
        axis = np.cross(a, tmp)
        axis /= (np.linalg.norm(axis) + 1e-16)
        return _rodrigues(axis, np.pi)
    axis = v / s
    theta = np.arctan2(s, c)
    return _rodrigues(axis, theta)


def _orth_proj_perp(u: np.ndarray) -> np.ndarray:
    """
    3×3 projector onto the plane perpendicular to vector `u`.
    """
    u = np.asarray(u, float)
    n = np.linalg.norm(u)
    if n < 1e-16:
        return np.eye(3)
    u = u / n
    return np.eye(3) - np.outer(u, u)


# =============================================================================
# Geometry utilities
# =============================================================================

def _coords3d(geom) -> np.ndarray:
    """
    Return (N, 3) coordinates in bohr (float).
    """
    return np.array(geom.coords3d, float)


def _rmsd(A: np.ndarray, B: np.ndarray) -> float:
    """
    Compute RMSD in Å (inputs are in bohr; conversion is applied internally).
    """
    A = np.asarray(A, float) * BOHR2ANG
    B = np.asarray(B, float) * BOHR2ANG
    return float(np.sqrt(np.mean(np.sum((A - B) ** 2, axis=1)))) if len(A) else float("nan")


def _set_all_coords_disabling_freeze(geom, coords3d_bohr: np.ndarray) -> None:
    """
    Temporarily clear `freeze_atoms`, set all coordinates (bohr), then restore.
    """
    old = np.array(getattr(geom, "freeze_atoms", []), int)
    try:
        geom.freeze_atoms = np.array([], int)
        geom.set_coords(np.asarray(coords3d_bohr, float).reshape(-1), cartesian=True)
    finally:
        geom.freeze_atoms = old


def _attach_calc_if_needed(geom, shared_calc=None) -> None:
    """
    Ensure a calculator is present; prefer `shared_calc` when available.

    ML/MM-only workflows mean no automatic fallback is provided—callers must either supply
    `shared_calc` (typically `mlmm(...)`) or attach a calculator beforehand.
    """
    try:
        has_calc = getattr(geom, "calculator", None) is not None
    except Exception:
        has_calc = False
    if shared_calc is not None:
        geom.set_calculator(shared_calc)
    elif not has_calc:
        raise RuntimeError(
            "align_freeze_atoms requires a calculator. Provide shared_calc=mlmm(...) or attach one to each geometry."
        )


def _freeze_union(g_ref, g_mob, n_atoms: Optional[int] = None) -> List[int]:
    """
    Union of `freeze_atoms` from `g_ref` and `g_mob` (0-based).
    If `n_atoms` is given, out-of-range indices are removed. Returns [] if empty.
    """
    fa0 = getattr(g_ref, "freeze_atoms", np.array([], int))
    fa1 = getattr(g_mob, "freeze_atoms", np.array([], int))
    cand = sorted(set(int(i) for i in list(fa0) + list(fa1)))
    if n_atoms is None:
        return cand
    good = [i for i in cand if 0 <= i < int(n_atoms)]
    return good


# =============================================================================
# Rigid alignment: special handling for freeze=1/2 → Kabsch otherwise
# =============================================================================

def align_second_to_first_kabsch_inplace(g_ref, g_mob,
                                         *, verbose: bool = True) -> Dict[str, Any]:
    """
    Rigidly align `g_mob` to `g_ref` (update coordinates of `g_mob` in place).

    Returns a dict with keys: {before_A, after_A, n_used, mode}
      - mode: "one_anchor" | "two_anchor" | "kabsch"

    Behavior:
      * freeze=1 atom: minimize the all-atom RMSD using rotations about the
        single anchor point only.
      * freeze=2 atoms: align the axis defined by the two anchors and optimize
        the rotation around that axis.
      * otherwise: solve Kabsch with the union of `freeze_atoms` (or all atoms
        if empty) and apply the resulting transform to all atoms.
    """
    P = _coords3d(g_ref)  # bohr
    Q = _coords3d(g_mob)  # bohr
    if P.shape != Q.shape:
        raise ValueError(f"Different atom counts: {P.shape[0]} vs {Q.shape[0]}")
    N = P.shape[0]
    idx = _freeze_union(g_ref, g_mob, n_atoms=N)

    def _set_all(Q_new: np.ndarray) -> None:
        _set_all_coords_disabling_freeze(g_mob, Q_new)

    mode = "kabsch"

    # ---- 1 anchor ----
    if len(idx) == 1:
        i = idx[0]
        before = _rmsd(P, Q)  # all-atom RMSD (design: constrained rotation minimizes all-atom)
        p0, q0 = P[i].copy(), Q[i].copy()
        Q_shift = Q + (p0 - q0)                # match the anchor point
        P_rel, Q_rel = P - p0, Q_shift - p0
        U, _, Vt = np.linalg.svd(P_rel.T @ Q_rel)
        R = Vt.T @ U.T
        if np.linalg.det(R) < 0.0:
            Vt[-1] *= -1.0
            R = Vt.T @ U.T
        Q_aln = (Q_rel @ R) + p0
        _set_all(Q_aln)
        after = _rmsd(P, Q_aln)
        mode = "one_anchor"
        if verbose:
            click.echo(f"[align] one-anchor: RMSD {before:.6f} Å → {after:.6f} Å (idx={i})")
        return {"before_A": before, "after_A": after, "n_used": 1, "mode": mode}

    # ---- 2 anchors ----
    if len(idx) == 2:
        before = _rmsd(P, Q)  # all-atom RMSD (design: constrained rotation minimizes all-atom)
        i0, i1 = idx[0], idx[1]
        p0, p1, q0, q1 = P[i0].copy(), P[i1].copy(), Q[i0].copy(), Q[i1].copy()
        pm, qm = 0.5 * (p0 + p1), 0.5 * (q0 + q1)
        vP, vQ = p1 - p0, q1 - q0

        if np.linalg.norm(vP) < 1e-16 or np.linalg.norm(vQ) < 1e-16:
            # Fallback: Kabsch
            pass
        else:
            Q0 = Q + (pm - qm)                      # match midpoints
            R_align = _rotation_align_vectors(vQ, vP)
            Q0 = ((Q0 - pm) @ R_align.T) + pm       # align axis direction (right-multiply → .T)

            u = vP / (np.linalg.norm(vP) + 1e-16)
            c = pm
            P_perp = _orth_proj_perp(u)
            A = (P - c) @ P_perp.T
            B = (Q0 - c) @ P_perp.T
            cross_u_B = np.cross(u, B)
            s1 = float(np.sum(A * B))
            s2 = float(np.sum(A * cross_u_B))
            theta = np.arctan2(s2, s1) if (abs(s1) + abs(s2)) > 1e-16 else 0.0

            R_axis = _rodrigues(u, theta)
            Q1 = ((Q0 - c) @ R_axis.T) + c
            _set_all(Q1)
            after = _rmsd(P, Q1)
            mode = "two_anchor"
            if verbose:
                click.echo(f"[align] two-anchors: RMSD {before:.6f} Å → {after:.6f} Å (idx=({i0},{i1}))")
            return {"before_A": before, "after_A": after, "n_used": 2, "mode": mode}

    # ---- Default: Kabsch (selected freeze atoms or all atoms) ----
    if len(idx) > 0:
        use = np.zeros(N, bool)
        for k in idx:
            if 0 <= k < N:
                use[k] = True
    else:
        use = np.ones(N, bool)

    P_sel, Q_sel = P[use], Q[use]
    n_used = int(P_sel.shape[0])

    # *** CHANGED: evaluate RMSD on the same selection used for Kabsch ***
    before_sel = _rmsd(P_sel, Q_sel)

    R, t = kabsch_R_t(P_sel, Q_sel)
    Q_aln = (Q @ R) + t
    _set_all(Q_aln)

    after_sel = _rmsd(P_sel, Q_aln[use])

    if verbose:
        click.echo(f"[align] kabsch:     RMSD {before_sel:.6f} Å → {after_sel:.6f} Å (used {n_used})")

    return {"before_A": before_sel, "after_A": after_sel, "n_used": n_used, "mode": mode}


# =============================================================================
# Scan + relaxation (stepwise matching of `freeze_atoms` to the reference)
# =============================================================================

def scan_freeze_atoms_toward_target_inplace(
    g_ref,
    g_mob,
    *,
    step_A: float = 0.1,
    per_step_cycles: int = 50,
    final_cycles: int = 200,
    max_steps: int = 1000,
    thresh: str = "gau",
    shared_calc=None,
    out_dir: Path = Path("./result_align_refine/"),
    verbose: bool = True,
) -> Dict[str, Any]:
    """
    Move the `freeze_atoms` of `g_mob` toward the reference (`g_ref`) by `step_A` Å
    per iteration. At each step, keep the frozen atoms fixed and run a short LBFGS
    relaxation on the remaining atoms. Finally, enforce exact coincidence of the
    frozen atoms and perform a finishing relaxation. Updates `g_mob` in place.

    Returns:
        dict(max_remaining_A, n_steps, converged)
    """
    P = _coords3d(g_ref)  # bohr
    Q = _coords3d(g_mob)  # bohr
    if P.shape != Q.shape:
        raise ValueError(f"Different atom counts: {P.shape[0]} vs {Q.shape[0]}")
    N = P.shape[0]

    original_freeze = np.array(getattr(g_mob, "freeze_atoms", []), int)
    try:
        idx = _freeze_union(g_ref, g_mob, n_atoms=N)

        if len(idx) == 0:
            if verbose:
                click.echo("[scan] freeze_atoms list is empty. Skipping scan and relaxation.")
            return {"max_remaining_A": 0.0, "n_steps": 0, "converged": True}

        # Attach a calculator if needed
        _attach_calc_if_needed(g_mob, shared_calc)

        out_dir = Path(out_dir)
        out_dir.mkdir(parents=True, exist_ok=True)

        step_bohr = float(step_A) / BOHR2ANG
        eps = 1e-12
        n_steps_done = 0
        converged = False
        max_remaining_A = None

        for istep in range(1, max_steps + 1):
            Q = _coords3d(g_mob)
            d = P[idx] - Q[idx]                    # bohr
            rem_bohr = np.linalg.norm(d, axis=1)
            max_rem_bohr = float(rem_bohr.max()) if len(rem_bohr) else 0.0
            max_remaining_A = max_rem_bohr * BOHR2ANG
            if verbose:
                click.echo(f"[scan] step {istep:03d}: max remaining = {max_remaining_A:.6f} Å")

            if max_rem_bohr <= step_bohr + 1e-12:
                # Final step: enforce exact coincidence
                Q_new = Q.copy()
                Q_new[idx] = P[idx]
                _set_all_coords_disabling_freeze(g_mob, Q_new)
                try:
                    # Finishing relaxation
                    g_mob.freeze_atoms = np.array(idx, int)
                    LBFGS(
                        g_mob,
                        out_dir=str(out_dir),
                        max_cycles=int(final_cycles),
                        print_every=100,
                        thresh=thresh,
                        dump=False,
                    ).run()
                except (ZeroStepLength, OptimizationError) as e:
                    if verbose:
                        click.echo(f"[scan] WARNING: Exception occurred in final relaxation: {e} (continue...)", err=True)
                g_mob.freeze_atoms = np.array([], int)
                converged = True
                n_steps_done = istep
                break

            # Take one step forward toward the target
            move = np.zeros_like(d)
            sel = rem_bohr > eps
            move[sel] = (d[sel] / rem_bohr[sel, None]) * step_bohr
            Q_next = Q.copy()
            Q_next[idx] = Q[idx] + move

            # Update coordinates → short relaxation with frozen atoms fixed
            _set_all_coords_disabling_freeze(g_mob, Q_next)
            try:
                g_mob.freeze_atoms = np.array(idx, int)
                LBFGS(
                    g_mob,
                    out_dir=str(out_dir),
                    max_cycles=int(per_step_cycles),
                    print_every=100,
                    thresh=thresh,
                    dump=False,
                ).run()
            except (ZeroStepLength, OptimizationError) as e:
                if verbose:
                    click.echo(f"[scan] WARNING: Exception occurred in relaxation: {e} (continue...)", err=True)
            finally:
                g_mob.freeze_atoms = np.array([], int)

            n_steps_done = istep
        else:
            if verbose:
                click.echo(f"[scan] WARNING: Reached max_steps={max_steps}.", err=True)

        return {"max_remaining_A": float(max_remaining_A or 0.0),
                "n_steps": int(n_steps_done),
                "converged": bool(converged)}
    finally:
        # Always restore the original freeze_atoms state
        g_mob.freeze_atoms = original_freeze


# =============================================================================
# High-level API: pair / sequence
# =============================================================================

def align_and_refine_pair_inplace(
    g_ref,
    g_mob,
    *,
    shared_calc=None,
    out_dir: Path = Path("./result_align_refine/"),
    step_A: float = 0.1,
    per_step_cycles: int = 50,
    final_cycles: int = 200,
    max_steps: int = 1000,
    thresh: str = "gau",
    verbose: bool = True,
) -> Dict[str, Any]:
    """
    For a pair (g_ref, g_mob), perform:
      (1) rigid alignment (special cases for freeze=1/2, otherwise Kabsch), then
      (2) scan + relaxation (stepwise matching of `freeze_atoms` to the reference).
    Updates `g_mob` in place.

    Returns:
        {
          "align": {before_A, after_A, n_used, mode},
          "scan":  {max_remaining_A, n_steps, converged}
        }
    """
    # Rigid alignment
    align_res = align_second_to_first_kabsch_inplace(g_ref, g_mob, verbose=verbose)
    # Scan + relaxation
    scan_res = scan_freeze_atoms_toward_target_inplace(
        g_ref, g_mob,
        step_A=step_A,
        per_step_cycles=per_step_cycles,
        final_cycles=final_cycles,
        max_steps=max_steps,
        thresh=thresh,
        shared_calc=shared_calc,
        out_dir=out_dir,
        verbose=verbose,
    )
    return {"align": align_res, "scan": scan_res}


def align_and_refine_sequence_inplace(
    geoms: Sequence[Any],
    *,
    shared_calc=None,
    out_dir: Path = Path("./result_align_refine/"),
    step_A: float = 0.1,
    per_step_cycles: int = 1000,
    final_cycles: int = 1000,
    max_steps: int = 10000,
    thresh: str = "gau",
    verbose: bool = True,
) -> List[Dict[str, Any]]:
    """
    For a list [g0, g1, g2, ...], apply `align_and_refine_pair_inplace` in order:
    (g0←g1), (g1←g2), ... i.e., each g_{i+1} is aligned/refined to g_i.
    Returns a list of per-pair result dicts.

    Intended to be used after pre-optimization in `path_search.py`.
    """
    geoms = list(geoms)
    if len(geoms) <= 1:
        return []

    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    results: List[Dict[str, Any]] = []
    for i in range(len(geoms) - 1):
        g_ref = geoms[i]
        g_mob = geoms[i + 1]
        pair_out = out_dir / f"pair_{i:02d}"
        pair_out.mkdir(parents=True, exist_ok=True)

        if verbose:
            click.echo(f"\n[align+scan] Pair {i:02d}: image {i} (ref) ← image {i+1} (mobile)")

        res = align_and_refine_pair_inplace(
            g_ref, g_mob,
            shared_calc=shared_calc,
            out_dir=pair_out,
            step_A=step_A,
            per_step_cycles=per_step_cycles,
            final_cycles=final_cycles,
            max_steps=max_steps,
            thresh=thresh,
            verbose=verbose,
        )
        results.append(res)

    return results


__all__ = [
    "align_second_to_first_kabsch_inplace",
    "scan_freeze_atoms_toward_target_inplace",
    "align_and_refine_pair_inplace",
    "align_and_refine_sequence_inplace",
    "kabsch_R_t",
]
