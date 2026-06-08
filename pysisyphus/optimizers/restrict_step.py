import numpy as np


def scale_by_max_step(steps, max_step):
    steps_max = np.abs(steps).max()
    if steps_max > max_step:
        steps *= max_step / steps_max
    return steps


def get_scale_max(max_element):
    def scale_max(step):
        step_max = np.abs(step).max()
        if step_max > max_element:
            step *= max_element / step_max
        return step
    return scale_max


def restrict_step(steps, max_step):
    too_big = np.abs(steps) > max_step
    signs = np.sign(steps[too_big])
    steps[too_big] = signs * max_step
    return steps


# ---------------------------------------------------------------------------
# Sella backport 4: per-coord-type weighted L_inf trust check
#
# Source: sella/optimize/restricted_step.py:177-212 (MaxInternalStep class).
# Sella applies a weight per coordinate kind (bond / angle / dihedral / other
# / cartesian / rotation) before computing the L_inf norm. The default
# weights wb=1.0 / wa=0.5 / wd=0.25 give dihedrals 4x less trust budget
# than bonds, which fixes the dominant failure mode of pysis's plain L2
# norm where a single flexible dihedral consumes the whole step.
#
# This is a pure helper — caller (e.g., an opt-in branch inside RSI-RFO)
# decides when to substitute it for `np.linalg.norm`. Not yet wired into
# any optimizer dispatch.
# ---------------------------------------------------------------------------

import numpy as np
from pysisyphus.intcoords.PrimTypes import PrimTypes as _PT

_BOND_TYPES = {
    _PT.BOND, _PT.AUX_BOND, _PT.HYDROGEN_BOND,
    _PT.INTERFRAG_BOND, _PT.AUX_INTERFRAG_BOND,
}
_ANGLE_TYPES = {
    _PT.BEND, _PT.BEND2, _PT.LINEAR_BEND, _PT.LINEAR_BEND_COMPLEMENT,
}
_DIHEDRAL_TYPES = {
    _PT.PROPER_DIHEDRAL, _PT.PROPER_DIHEDRAL2, _PT.IMPROPER_DIHEDRAL,
    _PT.OUT_OF_PLANE, _PT.DUMMY_TORSION, _PT.DUMMY_IMPROPER,
    _PT.ROBUST_TORSION1, _PT.ROBUST_TORSION2,
}
_CARTESIAN_TYPES = {
    _PT.CARTESIAN, _PT.CARTESIAN_X, _PT.CARTESIAN_Y, _PT.CARTESIAN_Z,
    _PT.TRANSLATION, _PT.TRANSLATION_X, _PT.TRANSLATION_Y, _PT.TRANSLATION_Z,
    _PT.ROTATION, _PT.ROTATION_A, _PT.ROTATION_B, _PT.ROTATION_C,
}


def per_coord_type_weights(
    typed_prims,
    *,
    wb: float = 1.0,
    wa: float = 0.5,
    wd: float = 0.25,
    wo: float = 0.5,
    wx: float = 1.0,
):
    """Build a per-coordinate weight vector from a `typed_prims` list.

    Maps each pysisyphus PrimType into one of {bond, angle, dihedral, cart,
    other} and returns the corresponding weight. Anything unrecognised gets
    `wo` (other). Default weights match Sella's defaults for saddle search.

    Parameters
    ----------
    typed_prims : iterable of (PrimType, *atom_indices)
        Output of `RedundantCoords.typed_prims`.
    wb, wa, wd, wo, wx : float
        Weights for bonds / angles / dihedrals / other / cartesian-rotation.

    Returns
    -------
    weights : ndarray, shape (n_internal,)
    """
    out = []
    for entry in typed_prims:
        pt = entry[0]
        if pt in _BOND_TYPES:
            out.append(wb)
        elif pt in _ANGLE_TYPES:
            out.append(wa)
        elif pt in _DIHEDRAL_TYPES:
            out.append(wd)
        elif pt in _CARTESIAN_TYPES:
            out.append(wx)
        else:
            out.append(wo)
    return np.asarray(out, dtype=float)


def weighted_max_internal_step(
    step_internal,
    typed_prims,
    *,
    wb: float = 1.0,
    wa: float = 0.5,
    wd: float = 0.25,
    wo: float = 0.5,
    wx: float = 1.0,
) -> float:
    """Sella's per-coord-type weighted L_inf norm.

    Replaces a plain `np.linalg.norm(step)` (= L2) trust check inside an
    internal-coord RS-RFO loop with a weighted L_inf that discounts
    dihedrals before they can use up the entire trust budget.

    Parameters
    ----------
    step_internal : ndarray, shape (n_internal,)
        Step vector in redundant internal coords.
    typed_prims : list of (PrimType, *atom_indices)
        From `RedundantCoords.typed_prims`. len == len(step_internal).
    wb, wa, wd, wo, wx : float
        See `per_coord_type_weights`.

    Returns
    -------
    weighted_max : float
        max_i (|step_i| * weight(prim_type_i))
    """
    step = np.asarray(step_internal, dtype=float).ravel()
    if len(step) != len(typed_prims):
        raise ValueError(
            f"step has {len(step)} entries but typed_prims has "
            f"{len(typed_prims)}; they must match"
        )
    weights = per_coord_type_weights(
        typed_prims, wb=wb, wa=wa, wd=wd, wo=wo, wx=wx
    )
    return float(np.max(np.abs(step) * weights))
