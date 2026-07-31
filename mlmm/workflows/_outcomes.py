"""Product-local scientific-outcome vocabulary for mlmm_toolkit.

The vocabulary separates three orthogonal questions that the single legacy
``status`` field used to conflate:

* did a command or stage actually run?       -> ``executed``
* did the underlying engine converge?        -> ``converged`` (tri-state)
* may a scientific consumer use the result?  -> ``usable``

Fail-closed invariant: no artifact existence and no
finite fallback promotes a *required* nonconverged or missing scientific leaf to
success.  ``converged`` is tri-state: ``True`` means the engine reported
convergence, ``False`` means it explicitly did not, and ``None`` means the
convergence is unknown (a missing signal is treated as *not* converged).
Artifacts remain reportable even when the leaf that produced them is unusable.

These types are additive.  They never replace or rename the legacy public
``status``/``schema_version`` fields; they are serialized alongside them so a
forward-compatible consumers can read explicit outcomes while legacy
consumers that read only ``status`` are unaffected.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Any, Iterable, List, Mapping, Optional, Sequence, Tuple


def _finite(value: Any) -> bool:
    """Return True only for a finite real number (not None, NaN, or inf)."""

    if value is None or isinstance(value, bool):
        return False
    try:
        return math.isfinite(float(value))
    except (TypeError, ValueError):
        return False


def _normalize_bool(value: Any) -> Optional[bool]:
    """Return a Python boolean for explicit Python/NumPy boolean scalars."""

    if value is True or value is False:
        return value
    # NumPy 2 exposes ``np.bool_`` as ``numpy.bool`` while older releases use
    # ``numpy.bool_``. Accept only those scalar types, not arbitrary truthy
    # objects.
    value_type = type(value)
    if (
        value_type.__module__ == "numpy"
        and value_type.__name__ in {"bool", "bool_"}
    ):
        try:
            return bool(value)
        except (TypeError, ValueError):
            return None
    return None


def _is_true(value: Any) -> bool:
    """Explicit-True test that is robust to NumPy booleans.

    Convergence is fail-closed: Python and NumPy boolean ``True`` count as
    converged. Missing or non-boolean truthy values such as ``1`` and
    ``"true"`` do not.
    """

    return _normalize_bool(value) is True


@dataclass(frozen=True)
class LeafOutcome:
    """One scientific leaf (a stage / engine run / endpoint / segment).

    ``usable`` is the explicit scientific-consumption gate.  Prefer the
    :func:`make_leaf` factory, which computes ``usable`` from the fail-closed
    rule; the raw constructor is retained for tests and for callers that carry
    an engine-specific usability concept.
    """

    stage: str
    item_id: str
    required: bool = True
    executed: Optional[bool] = None
    converged: Optional[bool] = None
    usable: bool = False
    reason: str = ""
    artifacts: Tuple[str, ...] = ()

    def to_dict(self) -> dict:
        return {
            "stage": self.stage,
            "item_id": self.item_id,
            "required": bool(self.required),
            "executed": self.executed,
            "converged": self.converged,
            "usable": bool(self.usable),
            "reason": self.reason,
            "artifacts": list(self.artifacts),
        }


@dataclass(frozen=True)
class ScanPointOutcome:
    """One attempted grid/scan point.

    A point is ``seed_eligible`` only when its optimizer explicitly converged,
    its unbiased energy is finite, and its geometry artifact was written.  Failed
    points are retained (for CSV/diagnostics) but must never become a seed, a
    reference minimum, an interpolation input, or the reported minimum energy.
    """

    point_id: str
    executed: Optional[bool] = None
    converged: Optional[bool] = None
    energy_valid: bool = False
    artifact_written: bool = False
    seed_eligible: bool = False
    reason: str = ""

    def to_dict(self) -> dict:
        return {
            "point_id": self.point_id,
            "executed": self.executed,
            "converged": self.converged,
            "energy_valid": bool(self.energy_valid),
            "artifact_written": bool(self.artifact_written),
            "seed_eligible": bool(self.seed_eligible),
            "reason": self.reason,
        }


@dataclass(frozen=True)
class AggregateTruth:
    """The single aggregate truth for a multi-leaf workflow.

    ``execution_status`` answers whether the constituent commands ran;
    ``scientific_status`` answers whether the science is complete and usable.
    """

    execution_status: str
    scientific_status: str
    status_reasons: Tuple[str, ...] = ()
    expected_item_ids: Tuple[str, ...] = ()
    observed_item_ids: Tuple[str, ...] = ()

    def to_dict(self) -> dict:
        return {
            "execution_status": self.execution_status,
            "scientific_status": self.scientific_status,
            "status_reasons": list(self.status_reasons),
            "expected_item_ids": list(self.expected_item_ids),
            "observed_item_ids": list(self.observed_item_ids),
        }


def make_leaf(
    stage: str,
    item_id: str,
    *,
    required: bool = True,
    executed: Optional[bool],
    converged: Optional[bool],
    energy_valid: bool = True,
    artifacts: Sequence[str] = (),
    reason: str = "",
) -> LeafOutcome:
    """Construct a :class:`LeafOutcome`, computing ``usable`` fail-closed.

    ``usable`` is True only when the stage executed, the engine explicitly
    converged (``converged is True``), and any required numeric result is finite.
    Anything else is unusable and carries a diagnostic ``reason``.  ``converged``
    left as ``None`` (a missing/unknown signal) is *not* usable.
    """

    executed = _normalize_bool(executed)
    converged = _normalize_bool(converged)
    usable = executed is True and converged is True and bool(energy_valid)
    if not reason:
        if executed is False:
            reason = "not_executed"
        elif converged is not True:
            reason = "not_converged" if converged is False else "convergence_unknown"
        elif not energy_valid:
            reason = "energy_invalid"
        else:
            reason = "ok"
    return LeafOutcome(
        stage=stage,
        item_id=item_id,
        required=bool(required),
        executed=executed,
        converged=converged,
        usable=bool(usable),
        reason=reason,
        artifacts=tuple(str(a) for a in artifacts),
    )


def make_scan_point(
    point_id: str,
    *,
    executed: Optional[bool],
    converged: Optional[bool],
    energy: Any = None,
    energy_valid: Optional[bool] = None,
    artifact_written: bool = False,
    reason: str = "",
) -> ScanPointOutcome:
    """Construct a :class:`ScanPointOutcome` with fail-closed seed eligibility."""

    executed = _normalize_bool(executed)
    converged = _normalize_bool(converged)
    if energy_valid is None:
        energy_valid = _finite(energy)
    seed_eligible = _is_true(converged) and bool(energy_valid) and bool(artifact_written)
    if not reason:
        if not _is_true(converged):
            reason = "not_converged" if converged is False else "convergence_unknown"
        elif not energy_valid:
            reason = "energy_invalid"
        elif not artifact_written:
            reason = "artifact_missing"
        else:
            reason = "ok"
    return ScanPointOutcome(
        point_id=str(point_id),
        executed=executed,
        converged=converged,
        energy_valid=bool(energy_valid),
        artifact_written=bool(artifact_written),
        seed_eligible=bool(seed_eligible),
        reason=reason,
    )


def eligible_points(points: Iterable[ScanPointOutcome]) -> List[ScanPointOutcome]:
    """Return only the seed-eligible points (converged, finite, artifact written)."""

    return [p for p in points if p.seed_eligible]


def optimizer_converged_bit(optimizer: Any, attr: str = "is_converged") -> Optional[bool]:
    """Read a pysisyphus optimizer's convergence as a fail-closed tri-state bit.

    Returns ``True``/``False`` only when the optimizer exposes an *explicit*
    boolean ``is_converged``. NumPy boolean scalars are normalized to Python
    booleans; anything else (a missing attribute, a non-boolean value such as
    ``1``) collapses to ``None`` (unknown). A "normal return" is not
    convergence: the caller must read this bit rather than assume a non-raising
    ``run()`` converged.
    """

    return _normalize_bool(getattr(optimizer, attr, None))


def combine_step_convergence(convs: Iterable[Optional[bool]]) -> Optional[bool]:
    """Fold per-step convergence so an early failure is never hidden.

    Any explicit ``False`` -> ``False`` (a middle step failed, even if the final
    step converged); otherwise any ``None`` -> ``None`` (unknown, fail-closed);
    all ``True`` -> ``True``; no steps at all -> ``None``.
    """

    convs = [_normalize_bool(conv) for conv in convs]
    if any(c is False for c in convs):
        return False
    if any(c is None for c in convs):
        return None
    return True if convs else None


def irc_hessian_cache_eligible(obj: Any, converged_attr: str) -> bool:
    """Return whether an IRC endpoint Hessian is eligible for caching.

    A direction's endpoint Hessian may be cached only when that direction
    *explicitly* converged: ``getattr(obj, converged_attr, None) is True``.  A
    nonconverged (max-cycle) or never-run direction, or a non-boolean signal,
    fails closed so no stale/never-converged Hessian seeds a downstream RFO.
    """

    return _normalize_bool(getattr(obj, converged_attr, None)) is True


def irc_direction_leaves(
    directions: Iterable[Tuple[str, bool, Any, Any, Sequence[str]]],
) -> Tuple[List[LeafOutcome], List[str]]:
    """Build one :class:`LeafOutcome` per IRC direction plus the expected IDs.

    ``directions`` yields ``(name, requested, converged, n_frames, artifacts)``.
    A requested direction is a required leaf that is usable only when it
    explicitly converged and produced at least one frame; its endpoint
    trajectory is retained as an artifact regardless.  A disabled direction is
    an optional (not-required) leaf that contributes no failure — a one-sided
    IRC request is a legitimate success.
    """

    leaves: List[LeafOutcome] = []
    expected: List[str] = []
    for name, requested, converged, n_frames, artifacts in directions:
        if not requested:
            leaves.append(
                make_leaf(
                    "irc",
                    name,
                    required=False,
                    executed=False,
                    converged=None,
                    reason="direction_disabled",
                )
            )
            continue
        expected.append(name)
        try:
            _n = int(n_frames)
        except (TypeError, ValueError):
            _n = 0
        leaves.append(
            make_leaf(
                "irc",
                name,
                required=True,
                executed=True,
                converged=_normalize_bool(converged),
                energy_valid=_n > 0,
                artifacts=list(artifacts),
            )
        )
    return leaves, expected


def ipopt_status_to_converged(status: Any) -> Tuple[Optional[bool], str]:
    """Map an IPOPT return status to a fail-closed (converged, reason) pair.

    IPOPT ``status == 0`` is Solve_Succeeded and ``1`` is
    Solved_To_Acceptable_Level; both count as converged.  Any other integer
    (max-iter, infeasible, ...) is *not* convergence.  A missing/unreadable
    status returns ``(None, "convergence_unknown")`` so DMF artifact existence
    never promotes a nonconverged solve.
    """

    if status is None:
        return None, "convergence_unknown"
    try:
        code = int(status)
    except (TypeError, ValueError):
        return None, "convergence_unknown"
    converged = code in (0, 1)
    return converged, ("ipopt_converged" if converged else f"ipopt_status_{code}")


def seed_eligible_mask(
    records: Sequence[Mapping[str, Any]],
    *,
    converged_key: str = "bias_converged",
    energy_key: str = "energy_hartree",
) -> List[bool]:
    """Boolean mask over scan ``records`` marking seed-eligible rows.

    A row is eligible only when its recorded convergence bit is exactly
    ``True``, its energy is finite, and its geometry artifact was written.
    Legacy rows lacking any of these provenance fields are ineligible.
    """

    mask: List[bool] = []
    for rec in records:
        conv = rec.get(converged_key)
        mask.append(
            _is_true(conv)
            and _finite(rec.get(energy_key))
            and _is_true(rec.get("artifact_written"))
        )
    return mask


def scan_scientific_status(points: Sequence[ScanPointOutcome]) -> Tuple[str, Tuple[str, ...]]:
    """Map a set of scan points to (scientific_status, reasons).

    Zero usable points is ``failed``; a mixture of usable and failed points is
    ``partial``; every point usable is ``success``.
    """

    attempted = len(points)
    usable = sum(1 for p in points if p.seed_eligible)
    reasons: List[str] = []
    if attempted == 0:
        return "failed", ("no_points_attempted",)
    failed_ids = [p.point_id for p in points if not p.seed_eligible]
    if usable == 0:
        reasons.append("no_usable_points")
        status = "failed"
    elif failed_ids:
        reasons.append(f"unusable_points:{len(failed_ids)}")
        status = "partial"
    else:
        status = "success"
    return status, tuple(reasons)


def aggregate_workflow_truth(
    leaves: Iterable[LeafOutcome],
    expected_item_ids: Sequence[Any] = (),
) -> AggregateTruth:
    """The single aggregate truth mapper for a multi-leaf workflow.

    Rules:

    * every expected *required* leaf usable and no expected ID missing -> success
    * at least one usable required leaf or usable diagnostic artifact, plus any
      missing / unusable required leaf -> partial
    * no usable required output, or an execution failure with nothing usable ->
      failed

    A required leaf whose ``converged`` is not ``True`` or whose ``usable`` is
    ``False`` cannot count toward completeness, regardless of any artifact it
    produced.
    """

    leaves = list(leaves)
    expected = tuple(dict.fromkeys(str(x) for x in expected_item_ids))
    observed = tuple(dict.fromkeys(leaf.item_id for leaf in leaves))
    required = [leaf for leaf in leaves if leaf.required]
    usable_required = [leaf for leaf in required if leaf.usable]
    observed_required_ids = {leaf.item_id for leaf in required}

    missing = [item for item in expected if item not in observed_required_ids]
    exec_failed = any(_normalize_bool(leaf.executed) is False for leaf in required)

    reasons: List[str] = []
    for leaf in required:
        if not leaf.usable:
            reasons.append(f"{leaf.stage}:{leaf.item_id}:{leaf.reason or 'unusable'}")
    for item in missing:
        reasons.append(f"missing:{item}")
    reasons = list(dict.fromkeys(reasons))

    execution_status = "failed" if exec_failed else "completed"

    has_usable_required = bool(usable_required)
    has_any_artifact = any(leaf.artifacts for leaf in leaves)
    all_required_usable = bool(required) and all(leaf.usable for leaf in required)

    if all_required_usable and not missing:
        scientific = "success"
    elif has_usable_required or has_any_artifact:
        scientific = "partial"
    else:
        scientific = "failed"

    if exec_failed and not has_usable_required:
        scientific = "failed"

    return AggregateTruth(
        execution_status=execution_status,
        scientific_status=scientific,
        status_reasons=tuple(reasons),
        expected_item_ids=expected,
        observed_item_ids=observed,
    )


def attach_outcomes(
    data: dict,
    *,
    truth: Optional[AggregateTruth] = None,
    stage_outcomes: Optional[Sequence[LeafOutcome]] = None,
    point_outcomes: Optional[Sequence[ScanPointOutcome]] = None,
    scientific_status: Optional[str] = None,
    scientific_status_reasons: Optional[Sequence[str]] = None,
    write_reasons_key: str = "scientific_status_reasons",
    write_status_reasons: bool = True,
) -> dict:
    """Add the outcome fields to a result/summary dict in place.

    This never touches the legacy ``status`` / ``schema_version`` keys.  It adds
    ``execution_status``, ``scientific_status``, expected/observed ID lists, and
    serialized ``stage_outcomes`` / ``point_outcomes``.  Reasons are written to a
    distinct key (default ``scientific_status_reasons``) so a workflow that
    already owns a ``status_reasons`` field is not clobbered.
    """

    if truth is not None:
        data["execution_status"] = truth.execution_status
        data["scientific_status"] = truth.scientific_status
        data["expected_item_ids"] = list(truth.expected_item_ids)
        data["observed_item_ids"] = list(truth.observed_item_ids)
        if write_status_reasons and truth.status_reasons:
            data[write_reasons_key] = list(truth.status_reasons)
    if scientific_status is not None:
        data["scientific_status"] = scientific_status
    if scientific_status_reasons is not None and write_status_reasons:
        reasons = list(scientific_status_reasons)
        if reasons:
            data[write_reasons_key] = reasons
    if stage_outcomes is not None:
        data["stage_outcomes"] = [o.to_dict() for o in stage_outcomes]
    if point_outcomes is not None:
        data["point_outcomes"] = [o.to_dict() for o in point_outcomes]
    return data


__all__ = [
    "LeafOutcome",
    "ScanPointOutcome",
    "AggregateTruth",
    "make_leaf",
    "make_scan_point",
    "eligible_points",
    "seed_eligible_mask",
    "scan_scientific_status",
    "aggregate_workflow_truth",
    "attach_outcomes",
    "optimizer_converged_bit",
    "combine_step_convergence",
    "irc_hessian_cache_eligible",
    "irc_direction_leaves",
    "ipopt_status_to_converged",
]
