"""Product-local immutable microiteration partition + outcome vocabulary (C9).

This private helper is shared by ``mlmm.workflows.opt`` and
``mlmm.workflows.tsopt`` so both drivers resolve *one* immutable macro/micro
partition and serialize *one* truthful optimizer outcome.  It is intentionally
product-local (an MLMM-private minor helper), not a PEScape public
``MicroiterationPartition``/``ActiveDofMap``; p2r has no ML/MM layer partition
and does not import it.

Three invariants this module enforces (the point of C9):

* **Immutable partition.** The user's original freeze mask is preserved in BOTH
  the macro and the micro phase (it is unioned into each phase freeze and never
  appears in either phase active) and is restored exactly on every exit path.
  A user-frozen atom therefore never moves in either phase.
* **Loud partition failure.** A partition that cannot be resolved from the
  accepted calculator core raises :class:`PartitionError`; it is never a
  swallowed empty set that a caller mistakes for a valid empty-ML region.  A
  valid empty-ML region is a *documented* fallback to the ordinary optimizer,
  not an unchanged geometry returned as a completed result.
* **Truthful macro/micro outcome.** Aggregate convergence requires an explicitly
  converged macro state *and* an explicitly converged latest required micro
  relaxation.  A micro plateau/stall/max-cycle exhaustion (or a missing
  convergence signal) fails closed and never reads as macro convergence.

The :class:`OptimizerOutcome` field names (``status``/``executed``/
``converged``/``cycles``/``max_cycles``/``stalled``/``stop_reason``) are kept
field-isomorphic with the C6/C7 vocabulary so a later PEScape migration can map
the frozen minor vectors into its public conformance without exposing these
private classes as a prematurely frozen API.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Dict, Iterable, List, Optional, Sequence, Tuple

from mlmm.core.utils import optimizer_cycle_count, optimizer_terminal_status
from mlmm.workflows._outcomes import optimizer_converged_bit


class PartitionError(RuntimeError):
    """A microiteration partition could not be resolved from the accepted core.

    Raised (never swallowed) when the calculator core does not expose a usable
    layer partition or when the resolved indices are out of range / overlapping.
    A caller must let this propagate to the ordinary error envelope; it must not
    return an unchanged geometry as if it were an optimized result.
    """


# ---------------------------------------------------------------------------
# Immutable macro/micro partition
# ---------------------------------------------------------------------------


def _clean_indices(values: Iterable[Any], n_atoms: int, *, label: str) -> Tuple[int, ...]:
    """Return validated, de-duplicated 0-based atom indices in ascending order."""

    out: List[int] = []
    seen: set[int] = set()
    for v in values:
        try:
            idx = int(v)
        except (TypeError, ValueError) as exc:
            raise PartitionError(f"{label}: non-integer atom index {v!r}.") from exc
        if idx < 0 or idx >= n_atoms:
            raise PartitionError(
                f"{label}: atom index {idx} out of range [0, {n_atoms})."
            )
        if idx not in seen:
            seen.add(idx)
            out.append(idx)
    return tuple(sorted(out))


@dataclass(frozen=True)
class MicroiterationPartition:
    """One immutable macro/micro partition shared by opt and tsopt.

    Every field is an ordered tuple derived by iterating ``range(n_atoms)`` (never
    by exposing set iteration order).  For each phase the active/freeze pair is
    disjoint and covers all atoms, and ``original_freeze`` is a subset of both
    ``macro_freeze_atoms`` and ``micro_freeze_atoms``.
    """

    n_atoms: int
    original_freeze: Tuple[int, ...]
    macro_active_atoms: Tuple[int, ...]
    macro_freeze_atoms: Tuple[int, ...]
    micro_active_atoms: Tuple[int, ...]
    micro_freeze_atoms: Tuple[int, ...]
    ml_atoms: Tuple[int, ...]
    link_parent_atoms: Tuple[int, ...]
    movable_mm_atoms: Tuple[int, ...]

    @property
    def has_macro_active(self) -> bool:
        """True when the macro phase has at least one movable (ML/link) atom."""

        return len(self.macro_active_atoms) > 0

    @property
    def has_micro_active(self) -> bool:
        """True when the micro phase has at least one movable MM atom."""

        return len(self.micro_active_atoms) > 0

    def provenance(self) -> Dict[str, Any]:
        """A small JSON-friendly provenance record of the partition sizes."""

        return {
            "n_atoms": int(self.n_atoms),
            "n_original_freeze": len(self.original_freeze),
            "n_ml_atoms": len(self.ml_atoms),
            "n_link_parent_atoms": len(self.link_parent_atoms),
            "n_movable_mm_atoms": len(self.movable_mm_atoms),
            "n_macro_active": len(self.macro_active_atoms),
            "n_micro_active": len(self.micro_active_atoms),
        }


def build_partition(
    n_atoms: int,
    *,
    ml: Iterable[Any],
    link_parents: Iterable[Any],
    hess_mm: Iterable[Any],
    movable_mm: Iterable[Any],
    frozen_mm: Iterable[Any],
    original_freeze: Iterable[Any],
) -> MicroiterationPartition:
    """Build the single immutable macro/micro partition.

    * macro candidate  = ML u link-parent
    * micro candidate  = (HessMM u MovableMM) \\ link-parent
    * each phase active = candidate minus the user's original freeze
    * each phase freeze = the ordered complement of that phase's active

    The user's ``original_freeze`` is therefore removed from both actives and
    present in both freezes; a frozen atom never moves in either phase.
    """

    if int(n_atoms) <= 0:
        raise PartitionError(f"n_atoms must be positive (got {n_atoms}).")
    n_atoms = int(n_atoms)

    ml_t = _clean_indices(ml, n_atoms, label="ml")
    link_t = _clean_indices(link_parents, n_atoms, label="link_parents")
    hess_mm_t = _clean_indices(hess_mm, n_atoms, label="hess_mm")
    movable_mm_t = _clean_indices(movable_mm, n_atoms, label="movable_mm")
    frozen_t = _clean_indices(original_freeze, n_atoms, label="original_freeze")
    frozen_layer_t = _clean_indices(frozen_mm, n_atoms, label="frozen_mm")

    # ``original_freeze`` is the exact accepted geometry mask (the restore
    # target). ``effective_freeze`` additionally holds frozen-layer atoms so they
    # never move in either phase even when a caller passed them only through the
    # calculator layer; ``original_freeze`` stays a subset of it, so restoring
    # the exact geometry mask is always well defined.
    original_set = set(frozen_t)
    effective_set = original_set | set(frozen_layer_t)
    link_set = set(link_t)

    macro_candidate = set(ml_t) | link_set
    micro_candidate = (set(hess_mm_t) | set(movable_mm_t)) - link_set

    macro_active = tuple(i for i in range(n_atoms) if i in macro_candidate and i not in effective_set)
    macro_freeze = tuple(i for i in range(n_atoms) if i not in set(macro_active))
    micro_active = tuple(i for i in range(n_atoms) if i in micro_candidate and i not in effective_set)
    micro_freeze = tuple(i for i in range(n_atoms) if i not in set(micro_active))

    # Disjoint + covering invariant for both phases (a validated partition), and
    # the user freeze mask preserved in each phase freeze.
    for phase, active, freeze in (
        ("macro", macro_active, macro_freeze),
        ("micro", micro_active, micro_freeze),
    ):
        if set(active) & set(freeze):
            raise PartitionError(f"{phase} active/freeze overlap.")
        if len(active) + len(freeze) != n_atoms:
            raise PartitionError(f"{phase} active/freeze does not cover all atoms.")
        if not original_set.issubset(set(freeze)):
            raise PartitionError(f"{phase} freeze does not preserve the user freeze mask.")
        if not effective_set.issubset(set(freeze)):
            raise PartitionError(f"{phase} freeze does not preserve frozen-layer atoms.")

    return MicroiterationPartition(
        n_atoms=n_atoms,
        original_freeze=frozen_t,
        macro_active_atoms=macro_active,
        macro_freeze_atoms=macro_freeze,
        micro_active_atoms=micro_active,
        micro_freeze_atoms=micro_freeze,
        ml_atoms=ml_t,
        link_parent_atoms=link_t,
        movable_mm_atoms=tuple(i for i in range(n_atoms) if i in micro_candidate),
    )


def _core_index_attr(core: Any, attr: str) -> List[int]:
    """Read an index list attribute from the accepted calculator core."""

    return list(getattr(core, attr, []) or [])


def resolve_partition_from_core(
    core: Any,
    n_atoms: int,
    original_freeze: Iterable[Any],
) -> MicroiterationPartition:
    """Resolve the partition strictly from the *already constructed* core (M44).

    Unlike the historical ``_collect_layer_atom_sets`` (which built a second
    calculator and swallowed every exception into four empty sets), this reads
    the layer indices off the accepted core and lets any failure surface as a
    :class:`PartitionError`.  A swallowed construction error can no longer
    masquerade as a valid empty-ML region.
    """

    if core is None:
        raise PartitionError("No calculator core available to resolve the partition.")
    try:
        ml = _core_index_attr(core, "ml_indices")
        hess_mm = _core_index_attr(core, "hess_mm_indices")
        movable_mm = _core_index_attr(core, "movable_mm_indices")
        frozen_mm = _core_index_attr(core, "frozen_layer_indices")
        link_parents = [
            int(mm_1) - 1 for (_ml_1, mm_1) in (getattr(core, "mlmm_links", []) or [])
        ]
    except PartitionError:
        raise
    except Exception as exc:  # noqa: BLE001 - re-raised as a loud PartitionError
        raise PartitionError(
            f"Failed to read the layer partition from the calculator core: {exc}"
        ) from exc
    return build_partition(
        n_atoms,
        ml=ml,
        link_parents=link_parents,
        hess_mm=hess_mm,
        movable_mm=movable_mm,
        frozen_mm=frozen_mm,
        original_freeze=original_freeze,
    )


def micro_reached_force_equilibrium(optimizer: Any) -> bool:
    """True when a stalled MM relaxation nonetheless met its force criteria.

    The micro stage exists to bring the MM subsystem to equilibrium before the
    macro step reads curvature, and equilibrium is a condition on forces. An
    energy plateau whose forces are already under their configured thresholds
    is that equilibrium; the step criteria that remain unmet describe how far
    the optimizer would still travel, which is not what the macro step needs
    from it. Reads only thresholds the optimizer already carries -- this
    introduces no new tolerance.
    """
    conv = getattr(optimizer, "convergence", None) or {}
    max_thresh = conv.get("max_force_thresh")
    if max_thresh is None:
        return False
    max_forces = getattr(optimizer, "max_forces", None) or []
    if not max_forces or float(max_forces[-1]) > float(max_thresh):
        return False
    rms_thresh = conv.get("rms_force_thresh")
    rms_forces = getattr(optimizer, "rms_forces", None) or []
    if rms_thresh is not None and rms_forces:
        return float(rms_forces[-1]) <= float(rms_thresh)
    return True


# ---------------------------------------------------------------------------
# One field-isomorphic optimizer outcome (C7) + nested micro outcome (C9)
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class OptimizerOutcome:
    """One optimizer's terminal state, field-isomorphic across every path.

    ``status`` uses the C7 vocabulary (``converged`` / ``not_converged`` /
    ``stalled``); ``converged`` is the fail-closed tri-state bit
    (``True``/``False``/``None``).  ``cycles`` counts *executed* optimizer
    cycles (``cur_cycle + 1``), never the configured budget.
    """

    status: str
    executed: bool
    converged: Optional[bool]
    cycles: Optional[int]
    max_cycles: Optional[int] = None
    stalled: bool = False
    stop_reason: Optional[str] = None

    @classmethod
    def from_optimizer(
        cls,
        optimizer: Any,
        *,
        max_cycles: Optional[int] = None,
        executed: bool = True,
    ) -> "OptimizerOutcome":
        """Build from a live pysisyphus optimizer (or product-local runner)."""

        status = optimizer_terminal_status(optimizer)
        converged = optimizer_converged_bit(optimizer)
        stalled = bool(getattr(optimizer, "is_stalled", False))
        if stalled:
            # A stall is explicitly not-converged (never None/unknown).
            converged = False
        return cls(
            status=status,
            executed=bool(executed),
            converged=converged,
            cycles=optimizer_cycle_count(optimizer),
            max_cycles=(int(max_cycles) if max_cycles is not None else None),
            stalled=stalled,
            stop_reason=(getattr(optimizer, "stop_reason", None) or None),
        )

    @classmethod
    def not_executed(
        cls, *, max_cycles: Optional[int] = None, reason: str = "not_executed"
    ) -> "OptimizerOutcome":
        """A leaf that never ran (fail-closed: unknown convergence, not usable)."""

        return cls(
            status="not_converged",
            executed=False,
            converged=None,
            cycles=None,
            max_cycles=(int(max_cycles) if max_cycles is not None else None),
            stalled=False,
            stop_reason=reason,
        )

    @classmethod
    def vacuous_success(cls, *, reason: str = "no_micro_active_dofs") -> "OptimizerOutcome":
        """A zero-cycle vacuous micro success (no movable MM coordinate to relax).

        Per the C9 blueprint a validated partition that leaves no micro-active
        DOF is *not* an error: it is a zero-cycle converged micro leaf.  No LBFGS
        is constructed with every atom frozen.
        """

        return cls(
            status="converged",
            executed=True,
            converged=True,
            cycles=0,
            max_cycles=None,
            stalled=False,
            stop_reason=None,
        )

    def to_dict(self) -> Dict[str, Any]:
        return {
            "status": self.status,
            "executed": bool(self.executed),
            "converged": self.converged,
            "cycles": self.cycles,
            "max_cycles": self.max_cycles,
            "stalled": bool(self.stalled),
            "stop_reason": self.stop_reason,
        }


@dataclass(frozen=True)
class MicroiterationOutcome:
    """The single terminal outcome of a macro/micro microiteration run (C9).

    ``aggregate`` is the one outcome a consumer reads; ``macro`` and the ordered
    ``micro_attempts`` retain both leaf truths so a failed micro is never
    flattened into a bare boolean.  ``macro_cycles`` counts executed macro
    evaluations; ``micro_cycles`` is the sum of executed micro optimizer cycles.
    """

    aggregate: OptimizerOutcome
    macro: OptimizerOutcome
    micro_attempts: Tuple[OptimizerOutcome, ...]
    macro_cycles: int
    micro_cycles: int
    partition: Optional[MicroiterationPartition] = None
    requested: bool = True
    used: bool = True
    fallback_reason: Optional[str] = None

    def to_result_object(self) -> Dict[str, Any]:
        """The additive ``microiteration`` object serialized into result.json."""

        obj: Dict[str, Any] = {
            "requested": bool(self.requested),
            "used": bool(self.used),
            "macro_cycles": int(self.macro_cycles),
            "micro_cycles": int(self.micro_cycles),
            "aggregate": self.aggregate.to_dict(),
            "macro": self.macro.to_dict(),
            "micro_attempts": [m.to_dict() for m in self.micro_attempts],
        }
        if self.fallback_reason:
            obj["fallback_reason"] = self.fallback_reason
        if self.partition is not None:
            obj["partition"] = self.partition.provenance()
        return obj


def build_aggregate(
    macro: OptimizerOutcome,
    micro_attempts: Sequence[OptimizerOutcome],
    *,
    max_cycles: Optional[int] = None,
) -> OptimizerOutcome:
    """Fold macro + micro leaves into one fail-closed aggregate outcome (M46).

    Aggregate is ``converged`` only when the macro is explicitly converged AND
    the latest *required* micro relaxation is explicitly converged.  A macro or
    latest-micro stall yields ``stalled`` (carrying its reason); a macro that
    ran out of cycles, or a latest micro that did not explicitly converge
    (``False`` / unknown / max-cycle), yields ``not_converged`` with the exact
    child reason.  A normal Python return from ``run()`` is never, by itself,
    treated as convergence.
    """

    latest = micro_attempts[-1] if micro_attempts else None

    if macro.stalled:
        return OptimizerOutcome(
            status="stalled",
            executed=True,
            converged=False,
            cycles=macro.cycles,
            max_cycles=max_cycles,
            stalled=True,
            stop_reason=macro.stop_reason or "macro energy plateau (not converged)",
        )
    if latest is not None and latest.stalled:
        return OptimizerOutcome(
            status="stalled",
            executed=True,
            converged=False,
            cycles=macro.cycles,
            max_cycles=max_cycles,
            stalled=True,
            stop_reason=latest.stop_reason
            or "energy plateau in the latest micro (MM) relaxation",
        )

    macro_conv = macro.converged is True
    latest_micro_conv = True if latest is None else (latest.converged is True)

    if macro_conv and latest_micro_conv:
        return OptimizerOutcome(
            status="converged",
            executed=True,
            converged=True,
            cycles=macro.cycles,
            max_cycles=max_cycles,
            stalled=False,
            stop_reason=None,
        )

    if not macro_conv:
        reason = macro.stop_reason or "macro_not_converged"
    else:
        # macro converged but the latest micro did not.
        if latest is not None and latest.converged is False:
            reason = latest.stop_reason or "micro_max_cycles"
        else:
            reason = (latest.stop_reason if latest is not None else None) or "micro_convergence_unknown"
    return OptimizerOutcome(
        status="not_converged",
        executed=True,
        converged=False,
        cycles=macro.cycles,
        max_cycles=max_cycles,
        stalled=False,
        stop_reason=reason,
    )


# ---------------------------------------------------------------------------
# IRC / Hessian device policy (M43)
# ---------------------------------------------------------------------------


def resolve_hessian_device(requested: str, cuda_available: bool) -> Tuple[str, str]:
    """Calibrated GPU-first device policy for the IRC integration Hessian (M43).

    * ``auto``  -> GPU-first: ``cuda`` when available, else ``cpu`` (logged).
    * ``cuda``  -> stays ``cuda``; raises :class:`ValueError` when CUDA is
      unavailable.  An explicit CUDA request is never silently moved to CPU.
    * ``cpu``   -> ``cpu`` (an explicit, calibrated offload choice).

    Returns ``(effective_device, reason)``.  Large tensors stay GPU-resident
    under the default (auto) policy; the CPU fallback is an explicit, logged
    user choice, not a silent inefficiency response.
    """

    req = (requested or "auto").strip().lower()
    if req == "cpu":
        return "cpu", "explicit_cpu"
    if req == "cuda":
        if not cuda_available:
            raise ValueError(
                "--hess-device cuda was requested but no CUDA device is available; "
                "an explicit CUDA request is never silently moved to CPU. "
                "Use --hess-device cpu (or auto) instead."
            )
        return "cuda", "explicit_cuda"
    if req == "auto":
        if cuda_available:
            return "cuda", "auto_gpu_first"
        return "cpu", "auto_no_cuda"
    raise ValueError(f"Unknown --hess-device value {requested!r}.")


__all__ = [
    "PartitionError",
    "MicroiterationPartition",
    "build_partition",
    "resolve_partition_from_core",
    "OptimizerOutcome",
    "MicroiterationOutcome",
    "build_aggregate",
    "resolve_hessian_device",
]
