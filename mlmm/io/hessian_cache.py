"""In-process Hessian cache for the ``all`` workflow.

Stores Hessian matrices as CPU tensors (or numpy arrays) in a module-level
dict so that subsequent CLI stages executed via ``_run_cli_main()`` can
reuse them without recomputation.

Each entry may carry ``active_dofs`` — a list of DOF indices that the stored
(partial) Hessian spans.  Consumers set ``geometry.within_partial_hessian``
before assigning partial Hessians to ``geometry.cart_hessian``.

Ownership and reuse semantics
-----------------------------
* **Defensive ownership.** ``store()`` copies supported array and container
  payloads in,
  and every read (``load`` / ``load_matching``) returns a *fresh* defensive
  snapshot.  A consumer may mass-weight, project, or otherwise mutate what it
  receives without corrupting the raw Cartesian artifact retained in the cache.
* **Reuse identity.** A cached Hessian is only reused when the captured run,
  system, evaluator, potential, active-space, constraint, and artifact fields
  match (coordinates use the established bohr round-trip tolerance).
  ``build_identity()`` records those fields as a
  private ``hessian-cache-identity/v1`` token and ``load_matching()`` is the
  single reuse chokepoint that requires an exact match.  A legacy
  coordinate-only entry (no identity) is never reused through
  ``load_matching``.

This is an implementation detail of mlmm_toolkit v0.3.3.
"""

import hashlib
import os
from collections.abc import Mapping
from typing import Any, Dict, Optional, Sequence

import numpy as np
import torch

from mlmm.backends.methods import normalize_mm_hessian_mode
from mlmm.core.result_commit import MLMM_RUN_ID_ENV as RUN_ID_ENV

_cache: Dict[str, Any] = {}

# Internal cache-token schema.
IDENTITY_SCHEMA = "hessian-cache-identity/v1"

# Role -> artifact method/source contract.  Producers stamp this; consumers
# expect it for the same role, so both sides derive the same identity.
_ROLE_SOURCE: Dict[str, str] = {
    "ts": "tsopt_exact",
    "irc_left": "irc_endpoint_quasi_newton",
    "irc_right": "irc_endpoint_quasi_newton",
    "irc_endpoint": "irc_endpoint_quasi_newton",
}

# MLMM stores the ML model name and inference precision under BACKEND-PREFIXED
# config keys, never the generic ``model`` / ``precision`` keys (see
# ``mlmm.core.defaults.MLMM_CALC_KW``).  These maps let the cache identity
# select the effective model + precision by ``backend`` so that two runs with a
# different ML model or a different fp32-vs-fp64 precision produce DIFFERENT
# identities.  Mirrors the canonical resolver in ``mlmm.core.utils``
# (``calculator_provenance``); kept local so this isolated cache module carries
# no extra import.
_BACKEND_MODEL_KEY: Dict[str, str] = {
    "uma": "uma_model",
    "orb": "orb_model",
    "mace": "mace_model",
    "aimnet2": "aimnet2_model",
}
_BACKEND_PRECISION_KEY: Dict[str, str] = {
    "uma": "uma_precision",
    "orb": "orb_precision",
    "mace": "mace_dtype",
    # aimnet2 exposes no precision toggle (always fp32); an explicit
    # ``aimnet2_precision`` is honored if a config ever sets one.
    "aimnet2": "aimnet2_precision",
}


# ---------------------------------------------------------------------------
# Defensive copy / canonicalization helpers
# ---------------------------------------------------------------------------
def _clone_value(value: Any) -> Any:
    """Recursively copy tensors / arrays / containers; leave scalars shared."""

    if isinstance(value, torch.Tensor):
        return value.detach().clone()
    if isinstance(value, np.ndarray):
        return np.array(value, copy=True)
    if isinstance(value, Mapping):
        return {key: _clone_value(val) for key, val in value.items()}
    if isinstance(value, list):
        return [_clone_value(val) for val in value]
    if isinstance(value, tuple):
        return tuple(_clone_value(val) for val in value)
    return value


def _canon(value: Any) -> Any:
    """Normalize to comparable Python primitives (arrays -> nested lists)."""

    if isinstance(value, torch.Tensor):
        value = value.detach().cpu().numpy()
    if isinstance(value, np.ndarray):
        return [_canon(item) for item in value.tolist()]
    if isinstance(value, Mapping):
        return {str(key): _canon(val) for key, val in value.items()}
    if isinstance(value, (list, tuple)):
        return [_canon(item) for item in value]
    if isinstance(value, np.integer):
        return int(value)
    if isinstance(value, np.floating):
        return float(value)
    if isinstance(value, (bool, int, float, str)) or value is None:
        return value
    return str(value)


def _as_list(value: Any) -> list:
    if value is None:
        return []
    if isinstance(value, torch.Tensor):
        value = value.detach().cpu().numpy()
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, (list, tuple, set, frozenset)):
        return list(value)
    return [value]


def _int_list(value: Any) -> list:
    return [int(item) for item in _as_list(value)]


def _atom_token(value: Any) -> Any:
    try:
        return int(value)
    except (TypeError, ValueError):
        return str(value)


def _norm_str(value: Any) -> Optional[str]:
    if value is None:
        return None
    return str(value)


def _norm_int(value: Any) -> Optional[int]:
    if value is None:
        return None
    try:
        return int(value)
    except (TypeError, ValueError):
        return None


def _file_digest(path: Any) -> Optional[str]:
    """Return the SHA-256 of a file's content, or *None* when unreadable."""

    try:
        hasher = hashlib.sha256()
        with open(os.fspath(path), "rb") as handle:
            for block in iter(lambda: handle.read(1 << 20), b""):
                hasher.update(block)
        return hasher.hexdigest()
    except (OSError, TypeError):
        return None


def _effective_model_precision(calc_cfg: Mapping) -> tuple:
    """Resolve the effective ML model + precision from the backend-specific keys.

    The identity must differ when the ML model or the inference precision
    differs (fp32 vs fp64 changes the finite-difference Hessian and its
    imaginary-mode count), so the value is selected by ``calc_cfg["backend"]``:
    ``uma`` -> ``uma_model`` / ``uma_precision``; ``orb`` -> ``orb_model`` /
    ``orb_precision``; ``mace`` -> ``mace_model`` / ``mace_dtype``; ``aimnet2``
    -> ``aimnet2_model`` (+ any ``aimnet2_precision``).  A custom ASE
    calculator (``backend='custom'``) carries no MLIP model variant, so its
    ``calc_file`` / ``calc_factory`` reference stands in as the model identity.
    Any absent key falls back to *None*.
    """

    backend = calc_cfg.get("backend")
    backend = str(backend).strip().lower() if backend is not None else None
    if backend == "custom":
        calc_file = calc_cfg.get("calc_file")
        factory = calc_cfg.get("calc_factory") or "get_calculator"
        if calc_file:
            model = f"{calc_file}:{factory}"
        else:
            model = str(factory)
        return model, None
    model_key = _BACKEND_MODEL_KEY.get(backend) if backend is not None else None
    precision_key = (
        _BACKEND_PRECISION_KEY.get(backend) if backend is not None else None
    )
    model = calc_cfg.get(model_key) if model_key is not None else None
    precision = calc_cfg.get(precision_key) if precision_key is not None else None
    return model, precision


def _potential_identity(calc_cfg: Mapping) -> Dict[str, Any]:
    """Composition of the evaluated ML/MM potential beyond the bare backend.

    Captures the QM/MM composition fields whose change makes a stored Hessian an
    invalid substitute: MM backend, link method, embedding, CMAP, the QM-region
    charge/spin, the region/layer definition, and the *content* of the topology
    / region files (so replacing a parm7's bytes at the same path rejects
    reuse. The topology-content identity is hashed here without importing
    hessian_ff).
    """

    potential: Dict[str, Any] = {}
    # Scalar QM/MM composition knobs.
    for key in (
        "mm_backend",
        "link_atom_method",
        "embedcharge",
        "embedcharge_cutoff",
        "use_cmap",
        "model_charge",
        "model_mult",
        "movable_cutoff",
        "use_bfactor_layers",
        "hessian_calc_mode",
        "H_double",
        "hess_cutoff",
        "return_partial_hessian",
        "symmetrize_hessian",
    ):
        val = calc_cfg.get(key)
        if val is None:
            continue
        if isinstance(val, str) and val.strip().lower() in ("", "none"):
            continue
        if key in {"hess_cutoff", "movable_cutoff"}:
            try:
                if np.isposinf(float(val)):
                    potential[key] = "unbounded"
                    continue
            except (TypeError, ValueError):
                pass
        potential[key] = _canon(val)
    backend = str(calc_cfg.get("backend") or "").strip().lower()
    if backend == "uma" and calc_cfg.get("uma_task_name") is not None:
        potential["uma_task_name"] = str(calc_cfg["uma_task_name"])
    model, _precision = _effective_model_precision(calc_cfg)
    model_digest = _file_digest(model)
    if model_digest is not None:
        potential["model_sha256"] = model_digest
    effective_mm_mode = normalize_mm_hessian_mode(
        calc_cfg.get("mm_hessian_mode"),
        mm_fd=bool(calc_cfg.get("mm_fd", True)),
    )
    potential["mm_hessian_mode"] = effective_mm_mode
    if effective_mm_mode == "finite_difference":
        potential["mm_fd_delta"] = float(calc_cfg.get("mm_fd_delta", 1.0e-3))
    # Topology / region source files — identity is the content digest so a
    # same-path byte replacement rejects reuse (topology content identity).
    for key in ("real_parm7", "model_pdb", "input_pdb"):
        path = calc_cfg.get(key)
        if path:
            digest = _file_digest(path)
            if digest is None:
                # Preserve a fail-closed identity when the file cannot be read.
                potential[key] = str(path)
            else:
                potential[f"{key}_sha256"] = digest
    calc_file = calc_cfg.get("calc_file")
    if calc_file:
        potential["calc_file"] = str(calc_file)
        digest = _file_digest(calc_file)
        if digest is not None:
            potential["calc_file_sha256"] = digest
    # Explicit ML/MM Hessian region atom list when threaded through the config.
    hess_mm = calc_cfg.get("hess_mm_atoms")
    if hess_mm is not None:
        potential["hess_mm_atoms"] = _int_list(hess_mm)
    # Explicit MLMM region partition + link map.  These are NOT subsumed by the
    # captured freeze set: the workflow only folds the B-factor
    # ``frozen_indices`` into ``calc_cfg["freeze_atoms"]`` (via
    # ``apply_layer_freeze_constraints``); an explicit ``frozen_mm_atoms`` /
    # ``movable_mm_atoms`` is consumed only inside the calculator's
    # ``_compute_layer_indices``, and ``link_mlmm`` is never reflected in the
    # freeze set at all.  Capture them (sorted int lists / a stable link repr)
    # so a region or link-boundary change at matching coordinates rejects reuse.
    for atoms_key in ("movable_mm_atoms", "frozen_mm_atoms"):
        atoms_val = calc_cfg.get(atoms_key)
        if atoms_val is not None:
            potential[atoms_key] = sorted(_int_list(atoms_val))
    link = calc_cfg.get("link_mlmm")
    if link is not None:
        canon_link = _canon(link)
        if isinstance(canon_link, list):
            try:
                canon_link = sorted(canon_link)
            except TypeError:
                pass
        potential["link_mlmm"] = canon_link
    return potential


# ---------------------------------------------------------------------------
# Identity construction and comparison
# ---------------------------------------------------------------------------
def build_identity(
    *,
    atoms,
    cart_coords,
    run_id: Optional[str] = None,
    backend=None,
    model=None,
    precision=None,
    charge=None,
    spin=None,
    potential: Optional[Mapping] = None,
    active_atoms=None,
    active_dofs=None,
    constraints: Optional[Mapping] = None,
    units: str = "bohr",
    representation: str = "raw_cartesian",
    method=None,
    source=None,
) -> Dict[str, Any]:
    """Assemble a private ``hessian-cache-identity/v1`` token.

    ``run_id`` defaults to the current in-process invocation ID.  Every field
    except ``cart_coords`` participates in an exact match; the coordinates use
    the established bohr round-trip tolerance.
    """

    if run_id is None:
        run_id = os.environ.get(RUN_ID_ENV)
    coords = cart_coords
    if isinstance(coords, torch.Tensor):
        coords = coords.detach().cpu().numpy()
    coords = np.asarray(coords, dtype=float).reshape(-1).copy()
    return {
        "schema": IDENTITY_SCHEMA,
        "run_id": run_id,
        "system": {
            "atoms": [_atom_token(a) for a in _as_list(atoms)],
            "cart_coords_bohr": coords,
        },
        "evaluator": {
            "backend": _norm_str(backend),
            "model": _norm_str(model),
            "precision": _norm_str(precision),
            "charge": _norm_int(charge),
            "spin": _norm_int(spin),
            "potential": _canon(potential or {}),
        },
        "active_space": {
            "active_atoms": _int_list(active_atoms),
            "active_dofs": _int_list(active_dofs),
            "constraints": _canon(constraints or {}),
        },
        "artifact": {
            "units": str(units),
            "representation": str(representation),
            "method": _norm_str(method),
            "source": _norm_str(source),
        },
    }


def identity_from_context(
    geometry,
    calc_cfg: Optional[Mapping],
    *,
    role: Optional[str] = None,
    source: Optional[str] = None,
    cart_coords=None,
) -> Dict[str, Any]:
    """Build an identity from a geometry + resolved evaluator config.

    Producers and consumers of a given role call this with the same inputs so
    that a compatible artifact produces a matching token.  The active space is
    derived deterministically from the freeze set so that both sides agree
    without threading partial-Hessian bookkeeping through the identity.
    """

    calc_cfg = dict(calc_cfg or {})
    if cart_coords is None:
        cart_coords = geometry.cart_coords
    atoms = getattr(geometry, "atomic_numbers", None)
    if atoms is None:
        atoms = getattr(geometry, "atoms", [])
    atoms_list = _as_list(atoms)
    if source is None and role is not None:
        source = _ROLE_SOURCE.get(role)

    freeze = calc_cfg.get("freeze_atoms")
    if freeze is None:
        freeze = getattr(geometry, "freeze_atoms", None)
    frozen = set(_int_list(freeze))
    n_atoms = len(atoms_list)
    derived_active_atoms = [i for i in range(n_atoms) if i not in frozen]
    derived_active_dofs = [
        3 * a + k for a in derived_active_atoms for k in range(3)
    ]

    # MLMM never sets the generic ``model`` / ``precision`` keys; resolve the
    # effective ML model + precision from the backend-specific keys so the
    # identity distinguishes two ML models and fp32 vs fp64.
    model, precision = _effective_model_precision(calc_cfg)

    return build_identity(
        atoms=atoms_list,
        cart_coords=cart_coords,
        backend=calc_cfg.get("backend"),
        model=model,
        precision=precision,
        charge=calc_cfg.get("charge", calc_cfg.get("model_charge")),
        spin=calc_cfg.get("spin", calc_cfg.get("model_mult")),
        potential=_potential_identity(calc_cfg),
        active_atoms=derived_active_atoms,
        active_dofs=derived_active_dofs,
        constraints={"freeze_atoms": sorted(frozen)},
        source=source,
        method=source,
    )


def persistent_identity_from_context(
    geometry,
    calc_cfg: Optional[Mapping],
    *,
    method: str = "mlmm_exact_hessian",
) -> Dict[str, Any]:
    """Return the cross-process PES identity, excluding run/coordinate fields."""

    identity = identity_from_context(
        geometry,
        calc_cfg,
        source=method,
    )
    persistent = _canon(identity)
    persistent.pop("run_id", None)
    system = persistent.get("system")
    if isinstance(system, Mapping):
        system = dict(system)
        system.pop("cart_coords_bohr", None)
        persistent["system"] = system
    return persistent


def _identity_coords(identity: Mapping) -> Optional[np.ndarray]:
    system = identity.get("system")
    if not isinstance(system, Mapping):
        return None
    coords = system.get("cart_coords_bohr")
    if coords is None:
        return None
    if isinstance(coords, torch.Tensor):
        coords = coords.detach().cpu().numpy()
    return np.asarray(coords, dtype=float).reshape(-1)


def _identity_without_coords(identity: Mapping) -> Dict[str, Any]:
    canon = _canon(identity)
    system = canon.get("system")
    if isinstance(system, Mapping):
        system = dict(system)
        system.pop("cart_coords_bohr", None)
        canon["system"] = system
    return canon


def identities_match(
    stored: Any,
    expected: Any,
    *,
    atol: float = 1.0e-5,
) -> bool:
    """Whether two identity tokens describe an interchangeable artifact."""

    if not isinstance(stored, Mapping) or not isinstance(expected, Mapping):
        return False
    if stored.get("schema") != IDENTITY_SCHEMA:
        return False
    if expected.get("schema") != IDENTITY_SCHEMA:
        return False
    stored_run = stored.get("run_id")
    expected_run = expected.get("run_id")
    if not stored_run or not expected_run or stored_run != expected_run:
        return False

    stored_coords = _identity_coords(stored)
    expected_coords = _identity_coords(expected)
    if stored_coords is None or expected_coords is None:
        return False
    if stored_coords.shape != expected_coords.shape:
        return False
    if not bool(
        np.allclose(stored_coords, expected_coords, rtol=0.0, atol=float(atol))
    ):
        return False

    return _identity_without_coords(stored) == _identity_without_coords(expected)


def reconcile_active_hessian(
    entry: Mapping[str, Any],
    required_dofs: Sequence[int],
    *,
    full_n_dof: int,
) -> Optional[torch.Tensor]:
    """Return a cached Hessian in exactly ``required_dofs`` order.

    Compact matrices without ordered DOF metadata are ambiguous and fail
    closed.  A full ``3N`` matrix may omit metadata because its row/column
    order is the canonical Cartesian order.
    """

    raw = entry.get("hessian")
    if raw is None:
        return None
    hessian = (
        raw.detach().clone()
        if isinstance(raw, torch.Tensor)
        else torch.as_tensor(raw, dtype=torch.float64).clone()
    )
    if hessian.ndim != 2 or hessian.shape[0] != hessian.shape[1]:
        return None

    required = [int(dof) for dof in required_dofs]
    if len(set(required)) != len(required):
        return None
    if any(dof < 0 or dof >= int(full_n_dof) for dof in required):
        return None

    active_raw = entry.get("active_dofs")
    if active_raw is None:
        if hessian.shape[0] != int(full_n_dof):
            return None
        active = list(range(int(full_n_dof)))
    else:
        try:
            active = [int(dof) for dof in active_raw]
        except (TypeError, ValueError):
            return None
        if len(active) != hessian.shape[0] or len(set(active)) != len(active):
            return None
        if any(dof < 0 or dof >= int(full_n_dof) for dof in active):
            return None

    positions = {dof: index for index, dof in enumerate(active)}
    if any(dof not in positions for dof in required):
        return None
    local = torch.as_tensor(
        [positions[dof] for dof in required],
        dtype=torch.long,
        device=hessian.device,
    )
    from pysisyphus._array import active_square

    return active_square(hessian, local)


# ---------------------------------------------------------------------------
# Store / load
# ---------------------------------------------------------------------------
def _store_identity(identity: Mapping) -> Dict[str, Any]:
    ident = _clone_value(dict(identity))
    system = ident.get("system")
    if isinstance(system, Mapping):
        system = dict(system)
        coords = system.get("cart_coords_bohr")
        if coords is not None:
            if isinstance(coords, torch.Tensor):
                coords = coords.detach().cpu().numpy()
            system["cart_coords_bohr"] = (
                np.asarray(coords, dtype=float).reshape(-1).copy()
            )
        ident["system"] = system
    return ident


def store(
    key: str,
    H,
    active_dofs: Optional[Sequence[int]] = None,
    meta: Optional[dict] = None,
    identity: Optional[Mapping] = None,
) -> None:
    """Cache a (possibly partial) Hessian on CPU with an optional identity."""

    if isinstance(H, torch.Tensor):
        h_cpu = H.detach().cpu().clone()
    else:
        h_cpu = np.array(H, copy=True)
    meta_copy = _clone_value(dict(meta or {}))
    if "cart_coords" in meta_copy:
        coords = meta_copy["cart_coords"]
        if isinstance(coords, torch.Tensor):
            coords = coords.detach().cpu().numpy()
        meta_copy["cart_coords"] = np.asarray(coords, dtype=float).copy()
    entry: Dict[str, Any] = {
        "hessian": h_cpu,
        "active_dofs": list(active_dofs) if active_dofs is not None else None,
        "meta": meta_copy,
    }
    if identity is not None:
        entry["identity"] = _store_identity(identity)
    _cache[key] = entry


def _snapshot(entry: Mapping) -> Dict[str, Any]:
    """Return a defensive copy of a cache entry."""

    snap: Dict[str, Any] = {
        "hessian": _clone_value(entry.get("hessian")),
        "active_dofs": (
            list(entry["active_dofs"])
            if entry.get("active_dofs") is not None
            else None
        ),
        "meta": _clone_value(entry.get("meta") or {}),
    }
    if entry.get("identity") is not None:
        snap["identity"] = _clone_value(entry["identity"])
    return snap


def load(key: str) -> Optional[Dict[str, Any]]:
    """Return a defensive snapshot of a cached entry, or *None*.

    The returned mapping is safe to mutate; the retained cache artifact stays
    raw, Cartesian, unweighted, and unprojected.
    """

    entry = _cache.get(key)
    if entry is None:
        return None
    return _snapshot(entry)


def load_matching(
    role: str,
    expected: Mapping,
    *,
    atol: float = 1.0e-5,
) -> Optional[Dict[str, Any]]:
    """Single reuse chokepoint: return a snapshot only on a full identity match.

    A cached entry without an identity token (legacy coordinate-only) is never
    reused here.  Missing/mismatched run ID, evaluator, system, active space,
    constraints, potential, or representation all reject.
    """

    entry = _cache.get(role)
    if entry is None:
        return None
    stored_identity = entry.get("identity")
    if stored_identity is None:
        return None
    if not identities_match(stored_identity, expected, atol=atol):
        return None
    return _snapshot(entry)


def matches_cart_coords(
    entry: Dict[str, Any],
    cart_coords,
    *,
    atol: float = 1.0e-5,
) -> bool:
    """Whether an entry belongs to ``cart_coords``.

    Entries without a coordinate fingerprint are not safe for cross-stage
    reuse. Coordinates are in bohr and the tolerance accommodates serialized
    structure round-tripping.  Coordinate identity alone is *not* sufficient
    for reuse — see ``load_matching`` for the full-identity chokepoint.
    """
    cached_coords = entry.get("meta", {}).get("cart_coords")
    if cached_coords is None:
        return False
    if isinstance(cart_coords, torch.Tensor):
        cart_coords = cart_coords.detach().cpu().numpy()
    cached_arr = np.asarray(cached_coords, dtype=float).reshape(-1)
    current_arr = np.asarray(cart_coords, dtype=float).reshape(-1)
    return cached_arr.shape == current_arr.shape and bool(
        np.allclose(cached_arr, current_arr, rtol=0.0, atol=float(atol))
    )


def discard(key: str) -> None:
    """Drop one cached entry if present."""
    _cache.pop(key, None)


def clear() -> None:
    """Drop all cached Hessians."""
    _cache.clear()
