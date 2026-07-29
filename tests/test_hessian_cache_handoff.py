"""Regression tests for geometry-safe in-process Hessian handoff."""

from __future__ import annotations

import numpy as np
import pytest
import torch

from mlmm.io import hessian_cache


def setup_function() -> None:
    hessian_cache.clear()


def teardown_function() -> None:
    hessian_cache.clear()


def test_endpoint_cache_matches_xyz_round_trip_coordinates() -> None:
    coords = np.array([0.0, 1.234567890, -2.345678901])
    hessian_cache.store(
        "irc_endpoint",
        np.eye(3),
        meta={"cart_coords": coords, "irc_direction": "forward"},
    )
    entry = hessian_cache.load("irc_endpoint")

    # XYZ serialization changes coordinates slightly; the endpoint cache is
    # still valid within the explicit bohr tolerance.
    round_tripped = coords + np.array([0.0, 2.0e-6, -2.0e-6])
    assert entry is not None
    assert hessian_cache.matches_cart_coords(entry, round_tripped)


def test_endpoint_cache_rejects_swapped_or_stale_geometry() -> None:
    hessian_cache.store(
        "irc_endpoint",
        np.eye(3),
        meta={"cart_coords": np.array([0.0, 0.0, 0.0])},
    )
    entry = hessian_cache.load("irc_endpoint")

    assert entry is not None
    assert not hessian_cache.matches_cart_coords(
        entry, np.array([0.0, 0.0, 1.0e-3])
    )
    assert not hessian_cache.matches_cart_coords(entry, np.zeros(6))


def test_cache_without_coordinate_identity_is_not_reused() -> None:
    hessian_cache.store("ts", np.eye(3))
    entry = hessian_cache.load("ts")

    assert entry is not None
    assert not hessian_cache.matches_cart_coords(entry, np.zeros(3))


def test_store_copies_tensor_and_coordinate_metadata() -> None:
    H = torch.eye(3, dtype=torch.float64)
    coords = torch.tensor([1.0, 2.0, 3.0], dtype=torch.float64)
    hessian_cache.store(
        "irc_endpoint",
        H,
        active_dofs=[0, 1, 2],
        meta={"cart_coords": coords},
    )
    H.add_(5.0)
    coords.add_(5.0)

    entry = hessian_cache.load("irc_endpoint")
    assert entry is not None
    torch.testing.assert_close(entry["hessian"], torch.eye(3, dtype=torch.float64))
    np.testing.assert_allclose(entry["meta"]["cart_coords"], [1.0, 2.0, 3.0])


def test_discard_prevents_missing_endpoint_from_reusing_previous_seed() -> None:
    hessian_cache.store("irc_endpoint", np.eye(3))
    hessian_cache.discard("irc_endpoint")
    assert hessian_cache.load("irc_endpoint") is None


# ---------------------------------------------------------------------------
# — defensive in-process ownership
# ---------------------------------------------------------------------------
def test_load_returns_independent_tensor_snapshots() -> None:
    hessian_cache.store(
        "ts",
        torch.eye(3, dtype=torch.float64),
        meta={"cart_coords": np.zeros(3)},
    )
    first = hessian_cache.load("ts")
    second = hessian_cache.load("ts")

    assert first is not None and second is not None
    assert first["hessian"].data_ptr() != second["hessian"].data_ptr()

    # Mutating a loaded snapshot must not corrupt the retained raw artifact.
    first["hessian"].add_(5.0)
    again = hessian_cache.load("ts")
    torch.testing.assert_close(again["hessian"], torch.eye(3, dtype=torch.float64))


def test_load_snapshot_isolates_numpy_meta_and_active_dofs() -> None:
    hessian_cache.store(
        "ts",
        np.eye(3),
        active_dofs=[0, 1, 2],
        meta={"cart_coords": np.zeros(3), "source": "tsopt_exact"},
    )
    snap = hessian_cache.load("ts")
    snap["meta"]["cart_coords"][0] = 99.0
    snap["active_dofs"].append(999)
    # Mutating a loaded numpy Hessian through torch.as_tensor must not persist.
    torch.as_tensor(snap["hessian"]).mul_(0.0)

    fresh = hessian_cache.load("ts")
    assert fresh["meta"]["cart_coords"][0] == 0.0
    assert fresh["active_dofs"] == [0, 1, 2]
    np.testing.assert_array_equal(fresh["hessian"], np.eye(3))


# ---------------------------------------------------------------------------
# — complete reuse identity
# ---------------------------------------------------------------------------
def _identity(
    *,
    run="run-A",
    backend="uma",
    model="m",
    precision="fp64",
    charge=0,
    spin=1,
    active_atoms=(0, 1),
    active_dofs=(0, 1, 2, 3, 4, 5),
    potential=None,
    constraints=None,
    source="tsopt_exact",
    coords=None,
):
    if coords is None:
        coords = np.zeros(6)
    return hessian_cache.build_identity(
        atoms=[1, 1],
        cart_coords=coords,
        run_id=run,
        backend=backend,
        model=model,
        precision=precision,
        charge=charge,
        spin=spin,
        potential=potential or {},
        active_atoms=list(active_atoms),
        active_dofs=list(active_dofs),
        constraints=constraints or {"freeze_atoms": []},
        source=source,
        method=source,
    )


def test_load_matching_accepts_only_full_identity() -> None:
    hessian_cache.store("ts", np.eye(6), identity=_identity())

    # Exact identity, exact coordinates.
    assert hessian_cache.load_matching("ts", _identity()) is not None
    # Coordinates within the bohr round-trip tolerance.
    assert hessian_cache.load_matching("ts", _identity(coords=np.full(6, 2.0e-6))) is not None
    # Every single-field change rejects.
    assert hessian_cache.load_matching("ts", _identity(coords=np.full(6, 1.0e-3))) is None
    assert hessian_cache.load_matching("ts", _identity(backend="orb")) is None
    assert hessian_cache.load_matching("ts", _identity(model="other")) is None
    assert hessian_cache.load_matching("ts", _identity(precision="fp32")) is None
    assert hessian_cache.load_matching("ts", _identity(charge=-1)) is None
    assert hessian_cache.load_matching("ts", _identity(spin=3)) is None
    assert hessian_cache.load_matching("ts", _identity(active_dofs=(0, 1, 2))) is None
    assert hessian_cache.load_matching("ts", _identity(active_atoms=(0,))) is None
    assert hessian_cache.load_matching(
        "ts", _identity(constraints={"freeze_atoms": [0]})
    ) is None
    assert hessian_cache.load_matching(
        "ts", _identity(potential={"mm_backend": "openmm"})
    ) is None
    assert hessian_cache.load_matching("ts", _identity(source="irc_endpoint_quasi_newton")) is None
    assert hessian_cache.load_matching("ts", _identity(run="run-B")) is None


def test_load_matching_returns_defensive_snapshot() -> None:
    hessian_cache.store("ts", np.eye(6), active_dofs=[0, 1, 2, 3, 4, 5], identity=_identity())
    snap = hessian_cache.load_matching("ts", _identity())
    assert snap is not None
    snap["hessian"][0, 0] = 42.0
    snap["active_dofs"].append(7)
    again = hessian_cache.load_matching("ts", _identity())
    np.testing.assert_array_equal(again["hessian"], np.eye(6))
    assert again["active_dofs"] == [0, 1, 2, 3, 4, 5]


def test_legacy_coordinate_only_entry_is_never_reused_by_load_matching() -> None:
    # A legacy entry carries no identity token.
    hessian_cache.store("ts", np.eye(6), meta={"cart_coords": np.zeros(6)})
    assert hessian_cache.load_matching("ts", _identity()) is None
    # The coordinate-only load path still returns a snapshot for legacy callers.
    assert hessian_cache.load("ts") is not None


def test_missing_run_id_never_reuses() -> None:
    hessian_cache.store("ts", np.eye(6), identity=_identity(run=None))
    assert hessian_cache.load_matching("ts", _identity(run=None)) is None
    assert hessian_cache.load_matching("ts", _identity(run="run-A")) is None


def test_identity_from_context_round_trips_through_cache(monkeypatch) -> None:
    from mlmm.core.result_commit import MLMM_RUN_ID_ENV as RUN_ID_ENV

    monkeypatch.setenv(RUN_ID_ENV, "run-X")

    class _Geom:
        atomic_numbers = np.array([1, 1, 8])
        cart_coords = np.arange(9, dtype=float)
        freeze_atoms = np.array([0])

    # Realistic mlmm calc_cfg: the ML model + precision live under the
    # backend-prefixed keys (uma_model / uma_precision), never the generic
    # model / precision keys.
    calc_cfg = {
        "backend": "uma",
        "uma_model": "uma-s-1p2",
        "uma_precision": "fp64",
        "charge": 0,
        "spin": 1,
        "freeze_atoms": [0],
    }
    hessian_cache.store(
        "ts",
        np.eye(6),
        identity=hessian_cache.identity_from_context(_Geom(), calc_cfg, role="ts"),
    )
    # Same evaluator context reuses.
    assert hessian_cache.load_matching(
        "ts", hessian_cache.identity_from_context(_Geom(), calc_cfg, role="ts")
    ) is not None
    # A different evaluator does not.
    other = dict(calc_cfg, backend="orb")
    assert hessian_cache.load_matching(
        "ts", hessian_cache.identity_from_context(_Geom(), other, role="ts")
    ) is None


@pytest.mark.parametrize(
    "change",
    [
        {"hessian_calc_mode": "Analytical"},
        {"mm_hessian_mode": "none"},
        {"mm_fd_delta": 2.0e-3},
        {"symmetrize_hessian": False},
        {"H_double": False},
        {"return_partial_hessian": False},
    ],
)
def test_identity_rejects_a_different_hessian_construction_policy(
    change,
    monkeypatch,
) -> None:
    from mlmm.core.result_commit import MLMM_RUN_ID_ENV as RUN_ID_ENV

    monkeypatch.setenv(RUN_ID_ENV, "run-hessian-policy")

    class _Geom:
        atoms = ["H", "H"]
        atomic_numbers = [1, 1]
        cart_coords = np.zeros(6)
        freeze_atoms = []

    base = {
        "backend": "uma",
        "uma_model": "uma-s-1p2",
        "uma_precision": "fp32",
        "model_charge": 0,
        "model_mult": 1,
        "hessian_calc_mode": "FiniteDifference",
        "mm_hessian_mode": "finite_difference",
        "mm_fd_delta": 1.0e-3,
        "symmetrize_hessian": True,
        "H_double": True,
        "return_partial_hessian": True,
    }

    base_identity = hessian_cache.identity_from_context(_Geom(), base, role="ts")
    other_identity = hessian_cache.identity_from_context(
        _Geom(), {**base, **change}, role="ts"
    )
    assert (
        base_identity["evaluator"]["potential"]
        != other_identity["evaluator"]["potential"]
    )
    hessian_cache.store("ts", np.eye(6), identity=base_identity)
    assert hessian_cache.load_matching("ts", other_identity) is None


def test_mm_hessian_mode_aliases_share_one_cache_identity() -> None:
    class _Geom:
        atomic_numbers = [1]
        cart_coords = np.zeros(3)
        freeze_atoms = []

    base = {
        "backend": "uma",
        "mm_hessian_mode": "finite_difference",
    }
    alias = {**base, "mm_hessian_mode": "FD"}

    canonical = hessian_cache.identity_from_context(_Geom(), base, role="ts")
    aliased = hessian_cache.identity_from_context(_Geom(), alias, role="ts")
    assert (
        canonical["evaluator"]["potential"]
        == aliased["evaluator"]["potential"]
    )


def test_persistent_identity_canonicalizes_unbounded_hessian_cutoff() -> None:
    """The freq runtime's +inf sentinel is the default all-movable region."""

    class _Geom:
        atomic_numbers = np.array([1, 8])
        cart_coords = np.zeros(6)
        freeze_atoms = np.array([], dtype=int)

    base = {
        "backend": "uma",
        "uma_model": "uma-s-1p2",
        "uma_precision": "fp32",
        "charge": 0,
        "spin": 1,
    }

    default = hessian_cache.persistent_identity_from_context(
        _Geom(), {**base, "hess_cutoff": None}
    )
    runtime = hessian_cache.persistent_identity_from_context(
        _Geom(), {**base, "hess_cutoff": float("inf")}
    )
    finite = hessian_cache.persistent_identity_from_context(
        _Geom(), {**base, "hess_cutoff": 6.0}
    )

    assert runtime == default
    assert finite != default
    assert finite["evaluator"]["potential"]["hess_cutoff"] == 6.0


def test_persistent_identity_uses_region_file_content_not_location(tmp_path) -> None:
    """Generated region files with identical bytes identify the same PES."""

    class _Geom:
        atomic_numbers = np.array([1, 8])
        cart_coords = np.zeros(6)
        freeze_atoms = np.array([], dtype=int)

    first = tmp_path / "freq" / "model_from_bfactor.pdb"
    second = tmp_path / "irc" / "model_from_bfactor.pdb"
    first.parent.mkdir()
    second.parent.mkdir()
    first.write_bytes(b"IDENTICAL-REGION-CONTENT")
    second.write_bytes(first.read_bytes())
    base = {
        "backend": "uma",
        "uma_model": "uma-s-1p2",
        "uma_precision": "fp32",
        "charge": 0,
        "spin": 1,
    }

    produced = hessian_cache.persistent_identity_from_context(
        _Geom(), {**base, "model_pdb": str(first)}
    )
    consumed = hessian_cache.persistent_identity_from_context(
        _Geom(), {**base, "model_pdb": str(second)}
    )

    assert produced == consumed
    second.write_bytes(b"DIFFERENT-REGION-CONTENT")
    changed = hessian_cache.persistent_identity_from_context(
        _Geom(), {**base, "model_pdb": str(second)}
    )
    assert changed != produced


def test_reconcile_active_hessian_extracts_required_dofs_in_order() -> None:
    source = torch.arange(36, dtype=torch.float64).reshape(6, 6)
    entry = {
        "hessian": source,
        "active_dofs": [3, 4, 5, 0, 1, 2],
    }

    actual = hessian_cache.reconcile_active_hessian(
        entry,
        [0, 1, 2],
        full_n_dof=6,
    )

    expected = source[torch.tensor([3, 4, 5])][:, torch.tensor([3, 4, 5])]
    torch.testing.assert_close(actual, expected)


def test_reconcile_active_hessian_rejects_wrong_or_ambiguous_basis() -> None:
    wrong_same_shape = {
        "hessian": torch.eye(3),
        "active_dofs": [0, 1, 2],
    }
    no_metadata = {"hessian": torch.eye(3), "active_dofs": None}

    assert hessian_cache.reconcile_active_hessian(
        wrong_same_shape,
        [3, 4, 5],
        full_n_dof=6,
    ) is None
    assert hessian_cache.reconcile_active_hessian(
        no_metadata,
        [3, 4, 5],
        full_n_dof=6,
    ) is None


def test_identity_from_context_rejects_backend_specific_model_and_precision() -> None:
    """The cache identity resolves the ML model + precision from the
    BACKEND-PREFIXED keys (uma_model / uma_precision, ...), so changing the ML
    model or fp32-vs-fp64 precision rejects reuse; a genuinely matching context
    still reuses.  The generic ``model`` / ``precision`` keys (which mlmm never
    sets) must NOT participate."""

    class _Geom:
        atomic_numbers = np.array([1, 1, 8])
        cart_coords = np.arange(9, dtype=float)
        freeze_atoms = np.array([], dtype=int)

    import os

    os.environ["MLMM_RUN_ID"] = "run-MP"
    try:
        base_cfg = {
            "backend": "uma",
            "uma_model": "uma-s-1p2",
            "uma_precision": "fp32",
            "charge": 0,
            "spin": 1,
            "freeze_atoms": [],
        }
        ident = hessian_cache.identity_from_context(_Geom(), base_cfg, role="ts")
        hessian_cache.store("ts", np.eye(9), identity=ident)

        # Exact same context reuses.
        assert hessian_cache.load_matching(
            "ts", hessian_cache.identity_from_context(_Geom(), base_cfg, role="ts")
        ) is not None

        # fp32 -> fp64 precision change rejects (scientifically load-bearing).
        assert hessian_cache.load_matching(
            "ts",
            hessian_cache.identity_from_context(
                _Geom(), dict(base_cfg, uma_precision="fp64"), role="ts"
            ),
        ) is None

        # A different ML model rejects.
        assert hessian_cache.load_matching(
            "ts",
            hessian_cache.identity_from_context(
                _Geom(), dict(base_cfg, uma_model="uma-s-other"), role="ts"
            ),
        ) is None

        # The generic model/precision keys are inert: setting them to bogus
        # values (while the backend keys match) must still reuse.
        assert hessian_cache.load_matching(
            "ts",
            hessian_cache.identity_from_context(
                _Geom(),
                dict(base_cfg, model="ignored", precision="ignored"),
                role="ts",
            ),
        ) is not None

        # The resolved identity actually carries the effective model/precision.
        assert ident["evaluator"]["model"] == "uma-s-1p2"
        assert ident["evaluator"]["precision"] == "fp32"
    finally:
        os.environ.pop("MLMM_RUN_ID", None)


def test_identity_from_context_custom_backend_uses_calc_file_as_model() -> None:
    """A custom ASE calculator has no MLIP model variant; its
    calc_file / calc_factory reference stands in as the model identity, so two
    different custom calculators reject reuse."""

    class _Geom:
        atomic_numbers = np.array([1, 1, 8])
        cart_coords = np.zeros(9, dtype=float)
        freeze_atoms = np.array([], dtype=int)

    import os

    os.environ["MLMM_RUN_ID"] = "run-C"
    try:
        cfg = {
            "backend": "custom",
            "calc_file": "/tmp/my_calc.py",
            "calc_factory": "get_calculator",
            "charge": 0,
            "spin": 1,
            "freeze_atoms": [],
        }
        ident = hessian_cache.identity_from_context(_Geom(), cfg, role="ts")
        hessian_cache.store("ts", np.eye(9), identity=ident)
        assert ident["evaluator"]["model"] == "/tmp/my_calc.py:get_calculator"

        assert hessian_cache.load_matching(
            "ts", hessian_cache.identity_from_context(_Geom(), cfg, role="ts")
        ) is not None
        # A different calc_file rejects.
        assert hessian_cache.load_matching(
            "ts",
            hessian_cache.identity_from_context(
                _Geom(), dict(cfg, calc_file="/tmp/other_calc.py"), role="ts"
            ),
        ) is None
        # A different factory rejects.
        assert hessian_cache.load_matching(
            "ts",
            hessian_cache.identity_from_context(
                _Geom(), dict(cfg, calc_factory="build"), role="ts"
            ),
        ) is None
    finally:
        os.environ.pop("MLMM_RUN_ID", None)


def test_custom_calculator_content_change_rejects_same_path_cache(
    monkeypatch,
    tmp_path,
) -> None:
    class _Geom:
        atomic_numbers = np.array([1])
        cart_coords = np.zeros(3, dtype=float)
        freeze_atoms = np.array([], dtype=int)

    calc_file = tmp_path / "calculator.py"
    calc_file.write_text("VALUE = 1\n", encoding="utf-8")
    cfg = {
        "backend": "custom",
        "calc_file": str(calc_file),
        "calc_factory": "get_calculator",
        "charge": 0,
        "spin": 1,
        "freeze_atoms": [],
    }
    monkeypatch.setenv("MLMM_RUN_ID", "run-custom-content")

    first = hessian_cache.identity_from_context(_Geom(), cfg, role="ts")
    hessian_cache.store("ts", np.eye(3), identity=first)
    assert hessian_cache.load_matching(
        "ts", hessian_cache.identity_from_context(_Geom(), cfg, role="ts")
    ) is not None

    calc_file.write_text("VALUE = 2\n", encoding="utf-8")
    changed = hessian_cache.identity_from_context(_Geom(), cfg, role="ts")

    assert first["evaluator"]["potential"]["calc_file_sha256"] != (
        changed["evaluator"]["potential"]["calc_file_sha256"]
    )
    assert hessian_cache.load_matching("ts", changed) is None


def test_mlmm_potential_identity_rejects_parm7_link_embed_region_changes(tmp_path) -> None:
    """The mlmm potential identity rejects a topology-content change,
    a link-method change, an embedding change, and a region-map change."""

    parm7 = tmp_path / "system.parm7"
    parm7.write_bytes(b"ORIGINAL-PRMTOP-CONTENT")

    class _Geom:
        atomic_numbers = np.array([1, 1, 8])
        cart_coords = np.zeros(9, dtype=float)
        freeze_atoms = np.array([], dtype=int)

    base_cfg = {
        "backend": "uma",
        "charge": 0,
        "spin": 1,
        "mm_backend": "openmm",
        "link_atom_method": "ratio",
        "embedcharge": True,
        "use_cmap": True,
        "real_parm7": str(parm7),
        "hess_mm_atoms": [3, 4, 5],
        "freeze_atoms": [],
    }

    monkeyrun = "run-P"
    import os

    os.environ["MLMM_RUN_ID"] = monkeyrun
    try:
        ident = hessian_cache.identity_from_context(_Geom(), base_cfg, role="ts")
        hessian_cache.store("ts", np.eye(9), identity=ident)

        # Exact same context reuses.
        assert hessian_cache.load_matching(
            "ts", hessian_cache.identity_from_context(_Geom(), base_cfg, role="ts")
        ) is not None

        # Replace the parm7's BYTES at the same path -> reject.
        parm7.write_bytes(b"MUTATED-PRMTOP-CONTENT-DIFFERENT")
        assert hessian_cache.load_matching(
            "ts", hessian_cache.identity_from_context(_Geom(), base_cfg, role="ts")
        ) is None
        # Restore original bytes -> reuse again (content, not path, is authoritative).
        parm7.write_bytes(b"ORIGINAL-PRMTOP-CONTENT")
        assert hessian_cache.load_matching(
            "ts", hessian_cache.identity_from_context(_Geom(), base_cfg, role="ts")
        ) is not None

        # Link-method / embedding / region-map changes each reject.
        assert hessian_cache.load_matching(
            "ts",
            hessian_cache.identity_from_context(
                _Geom(), dict(base_cfg, link_atom_method="scaled"), role="ts"
            ),
        ) is None
        assert hessian_cache.load_matching(
            "ts",
            hessian_cache.identity_from_context(
                _Geom(), dict(base_cfg, embedcharge=False), role="ts"
            ),
        ) is None
        assert hessian_cache.load_matching(
            "ts",
            hessian_cache.identity_from_context(
                _Geom(), dict(base_cfg, hess_mm_atoms=[3, 4]), role="ts"
            ),
        ) is None
    finally:
        os.environ.pop("MLMM_RUN_ID", None)


def test_mlmm_potential_identity_rejects_explicit_region_and_link_changes() -> None:
    """The explicit region partition (movable_mm_atoms / frozen_mm_atoms) and
    the link-atom map (link_mlmm) are NOT folded into the captured freeze set,
    so a change to any of them at matching coordinates + run must reject reuse.
    A genuinely matching context still reuses."""

    class _Geom:
        atomic_numbers = np.array([1, 1, 8, 8, 8, 8])
        cart_coords = np.zeros(18, dtype=float)
        freeze_atoms = np.array([], dtype=int)

    import os

    os.environ["MLMM_RUN_ID"] = "run-RL"
    try:
        base_cfg = {
            "backend": "uma",
            "uma_model": "uma-s-1p2",
            "uma_precision": "fp32",
            "charge": 0,
            "spin": 1,
            "movable_mm_atoms": [3, 4],
            "frozen_mm_atoms": [5],
            "link_mlmm": [("A:CYS100:CB", "A:CYS100:CA")],
            "freeze_atoms": [],
        }
        ident = hessian_cache.identity_from_context(_Geom(), base_cfg, role="ts")
        hessian_cache.store("ts", np.eye(18), identity=ident)

        # Exact same context reuses.
        assert hessian_cache.load_matching(
            "ts", hessian_cache.identity_from_context(_Geom(), base_cfg, role="ts")
        ) is not None

        # A movable-region change rejects.
        assert hessian_cache.load_matching(
            "ts",
            hessian_cache.identity_from_context(
                _Geom(), dict(base_cfg, movable_mm_atoms=[3]), role="ts"
            ),
        ) is None
        # A frozen-region change rejects.
        assert hessian_cache.load_matching(
            "ts",
            hessian_cache.identity_from_context(
                _Geom(), dict(base_cfg, frozen_mm_atoms=[4, 5]), role="ts"
            ),
        ) is None
        # A link-boundary change rejects.
        assert hessian_cache.load_matching(
            "ts",
            hessian_cache.identity_from_context(
                _Geom(),
                dict(base_cfg, link_mlmm=[("A:CYS100:SG", "A:CYS100:CB")]),
                role="ts",
            ),
        ) is None

        # Region-partition order is canonical: reordering movable_mm_atoms reuses.
        assert hessian_cache.load_matching(
            "ts",
            hessian_cache.identity_from_context(
                _Geom(), dict(base_cfg, movable_mm_atoms=[4, 3]), role="ts"
            ),
        ) is not None
    finally:
        os.environ.pop("MLMM_RUN_ID", None)
