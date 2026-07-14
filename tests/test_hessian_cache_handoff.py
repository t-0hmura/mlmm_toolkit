"""Regression tests for geometry-safe in-process Hessian handoff."""

from __future__ import annotations

import numpy as np
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
