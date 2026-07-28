"""Regression tests for authoritative final-frequency Hessian bases and the
HessianDimer cycle/coordinate guards."""

from __future__ import annotations

from types import SimpleNamespace

import click
import numpy as np
import pytest
import torch

from mlmm.workflows.freq import _reconcile_hessian_analysis_basis
from mlmm.workflows import tsopt
from mlmm.workflows.tsopt import HessianDimer


def _geometry(n_atoms: int, coverage: list[int]) -> SimpleNamespace:
    active_dofs = [
        3 * atom + axis for atom in coverage for axis in range(3)
    ]
    return SimpleNamespace(
        atomic_numbers=np.ones(n_atoms, dtype=int),
        within_partial_hessian={
            "active_n_dof": len(active_dofs),
            "full_n_dof": 3 * n_atoms,
            "active_atoms": coverage,
            "active_dofs": active_dofs,
        },
        _hess_active_atoms_last=np.asarray(coverage, dtype=int),
        _hess_active_dofs_last=np.asarray(active_dofs, dtype=int),
    )


def test_compact_superset_is_sliced_in_declared_coverage_order() -> None:
    geometry = _geometry(3, [2, 0, 1])
    source = torch.arange(81, dtype=torch.float64).reshape(9, 9)

    actual, requested, coverage, storage = _reconcile_hessian_analysis_basis(
        source, geometry, [1, 0],
    )

    local = torch.tensor([3, 4, 5, 6, 7, 8])
    expected = source.index_select(0, local).index_select(1, local)
    torch.testing.assert_close(actual, expected)
    assert requested == [0, 1]
    assert coverage == [2, 0, 1]
    assert storage == "compact"


def test_matching_basis_reuses_the_original_hessian() -> None:
    geometry = _geometry(3, [0, 1, 2])
    source = torch.eye(9, dtype=torch.float64)

    actual, requested, coverage, storage = _reconcile_hessian_analysis_basis(
        source, geometry, [0, 1, 2],
    )

    assert actual is source
    assert requested == [0, 1, 2]
    assert coverage == [0, 1, 2]
    assert storage == "compact"


def test_full_shape_does_not_claim_curvature_outside_coverage() -> None:
    geometry = _geometry(3, [0, 1])
    zero_expanded = torch.zeros((9, 9), dtype=torch.float64)

    with pytest.raises(click.ClickException, match="missing 1-based atom indices: 3"):
        _reconcile_hessian_analysis_basis(zero_expanded, geometry, [2])


def test_compact_hessian_without_ordered_metadata_fails_closed() -> None:
    geometry = SimpleNamespace(
        atomic_numbers=np.ones(3, dtype=int),
        within_partial_hessian=None,
    )
    with pytest.raises(click.ClickException, match="without ordered coverage"):
        _reconcile_hessian_analysis_basis(
            torch.eye(6, dtype=torch.float64), geometry, [0, 1],
        )


def test_tsopt_rejects_configured_basis_gap_before_optimization(
    monkeypatch,
) -> None:
    layers = {
        "ml": {0, 1},
        "hess_mm": {2},
        "movable_mm": {3},
        "frozen_mm": set(),
    }

    def resolve(_cfg, _n_atoms, mode):
        if mode == "ml-only":
            return {0, 1}, layers
        return {0, 1, 2, 3}, layers

    monkeypatch.setattr(tsopt, "_resolve_active_atom_indices", resolve)

    with pytest.raises(
        click.ClickException,
        match=r"configured Hessian coverage.*indices: 4",
    ):
        tsopt._resolve_validated_hessian_analysis_atoms(
            {}, 4, "partial", [], validate_coverage=True
        )

    assert tsopt._resolve_validated_hessian_analysis_atoms(
        {}, 4, "ml-only", [], validate_coverage=True
    ) == [0, 1]
    assert tsopt._resolve_validated_hessian_analysis_atoms(
        {}, 4, "partial", [], validate_coverage=False
    ) == [0, 1, 2, 3]


def test_dimer_reserves_the_strict_threshold_cycle() -> None:
    runner = HessianDimer.__new__(HessianDimer)
    runner.max_total_cycles = 1
    runner._cycles_spent = 0
    runner.update_interval_hessian = 10
    runner.is_stalled = False
    runner.geom = SimpleNamespace(cart_coords=np.zeros(3))
    calls: list[tuple[str, int]] = []

    def segment(threshold: str, steps: int) -> tuple[int, bool]:
        calls.append((threshold, steps))
        return 1, True

    runner._dimer_segment = segment
    assert runner._dimer_loop("gau_loose", reserve_cycles=1) == (0, False, False)
    assert runner._dimer_loop("baker") == (1, True, True)
    assert calls == [("baker", 1)]


def test_direct_dimer_rejects_internal_coordinates(tmp_path) -> None:
    with pytest.raises(ValueError, match="coord_type must be 'cart'"):
        HessianDimer(
            "unused.xyz",
            out_dir=str(tmp_path / "ts"),
            geom_kwargs={"coord_type": "dlc"},
        )


def test_dimer_hessian_cache_separates_guidance_and_exact_curvature(
    monkeypatch,
) -> None:
    runner = HessianDimer.__new__(HessianDimer)
    runner.geom = SimpleNamespace(
        atoms=("H",),
        atomic_numbers=np.array([1]),
        cart_coords=np.zeros(3),
        freeze_atoms=[],
    )
    runner.device = torch.device("cpu")
    runner._raw_hessian_cache_cpu = None
    runner._raw_hessian_coords_cpu = None
    runner._raw_hessian_identity = None
    runner._compact_hessian_to_computed_coverage = lambda value: value
    calls = []

    def calculate(geometry, kwargs, device):
        calls.append(kwargs["mm_hessian_mode"])
        value = 1.0 if kwargs["mm_hessian_mode"] == "none" else 2.0
        return torch.eye(3, dtype=torch.float64) * value

    monkeypatch.setattr(tsopt, "_calc_full_hessian_torch", calculate)
    guide = {
        "backend": "uma",
        "mm_hessian_mode": "none",
        "return_partial_hessian": False,
    }
    exact = {
        "backend": "uma",
        "mm_hessian_mode": "analytical",
        "return_partial_hessian": False,
    }

    first = runner._calc_full_hessian_cached(guide, allow_reuse=False)
    second = runner._calc_full_hessian_cached(exact, allow_reuse=True)
    third = runner._calc_full_hessian_cached(exact, allow_reuse=True)

    torch.testing.assert_close(first, torch.eye(3, dtype=torch.float64))
    torch.testing.assert_close(second, torch.eye(3, dtype=torch.float64) * 2.0)
    torch.testing.assert_close(third, second)
    assert calls == ["none", "analytical"]
