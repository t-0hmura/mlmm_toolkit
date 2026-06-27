"""Unit tests for mlmm.domain.bond_changes pure helper functions."""

from __future__ import annotations

from types import SimpleNamespace

import numpy as np
import pytest


class TestBondStr:
    def test_one_based(self):
        from mlmm.domain.bond_changes import _bond_str

        result = _bond_str(0, 1, ["C", "O"], one_based=True)
        assert result == "C1-O2"

    def test_zero_based(self):
        from mlmm.domain.bond_changes import _bond_str

        result = _bond_str(0, 1, ["C", "O"], one_based=False)
        assert result == "C0-O1"

    def test_same_element(self):
        from mlmm.domain.bond_changes import _bond_str

        result = _bond_str(2, 5, ["H", "H", "C", "N", "O", "C"], one_based=True)
        assert result == "C3-C6"


class TestElementArrays:
    def test_basic(self):
        from mlmm.domain.bond_changes import _element_arrays

        elems, radii = _element_arrays(["C", "O", "N"])
        assert elems == ["C", "O", "N"]
        assert len(radii) == 3
        # Carbon covalent radius should be positive
        assert radii[0] > 0

    def test_lowercase_normalization(self):
        from mlmm.domain.bond_changes import _element_arrays

        elems, radii = _element_arrays(["c", "o"])
        assert elems == ["C", "O"]

    def test_empty(self):
        from mlmm.domain.bond_changes import _element_arrays

        elems, radii = _element_arrays([])
        assert elems == []
        assert len(radii) == 0


class TestResolveDevice:
    def test_cpu(self):
        import torch
        from mlmm.domain.bond_changes import _resolve_device

        dev = _resolve_device("cpu")
        assert dev == torch.device("cpu")

    def test_auto_returns_device(self):
        import torch
        from mlmm.domain.bond_changes import _resolve_device

        dev = _resolve_device("auto")
        assert isinstance(dev, torch.device)


# Bohr; well inside the C-C covalent cutoff (~3.27 Bohr after the default margin)
_BOND_BOHR = 2.8


def _line_geom(coords):
    """Minimal geom stub: all-carbon chain with explicit Bohr coordinates."""
    coords = np.asarray(coords, dtype=float)
    return SimpleNamespace(atoms=["C"] * len(coords), coords3d=coords)


def _as_int_pairs(pairs):
    # compare_structures returns numpy-int tuples; normalize for comparison
    return {tuple(int(x) for x in p) for p in pairs}


class TestCompareStructures:
    def test_chunk_boundary_formed_and_broken(self):
        """Row-chunked detection must be correct across the 1024-row block
        boundary: one bond forms in the first block, one breaks in the second
        (exercises the per-block upper-triangle offset and the ``+= i0`` global
        index correction)."""
        from mlmm.domain.bond_changes import compare_structures

        N = 1026
        base = np.zeros((N, 3), dtype=float)
        base[:, 0] = np.arange(N) * 20.0  # 20 Bohr spacing -> no bonds anywhere

        r1 = base.copy()
        r2 = base.copy()
        # formed (2, 3): apart in r1, bonded in r2  -> first row block
        r2[3, 0] = base[2, 0] + _BOND_BOHR
        # broken (1024, 1025): bonded in r1, apart in r2  -> second row block
        r1[1025, 0] = base[1024, 0] + _BOND_BOHR

        res = compare_structures(_line_geom(r1), _line_geom(r2), device="cpu")

        assert _as_int_pairs(res.formed_covalent) == {(2, 3)}
        assert _as_int_pairs(res.broken_covalent) == {(1024, 1025)}
        # memory-bounded path: full N x N matrices dropped, sparse lengths kept
        assert res.distances_1 is None and res.distances_2 is None
        lengths = {tuple(int(x) for x in k): v for k, v in res.changed_lengths.items()}
        assert lengths[(2, 3)][1] == pytest.approx(_BOND_BOHR, abs=1e-6)      # d2 (formed)
        assert lengths[(1024, 1025)][0] == pytest.approx(_BOND_BOHR, abs=1e-6)  # d1 (broken)

    def test_summarize_uses_sparse_lengths(self):
        """summarize_changes still prints the ``D1 Å --> D2 Å`` lengths from the
        sparse map even though the full distance matrices are gone."""
        from mlmm.domain.bond_changes import compare_structures, summarize_changes

        r1 = [[0.0, 0.0, 0.0], [20.0, 0.0, 0.0]]            # apart
        r2 = [[0.0, 0.0, 0.0], [_BOND_BOHR, 0.0, 0.0]]      # bonded -> formed
        res = compare_structures(_line_geom(r1), _line_geom(r2), device="cpu")

        assert _as_int_pairs(res.formed_covalent) == {(0, 1)}
        text = summarize_changes(_line_geom(r2), res, one_based=True)
        assert "C1-C2" in text
        assert "Å -->" in text
