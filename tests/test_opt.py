"""Unit tests for lightweight helper functions in mlmm.opt."""

from __future__ import annotations

import importlib
import sys
import types

import numpy as np
import pytest

pytestmark = pytest.mark.skipif(
    sys.version_info < (3, 11),
    reason="mlmm requires Python >= 3.11",
)


@pytest.fixture(scope="module")
def opt_module():
    """Import mlmm.opt with mlmm_calc stubbed out."""
    stub = types.ModuleType("mlmm.mlmm_calc")

    def _stub_mlmm(*_args, **_kwargs):
        raise RuntimeError("stub mlmm should not be called in helper tests")

    stub.mlmm = _stub_mlmm
    stub.mlmm_mm_only = _stub_mlmm
    sys.modules["mlmm.mlmm_calc"] = stub
    sys.modules.pop("mlmm.opt", None)
    mod = importlib.import_module("mlmm.opt")
    return mod


def test_parse_freeze_atoms_valid(opt_module):
    parsed = opt_module._parse_freeze_atoms("3,1,3,2")
    assert parsed == [0, 1, 2]


def test_parse_freeze_atoms_rejects_invalid(opt_module):
    with pytest.raises(Exception, match="Invalid integer"):
        opt_module._parse_freeze_atoms("1,x,3")

    with pytest.raises(Exception, match="1-based positive"):
        opt_module._parse_freeze_atoms("0,2")


def test_normalize_geom_freeze(opt_module):
    assert opt_module._normalize_geom_freeze("3,1,2") == [1, 2, 3]
    assert opt_module._normalize_geom_freeze([5, 3, 4]) == [3, 4, 5]

    with pytest.raises(Exception, match="must contain integers"):
        opt_module._normalize_geom_freeze("1,a")


def test_parse_dist_freeze_valid_and_index_conversion(opt_module):
    parsed = opt_module._parse_dist_freeze_args(["[(1,2,1.5),(2,3)]"], one_based=True, atom_meta=None)
    assert parsed == [(0, 1, 1.5), (1, 2, None)]


def test_parse_dist_freeze_rejects_invalid_target(opt_module):
    with pytest.raises(Exception, match="Target distance must be > 0"):
        opt_module._parse_dist_freeze_args(["(1,2,0.0)"], one_based=True, atom_meta=None)


def test_resolve_dist_freeze_targets_uses_current_distance_when_none(opt_module):
    class DummyGeom:
        # coords3d in Bohr; atom distance is exactly 1.0 Å
        coords3d = np.array([
            [0.0, 0.0, 0.0],
            [1.0 / opt_module.BOHR2ANG, 0.0, 0.0],
        ])

    resolved = opt_module._resolve_dist_freeze_targets(DummyGeom(), [(0, 1, None), (0, 1, 2.5)])
    assert abs(resolved[0][2] - 1.0) < 1e-10
    assert resolved[1] == (0, 1, 2.5)


def test_pdb_key_parser_handles_missing_columns(opt_module):
    full_line = (
        "ATOM      1  CA  ALA A  10      10.000  11.000  12.000  1.00 20.00           C\n"
    )
    key_full, key_simple = opt_module._pdb_keys_from_line(full_line)
    assert key_full == ("A", 10, "", "ALA", "CA", "")
    assert key_simple == ("A", 10, "", "CA")

    short_line = "ATOM      1  N   ALA\n"
    key_full2, key_simple2 = opt_module._pdb_keys_from_line(short_line)
    assert key_full2[1] == -(10 ** 9)
    assert key_simple2[1] == -(10 ** 9)


def test_format_with_bfactor_pads_short_lines(opt_module):
    short_line = "ATOM      1  N   ALA A   1\n"
    out = opt_module._format_with_bfactor(short_line, 50.0)

    assert len(out) >= 67
    assert out[60:66] == " 50.00"
