import numpy as np
import pytest

from mlmm.workflows import freq


def test_partial_active_dof_includes_hessian_target_mm(monkeypatch):
    monkeypatch.setattr(
        freq,
        "_collect_layer_atom_sets",
        lambda _cfg: {
            "ml": {0},
            "hess_mm": {1, 3},
            "movable_mm": {2},
            "frozen_mm": {4},
        },
    )

    active, _ = freq._resolve_active_atom_indices({}, 5, "partial")

    assert active == {0, 1, 2, 3}


def test_partial_hessian_metadata_must_match_shape():
    class DummyGeom:
        pass

    geom = DummyGeom()
    geom.within_partial_hessian = {
        "active_atoms": np.array([0, 2, 4]),
        "active_dofs": np.array([0, 1, 2, 6, 7, 8, 12, 13, 14]),
    }

    assert freq._active_atoms_from_partial_hessian_metadata(geom, 9) == [0, 2, 4]
    assert freq._active_atoms_from_partial_hessian_metadata(geom, 12) is None


def test_movable_cutoff_becomes_the_default_hessian_cutoff():
    calc_cfg = {
        "use_bfactor_layers": False,
        "movable_cutoff": 8.0,
        "hess_cutoff": None,
    }

    assert freq._align_three_layer_hessian_targets(calc_cfg) is True
    assert calc_cfg["hess_cutoff"] == 8.0


def test_hessian_basis_gap_is_rejected_before_evaluation():
    with pytest.raises(Exception, match="before evaluation"):
        freq._validate_hessian_basis_coverage(
            {0, 1, 2},
            {
                "ml": {0},
                "hess_mm": {1},
                "movable_mm": {2},
                "frozen_mm": set(),
            },
            set(),
        )
