import numpy as np

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

