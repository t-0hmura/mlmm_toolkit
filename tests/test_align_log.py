"""Regression tests for alignment-selection provenance in user logs."""

import numpy as np
import pytest

from pysisyphus.Geometry import Geometry
import mlmm.workflows.align_freeze as align_freeze
from mlmm.workflows.align_freeze import align_second_to_first_kabsch_inplace
from mlmm.workflows.align_freeze import (
    alignment_failed_pair_indices,
    scan_freeze_atoms_toward_target_inplace,
)


def test_alignment_failure_indices_read_nested_scan_outcome():
    results = [
        {"align": {}, "scan": {"converged": True}},
        {"align": {}, "scan": {"converged": False}},
        {"align": {}, "scan": {}},
        {},
    ]

    assert alignment_failed_pair_indices(results) == [1, 2, 3]


def test_final_relaxation_nonconvergence_is_not_alignment_success(
    tmp_path, monkeypatch
):
    class _NonconvergedLBFGS:
        is_converged = False

        def __init__(self, *_args, **_kwargs):
            pass

        def run(self):
            return None

    ref = Geometry(["H"], np.zeros(3), coord_type="cart")
    mob = Geometry(["H"], np.array([0.01, 0.0, 0.0]), coord_type="cart")
    ref.freeze_atoms = np.array([0], dtype=int)
    mob.freeze_atoms = np.array([0], dtype=int)
    monkeypatch.setattr(align_freeze, "LBFGS", _NonconvergedLBFGS)
    monkeypatch.setattr(
        align_freeze, "_attach_calc_if_needed", lambda *_args, **_kwargs: None
    )

    outcome = scan_freeze_atoms_toward_target_inplace(
        ref,
        mob,
        step_A=0.1,
        out_dir=tmp_path,
        verbose=False,
    )

    assert outcome["max_remaining_A"] == 0.0
    assert outcome["converged"] is False


@pytest.mark.parametrize(
    ("freeze_atoms", "expected"),
    [
        ([], "used 3 atoms"),
        ([0, 1, 2], "used 3 freeze atoms"),
    ],
)
def test_kabsch_log_names_the_alignment_selection(
    monkeypatch, freeze_atoms, expected
):
    coords = np.array(
        [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]]
    )
    ref = Geometry(["H", "H", "H"], coords.reshape(-1), coord_type="cart")
    mob = Geometry(
        ["H", "H", "H"],
        (coords + np.array([0.2, -0.1, 0.3])).reshape(-1),
        coord_type="cart",
    )
    ref.freeze_atoms = np.asarray(freeze_atoms, dtype=int)
    mob.freeze_atoms = np.asarray(freeze_atoms, dtype=int)
    messages = []
    monkeypatch.setattr(
        align_freeze,
        "emit",
        lambda message, **kwargs: messages.append(message),
    )

    align_second_to_first_kabsch_inplace(ref, mob, verbose=True)

    assert expected in messages[-1]


def test_final_step_log_states_the_achieved_coincidence(tmp_path, monkeypatch):
    """The step line reports the pre-snap gap, so the snap itself must be logged."""

    class _ConvergedLBFGS:
        is_converged = True

        def __init__(self, *_args, **_kwargs):
            pass

        def run(self):
            return None

    ref = Geometry(["H"], np.zeros(3), coord_type="cart")
    mob = Geometry(["H"], np.array([0.05, 0.0, 0.0]), coord_type="cart")
    ref.freeze_atoms = np.array([0], dtype=int)
    mob.freeze_atoms = np.array([0], dtype=int)
    monkeypatch.setattr(align_freeze, "LBFGS", _ConvergedLBFGS)
    monkeypatch.setattr(
        align_freeze, "_attach_calc_if_needed", lambda *_args, **_kwargs: None
    )
    messages = []
    monkeypatch.setattr(
        align_freeze, "emit", lambda message, **kwargs: messages.append(message)
    )

    outcome = scan_freeze_atoms_toward_target_inplace(
        ref, mob, step_A=0.1, out_dir=tmp_path, verbose=True,
    )

    snap = [m for m in messages if "anchors set to the reference" in m]
    assert len(snap) == 1, messages
    assert "-> 0.000000" in snap[0]
    assert "frozen" in snap[0]
    assert outcome["max_remaining_A"] == 0.0
    assert np.allclose(mob.cart_coords, ref.cart_coords)
