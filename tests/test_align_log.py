"""Regression tests for alignment-selection provenance in user logs."""

import numpy as np
import pytest

from pysisyphus.Geometry import Geometry
import mlmm.workflows.align_freeze as align_freeze
from mlmm.workflows.align_freeze import align_second_to_first_kabsch_inplace


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
