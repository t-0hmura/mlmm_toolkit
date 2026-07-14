"""Tests for carrying the MEP reaction direction into TS optimization."""

from __future__ import annotations

import numpy as np
from ase import Atoms
from ase.io import write

from mlmm.workflows.all import _ensure_hei_path_tangent
from mlmm.workflows.path_search import _write_xyz_trj_with_energy_from_ase


def test_hei_path_tangent_uses_neighbors_of_matching_image(tmp_path) -> None:
    images = [
        Atoms("H2", positions=positions)
        for positions in (
            [[0.0, 0.0, 0.0], [0.7, -0.2, 0.0]],
            [[0.0, 0.0, 0.0], [0.7, 0.0, 0.0]],
            [[0.0, 0.0, 0.0], [1.0, 0.3, 0.0]],
        )
    ]
    mep = tmp_path / "mep_trj.xyz"
    hei = tmp_path / "hei.xyz"
    mode = tmp_path / "hei_mode.txt"
    write(mep, images)
    write(hei, images[1])

    assert _ensure_hei_path_tangent(mep, hei, mode) == mode
    tangent = np.loadtxt(mode)
    incoming = (images[1].positions - images[0].positions).reshape(-1)
    outgoing = (images[2].positions - images[1].positions).reshape(-1)
    incoming /= np.linalg.norm(incoming)
    outgoing /= np.linalg.norm(outgoing)
    expected = incoming + outgoing
    expected /= np.linalg.norm(expected)
    np.testing.assert_allclose(tangent, expected)


def test_hei_path_tangent_uses_energy_upwinding_when_available(tmp_path) -> None:
    images = [
        Atoms("H", positions=positions)
        for positions in (
            [[-2.0, 0.0, 0.0]],
            [[0.0, 0.0, 0.0]],
            [[0.0, 1.0, 0.0]],
        )
    ]
    mep = tmp_path / "mep_trj.xyz"
    hei = tmp_path / "hei.xyz"
    mode = tmp_path / "hei_mode.txt"
    _write_xyz_trj_with_energy_from_ase(images, [0.0, 2.0, 0.5], mep)
    write(hei, images[1])

    assert _ensure_hei_path_tangent(mep, hei, mode) == mode
    outgoing = (images[2].positions - images[1].positions).reshape(-1)
    incoming = (images[1].positions - images[0].positions).reshape(-1)
    expected = 2.0 * outgoing + 1.5 * incoming
    expected /= np.linalg.norm(expected)
    np.testing.assert_allclose(np.loadtxt(mode), expected)
