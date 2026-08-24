"""Tests for carrying the MEP reaction direction into TS optimization."""

from __future__ import annotations

import numpy as np
import click
import pytest
from ase import Atoms
from ase.io import write

from mlmm.workflows.all import (
    _ensure_hei_path_tangent,
    _required_xyz_block_energies,
    cli,
)
from mlmm.workflows.path_search import _write_xyz_trj_with_energy_from_ase
from mlmm.io.path_mode_cache import (
    build_path_mode_candidates,
    load_path_mode_cache,
    write_path_mode_cache,
)


def _replace_xyz_comments(path, comments) -> None:
    lines = path.read_text(encoding="utf-8").splitlines()
    cursor = 0
    for comment in comments:
        atom_count = int(lines[cursor])
        lines[cursor + 1] = comment
        cursor += atom_count + 2
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def test_mep_tangent_handoff_is_default_on_and_can_be_disabled() -> None:
    option = next(param for param in cli.params if param.name == "tsopt_from_mep_tan")
    assert option.default is True
    assert "--no-tsopt-from-mep-tan" in option.opts + option.secondary_opts
    assert all(param.name != "use_mep_tangent" for param in cli.params)
    assert all(param.name != "freq_symmetry_number" for param in cli.params)


def test_disabled_mep_tangent_does_not_build_a_reference_mode(
    tmp_path, monkeypatch
) -> None:
    from mlmm.workflows import all as all_workflow

    monkeypatch.setattr(
        all_workflow,
        "_ensure_hei_path_tangent",
        lambda *args: pytest.fail("disabled tangent builder was called"),
    )
    path, status = all_workflow._select_hei_reference_mode(
        False,
        tmp_path / "mep.xyz",
        tmp_path / "hei.xyz",
        tmp_path / "mode.txt",
    )
    assert path is None
    assert status["status"] == "disabled"


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
    _replace_xyz_comments(
        mep,
        [
            "mode 1 -120.0 cm-1 frame=0",
            "mode 1 -80.0 cm-1 frame=1",
            "mode 1 -40.0 cm-1 frame=2",
        ],
    )
    write(hei, images[1])

    cache_path, status = _ensure_hei_path_tangent(mep, hei, mode)
    assert cache_path == mode.with_suffix(".npz")
    assert status["status"] == "created"
    assert status["candidate_count"] >= 1
    tangent = np.loadtxt(mode.with_suffix(".txt"))
    incoming = (images[1].positions - images[0].positions).reshape(-1)
    outgoing = (images[2].positions - images[1].positions).reshape(-1)
    incoming /= np.linalg.norm(incoming)
    outgoing /= np.linalg.norm(outgoing)
    expected = incoming + outgoing
    expected /= np.linalg.norm(expected)
    np.testing.assert_allclose(tangent, expected)


def test_required_segment_energies_fail_on_one_ambiguous_frame(tmp_path) -> None:
    path = tmp_path / "segment.xyz"
    blocks = [
        ["1", "-1.0", "H 0 0 0"],
        ["1", "mode 1 -100.0 cm-1", "H 0 0 0"],
        ["1", "E=-0.5 Ha", "H 0 0 0"],
    ]
    with pytest.raises(click.ClickException, match=r"frame 2 .*ambiguous"):
        _required_xyz_block_energies(
            blocks, path=path, context="path-opt segment 1"
        )
    assert _required_xyz_block_energies(
        [blocks[0], blocks[2]], path=path, context="path-opt segment 1"
    ) == [-1.0, -0.5]


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

    cache_path, status = _ensure_hei_path_tangent(mep, hei, mode)
    assert cache_path == mode.with_suffix(".npz")
    assert status["status"] == "created"
    outgoing = (images[2].positions - images[1].positions).reshape(-1)
    incoming = (images[1].positions - images[0].positions).reshape(-1)
    expected = 2.0 * outgoing + 1.5 * incoming
    expected /= np.linalg.norm(expected)
    np.testing.assert_allclose(np.loadtxt(mode.with_suffix(".txt")), expected)


def test_path_mode_cache_is_cpu_only_reusable_and_rejects_stale_inputs(tmp_path) -> None:
    trajectory = tmp_path / "mep.xyz"
    hei = tmp_path / "hei.xyz"
    mode = tmp_path / "hei_mode.txt"
    images = [Atoms("H", positions=[[0.0, 0.0, 0.0]]), Atoms("H", positions=[[1.0, 0.0, 0.0]])]
    write(trajectory, images)
    write(hei, images[1])
    cache = write_path_mode_cache(
        mode.with_suffix(".npz"),
        [im.positions for im in images],
        1,
        energies=[0.0, -1.25],
        trajectory_path=trajectory,
        hei_path=hei,
        atom_numbers=[1],
        primary_text_path=mode,
    )
    assert cache is not None
    assert cache.path == mode.with_suffix(".npz")
    assert cache.path.with_suffix(".json").exists()
    loaded = load_path_mode_cache(
        cache.path,
        trajectory_path=trajectory,
        hei_path=hei,
        expected_size=3,
        atom_numbers=[1],
    )
    assert loaded.metadata["image_index"] == 1
    np.testing.assert_allclose(loaded.primary, [1.0, 0.0, 0.0])
    trajectory.write_text(trajectory.read_text(encoding="utf-8") + "\n", encoding="utf-8")
    with pytest.raises(ValueError, match="trajectory fingerprint"):
        load_path_mode_cache(
            cache.path, trajectory_path=trajectory, hei_path=hei
        )
