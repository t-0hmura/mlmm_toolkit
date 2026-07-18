"""Tests for geometry-identified Hessian file handoff."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

from mlmm.core import result_commit
from mlmm.core.result_commit import ResultCommitError
from mlmm.io.hessian_file import load_hessian_file, save_hessian_file
from mlmm.workflows.freq import _record_hessian_result_path


def test_partial_hessian_round_trip_preserves_identity_and_active_dofs(tmp_path) -> None:
    path = tmp_path / "partial.npz"
    coords = np.arange(9, dtype=float) / 10.0
    numbers = np.array([6, 1, 8])
    active = [0, 1, 2, 6, 7, 8]
    save_hessian_file(
        path,
        hessian=np.eye(6),
        energy_ha=-1.25,
        cart_coords_bohr=coords,
        atomic_numbers=numbers,
        partial_metadata={
            "active_dofs": active,
            "active_n_dof": 6,
            "full_n_dof": 9,
        },
    )

    loaded = load_hessian_file(
        path,
        cart_coords_bohr=coords + 1.0e-5,
        atomic_numbers=numbers,
        expected_active_dofs=active,
    )

    np.testing.assert_allclose(loaded["hessian"], np.eye(6))
    assert loaded["energy_ha"] == pytest.approx(-1.25)
    assert loaded["partial_metadata"] == {
        "active_n_dof": 6,
        "full_n_dof": 9,
        "active_dofs": active,
        "active_atoms": [0, 2],
    }


@pytest.mark.parametrize(
    ("coords", "numbers", "message"),
    [
        (np.array([0.0, 0.0, 0.1]), np.array([1]), "coordinates"),
        (np.zeros(3), np.array([8]), "atomic numbers/order"),
    ],
)
def test_hessian_file_rejects_wrong_geometry(tmp_path, coords, numbers, message) -> None:
    path = tmp_path / "full.npz"
    save_hessian_file(
        path,
        hessian=np.eye(3),
        energy_ha=0.0,
        cart_coords_bohr=np.zeros(3),
        atomic_numbers=np.array([1]),
    )
    with pytest.raises(ValueError, match=message):
        load_hessian_file(path, cart_coords_bohr=coords, atomic_numbers=numbers)


def test_hessian_file_rejects_legacy_unidentified_npz(tmp_path) -> None:
    path = tmp_path / "legacy.npz"
    np.savez_compressed(path, hessian=np.eye(3), energy_ha=0.0)
    with pytest.raises(ValueError, match="lacks geometry identity metadata"):
        load_hessian_file(
            path,
            cart_coords_bohr=np.zeros(3),
            atomic_numbers=np.array([1]),
        )


def test_hessian_file_rejects_inconsistent_active_metadata(tmp_path) -> None:
    path = tmp_path / "bad-active.npz"
    np.savez_compressed(
        path,
        schema_version=np.int64(1),
        hessian=np.eye(3),
        energy_ha=0.0,
        cart_coords_bohr=np.zeros(6),
        atomic_numbers=np.array([1, 1]),
        wph_active_dofs=np.array([0, 1]),
        wph_active_n_dof=np.int64(2),
        wph_full_n_dof=np.int64(6),
    )
    with pytest.raises(ValueError, match="inconsistent active-DOF metadata"):
        load_hessian_file(
            path,
            cart_coords_bohr=np.zeros(6),
            atomic_numbers=np.array([1, 1]),
        )


def test_hessian_file_rejects_different_current_active_basis(tmp_path) -> None:
    path = tmp_path / "partial.npz"
    coords = np.zeros(9)
    numbers = np.array([6, 1, 8])
    save_hessian_file(
        path,
        hessian=np.eye(6),
        energy_ha=0.0,
        cart_coords_bohr=coords,
        atomic_numbers=numbers,
        partial_metadata={
            "active_dofs": [0, 1, 2, 6, 7, 8],
            "active_n_dof": 6,
            "full_n_dof": 9,
        },
    )

    with pytest.raises(ValueError, match="active-DOF basis"):
        load_hessian_file(
            path,
            cart_coords_bohr=coords,
            atomic_numbers=numbers,
            expected_active_dofs=[0, 1, 2, 3, 4, 5],
        )


@pytest.mark.parametrize("name", ["hessian", "hessian.bin", "hessian.npz"])
def test_hessian_save_uses_exact_requested_path(tmp_path, name: str) -> None:
    requested = tmp_path / name
    coords = np.arange(6, dtype=float) / 10.0
    numbers = np.array([1, 8])
    returned = save_hessian_file(
        requested,
        hessian=np.eye(6),
        energy_ha=-2.5,
        cart_coords_bohr=coords,
        atomic_numbers=numbers,
    )

    assert returned == requested
    assert requested.exists()
    assert not Path(str(requested) + ".npz").exists()
    loaded = load_hessian_file(
        requested,
        cart_coords_bohr=coords,
        atomic_numbers=numbers,
    )
    np.testing.assert_allclose(loaded["hessian"], np.eye(6))
    assert loaded["energy_ha"] == pytest.approx(-2.5)

    files = _record_hessian_result_path({"frequencies_txt": "frequencies_cm-1.txt"}, returned)
    assert files["hessian_npz"] == str(requested)
    assert Path(files["hessian_npz"]).exists()


def test_hessian_atomic_publish_failure_preserves_old_exact_artifact(
    tmp_path, monkeypatch: pytest.MonkeyPatch
) -> None:
    requested = tmp_path / "hessian.bin"
    old = b"old-hessian"
    requested.write_bytes(old)

    def fail_replace(source, destination):
        raise OSError("injected")

    monkeypatch.setattr(result_commit.os, "replace", fail_replace)
    with pytest.raises(ResultCommitError, match="publish"):
        save_hessian_file(
            requested,
            hessian=np.eye(3),
            energy_ha=0.0,
            cart_coords_bohr=np.zeros(3),
            atomic_numbers=np.array([1]),
        )
    assert requested.read_bytes() == old
    assert not Path(str(requested) + ".npz").exists()
    assert list(tmp_path.glob(".*.tmp")) == []
