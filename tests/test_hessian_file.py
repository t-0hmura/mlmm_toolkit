"""Tests for geometry-identified Hessian file handoff."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

from mlmm.core import result_commit
from mlmm.core.result_commit import ResultCommitError
from mlmm.io.hessian_file import load_hessian_file, save_hessian_file
from mlmm.workflows.freq import _record_hessian_result_path


PES_IDENTITY = {
    "schema": "hessian-cache-identity/v1",
    "system": {"atoms": [6, 1, 8]},
    "evaluator": {
        "backend": "uma",
        "model": "uma-s-1p1",
        "precision": "float64",
        "potential": {"mm_backend": "hessian_ff"},
    },
}
SAVE_STATE = {
    "model_charge": 0,
    "model_mult": 1,
    "potential_identity": PES_IDENTITY,
}
LOAD_STATE = {
    "expected_model_charge": 0,
    "expected_model_mult": 1,
    "expected_potential_identity": PES_IDENTITY,
}


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
        **SAVE_STATE,
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
        **LOAD_STATE,
    )

    np.testing.assert_allclose(loaded["hessian"], np.eye(6))
    assert loaded["energy_ha"] == pytest.approx(-1.25)
    assert loaded["schema_version"] == 3
    assert loaded["model_charge"] == 0
    assert loaded["model_mult"] == 1
    assert loaded["electronic_state_verified"] is True
    assert loaded["potential_identity_verified"] is True
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
        **SAVE_STATE,
    )
    with pytest.raises(ValueError, match=message):
        load_hessian_file(
            path,
            cart_coords_bohr=coords,
            atomic_numbers=numbers,
            **LOAD_STATE,
        )


def test_hessian_file_rejects_legacy_unidentified_npz(tmp_path) -> None:
    path = tmp_path / "legacy.npz"
    np.savez_compressed(path, hessian=np.eye(3), energy_ha=0.0)
    with pytest.raises(ValueError, match="lacks geometry identity metadata"):
        load_hessian_file(
            path,
            cart_coords_bohr=np.zeros(3),
            atomic_numbers=np.array([1]),
            **LOAD_STATE,
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
            allow_unverified_state=True,
            allow_unverified_pes=True,
            **LOAD_STATE,
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
        **SAVE_STATE,
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
            **LOAD_STATE,
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
        **SAVE_STATE,
    )

    assert returned == requested
    assert requested.exists()
    assert not Path(str(requested) + ".npz").exists()
    loaded = load_hessian_file(
        requested,
        cart_coords_bohr=coords,
        atomic_numbers=numbers,
        **LOAD_STATE,
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
            **SAVE_STATE,
        )
    assert requested.read_bytes() == old
    assert not Path(str(requested) + ".npz").exists()
    assert list(tmp_path.glob(".*.tmp")) == []


@pytest.mark.parametrize(
    ("expected_charge", "expected_mult"),
    [(1, 1), (0, 3)],
)
def test_hessian_file_rejects_electronic_state_mismatch_even_with_override(
    tmp_path, expected_charge, expected_mult
) -> None:
    path = tmp_path / "state.npz"
    save_hessian_file(
        path,
        hessian=np.eye(3),
        energy_ha=0.0,
        cart_coords_bohr=np.zeros(3),
        atomic_numbers=np.array([1]),
        **SAVE_STATE,
    )
    with pytest.raises(ValueError, match="electronic state"):
        load_hessian_file(
            path,
            cart_coords_bohr=np.zeros(3),
            atomic_numbers=np.array([1]),
            expected_model_charge=expected_charge,
            expected_model_mult=expected_mult,
            allow_unverified_state=True,
        )


def test_schema_one_requires_explicit_unverified_state_opt_in(tmp_path) -> None:
    path = tmp_path / "schema-one.npz"
    np.savez_compressed(
        path,
        schema_version=np.int64(1),
        hessian=np.eye(3),
        energy_ha=0.0,
        cart_coords_bohr=np.zeros(3),
        atomic_numbers=np.array([1]),
    )
    with pytest.raises(ValueError, match="does not identify model charge"):
        load_hessian_file(
            path,
            cart_coords_bohr=np.zeros(3),
            atomic_numbers=np.array([1]),
            **LOAD_STATE,
        )
    loaded = load_hessian_file(
        path,
        cart_coords_bohr=np.zeros(3),
        atomic_numbers=np.array([1]),
        allow_unverified_state=True,
        allow_unverified_pes=True,
        **LOAD_STATE,
    )
    assert loaded["schema_version"] == 1
    assert loaded["model_charge"] is None
    assert loaded["model_mult"] is None
    assert loaded["electronic_state_verified"] is False
    assert loaded["potential_identity_verified"] is False


@pytest.mark.parametrize("missing", ["model_charge", "model_mult"])
def test_schema_two_requires_complete_electronic_state_metadata(tmp_path, missing) -> None:
    payload = {
        "schema_version": np.int64(2),
        "hessian": np.eye(3),
        "energy_ha": 0.0,
        "cart_coords_bohr": np.zeros(3),
        "atomic_numbers": np.array([1]),
        "model_charge": np.int64(0),
        "model_mult": np.int64(1),
    }
    payload.pop(missing)
    path = tmp_path / f"missing-{missing}.npz"
    np.savez_compressed(path, **payload)
    with pytest.raises(ValueError, match="lacks electronic-state metadata"):
        load_hessian_file(
            path,
            cart_coords_bohr=np.zeros(3),
            atomic_numbers=np.array([1]),
            allow_unverified_state=True,
            **LOAD_STATE,
        )


@pytest.mark.parametrize(
    ("field", "value", "message"),
    [
        ("model_charge", 0.5, "scalar integer"),
        ("model_mult", 0, "must be >= 1"),
    ],
)
def test_hessian_file_rejects_invalid_state_metadata(tmp_path, field, value, message) -> None:
    kwargs = dict(SAVE_STATE)
    kwargs[field] = value
    with pytest.raises(ValueError, match=message):
        save_hessian_file(
            tmp_path / "invalid.npz",
            hessian=np.eye(3),
            energy_ha=0.0,
            cart_coords_bohr=np.zeros(3),
            atomic_numbers=np.array([1]),
            **kwargs,
        )


@pytest.mark.parametrize("value", [float("nan"), float("inf"), float("-inf")])
def test_hessian_file_rejects_nonfinite_pes_identity(tmp_path, value) -> None:
    identity = {
        **PES_IDENTITY,
        "evaluator": {
            **PES_IDENTITY["evaluator"],
            "potential": {"unexpected_nonfinite": value},
        },
    }
    with pytest.raises(ValueError, match="finite JSON-compatible values"):
        save_hessian_file(
            tmp_path / "nonfinite-identity.npz",
            hessian=np.eye(3),
            energy_ha=0.0,
            cart_coords_bohr=np.zeros(3),
            atomic_numbers=np.array([1]),
            **{**SAVE_STATE, "potential_identity": identity},
        )


@pytest.mark.parametrize("value", [float("nan"), float("inf"), float("-inf")])
def test_hessian_file_rejects_nonfinite_energy_on_save(tmp_path, value) -> None:
    with pytest.raises(ValueError, match="energy_ha must be a finite scalar"):
        save_hessian_file(
            tmp_path / "nonfinite-energy.npz",
            hessian=np.eye(3),
            energy_ha=value,
            cart_coords_bohr=np.zeros(3),
            atomic_numbers=np.array([1]),
            **SAVE_STATE,
        )


@pytest.mark.parametrize("value", [float("nan"), float("inf"), float("-inf")])
def test_hessian_file_rejects_nonfinite_energy_on_load(tmp_path, value) -> None:
    path = tmp_path / "nonfinite-energy.npz"
    np.savez_compressed(
        path,
        schema_version=np.int64(3),
        hessian=np.eye(3),
        energy_ha=value,
        cart_coords_bohr=np.zeros(3),
        atomic_numbers=np.array([1]),
        model_charge=np.int64(0),
        model_mult=np.int64(1),
        potential_identity_json=np.str_(
            '{"evaluator":{"model":"uma-s-1p1","potential":{"mm_backend":'
            '"hessian_ff"},"precision":"float64","backend":"uma"},'
            '"schema":"hessian-cache-identity/v1","system":{"atoms":[6,1,8]}}'
        ),
    )

    with pytest.raises(ValueError, match="energy_ha must be a finite scalar"):
        load_hessian_file(
            path,
            cart_coords_bohr=np.zeros(3),
            atomic_numbers=np.array([1]),
            **LOAD_STATE,
        )


def test_hessian_file_rejects_fractional_schema_version(tmp_path) -> None:
    path = tmp_path / "fractional-schema.npz"
    np.savez_compressed(
        path,
        schema_version=np.float64(3.5),
        hessian=np.eye(3),
        energy_ha=0.0,
        cart_coords_bohr=np.zeros(3),
        atomic_numbers=np.array([1]),
    )

    with pytest.raises(ValueError, match="schema_version must be a scalar integer"):
        load_hessian_file(
            path,
            cart_coords_bohr=np.zeros(3),
            atomic_numbers=np.array([1]),
            **LOAD_STATE,
        )


def test_hessian_file_rejects_pes_identity_mismatch(tmp_path) -> None:
    path = tmp_path / "pes.npz"
    save_hessian_file(
        path,
        hessian=np.eye(3),
        energy_ha=0.0,
        cart_coords_bohr=np.zeros(3),
        atomic_numbers=np.array([1]),
        **SAVE_STATE,
    )

    different = {
        **PES_IDENTITY,
        "evaluator": {
            **PES_IDENTITY["evaluator"],
            "precision": "float32",
        },
    }
    with pytest.raises(ValueError, match="PES identity"):
        load_hessian_file(
            path,
            cart_coords_bohr=np.zeros(3),
            atomic_numbers=np.array([1]),
            expected_model_charge=0,
            expected_model_mult=1,
            expected_potential_identity=different,
        )
