"""Tests for the plain ``.npy`` Hessian files of ``--dump-hess`` / ``--read-hess``."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest
from click.testing import CliRunner

from mlmm.core import result_commit
from mlmm.core.result_commit import ResultCommitError
from mlmm.io.hessian_file import load_hessian_file, save_hessian_file
from mlmm.workflows.freq import _record_hessian_result_path

SMOKE = Path(__file__).resolve().parents[1] / "tests" / "smoke"


def test_full_hessian_round_trips_as_a_plain_npy_array(tmp_path) -> None:
    path = tmp_path / "h.npy"
    hess = np.diag(np.arange(1.0, 10.0))
    save_hessian_file(path, hess)

    np.testing.assert_array_equal(np.load(path), hess)  # readable by any NumPy user
    np.testing.assert_array_equal(load_hessian_file(path, n_atoms=3), hess)


def test_full_hessian_is_restricted_to_the_movable_atoms(tmp_path) -> None:
    path = tmp_path / "full.npy"
    save_hessian_file(path, np.diag(np.arange(1.0, 10.0)))

    loaded = load_hessian_file(path, n_atoms=3, active_dofs=[6, 7, 8])

    np.testing.assert_array_equal(loaded, np.diag([7.0, 8.0, 9.0]))


def test_movable_atom_block_is_used_as_is(tmp_path) -> None:
    path = tmp_path / "partial.npy"
    np.save(path, np.eye(6))

    loaded = load_hessian_file(path, n_atoms=3, active_dofs=[0, 1, 2, 6, 7, 8])

    np.testing.assert_array_equal(loaded, np.eye(6))


@pytest.mark.parametrize(
    ("active", "message"),
    [
        (None, r"is 6x6; expected 9x9 \(all atoms\)\.$"),
        ([0, 1, 2], r"is 6x6; expected 9x9 \(all atoms\) or 3x3 \(movable atoms only\)"),
    ],
)
def test_wrong_size_is_rejected(tmp_path, active, message) -> None:
    path = tmp_path / "h.npy"
    np.save(path, np.eye(6))
    with pytest.raises(ValueError, match=message):
        load_hessian_file(path, n_atoms=3, active_dofs=active)


@pytest.mark.parametrize(
    ("array", "message"),
    [
        (np.triu(np.ones((3, 3))), "is not symmetric"),
        (np.full((3, 3), np.nan), "non-finite"),
        (np.zeros(9), "real square matrix"),
        (np.eye(3, dtype=complex), "real square matrix"),
    ],
)
def test_invalid_matrix_is_rejected(tmp_path, array, message) -> None:
    path = tmp_path / "h.npy"
    np.save(path, array)
    with pytest.raises(ValueError, match=message):
        load_hessian_file(path, n_atoms=1)


def test_npz_archive_is_rejected(tmp_path) -> None:
    path = tmp_path / "h.npz"
    np.savez(path, hessian=np.eye(3))
    with pytest.raises(ValueError, match="is not a .npy array"):
        load_hessian_file(path, n_atoms=1)


@pytest.mark.parametrize("content", [b"", b"not numpy\n", b"PK\x03\x04broken"])
def test_unreadable_file_is_rejected(tmp_path, content) -> None:
    path = tmp_path / "h.npy"
    path.write_bytes(content)
    with pytest.raises(ValueError, match="Cannot read"):
        load_hessian_file(path, n_atoms=1)


@pytest.mark.parametrize("array", [np.ones((2, 3)), np.array([[np.inf]])])
def test_save_rejects_a_non_square_or_non_finite_matrix(tmp_path, array) -> None:
    with pytest.raises(ValueError, match="finite square matrix"):
        save_hessian_file(tmp_path / "h.npy", array)


@pytest.mark.parametrize("name", ["hessian", "hessian.bin", "hessian.npy"])
def test_save_uses_the_exact_requested_path(tmp_path, name: str) -> None:
    requested = tmp_path / name

    assert save_hessian_file(requested, np.eye(3)) == requested
    assert [path.name for path in tmp_path.iterdir()] == [name]
    np.testing.assert_array_equal(load_hessian_file(requested, n_atoms=1), np.eye(3))

    files = _record_hessian_result_path({"frequencies_txt": "frequencies_cm-1.txt"}, requested)
    assert files["hessian_npy"] == str(requested.resolve())


def test_publish_failure_preserves_the_old_file(tmp_path, monkeypatch) -> None:
    requested = tmp_path / "hessian.npy"
    requested.write_bytes(b"old-hessian")

    def fail_replace(source, destination):
        raise OSError("injected")

    monkeypatch.setattr(result_commit.os, "replace", fail_replace)
    with pytest.raises(ResultCommitError, match="publish"):
        save_hessian_file(requested, np.eye(3))
    assert requested.read_bytes() == b"old-hessian"
    assert [path.name for path in tmp_path.iterdir()] == [requested.name]


def _smoke_args(out: Path) -> list[str]:
    return [
        "-i", str(SMOKE / "p_complex_layered.pdb"),
        "--parm", str(SMOKE / "p_complex.parm7"),
        "-q", "-1", "-m", "1",
        "--out-dir", str(out),
    ]


@pytest.mark.parametrize(
    ("command", "option", "message"),
    [
        ("freq", "--read-hess", "collides with a reserved frequency output"),
        ("freq", "--dump-hess", "collides with a reserved frequency output"),
        ("tsopt", "--dump-hess", "collides with a reserved TSOPT output"),
    ],
)
def test_hessian_file_at_a_reserved_output_path_is_rejected(
    tmp_path, command, option, message
) -> None:
    import importlib

    module = importlib.import_module(f"mlmm.workflows.{command}")
    out = tmp_path / "out"
    out.mkdir()
    reserved = out / "result.json"
    reserved.write_bytes(b"existing")

    result = CliRunner().invoke(module.cli, _smoke_args(out) + [option, str(reserved)])

    assert result.exit_code == 2, result.output
    assert message in result.output
    assert reserved.read_bytes() == b"existing"


@pytest.mark.parametrize("dry_run", [[], ["--dry-run"]])
@pytest.mark.parametrize(
    ("command", "section", "extra"),
    [("irc", "irc", []), ("tsopt", "rsirfo", ["--opt-mode", "hess"])],
)
def test_read_hess_needs_the_calc_hessian_init(
    tmp_path, command, section, extra, dry_run
) -> None:
    import importlib

    module = importlib.import_module(f"mlmm.workflows.{command}")
    config = tmp_path / "config.yaml"
    config.write_text(f"{section}:\n  hessian_init: unit\n", encoding="utf-8")
    hess = tmp_path / "h.npy"
    np.save(hess, np.eye(3))

    result = CliRunner().invoke(
        module.cli,
        _smoke_args(tmp_path / "out")
        + extra
        + dry_run
        + ["--config", str(config), "--read-hess", str(hess)],
    )

    assert result.exit_code == 2, result.output
    assert "--read-hess needs hessian_init: calc." in result.output


def test_tsopt_final_hessian_round_trips_through_dump_and_read(
    tmp_path, capsys, monkeypatch
) -> None:
    from types import SimpleNamespace

    import torch

    from mlmm.core.result_commit import MLMM_RUN_ID_ENV
    from mlmm.io import hessian_cache
    from mlmm.workflows import tsopt

    monkeypatch.delenv(MLMM_RUN_ID_ENV, raising=False)  # a standalone tsopt run
    geom = SimpleNamespace(
        atomic_numbers=np.array([6, 1, 8]),
        atoms=["C", "H", "O"],
        cart_coords=np.arange(9, dtype=float) / 10.0,
        freeze_atoms=[1],
    )
    calc = {"backend": "uma", "model_charge": 0, "model_mult": 1, "freeze_atoms": [1]}
    block = np.diag(np.arange(1.0, 7.0))
    hessian_cache.clear()
    try:
        hessian_cache.store(
            "ts",
            torch.as_tensor(block),
            active_dofs=[0, 1, 2, 6, 7, 8],
            meta={"cart_coords": geom.cart_coords},
            identity=hessian_cache.identity_from_context(geom, calc, role="ts"),
        )
        written = tsopt._dump_terminal_hessian(tmp_path / "ts.npy", geom, calc)
        np.testing.assert_array_equal(np.load(written), block)
        loaded = tsopt._load_initial_hessian_file(written, geom, calc)
        np.testing.assert_array_equal(loaded["hessian"], block)
        assert loaded["active_dofs"] == [0, 1, 2, 6, 7, 8]

        geom.cart_coords = geom.cart_coords + 0.01
        assert tsopt._dump_terminal_hessian(tmp_path / "moved.npy", geom, calc) is None
        assert "--dump-hess file was not written" in capsys.readouterr().err
        assert not (tmp_path / "moved.npy").exists()
    finally:
        hessian_cache.clear()
