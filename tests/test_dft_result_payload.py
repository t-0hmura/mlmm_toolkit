"""Tests for effective ML/MM DFT result provenance and terminal control."""

from __future__ import annotations

from pathlib import Path

import pytest

from mlmm.core.result_commit import ResultCommitError
from mlmm.workflows.dft import (
    _apply_explicit_dft_overrides,
    _build_dft_result_payload,
    _finalize_dft_result,
)


def _payload(*, converged: bool, engine: str = "pyscf(cpu)"):
    return _build_dft_result_payload(
        converged=converged,
        energy_hartree=-10.0,
        energy_kcal_per_mol=-6275.0,
        xc="wb97m-v",
        basis="def2-tzvpd",
        engine_label=engine,
        using_gpu="gpu4pyscf" in engine,
        using_lowmem="lowmem" in engine,
        dft_kw={"grid_level": 7, "conv_tol": 2.0e-11, "max_cycle": 17},
        calc_kw={
            "backend": "orb",
            "orb_model": "orb-v3-conservative-inf-omat",
            "orb_precision": "fp64",
            "model_charge": -1,
            "model_mult": 2,
            "mm_backend": "hessian_ff",
            "link_atom_method": "fixed",
            "use_cmap": False,
        },
        n_atoms=12,
        input_path=Path("relative/input.pdb"),
        charges={"mulliken": [0.1]},
        spin_densities={"mulliken": [0.2]},
    )


def test_yaml_only_dft_values_survive_when_cli_is_omitted() -> None:
    yaml_values = {"grid_level": 7, "conv_tol": 2.0e-11, "max_cycle": 17}
    resolved = _apply_explicit_dft_overrides(
        yaml_values,
        is_param_explicit=lambda name: False,
        conv_tol=1.0e-8,
        max_cycle=100,
        grid_level=3,
        out_dir=Path("default"),
        lowmem=False,
    )
    assert resolved == yaml_values
    assert resolved is not yaml_values


def test_explicit_cli_dft_values_override_yaml() -> None:
    explicit = {"conv_tol", "max_cycle", "grid_level", "out_dir", "lowmem"}
    resolved = _apply_explicit_dft_overrides(
        {"grid_level": 7, "conv_tol": 2.0e-11, "max_cycle": 17},
        is_param_explicit=lambda name: name in explicit,
        conv_tol=4.0e-9,
        max_cycle=23,
        grid_level=5,
        out_dir=Path("explicit"),
        lowmem=True,
    )
    assert resolved["conv_tol"] == pytest.approx(4.0e-9)
    assert resolved["max_cycle"] == 23
    assert resolved["grid_level"] == 5
    assert resolved["out_dir"] == "explicit"
    assert resolved["lowmem"] is True


@pytest.mark.parametrize(
    ("converged", "status"),
    [(True, "converged"), (False, "not_converged")],
)
@pytest.mark.parametrize(
    "engine",
    ["pyscf(cpu)", "gpu4pyscf", "gpu4pyscf(rks_lowmem)"],
)
def test_payload_preserves_legacy_keys_and_records_effective_values(
    converged: bool, status: str, engine: str
) -> None:
    payload = _payload(converged=converged, engine=engine)
    legacy_keys = {
        "converged",
        "energy_hartree",
        "energy_kcal_per_mol",
        "xc_functional",
        "basis_set",
        "engine",
        "used_gpu",
        "used_lowmem",
        "mlip_backend",
        "mlip_model",
        "mlip_precision",
        "mm_backend",
        "link_atom_method",
        "use_cmap",
        "charge",
        "spin",
        "n_atoms",
        "grid_level",
        "conv_tol",
        "input_file",
        "charges",
        "spin_densities",
        "files",
    }
    assert legacy_keys <= payload.keys()
    assert payload["status"] == status
    assert payload["converged"] is converged
    assert payload["engine"] == engine
    assert payload["grid_level"] == 7
    assert payload["conv_tol"] == pytest.approx(2.0e-11)
    assert payload["max_cycle"] == 17


def test_nonconverged_payload_commits_before_exit_three(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    events = []

    def fake_write(out_dir, payload, **kwargs):
        events.append(("write", payload["status"], Path(out_dir)))
        return Path(out_dir) / "result.json"

    monkeypatch.setattr("mlmm.core.utils.write_result_json", fake_write)
    with pytest.raises(SystemExit) as exc_info:
        _finalize_dft_result(
            out_json=True,
            out_dir=tmp_path,
            payload=_payload(converged=False),
            elapsed_seconds=1.0,
        )
    assert exc_info.value.code == 3
    assert events == [("write", "not_converged", tmp_path)]


def test_result_write_failure_is_not_reclassified_as_scf_nonconvergence(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    def fail_write(*args, **kwargs):
        raise ResultCommitError("publish", tmp_path / "result.json", OSError("injected"))

    monkeypatch.setattr("mlmm.core.utils.write_result_json", fail_write)
    with pytest.raises(ResultCommitError, match="publish"):
        _finalize_dft_result(
            out_json=True,
            out_dir=tmp_path,
            payload=_payload(converged=False),
            elapsed_seconds=1.0,
        )
