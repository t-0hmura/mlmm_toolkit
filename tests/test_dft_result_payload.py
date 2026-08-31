"""Tests for effective ML/MM DFT result provenance and terminal control."""

from __future__ import annotations

from pathlib import Path

import click
import pytest
from click.testing import CliRunner

from mlmm.core.result_commit import ResultCommitError
from mlmm.backends.pyscf_dft import _subtract_packed_mm_hcore
from mlmm.workflows.dft import (
    _apply_explicit_dft_overrides,
    _build_dft_result_payload,
    _finalize_dft_result,
)


def test_standalone_dft_accepts_pyscf_object_defaults_without_false_conflict(
    tmp_path: Path,
) -> None:
    from mlmm.cli import cli

    config = tmp_path / "pyscf.yaml"
    config.write_text(
        "dft:\n"
        "  pyscf:\n"
        "    mf:\n"
        "      conv_tol: 2.0e-8\n"
        "    grids:\n"
        "      level: 1\n"
        "    mol:\n"
        "      max_memory: 4096\n",
        encoding="utf-8",
    )
    smoke = Path(__file__).resolve().parent / "smoke"

    result = CliRunner().invoke(
        cli,
        [
            "dft", "-v", "3", "-i", str(smoke / "r_complex_layered.pdb"),
            "--parm", str(smoke / "p_complex.parm7"), "-q", "-1", "-m", "1",
            "--engine", "cpu", "--config", str(config), "--show-config",
            "--dry-run", "--out-dir", str(tmp_path / "result"),
        ],
    )

    assert result.exit_code == 0, result.output
    assert "conv_tol: 2.0e-08" in result.output
    assert "grid_level: 1" in result.output
    assert "memory_mb: 4096" in result.output


def test_standalone_dft_honors_explicit_func_basis(tmp_path: Path) -> None:
    from mlmm.cli import cli

    smoke = Path(__file__).resolve().parent / "smoke"
    result = CliRunner().invoke(
        cli,
        [
            "dft",
            "-v",
            "3",
            "-i",
            str(smoke / "r_complex_layered.pdb"),
            "--parm",
            str(smoke / "p_complex.parm7"),
            "-q",
            "-1",
            "-m",
            "1",
            "--func-basis",
            "hf/sto-3g",
            "--engine",
            "cpu",
            "--dry-run",
            "--out-dir",
            str(tmp_path / "result"),
        ],
    )

    assert result.exit_code == 0, result.output
    assert "xc: hf" in result.output
    assert "basis: sto-3g" in result.output


def test_leaf_dft_checkpoint_uses_yaml_effective_output_directory(
    tmp_path: Path,
) -> None:
    from mlmm.core.dft_settings import (
        DFT_CLI_META_KEY,
        finalize_dft_calculator_config,
    )

    ctx = click.Context(click.Command("sp"), info_name="sp")
    ctx.params["out_dir"] = tmp_path / "click-default"
    ctx.meta[DFT_CLI_META_KEY] = {"save_scf_checkpoint": True}
    calc_cfg = {"backend": "dft", "model_charge": 0, "model_mult": 1}
    effective_out = tmp_path / "yaml-output"

    finalize_dft_calculator_config(ctx, calc_cfg, output_dir=effective_out)

    assert calc_cfg["dft_settings"]["checkpoint_path"] == str(
        effective_out / "_work" / "dft_scf" / "state.chk"
    )


@pytest.mark.parametrize(
    "dft_config",
    [
        {"lowmem": "false"},
        {"density_fit": "false"},
        {"save_scf_checkpoint": "false"},
        {"embedcharge": "false"},
        {"pyscf": {"density_fit": {"enabled": "false"}}},
    ],
)
def test_dft_yaml_booleans_require_yaml_boolean_type(dft_config) -> None:
    from mlmm.core.dft_settings import resolve_dft_settings

    with pytest.raises(click.BadParameter, match="must be true or false"):
        resolve_dft_settings({"backend": "dft", "dft": dft_config})


@pytest.mark.parametrize("pbs_var", ["PBS_NP", "PBS_NUM_PPN"])
def test_dft_resources_honor_pbs_cpu_counts(monkeypatch, pbs_var) -> None:
    from mlmm.core import dft_settings

    for name in (
        "OMP_NUM_THREADS",
        "SLURM_CPUS_PER_TASK",
        "NSLOTS",
        "PBS_NP",
        "PBS_NUM_PPN",
        "PBS_NODEFILE",
    ):
        monkeypatch.delenv(name, raising=False)
    monkeypatch.setenv(pbs_var, "3")
    monkeypatch.setattr(dft_settings, "_affinity_count", lambda: None)

    assert dft_settings.resolve_dft_settings({"backend": "dft"}).nprocs == 3


def test_dft_resources_honor_pbs_nodefile(monkeypatch, tmp_path: Path) -> None:
    from mlmm.core import dft_settings

    for name in (
        "OMP_NUM_THREADS",
        "SLURM_CPUS_PER_TASK",
        "NSLOTS",
        "PBS_NP",
        "PBS_NUM_PPN",
        "PBS_NODEFILE",
    ):
        monkeypatch.delenv(name, raising=False)
    nodefile = tmp_path / "pbs_nodes"
    nodefile.write_text("node02\nnode02\nnode03\n", encoding="utf-8")
    monkeypatch.setenv("PBS_NODEFILE", str(nodefile))
    monkeypatch.setattr(dft_settings, "_affinity_count", lambda: None)

    assert dft_settings.resolve_dft_settings({"backend": "dft"}).nprocs == 3


def test_gpu_lowmem_embedding_preserves_packed_hcore_layout_and_sign() -> None:
    import numpy as np

    class DeviceArray:
        def __init__(self, value):
            self.value = value

        def get(self):
            return self.value

    def pack_lower_triangle(square):
        return DeviceArray(square[np.tril_indices(square.shape[0])])

    hcore = np.array([10.0, 20.0, 30.0])
    mm_potential = np.array([[1.0, 2.0], [2.0, 4.0]])

    combined = _subtract_packed_mm_hcore(
        hcore, mm_potential, pack_lower_triangle
    )

    assert np.array_equal(combined, np.array([9.0, 18.0, 26.0]))


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
        dft_kw={
            "grid_level": 7,
            "conv_tol": 2.0e-11,
            "max_cycle": 17,
            "lowmem": True,
            "memory_mode": "direct_jk",
            "nprocs": 8,
            "nprocs_source": "explicit",
            "memory_mb": 64000,
            "memory_source": "explicit",
        },
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
    assert payload["lowmem_requested"] is True
    assert payload["dft_resources"]["memory_mode"] == "direct_jk"
    assert payload["dft_resources"]["memory_mb"] == 64000
    assert payload["mlip_backend"] == "dft"
    assert payload["mlip_model"] is None
    assert payload["mlip_model_label"] is None
    assert payload["mlip_task"] is None
    assert payload["mlip_precision"] is None


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
