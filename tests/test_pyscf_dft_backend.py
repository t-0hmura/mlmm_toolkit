"""Small CPU contracts for the native PySCF ML/MM high-level backend."""

from __future__ import annotations

import json
import shutil
import numpy as np
import pytest
import torch
from ase import Atoms


pytest.importorskip("pyscf")


def test_dft_rejects_redundant_backend_selector() -> None:
    from click.testing import CliRunner
    from mlmm.cli import cli

    result = CliRunner().invoke(cli, ["dft", "-b", "dft"])

    assert result.exit_code == 2
    assert "already selects the DFT evaluator" in result.output
    assert "sp -b dft" in result.output


def test_calculator_leaf_help_describes_native_dft_embedding() -> None:
    from click.testing import CliRunner
    from mlmm.cli import cli

    result = CliRunner().invoke(cli, ["sp", "--help-advanced"])

    assert result.exit_code == 0, result.output
    assert "xTB point-charge delta" in result.output
    assert "native PySCF MM point charges" in result.output


def test_dft_resource_defaults_are_lowmem_and_explicit_values_normalize(
    monkeypatch,
) -> None:
    from mlmm.core.dft_settings import resolve_dft_settings

    monkeypatch.setenv("OMP_NUM_THREADS", "6")
    settings = resolve_dft_settings({"backend": "dft"})
    explicit = resolve_dft_settings({
        "backend": "dft",
        "dft": {"nprocs": 3, "memory": "64GB"},
    })

    assert settings.func_basis == "wb97m-v/def2-svp"
    assert settings.lowmem is True
    assert settings.density_fit is False
    assert settings.nprocs <= 6
    assert explicit.nprocs == 3
    assert explicit.memory_mb == 64000


def test_dft_resources_are_provenance_but_not_checkpoint_identity() -> None:
    from mlmm.core.dft_settings import resolve_dft_settings
    from mlmm.core.utils import calculator_provenance

    small = resolve_dft_settings({
        "backend": "dft", "dft": {"nprocs": 2, "memory": "4GB"}
    })
    large = resolve_dft_settings({
        "backend": "dft", "dft": {"nprocs": 16, "memory": "64GB"}
    })
    provenance = calculator_provenance({
        "backend": "dft", "dft_settings": large.to_dict()
    })

    assert small.scientific_identity() == large.scientific_identity()
    assert provenance["dft_resources"] == {
        "memory_mode": "gpu4pyscf_rks_lowmem",
        "nprocs": 16,
        "nprocs_source": "explicit",
        "memory_mb": 64000,
        "memory_source": "explicit",
    }


@pytest.mark.parametrize(
    ("dft", "multiplicity", "expected"),
    [
        ({}, 1, "gpu4pyscf_rks_lowmem"),
        ({"engine": "cpu"}, 1, "direct_jk"),
        ({}, 2, "direct_jk"),
        ({"lowmem": False}, 1, "density_fit"),
    ],
)
def test_dft_memory_mode_tracks_the_effective_driver(
    dft, multiplicity, expected
) -> None:
    from mlmm.core.dft_settings import resolve_dft_settings

    calc = {"backend": "dft", "model_mult": multiplicity, "dft": dft}
    assert resolve_dft_settings(calc).memory_mode == expected


@pytest.mark.parametrize(
    "key", ["charge", "multiplicity", "embedcharge", "embedcharge_cutoff"]
)
def test_calculator_dft_mapping_rejects_top_level_owned_state(key) -> None:
    import click

    from mlmm.core.dft_settings import resolve_dft_settings

    with pytest.raises(click.BadParameter, match=f"calc.dft.{key}"):
        resolve_dft_settings({"backend": "dft", "dft": {key: 1}})


def test_dft_yaml_nprocs_reports_a_click_validation_error() -> None:
    import click

    from mlmm.core.dft_settings import resolve_dft_settings

    with pytest.raises(click.BadParameter, match="positive integer"):
        resolve_dft_settings({"backend": "dft", "dft": {"nprocs": "many"}})


@pytest.mark.parametrize(
    "pyscf_config, field",
    [
        ({"mol": {"basis": "sto-3g"}}, "basis"),
        ({"mf": {"xc": "pbe"}}, "xc"),
    ],
)
def test_dft_rejects_resolver_owned_pyscf_method_fields(
    pyscf_config, field
) -> None:
    import click

    from mlmm.core.dft_settings import resolve_dft_settings

    with pytest.raises(click.BadParameter, match=field):
        resolve_dft_settings(
            {"backend": "dft", "dft": {"pyscf": pyscf_config}}
        )


@pytest.mark.parametrize("conv_tol", [float("nan"), float("inf"), float("-inf")])
def test_dft_rejects_nonfinite_convergence_tolerance(conv_tol) -> None:
    import click

    from mlmm.core.dft_settings import resolve_dft_settings

    with pytest.raises(click.BadParameter, match="conv_tol"):
        resolve_dft_settings(
            {"backend": "dft", "dft": {"conv_tol": conv_tol}}
        )


def test_mlmm_factory_binds_canonical_charge_and_multiplicity() -> None:
    from mlmm.backends.mlmm_calc import _create_ml_backend

    backend = _create_ml_backend(
        "dft",
        dft_settings=_settings(charge=0, multiplicity=1),
        model_charge=-1,
        model_mult=2,
        ml_device=torch.device("cpu"),
    )

    assert backend.settings.charge == -1
    assert backend.settings.multiplicity == 2


def test_dft_scientific_identity_tracks_active_embedding_cutoff() -> None:
    from mlmm.core.dft_settings import resolve_dft_settings

    disabled = resolve_dft_settings(
        {"backend": "dft", "embedcharge": False, "embedcharge_cutoff": 4.0}
    )
    enabled = resolve_dft_settings(
        {"backend": "dft", "embedcharge": True, "embedcharge_cutoff": 4.0}
    )

    assert disabled.scientific_identity()["embedcharge_cutoff"] is None
    assert enabled.scientific_identity()["embedcharge_cutoff"] == 4.0


@pytest.mark.parametrize(
    "calc_cfg, message",
    [
        ({"model_charge": 0.5}, "model_charge"),
        ({"model_mult": 0}, "model_mult"),
        ({"embedcharge_cutoff": float("nan")}, "embedcharge_cutoff"),
        ({"dft": {"nprocs": 1.5}}, "DFT nprocs"),
    ],
)
def test_dft_rejects_invalid_canonical_charge_spin_and_cutoff(
    calc_cfg, message
) -> None:
    import click

    from mlmm.core.dft_settings import resolve_dft_settings

    with pytest.raises(click.BadParameter, match=message):
        resolve_dft_settings({"backend": "dft", **calc_cfg})


@pytest.mark.parametrize(
    "key",
    [
        "typo_setting",
        "charge",
        "multiplicity",
        "save_scf_checkpoint",
        "checkpoint_path",
        "nprocs_source",
        "memory_source",
        "embedcharge",
        "embedcharge_cutoff",
    ],
)
def test_standalone_dft_mapping_fails_closed_on_unowned_key(key) -> None:
    import click

    from mlmm.core.dft_settings import standalone_dft_settings_mapping

    with pytest.raises(click.BadParameter, match=f"dft.{key}"):
        standalone_dft_settings_mapping({"out_dir": "result", key: 1})


def _settings(**overrides):
    values = {
        "func_basis": "hf/sto-3g",
        "engine": "cpu",
        "charge": 0,
        "multiplicity": 1,
        "density_fit": False,
        "save_scf_checkpoint": False,
    }
    values.update(overrides)
    return values


def test_he_cache_and_analytical_hessian() -> None:
    from mlmm.backends.pyscf_dft import create_dft_backend

    backend = create_dft_backend(_settings())
    atoms = Atoms("He", positions=[[0.0, 0.0, 0.0]])
    energy, forces, _ = backend.eval(atoms, need_grad=True)
    backend.eval(atoms, need_grad=True)
    hessian = backend.hessian_analytical(
        atoms, 1, dtype=torch.float64
    )

    assert energy == pytest.approx(-2.80778395754 * 27.211386245988, abs=1.0e-8)
    assert np.linalg.norm(forces) < 1.0e-9
    assert hessian.shape == (1, 3, 1, 3)
    assert len(backend.session.metrics) == 1
    assert backend.session._last_good is None


def test_same_coordinate_energy_then_force_runs_one_scf() -> None:
    from mlmm.backends.pyscf_dft import create_dft_backend

    backend = create_dft_backend(_settings())
    atoms = Atoms("He", positions=[[0.0, 0.0, 0.0]])

    backend.energy(atoms)
    backend.eval(atoms, need_grad=True)

    assert len(backend.session.metrics) == 1


@pytest.mark.parametrize("from_checkpoint", [False, True])
def test_gpu_lowmem_passes_reused_density_to_rebuilt_method(
    monkeypatch, from_checkpoint
) -> None:
    from mlmm.backends.pyscf_dft import PySCFDFTSession
    from mlmm.core.dft_settings import resolve_dft_settings

    settings = resolve_dft_settings({
        "backend": "dft",
        "dft": {"func_basis": "lda/sto-3g", "engine": "gpu"},
    })
    density = object()

    class FakeMethod:
        converged = True

        def __init__(self):
            self.dm0 = None

        def make_rdm1(self):
            return density

        def kernel(self, dm0=None):
            self.dm0 = dm0
            return -1.0

    session = PySCFDFTSession(settings)
    rebuilt = FakeMethod()
    if from_checkpoint:
        session._pending_checkpoint = True
        session._last_good = {"loaded": True}
    else:
        session._scanner = FakeMethod()
        session._using_rks_lowmem = True

    def build_scanner(mol):
        session._scanner = rebuilt
        session._using_rks_lowmem = True

    monkeypatch.setattr(session, "_build_scanner", build_scanner)
    monkeypatch.setattr(session, "_update_mm_mol", lambda: None)

    _, guess_source = session._run_scf(object())

    assert rebuilt.dm0 is density
    assert guess_source == (
        "checkpoint" if from_checkpoint else "previous_density"
    )


def test_native_embedding_returns_qm_and_mm_site_forces() -> None:
    from mlmm.backends.pyscf_dft import create_dft_backend

    backend = create_dft_backend(
        _settings(func_basis="lda/sto-3g", embedcharge=True)
    )
    atoms = Atoms("H2", positions=[[0.0, 0.0, -0.7], [0.0, 0.0, 0.7]])
    backend.set_embedding(
        np.array([[2.0, 0.0, 0.0]]), np.array([0.2]), [7]
    )
    energy, qm_forces, _ = backend.eval(atoms, need_grad=True)
    indices, mm_forces = backend.mm_forces()

    assert np.isfinite(energy)
    assert qm_forces.shape == (2, 3)
    assert indices == [7]
    assert mm_forces is not None and mm_forces.shape == (1, 3)
    assert np.allclose(qm_forces.sum(axis=0) + mm_forces.sum(axis=0), 0.0, atol=1.0e-7)


def test_open_shell_embedding_uses_the_total_spin_density() -> None:
    from mlmm.backends.pyscf_dft import create_dft_backend

    backend = create_dft_backend(
        _settings(func_basis="hf/sto-3g", multiplicity=2, embedcharge=True)
    )
    atoms = Atoms("H", positions=[[0.0, 0.0, 0.0]])
    backend.set_embedding(
        np.array([[2.0, 0.0, 0.0]]), np.array([0.2]), [4]
    )
    energy, qm_forces, _ = backend.eval(atoms, need_grad=True)
    indices, mm_forces = backend.mm_forces()

    assert np.isfinite(energy)
    assert indices == [4]
    assert mm_forces is not None and mm_forces.shape == (1, 3)
    assert np.allclose(qm_forces.sum(axis=0) + mm_forces.sum(axis=0), 0.0, atol=1.0e-7)


def test_checkpoint_round_trip(tmp_path) -> None:
    from mlmm.backends.pyscf_dft import create_dft_backend

    checkpoint = tmp_path / "state.chk"
    settings = _settings(
        save_scf_checkpoint=True, checkpoint_path=str(checkpoint)
    )
    atoms = Atoms("He", positions=[[0.0, 0.0, 0.0]])
    first = create_dft_backend(settings)
    first.eval(atoms, need_grad=True)
    first.save_scf_checkpoint(checkpoint, atoms)
    restored = create_dft_backend(settings)
    assert restored.load_scf_checkpoint(checkpoint, atoms)
    restored.eval(atoms, need_grad=True)

    assert checkpoint.is_file()
    assert checkpoint.with_suffix(".chk.json").is_file()
    assert restored.session.metrics[0]["guess_source"] == "checkpoint"

    mismatched = create_dft_backend(settings)
    shifted = Atoms("He", positions=[[0.1, 0.0, 0.0]])
    assert not mismatched.load_scf_checkpoint(checkpoint, shifted)
    assert mismatched.session.checkpoint_status["reason"] == "coordinates_angstrom_mismatch"


def test_last_good_checkpoint_survives_lost_scanner(tmp_path) -> None:
    from mlmm.backends.pyscf_dft import create_dft_backend

    checkpoint = tmp_path / "last-good.chk"
    backend = create_dft_backend(_settings(save_scf_checkpoint=True))
    atoms = Atoms("He", positions=[[0.0, 0.0, 0.0]])
    backend.energy(atoms)
    backend.session._scanner = None

    backend.save_scf_checkpoint(checkpoint, atoms)

    assert checkpoint.is_file()
    assert checkpoint.with_suffix(".chk.json").is_file()


def test_checkpoint_rejects_mixed_binary_and_metadata_generations(tmp_path) -> None:
    from mlmm.backends.pyscf_dft import create_dft_backend

    atoms = Atoms("He", positions=[[0.0, 0.0, 0.0]])
    first_path = tmp_path / "first.chk"
    second_path = tmp_path / "second.chk"
    for path in (first_path, second_path):
        backend = create_dft_backend(_settings(save_scf_checkpoint=True))
        backend.energy(atoms)
        backend.save_scf_checkpoint(path, atoms)
        backend.close()

    shutil.copyfile(
        second_path.with_suffix(".chk.json"),
        first_path.with_suffix(".chk.json"),
    )
    metadata_path = first_path.with_suffix(".chk.json")
    metadata = json.loads(metadata_path.read_text(encoding="utf-8"))
    metadata["schema"] = 1
    metadata_path.write_text(json.dumps(metadata), encoding="utf-8")
    restored = create_dft_backend(_settings(save_scf_checkpoint=True))

    assert not restored.load_scf_checkpoint(first_path, atoms)
    assert restored.session.checkpoint_status["reason"] == "schema_mismatch"
    metadata["schema"] = 2
    metadata_path.write_text(json.dumps(metadata), encoding="utf-8")
    assert not restored.load_scf_checkpoint(first_path, atoms)
    assert restored.session.checkpoint_status["reason"] == "generation_mismatch"


def test_settings_forward_pyscf_objects_and_reject_analytical_embedding() -> None:
    import click
    from click.testing import CliRunner

    from mlmm.core.dft_settings import (
        finalize_dft_calculator_config,
        resolve_dft_settings,
    )

    settings = resolve_dft_settings(
        {
            "backend": "dft",
            "dft": {
                "func_basis": "hf/sto-3g",
                "engine": "cpu",
                "pyscf": {
                    "mf": {"max_cycle": 42},
                },
            },
        }
    )
    assert settings.max_cycle == 42

    @click.command()
    @click.pass_context
    def command(ctx):
        finalize_dft_calculator_config(
            ctx,
            {
                "backend": "dft",
                "embedcharge": True,
                "hessian_calc_mode": "Analytical",
                "dft": {"func_basis": "hf/sto-3g", "engine": "cpu"},
            },
        )

    result = CliRunner().invoke(command)
    assert result.exit_code == 2
    assert "complete QM-MM response" in result.output


def test_leaf_checkpoint_default_is_output_local(tmp_path) -> None:
    import click

    from mlmm.core.dft_settings import (
        DFT_CLI_META_KEY,
        finalize_dft_calculator_config,
    )

    ctx = click.Context(click.Command("sp"), info_name="sp")
    ctx.params["out_dir"] = tmp_path / "result"
    ctx.meta[DFT_CLI_META_KEY] = {"save_scf_checkpoint": True}
    calc_cfg = {"backend": "dft", "model_charge": 0, "model_mult": 1}

    finalize_dft_calculator_config(ctx, calc_cfg)

    assert calc_cfg["dft_settings"]["checkpoint_path"] == str(
        tmp_path / "result" / "_work" / "dft_scf" / "state.chk"
    )


@pytest.mark.parametrize(
    "dft_config, message",
    [
        ({"solvent": "water"}, "calc.dft.solvent"),
        ({"solvent_model": "pcm"}, "calc.dft.solvent_model"),
        ({"pyscf": {"with_solvent": {}}}, "with_solvent"),
    ],
)
def test_implicit_solvent_settings_are_rejected(
    dft_config, message
) -> None:
    import click

    from mlmm.core.dft_settings import resolve_dft_settings

    with pytest.raises(click.BadParameter, match=message):
        resolve_dft_settings(
            {
                "backend": "dft",
                "dft": dft_config,
            }
        )


def test_embedding_charge_indices_are_fixed_for_the_pes_lifetime() -> None:
    from types import SimpleNamespace

    from mlmm.backends.mlmm_calc import MLMMCore

    recorded = []
    core = MLMMCore.__new__(MLMMCore)
    core.backend_name = "dft"
    core.embedcharge = True
    core.embedcharge_cutoff = 1.5
    core.selection_indices = [0]
    core._dft_embedding_indices = None
    core._ml_backend = SimpleNamespace(
        set_embedding=lambda coords, charges, indices: recorded.append(list(indices))
    )
    core._get_mm_charges = lambda indices: np.zeros(len(indices))
    core.print_timing = False

    initial = Atoms("HHH", positions=[[0, 0, 0], [1, 0, 0], [3, 0, 0]])
    moved = Atoms("HHH", positions=[[0, 0, 0], [3, 0, 0], [1, 0, 0]])
    core._configure_native_dft_embedding(initial)
    core._configure_native_dft_embedding(moved)

    assert recorded == [[1], [1]]
