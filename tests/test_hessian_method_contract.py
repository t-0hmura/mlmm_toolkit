"""Strict ML/MM constructor and YAML method-enum contracts."""

from __future__ import annotations

from pathlib import Path

import pytest
import yaml
from click.testing import CliRunner

from mlmm.backends import normalize_calculator_methods
from mlmm.backends.mlmm_calc import (
    MLMMCore,
    _gather_atom_hessian_square,
    normalize_hessian_calc_mode,
    normalize_link_atom_method,
    normalize_mm_hessian_mode,
)
from mlmm.cli import cli as root_cli


def test_atom_hessian_square_gather_matches_chained_reference() -> None:
    import torch

    generator = torch.Generator().manual_seed(123)
    hessian = torch.randn(
        8, 3, 8, 3, generator=generator, dtype=torch.float64
    )
    indices = torch.tensor([6, 1, 4], dtype=torch.long)

    expected = hessian.index_select(0, indices).index_select(2, indices)
    actual = _gather_atom_hessian_square(hessian, indices)

    assert torch.equal(actual, expected)


@pytest.mark.parametrize(
    ("raw", "expected"),
    [
        ("FiniteDifference", "FiniteDifference"),
        ("finitedifference", "FiniteDifference"),
        ("ANALYTICAL", "Analytical"),
        (None, "FiniteDifference"),
    ],
)
def test_hessian_mode_normalizes_only_the_supported_vocabulary(raw, expected):
    assert normalize_hessian_calc_mode(raw) == expected


@pytest.mark.parametrize(
    ("raw", "expected"),
    [("SCALED", "scaled"), ("fixed", "fixed"), (None, "scaled")],
)
def test_link_atom_method_normalizes_only_the_supported_vocabulary(raw, expected):
    assert normalize_link_atom_method(raw) == expected


@pytest.mark.parametrize(
    ("raw", "mm_fd", "expected"),
    [
        (None, True, "finite_difference"),
        (None, False, "analytical"),
        ("finite-difference", False, "finite_difference"),
        ("analytical", True, "analytical"),
        ("none", True, "none"),
    ],
)
def test_mm_hessian_mode_separates_method_from_high_only_policy(
    raw, mm_fd, expected
):
    assert normalize_mm_hessian_mode(raw, mm_fd=mm_fd) == expected


def test_hessianff_finite_difference_dispatches_numerical_engine(
    monkeypatch,
) -> None:
    import numpy as np
    from ase import Atoms

    from mlmm.backends.mlmm_calc import hessianffCalculator
    from mlmm.io import hessian_calc as hessian_calc_module

    calls = []

    def fake_hessian_calc(atoms, calculator, **kwargs):
        calls.append((atoms, calculator, kwargs))
        return 2.0 * np.eye(3 * len(atoms))

    monkeypatch.setattr(
        hessian_calc_module,
        "hessian_calc",
        fake_hessian_calc,
    )
    calculator = hessianffCalculator.__new__(hessianffCalculator)
    atoms = Atoms("H", positions=[[0.0, 0.0, 0.0]])

    hessian, active = calculator.finite_difference_hessian(
        atoms,
        delta=0.002,
        return_partial_hessian=False,
    )

    np.testing.assert_allclose(hessian, 2.0 * np.eye(3))
    assert active is None
    assert calls[0][1] is calculator
    assert calls[0][2]["delta"] == pytest.approx(0.002)


@pytest.mark.parametrize(
    ("keyword", "value", "message"),
    [
        ("hessian_calc_mode", "finite-difference", "hessian_calc_mode"),
        ("link_atom_method", "legacy", "link_atom_method"),
    ],
)
def test_direct_core_rejects_invalid_method_before_workspace_io(
    keyword: str, value: str, message: str, monkeypatch
) -> None:
    def workspace_must_not_be_created(*_args, **_kwargs):
        raise AssertionError("temporary workspace was created before enum validation")

    monkeypatch.setattr(
        "mlmm.backends.mlmm_calc.tempfile.TemporaryDirectory",
        workspace_must_not_be_created,
    )
    kwargs = {
        "input_pdb": "does-not-exist.pdb",
        "real_parm7": "does-not-exist.parm7",
        "model_pdb": "does-not-exist-model.pdb",
        keyword: value,
    }
    with pytest.raises(ValueError, match=message):
        MLMMCore(**kwargs)


@pytest.mark.parametrize(
    ("yaml_text", "message"),
    [
        ("calc:\n  hessian_calc_mode: typo\n", "hessian_calc_mode"),
        ("calc:\n  link_atom_method: typo\n", "link_atom_method"),
    ],
)
def test_resolved_yaml_config_uses_the_same_strict_normalizers(
    yaml_text: str, message: str
) -> None:
    calc_cfg = yaml.safe_load(yaml_text)["calc"]
    with pytest.raises(ValueError, match=message):
        normalize_calculator_methods(calc_cfg)


@pytest.mark.parametrize(
    ("key", "value", "message"),
    [
        ("hessian_calc_mode", "typo", "hessian_calc_mode"),
        ("link_atom_method", "typo", "link_atom_method"),
    ],
)
def test_opt_dry_run_rejects_invalid_yaml_method(
    tmp_path: Path, key: str, value: str, message: str
) -> None:
    repo = Path(__file__).resolve().parents[1]
    in_pdb = repo / "examples" / "toy_system" / "p_toy.pdb"
    parm = repo / "examples" / "toy_system" / "p_toy.parm7"
    if not (in_pdb.exists() and parm.exists()):
        pytest.skip("toy_system example inputs not present")
    config = tmp_path / "invalid-method.yaml"
    config.write_text(
        f"calc:\n  {key}: {value}\n",
        encoding="utf-8",
    )

    result = CliRunner().invoke(
        root_cli,
        [
            "opt",
            "-i",
            str(in_pdb),
            "--parm",
            str(parm),
            "-q",
            "0",
            "--detect-layer",
            "--config",
            str(config),
            "--dry-run",
            "--out-dir",
            str(tmp_path / "out"),
        ],
    )

    assert result.exit_code != 0
    assert message in result.output
    assert "[Dry run] --dry-run completed." not in result.output


@pytest.mark.parametrize(
    ("module_name", "command_args"),
    [
        (
            "mlmm.workflows.path_opt",
            lambda smoke, out: [
                "-i",
                str(smoke / "r_complex_layered.pdb"),
                str(smoke / "p_complex_layered.pdb"),
                "--parm",
                str(smoke / "p_complex.parm7"),
                "-q",
                "-1",
                "--dry-run",
                "--out-dir",
                str(out),
            ],
        ),
        (
            "mlmm.workflows.path_search",
            lambda smoke, out: [
                "-i",
                str(smoke / "r_complex_layered.pdb"),
                "-i",
                str(smoke / "p_complex_layered.pdb"),
                "--parm",
                str(smoke / "p_complex.parm7"),
                "-q",
                "-1",
                "--dry-run",
                "--out-dir",
                str(out),
            ],
        ),
    ],
)
@pytest.mark.parametrize(
    ("key", "message"),
    [
        ("hessian_calc_mode", "hessian_calc_mode"),
        ("link_atom_method", "link_atom_method"),
    ],
)
def test_path_dry_run_rejects_invalid_method_from_final_override(
    tmp_path: Path,
    monkeypatch,
    module_name: str,
    command_args,
    key: str,
    message: str,
) -> None:
    import importlib

    repo = Path(__file__).resolve().parents[1]
    smoke = repo / "tests" / "smoke"
    required = [
        smoke / "r_complex_layered.pdb",
        smoke / "p_complex_layered.pdb",
        smoke / "p_complex.parm7",
    ]
    if not all(path.exists() for path in required):
        pytest.skip("smoke inputs not present")

    override = tmp_path / f"{key}-override.yaml"
    override.write_text(f"calc:\n  {key}: typo\n", encoding="utf-8")
    module = importlib.import_module(module_name)
    monkeypatch.setattr(
        module,
        "resolve_yaml_sources",
        lambda config_yaml, override_yaml, args_yaml_legacy: (
            config_yaml,
            override,
            False,
        ),
    )

    result = CliRunner().invoke(
        module.cli,
        command_args(smoke, tmp_path / f"out-{key}"),
    )

    assert result.exit_code != 0
    assert message in result.output
    assert "[Dry run] --dry-run completed." not in result.output
