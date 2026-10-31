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


@pytest.mark.parametrize("backend", ["hessian_ff", "openmm"])
@pytest.mark.parametrize("fixed", [[], [4, 1], [0, 1, 2, 3, 4, 5]])
def test_mm_finite_difference_partial_equals_restricted_full(backend, fixed) -> None:
    import numpy as np
    from ase import Atoms
    from ase.constraints import FixAtoms

    from mlmm.backends.mlmm_calc import OpenMMCalculator, hessianffCalculator

    n_atoms = 6
    vsites = [2] if backend == "openmm" else []
    base = hessianffCalculator if backend == "hessian_ff" else OpenMMCalculator

    class _Quadratic(base):
        """No topology: forces of E = x.K.x/2 with a fixed dense symmetric K."""

        def __init__(self) -> None:
            a = np.random.default_rng(5).standard_normal((3 * n_atoms, 3 * n_atoms))
            self.k = a + a.T
            self._vsite_idx = list(vsites)

        def get_forces(self, atoms=None):
            return -(self.k @ atoms.get_positions().reshape(-1)).reshape(-1, 3)

    calculator = _Quadratic()
    positions = np.random.default_rng(1).standard_normal((n_atoms, 3))
    atoms = Atoms("H6", positions=positions)
    if fixed:
        atoms.set_constraint(FixAtoms(indices=fixed))

    full, no_atoms = calculator.finite_difference_hessian(
        atoms, delta=1.0e-3, return_partial_hessian=False
    )
    partial, active = calculator.finite_difference_hessian(
        atoms, delta=1.0e-3, return_partial_hessian=True
    )

    expected = [i for i in range(n_atoms) if i not in set(fixed) | set(vsites)]
    assert no_atoms is None
    assert full.shape == (3 * n_atoms, 3 * n_atoms)
    np.testing.assert_array_equal(active, expected)
    assert active.dtype == np.int64
    idx3 = (3 * np.asarray(expected, dtype=np.int64)[:, None] + np.arange(3)).ravel()
    assert partial.shape == (idx3.size, idx3.size)
    assert partial.tobytes() == full[np.ix_(idx3, idx3)].tobytes()


@pytest.mark.parametrize("active", [[0, 2, 5], [5, 1, 3], [], [0, 1, 2, 3, 4, 5]])
@pytest.mark.parametrize("dtype", ["float64", "float32"])
def test_expand_partial_hessian_matches_block_loop(active, dtype) -> None:
    import numpy as np

    from mlmm.backends.mlmm_calc import _expand_partial_hessian

    n_atoms = 6
    n_dof = 3 * len(active)
    h_sub = np.random.default_rng(11).standard_normal((n_dof, n_dof))
    expected = np.zeros((3 * n_atoms, 3 * n_atoms), dtype=dtype)
    for i_local, i_atom in enumerate(active):
        for j_local, j_atom in enumerate(active):
            expected[3 * i_atom:3 * i_atom + 3, 3 * j_atom:3 * j_atom + 3] = h_sub[
                3 * i_local:3 * i_local + 3, 3 * j_local:3 * j_local + 3
            ]

    actual = _expand_partial_hessian(
        h_sub, np.asarray(active, dtype=int), n_atoms, dtype=np.dtype(dtype)
    )

    assert actual.dtype == expected.dtype
    assert actual.shape == expected.shape
    assert actual.tobytes() == expected.tobytes()


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
def test_path_dry_run_rejects_invalid_method_from_config(
    tmp_path: Path,
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

    config = tmp_path / f"{key}-config.yaml"
    config.write_text(f"calc:\n  {key}: typo\n", encoding="utf-8")
    module = importlib.import_module(module_name)

    result = CliRunner().invoke(
        module.cli,
        [*command_args(smoke, tmp_path / f"out-{key}"), "--config", str(config)],
    )

    assert result.exit_code != 0
    assert message in result.output
    assert "[Dry run] --dry-run completed." not in result.output
