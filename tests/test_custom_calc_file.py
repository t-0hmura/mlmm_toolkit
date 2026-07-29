"""Smoke + unit tests for the ``--calc-file`` custom ML-region backend.

Exercises loading an arbitrary ASE Calculator from a user Python file and using
it as the ML-region engine of the ML/MM ONIOM coupling (the R1 reviewer point:
couple GFN-xTB / DFTB+ / any ASE engine). Uses a dependency-free, element-
agnostic toy harmonic calculator so the test needs no MLIP weights or GPU.
"""

from __future__ import annotations

import textwrap
from pathlib import Path

import numpy as np
import pytest
import torch

from mlmm.backends import apply_calc_file_to_calc_cfg

# Toy ASE calculator: V = 0.5 * sum(pos**2) eV, F = -pos eV/Ang.
TOY_CALC = textwrap.dedent(
    '''
    import numpy as np
    from ase.calculators.calculator import Calculator, all_changes


    class ToyHarmonic(Calculator):
        implemented_properties = ["energy", "forces"]

        def calculate(self, atoms=None, properties=("energy",), system_changes=all_changes):
            super().calculate(atoms, properties, system_changes)
            pos = atoms.get_positions()
            self.results["energy"] = 0.5 * float(np.sum(pos ** 2))
            self.results["forces"] = -pos


    def get_calculator(charge=0, spin=1, device="auto", **kwargs):
        return ToyHarmonic()
    '''
)


def _write(path: Path, text: str) -> Path:
    path.write_text(text, encoding="utf-8")
    return path


def test_load_ase_calculator(tmp_path: Path) -> None:
    from mlmm.backends.custom import load_ase_calculator

    calc_file = _write(tmp_path / "toy.py", TOY_CALC)
    ase_calc = load_ase_calculator(str(calc_file))
    assert hasattr(ase_calc, "get_potential_energy")
    assert hasattr(ase_calc, "get_forces")


def test_calculator_class_export_is_instantiated(tmp_path: Path) -> None:
    from ase.calculators.calculator import Calculator

    from mlmm.backends.custom import load_ase_calculator

    calc_file = _write(
        tmp_path / "class_calc.py",
        "from ase.calculators.emt import EMT\nget_calculator = EMT\n",
    )

    assert isinstance(load_ase_calculator(str(calc_file)), Calculator)


def test_load_ase_calculator_errors(tmp_path: Path) -> None:
    from mlmm.backends.custom import load_ase_calculator

    missing = _write(tmp_path / "no_factory.py", "x = 1\n")
    with pytest.raises(ValueError):
        load_ase_calculator(str(missing))

    not_a_calc = _write(tmp_path / "bad.py", "def get_calculator(**kw):\n    return 42\n")
    with pytest.raises(ValueError):
        load_ase_calculator(str(not_a_calc))


def test_custom_backend_eval(tmp_path: Path) -> None:
    from ase import Atoms

    from mlmm.backends.mlmm_calc import _CustomBackend

    calc_file = _write(tmp_path / "toy.py", TOY_CALC)
    backend = _CustomBackend(
        calc_file=str(calc_file),
        model_charge=0,
        model_mult=1,
        ml_device=torch.device("cpu"),
    )
    coord = np.array([[0.0, 0.0, 0.0], [0.0, 0.757, 0.587], [0.0, -0.757, 0.587]])
    atoms = Atoms(symbols=["O", "H", "H"], positions=coord)
    energy, forces, _ = backend.eval(atoms, need_grad=True)
    assert abs(energy - 0.5 * float(np.sum(coord ** 2))) < 1e-9
    assert np.allclose(forces, -coord)


def test_apply_calc_file_switches_backend() -> None:
    cfg = {"backend": "uma", "uma_model": "uma-s-1p1"}
    apply_calc_file_to_calc_cfg(cfg, "/path/to/toy.py", "get_calculator")
    assert cfg["backend"] == "custom"
    assert cfg["calc_file"] == "/path/to/toy.py"
    assert cfg["calc_factory"] == "get_calculator"
    assert "uma_model" not in cfg  # per-backend model defaults dropped

    # No calc-file -> the --backend selection is untouched.
    cfg2 = {"backend": "uma"}
    apply_calc_file_to_calc_cfg(cfg2, None, None)
    assert cfg2["backend"] == "uma"


def test_custom_calculator_provenance_uses_loader_default_factory() -> None:
    from mlmm.core.utils import calculator_provenance

    provenance = calculator_provenance({
        "backend": "custom",
        "calc_file": "/tmp/toy.py",
    })

    assert provenance["mlip_backend"] == "custom"
    assert provenance["mlip_model"] == "toy.py:get_calculator"
    assert provenance["mlip_precision"] is None


@pytest.mark.parametrize(
    ("backend", "key", "value", "expected"),
    [
        ("uma", "uma_precision", "fp64", "fp64"),
        ("orb", "orb_precision", "float32-high", "fp32"),
        ("orb", "orb_precision", "float32-highest", "fp32"),
        ("mace", "mace_dtype", "float32", "fp32"),
        ("aimnet2", None, "fp32", "fp32"),
    ],
)
def test_calculator_provenance_records_effective_precision(
    backend, key, value, expected
) -> None:
    from mlmm.core.utils import calculator_provenance

    cfg = {"backend": backend}
    if key is not None:
        cfg[key] = value
    assert calculator_provenance(cfg)["mlip_precision"] == expected


def test_sp_yaml_custom_factory_is_not_overwritten_by_cli_default(
    tmp_path: Path,
) -> None:
    from click.testing import CliRunner

    from mlmm.cli import cli as root_cli

    calc_file = _write(
        tmp_path / "toy_build.py", TOY_CALC.replace("def get_calculator", "def build")
    )
    input_pdb = _write(
        tmp_path / "carbon.pdb",
        "HETATM    1  C1  MOL A   1       0.000   0.000   0.000  1.00  0.00           C\nEND\n",
    )
    parm = _write(tmp_path / "empty.parm7", "placeholder\n")
    config = _write(
        tmp_path / "config.yaml",
        f"calc:\n  calc_file: {calc_file}\n  calc_factory: build\n",
    )

    result = CliRunner().invoke(
        root_cli,
        [
            "sp",
            "-i",
            str(input_pdb),
            "--parm",
            str(parm),
            "-q",
            "0",
            "-m",
            "1",
            "--config",
            str(config),
            "--show-config",
        ],
        catch_exceptions=False,
    )
    assert result.exit_code == 0, result.output
    assert "backend: custom" in result.output
    assert "calc_factory: build" in result.output
