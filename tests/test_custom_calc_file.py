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
