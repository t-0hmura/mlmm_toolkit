from __future__ import annotations

import torch

from hessian_ff.loaders import load_coords


class _DummyASEAtoms:
    def __init__(self, positions):
        self._positions = positions

    def get_positions(self):
        return self._positions


def test_load_coords_xyz(tmp_path) -> None:
    xyz_path = tmp_path / "mini.xyz"
    xyz_path.write_text(
        "2\n"
        "comment\n"
        "H 0.0 1.0 2.0\n"
        "O 3.0 4.0 5.0\n",
        encoding="utf-8",
    )
    x = load_coords(xyz_path, natom=2, dtype=torch.float64, device="cpu")
    assert tuple(x.shape) == (2, 3)
    assert x.dtype == torch.float64
    assert float(x[1, 2]) == 5.0


def test_load_coords_ase_like() -> None:
    atoms = _DummyASEAtoms([[0.0, 0.1, 0.2], [1.0, 1.1, 1.2]])
    x = load_coords(atoms, natom=2, dtype=torch.float32, device="cpu")
    assert tuple(x.shape) == (2, 3)
    assert x.dtype == torch.float32
    assert abs(float(x[0, 1]) - 0.1) < 1.0e-6


def test_load_coords_tensor() -> None:
    x0 = torch.tensor([[1.0, 2.0, 3.0]], dtype=torch.float64)
    x = load_coords(x0, natom=1, dtype=torch.float32, device="cpu")
    assert tuple(x.shape) == (1, 3)
    assert x.dtype == torch.float32
    assert float(x[0, 0]) == 1.0
