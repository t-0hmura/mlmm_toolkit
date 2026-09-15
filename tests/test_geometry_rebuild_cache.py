"""Actual Geometry cache regressions; run separately against each product."""

from pathlib import Path
import runpy

import numpy as np
import pytest
import pysisyphus

from pysisyphus.Geometry import Geometry
from pysisyphus.calculators.Calculator import Calculator
from pysisyphus.intcoords.exceptions import (
    NeedNewInternalsException,
    RebuiltInternalsException,
)
from pysisyphus.optimizers.guess_hessians import get_guess_hessian


# Reuse the fixture from the product selected by the caller's import path.
_fixture = runpy.run_path(str(
    Path(pysisyphus.__file__).resolve().parents[1]
    / "tests" / "test_freeze_noncart_geometry.py"
))
ATOMS, COORDS, FREEZE = (_fixture[key] for key in ("ATOMS", "COORDS", "FREEZE"))


class CountingQuarticCalculator(Calculator):
    def __init__(self, out_dir):
        super().__init__(out_dir=out_dir, check_mem=False)
        self.calls = []

    @staticmethod
    def exact(coords):
        x = np.asarray(coords, dtype=float).reshape(-1)
        a = 0.05
        return {
            "energy": float(np.sum(0.5 * x**2 + 0.25 * a * x**4)),
            "forces": -(x + a * x**3),
            "hessian": np.diag(1.0 + 3.0 * a * x**2),
        }

    def _evaluate(self, kind, coords):
        self.calls.append((kind, np.asarray(coords, dtype=float).copy()))
        values = self.exact(coords)
        keys = {
            "energy": ("energy",),
            "forces": ("energy", "forces"),
            "hessian": ("energy", "forces", "hessian"),
        }[kind]
        return {key: values[key] for key in keys}

    def get_energy(self, atoms, coords, **kwargs):
        return self._evaluate("energy", coords)

    def get_forces(self, atoms, coords, **kwargs):
        return self._evaluate("forces", coords)

    def get_hessian(self, atoms, coords, **kwargs):
        return self._evaluate("hessian", coords)


def _seed_geometry(coord_type, tmp_path):
    geom = Geometry(ATOMS, COORDS.copy(), coord_type=coord_type, freeze_atoms=FREEZE)
    calc = CountingQuarticCalculator(tmp_path)
    geom.set_calculator(calc)
    _ = geom.energy, geom.cart_forces, geom.cart_hessian
    assert [kind for kind, _ in calc.calls] == ["energy", "forces", "hessian"]
    return geom, calc


def _assert_current_results(geom, calc):
    expected = calc.exact(geom.cart_coords)
    expected["forces"].reshape(-1, 3)[FREEZE] = 0.0
    guess, label = get_guess_hessian(geom, "calc")
    assert label == "calculated exact"
    assert geom.energy == pytest.approx(expected["energy"])
    np.testing.assert_allclose(geom.cart_forces, expected["forces"])
    np.testing.assert_allclose(geom.cart_hessian, expected["hessian"])
    int_gradient = geom.internal.transform_forces(-expected["forces"])
    expected_guess = geom.internal.transform_hessian(expected["hessian"], int_gradient)
    np.testing.assert_allclose(guess, expected_guess, atol=1e-10, rtol=1e-10)
    np.testing.assert_allclose(geom.internal.coords3d, geom.coords3d, atol=1e-12)


@pytest.mark.parametrize("coord_type", ["redund", "dlc", "tric"])
def test_rebuild_refreshes_committed_geometry_results(coord_type, tmp_path, monkeypatch):
    geom, calc = _seed_geometry(coord_type, tmp_path)
    frozen_ref = geom.coords3d[FREEZE].copy()
    committed = geom.coords3d.copy()
    committed[1] += 0.1

    def raise_rebuild(_int_step, update_constraints=False):
        trial = committed.copy()
        trial[FREEZE] += 5.0
        raise NeedNewInternalsException(trial, invalid_inds=())

    monkeypatch.setattr(geom.internal, "transform_int_step", raise_rebuild)
    with pytest.raises(RebuiltInternalsException):
        geom.coords = geom.coords.copy() + 1e-3

    np.testing.assert_allclose(geom.coords3d, committed, atol=1e-12)
    np.testing.assert_allclose(geom.coords3d[FREEZE], frozen_ref, atol=1e-12)
    _assert_current_results(geom, calc)
    assert len(calc.calls) == 4
    assert calc.calls[-1][0] == "hessian"
    np.testing.assert_array_equal(calc.calls[-1][1], geom.cart_coords)


@pytest.mark.parametrize("coord_type", ["redund", "dlc", "tric"])
def test_normal_internal_setter_refreshes_results(coord_type, tmp_path):
    geom, calc = _seed_geometry(coord_type, tmp_path)
    previous = geom.cart_coords.copy()
    frozen_ref = geom.coords3d[FREEZE].copy()
    geom.coords = geom.coords.copy() + 1e-3
    assert not np.array_equal(geom.cart_coords, previous)
    np.testing.assert_allclose(geom.coords3d[FREEZE], frozen_ref, atol=1e-12)
    _assert_current_results(geom, calc)
    assert len(calc.calls) == 4
    assert calc.calls[-1][0] == "hessian"
    np.testing.assert_array_equal(calc.calls[-1][1], geom.cart_coords)


@pytest.mark.parametrize("coord_type", ["redund", "dlc", "tric"])
def test_unchanged_reads_reuse_cached_results(coord_type, tmp_path):
    geom, calc = _seed_geometry(coord_type, tmp_path)
    _assert_current_results(geom, calc)
    _assert_current_results(geom, calc)
    assert len(calc.calls) == 3
