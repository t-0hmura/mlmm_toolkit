"""DMF core ownership, solver outcome, and frozen-anchor regressions."""

from __future__ import annotations

from types import ModuleType, SimpleNamespace
from typing import Any
import sys

import numpy as np
import pytest

from mlmm.workflows import path_opt, path_search


class _Geom:
    def __init__(self, positions, *, freeze_atoms=()):
        self.positions = np.asarray(positions, dtype=float)
        self.freeze_atoms = np.asarray(freeze_atoms, dtype=int)
        self.coord_type = "cart"

    def as_xyz(self) -> str:
        lines = [str(len(self.positions)), "test geometry"]
        lines.extend(
            f"H {x:.12f} {y:.12f} {z:.12f}" for x, y, z in self.positions
        )
        return "\n".join(lines) + "\n"


class _SharedCalculator:
    def __init__(self):
        self.core = SimpleNamespace(model_charge=-1, model_mult=2)
        self.energy_calls = 0

    def get_energy(self, elements, coords_bohr):
        self.energy_calls += 1
        coords = np.asarray(coords_bohr, dtype=float)
        return {"energy": float(np.sum(coords * coords))}


def _install_fake_dmf(monkeypatch, *, status: int, reason: bytes = b"status"):
    module = ModuleType("dmf")

    class _Interpolated:
        def __init__(self, images):
            self.images = [image.copy() for image in images]
            self.coefs = np.ones((1,), dtype=float)

    class _DirectMaxFlux:
        def __init__(self, images, **kwargs):
            first = images[0].copy()
            last = images[-1].copy()
            middle = first.copy()
            middle.set_positions(
                0.5 * (first.get_positions() + last.get_positions())
            )
            self.images = [first, middle, last]
            self.options = []

        def add_ipopt_options(self, options):
            self.options.append(dict(options))

        def solve(self, tol):
            return np.zeros((1,), dtype=float), {
                "status": status,
                "status_msg": reason,
            }

    def _interpolate(images, **kwargs):
        return _Interpolated(images)

    module.DirectMaxFlux = _DirectMaxFlux
    module.interpolate_fbenm = _interpolate
    monkeypatch.setitem(sys.modules, "dmf", module)


def test_solver_status_zero_is_the_only_converged_outcome() -> None:
    assert path_opt._dmf_solver_outcome(
        (object(), {"status": 0, "status_msg": b"Solve succeeded"})
    ) == (True, 0, "Solve succeeded")
    assert path_opt._dmf_solver_outcome(
        (object(), {"status": -1, "status_msg": "Maximum iterations exceeded"})
    ) == (False, -1, "Maximum iterations exceeded")
    assert path_opt._dmf_solver_outcome(object()) == (
        False,
        None,
        "IPOPT status was not reported.",
    )


def test_dmf_interpolation_cache_is_released_without_emptying_mid_phase() -> None:
    calls = []

    class Stage:
        def release_device_cache(self, *, empty_cache):
            calls.append(empty_cache)

    path_opt._release_dmf_interpolation_cache(Stage())
    path_opt._release_dmf_interpolation_cache(object())  # pydmf 1.2 fallback

    assert calls == [False]


def test_torch_dmf_runtime_options_disable_unused_history_and_preserve_precision(
    monkeypatch,
) -> None:
    monkeypatch.setattr(path_opt.torch.cuda, "is_available", lambda: True)
    assert path_opt._torch_dmf_runtime_kwargs(
        "cpu", {"keep_history": True}, {}, {}
    ) == {}
    assert path_opt._torch_dmf_runtime_kwargs(
        "gpu",
        {"device": "cpu", "dtype": "float64"},
        {"device": "cuda", "dtype": "float32"},
        {"dtype": "float64"},
    ) == {"keep_history": False, "device": "cuda", "dtype": "float64"}
    assert path_opt._torch_dmf_runtime_kwargs(
        "gpu", {"keep_history": True}, {}, {}
    ) == {"keep_history": True, "device": "cuda"}
    assert path_opt._torch_dmf_runtime_kwargs(
        "gpu", {"keep_history": True}, {}, {}, supports_keep_history=False
    ) == {"device": "cuda"}


def test_torch_dmf_gpu_requires_a_visible_or_explicit_device(monkeypatch) -> None:
    monkeypatch.setattr(path_opt.torch.cuda, "is_available", lambda: False)

    with pytest.raises(RuntimeError, match="requires a visible CUDA device"):
        path_opt._torch_dmf_runtime_kwargs("gpu", {}, {}, {})

    assert path_opt._torch_dmf_runtime_kwargs(
        "gpu", {"device": "cpu"}, {}, {}
    ) == {"keep_history": False, "device": "cpu"}


def test_nonconverged_result_payload_does_not_become_completed() -> None:
    result = path_opt.DMFMepResult(
        images=(object(), object(), object()),
        energies=(-10.0, -9.5, -10.2),
        hei_idx=1,
        converged=False,
        ipopt_status=-1,
        reason="Maximum iterations exceeded",
    )
    payload = path_opt._build_dmf_result_data(
        result,
        {"backend": "uma", "model_charge": 0, "model_mult": 1},
    )
    assert payload["status"] == "not_converged"
    assert payload["converged"] is False
    assert payload["ipopt_status"] == -1
    assert payload["reason"] == "Maximum iterations exceeded"
    assert payload["image_energies_hartree"] == [-10.0, -9.5, -10.2]
    assert payload["hei_index"] == 1


def test_path_opt_dmf_reuses_one_shared_core_for_solve_and_final_energy(
    tmp_path, monkeypatch
) -> None:
    _install_fake_dmf(
        monkeypatch,
        status=-1,
        reason=b"Maximum Number of Iterations Exceeded",
    )
    seen_cores = []
    original_ase_wrapper = path_opt.MLMMASECalculator

    def _capture_wrapper(*, core):
        seen_cores.append(core)
        return original_ase_wrapper(core=core)

    monkeypatch.setattr(path_opt, "MLMMASECalculator", _capture_wrapper)
    monkeypatch.setattr(
        path_opt,
        "mlmm",
        lambda **kwargs: pytest.fail("DMF constructed a second heavy mlmm core"),
    )

    shared = _SharedCalculator()
    result = path_opt._run_dmf_mep(
        [
            _Geom([[0.0, 0.0, 0.0], [0.0, 0.0, 1.0]], freeze_atoms=[0]),
            _Geom([[1.0, 0.0, 0.0], [0.0, 0.0, 1.5]], freeze_atoms=[0]),
        ],
        shared,
        tmp_path,
        [tmp_path / "input.xyz"],
        3,
        [0],
        dmf_cfg={"backend": "cpu", "ipopt_options": {}},
    )

    assert seen_cores == [shared.core]
    assert shared.energy_calls == len(result.images) == 3
    assert result.converged is False
    assert result.ipopt_status == -1
    assert result.reason == "Maximum Number of Iterations Exceeded"
    assert all(image.calc is None for image in result.images)
    assert (tmp_path / "final_geometries_trj.xyz").exists()
    assert (tmp_path / "hei.xyz").exists()


def test_shared_frozen_reference_is_copied_from_first_image_and_read_only() -> None:
    from ase import Atoms

    first = Atoms("HH", positions=[[0.0, 0.0, 0.0], [0.0, 0.0, 1.0]])
    second = Atoms("HH", positions=[[4.0, 0.0, 0.0], [0.0, 0.0, 2.0]])
    reference = path_opt._shared_frozen_reference([first, second], [0])

    assert reference is not None
    np.testing.assert_allclose(reference, [[0.0, 0.0, 0.0]])
    first.positions[0, 0] = 9.0
    np.testing.assert_allclose(reference, [[0.0, 0.0, 0.0]])
    assert reference.flags.writeable is False
    with pytest.raises(ValueError, match="outside the path image bounds"):
        path_opt._shared_frozen_reference([first, second], [2])


def test_recursive_dmf_reuses_first_image_anchor_for_every_restraint(
    tmp_path, monkeypatch
) -> None:
    _install_fake_dmf(monkeypatch, status=0, reason=b"Solve succeeded")
    captured = []
    from ase.calculators.calculator import Calculator
    from mlmm.workflows import restraints

    class _CapturingHarmonic(Calculator):
        implemented_properties = ["energy", "forces"]

        def __init__(self, *, indices, ref_positions, k_fix):
            super().__init__()
            captured.append(np.asarray(ref_positions, dtype=float).copy())

    monkeypatch.setattr(restraints, "HarmonicFixAtoms", _CapturingHarmonic)
    monkeypatch.setattr(path_search, "run_trj2fig", lambda *args, **kwargs: None)

    shared = _SharedCalculator()
    result = path_search._run_dmf_between(
        _Geom([[0.0, 0.0, 0.0], [0.0, 0.0, 1.0]], freeze_atoms=[0]),
        _Geom([[3.0, 0.0, 0.0], [0.0, 0.0, 2.0]], freeze_atoms=[0]),
        shared,
        {"model_charge": -1, "model_mult": 2},
        tmp_path,
        "seg",
        None,
        3,
        {"backend": "cpu", "ipopt_options": {}},
    )

    assert len(result.images) == len(captured) == 3
    for reference in captured:
        np.testing.assert_allclose(reference, [[0.0, 0.0, 0.0]])
    assert not np.allclose(captured[-1], [[3.0, 0.0, 0.0]])
