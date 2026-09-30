"""A restart array without its backend spec is a corrupt checkpoint, not NumPy/CPU."""

from __future__ import annotations

import numpy as np
import pytest
import yaml

from pysisyphus.Geometry import Geometry
from pysisyphus.calculators.Calculator import Calculator
from pysisyphus.optimizers import checkpoint
from pysisyphus.optimizers.RFOptimizer import RFOptimizer


class _QuadraticCalculator(Calculator):
    def __init__(self, out_dir):
        super().__init__(out_dir=out_dir, check_mem=False)

    @staticmethod
    def _results(coords):
        coords = np.asarray(coords, dtype=float)
        return {"energy": float(coords @ coords), "forces": -2.0 * coords}

    def get_energy(self, atoms, coords, **prepare_kwargs):
        return self._results(coords)

    def get_forces(self, atoms, coords, **prepare_kwargs):
        return self._results(coords)

    def get_hessian(self, atoms, coords, **prepare_kwargs):
        results = self._results(coords)
        results["hessian"] = 2.0 * np.eye(len(coords))
        return results


def _rfo(out_dir, start, **kwargs):
    out_dir.mkdir(parents=True, exist_ok=True)
    geom = Geometry(["H"], np.array(start, dtype=float), coord_type="cart")
    geom.set_calculator(_QuadraticCalculator(out_dir))
    return RFOptimizer(
        geom, hessian_init="unit", trust_radius=0.1, trust_min=1e-4,
        trust_max=0.1, line_search=False, gdiis=False, out_dir=out_dir,
        dump=False, thresh="gau_loose", **kwargs,
    )


def _restart_info(tmp_path):
    opt = _rfo(tmp_path / "a", [0.5, 0.0, 0.0], max_cycles=2)
    opt.run()
    info = opt._get_opt_restart_info()
    info["_sy_buffer_S"] = [[1.0, 0.0, 0.0]]
    info["_sy_buffer_Y"] = [[2.0, 0.0, 0.0]]
    info["_sy_buffer_S_backend"] = [opt._restart_backend(np.ones(3))]
    info["_sy_buffer_Y_backend"] = [opt._restart_backend(np.ones(3))]
    info["_prev_eigvec_min"] = [0.0, 1.0, 0.0]
    info["_prev_eigvec_min_backend"] = opt._restart_backend(np.ones(3))
    return info


def test_complete_backend_specs_restore(tmp_path) -> None:
    info = _restart_info(tmp_path)
    opt = _rfo(tmp_path / "b", [0.5, 0.0, 0.0], max_cycles=2)

    opt._set_opt_restart_info(info)

    assert isinstance(opt.H, np.ndarray)
    np.testing.assert_allclose(opt._sy_buffer_S[0], [1.0, 0.0, 0.0])
    np.testing.assert_allclose(opt._prev_eigvec_min, [0.0, 1.0, 0.0])


@pytest.mark.parametrize(
    "corrupt, label",
    [
        (lambda info: info.__setitem__("H_backend", None), "H_backend"),
        (lambda info: info["_sy_buffer_S_backend"].__setitem__(0, None),
         r"_sy_buffer_S_backend\[0\]"),
        (lambda info: info["_sy_buffer_Y_backend"].__setitem__(0, "numpy"),
         r"_sy_buffer_Y_backend\[0\]"),
        (lambda info: info.__setitem__("_prev_eigvec_min_backend", None),
         "_prev_eigvec_min_backend"),
    ],
    ids=["H", "sy_buffer_S", "sy_buffer_Y", "prev_eigvec_min"],
)
def test_missing_backend_spec_is_rejected(tmp_path, corrupt, label) -> None:
    info = _restart_info(tmp_path)
    corrupt(info)
    opt = _rfo(tmp_path / "b", [0.5, 0.0, 0.0], max_cycles=2)

    with pytest.raises(ValueError, match=rf"Corrupt checkpoint: {label} must be a mapping"):
        opt._set_opt_restart_info(info)


def test_checkpoint_file_without_hessian_backend_is_rejected(tmp_path) -> None:
    opt_a = _rfo(tmp_path / "a", [0.5, 0.0, 0.0], max_cycles=2)
    opt_a.run()
    ck = tmp_path / "restart.yaml"
    checkpoint.save_checkpoint(opt_a, ck)
    payload = yaml.safe_load(ck.read_text())
    payload["restart_info"]["H_backend"] = None
    ck.write_text(yaml.safe_dump(payload))

    opt_b = _rfo(tmp_path / "b", [9.9, 9.9, 9.9], max_cycles=200)
    with pytest.raises(
        checkpoint.CheckpointValidationError, match="H_backend must be a mapping"
    ):
        checkpoint.load_and_apply(opt_b, ck)

    # Rejected before set_restart_info could overwrite the base history.
    assert opt_b.coords == [] and opt_b.energies == []


@pytest.mark.parametrize(
    "corrupt, message",
    [
        (lambda info: None, None),
        (lambda info: info.__setitem__("H_backend", None), "H_backend must be a mapping"),
        (lambda info: info["_sy_buffer_S_backend"].__setitem__(0, None),
         r"_sy_buffer_S_backend\[0\] must be a mapping"),
        (lambda info: info.__setitem__("_sy_buffer_Y_backend", []),
         "_sy_buffer_Y_backend must hold one spec per _sy_buffer_Y entry"),
        (lambda info: info.__setitem__("_prev_eigvec_min_backend", None),
         "_prev_eigvec_min_backend must be a mapping"),
    ],
    ids=["complete", "H", "sy_buffer_S", "sy_buffer_Y_length", "prev_eigvec_min"],
)
def test_validate_payload_checks_backend_specs(tmp_path, corrupt, message) -> None:
    opt_a = _rfo(tmp_path / "a", [0.5, 0.0, 0.0], max_cycles=2)
    opt_a.run()
    ck = tmp_path / "restart.yaml"
    checkpoint.save_checkpoint(opt_a, ck)
    payload = checkpoint.load_payload(ck)
    info = payload["restart_info"]
    info["_sy_buffer_S"] = [[1.0, 0.0, 0.0]]
    info["_sy_buffer_Y"] = [[2.0, 0.0, 0.0]]
    info["_sy_buffer_S_backend"] = [opt_a._restart_backend(np.ones(3))]
    info["_sy_buffer_Y_backend"] = [opt_a._restart_backend(np.ones(3))]
    info["_prev_eigvec_min"] = [0.0, 1.0, 0.0]
    info["_prev_eigvec_min_backend"] = opt_a._restart_backend(np.ones(3))
    corrupt(info)
    opt_b = _rfo(tmp_path / "b", [9.9, 9.9, 9.9], max_cycles=200)

    if message is None:
        checkpoint.validate_payload(payload, opt_b)
    else:
        with pytest.raises(checkpoint.CheckpointValidationError, match=message):
            checkpoint.validate_payload(payload, opt_b)
