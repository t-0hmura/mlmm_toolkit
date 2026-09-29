"""``--opt-mode grad|hess`` on scan, scan2d and scan3d.

``hess`` relaxes every scan step with a standard RFO whose initial Hessian is
the exact Hessian of the PES attached to the geometry; ``grad`` keeps LBFGS.
"""

from __future__ import annotations

import json
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest
from click.testing import CliRunner
from pysisyphus.constants import ANG2BOHR

SCAN_COMMANDS = ("scan", "scan2d", "scan3d")


def _has_option_header(output: str, option_prefix: str) -> bool:
    for line in output.splitlines():
        stripped = line.lstrip()
        if len(line) - len(stripped) > 2 or not stripped.startswith(option_prefix):
            continue
        tail = stripped[len(option_prefix):]
        if (not tail) or tail[0].isspace() or tail[0] in {",", "/"}:
            return True
    return False


def _command(name: str):
    from mlmm.workflows import scan, scan2d, scan3d

    return {"scan": scan, "scan2d": scan2d, "scan3d": scan3d}[name].cli


class _QuadraticCalculator:
    """E = 0.5 * k * |x - x0|^2 (Hartree, Bohr) with an exact Hessian.

    Without ``x0`` energy and forces are zero everywhere; only the Hessian is k * I.
    """

    freeze_atoms: list = []
    analytical_2d = False

    def __init__(self, x0=None, k: float = 0.5):
        self.x0 = None if x0 is None else np.asarray(x0, dtype=float).reshape(-1)
        self.k = float(k)
        self.hessian_calls = 0

    def _delta(self, positions):
        x = np.asarray(positions, dtype=float).reshape(-1)
        return x - (x if self.x0 is None else self.x0)

    def get_energy(self, atoms, positions, **kwargs):
        d = self._delta(positions)
        return {"energy": 0.5 * self.k * float(d @ d)}

    def get_forces(self, atoms, positions, **kwargs):
        d = self._delta(positions)
        return {"energy": 0.5 * self.k * float(d @ d), "forces": -self.k * d}

    def get_hessian(self, atoms, positions, **kwargs):
        self.hessian_calls += 1
        d = self._delta(positions)
        return {
            "energy": 0.5 * self.k * float(d @ d),
            "forces": -self.k * d,
            "hessian": self.k * np.eye(d.size),
        }


def _three_atom_geometry(calc):
    from pysisyphus.Geometry import Geometry

    coords = np.array([[0.0, 0.0, 0.0], [1.5, 0.0, 0.0], [0.0, 2.0, 0.0]]) * ANG2BOHR
    geom = Geometry(["C", "C", "C"], coords.reshape(-1), coord_type="cart")
    geom.set_calculator(calc)
    return geom


def _recording_rfo(monkeypatch):
    """Replace the RFO class used by the scan helper; record each construction.

    ``run`` places each restrained second atom at its target along x.
    """
    from mlmm.workflows import scan_common

    calls = []

    def factory(geometry, **kwargs):
        seeded = geometry._hessian
        if hasattr(seeded, "detach"):
            # On a CUDA node the seeded ML/MM Hessian is a device tensor.
            seeded = seeded.detach().cpu().numpy()
        calls.append({
            "prefix": kwargs.get("prefix"),
            "kwargs": dict(kwargs),
            "seeded": None if seeded is None else np.asarray(seeded).copy(),
        })

        def run():
            coords = np.array(geometry.coords3d, dtype=float)
            for first, second, target in getattr(geometry.calculator, "_restraints", []):
                coords[second] = coords[first] + np.array([float(target) * ANG2BOHR, 0.0, 0.0])
            geometry.coords3d = coords

        return SimpleNamespace(run=run, is_converged=True, is_stalled=False, stop_reason="")

    monkeypatch.setattr(scan_common, "RFOptimizer", factory)
    return calls


@pytest.mark.parametrize("command", SCAN_COMMANDS)
def test_opt_mode_is_advanced_only_and_defaults_to_grad(command) -> None:
    import click

    from mlmm.cli import cli as root_cli

    runner = CliRunner()
    primary = runner.invoke(root_cli, [command, "--help"])
    assert primary.exit_code == 0, primary.output
    assert not _has_option_header(primary.output, "--opt-mode"), primary.output
    advanced = runner.invoke(root_cli, [command, "--help-advanced"])
    assert advanced.exit_code == 0, advanced.output
    assert _has_option_header(advanced.output, "--opt-mode"), advanced.output

    param, = [p for p in _command(command).params if p.name == "opt_mode"]
    assert param.help == "Relaxation mode: grad (=LBFGS) or hess (=RFO)."
    assert param.default == "grad"
    assert isinstance(param.type, click.Choice)
    assert list(param.type.choices) == ["grad", "hess"]
    assert param.type.case_sensitive is False


def test_normalize_scan_opt_mode_maps_aliases() -> None:
    import click

    from mlmm.workflows.scan_common import normalize_scan_opt_mode

    for raw in ("grad", "GRAD", "lbfgs"):
        assert normalize_scan_opt_mode(raw) == "lbfgs"
    for raw in ("hess", "Hess", "rfo"):
        assert normalize_scan_opt_mode(raw) == "rfo"
    with pytest.raises(click.BadParameter):
        normalize_scan_opt_mode("dimer")


@pytest.mark.parametrize(
    "yaml_cfg",
    [
        {"rfo": {"trust_max": 0.07, "hessian_update": "bfgs"}},
        {"opt": {"rfo": {"trust_max": 0.07, "hessian_update": "bfgs"}}},
    ],
    ids=["top-level", "nested-under-opt"],
)
def test_rfo_yaml_section_reaches_scan_rfo(yaml_cfg) -> None:
    from mlmm.workflows import scan2d
    from mlmm.workflows.scan_common import resolve_scan_optimizer_configs

    common = dict(
        opt_defaults=scan2d.OPT_BASE_KW,
        lbfgs_defaults=scan2d.LBFGS_KW,
        thresh="baker",
        relax_max_cycles=10000,
        is_param_explicit=lambda name: False,
    )
    opt_cfg, sopt_cfg = resolve_scan_optimizer_configs(
        yaml_cfg, rfo_defaults=scan2d.RFO_KW, kind="rfo", **common,
    )
    assert sopt_cfg["trust_max"] == 0.07
    assert sopt_cfg["hessian_update"] == "bfgs"
    assert sopt_cfg["hessian_init"] == "calc"
    assert sopt_cfg["thresh"] == opt_cfg["thresh"] == "baker"
    assert "rfo" not in opt_cfg and "lbfgs" not in opt_cfg

    # The default optimizer is still LBFGS, and an rfo block never leaks into it.
    opt_cfg, sopt_cfg = resolve_scan_optimizer_configs(yaml_cfg, **common)
    assert "keep_last" in sopt_cfg and "trust_max" not in sopt_cfg
    assert "rfo" not in opt_cfg


def test_build_scan_rfo_kwargs_caps_trust_radii_by_scan_step(tmp_path) -> None:
    from mlmm.workflows import scan
    from mlmm.workflows.scan_common import build_scan_rfo_kwargs

    rfo_cfg = dict(scan.RFO_KW, trust_radius=0.30, trust_max=0.50)
    opt_cfg = {"thresh": "baker", "max_cycles": 7}
    args = build_scan_rfo_kwargs(
        rfo_cfg, opt_cfg, max_step_bohr=0.2, out_dir=tmp_path, prefix="scan_s0001",
    )
    assert args["trust_radius"] == args["trust_max"] == 0.2
    assert args["thresh"] == "baker" and args["max_cycles"] == 7
    assert args["out_dir"] == str(tmp_path) and args["prefix"] == "scan_s0001"
    assert "flatten_enabled" not in args

    args = build_scan_rfo_kwargs(
        dict(scan.RFO_KW), {}, max_step_bohr=1.0, out_dir=tmp_path, prefix="p",
    )
    assert args["trust_radius"] == scan.RFO_KW["trust_radius"]
    assert args["trust_max"] == scan.RFO_KW["trust_max"]


def test_make_scan_rfo_seeds_exact_hessian_of_attached_pes(tmp_path, monkeypatch) -> None:
    from mlmm.workflows import scan
    from mlmm.workflows.scan_common import make_scan_rfo

    calls = _recording_rfo(monkeypatch)
    calc = _QuadraticCalculator(k=0.25)
    geom = _three_atom_geometry(calc)
    make_scan_rfo(
        geom, dict(scan.RFO_KW), {}, max_step_bohr=0.2,
        out_dir=tmp_path, prefix="p", calc_cfg={"ml_device": "cpu"},
    )
    assert calc.hessian_calls == 1
    np.testing.assert_allclose(calls[-1]["seeded"], 0.25 * np.eye(9))

    # A non-exact initial Hessian is left to the optimizer.
    geom = _three_atom_geometry(calc)
    make_scan_rfo(
        geom, dict(scan.RFO_KW, hessian_init="unit"), {}, max_step_bohr=0.2,
        out_dir=tmp_path, prefix="p", calc_cfg={"ml_device": "cpu"},
    )
    assert calc.hessian_calls == 1
    assert calls[-1]["seeded"] is None


def test_make_scan_rfo_runs_real_rfo_on_seeded_hessian(tmp_path) -> None:
    from pysisyphus.optimizers.RFOptimizer import RFOptimizer

    from mlmm.workflows import scan2d
    from mlmm.workflows.scan_common import make_scan_rfo, resolve_scan_optimizer_configs

    opt_cfg, rfo_cfg = resolve_scan_optimizer_configs(
        {},
        opt_defaults=scan2d.OPT_BASE_KW,
        lbfgs_defaults=scan2d.LBFGS_KW,
        rfo_defaults=scan2d.RFO_KW,
        kind="rfo",
        thresh="baker",
        relax_max_cycles=50,
        is_param_explicit=lambda name: name == "relax_max_cycles",
    )
    probe = _three_atom_geometry(_QuadraticCalculator())
    x0 = probe.cart_coords.copy()
    x0[3] += 0.05  # minimum 0.05 Bohr away along x of atom 2
    calc = _QuadraticCalculator(x0=x0, k=0.5)
    geom = _three_atom_geometry(calc)
    optimizer = make_scan_rfo(
        geom, rfo_cfg, opt_cfg, max_step_bohr=0.2,
        out_dir=tmp_path, prefix="rfo", calc_cfg={"ml_device": "cpu"},
    )
    assert isinstance(optimizer, RFOptimizer)
    optimizer.run()
    assert optimizer.is_converged is True
    np.testing.assert_allclose(geom.cart_coords, x0, atol=1e-4)
    # The seeded Hessian is used as is, not recalculated at the first cycle.
    assert calc.hessian_calls == 1


def _write_carbons(path: Path, coords) -> None:
    path.write_text("".join(
        f"HETATM{i:5d}  C{i:<2d} MOL A   1    "
        f"{x:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00           C  \n"
        for i, (x, y, z) in enumerate(coords, 1)
    ) + "END\n")


def _run_scan(tmp_path, monkeypatch, *extra):
    from mlmm.cli import cli as root_cli
    from mlmm.workflows import scan as workflow

    lbfgs_calls = []

    def lbfgs(geometry, *args, **kwargs):
        lbfgs_calls.append(kwargs.get("prefix"))
        return SimpleNamespace(
            run=lambda: None, is_converged=True, is_stalled=False, stop_reason="",
        )

    rfo_calls = _recording_rfo(monkeypatch)
    calc = _QuadraticCalculator()
    monkeypatch.setattr(workflow, "mlmm", lambda **kwargs: calc)
    monkeypatch.setattr(workflow, "LBFGS", lbfgs)
    monkeypatch.setattr(workflow, "_has_bond_change", lambda *args, **kwargs: (False, ""))
    source = tmp_path / "three_carbons.pdb"
    _write_carbons(source, [(0.0, 0.0, 0.0), (1.5, 0.0, 0.0), (0.0, 2.0, 0.0)])
    output = tmp_path / "scan"
    result = CliRunner().invoke(root_cli, [
        "scan", "-i", str(source), "-q", "0", "-m", "1",
        "--parm", str(source), "--model-indices", "1-3", "--no-detect-layer",
        "--preopt", "--endopt", "--no-convert-files", "--out-json", "--one-based",
        "--scan-lists", "[(1,2,1.6)]", "--max-step-size", "0.1",
        "-o", str(output), *extra,
    ])
    assert result.exit_code == 0, result.output + repr(result.exception)
    payload = json.loads((output / "result.json").read_text())
    return result, payload, lbfgs_calls, rfo_calls, calc


def test_scan_default_relaxes_with_lbfgs(tmp_path, monkeypatch) -> None:
    result, payload, lbfgs_calls, rfo_calls, calc = _run_scan(tmp_path, monkeypatch)
    assert rfo_calls == []
    assert lbfgs_calls[0] == "preopt" and lbfgs_calls[-1] == "endopt"
    assert "relaxation (lbfgs)" in result.output
    assert payload["scan_opt_mode"] == "grad"
    assert payload["scan_optimizer"] == "lbfgs"
    assert calc.hessian_calls == 0


@pytest.mark.parametrize(
    "yaml_text",
    [
        "rfo:\n  trust_max: 0.05\n  hessian_update: bfgs\n",
        "opt:\n  rfo:\n    trust_max: 0.05\n    hessian_update: bfgs\n",
    ],
    ids=["top-level", "nested-under-opt"],
)
def test_scan_hess_relaxes_every_stage_with_seeded_rfo(tmp_path, monkeypatch, yaml_text) -> None:
    config = tmp_path / "scan.yaml"
    config.write_text(yaml_text)
    result, payload, lbfgs_calls, rfo_calls, calc = _run_scan(
        tmp_path, monkeypatch, "--opt-mode", "HESS", "--config", str(config),
    )
    assert lbfgs_calls == []
    prefixes = [call["prefix"] for call in rfo_calls]
    assert prefixes[0] == "preopt" and prefixes[-1] == "endopt"
    assert len(prefixes) >= 3
    assert all(p.startswith("scan_s") for p in prefixes[1:-1])
    for call in rfo_calls:
        assert call["seeded"] is not None and call["seeded"].shape == (9, 9)
        assert call["kwargs"]["trust_max"] == 0.05
        assert call["kwargs"]["hessian_update"] == "bfgs"
        assert call["kwargs"]["trust_radius"] == pytest.approx(0.10)
        assert "rfo" not in call["kwargs"]
    assert calc.hessian_calls == len(rfo_calls)
    assert "relaxation (rfo)" in result.output
    assert payload["scan_opt_mode"] == "hess"
    assert payload["scan_optimizer"] == "rfo"


@pytest.mark.parametrize(
    ("command", "scan_lists", "n_points"),
    [
        ("scan2d", "[(1,2,1.500,1.508),(2,3,1.500,1.504)]", 6),
        ("scan3d", "[(1,2,1.500,1.508),(2,3,1.500,1.504),(3,4,1.500,1.504)]", 12),
    ],
)
def test_grid_scan_hess_relaxes_every_point_with_seeded_rfo(
    tmp_path, monkeypatch, command, scan_lists, n_points,
) -> None:
    import re

    from mlmm.cli import cli as root_cli
    from mlmm.workflows import scan2d, scan3d

    module = scan2d if command == "scan2d" else scan3d
    rfo_calls = _recording_rfo(monkeypatch)
    lbfgs_calls = []
    monkeypatch.setattr(module, "_make_lbfgs", lambda *a, **k: lbfgs_calls.append(k))
    calc = _QuadraticCalculator()
    monkeypatch.setattr(module, "mlmm", lambda **kwargs: calc)
    monkeypatch.setattr(scan2d, "write_plotly_image", lambda *args, **kwargs: None)
    source = tmp_path / "four_carbons.pdb"
    _write_carbons(source, [(1.5 * i, 0.0, 0.0) for i in range(4)])
    config = tmp_path / "scan.yaml"
    config.write_text("rfo:\n  hessian_update: bfgs\n")
    out_dir = tmp_path / "out"
    result = CliRunner().invoke(root_cli, [
        command, "-i", str(source), "--parm", str(source),
        "--model-indices", "1-4", "--no-detect-layer", "-q", "0", "-m", "1",
        "--scan-lists", scan_lists, "--max-step-size", "0.004", "--preopt",
        "--opt-mode", "hess", "--config", str(config), "--out-dir", str(out_dir),
    ])
    assert result.exit_code == 0, result.output + repr(result.exception)
    assert lbfgs_calls == []
    prefixes = [call["prefix"] for call in rfo_calls]
    assert prefixes[0] == "preopt"
    dims = 2 if command == "scan2d" else 3
    point = "_".join(rf"d{n}_\d{{3}}" for n in range(1, dims + 1))
    assert len([p for p in prefixes if re.fullmatch(point, p)]) == n_points
    for call in rfo_calls:
        assert call["seeded"] is not None and call["seeded"].shape == (12, 12)
        assert call["kwargs"]["hessian_update"] == "bfgs"
        assert call["kwargs"]["thresh"] == "baker"
        assert call["kwargs"]["trust_max"] <= 0.004 * ANG2BOHR + 1e-12
    assert calc.hessian_calls == len(rfo_calls)


@pytest.mark.parametrize(
    ("command", "scan_lists"),
    [
        ("scan", "[(1,2,1.6)]"),
        ("scan2d", "[(1,2,1.500,1.508),(2,3,1.500,1.504)]"),
        ("scan3d", "[(1,2,1.500,1.508),(2,3,1.500,1.504),(3,4,1.500,1.504)]"),
    ],
)
def test_dry_run_accepts_hess_without_building_an_optimizer(
    tmp_path, monkeypatch, command, scan_lists,
) -> None:
    from mlmm.cli import cli as root_cli
    from mlmm.core.utils import DRY_RUN_COMPLETE_MESSAGE

    rfo_calls = _recording_rfo(monkeypatch)
    source = tmp_path / "four_carbons.pdb"
    _write_carbons(source, [(1.5 * i, 0.0, 0.0) for i in range(4)])
    result = CliRunner().invoke(root_cli, [
        command, "-i", str(source), "--parm", str(source), "--model-pdb", str(source),
        "-q", "0", "-m", "1", "--scan-lists", scan_lists, "--opt-mode", "hess",
        "--out-dir", str(tmp_path / "out"), "--dry-run",
    ])
    assert result.exit_code == 0, result.output + repr(result.exception)
    assert result.output.rstrip().splitlines()[-1] == DRY_RUN_COMPLETE_MESSAGE
    assert rfo_calls == []
