"""Hessian TS optimizer construction in the microiteration path."""

from __future__ import annotations

import sys
from copy import deepcopy
from pathlib import Path

import numpy as np
import pytest
from click.testing import CliRunner

pytestmark = pytest.mark.skipif(
    sys.version_info < (3, 11),
    reason="mlmm requires Python >= 3.11",
)


@pytest.fixture(scope="module")
def tsopt_mod():
    return pytest.importorskip("mlmm.workflows.tsopt")


def _tiny_geom():
    """A minimal 3-atom geometry (bohr) — enough to construct an optimizer."""
    from pysisyphus.Geometry import Geometry

    ang2bohr = 1.0 / 0.529177210903
    coords = np.array([0, 0, 0, 0.96, 0, 0, -0.24, 0.93, 0], dtype=float) * ang2bohr
    return Geometry(atoms=("O", "H", "H"), coords=coords)


def test_tsopt_class_map_covers_all_hessian_modes(tsopt_mod):
    assert set(tsopt_mod.TSOPT_CLASS_MAP) == {"rsirfo", "rsprfo", "trim"}


@pytest.mark.parametrize("mode", ["rsirfo", "rsprfo", "trim"])
def test_microiter_macro_optimizer_builds(tsopt_mod, mode, tmp_path):
    from mlmm.core.defaults import RSIRFO_KW

    kw = tsopt_mod._build_rsirfo_kwargs(
        dict(RSIRFO_KW),
        max_cycles=1,
        out_dir=tmp_path,
        macro_thresh="baker",
        mode=mode,
    )
    if mode == "rsirfo":
        assert "min_line_search" not in kw
        assert "max_line_search" not in kw
    elif mode == "rsprfo":
        assert kw["min_line_search"] is False
        assert kw["max_line_search"] is False
    else:
        assert "min_line_search" not in kw
        assert "max_line_search" not in kw

    opt = tsopt_mod.TSOPT_CLASS_MAP[mode](_tiny_geom(), **kw)
    assert type(opt).__name__ in ("RSIRFOptimizer", "RSPRFOptimizer", "TRIM")


def test_rsprfo_honors_explicit_line_search_values(tsopt_mod, tmp_path):
    kw = tsopt_mod._build_rsirfo_kwargs(
        {"min_line_search": True, "max_line_search": True},
        max_cycles=1,
        out_dir=tmp_path,
        mode="rsprfo",
    )

    assert kw["min_line_search"] is True
    assert kw["max_line_search"] is True


def test_hessian_dimer_defaults_are_isolated_and_bounded_by_runner(tsopt_mod):
    assert tsopt_mod.hessian_dimer_KW["dimer"]["write_orientations"] is False
    assert "max_cycles" not in tsopt_mod.hessian_dimer_KW["lbfgs"]


@pytest.mark.parametrize(
    ("yaml_section", "yaml_out", "cli_out", "expected"),
    [
        (None, None, None, "result_tsopt"),
        ("opt", "yaml_ts", None, "yaml_ts"),
        ("rsirfo", "rsirfo_ts", None, "rsirfo_ts"),
        ("opt", "yaml_ts", "cli_ts", "cli_ts"),
    ],
)
def test_tsopt_output_dir_precedence(
    tsopt_mod, monkeypatch, tmp_path, yaml_section, yaml_out, cli_out, expected,
):
    source_repo = Path(__file__).resolve().parents[1]
    smoke = source_repo / "tests" / "smoke"
    source = smoke / "p_complex_layered.pdb"
    parm = smoke / "p_complex.parm7"
    if not source.is_file() or not parm.is_file():
        pytest.skip("smoke inputs are not present")

    monkeypatch.chdir(tmp_path)
    args = [
        "-i", str(source), "--parm", str(parm), "-q", "-1", "-m", "1",
        "--opt-mode", "hess", "--dry-run",
    ]
    if yaml_out is not None:
        config = tmp_path / "tsopt.yaml"
        config.write_text(
            f"{yaml_section}:\n  out_dir: {yaml_out}\n", encoding="utf-8"
        )
        args.extend(["--config", str(config)])
    if cli_out is not None:
        args.extend(["--out-dir", cli_out])

    dry_run_plan = {}

    def capture_pretty_block(title, content, **_kwargs):
        if title == "dry_run_plan":
            dry_run_plan.update(content)
        return ""

    monkeypatch.setattr(tsopt_mod, "pretty_block", capture_pretty_block)
    result = CliRunner().invoke(tsopt_mod.cli, args)

    assert result.exit_code == 0, result.output
    assert dry_run_plan["output_dir"] == str(tmp_path / expected)


@pytest.mark.parametrize(
    ("mode", "config_text", "message"),
    [
        (
            "hess",
            "opt:\n  max_cycles: 10\nrsirfo:\n  max_cycles: 20\n",
            "opt.max_cycles and rsirfo.max_cycles conflict",
        ),
        (
            "hess",
            "opt:\n  out_dir: root-out\nrsirfo:\n  out_dir: downstream-out\n",
            "opt.out_dir and rsirfo.out_dir conflict",
        ),
        (
            "grad",
            "opt:\n  thresh: gau\nhessian_dimer:\n  thresh: baker\n",
            "opt.thresh and hessian_dimer.thresh conflict",
        ),
        (
            "grad",
            "hessian_dimer:\n  lbfgs:\n    max_cycles: 5\n",
            "hessian_dimer.lbfgs.max_cycles is not configurable",
        ),
    ],
)
def test_tsopt_rejects_ambiguous_shared_optimizer_config(
    tsopt_mod, monkeypatch, tmp_path, mode, config_text, message,
):
    source_repo = Path(__file__).resolve().parents[1]
    source = source_repo / "tests" / "smoke" / "p_complex_layered.pdb"
    parm = source_repo / "tests" / "smoke" / "p_complex.parm7"
    if not source.is_file() or not parm.is_file():
        pytest.skip("smoke inputs are not present")

    config = tmp_path / "tsopt.yaml"
    config.write_text(config_text, encoding="utf-8")
    monkeypatch.chdir(tmp_path)
    result = CliRunner().invoke(
        tsopt_mod.cli,
        [
            "-i", str(source), "--parm", str(parm), "-q", "-1", "-m", "1",
            "--opt-mode", mode, "--dry-run", "--config", str(config),
        ],
    )

    assert result.exit_code != 0
    assert message in result.output


def test_hessian_dimer_nested_yaml_does_not_leak_between_invocations(
    tsopt_mod, monkeypatch, tmp_path,
):
    source_repo = Path(__file__).resolve().parents[1]
    source = source_repo / "tests" / "smoke" / "p_complex_layered.pdb"
    parm = source_repo / "tests" / "smoke" / "p_complex.parm7"
    if not source.is_file() or not parm.is_file():
        pytest.skip("smoke inputs are not present")

    defaults_before = deepcopy(tsopt_mod.hessian_dimer_KW)
    config = tmp_path / "tsopt.yaml"
    config.write_text(
        "hessian_dimer:\n"
        "  dimer:\n"
        "    write_orientations: true\n"
        "  lbfgs:\n"
        "    line_search: false\n",
        encoding="utf-8",
    )
    common = [
        "-i", str(source), "--parm", str(parm), "-q", "-1", "-m", "1",
        "--opt-mode", "grad", "--dry-run",
    ]
    monkeypatch.chdir(tmp_path)

    first = CliRunner().invoke(tsopt_mod.cli, [*common, "--config", str(config)])
    second = CliRunner().invoke(tsopt_mod.cli, common)

    assert first.exit_code == 0, first.output
    assert second.exit_code == 0, second.output
    assert tsopt_mod.hessian_dimer_KW == defaults_before


def test_shared_optimizer_value_rejects_explicit_conflict(tsopt_mod):
    opt_cfg = {"max_cycles": 10}
    rsirfo_cfg = {"max_cycles": 20}
    with pytest.raises(tsopt_mod.click.BadParameter, match="opt.max_cycles"):
        tsopt_mod._resolve_shared_optimizer_value(
            opt_cfg,
            rsirfo_cfg,
            "max_cycles",
            opt_explicit=True,
            downstream_explicit=True,
            downstream_default=300,
            downstream_section="rsirfo",
        )


@pytest.mark.parametrize(
    "opt_explicit, downstream_explicit, expected",
    [(True, False, 10), (False, True, 20), (False, False, 300)],
)
def test_shared_optimizer_value_precedence(
    tsopt_mod, opt_explicit, downstream_explicit, expected,
):
    opt_cfg = {"max_cycles": 10}
    rsirfo_cfg = {"max_cycles": 20}
    tsopt_mod._resolve_shared_optimizer_value(
        opt_cfg,
        rsirfo_cfg,
        "max_cycles",
        opt_explicit=opt_explicit,
        downstream_explicit=downstream_explicit,
        downstream_default=300,
        downstream_section="rsirfo",
    )
    assert opt_cfg["max_cycles"] == expected
    assert rsirfo_cfg["max_cycles"] == expected


@pytest.mark.parametrize("mode", ["rsirfo", "trim"])
def test_unused_line_search_values_are_removed(tsopt_mod, mode, tmp_path):
    kw = tsopt_mod._build_rsirfo_kwargs(
        {"min_line_search": True, "max_line_search": True},
        max_cycles=1,
        out_dir=tmp_path,
        mode=mode,
    )

    assert "min_line_search" not in kw
    assert "max_line_search" not in kw


def test_hessian_ts_kwargs_require_one_root(tsopt_mod, tmp_path):
    with pytest.raises(tsopt_mod.click.BadParameter, match="exactly one root"):
        tsopt_mod._build_rsirfo_kwargs(
            {"roots": [0, 1]},
            max_cycles=1,
            out_dir=tmp_path,
            mode="rsprfo",
        )


def test_ts_kwargs_drop_ordinary_rfo_overlap_tracking(tsopt_mod, tmp_path):
    kw = tsopt_mod._build_rsirfo_kwargs(
        {"rfo_overlaps": True},
        max_cycles=1,
        out_dir=tmp_path,
        mode="rsprfo",
    )

    assert "rfo_overlaps" not in kw
