"""Cartesian TS defaults preserve explicit YAML meanings and reach CLI owners."""

from copy import deepcopy

import pytest

from mlmm.core.defaults import GEOM_KW_DEFAULT, OPT_BASE_KW, RSIRFO_KW
from mlmm.core.utils import apply_yaml_overrides
from mlmm.workflows import tsopt


def resolve(base=None, override=None, *, kind="rsprfo", geom=None):
    base, override = deepcopy(base or {}), deepcopy(override or {})
    raw_before = deepcopy((base, override))
    opt, rs = dict(OPT_BASE_KW), deepcopy(RSIRFO_KW)
    geometry = deepcopy(GEOM_KW_DEFAULT)
    for layer in (base, override):
        apply_yaml_overrides(layer, [(opt, (("opt",),)), (rs, (("rsirfo",),)),
                                     (geometry, (("geom",),))])
    if geom:
        geometry.update(geom)
    assert (base, override) == raw_before
    return tsopt._build_rsirfo_kwargs(rs, max_cycles=1, out_dir=".", mode=kind, opt_cfg=opt)


def assert_baseline(kwargs):
    assert kwargs.get("trust_norm", "l2") == "l2"
    assert kwargs["trust_radius"] == kwargs["trust_max"] == 0.1
    assert kwargs["trust_min"] == 1e-4


@pytest.mark.parametrize("coord", ["cart", "cartesian"])
def test_unconfigured_cartesian_defaults_preserve_baseline(coord):
    before = deepcopy(RSIRFO_KW)
    kwargs = resolve(geom={"coord_type": coord})
    assert_baseline(kwargs)
    assert kwargs["hessian_update"] == "bofill"
    assert RSIRFO_KW == before  # other TS owners still share the old dictionary


@pytest.mark.parametrize("section", ["opt", "rsirfo"])
@pytest.mark.parametrize("layer", ["base", "override"])
def test_explicit_old_values_survive_even_when_equal_to_defaults(section, layer):
    config = {section: {"trust_norm": "l2", "trust_radius": 0.1,
                        "trust_min": 1e-4, "trust_max": 0.1,
                        "hessian_update": "bofill"}}
    kwargs = resolve(**{layer: config})
    assert kwargs["trust_norm"] == "l2"
    assert kwargs["trust_radius"] == kwargs["trust_max"] == 0.1
    assert kwargs["trust_min"] == 1e-4
    assert kwargs["hessian_update"] == "bofill"


@pytest.mark.parametrize("section", ["opt", "rsirfo"])
@pytest.mark.parametrize("key,value", [
    ("trust_radius", 0.1), ("trust_min", 1e-4), ("trust_max", 0.1),
    ("trust_radius", 0.025),
])
def test_radius_alone_keeps_legacy_global_l2_meaning(section, key, value):
    kwargs = resolve({section: {key: value}})
    assert kwargs.get("trust_norm", "l2") == "l2"
    assert kwargs[key] == value
    assert kwargs["trust_radius"] == (value if key == "trust_radius" else 0.1)
    assert kwargs["trust_max"] == 0.1
    assert kwargs["hessian_update"] == "bofill"


@pytest.mark.parametrize("section", ["opt", "rsirfo"])
def test_explicit_l2_keeps_baseline_radius_and_update(section):
    kwargs = resolve({section: {"trust_norm": "l2"}})
    assert kwargs["trust_norm"] == "l2"
    assert kwargs["trust_radius"] == kwargs["trust_max"] == 0.1
    assert kwargs["hessian_update"] == "bofill"


@pytest.mark.parametrize("section", ["opt", "rsirfo"])
def test_explicit_bofill_keeps_baseline_trust(section):
    kwargs = resolve({section: {"hessian_update": "bofill"}})
    assert_baseline(kwargs)
    assert kwargs["hessian_update"] == "bofill"


@pytest.mark.parametrize("section", ["opt", "rsirfo"])
def test_explicit_max_atom_preserves_default_and_explicit_radii(section):
    kwargs = resolve({section: {"trust_norm": "max_atom"}})
    assert kwargs["trust_norm"] == "max_atom"
    assert kwargs["trust_radius"] == kwargs["trust_max"] == 0.1
    kwargs = resolve({section: {"trust_norm": "max_atom", "trust_radius": 0.025,
                               "trust_min": 2e-4, "trust_max": 0.075}})
    assert kwargs["trust_norm"] == "max_atom"
    assert kwargs["trust_radius"] == 0.025
    assert kwargs["trust_min"] == 2e-4
    assert kwargs["trust_max"] == 0.075


def test_override_layer_and_ignored_section_keep_explicit_provenance():
    base = {"rsirfo": {"trust_radius": 0.025, "hessian_update": "bofill"}}
    overridden = resolve(base, {"rsirfo": {"trust_radius": 0.1}})
    assert overridden["trust_radius"] == 0.1
    assert overridden.get("trust_norm", "l2") == "l2"
    assert overridden["hessian_update"] == "bofill"
    ignored = resolve(base, {"rsirfo": None})
    assert ignored["trust_radius"] == 0.025
    assert ignored.get("trust_norm", "l2") == "l2"
    assert ignored["hessian_update"] == "bofill"


def test_explicit_atomic_override_preserves_base_radius():
    kwargs = resolve({"rsirfo": {"trust_radius": 0.1}},
                     {"rsirfo": {"trust_norm": "max_atom"}})
    assert kwargs["trust_norm"] == "max_atom"
    assert kwargs["trust_radius"] == 0.1
    assert kwargs["trust_max"] == 0.1


@pytest.mark.parametrize("kind,geom,extra", [
    ("rsirfo", {}, {}), ("trim", {}, {}), ("dimer", {}, {}),
    ("rsprfo", {"coord_type": "dlc"}, {}),
    ("rsprfo", {"coord_type": "hdlc"}, {}),
    ("rsprfo", {"coord_type": "mwcartesian"}, {}),
    ("rsprfo", {"coord_type": "cartesian", "coord_kwargs": {"mass_weighted": True}}, {}),
    ("rsprfo", {}, {"opt": {"weighted_trust": True}}),
    ("rsprfo", {}, {"rsirfo": {"weighted_trust": True}}),
])
def test_other_methods_coordinates_and_weighted_trust_keep_existing_defaults(kind, geom, extra):
    kwargs = resolve(extra, kind=kind, geom=geom)
    assert kwargs.get("trust_norm", "l2") == "l2"
    assert kwargs["trust_radius"] == kwargs["trust_max"] == 0.1
    assert kwargs["hessian_update"] == "bofill"


def test_existing_opt_rsirfo_priority_is_retained():
    kwargs = resolve({"opt": {"trust_norm": "max_atom", "weighted_trust": False},
                      "rsirfo": {"trust_norm": "l2", "weighted_trust": True}})
    assert kwargs["trust_norm"] == "max_atom"
    assert kwargs["weighted_trust"] is False
    assert kwargs["trust_radius"] == 0.1
    assert kwargs["hessian_update"] == "bofill"


@pytest.fixture
def capture_cli(monkeypatch, tmp_path):
    import json
    from types import SimpleNamespace
    import numpy as np
    import torch
    from click.testing import CliRunner
    from pysisyphus.Geometry import Geometry

    class ConstructorReached(RuntimeError):
        pass

    captured = {}

    class CaptureOptimizer:
        def __init__(self, geometry, **kwargs):
            captured.update(kwargs)
            captured["coord_type"] = geometry.coord_type
            raise ConstructorReached("ConstructorReached: no model evaluation requested.")

    def load_geometry(_path, **kwargs):
        return Geometry(("O", "H", "H"), np.array([0., 0., 0., 1.8, 0., 0., -.4, 1.7, 0.]),
                        coord_type=kwargs.get("coord_type", "cart"))

    source, parm = tmp_path / "input.pdb", tmp_path / "input.parm7"
    source.write_text("Prepared-structure parsing is replaced for constructor capture.\n")
    parm.write_text("No MM engine or topology parser is constructed.\n")
    prepared = SimpleNamespace(original_path=source, source_path=source,
                               geom_path=source, cleanup=lambda: None)
    monkeypatch.chdir(tmp_path)
    monkeypatch.setattr(tsopt, "prepare_input_structure", lambda *_a, **_k: prepared)
    monkeypatch.setattr(tsopt, "resolve_charge_spin_or_raise", lambda *_a, **_k: (0, 1))
    monkeypatch.setattr(tsopt, "resolve_ml_layer_assignment", lambda **_k: (source, None))
    monkeypatch.setattr(tsopt, "geom_loader", load_geometry)
    monkeypatch.setattr(tsopt, "mlmm", lambda **_k:
                        SimpleNamespace(core=SimpleNamespace(hess_active_atoms=[0, 1, 2])))
    monkeypatch.setattr(tsopt, "_calc_full_hessian_torch", lambda *_a, **_k:
                        torch.eye(9, dtype=torch.float64))
    monkeypatch.setattr(tsopt, "_torch_device", lambda *_a: torch.device("cpu"))
    monkeypatch.setitem(tsopt.TSOPT_CLASS_MAP, "rsprfo", CaptureOptimizer)

    def capture_micro(geometry, _calc, rs_cfg, _lbfgs, opt_cfg, _micro, out_dir, **kwargs):
        # Capture the resolved dictionary at the real microdriver boundary.
        # Its existing kwargs builder is shared with the direct constructor.
        opts = tsopt._build_rsirfo_kwargs(
            rs_cfg, max_cycles=opt_cfg["max_cycles"], out_dir=out_dir,
            macro_thresh=kwargs.get("thresh") or rs_cfg["thresh"],
            mode=kwargs["mode"], opt_cfg=opt_cfg,
            dump=kwargs["dump"], reference_mode=kwargs.get("reference_mode"),
            flatten_enabled=kwargs["flatten_enabled"],
        )
        captured["micro_boundary"] = True
        CaptureOptimizer(geometry, **opts)

    monkeypatch.setattr(tsopt, "_run_microiter_tsopt", capture_micro)

    def invoke(config=None, *, mode=None, coord=None, micro=False):
        args = ["-i", str(source), "--parm", str(parm), "-q", "0", "-m", "1",
                "-o", str(tmp_path / "out"), "--no-flatten", "--no-dump",
                "--no-convert-files", "--microiter" if micro else "--no-microiter"]
        if config is not None:
            path = tmp_path / "config.yaml"
            path.write_text(json.dumps(config))
            args += ["--config", str(path)]
        if mode is not None:
            args += ["--opt-mode", mode]
        if coord is not None:
            args += ["--coord-type", coord]
        result = CliRunner().invoke(tsopt.cli, args)
        assert result.exit_code == 1 and "ConstructorReached" in result.output, result.output
        assert captured, result.output
        return captured

    return invoke


@pytest.mark.parametrize("mode", [None, "hess", "rsprfo"])
def test_real_default_entry_reaches_constructor(capture_cli, mode):
    kwargs = capture_cli(mode=mode)
    assert_baseline(kwargs)
    assert kwargs["hessian_update"] == "bofill"


@pytest.mark.parametrize("section", ["opt", "rsirfo"])
def test_real_entry_preserves_explicit_legacy_choices(capture_cli, section):
    kwargs = capture_cli({section: {"trust_norm": "l2", "trust_radius": 0.1,
                                   "trust_max": 0.1, "hessian_update": "bofill"}})
    assert kwargs["trust_norm"] == "l2"
    assert kwargs["trust_radius"] == kwargs["trust_max"] == 0.1
    assert kwargs["hessian_update"] == "bofill"


def test_real_entry_uses_final_cli_coordinate_choice(capture_cli):
    kwargs = capture_cli({"geom": {"coord_type": "dlc"}}, coord="cart")
    assert kwargs["coord_type"] == "cart"
    assert_baseline(kwargs)


@pytest.mark.parametrize("section", ["opt", "rsirfo"])
def test_real_entry_radius_alone_keeps_l2(capture_cli, section):
    kwargs = capture_cli({section: {"trust_radius": 0.1}})
    assert kwargs.get("trust_norm", "l2") == "l2"
    assert kwargs["trust_radius"] == kwargs["trust_max"] == 0.1
    assert kwargs["hessian_update"] == "bofill"


def test_real_default_entry_reaches_microdriver(capture_cli):
    kwargs = capture_cli(micro=True)
    assert kwargs["micro_boundary"] is True
    assert_baseline(kwargs)
    assert kwargs["hessian_update"] == "bofill"


@pytest.mark.parametrize("section", ["opt", "rsirfo"])
def test_real_entry_preserves_explicit_atomic_norm(capture_cli, section):
    kwargs = capture_cli({section: {"trust_norm": "max_atom", "trust_radius": 0.025,
                                   "trust_max": 0.075, "hessian_update": "ts_bfgs"}})
    assert kwargs["trust_norm"] == "max_atom"
    assert kwargs["trust_radius"] == 0.025
    assert kwargs["trust_max"] == 0.075
    assert kwargs["hessian_update"] == "ts_bfgs"
