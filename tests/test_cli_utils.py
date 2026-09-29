"""Unit tests for mlmm.cli.decorators."""

import sys

import pytest
import yaml

# mlmm.cli.decorators uses ``from __future__ import annotations`` (Python 3.7+)
# and other modern features.  Skip the entire module on older interpreters.
pytestmark = pytest.mark.skipif(
    sys.version_info < (3, 11),
    reason="mlmm.cli.decorators requires Python >= 3.11",
)


@pytest.fixture
def _parse_bool():
    from mlmm.cli.decorators import parse_bool
    return parse_bool


class TestParseBool:
    @pytest.mark.parametrize("value", ["true", "True", "TRUE", "1", "yes", "y", "t"])
    def test_true_values(self, value, _parse_bool):
        assert _parse_bool(value) is True

    @pytest.mark.parametrize("value", ["false", "False", "FALSE", "0", "no", "n", "f"])
    def test_false_values(self, value, _parse_bool):
        assert _parse_bool(value) is False

    def test_none_raises(self, _parse_bool):
        with pytest.raises(ValueError, match="None"):
            _parse_bool(None)

    def test_invalid_raises(self, _parse_bool):
        with pytest.raises(ValueError, match="Invalid boolean"):
            _parse_bool("maybe")


@pytest.fixture
def _load_merged_yaml_cfg():
    from mlmm.cli.decorators import load_merged_yaml_cfg
    return load_merged_yaml_cfg


class TestLoadMergedYamlCfg:
    def test_merge_keeps_source_layers_and_effective_tree_independent(
        self, tmp_path, _load_merged_yaml_cfg
    ):
        cfg_file = tmp_path / "config.yaml"
        cfg_file.write_text(
            yaml.dump({"opt": {"max_cycles": 100, "nested": {"left": 1}}})
        )
        ovr_file = tmp_path / "override.yaml"
        ovr_file.write_text(
            yaml.dump({"opt": {"max_cycles": 50, "nested": {"right": 2}}})
        )

        merged, cfg, ovr = _load_merged_yaml_cfg(cfg_file, ovr_file)

        assert cfg["opt"] == {"max_cycles": 100, "nested": {"left": 1}}
        assert ovr["opt"] == {"max_cycles": 50, "nested": {"right": 2}}
        assert merged["opt"] == {
            "max_cycles": 50,
            "nested": {"left": 1, "right": 2},
        }

        merged["opt"]["nested"]["left"] = 99
        assert cfg["opt"]["nested"]["left"] == 1
        assert "left" not in ovr["opt"]["nested"]

    def test_unknown_top_level_section_warns(self, tmp_path, capsys, _load_merged_yaml_cfg):
        cfg_file = tmp_path / "config.yaml"
        cfg_file.write_text("clac: {charge: 0}\n1: {}\nopt: {max_cycles: 10}\n")
        _load_merged_yaml_cfg(cfg_file, None)
        err = capsys.readouterr().err
        assert "YAML section(s) 1, clac are not recognized" in err

    def test_shared_config_sections_do_not_warn(self, tmp_path, capsys, _load_merged_yaml_cfg):
        sections = (
            "geom", "mlmm", "opt", "lbfgs", "rfo", "rsirfo", "stopt", "gs", "dmf",
            "sp", "freq", "thermo", "microiter", "hessian_dimer", "irc", "dft", "bond",
            "bias", "search",
        )
        cfg_file = tmp_path / "config.yaml"
        cfg_file.write_text(yaml.dump({name: {} for name in sections}))
        _load_merged_yaml_cfg(cfg_file, None)
        assert capsys.readouterr().err == ""
