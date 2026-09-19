import click
from click.core import ParameterSource
from click.testing import CliRunner
import pytest

from mlmm.cli.app import cli as root_cli
from mlmm.cli.decorators import (
    canonicalize_calculator_section,
    resolve_model_indices_setting,
)
from mlmm.core.utils import apply_yaml_overrides


def _option(command: str, name: str) -> click.Option:
    cmd = root_cli.get_command(click.Context(root_cli), command)
    assert cmd is not None
    return next(param for param in cmd.params if param.name == name)


class _Context:
    def __init__(self, **sources):
        self.sources = sources

    def get_parameter_source(self, name):
        return self.sources.get(name, ParameterSource.DEFAULT)


def test_yaml_model_indices_base_and_fixed_cli_base() -> None:
    spec, one_based = resolve_model_indices_setting(
        _Context(), {"calc": {"model_indices": [0, 2], "model_indices_base": 0}},
        None, True,
    )
    assert (spec, one_based) == ("0,2", False)

    spec, one_based = resolve_model_indices_setting(
        _Context(model_indices_str=ParameterSource.COMMANDLINE),
        {"calc": {"model_indices": [0, 2], "model_indices_base": 0}},
        "1,3", False,
    )
    assert (spec, one_based) == ("1,3", True)

    spec, one_based = resolve_model_indices_setting(
        _Context(
            model_indices_str=ParameterSource.COMMANDLINE,
            model_indices_one_based=ParameterSource.COMMANDLINE,
        ), {}, "0,2", False,
    )
    assert (spec, one_based) == ("0,2", False)


def test_calc_and_legacy_section_conflicts_are_rejected() -> None:
    merged = canonicalize_calculator_section(
        {"calc": {"backend": "uma"}, "mlmm": {"workers": 2}}
    )
    assert merged["calc"] == {"backend": "uma", "workers": 2}
    with pytest.raises(click.BadParameter, match="Conflicting YAML values"):
        canonicalize_calculator_section(
            {"calc": {"backend": "uma"}, "mlmm": {"backend": "orb"}}
        )


def test_calc_and_mlmm_sections_merge_in_runtime_override_path() -> None:
    target = {}
    apply_yaml_overrides(
        {"calc": {"backend": "uma"}, "mlmm": {"workers": 2}},
        [(target, (("calc",), ("mlmm",)))],
    )
    assert target == {"backend": "uma", "workers": 2}
    with pytest.raises(click.BadParameter, match="Conflicting YAML values"):
        apply_yaml_overrides(
            {"calc": {"backend": "uma"}, "mlmm": {"backend": "orb"}},
            [(target, (("calc",), ("mlmm",)))],
        )


def test_canonical_names_keep_compatibility_aliases() -> None:
    assert _option("all", "workers").opts == ["--uma-workers", "--workers"]
    assert _option("all", "parm7_override").opts == ["--parm7", "--parm"]
    assert _option("all", "dft_func_basis").opts == ["--func-basis", "--dft-func-basis"]
    assert _option("tsopt", "hess_cutoff").opts == ["--hessian-cutoff", "--radius-hessian", "--hess-cutoff"]
    assert _option("tsopt", "model_indices_one_based").hidden is True


def test_conflicting_alias_values_are_rejected_before_execution() -> None:
    result = CliRunner().invoke(
        root_cli,
        ["dft", "--scf-tol", "1e-8", "--conv-tol", "1e-7"],
    )
    assert result.exit_code == 2
    assert "Conflicting values were supplied through aliases" in result.output
