"""Pins the explicit FAIR-Chem InferenceSettings behind each UMA precision.

fairchem-core 2.24 turned its named default into a TF32 mode, so every UMA
route builds its settings from fields: fp32 must stay IEEE FP32 whatever the
installed fairchem release names its presets, and tf32 is the only TF32 route.
"""

from __future__ import annotations

from dataclasses import dataclass
from types import SimpleNamespace

import pytest

import mlmm.backends.mlmm_calc as mlmm_calc


@dataclass
class _Settings:
    tf32: bool = False
    activation_checkpointing: bool = True
    merge_mole: bool = False
    compile: bool = False
    external_graph_gen: bool = False
    internal_graph_gen_version: int = 2
    base_precision_dtype: str = "float32"
    nitro: bool = False


@dataclass
class _SettingsWithoutNitro:
    tf32: bool = False
    activation_checkpointing: bool = True
    merge_mole: bool = False
    compile: bool = False
    external_graph_gen: bool = False
    internal_graph_gen_version: int = 2
    base_precision_dtype: str = "float32"


@pytest.fixture
def settings_cls(monkeypatch):
    monkeypatch.setattr(mlmm_calc, "_UMAInferenceSettings", _Settings)
    return _Settings


def test_fp32_is_the_ieee_fast_path(settings_cls) -> None:
    s = mlmm_calc._uma_inference_settings("fp32", analytical_hessian=False, parallel=False)
    assert s == settings_cls(
        tf32=False, activation_checkpointing=False, merge_mole=True, compile=True,
        external_graph_gen=False, internal_graph_gen_version=2,
    )


@pytest.mark.usefixtures("settings_cls")
def test_tf32_enables_tf32_and_nitro_on_one_worker() -> None:
    s = mlmm_calc._uma_inference_settings("tf32", analytical_hessian=False, parallel=False)
    assert s.tf32 and s.nitro and s.merge_mole and s.compile


@pytest.mark.usefixtures("settings_cls")
def test_tf32_on_workers_skips_nitro() -> None:
    s = mlmm_calc._uma_inference_settings("tf32", analytical_hessian=False, parallel=True)
    assert s.tf32 and not s.nitro


def test_tf32_without_nitro_field(monkeypatch) -> None:
    monkeypatch.setattr(mlmm_calc, "_UMAInferenceSettings", _SettingsWithoutNitro)
    s = mlmm_calc._uma_inference_settings("tf32", analytical_hessian=False, parallel=False)
    assert s.tf32 and not hasattr(s, "nitro")


@pytest.mark.parametrize(
    ("precision", "dtype", "tf32"),
    [("fp64", "float64", False), ("fp32", "float32", False), ("tf32", "float32", True)],
)
def test_analytical_and_fp64_routes_are_not_compiled(settings_cls, precision, dtype, tf32) -> None:
    analytical = precision != "fp64"
    s = mlmm_calc._uma_inference_settings(precision, analytical_hessian=analytical, parallel=False)
    assert s == settings_cls(tf32=tf32, compile=False, base_precision_dtype=dtype)


@pytest.mark.parametrize(("precision", "loaded"), [("fp32", True), ("tf32", False), ("fp64", True)])
def test_verify_rejects_a_predictor_in_the_other_mode(precision, loaded) -> None:
    predictor = SimpleNamespace(inference_settings=SimpleNamespace(tf32=loaded))
    with pytest.raises(RuntimeError, match="tf32"):
        mlmm_calc._verify_uma_precision(predictor, precision)


@pytest.mark.parametrize(("precision", "loaded"), [("fp32", False), ("tf32", True), ("fp64", False)])
def test_verify_accepts_the_requested_mode(precision, loaded) -> None:
    predictor = SimpleNamespace(inference_settings=SimpleNamespace(tf32=loaded))
    mlmm_calc._verify_uma_precision(predictor, precision)


def test_cli_precision_choice_lists_tf32() -> None:
    import click

    from mlmm.cli.common_options import add_precision_option

    @click.command()
    @add_precision_option()
    def cmd(precision):
        click.echo(precision)

    choices = next(p for p in cmd.params if p.name == "precision").type.choices
    assert list(choices) == ["tf32", "fp32", "fp64"]
