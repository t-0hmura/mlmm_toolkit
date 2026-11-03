"""Pins the explicit FAIR-Chem InferenceSettings behind each UMA precision.

fairchem-core 2.24 turned its named default into a TF32 mode, so every UMA
route builds its settings from fields: fp32 must stay IEEE FP32 whatever the
installed fairchem release names its presets, and tf32 is the only TF32 route.
Every route except the analytical Hessian runs the fast path (fp64 included),
compiled only when --compile asks for it.
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


@pytest.mark.parametrize("compiled", [False, True], ids=["uncompiled", "compiled"])
@pytest.mark.parametrize(
    ("precision", "dtype"), [("fp32", "float32"), ("fp64", "float64")],
)
def test_fp32_and_fp64_run_the_ieee_fast_path(settings_cls, precision, dtype, compiled) -> None:
    s = mlmm_calc._uma_inference_settings(
        precision, analytical_hessian=False, parallel=False, compile_model=compiled
    )
    assert s == settings_cls(
        tf32=False, activation_checkpointing=False, merge_mole=True, compile=compiled,
        external_graph_gen=False, internal_graph_gen_version=2,
        base_precision_dtype=dtype,
    )


def test_fast_path_is_uncompiled_by_default(settings_cls) -> None:
    for precision in ("tf32", "fp32", "fp64"):
        s = mlmm_calc._uma_inference_settings(precision, analytical_hessian=False, parallel=False)
        assert s.compile is False and s.merge_mole and not s.activation_checkpointing


@pytest.mark.parametrize("compiled", [False, True], ids=["uncompiled", "compiled"])
def test_tf32_enables_tf32_and_nitro_on_one_worker(settings_cls, compiled) -> None:
    s = mlmm_calc._uma_inference_settings(
        "tf32", analytical_hessian=False, parallel=False, compile_model=compiled
    )
    assert s == settings_cls(
        tf32=True, activation_checkpointing=False, merge_mole=True, compile=compiled,
        external_graph_gen=False, internal_graph_gen_version=2, nitro=True,
    )


@pytest.mark.usefixtures("settings_cls")
def test_tf32_on_workers_skips_nitro() -> None:
    s = mlmm_calc._uma_inference_settings("tf32", analytical_hessian=False, parallel=True)
    assert s.tf32 and not s.nitro


def test_tf32_without_nitro_field(monkeypatch) -> None:
    monkeypatch.setattr(mlmm_calc, "_UMAInferenceSettings", _SettingsWithoutNitro)
    s = mlmm_calc._uma_inference_settings("tf32", analytical_hessian=False, parallel=False)
    assert s.tf32 and not hasattr(s, "nitro")


@pytest.mark.parametrize("compiled", [False, True], ids=["uncompiled", "compile-requested"])
@pytest.mark.parametrize(
    ("precision", "dtype", "tf32"),
    [("fp64", "float64", False), ("fp32", "float32", False), ("tf32", "float32", True)],
)
def test_analytical_route_keeps_fairchem_defaults_uncompiled(
    settings_cls, precision, dtype, tf32, compiled
) -> None:
    s = mlmm_calc._uma_inference_settings(
        precision, analytical_hessian=True, parallel=False, compile_model=compiled
    )
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
