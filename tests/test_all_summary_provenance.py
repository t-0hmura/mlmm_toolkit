"""Summary-status and MLIP provenance regression tests."""

from __future__ import annotations

from mlmm.workflows.all import (
    _derive_pipeline_status,
    _enrich_summary,
    _resolve_mlip_provenance,
    _tsopt_continuation_decision,
)
import pytest


def test_resolved_mlip_provenance_uses_backend_defaults() -> None:
    assert _resolve_mlip_provenance(
        backend="orb",
        backend_model=None,
        calc_file=None,
        calc_factory="get_calculator",
    ) == ("orb", "orb_v3_conservative_omol", "fp64")
    assert _resolve_mlip_provenance(
        backend=None,
        backend_model=None,
        calc_file="/tmp/custom_calc.py",
        calc_factory="make_calc",
    ) == ("custom", "custom_calc.py:make_calc", None)


def test_resolved_mlip_provenance_uses_yaml_custom_calculator() -> None:
    assert _resolve_mlip_provenance(
        backend=None,
        backend_model=None,
        calc_file=None,
        calc_factory=None,
        merged_yaml_cfg={
            "calc": {"calc_file": "/tmp/custom.py", "calc_factory": "build"}
        },
    ) == ("custom", "custom.py:build", None)


def test_resolved_mlip_provenance_defaults_custom_factory() -> None:
    assert _resolve_mlip_provenance(
        backend=None,
        backend_model=None,
        calc_file="/tmp/custom.py",
        calc_factory=None,
    ) == ("custom", "custom.py:get_calculator", None)


@pytest.mark.parametrize(
    ("backend", "precision", "expected"),
    [
        ("uma", "fp64", "fp64"),
        ("orb", "fp32", "fp32"),
        ("mace", "fp32", "fp32"),
        ("mace", None, "fp64"),
    ],
)
def test_resolved_mlip_provenance_records_effective_precision(
    backend: str, precision: str | None, expected: str
) -> None:
    assert _resolve_mlip_provenance(
        backend=backend,
        backend_model=None,
        calc_file=None,
        calc_factory=None,
        precision=precision,
    )[2] == expected


def test_pipeline_status_keeps_frequency_count_diagnostic() -> None:
    status, reasons = _derive_pipeline_status(
        {"segments": [{"index": 1}], "energy_diagrams": [{"name": "MEP"}]},
        post_segments=[
            {
                "index": 1,
                "mlip": {"barrier_kcal": 10.0},
                "gibbs_mlip": {"barrier_kcal": 9.0},
                "irc_traj": "irc.xyz",
                "ts_imag": {"n_imag": 0},
            }
        ],
        config={"tsopt": True, "thermo": True, "dft": False},
    )

    assert status == "success"
    assert reasons == []


def test_all_stops_before_irc_for_zero_imaginary_modes() -> None:
    decision = _tsopt_continuation_decision(
        {
            "optimization_status": "converged",
            "hessian_status": "completed",
            "n_imaginary_modes": 0,
        },
        skip_final_freq=False,
    )
    assert decision["continue_irc"] is False
    assert decision["reason"] == "no_imaginary_reaction_mode"


def test_all_continues_diagnostic_irc_for_converged_higher_order_saddle() -> None:
    decision = _tsopt_continuation_decision(
        {
            "optimization_status": "converged",
            "hessian_status": "completed",
            "n_imaginary_modes": 2,
            "imaginary_frequencies_cm": [-400.0, -100.0],
            "reaction_mode_index": 1,
            "reaction_mode_frequency_cm": -100.0,
        },
        skip_final_freq=False,
    )
    assert decision["continue_irc"] is True
    assert decision["reason"] == "higher_order_saddle"


def test_all_stops_before_irc_when_terminal_frequency_is_skipped() -> None:
    decision = _tsopt_continuation_decision(
        {
            "optimization_status": "converged",
            "hessian_status": "skipped",
            "n_imaginary_modes": None,
        },
        skip_final_freq=True,
    )
    assert decision["continue_irc"] is False
    assert decision["reason"] == "terminal_hessian_explicitly_skipped"


def test_all_rejects_nonnegative_recorded_reaction_root() -> None:
    decision = _tsopt_continuation_decision(
        {
            "optimization_status": "converged",
            "hessian_status": "completed",
            "n_imaginary_modes": 2,
            "imaginary_frequencies_cm": [-400.0, -100.0],
            "reaction_mode_index": 9,
            "reaction_mode_frequency_cm": 25.0,
        },
        skip_final_freq=False,
    )
    assert decision["continue_irc"] is True
    assert decision["reaction_mode_index"] == 0
    assert decision["reaction_mode_frequency_cm"] == -400.0
    assert decision["reaction_mode_fallback"] is True

def test_enriched_rate_limit_uses_refined_barrier(tmp_path) -> None:
    summary = {
        "out_dir": str(tmp_path / "_work"),
        "segments": [{"index": 1, "kind": "seg", "barrier_kcal": 5.0}],
        "energy_diagrams": [
            {
                "name": "energy_diagram_MLIP_all",
                "energies_kcal": [0.0, 12.0, -1.0],
            }
        ],
    }
    post = [
        {
            "index": 1,
            "mlip": {"barrier_kcal": 12.0},
            "irc_traj": "irc.xyz",
            "tsopt": {"continue_irc": True},
            "irc": {"forward_status": "stopped"},
            "ts_imag": {"n_imag": 1},
            "endpoint_opt": {"reactant_converged": True},
        }
    ]
    _enrich_summary(
        summary,
        version="",
        pipeline_mode="path-opt",
        mlip_backend="orb",
        mlip_model="orb_v3_conservative_omol",
        mlip_precision="fp64",
        charge=0,
        spin=1,
        post_segments=post,
        config={
            "tsopt": True,
            "thermo": False,
            "dft": False,
            "path_opt_mode": "grad",
            "post_opt_mode": "hess",
            "ts_opt_mode": "hess",
            "endpoint_opt_mode": "grad",
            "mep_mode": "gsm",
        },
        out_dir=tmp_path,
    )

    assert summary["status"] == "success"
    assert summary["mlip_precision"] == "fp64"
    assert summary["rate_limiting_step"] == {
        "segment": 1,
        "barrier_kcal": 12.0,
        "method": "MLIP",
        "mep_barrier_kcal": 5.0,
    }
    assert summary["mlip_backend"] == "orb"
    assert summary["mlip_model"] == "orb_v3_conservative_omol"
    assert summary["mlip_model_label"] == "ORB-v3-conservative-OMol"
    assert any(ref["method"] == "mlmm-toolkit" for ref in summary["references"])
    assert any(ref["method"] == "Orb-v3" for ref in summary["references"])
    assert any(
        ref["method"] == "Growing String Method (GSM)"
        for ref in summary["references"]
    )
    assert any(
        ref["method"] == "RS-P-RFO" for ref in summary["references"]
    )
    assert any(
        ref["method"] == "Limited-memory BFGS (L-BFGS)"
        for ref in summary["references"]
    )


def test_enriched_references_only_add_omol25_for_omol(tmp_path) -> None:
    summary = {"segments": [], "energy_diagrams": []}
    _enrich_summary(
        summary,
        version="",
        pipeline_mode="path-opt",
        mlip_backend="uma",
        mlip_model="uma-s-1p2",
        charge=0,
        spin=1,
        calculator_config={
            "backend": "uma",
            "uma_model": "uma-s-1p2",
            "uma_task_name": "non-omol",
        },
        out_dir=tmp_path,
    )

    methods = [reference["method"] for reference in summary["references"]]
    assert summary["mlip_task"] == "non-omol"
    assert "UMA" in methods
    assert "OMol25" not in methods


def test_ts_only_summary_does_not_assign_reaction_direction(tmp_path) -> None:
    summary = {
        "out_dir": str(tmp_path),
        "segments": [{
            "index": 1,
            "kind": "tsopt",
            "barrier_from_endpoint_1_kcal": 8.0,
            "barrier_from_endpoint_2_kcal": 9.0,
        }],
        "energy_diagrams": [{
            "name": "energy_diagram_MLIP_all",
            "labels": ["E1", "TS", "E2"],
            "energies_kcal": [0.0, 8.0, -1.0],
        }],
    }
    _enrich_summary(
        summary,
        version="",
        pipeline_mode="tsopt-only",
        mlip_backend="orb",
        charge=0,
        spin=1,
        out_dir=tmp_path,
        config={"tsopt": False},
    )

    assert "rate_limiting_step" not in summary
    assert "overall_reaction_energy_kcal" not in summary
    assert "overall_reaction_energy_method" not in summary
