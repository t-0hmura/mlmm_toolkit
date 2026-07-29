"""Summary-status and MLIP provenance regression tests."""

from __future__ import annotations

from mlmm.workflows.all import (
    _derive_pipeline_status,
    _enrich_summary,
    _resolve_mlip_provenance,
    _validate_tsopt_result_payload,
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


def test_pipeline_status_rejects_zero_imaginary_modes() -> None:
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

    assert status == "partial"
    assert any("n_imag=0" in reason for reason in reasons)


@pytest.mark.parametrize("n_imag", [0, 2])
def test_all_stops_before_irc_for_wrong_saddle_order(n_imag: int) -> None:
    with pytest.raises(Exception, match="IRC was not started"):
        _validate_tsopt_result_payload(
            {"status": "not_converged", "n_imaginary_modes": n_imag},
            skip_final_freq=False,
        )


def test_all_allows_only_explicitly_unverified_skip() -> None:
    _validate_tsopt_result_payload(
        {"status": "unverified", "n_imaginary_modes": None},
        skip_final_freq=True,
    )
    with pytest.raises(Exception, match="IRC was not started"):
        _validate_tsopt_result_payload(
            {"status": "unverified", "n_imaginary_modes": None},
            skip_final_freq=False,
        )


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
    assert any(ref["method"] == "mlmm-toolkit" for ref in summary["references"])
    assert any(
        ref["method"] == "Growing String Method (GSM)"
        for ref in summary["references"]
    )
    assert any(
        ref["method"] == "RS-I-RFO" for ref in summary["references"]
    )
    assert any(
        ref["method"] == "Limited-memory BFGS (L-BFGS)"
        for ref in summary["references"]
    )


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
