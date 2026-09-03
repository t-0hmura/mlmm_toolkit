import click
import pytest

from mlmm.workflows.all import (
    _derive_pipeline_status,
    _pipeline_aggregate_truth,
    _validate_postprocessing_dependencies,
)


@pytest.mark.parametrize(
    ("do_thermo", "do_dft"),
    [(True, False), (False, True), (True, True)],
)
def test_all_rejects_ts_labeled_postprocessing_without_tsopt(
    do_thermo, do_dft
):
    with pytest.raises(click.UsageError, match="require `--tsopt`"):
        _validate_postprocessing_dependencies(
            do_tsopt=False,
            do_thermo=do_thermo,
            do_dft=do_dft,
        )


@pytest.mark.parametrize(
    ("do_tsopt", "do_thermo", "do_dft"),
    [(False, False, False), (True, False, False), (True, True, True)],
)
def test_all_accepts_scientifically_defined_postprocessing_combinations(
    do_tsopt, do_thermo, do_dft
):
    _validate_postprocessing_dependencies(
        do_tsopt=do_tsopt,
        do_thermo=do_thermo,
        do_dft=do_dft,
    )


def test_missing_requested_post_segment_is_partial_and_unobserved():
    summary = {
        "segments": [
            {"index": 1, "kind": "seg", "converged": True},
            {"index": 2, "kind": "seg", "converged": True},
        ],
        "energy_diagrams": [{"name": "MEP"}],
    }
    post_segments = [
        {
            "index": 1,
            # A real `all --tsopt` record always carries the TS decision; the
            # aggregate fails closed without it.
            "tsopt": {"continue_irc": True},
            "mlip": {},
            "irc_traj": "seg_01/irc.trj",
            "irc": {"usable": True, "traj": "seg_01/irc.trj"},
            "endpoint_assignment": {"connectivity_validated": True},
            "endpoint_opt": {
                "reactant_converged": True,
                "product_converged": True,
                "connectivity_validated": True,
            },
            "ts_imag": {"n_imag": 1},
            "gibbs_mlip": {},
        }
    ]
    config = {"tsopt": True, "thermo": True, "dft": False}

    status, reasons = _derive_pipeline_status(
        summary,
        post_segments=post_segments,
        config=config,
    )
    truth = _pipeline_aggregate_truth(
        summary,
        post_segments=post_segments,
        config=config,
        legacy_status=status,
        legacy_reasons=reasons,
    )

    assert status == "partial"
    assert any("segment 2" in reason and "record is missing" in reason for reason in reasons)
    assert truth.scientific_status == "partial"
    assert truth.expected_item_ids == ("segment_1", "segment_2")
    assert truth.observed_item_ids == ("segment_1",)


def test_explicit_no_change_segment_does_not_require_postprocessing() -> None:
    summary = {
        "segments": [
            {
                "index": 1,
                "kind": "seg",
                "converged": True,
                "bond_changes": "(no covalent changes detected)",
            }
        ],
        "energy_diagrams": [{"name": "MEP"}],
    }
    config = {"tsopt": True, "thermo": False, "dft": False}

    status, reasons = _derive_pipeline_status(
        summary, post_segments=[], config=config
    )
    truth = _pipeline_aggregate_truth(
        summary,
        post_segments=[],
        config=config,
        legacy_status=status,
        legacy_reasons=reasons,
    )

    # The no-change segment still requires no per-segment record, which is what
    # this test guards. The headline verdict now says what actually happened:
    # TSOPT was requested and never ran, because nothing was reactive.
    assert status == "partial"
    assert any("no reactive segment" in reason for reason in reasons)
    assert truth.scientific_status == "partial"
    assert truth.expected_item_ids == ()


def test_failed_segment_dft_is_partial() -> None:
    summary = {
        "segments": [{"index": 1, "kind": "seg", "converged": True}],
        "energy_diagrams": [{"name": "MEP"}],
    }
    post = [{
        "index": 1,
        "mlip": {},
        "irc_traj": "seg_01/irc.trj",
        "ts_imag": {"n_imag": 1},
        "dft": {"status": "failed", "failed_states": ["TS"]},
    }]

    status, reasons = _derive_pipeline_status(
        summary,
        post_segments=post,
        config={"tsopt": True, "thermo": False, "dft": True},
    )

    assert status == "partial"
    assert "segment 1: DFT failed (TS)" in reasons
