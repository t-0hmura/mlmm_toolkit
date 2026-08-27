"""Regression tests for summary_log payload handling and rendering."""

from __future__ import annotations

import json
import sys
import time
from pathlib import Path

import pytest

pytestmark = pytest.mark.skipif(
    sys.version_info < (3, 11),
    reason="mlmm requires Python >= 3.11",
)


def test_normalize_summary_payload_sets_defaults():
    from mlmm.io.summary import normalize_summary_payload

    payload = normalize_summary_payload({"pipeline_mode": "path-search"})
    assert payload["pipeline_mode"] == "path-search"
    assert payload["root_out_dir"] == "-"
    assert payload["path_module_dir"] == "-"
    assert payload["segments"] == []
    assert payload["energy_diagrams"] == []


def test_write_summary_log_accepts_empty_payload(tmp_path: Path):
    from mlmm.io.summary import write_summary_log

    dest = tmp_path / "summary.log"
    write_summary_log(dest, {})

    text = dest.read_text(encoding="utf-8")
    assert "mlmm summary.log" in text
    assert "missing payload keys replaced with defaults" in text


def test_write_summary_log_renders_segment_section(tmp_path: Path):
    from mlmm.io.summary import write_summary_log

    out_root = tmp_path / "run"
    out_root.mkdir(parents=True, exist_ok=True)

    payload = {
        "root_out_dir": str(out_root),
        "path_module_dir": "path_search",
        "pipeline_mode": "path-search",
        "segments": [
            {
                "index": 1,
                "tag": "seg_01",
                "kind": "seg",
                "barrier_kcal": 12.3,
                "delta_kcal": -1.2,
                "bond_changes": "Broken: C1-O1",
            }
        ],
        "energy_diagrams": [],
    }

    dest = out_root / "summary.log"
    write_summary_log(dest, payload)

    text = dest.read_text(encoding="utf-8")
    assert "[2] Segment-level MEP summary" in text
    assert "Segment 01 [seg]" in text
    assert "Broken: C1-O1" in text


def test_write_summary_log_marks_non_successful_results_and_precision(
    tmp_path: Path,
) -> None:
    from mlmm.io.summary import write_summary_log

    dest = tmp_path / "summary.log"
    write_summary_log(
        dest,
        {
            "root_out_dir": str(tmp_path),
            "path_module_dir": "path_search",
            "pipeline_mode": "path-search",
            "segments": [{"index": 1, "barrier_kcal": 12.3}],
            "energy_diagrams": [],
            "mlip_precision": "fp64",
            "execution_status": "completed",
            "scientific_status": "partial",
            "scientific_status_reasons": ["segment 1 did not converge"],
        },
    )

    text = dest.read_text(encoding="utf-8")
    assert "MLIP precision      : fp64" in text
    assert "Execution status    : completed" in text
    assert "Scientific status   : partial" in text
    assert "RESULT WARNING      : Segment 1 did not converge." in text
    assert "Status reason" not in text
    assert text.index("RESULT WARNING") < text.index("ΔE‡")


@pytest.mark.parametrize(
    ("reason", "expected"),
    [
        (
            "all:segment_2:mep_not_converged",
            "Segment 2: MEP optimization did not converge. Review the MEP trajectory and convergence log.",
        ),
        (
            "all:segment_2:endpoint_opt:endpoint_2_converged",
            "Segment 2: the endpoint 2 optimization did not converge or could not be confirmed. "
            "Review the endpoint structure and optimizer log.",
        ),
        (
            "all:segment_3:irc:irc:forward:not_converged;irc:backward:energy_invalid",
            "Segment 3: Forward IRC did not converge. Review its trajectory and IRC log. "
            "Backward IRC did not produce a valid energy profile. Review its trajectory and IRC log.",
        ),
        (
            "path-search:endpoint_hei;engine_nonconverged",
            "No reactive segment was identified and the path-search engine did not converge. "
            "Review the path and path-search log.",
        ),
        (
            "missing:segment_4",
            "Segment 4: the expected segment result is missing. Review the workflow outputs.",
        ),
        (
            "all:segment_2:tsopt:ts_optimization_not_converged",
            "Segment 2: TS optimization did not converge. Review the TS trajectory.",
        ),
        (
            "segment 2: TS imaginary-mode validation found n_imag=2, expected 1",
            "Segment 2: TS imaginary-mode validation found n_imag=2. Consider --flatten --refine-path.",
        ),
    ],
)
def test_format_result_warning_explains_priority_status_codes(reason, expected):
    from mlmm.io.summary import format_result_warning

    assert format_result_warning(reason) == expected


def test_format_result_warning_omits_already_active_recovery_flags():
    from mlmm.io.summary import format_result_warning

    reason = "segment 2: TS imaginary-mode validation found n_imag=2, expected 1"
    assert format_result_warning(reason, refine_path=True) == (
        "Segment 2: TS imaginary-mode validation found n_imag=2. "
        "Consider --flatten."
    )
    assert format_result_warning(reason, refine_path=True, flatten=True) == (
        "Segment 2: TS imaginary-mode validation found n_imag=2."
    )


def test_summary_log_tree_lists_only_current_run_paths(tmp_path):
    from mlmm.io.summary import write_summary_log

    current = tmp_path / "segments" / "seg_01" / "structures" / "ts.pdb"
    stale = tmp_path / "segments" / "seg_02" / "structures" / "old.pdb"
    stale_diagram = tmp_path / "irc_plot_all.png"
    for path in (current, stale, stale_diagram):
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text("x", encoding="utf-8")
    dest = tmp_path / "summary.log"

    write_summary_log(
        dest,
        {
            "root_out_dir": str(tmp_path),
            "current_output_paths": ["segments/seg_01/structures/ts.pdb"],
        },
    )

    text = dest.read_text(encoding="utf-8")
    assert "seg_01/" in text
    assert "ts.pdb" in text
    assert "seg_02" not in text
    assert "irc_plot_all.png" not in text


def test_summary_log_renders_mlmm_provenance_and_labels(tmp_path):
    from mlmm.io.summary import write_summary_log

    dest = tmp_path / "summary.log"
    write_summary_log(
        dest,
        {
            "root_out_dir": str(tmp_path),
            "layer_counts": {"ml": 12, "movable": 20, "frozen": 30},
            "post_segments": [
                {
                    "index": 1,
                    "ts_imag": {
                        "n_imag": 1,
                        "nu_imag_max_cm": -430.0,
                        "frequency_zero_cutoff_cm": 5.0,
                    },
                    "thermo_symmetry": {
                        "R": {
                            "symmetry_number": 1,
                            "symmetry_number_source": "automatic",
                        },
                    },
                    "mlip": {
                        "labels": ["R", "TS", "P"],
                        "energies_au": [-3.0, -2.9, -3.1],
                        "energies_kcal": [0.0, 62.75, -62.75],
                    },
                }
            ],
        },
    )

    text = dest.read_text(encoding="utf-8")
    assert "ML (B=0)" in text and "12 atoms" in text
    assert "Movable MM (B=10)" in text and "20 atoms" in text
    assert "Frozen MM (B=20)" in text and "30 atoms" in text
    assert "zero cutoff  : 5.00 cm^-1" in text
    assert "Thermo symmetry   : R=1 (automatic)" in text
    assert "ML/MM energies (TSOPT+IRC)" in text
    assert "MLIP energies (TSOPT+IRC)" not in text


def test_write_summary_log_publish_failure_preserves_previous_file(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    from mlmm.core import result_commit
    from mlmm.core.result_commit import ResultCommitError
    from mlmm.io.summary import write_summary_log

    dest = tmp_path / "summary.log"
    dest.write_bytes(b"previous complete summary\n")

    def fail_replace(_staged, _destination):
        raise OSError("injected")

    monkeypatch.setattr(result_commit, "_replace_exact", fail_replace)
    with pytest.raises(ResultCommitError, match="publish"):
        write_summary_log(dest, {})

    assert dest.read_bytes() == b"previous complete summary\n"


def test_write_summary_log_ts_only_separates_model_dft_from_composite_gibbs(
    tmp_path: Path,
):
    from mlmm.io.summary import write_summary_log

    dest = tmp_path / "summary.log"
    write_summary_log(
        dest,
        {
            "root_out_dir": str(tmp_path),
            "path_module_dir": "-",
            "pipeline_mode": "tsopt-only",
            "n_images": 5,
            "n_segments": 1,
            "mep": {"n_images": 0, "n_segments": 0},
            "segments": [{
                "index": 1,
                "tag": "seg_01",
                "kind": "tsopt",
                "barrier_from_endpoint_1_kcal": 8.0,
                "barrier_from_endpoint_2_kcal": 9.0,
            }],
            "post_segments": [{
                "index": 1,
                "tag": "seg_01",
                "kind": "tsopt",
                "dft": {
                    "barrier_from_endpoint_1_kcal": 7.8,
                    "barrier_from_endpoint_2_kcal": 8.2,
                },
                "gibbs_dft_mlip": {
                    "barrier_from_endpoint_1_kcal": 9.1,
                    "barrier_from_endpoint_2_kcal": 9.5,
                },
            }],
            "energy_diagrams": [],
        },
    )

    text = dest.read_text(encoding="utf-8")
    assert "Number of IRC frames : 5" in text
    assert "Number of segments   : 1" in text
    assert "chemically unassigned endpoint" in text
    assert "MEP ΔE" not in text
    assert "model-region DFT ΔE‡ E1->TS" in text
    assert "model-region DFT ΔE‡ E2->TS" in text
    assert "DFT//MLIP/MM ΔG‡ E1->TS" in text
    assert "DFT//MLIP/MM ΔG‡ E2->TS" in text
    assert "DFT//MLIP/MM ΔE" not in text


def test_method_citations_follow_resolved_methods_and_match_stdout(
    tmp_path: Path, capsys
):
    from mlmm.io.summary import (
        emit_method_citations,
        format_method_citations,
        method_references,
        write_summary_log,
    )

    payload = {
        "root_out_dir": str(tmp_path),
        "path_module_dir": "path_search",
        "pipeline_mode": "path-search",
        "mep_mode": "gsm",
        "path_opt_mode": "grad",
        "post_opt_mode": "hess",
        "ts_opt_mode": "hess",
        "endpoint_opt_mode": "hess",
        "post_segments": [
            {
                "tsopt": {"continue_irc": True},
                "irc": {"forward_converged": True},
                "endpoint_opt": {"reactant_converged": True},
                "ts_imag": {"n_imag": 1},
                "thermo_symmetry": {"R": {"symmetry_number": 1}},
            }
        ],
        "segments": [],
        "energy_diagrams": [],
    }
    dest = tmp_path / "summary.log"

    write_summary_log(dest, payload)
    lines = format_method_citations(payload)
    references = method_references(payload)
    emit_method_citations(payload)

    text = dest.read_text(encoding="utf-8")
    stdout = capsys.readouterr().out
    block = "\n".join(lines)
    assert block in text
    assert text.rstrip().endswith(block)
    # Same citations, destination-appropriate header: the log keeps its section
    # index, stdout heads the block like every other console section.
    assert stdout.splitlines()[0] == "====== Citations & References ======"
    assert stdout.splitlines()[1:] == lines[1:]
    assert lines[0] == "[6] Methods and citations"
    assert "mlmm-toolkit:" in block
    assert "Growing String Method (GSM)" in block
    assert "RFO / P-RFO" in block
    assert "RS-P-RFO" in block
    assert "quasi-RRHO thermochemistry" in block
    assert "Direct Max Flux (DMF)" not in block
    assert all(set(ref) == {"method", "citation", "doi"} for ref in references)
    assert len({ref["doi"] for ref in references}) == len(references)
    assert lines[1] == "Please cite the software and methods used:"
    cursor = 2
    previous_method = None
    method_index = 0
    for reference in references:
        if reference["method"] != previous_method:
            method_index += 1
            assert lines[cursor] == f"({method_index}) {reference['method']}:"
            cursor += 1
            previous_method = reference["method"]
        assert lines[cursor] == f"- {reference['citation']}"
        cursor += 1
    assert cursor == len(lines)


def test_method_citations_use_actual_path_and_post_stages() -> None:
    from mlmm.io.summary import format_method_citations

    path_only = {
        "pipeline_mode": "path-search",
        "mep_mode": "dmf",
        "path_opt_mode": "grad",
        "post_opt_mode": "hess",
        "post_segments": [],
    }
    ts_only = {
        **path_only,
        "ts_opt_mode": "hess",
        "post_segments": [{"tsopt": {}}],
    }
    irc_only = {
        **path_only,
        "post_segments": [{"irc": {}}],
    }
    endpoint_only = {
        **path_only,
        "endpoint_opt_mode": "hess",
        "post_segments": [{"endpoint_opt": {}}],
    }
    complete = {
        **path_only,
        "ts_opt_mode": "hess",
        "endpoint_opt_mode": "hess",
        "post_segments": [
            {
                "tsopt": {},
                "irc": {},
                "endpoint_opt": {},
                "ts_imag": {"n_imag": 1},
            }
        ],
    }

    path_text = "\n".join(format_method_citations(path_only))
    ts_text = "\n".join(format_method_citations(ts_only))
    irc_text = "\n".join(format_method_citations(irc_only))
    endpoint_text = "\n".join(format_method_citations(endpoint_only))
    complete_text = "\n".join(format_method_citations(complete))

    assert "Limited-memory BFGS (L-BFGS)" in path_text
    assert "RFO / P-RFO" not in path_text
    assert "RS-P-RFO" not in path_text
    assert "Euler predictor-corrector IRC" not in path_text
    assert "quasi-RRHO thermochemistry" not in path_text
    assert "RS-P-RFO" in ts_text
    assert "Euler predictor-corrector IRC" not in ts_text
    assert "Euler predictor-corrector IRC" in irc_text
    assert "RFO / P-RFO" not in irc_text
    assert "RFO / P-RFO" in endpoint_text
    assert "RS-P-RFO" not in endpoint_text
    assert "Euler predictor-corrector IRC" not in endpoint_text
    assert "RS-P-RFO" in complete_text
    assert "Euler predictor-corrector IRC" in complete_text
    assert "quasi-RRHO thermochemistry" not in complete_text


@pytest.mark.parametrize(
    ("backend", "model", "task", "expected"),
    [
        ("uma", "uma-s-1p2", "omol", ["UMA", "OMol25"]),
        ("uma", "uma-s-1p2", None, ["UMA", "OMol25"]),
        ("uma", "uma-s-1p2", "non-omol", ["UMA"]),
        (
            "orb", "orb_v3_conservative_omol", None,
            ["Orb-v3", "OMol25"],
        ),
        (
            "mace", "MACE-OMOL-0", None,
            ["MACE", "MACE", "OMol25"],
        ),
        ("aimnet2", "aimnet2", None, ["AIMNet2"]),
    ],
)
def test_method_citations_include_the_executed_mlip_model(
    backend: str, model: str, task, expected: list[str],
) -> None:
    from mlmm.io.summary import method_references

    references = method_references({
        "mlip_backend": backend,
        "mlip_model": model,
        "mlip_task": task,
    })
    mlip_references = [
        reference for reference in references
        if reference["method"] in {
            "UMA", "Orb-v3", "MACE", "OMol25", "AIMNet2"
        }
    ]
    assert [reference["method"] for reference in mlip_references] == expected
    for reference in mlip_references:
        assert "et al." not in reference["citation"]
        assert "https://doi.org/" in reference["citation"]
    if backend == "aimnet2":
        assert mlip_references[0]["doi"] == "10.1039/D4SC08572H"
        assert "AIMNet2: A Neural Network Potential" in mlip_references[0][
            "citation"
        ]


def test_mace_citations_share_one_heading_and_keep_both_papers() -> None:
    from mlmm.io.summary import format_method_citations

    block = "\n".join(format_method_citations({
        "mlip_backend": "mace",
        "mlip_model": "MACE-OMOL-0",
    }))

    assert sum(line.endswith(" MACE:") for line in block.splitlines()) == 1
    assert "11423-11436" in block
    assert "56-67" in block
    assert "https://doi.org/10.1038/s42256-024-00956-x" in block
    assert "arXiv:2505.08762" in block


def test_orb_omat_model_does_not_claim_omol25_training_data() -> None:
    from mlmm.io.summary import method_references

    methods = [reference["method"] for reference in method_references({
        "mlip_backend": "orb",
        "mlip_model": "orb_v3_conservative_inf_omat",
    })]
    assert "OMol25" not in methods


def test_unknown_orb_family_does_not_claim_orb_v3() -> None:
    from mlmm.io.summary import method_references

    methods = [reference["method"] for reference in method_references({
        "mlip_backend": "orb",
        "mlip_model": "orb_v2",
    })]
    assert "Orb-v3" not in methods


def test_runtime_citations_match_the_p2r_paper_bibliography() -> None:
    from mlmm.io.summary import format_method_citations

    block = "\n".join(format_method_citations({
        "mlip_backend": "mace",
        "mlip_model": "MACE-OMOL-0",
        "pipeline_mode": "path-search",
        "mep_mode": "gsm",
        "post_segments": [{"tsopt": {}, "irc": {}, "endpoint_opt": {}}],
        "ts_opt_mode": "hess",
        "endpoint_opt_mode": "hess",
    }))

    assert "Levine, D. S.; Shuaibi, M." in block
    assert "Batatia, I.; Kovács, D. P." in block
    assert "https://doi.org/10.1021/ct400319w" in block
    assert "https://doi.org/10.1063/1.3514202" in block
    assert "https://doi.org/10.1063/1.1724823" in block


def test_dmf_and_split_ts_endpoint_references_follow_effective_settings() -> None:
    from mlmm.io.summary import format_method_citations

    base = {
        "pipeline_mode": "tsopt-only",
        "mep_mode": "dmf",
        "post_segments": [{"tsopt": {}, "irc": {}, "endpoint_opt": {}}],
        "ts_opt_mode": "hess",
        "endpoint_opt_mode": "grad",
    }
    ts_only = "\n".join(format_method_citations(base))
    correlated_path = "\n".join(
        format_method_citations(
            {
                **base,
                "pipeline_mode": "path-search",
                "path_opt_mode": "grad",
                "dmf_correlated": True,
            }
        )
    )

    assert "RS-P-RFO" in ts_only
    assert "Limited-memory BFGS (L-BFGS)" in ts_only
    assert "Euler predictor-corrector IRC" in ts_only
    assert "Correlated FB-ENM (CFB-ENM)" in correlated_path


def test_final_stdout_places_citations_immediately_before_elapsed(capsys) -> None:
    from mlmm.workflows.all import _emit_final_summary

    _emit_final_summary(
        None,
        time.time(),
        citation_payload={
            "pipeline_mode": "path-search",
            "mep_mode": "dmf",
            "path_opt_mode": "grad",
            "post_segments": [],
        },
    )

    lines = [line for line in capsys.readouterr().out.splitlines() if line]
    assert lines[-1].startswith("[time] Elapsed Time for Whole Pipeline")
    # stdout heads its blocks like every other section; the numbered form
    # belongs to summary.log alone.
    assert "====== Citations & References ======" in lines[:-1]
    assert not any(line.startswith("[6] ") for line in lines)


def test_final_stdout_explains_non_success_scientific_status(
    tmp_path, capsys
) -> None:
    from mlmm.workflows.all import _emit_final_summary

    (tmp_path / "summary.json").write_text(
        json.dumps(
            {
                "execution_status": "completed",
                "scientific_status": "partial",
                "scientific_status_reasons": ["IRC endpoint was not validated"],
            }
        ),
        encoding="utf-8",
    )

    _emit_final_summary(tmp_path, time.time())

    output = capsys.readouterr().out
    assert "Scientific status: partial" in output
    assert "RESULT WARNING: IRC endpoint was not validated." in output
    assert "Status reason:" not in output
    assert output.rstrip().splitlines()[-1].startswith(
        "[time] Elapsed Time for Whole Pipeline"
    )


def test_final_stdout_does_not_repeat_active_recovery_flags(tmp_path, capsys) -> None:
    from mlmm.workflows.all import _emit_final_summary

    (tmp_path / "summary.json").write_text(
        json.dumps(
            {
                "execution_status": "completed",
                "scientific_status": "partial",
                "scientific_status_reasons": [
                    "segment 2: TS imaginary-mode validation found n_imag=2, expected 1"
                ],
                "config": {"refine_path": True, "flatten": True},
            }
        ),
        encoding="utf-8",
    )

    _emit_final_summary(tmp_path, time.time())

    output = capsys.readouterr().out
    assert "TS imaginary-mode validation found n_imag=2." in output
    assert "Consider --" not in output


def test_citation_block_headers_match_their_destination() -> None:
    """`summary.log` numbers its sections, stdout uses `====== ... ======`. The
    two shared one renderer, so the log's section index leaked to the console.
    """
    from mlmm.io.summary import format_method_citations

    payload = {
        "pipeline_mode": "all",
        "ts_opt_mode": "rsirfo",
        "tsopt_executed": True,
    }

    log_block = format_method_citations(payload)
    assert log_block[0] == "[6] Methods and citations"
    assert "RS-I-RFO" in "\n".join(log_block)
    assert "RS-P-RFO" not in "\n".join(log_block)

    stdout_block = format_method_citations(
        payload, header="====== Citations & References ======"
    )
    assert stdout_block[0] == "====== Citations & References ======"
    # Only the header differs; the citations themselves are one source.
    assert log_block[1:] == stdout_block[1:]
