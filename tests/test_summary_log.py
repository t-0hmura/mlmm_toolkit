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
    assert "RESULT WARNING" in text
    assert "Status reason       : segment 1 did not converge" in text
    assert text.index("RESULT WARNING") < text.index("ΔE‡")


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
    assert "RS-I-RFO" in block
    assert "quasi-RRHO thermochemistry" in block
    assert "Direct Max Flux (DMF)" not in block
    assert all(set(ref) == {"method", "citation", "doi"} for ref in references)
    assert len({ref["doi"] for ref in references}) == len(references)


def test_method_citations_use_actual_path_and_post_stages() -> None:
    from mlmm.io.summary import format_method_citations

    path_only = {
        "pipeline_mode": "path-search",
        "mep_mode": "dmf",
        "path_opt_mode": "grad",
        "post_opt_mode": "hess",
        "post_segments": [],
    }
    mixed = {
        **path_only,
        "post_segments": [
            {"endpoint_opt": {}, "ts_imag": {"n_imag": 1}}
        ],
    }

    path_text = "\n".join(format_method_citations(path_only))
    mixed_text = "\n".join(format_method_citations(mixed))

    assert "Limited-memory BFGS (L-BFGS)" in path_text
    assert "RFO / P-RFO" not in path_text
    assert "RS-I-RFO" not in path_text
    assert "quasi-RRHO thermochemistry" not in path_text
    assert "Limited-memory BFGS (L-BFGS)" in mixed_text
    assert "RFO / P-RFO" in mixed_text
    assert "RS-I-RFO" in mixed_text
    assert "quasi-RRHO thermochemistry" not in mixed_text


def test_dmf_and_split_ts_endpoint_references_follow_effective_settings() -> None:
    from mlmm.io.summary import format_method_citations

    base = {
        "pipeline_mode": "tsopt-only",
        "mep_mode": "dmf",
        "post_segments": [{"endpoint_opt": {}}],
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

    assert "RS-I-RFO" in ts_only
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
    assert "RESULT WARNING:" in output
    assert "Status reason: IRC endpoint was not validated" in output
    assert output.rstrip().splitlines()[-1].startswith(
        "[time] Elapsed Time for Whole Pipeline"
    )


def test_citation_block_headers_match_their_destination() -> None:
    """`summary.log` numbers its sections, stdout uses `====== ... ======`. The
    two shared one renderer, so the log's section index leaked to the console.
    """
    from mlmm.io.summary import format_method_citations

    payload = {"pipeline_mode": "all", "ts_opt_mode": "rsirfo"}

    log_block = format_method_citations(payload)
    assert log_block[0] == "[6] Methods and citations"

    stdout_block = format_method_citations(
        payload, header="====== Citations & References ======"
    )
    assert stdout_block[0] == "====== Citations & References ======"
    # Only the header differs; the citations themselves are one source.
    assert log_block[1:] == stdout_block[1:]
