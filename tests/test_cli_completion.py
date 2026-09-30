"""Result publication and exit status agree across completion paths."""

import json

import click
import pytest
from click.testing import CliRunner

from mlmm.cli.completion import completion_guard, record_completion
from mlmm.core.utils import write_result_json


@pytest.mark.parametrize(
    ("terminal", "execution", "scientific"),
    [
        ("converged", "completed", "success"),
        ("not_converged", "completed", "failed"),
        ("stalled", "completed", "failed"),
        ("partial", "completed", "partial"),
        ("error", "failed", "failed"),
    ],
)
@pytest.mark.parametrize("out_json", [False, True])
def test_verdict_and_exit_code_do_not_depend_on_json(
    tmp_path, terminal, execution, scientific, out_json
):
    cleaned = []

    @click.command()
    def command():
        data = {"status": terminal}
        try:
            record_completion(data, command="opt")
            if out_json:
                write_result_json(tmp_path, data, command="opt")
        finally:
            cleaned.append(True)

    command.callback = completion_guard(command.callback)
    result = CliRunner().invoke(command)

    assert result.exit_code == (1 if scientific == "failed" else 0), result.output
    assert cleaned == [True]
    if out_json:
        primary = tmp_path / "result.json"
        payload = json.loads(primary.read_text())
        assert payload["execution_status"] == execution
        assert payload["scientific_status"] == scientific
        assert "status" not in payload
        assert primary.read_bytes() == (tmp_path / "summary.json").read_bytes()
        if terminal in {"converged", "not_converged", "stalled"}:
            assert payload["optimization_status"] == terminal
    else:
        assert not (tmp_path / "result.json").exists()


def test_execution_failure_with_partial_results_exits_one(tmp_path):
    @click.command()
    def command():
        write_result_json(
            tmp_path, {"execution_status": "failed", "scientific_status": "partial"},
            command="all",
        )

    command.callback = completion_guard(command.callback)
    result = CliRunner().invoke(command)
    assert result.exit_code == 1
    assert json.loads((tmp_path / "result.json").read_text())["scientific_status"] == "partial"


@pytest.mark.parametrize("code", [1, 2, 130])
def test_usage_and_interrupt_codes_remain_distinct_from_execution_errors(code):
    @click.command()
    def command():
        raise SystemExit(code)

    command.callback = completion_guard(command.callback)
    result = CliRunner().invoke(command)
    assert result.exit_code == code


def test_removed_scan_flag_is_rejected_for_each_scan_command(tmp_path):
    from mlmm.cli import cli

    source = tmp_path / "h2.xyz"
    source.write_text("2\nH2\nH 0 0 0\nH 0 0 2\n")
    parm = tmp_path / "system.parm7"
    parm.write_text("fixture")
    for command in ("scan", "scan2d", "scan3d"):
        for flag in ("--print-parsed", "--no-print-parsed"):
            result = CliRunner().invoke(cli, [command, "-i", str(source), "--parm7", str(parm), "-q", "0", "--dry-run", flag])
            assert result.exit_code == 2
            assert flag in result.output, result.output


def test_limited_smoke_checks_the_actual_exit_and_rejects_execution_failure():
    import importlib.util
    from pathlib import Path

    path = Path(__file__).with_name("smoke") / "run_limited.py"
    spec = importlib.util.spec_from_file_location("limited_smoke", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    success = {"execution_status": "completed", "scientific_status": "success"}
    partial = {"execution_status": "completed", "scientific_status": "partial"}
    unconverged = {"execution_status": "completed", "scientific_status": "failed"}
    assert module.check_completion(0, success) == "success"
    assert module.check_completion(0, partial) == "partial"
    assert module.check_completion(1, unconverged) == "failed"
    with pytest.raises(ValueError):
        module.check_completion(0, unconverged)
    with pytest.raises(ValueError):
        module.check_completion(1, {**partial, "execution_status": "failed"})
    with pytest.raises(ValueError):
        module.check_completion(2, unconverged)


def test_skill_checker_rejects_removed_status_and_wrong_axis_values(tmp_path, monkeypatch):
    import importlib.util
    from pathlib import Path

    checker_path = Path(__file__).parents[1] / ".github/scripts/check_skill_drift.py"
    spec = importlib.util.spec_from_file_location("result_field_checker", checker_path)
    checker = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(checker)
    monkeypatch.setattr(checker, "REPO_ROOT", tmp_path)
    skill = tmp_path / "skills/mlmm-cli/example.md"
    skill.parent.mkdir(parents=True)
    skill.write_text('{"execution_status": "completed", "scientific_status": "partial"}')
    assert checker._scan_file(skill, set()) == []
    for invalid in ('{"status": "completed"}', '{"execution_status": "success"}',
                    '{"scientific_status": "completed"}'):
        skill.write_text(invalid)
        assert len(checker._scan_file(skill, set())) == 1


@pytest.mark.parametrize("child_execution", ["completed", "failed"])
def test_parent_distinguishes_child_nonconvergence_from_caught_exception(tmp_path, child_execution):
    from mlmm.workflows.all import _run_cli_main

    @click.command()
    def child():
        write_result_json(
            tmp_path / "child",
            {"execution_status": child_execution, "scientific_status": "failed"},
            command="opt",
        )
        raise SystemExit(1)

    child.callback = completion_guard(child.callback)

    @click.command()
    def parent():
        _run_cli_main("opt", child, [], on_nonzero="warn", on_exception="warn")
        write_result_json(
            tmp_path / "parent",
            {"execution_status": "completed", "scientific_status": "partial"},
            command="all",
        )

    parent.callback = completion_guard(parent.callback)
    result = CliRunner().invoke(parent)
    assert result.exit_code == (1 if child_execution == "failed" else 0), result.output
    payload = json.loads((tmp_path / "parent/result.json").read_text())
    assert payload["execution_status"] == child_execution
    assert payload["scientific_status"] == "partial"


@pytest.mark.parametrize("kind", ["usage", "interrupt"])
def test_parent_propagates_child_usage_error_and_interrupt(kind):
    from mlmm.workflows.all import _run_cli_main

    @click.command()
    def child():
        if kind == "usage":
            raise click.BadParameter("invalid child configuration")
        raise SystemExit(130)

    @click.command()
    def parent():
        _run_cli_main("opt", child, [], on_nonzero="warn", on_exception="warn")

    result = CliRunner().invoke(parent)
    assert result.exit_code == (2 if kind == "usage" else 130), result.output


@pytest.mark.parametrize("segments", [[], [{"index": 1, "kind": "seg", "converged": False}]])
def test_all_unusable_results_do_not_imply_an_execution_exception(tmp_path, segments):
    from mlmm.workflows.all import _apply_pipeline_truth

    @click.command()
    def command():
        summary = {"segments": segments}
        _apply_pipeline_truth(summary, post_segments=None, config={}, legacy_status="failed")
        write_result_json(tmp_path, summary, command="all")

    command.callback = completion_guard(command.callback)
    result = CliRunner().invoke(command)
    assert result.exit_code == 1, result.output
    payload = json.loads((tmp_path / "result.json").read_text())
    assert payload["execution_status"] == "completed"
    assert payload["scientific_status"] == "failed"


def test_add_elem_info_runtime_io_error_exits_one(monkeypatch, tmp_path):
    from mlmm.domain import add_elem_info

    source = tmp_path / "input.pdb"
    source.write_text("END\n")
    def cannot_write(*args, **kwargs):
        raise OSError("disk write failed")
    monkeypatch.setattr(add_elem_info, "assign_elements", cannot_write)
    result = CliRunner().invoke(add_elem_info.cli, ["-i", str(source)])
    assert result.exit_code == 1, result.output
    assert "disk write failed" in result.output


def test_utility_input_errors_exit_two(tmp_path):
    from mlmm.cli import cli

    empty = tmp_path / "empty"
    empty.mkdir()
    result = CliRunner().invoke(cli, ["fix-altloc", "-i", str(empty)])
    assert result.exit_code == 2, result.output
    result = CliRunner().invoke(cli, ["bond-summary", "-i", str(tmp_path / "missing.xyz"), "-i", str(tmp_path / "other.xyz")])
    assert result.exit_code == 2, result.output
