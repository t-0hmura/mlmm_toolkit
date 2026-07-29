"""Keep the scheduled smoke script synchronized with the live Click CLI."""

from __future__ import annotations

import importlib.util
import subprocess
import sys
from pathlib import Path


def _load_docs_contract():
    path = (
        Path(__file__).parents[1]
        / ".github"
        / "scripts"
        / "docs_command_contract.py"
    )
    spec = importlib.util.spec_from_file_location("mlmm_docs_command_contract", path)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


contract = _load_docs_contract()
SMOKE_SCRIPT = Path(__file__).parents[1] / "tests" / "smoke" / "run.sh"


def _literal_smoke_commands(path: Path):
    return [
        command
        for command in contract.extract_shell_commands([path])
        if command.text.startswith("mlmm ")
    ]


def test_smoke_script_is_valid_shell_and_uses_live_click_options() -> None:
    syntax = subprocess.run(
        ["bash", "-n", str(SMOKE_SCRIPT)],
        text=True,
        capture_output=True,
        check=False,
    )
    assert syntax.returncode == 0, syntax.stderr

    commands = _literal_smoke_commands(SMOKE_SCRIPT)
    assert len(commands) >= 80
    assert contract.validate_option_names(commands, contract.load_root_cli()) == []


def test_smoke_contract_rejects_an_invented_option(tmp_path: Path) -> None:
    script = tmp_path / "run.sh"
    script.write_text(
        "mlmm path-search -i r.pdb p.pdb --invented-option out\n",
        encoding="utf-8",
    )
    errors = contract.validate_option_names(
        _literal_smoke_commands(script),
        contract.load_root_cli(),
    )
    assert len(errors) == 1
    assert "unknown option '--invented-option' for path-search" in errors[0]


def test_manual_all_smoke_disables_layer_detection() -> None:
    command = next(
        command.text
        for command in _literal_smoke_commands(SMOKE_SCRIPT)
        if "--parm \"${test19_parms[0]}\"" in command.text
    )
    assert "--model-pdb test19/ml_region.pdb" in command
    assert "--no-detect-layer" in command


def test_required_positive_lane_uses_release_settings_and_runs_last() -> None:
    command = next(
        command.text
        for command in _literal_smoke_commands(SMOKE_SCRIPT)
        if "--out-dir test19" in command.text
    )
    assert "--deterministic" in command
    assert "--no-refine-path" in command
    assert "--thresh gau" in command
    assert "--thresh-post baker" in command
    assert "--tsopt" in command
    assert "--thermo" in command
    assert "--dft" in command
    assert "--flatten" in command
    assert "--irc-never-stop" in command
    assert "--max-cycles" not in command
    assert "--tsopt-max-cycles" not in command

    script = SMOKE_SCRIPT.read_text(encoding="utf-8")
    assert script.index("--out-dir test74_flatten") < script.index("--out-dir test19")
    assert script.index("backend_hessian.out") < script.index("--out-dir test19")
    assert script.index("--out-dir test19") < script.index("--out-dir test20")
