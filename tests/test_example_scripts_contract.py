"""Contract for the advertised example shell scripts (M65).

Positive: the three README-advertised scripts exist, pass ``bash -n``, and their
47 invocations validate against the live CLI. Negative: a syntactically broken
script fails ``bash -n``, and a script carrying an invented option fails live
option validation.
"""

from __future__ import annotations

import subprocess
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
SCRIPTS = REPO_ROOT / ".github" / "scripts"
for _p in (str(REPO_ROOT), str(SCRIPTS)):
    if _p not in sys.path:
        sys.path.insert(0, _p)

import check_example_scripts as ces  # noqa: E402
import docs_command_contract as dc  # noqa: E402


def test_example_scripts_checker_passes_on_real_repo(monkeypatch) -> None:
    monkeypatch.setattr(sys, "argv", ["check_example_scripts.py"])
    assert ces.main() == 0


def test_advertised_scripts_exist_and_yield_expected_invocations() -> None:
    scripts = dc.public_shell_examples()
    assert [s.exists() for s in scripts] == [True, True, True]
    commands = dc.extract_shell_commands(scripts)
    # toy_system 32 + run_all 1 + run_stepwise 14 = 47.
    assert len(commands) == 47


def test_broken_shell_syntax_fails_bash_n(tmp_path: Path) -> None:
    bad = tmp_path / "broken.sh"
    bad.write_text("if [ -z ]; then\n  echo missing fi\n", encoding="utf-8")
    completed = subprocess.run(["bash", "-n", str(bad)], capture_output=True, text=True)
    assert completed.returncode != 0


def test_invented_option_in_example_fails_validation(tmp_path: Path) -> None:
    script = tmp_path / "run.sh"
    script.write_text(
        "mlmm all -i R.pdb P.pdb -c LIG --invented-option\n", encoding="utf-8"
    )
    commands = dc.extract_shell_commands([script])
    assert len(commands) == 1
    errors = dc.validate_option_names(commands, dc.load_root_cli())
    assert len(errors) == 1
    assert "--invented-option" in errors[0]
