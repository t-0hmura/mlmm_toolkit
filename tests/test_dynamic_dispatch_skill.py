"""Contract checks for the copyable PBS dynamic-dispatch recipe."""

from pathlib import Path
import re
import subprocess


SKILL = Path(__file__).parents[1] / "skills" / "mlmm-hpc" / "dynamic-dispatch.md"


def _dispatcher_script() -> str:
    text = SKILL.read_text(encoding="utf-8")
    section = text.split("## dispatcher.sh", 1)[1]
    return section.split("```bash", 1)[1].split("```", 1)[0]


def test_dispatcher_is_copyable_and_preserves_failure_contract() -> None:
    script = _dispatcher_script()
    assert '${PBS_JOBID:?PBS_JOBID is required' in script
    assert '_state.${JOB_ID}.txt' in script
    assert '_lock.${JOB_ID}' in script
    assert '_worker.${JOB_ID}.sh' in script
    assert "awk 'END {print NR}'" in script
    assert "failed=1" in script
    assert 'exit "${failed}"' in script

    copyable = re.sub(r"<[A-Z][A-Z0-9_:]*>", "placeholder", script)
    result = subprocess.run(
        ["bash", "-n"],
        input=copyable,
        text=True,
        capture_output=True,
        check=False,
    )
    assert result.returncode == 0, result.stderr
