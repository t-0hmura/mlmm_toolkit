"""all keeps the DFT child's output in dft/run.log once and shows its stderr."""

from __future__ import annotations

import subprocess
from types import SimpleNamespace

import pytest

from mlmm.workflows import all as workflow


@pytest.mark.parametrize("returncode", [0, 1])
def test_parent_dft_keeps_the_child_run_log_once(tmp_path, monkeypatch, returncode):
    child_stdout = "[mode] dft\nconverged SCF energy = -1.11\n"
    child_stderr = "warning: soft limit\n"
    error = "Error: PySCF import failed.\n" if returncode else ""
    echoes = []

    def run(cmd, **kwargs):
        # The real child tees its output to its own run.log; click prints an
        # error that ends the command after that log closes.
        out_dir = cmd[cmd.index("--out-dir") + 1]
        with open(f"{out_dir}/run.log", "a", encoding="utf-8") as handle:
            handle.write("[command] mlmm dft\n" + child_stdout + child_stderr)
        return SimpleNamespace(returncode=returncode, stdout=child_stdout,
                               stderr=child_stderr + error)

    monkeypatch.setattr(subprocess, "run", run)
    monkeypatch.setattr(workflow, "_echo", lambda text, **kw: echoes.append((text, kw)))
    result = workflow._run_dft_for_state(
        tmp_path / "ts.pdb", 0, 1, tmp_path / "real.parm7", tmp_path / "model.pdb",
        False, tmp_path / "dft", None, func_basis="hf/sto-3g",
    )

    log = (tmp_path / "dft/run.log").read_text()
    assert log.count("converged SCF energy") == 1
    assert log.count("warning: soft limit") == 1
    assert log.count("Error: PySCF import failed.") == (1 if returncode else 0)
    stderr_echoes = [kw for text, kw in echoes if "warning: soft limit" in text]
    assert stderr_echoes == [{"err": True, "narrative": True}]
    assert result["_dft_failed"] is True
