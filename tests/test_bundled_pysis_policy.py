"""Policies for the bundled pysisyphus integration."""

from __future__ import annotations

import inspect
import os
from pathlib import Path
import subprocess
import sys
from types import SimpleNamespace

import numpy as np
import pytest

from pysisyphus.irc.IRC import IRC
from pysisyphus.tr_projection import (
    TR_PROJECTION_MODES,
    normalize_tr_projection_mode,
)
from mlmm.workflows.irc import IRC_KW_DEFAULT
from mlmm.workflows.tsopt import (
    _finalize_dimer_saddle_status,
    _hessian_postprocessing_is_ready,
    _hessian_result_status,
    _heavy_ts_terminal_status,
)


def test_missing_optional_config_is_silent(tmp_path: Path) -> None:
    repo = Path(__file__).parents[1]
    env = os.environ.copy()
    env["HOME"] = str(tmp_path)
    env.pop("PYSISRC", None)
    env["PYTHONPATH"] = os.pathsep.join(
        part for part in (str(repo), env.get("PYTHONPATH", "")) if part
    )

    proc = subprocess.run(
        [sys.executable, "-c", "import pysisyphus.config"],
        cwd=repo,
        env=env,
        capture_output=True,
        text=True,
        check=False,
    )

    assert proc.returncode == 0, proc.stderr
    assert "Couldn't find configuration file" not in proc.stdout
    assert "Couldn't find configuration file" not in proc.stderr


def test_hessian_postprocessing_requires_numerical_convergence() -> None:
    assert not _hessian_postprocessing_is_ready(None)
    assert not _hessian_postprocessing_is_ready(SimpleNamespace())
    assert _hessian_postprocessing_is_ready(SimpleNamespace(is_converged=True))
    assert not _hessian_postprocessing_is_ready(SimpleNamespace(is_stalled=True))
    assert not _hessian_postprocessing_is_ready(
        SimpleNamespace(_last_exact_failure_reason="RuntimeError: failed")
    )


def test_hessian_result_status_distinguishes_skip_failure_and_completion() -> None:
    assert _hessian_result_status(
        n_imaginary=None,
        hessian_error=None,
        postprocessing_ready=False,
    ) == "skipped"
    assert _hessian_result_status(
        n_imaginary=None,
        hessian_error=None,
        postprocessing_ready=True,
        explicitly_skipped=True,
    ) == "skipped"
    assert _hessian_result_status(
        n_imaginary=None,
        hessian_error="RuntimeError: failed",
        postprocessing_ready=False,
    ) == "failed"
    assert _hessian_result_status(
        n_imaginary=1,
        hessian_error="RuntimeError: failed",
        postprocessing_ready=True,
    ) == "failed"
    assert _hessian_result_status(
        n_imaginary=1,
        hessian_error=None,
        postprocessing_ready=True,
    ) == "completed"


def test_missing_explicit_config_does_not_report_false_success(
    tmp_path: Path,
) -> None:
    repo = Path(__file__).parents[1]
    env = os.environ.copy()
    env["PYSISRC"] = str(tmp_path / "missing.pysisyphusrc")
    env["PYTHONPATH"] = os.pathsep.join(
        part for part in (str(repo), env.get("PYTHONPATH", "")) if part
    )

    proc = subprocess.run(
        [sys.executable, "-c", "import pysisyphus.config"],
        cwd=repo,
        env=env,
        capture_output=True,
        text=True,
        check=False,
    )

    assert proc.returncode == 0, proc.stderr
    assert "Read pysisyphus configuration" not in proc.stdout
    assert "Read pysisyphus configuration" not in proc.stderr


def test_irc_hdf5_is_opt_in_and_never_contains_a_hessian() -> None:
    assert inspect.signature(IRC).parameters["dump_every"].default is None
    assert "dump_every" in IRC_KW_DEFAULT
    assert IRC_KW_DEFAULT["dump_every"] is None
    assert IRC_KW_DEFAULT["dump_fn"] == "irc_data.h5"

    irc = IRC.__new__(IRC)
    irc.all_energies = [0.0]
    irc.all_coords = [np.zeros(3)]
    irc.all_gradients = [np.zeros(3)]
    irc.all_mw_coords = [np.zeros(3)]
    irc.all_mw_gradients = [np.zeros(3)]
    irc.ts_index = 0

    assert all(
        "hess" not in key.lower() for key in irc.get_full_irc_data()
    )


def test_the_removed_legacy_projection_is_rejected_not_silently_accepted() -> None:
    assert TR_PROJECTION_MODES == ("constrained",)
    assert normalize_tr_projection_mode(None) == "constrained"
    with pytest.raises(ValueError, match="Unknown TR projection mode"):
        normalize_tr_projection_mode("legacy-active")

    runner = SimpleNamespace(
        tr_projection="constrained",
        freeze_atoms=[0],
        is_converged=True,
        stop_reason="",
    )
    indices = _finalize_dimer_saddle_status(
        runner, np.array([-100.0, 25.0]), 5.0
    )

    assert indices.tolist() == [0]
    assert runner.n_imaginary_modes == 1
    assert runner.imaginary_frequencies_cm == [-100.0]
    assert runner.saddle_order_verified is True
    assert runner.is_converged is True

    assert _heavy_ts_terminal_status(
        optimizer_converged=True,
        n_imag=1,
        stalled=False,
    ) == "converged"
    assert _heavy_ts_terminal_status(
        optimizer_converged=True, n_imag=2, stalled=False
    ) == "converged"
