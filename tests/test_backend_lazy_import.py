from __future__ import annotations

import subprocess
import sys
import textwrap


def test_broken_fairchem_install_stops_only_uma_runs() -> None:
    # A fairchem import that fails with a non-ImportError (warp needing a newer
    # glibc raises RuntimeError) must not stop commands that never build UMA.
    code = textwrap.dedent(
        """
        import importlib.abc
        import sys

        class BrokenFairchem(importlib.abc.MetaPathFinder):
            def find_spec(self, name, path=None, target=None):
                if name == "fairchem" or name.startswith("fairchem."):
                    raise RuntimeError("Failed to initialize warp")
                return None

        sys.meta_path.insert(0, BrokenFairchem())
        import mlmm.cli.app
        import mlmm.workflows.all
        import mlmm.workflows.dft
        from mlmm.backends import mlmm_calc

        try:
            mlmm_calc._load_fairchem()
        except RuntimeError as exc:
            print("UMA:", exc)
        """
    )
    proc = subprocess.run([sys.executable, "-c", code], capture_output=True, text=True)
    assert proc.returncode == 0, proc.stderr
    assert "UMA: Failed to initialize warp" in proc.stdout
