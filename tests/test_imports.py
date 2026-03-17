"""Import smoke tests for mlmm."""

import sys

import pytest

# The cli and cli_utils modules use Python 3.10+ syntax (str | None)
# and 3.7+ features (from __future__ import annotations).
# Skip those tests on older interpreters.
_need_py311 = pytest.mark.skipif(
    sys.version_info < (3, 11),
    reason="mlmm.cli / cli_utils require Python >= 3.11",
)


def test_import_mlmm():
    import mlmm


def test_import_defaults():
    from mlmm.defaults import MLMM_CALC_KW, GEOM_KW_DEFAULT, OPT_BASE_KW


@_need_py311
def test_import_cli():
    from mlmm.cli import cli
    assert callable(cli)


@_need_py311
def test_import_cli_utils():
    from mlmm.cli_utils import parse_bool, run_cli
    assert callable(parse_bool)
    assert callable(run_cli)
