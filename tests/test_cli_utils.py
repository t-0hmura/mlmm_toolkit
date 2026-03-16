"""Unit tests for cli_utils.py."""

import sys

import pytest

# cli_utils uses ``from __future__ import annotations`` (Python 3.7+)
# and other modern features.  Skip the entire module on older interpreters.
pytestmark = pytest.mark.skipif(
    sys.version_info < (3, 11),
    reason="mlmm.cli_utils requires Python >= 3.11",
)


@pytest.fixture
def _parse_bool():
    from mlmm.cli_utils import parse_bool
    return parse_bool


@pytest.fixture
def _argparse_bool():
    from mlmm.cli_utils import argparse_bool
    return argparse_bool


class TestParseBool:
    @pytest.mark.parametrize("value", ["true", "True", "TRUE", "1", "yes", "y", "t"])
    def test_true_values(self, value, _parse_bool):
        assert _parse_bool(value) is True

    @pytest.mark.parametrize("value", ["false", "False", "FALSE", "0", "no", "n", "f"])
    def test_false_values(self, value, _parse_bool):
        assert _parse_bool(value) is False

    def test_none_raises(self, _parse_bool):
        with pytest.raises(ValueError, match="None"):
            _parse_bool(None)

    def test_invalid_raises(self, _parse_bool):
        with pytest.raises(ValueError, match="Invalid boolean"):
            _parse_bool("maybe")


class TestArgparseBool:
    def test_valid(self, _argparse_bool):
        assert _argparse_bool("true") is True
        assert _argparse_bool("false") is False

    def test_invalid_raises_argparse_error(self, _argparse_bool):
        import argparse
        with pytest.raises(argparse.ArgumentTypeError):
            _argparse_bool("invalid")
