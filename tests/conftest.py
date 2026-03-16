"""Shared fixtures for mlmm tests."""

import sys

import pytest

# Minimum Python version required by mlmm
REQUIRES_PY311 = pytest.mark.skipif(
    sys.version_info < (3, 11),
    reason="mlmm requires Python >= 3.11",
)
