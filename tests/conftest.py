"""Shared fixtures for mlmm_toolkit tests."""

import sys

import pytest

# Minimum Python version required by mlmm_toolkit
REQUIRES_PY311 = pytest.mark.skipif(
    sys.version_info < (3, 11),
    reason="mlmm_toolkit requires Python >= 3.11",
)
