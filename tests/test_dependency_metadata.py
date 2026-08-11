"""Release metadata contracts for dependencies with coupled version bounds."""

from pathlib import Path


def test_parmed_requirement_supports_the_numpy_2_runtime() -> None:
    pyproject = (Path(__file__).parents[1] / "pyproject.toml").read_text(encoding="utf-8")

    assert '"numpy>=2.0,<2.5"' in pyproject
    assert '"ParmEd>=4.3.1"' in pyproject
