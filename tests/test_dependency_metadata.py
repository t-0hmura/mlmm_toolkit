"""Release metadata contracts for dependencies with coupled version bounds."""

from pathlib import Path


def test_parmed_requirement_supports_the_numpy_2_runtime() -> None:
    pyproject = (Path(__file__).parents[1] / "pyproject.toml").read_text(encoding="utf-8")

    assert '"numpy>=2.0,<2.5"' in pyproject
    assert '"ParmEd>=4.3.1"' in pyproject


def test_no_deps_ci_jobs_install_the_declared_torch_minor() -> None:
    root = Path(__file__).parents[1]
    pyproject = (root / "pyproject.toml").read_text(encoding="utf-8")
    assert '"torch~=2.8.0"' in pyproject

    workflows = ("pytest.yml", "docs_quality.yml", "markers.yml")
    install = "pip install torch==2.8.0 --index-url https://download.pytorch.org/whl/cpu"
    for name in workflows:
        text = (root / ".github" / "workflows" / name).read_text(encoding="utf-8")
        assert install in text, name
        assert "pip install torch --index-url" not in text, name
