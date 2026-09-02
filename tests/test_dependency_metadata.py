"""Release metadata contracts for dependencies with coupled runtime behavior."""

from pathlib import Path
import tomllib


def test_parmed_requirement_supports_the_unbounded_numpy_runtime() -> None:
    pyproject = (Path(__file__).parents[1] / "pyproject.toml").read_text(encoding="utf-8")

    assert '"numpy"' in pyproject
    assert '"numpy>=' not in pyproject
    assert '"ParmEd>=4.3.1"' in pyproject


def test_no_deps_ci_jobs_install_the_current_fairchem_torch_minor() -> None:
    root = Path(__file__).parents[1]
    pyproject = (root / "pyproject.toml").read_text(encoding="utf-8")
    assert '"torch"' in pyproject
    assert '"torch~=' not in pyproject

    workflows = ("pytest.yml", "docs_quality.yml", "markers.yml")
    install = "pip install torch==2.13.0 --index-url https://download.pytorch.org/whl/cpu"
    for name in workflows:
        text = (root / ".github" / "workflows" / name).read_text(encoding="utf-8")
        assert install in text, name


def test_runtime_dependency_floors_match_consumed_apis() -> None:
    root = Path(__file__).parents[1]
    project = tomllib.loads(
        (root / "pyproject.toml").read_text(encoding="utf-8")
    )["project"]
    dependencies = set(project["dependencies"])
    extras = project["optional-dependencies"]

    assert "ParmEd>=4.3.1" in dependencies
    assert "plotly>=6.1.1" in dependencies
    assert extras["orb"] == ["orb-models>=0.7.0"]
    assert extras["aimnet"] == ["aimnet>=0.2.0"]
    assert extras["mcp"] == ["mcp[cli]>=1.29,<2"]
