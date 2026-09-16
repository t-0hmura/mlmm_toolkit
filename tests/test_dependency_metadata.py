"""Release metadata contracts for dependencies with coupled runtime behavior."""

from pathlib import Path
import tomllib

import pytest
from packaging.requirements import Requirement
from packaging.specifiers import SpecifierSet


def test_parmed_requirement_supports_the_selected_numpy_runtime() -> None:
    pyproject = (Path(__file__).parents[1] / "pyproject.toml").read_text(encoding="utf-8")

    requirements = {r.name: r for r in map(Requirement, tomllib.loads(pyproject)["project"]["dependencies"])}
    assert "2.4.6" in requirements["numpy"].specifier
    assert "2.5.0" not in requirements["numpy"].specifier
    assert '"ParmEd>=4.3.1"' in pyproject


def test_no_deps_ci_jobs_install_the_current_fairchem_torch_minor() -> None:
    root = Path(__file__).parents[1]
    pyproject = (root / "pyproject.toml").read_text(encoding="utf-8")
    requirements = {r.name: r for r in map(Requirement, tomllib.loads(pyproject)["project"]["dependencies"])}
    assert "2.13.0" in requirements["torch"].specifier
    assert "2.14.0" not in requirements["torch"].specifier
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
    assert extras["aimnet"] == ["aimnet>=0.2.0"]
    assert extras["mcp"] == ["mcp[cli]>=1.29,<2"]


@pytest.mark.parametrize("python_version", ["3.11", "3.12", "3.13"])
def test_orb_extra_selects_the_supported_api_for_each_python(python_version):
    project = tomllib.loads(
        (Path(__file__).resolve().parents[1] / "pyproject.toml").read_text(encoding="utf-8")
    )["project"]
    requirements = [Requirement(value) for value in project["optional-dependencies"]["orb"]]
    environment = {"python_version": python_version, "python_full_version": python_version + ".0"}
    eligible = [requirement for requirement in requirements
                if requirement.marker is None or requirement.marker.evaluate(environment)]
    assert len(eligible) == 1
    assert eligible[0].name == "orb-models"
    specifier = eligible[0].specifier
    if python_version == "3.11":
        assert "0.5.5" in specifier and "0.5.99" in specifier
        assert "0.5.4" not in specifier and "0.6.0" not in specifier and "0.7.0" not in specifier
    else:
        assert "0.7.0" in specifier and "0.8.0" in specifier
        assert "0.6.99" not in specifier


def test_version_tuple_matches_the_packaged_version() -> None:
    from packaging.version import Version
    from mlmm import _version

    assert _version.__version_tuple__ == _version.version_tuple == Version(_version.__version__).release
    assert _version.version == _version.__version__


@pytest.mark.parametrize(("version", "supported"), [
    ("3.9.23", False), ("3.10.18", False), ("3.11.0", True),
    ("3.11.13", True), ("3.12.0", True), ("3.12.13", True),
    ("3.13.0", True), ("3.14.0", True),
])
def test_python_requirement_keeps_core_floor_without_orb_upper_bound(version, supported):
    project = tomllib.loads(
        (Path(__file__).resolve().parents[1] / "pyproject.toml").read_text(encoding="utf-8")
    )["project"]
    assert (version in SpecifierSet(project["requires-python"])) is supported
