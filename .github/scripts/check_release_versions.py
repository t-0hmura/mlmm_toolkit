#!/usr/bin/env python3
"""Fail when release metadata, landing headers, license identity, the release/docs
workflows or the bundled Colab ref disagree."""

from __future__ import annotations

import argparse
import ast
import json
import re
import tomllib
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[2]
NOTEBOOK = REPO_ROOT / "examples" / "mlmm_colab.ipynb"
NOTEBOOK_REF = "mlmm_toolkit_version"
CHANGELOG = REPO_ROOT / "CHANGELOG.md"
WORKFLOWS = REPO_ROOT / ".github" / "workflows"
RELEASE_WORKFLOW = WORKFLOWS / "release.yml"
DOCS_WORKFLOW = WORKFLOWS / "docs.yml"
MANIFEST = REPO_ROOT / "MANIFEST.in"
PYPROJECT = REPO_ROOT / "pyproject.toml"
# The CPU-only CI jobs install torch without dependencies, so they must follow
# the torch minor that pyproject allows.
NO_DEPS_CI_WORKFLOWS = ("pytest.yml", "docs_quality.yml", "markers.yml")
CI_TORCH_INSTALL = (
    "pip install torch==2.13.0 --index-url https://download.pytorch.org/whl/cpu"
)
# The primary-source title that every citation surface repeats.
CITATION_SURFACES = (
    "CITATION.cff",
    "README.md",
    "docs/index.md",
    "docs/ja/index.md",
    "mlmm/io/summary.py",
    "examples/mlmm_colab.ipynb",
)
LANDING_PAGES = (REPO_ROOT / "docs" / "index.md", REPO_ROOT / "docs" / "ja" / "index.md")
LANDING_SUBSTITUTION = "{{ release }}"
_LANDING_HEADER_MARKERS = ("Version:", "バージョン:")
_VERSION_LITERAL_RE = re.compile(r"v?\d+\.\d+(?:\.\d+)?")
EXPECTED_LICENSE = "GPL-3.0-or-later"


def _cff_version() -> str:
    text = (REPO_ROOT / "CITATION.cff").read_text(encoding="utf-8")
    match = re.search(r"^version:\s*[\"']?([^\"'#\s]+)", text, re.MULTILINE)
    if match is None:
        raise ValueError("CITATION.cff has no top-level version field")
    return match.group(1)


def _docs_release() -> str:
    tree = ast.parse((REPO_ROOT / "docs" / "conf.py").read_text(encoding="utf-8"))
    for node in tree.body:
        if isinstance(node, ast.Assign) and any(
            isinstance(target, ast.Name) and target.id == "release"
            for target in node.targets
        ):
            value = ast.literal_eval(node.value)
            return str(value)
    raise ValueError("docs/conf.py has no literal release assignment")


def _notebook_version() -> str:
    notebook = json.loads(NOTEBOOK.read_text(encoding="utf-8"))
    pattern = re.compile(rf"^{NOTEBOOK_REF}\s*=\s*[\"']v?([^\"']+)[\"']", re.MULTILINE)
    for cell in notebook.get("cells", []):
        if cell.get("cell_type") != "code":
            continue
        source = cell.get("source", "")
        if isinstance(source, list):
            # nbformat also allows a list of lines; Colab exports that form.
            source = "".join(source)
        match = pattern.search(source)
        if match is not None:
            return match.group(1)
    raise ValueError(f"{NOTEBOOK.name} has no {NOTEBOOK_REF} assignment")


def _module_version() -> str:
    tree = ast.parse(
        (REPO_ROOT / "mlmm" / "_version.py").read_text(encoding="utf-8")
    )
    for node in tree.body:
        if not isinstance(node, ast.Assign):
            continue
        if any(isinstance(target, ast.Name) and target.id == "__version__" for target in node.targets):
            return str(ast.literal_eval(node.value))
    raise ValueError("mlmm/_version.py has no literal __version__")


def _cff_license() -> str:
    text = (REPO_ROOT / "CITATION.cff").read_text(encoding="utf-8")
    match = re.search(r"^license:\s*[\"']?([^\"'#\s]+)", text, re.MULTILINE)
    if match is None:
        raise ValueError("CITATION.cff has no top-level license field")
    return match.group(1)


def _pyproject_license() -> str:
    data = tomllib.loads((REPO_ROOT / "pyproject.toml").read_text(encoding="utf-8"))
    license_value = data.get("project", {}).get("license")
    if isinstance(license_value, dict):  # legacy table form
        license_value = license_value.get("text")
    if not isinstance(license_value, str) or not license_value:
        raise ValueError("pyproject.toml [project] has no SPDX license expression")
    return license_value


def _check_landing_pages() -> list[str]:
    errors: list[str] = []
    for path in LANDING_PAGES:
        try:
            rel = path.relative_to(REPO_ROOT)
        except ValueError:
            rel = path
        text = path.read_text(encoding="utf-8")
        header = next(
            (
                line
                for line in text.splitlines()
                if any(marker in line for marker in _LANDING_HEADER_MARKERS)
            ),
            None,
        )
        if header is None:
            errors.append(f"{rel}: no landing version header line found")
            continue
        if LANDING_SUBSTITUTION not in header:
            errors.append(
                f"{rel}: landing header must render the version via "
                f"'{LANDING_SUBSTITUTION}': {header.strip()!r}"
            )
        # Link targets such as the paper DOI are not version text.
        if _VERSION_LITERAL_RE.search(re.sub(r"\]\([^)]*\)", "]", header)):
            errors.append(
                f"{rel}: landing header must not hardcode a version literal: "
                f"{header.strip()!r}"
            )
    return errors


def _check_license() -> list[str]:
    errors: list[str] = []
    cff = _cff_license()
    pyproject = _pyproject_license()
    if pyproject != EXPECTED_LICENSE:
        errors.append(
            f"pyproject.toml license={pyproject!r}; expected {EXPECTED_LICENSE!r}"
        )
    if cff != pyproject:
        errors.append(
            f"CITATION.cff license={cff!r} != pyproject.toml license={pyproject!r}"
        )
    return errors


def _rel(path: Path) -> str:
    try:
        return str(path.relative_to(REPO_ROOT))
    except ValueError:
        return str(path)


def _text_errors(path: Path, *, present=(), absent=()) -> list[str]:
    """Required and forbidden fragments of one text file."""
    text = path.read_text(encoding="utf-8")
    errors = [f"{_rel(path)}: missing {fragment!r}" for fragment in present if fragment not in text]
    errors += [
        f"{_rel(path)}: must not contain {fragment!r}" for fragment in absent if fragment in text
    ]
    return errors


def _check_release_workflow() -> list[str]:
    """The release workflow builds a published tag, checks its title, attaches both
    distributions to the GitHub release and publishes to PyPI last."""
    publish, attach = "Publish to PyPI (Trusted Publishing)", "Attach distribution to GitHub Release"
    errors = _text_errors(
        RELEASE_WORKFLOW,
        present=(
            "release:\n    types:\n      - published",
            "workflow_dispatch:",
            "Existing published release tag to build and publish",
            "contents: write",
            "id-token: write",
            "dist/*.whl dist/*.tar.gz --clobber",
            "ref: ${{ env.RELEASE_TAG }}",
            'git merge-base --is-ancestor "$tag_commit" "$main_commit"',
            'python -m pip install --force-reinstall "${wheels[0]}"',
            "RELEASE_TAG: ${{ github.event_name == 'release' "
            "&& github.event.release.tag_name || inputs.tag }}",
            "RELEASE_NAME: ${{ github.event.release.name }}",
            '[[ "$RELEASE_NAME" != "$RELEASE_TAG" ]]',
            'gh release view "$RELEASE_TAG"',
            "--json isDraft --jq .isDraft",
            "--json isPrerelease --jq .isPrerelease",
            attach,
            publish,
        ),
        absent=("--no-deps", "attestations: false", "skip-existing: true"),
    )
    text = RELEASE_WORKFLOW.read_text(encoding="utf-8")
    name = _rel(RELEASE_WORKFLOW)
    for fragment in ("pypa/gh-action-pypi-publish@release/v1", 'gh release upload "$RELEASE_TAG"'):
        if text.count(fragment) != 1:
            errors.append(f"{name}: {fragment!r} must appear exactly once")
    if attach in text and publish in text and text.index(attach) > text.index(publish):
        errors.append(f"{name}: the release attachment must come before the PyPI publish")
    last_step = text.rsplit("      - name: ", 1)[-1]
    if publish in text and not last_step.startswith(publish):
        errors.append(f"{name}: the PyPI publish must be the last step")
    title_check = "Verify release title matches tag"
    if title_check in text and "- name: Checkout" in text:
        if text.index(title_check) > text.index("- name: Checkout"):
            errors.append(f"{name}: the title check must run before checkout")
    return errors


def _check_docs_deploy_workflow() -> list[str]:
    """Pages deploys only after Docs Quality succeeded on the same main commit;
    a manual dispatch stays an explicit override."""
    return _text_errors(
        DOCS_WORKFLOW,
        present=(
            'workflows: ["Docs Quality"]',
            "types: [completed]",
            "github.event.workflow_run.conclusion == 'success'",
            "github.event.workflow_run.event == 'push'",
            "github.event.workflow_run.head_branch == 'main'",
            "github.event.workflow_run.head_sha == github.sha",
            "github.event.workflow_run.head_sha || github.sha",
            "needs: build",
            "workflow_dispatch:",
            "github.event_name == 'workflow_dispatch'",
        ),
        absent=("\n  push:\n",),
    )


def _check_sdist_manifest() -> list[str]:
    """The sdist ships tests but prunes the local smoke products."""
    return _text_errors(
        MANIFEST,
        present=("recursive-include tests *", "\nprune tests/smoke\n"),
        absent=("recursive-exclude tests/smoke",),
    )


def _check_ci_torch_pin() -> list[str]:
    """The CPU-only CI jobs install the torch minor that the fairchem stack needs."""
    errors: list[str] = []
    for name in NO_DEPS_CI_WORKFLOWS:
        errors.extend(_text_errors(WORKFLOWS / name, present=(CI_TORCH_INSTALL,)))
    return errors


def _check_pyproject_metadata() -> list[str]:
    """Dependency floors, the ORB and DFT extras, the Python range and the torch,
    numpy and fairchem ranges that the runtime APIs require."""
    from packaging.markers import default_environment
    from packaging.requirements import Requirement
    from packaging.specifiers import SpecifierSet

    name = _rel(PYPROJECT)
    text = PYPROJECT.read_text(encoding="utf-8")
    project = tomllib.loads(text)["project"]
    dependencies = project["dependencies"]
    extras = project["optional-dependencies"]
    errors: list[str] = []

    for pin in ("ParmEd>=4.3.1", "pydmf>=1.2", "plotly>=6.1.1"):
        if pin not in dependencies:
            errors.append(f"{name}: dependencies must contain {pin!r}")
    if extras.get("aimnet") != ["aimnet>=0.2.0"]:
        errors.append(f"{name}: the aimnet extra must be ['aimnet>=0.2.0']")

    orb = [Requirement(value) for value in extras.get("orb", [])]
    if len(orb) != 1 or orb[0].marker is not None or orb[0].name != "orb-models":
        errors.append(f"{name}: the orb extra must be one unconditional orb-models requirement")
    else:
        specifier = orb[0].specifier
        if "0.7.0" not in specifier or "0.8.0" not in specifier:
            errors.append(f"{name}: the orb extra must allow orb-models 0.7.0 and later versions")
        if "0.5.5" in specifier or "0.6.99" in specifier:
            errors.append(f"{name}: the orb extra must exclude orb-models before 0.7.0")

    python = SpecifierSet(project["requires-python"])
    for version, supported in (
        ("3.10.18", False), ("3.11.13", False), ("3.12.0", True),
        ("3.12.13", True), ("3.13.0", True), ("3.14.0", False),
    ):
        if (version in python) is not supported:
            errors.append(f"{name}: requires-python must {'allow' if supported else 'exclude'} {version}")

    dft = {requirement.name: requirement for requirement in map(Requirement, extras.get("dft", []))}
    for package in ("gpu4pyscf-cuda13x", "cupy-cuda13x"):
        marker = dft[package].marker if package in dft else None
        if marker is None:
            errors.append(f"{name}: the dft extra needs a platform marker on {package}")
            continue
        for platform, expected in (("linux", True), ("darwin", False), ("win32", False)):
            environment = default_environment()
            environment.update({"sys_platform": platform, "platform_machine": "x86_64"})
            if marker.evaluate(environment) is not expected:
                errors.append(f"{name}: {package} must be {'installed' if expected else 'skipped'} on {platform}")
    cuda12 = {requirement.name for requirement in map(Requirement, extras.get("dft-cuda12", []))}
    if not {"gpu4pyscf-cuda12x", "cupy-cuda12x"} <= cuda12:
        errors.append(f"{name}: the dft-cuda12 extra must list gpu4pyscf-cuda12x and cupy-cuda12x")

    if '"torch~=' in text:
        errors.append(f"{name}: torch must not use a compatible-release pin")
    runtime = {requirement.name: requirement for requirement in map(Requirement, dependencies)}
    for package, allowed, excluded in (
        ("fairchem-core", "2.22.0", "3.0.0"),
        ("numpy", "2.4.6", "2.5.0"),
        ("torch", "2.13.0", "2.14.0"),
    ):
        if package not in runtime:
            errors.append(f"{name}: dependencies must contain {package}")
        elif allowed not in runtime[package].specifier or excluded in runtime[package].specifier:
            errors.append(f"{name}: {package} must allow {allowed} and exclude {excluded}")
    return errors


def _check_citation_title() -> list[str]:
    """Every citation surface carries the published title, not the earlier wording."""
    errors: list[str] = []
    for relative in CITATION_SURFACES:
        text = (REPO_ROOT / relative).read_text(encoding="utf-8")
        if "Towards Accelerated" not in text:
            errors.append(f"{relative}: missing the primary-source title 'Towards Accelerated'")
        if "Toward Accelerated" in text:
            errors.append(f"{relative}: must not use the earlier title 'Toward Accelerated'")
    return errors


def _cell_text(cell: dict) -> str:
    source = cell.get("source", "")
    return "".join(source) if isinstance(source, list) else str(source)


def _check_notebook_metadata() -> list[str]:
    """The Colab notebook is nbformat 4 on a GPU runtime image, keeps its three
    cells, links the repository and the paper, and every code cell parses."""
    notebook = json.loads(NOTEBOOK.read_text(encoding="utf-8"))
    name = NOTEBOOK.name
    metadata = notebook.get("metadata", {})
    cells = notebook.get("cells", [])
    errors: list[str] = []
    if notebook.get("nbformat") != 4:
        errors.append(f"{name}: nbformat must be 4")
    if metadata.get("accelerator") != "GPU":
        errors.append(f"{name}: metadata.accelerator must be 'GPU'")
    runtime_version = metadata.get("colab", {}).get("runtime_attributes", {}).get("runtime_version")
    if runtime_version != "2026.07":
        errors.append(f"{name}: the Colab runtime image must be pinned to 2026.07 (Python 3.12)")
    if len(cells) != 3:
        errors.append(f"{name}: expected 3 cells (introduction, Setup, app), found {len(cells)}")
    introduction = _cell_text(cells[0]) if cells else ""
    for fragment in (
        "[t-0hmura/mlmm_toolkit](https://github.com/t-0hmura/mlmm_toolkit)",
        "ML/MM ONIOM reaction paths on full solvated proteins",
        "**① Input → ② Setup → ③ Options → (Run) → ④ Results**",
        "[DOI: 10.26434/chemrxiv-2025-jft1k](https://doi.org/10.26434/chemrxiv-2025-jft1k)",
    ):
        if fragment not in introduction:
            errors.append(f"{name}: the introduction cell is missing {fragment!r}")
    for index, cell in enumerate(cells):
        if cell.get("cell_type") == "code":
            try:
                ast.parse(_cell_text(cell))
            except SyntaxError as exc:
                errors.append(f"{name}: code cell {index} does not parse: {exc}")
    return errors


def changelog_release_errors(text: str, expected: str) -> list[str]:
    errors: list[str] = []
    unreleased = re.search(
        r"^## \[Unreleased\][^\n]*\n(?P<body>.*?)(?=^## \[)",
        text,
        re.MULTILINE | re.DOTALL,
    )
    if unreleased is None:
        errors.append("CHANGELOG.md has no [Unreleased] section")
    else:
        body = unreleased.group("body").strip()
        if body not in {"", "_No changes yet._"}:
            errors.append("CHANGELOG.md [Unreleased] contains release payload")
    if not re.search(
        rf"^## \[{re.escape(expected)}\] — \d{{4}}-\d{{2}}-\d{{2}}$",
        text,
        re.MULTILINE,
    ):
        errors.append(
            f"CHANGELOG.md [{expected}] must use the actual YYYY-MM-DD release date"
        )
    return errors


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--expected-version")
    parser.add_argument("--release-mode", action="store_true")
    args = parser.parse_args()

    values = {
        "CITATION.cff": _cff_version(),
        "docs/conf.py": _docs_release(),
        "mlmm/_version.py": _module_version(),
        NOTEBOOK.name: _notebook_version(),
    }
    expected = str(args.expected_version or values["CITATION.cff"]).removeprefix("v")

    errors: list[str] = []
    mismatches = {name: value for name, value in values.items() if value != expected}
    if mismatches:
        details = ", ".join(f"{name}={value}" for name, value in mismatches.items())
        errors.append(f"expected version {expected}; mismatched: {details}")
    errors.extend(_check_landing_pages())
    errors.extend(_check_license())
    errors.extend(_check_release_workflow())
    errors.extend(_check_docs_deploy_workflow())
    errors.extend(_check_sdist_manifest())
    errors.extend(_check_ci_torch_pin())
    errors.extend(_check_pyproject_metadata())
    errors.extend(_check_citation_title())
    errors.extend(_check_notebook_metadata())
    if args.release_mode:
        errors.extend(changelog_release_errors(CHANGELOG.read_text(encoding="utf-8"), expected))

    if errors:
        print("[release-version] FAIL:")
        for error in errors:
            print(f"- {error}")
        return 1
    print(
        f"[release-version] OK: version {expected} ({', '.join(values)}); "
        f"landing headers use {LANDING_SUBSTITUTION}; license {EXPECTED_LICENSE}; "
        "release and docs workflows, sdist manifest, pyproject extras, citation title "
        "and Colab notebook metadata checked"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
