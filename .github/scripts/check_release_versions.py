#!/usr/bin/env python3
"""Fail when release metadata, landing headers, or license identity disagree."""

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
        if _VERSION_LITERAL_RE.search(header):
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
    if args.release_mode:
        errors.extend(changelog_release_errors(CHANGELOG.read_text(encoding="utf-8"), expected))

    if errors:
        print("[release-version] FAIL:")
        for error in errors:
            print(f"- {error}")
        return 1
    print(
        f"[release-version] OK: version {expected} ({', '.join(values)}); "
        f"landing headers use {LANDING_SUBSTITUTION}; license {EXPECTED_LICENSE}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
