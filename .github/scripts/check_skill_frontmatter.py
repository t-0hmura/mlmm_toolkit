#!/usr/bin/env python3
"""Validate every distributable skill's YAML frontmatter."""

from __future__ import annotations

import re
from pathlib import Path

import yaml


REPO_ROOT = Path(__file__).resolve().parents[2]
SKILLS_DIR = REPO_ROOT / "skills"
ALLOWED_KEYS = {"name", "description", "license", "allowed-tools", "metadata"}
NAME_RE = re.compile(r"^[a-z0-9]+(?:-[a-z0-9]+)*$")
FRONTMATTER_RE = re.compile(r"\A---\r?\n(.*?)\r?\n---(?:\r?\n|\Z)", re.DOTALL)


def _validate(path: Path) -> list[str]:
    errors: list[str] = []
    match = FRONTMATTER_RE.match(path.read_text(encoding="utf-8"))
    if match is None:
        return [f"{path.relative_to(REPO_ROOT)}: missing or malformed YAML frontmatter"]
    try:
        data = yaml.safe_load(match.group(1))
    except yaml.YAMLError as exc:
        return [f"{path.relative_to(REPO_ROOT)}: invalid YAML: {exc}"]
    if not isinstance(data, dict):
        return [f"{path.relative_to(REPO_ROOT)}: frontmatter must be a mapping"]

    unexpected = set(data) - ALLOWED_KEYS
    if unexpected:
        errors.append(f"unexpected keys: {', '.join(sorted(unexpected))}")
    name = data.get("name")
    if not isinstance(name, str) or not name.strip():
        errors.append("name must be a non-empty string")
    elif not NAME_RE.fullmatch(name) or len(name) > 64:
        errors.append("name must be <=64 characters of lowercase hyphen-case")
    elif name != path.parent.name:
        errors.append(f"name {name!r} must match directory {path.parent.name!r}")

    description = data.get("description")
    if not isinstance(description, str) or not description.strip():
        errors.append("description must be a non-empty string")
    elif len(description.strip()) > 1024:
        errors.append("description exceeds 1024 characters")
    elif "<" in description or ">" in description:
        errors.append("description cannot contain angle brackets")

    rel = path.relative_to(REPO_ROOT)
    return [f"{rel}: {message}" for message in errors]


def main() -> int:
    paths = sorted(SKILLS_DIR.glob("*/SKILL.md"))
    errors = [error for path in paths for error in _validate(path)]
    if errors:
        print(f"[skill-frontmatter] FAIL: {len(errors)} error(s)")
        for error in errors:
            print(error)
        return 1
    print(f"[skill-frontmatter] OK: {len(paths)} skills")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
