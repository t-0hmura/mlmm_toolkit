#!/usr/bin/env python3
"""Generate auto-derived CLI/YAML reference pages for docs."""

from __future__ import annotations

import argparse
import hashlib
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

import click
import yaml
from click.testing import CliRunner

REPO_ROOT = Path(__file__).resolve().parents[2]
DOCS_ROOT = REPO_ROOT / "docs"
REF_ROOT = DOCS_ROOT / "reference"
COMMANDS_ROOT = REF_ROOT / "commands"
YAML_REF_PATH = REF_ROOT / "yaml.md"
TOOL_NAME = "mlmm"
PACKAGE_NAME = "mlmm"

sys.path.insert(0, str(REPO_ROOT))

from mlmm.cli import cli as root_cli  # noqa: E402


_ALL_TEMPLATE = """# Starter config for `mlmm all`

calc:
  backend: uma              # ML backend: uma, orb, mace, aimnet2
  orb_model: orb_v3_conservative_omol  # ORB model name (when backend=orb)
  orb_precision: float64    # ORB precision default (when backend=orb; "float32-high" = TF32 matmul, also via --precision fp32; legacy "float32" alias accepted)
  mace_model: MACE-OMOL-0   # MACE model path or name (when backend=mace)
  mace_dtype: float64       # MACE dtype, e.g. float32 / float64 (when backend=mace)
  aimnet2_model: aimnet2    # AIMNet2 model name (when backend=aimnet2)
  embedcharge: false        # Enable xTB point-charge embedding correction
  embedcharge_step: 1.0e-3   # Numerical Hessian step for embedding correction (Å)
  xtb_cmd: xtb              # Path or command for the xTB executable
  xtb_acc: 0.2              # xTB SCF accuracy parameter
  xtb_workdir: tmp          # Working directory for xTB scratch files
  xtb_keep_files: false     # Keep xTB intermediate files after completion
  xtb_ncores: 4             # Number of CPU cores for xTB

extract:
  radius: 2.6
  radius_het2het: 0.0

path_search:
  max_nodes: 20
  max_cycles: 300

scan:
  max_step_size: 0.2
  bias_k: 300.0
  relax_max_cycles: 10000

tsopt:
  max_cycles: 10000

freq:
  max_write: 10
  amplitude_ang: 0.8
  n_frames: 20
  sort: value
  temperature: 298.15
  pressure_atm: 1.0

dft:
  func_basis: wb97m-v/def2-tzvpd
  max_cycle: 100
  conv_tol: 1.0e-9
  grid_level: 3
"""


# Ordered top-level sections this curated starter snapshot deliberately shows.
# The snapshot is intentionally NON-EXHAUSTIVE; the full schema lives in the
# hand-authored docs/yaml-reference.md.
_CURATED_SECTIONS: tuple[str, ...] = (
    "calc",
    "extract",
    "path_search",
    "scan",
    "tsopt",
    "freq",
    "dft",
)

# Value-free pointer from each curated scalar path to its single runtime owner.
# No expected value is hardcoded here: the generator resolves the live owner
# value and asserts the template scalar EQUALS it, so a stale snapshot value or
# an unused parity declaration fails generation. Owner kinds:
#   ("defaults", SYMBOL, KEY) -> mlmm.core.defaults.SYMBOL[KEY]
#   ("click", COMMAND, "--opt") -> the live `mlmm COMMAND --opt` Click default
_STARTER_OWNERS: dict[str, tuple[str, str, str]] = {
    "calc.backend": ("defaults", "MLMM_CALC_KW", "backend"),
    "calc.orb_model": ("defaults", "MLMM_CALC_KW", "orb_model"),
    "calc.orb_precision": ("defaults", "MLMM_CALC_KW", "orb_precision"),
    "calc.mace_model": ("defaults", "MLMM_CALC_KW", "mace_model"),
    "calc.mace_dtype": ("defaults", "MLMM_CALC_KW", "mace_dtype"),
    "calc.aimnet2_model": ("defaults", "MLMM_CALC_KW", "aimnet2_model"),
    "calc.embedcharge": ("defaults", "MLMM_CALC_KW", "embedcharge"),
    "calc.embedcharge_step": ("defaults", "MLMM_CALC_KW", "embedcharge_step"),
    "calc.xtb_cmd": ("defaults", "MLMM_CALC_KW", "xtb_cmd"),
    "calc.xtb_acc": ("defaults", "MLMM_CALC_KW", "xtb_acc"),
    "calc.xtb_workdir": ("defaults", "MLMM_CALC_KW", "xtb_workdir"),
    "calc.xtb_keep_files": ("defaults", "MLMM_CALC_KW", "xtb_keep_files"),
    "calc.xtb_ncores": ("defaults", "MLMM_CALC_KW", "xtb_ncores"),
    "extract.radius": ("click", "all", "--radius"),
    "extract.radius_het2het": ("click", "all", "--radius-het2het"),
    "path_search.max_nodes": ("defaults", "GS_KW", "max_nodes"),
    "path_search.max_cycles": ("defaults", "STOPT_KW", "max_cycles"),
    "scan.max_step_size": ("click", "scan", "--max-step-size"),
    "scan.bias_k": ("defaults", "BIAS_KW", "k"),
    "scan.relax_max_cycles": ("defaults", "OPT_BASE_KW", "max_cycles"),
    "tsopt.max_cycles": ("defaults", "OPT_BASE_KW", "max_cycles"),
    "freq.max_write": ("defaults", "FREQ_KW", "max_write"),
    "freq.amplitude_ang": ("defaults", "FREQ_KW", "amplitude_ang"),
    "freq.n_frames": ("defaults", "FREQ_KW", "n_frames"),
    "freq.sort": ("defaults", "FREQ_KW", "sort"),
    "freq.temperature": ("defaults", "THERMO_KW", "temperature"),
    "freq.pressure_atm": ("defaults", "THERMO_KW", "pressure_atm"),
    "dft.func_basis": ("defaults", "DFT_KW", "func_basis"),
    "dft.max_cycle": ("defaults", "DFT_KW", "max_cycle"),
    "dft.conv_tol": ("defaults", "DFT_KW", "conv_tol"),
    "dft.grid_level": ("defaults", "DFT_KW", "grid_level"),
}


@dataclass(frozen=True)
class RenderedFile:
    path: Path
    content: str


@dataclass(frozen=True)
class CommandDoc:
    name: str
    stem: str


def _collect_subcommands() -> list[str]:
    ctx = click.Context(root_cli)
    return sorted(root_cli.list_commands(ctx))


def _doc_stem(command_name: str) -> str:
    # Keep command spelling in headings/tables, but standardize filenames.
    return command_name.replace("-", "_")


def _collect_command_docs() -> list[CommandDoc]:
    docs: list[CommandDoc] = []
    used_stems: dict[str, str] = {}
    for name in _collect_subcommands():
        stem = _doc_stem(name)
        prev = used_stems.get(stem)
        if prev is not None and prev != name:
            raise RuntimeError(
                f"Command doc stem collision: '{prev}' and '{name}' both map to '{stem}'."
            )
        used_stems[stem] = name
        docs.append(CommandDoc(name=name, stem=stem))
    return docs


import re

_VERSION_LINE_RE = re.compile(r"^mlmm(?:-toolkit)? ver\. \S+\n", re.MULTILINE)
_PYSISRC_LINE_RE = re.compile(
    r"^Couldn't find configuration file\. Expected it at .*\n", re.MULTILINE
)


def _capture_help(command_name: str, *, advanced: bool) -> str:
    runner = CliRunner()
    args = [command_name, "--help-advanced"] if advanced else [command_name, "--help"]
    # The CLI start-header guard inspects sys.argv to detect a help/version
    # request. Under CliRunner the help flag lives in `args`, not sys.argv, so
    # without this the banner / [command] / [mode] lines leak into the captured
    # reference. Point sys.argv at the help invocation so the guard suppresses
    # them (real-invocation CLI behavior is unchanged).
    saved_argv = sys.argv
    sys.argv = [TOOL_NAME, *args]
    try:
        result = runner.invoke(
            root_cli,
            args,
            catch_exceptions=False,
            prog_name=TOOL_NAME,
        )
    finally:
        sys.argv = saved_argv
    if result.exit_code != 0:
        raise RuntimeError(
            f"Failed to collect help for '{TOOL_NAME} {command_name}' "
            f"(advanced={advanced}):\n{result.output}"
        )
    if "[Unavailable]" in result.output:
        raise RuntimeError(
            f"Cannot generate help for '{TOOL_NAME} {command_name}': the lazy "
            "subcommand could not be imported. Install the repository's "
            "development/runtime dependencies and retry."
        )
    # Strip version line and pysisyphus config warning for reproducibility
    text = _VERSION_LINE_RE.sub("", result.output)
    text = _PYSISRC_LINE_RE.sub("", text)
    return text.rstrip() + "\n"


def _render_command_page(command_name: str, help_text: str) -> str:
    return (
        f"# `{TOOL_NAME} {command_name}`\n\n"
        "```text\n"
        f"{help_text.rstrip()}\n"
        "```\n"
    )


def _render_commands_index(command_docs: list[CommandDoc]) -> str:
    toctree_body = "\n".join(doc.stem for doc in command_docs)
    rows = "\n".join(
        f"| `{TOOL_NAME} {doc.name}` | [{doc.name}]({doc.stem}.md) |" for doc in command_docs
    )
    return (
        "# CLI Command Reference\n\n"
        "These pages are the generated `--help` reference for each command. For usage "
        "guidance and worked examples, see the command's narrative guide in the main "
        "documentation (e.g. [all](../../all.md), [opt](../../opt.md)).\n\n"
        "```{toctree}\n"
        ":maxdepth: 1\n"
        ":hidden:\n\n"
        f"{toctree_body}\n"
        "```\n\n"
        "| Command | Page |\n"
        "|---|---|\n"
        f"{rows}\n"
    )


def _iter_scalar_items(data: dict, prefix: str = "") -> Iterable[tuple[str, object]]:
    for key, value in data.items():
        full_key = f"{prefix}.{key}" if prefix else str(key)
        if isinstance(value, dict):
            yield from _iter_scalar_items(value, prefix=full_key)
        else:
            yield full_key, value


def _scalar_equal(template_value: object, owner_value: object) -> bool:
    if isinstance(template_value, bool) or isinstance(owner_value, bool):
        return template_value is owner_value or template_value == owner_value
    return template_value == owner_value


def _resolve_owner_values(root_cli) -> dict[str, tuple[str, object]]:
    """Resolve each curated scalar's runtime owner label and live value."""
    from mlmm.core import defaults as mlmm_defaults  # noqa: E402

    ctx = root_cli.make_context(TOOL_NAME, [], resilient_parsing=True)
    click_cache: dict[str, dict[str, object]] = {}

    def _click_default(command_name: str, opt_name: str) -> object:
        cache = click_cache.get(command_name)
        if cache is None:
            command = root_cli.get_command(ctx, command_name)
            if command is None or (command.help or "").startswith("[Unavailable]"):
                raise RuntimeError(
                    f"live command {command_name!r} unavailable for owner resolution"
                )
            cache = {}
            for param in command.params:
                if isinstance(param, click.Option):
                    for opt in param.opts:
                        if opt.startswith("--"):
                            cache[opt] = param.default
            click_cache[command_name] = cache
        if opt_name not in cache:
            raise RuntimeError(f"owner option {opt_name!r} not found on {command_name!r}")
        return cache[opt_name]

    resolved: dict[str, tuple[str, object]] = {}
    try:
        for path, owner in _STARTER_OWNERS.items():
            kind = owner[0]
            if kind == "defaults":
                _, symbol, key = owner
                container = getattr(mlmm_defaults, symbol, None)
                if not isinstance(container, dict) or key not in container:
                    raise RuntimeError(f"owner {symbol}[{key!r}] missing in defaults")
                resolved[path] = (f'`{symbol}["{key}"]`', container[key])
            elif kind == "click":
                _, command_name, opt_name = owner
                resolved[path] = (
                    f"`mlmm {command_name} {opt_name}` default",
                    _click_default(command_name, opt_name),
                )
            else:  # pragma: no cover - guarded by the static owner map
                raise RuntimeError(f"unknown owner kind {kind!r} for {path!r}")
    finally:
        ctx.close()
    return resolved


def _validate_starter_snapshot(
    template_data: dict, owner_values: dict[str, tuple[str, object]]
) -> None:
    """Assert exact bidirectional scalar-path coverage and per-scalar parity."""
    template_items = dict(_iter_scalar_items(template_data))
    scalar_paths = set(template_items)
    errors: list[str] = []

    for path in sorted(scalar_paths):
        if path not in owner_values:
            errors.append(f"ownerless starter scalar: '{path}' declares no runtime owner")
    for path in sorted(owner_values):
        if path not in scalar_paths:
            errors.append(
                f"stale owner declaration: '{path}' is not present in the starter template"
            )
    for path in sorted(scalar_paths & set(owner_values)):
        label, owner_value = owner_values[path]
        template_value = template_items[path]
        if not _scalar_equal(template_value, owner_value):
            errors.append(
                f"starter value drift at '{path}': template={template_value!r} "
                f"!= runtime owner {label} = {owner_value!r}"
            )

    top = list(template_data.keys()) if isinstance(template_data, dict) else []
    if top != list(_CURATED_SECTIONS):
        errors.append(
            f"curated sections changed: {top} != {list(_CURATED_SECTIONS)}"
        )

    if errors:
        raise RuntimeError(
            "[yaml-reference] curated starter snapshot parity failed:\n"
            + "\n".join(errors)
        )


def _render_yaml_reference(owner_values: dict[str, tuple[str, object]]) -> str:
    template_data = yaml.safe_load(_ALL_TEMPLATE) or {}
    _validate_starter_snapshot(template_data, owner_values)

    top_keys = list(template_data.keys()) if isinstance(template_data, dict) else []
    top_rows = "\n".join(f"| `{k}` |" for k in top_keys)
    scalar_rows = "\n".join(
        f"| `{k}` | `{type(v).__name__}` | `{v!r}` | {owner_values[k][0]} |"
        for k, v in _iter_scalar_items(template_data)
    )
    digest = hashlib.sha256(_ALL_TEMPLATE.encode("utf-8")).hexdigest()[:12]

    scan_schema = """```yaml
# scan (1D staged)
one_based: false
stages:
  - - [1, 2, 1.65]
  - - [2, 3, 2.30]

# scan2d / scan3d
one_based: false
pairs:
  - [1, 2, 1.40, 2.20]
  - [2, 3, 1.20, 2.00]
  - [3, 4, 1.00, 1.80]  # required only for scan3d
```"""

    return (
        "# Curated `mlmm all` Starter Snapshot\n\n"
        "This page is a **curated, non-exhaustive** starter snapshot for "
        "`mlmm all`. It shows a common subset of keys whose values are pinned to "
        "(and equal) their runtime owners; it is **not** the full configuration "
        "schema. For every configurable section and option, see the "
        "[YAML Reference](../yaml-reference.md).\n\n"
        f"- Source template: `.github/scripts/generate_reference.py::_ALL_TEMPLATE`\n"
        f"- Template digest: `{digest}`\n\n"
        "## Included Sections\n\n"
        "| Section |\n"
        "|---|\n"
        f"{top_rows}\n\n"
        "## Starter Template\n\n"
        "```yaml\n"
        f"{_ALL_TEMPLATE.strip()}\n"
        "```\n\n"
        "## Scalar Defaults\n\n"
        "Each scalar is pinned to (and equals) the runtime owner shown.\n\n"
        "| Key | Type | Default | Runtime owner |\n"
        "|---|---|---|---|\n"
        f"{scalar_rows}\n\n"
        "## Scan Spec Shapes\n\n"
        "Accepted by `scan`, `scan2d`, and `scan3d` with `-s/--scan-lists`.\n\n"
        f"{scan_schema}\n"
    )


def _render() -> list[RenderedFile]:
    command_docs = _collect_command_docs()
    rendered: list[RenderedFile] = []

    for doc in command_docs:
        try:
            help_text = _capture_help(doc.name, advanced=True)
        except RuntimeError:
            help_text = _capture_help(doc.name, advanced=False)
        rendered.append(
            RenderedFile(
                path=COMMANDS_ROOT / f"{doc.stem}.md",
                content=_render_command_page(doc.name, help_text),
            )
        )

    rendered.append(
        RenderedFile(
            path=COMMANDS_ROOT / "index.md",
            content=_render_commands_index(command_docs),
        )
    )
    owner_values = _resolve_owner_values(root_cli)
    rendered.append(
        RenderedFile(path=YAML_REF_PATH, content=_render_yaml_reference(owner_values))
    )
    return rendered


def _check_or_write(rendered: list[RenderedFile], *, check: bool) -> int:
    expected_paths = {item.path.resolve() for item in rendered}
    existing_paths = set(COMMANDS_ROOT.glob("*.md"))
    stale_extra = sorted(p.resolve() for p in existing_paths if p.resolve() not in expected_paths)

    stale: list[Path] = []
    for item in rendered:
        current = item.path.read_text(encoding="utf-8") if item.path.exists() else None
        if current != item.content:
            stale.append(item.path.resolve())

    if check:
        if stale or stale_extra:
            print("Reference files are out of date. Run: python .github/scripts/generate_reference.py")
            for path in sorted(stale):
                print(f"  MISSING/DIFF: {path.relative_to(REPO_ROOT)}")
            for path in stale_extra:
                print(f"  EXTRA: {path.relative_to(REPO_ROOT)}")
            return 1
        print("Reference files are up to date.")
        return 0

    for item in rendered:
        item.path.parent.mkdir(parents=True, exist_ok=True)
        item.path.write_text(item.content, encoding="utf-8")

    for path in stale_extra:
        path.unlink()

    print(f"Generated {len(rendered)} files under {REF_ROOT.relative_to(REPO_ROOT)}")
    return 0


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true", help="Fail when generated files are stale.")
    args = parser.parse_args()
    rendered = _render()
    return _check_or_write(rendered, check=args.check)


if __name__ == "__main__":
    raise SystemExit(main())
