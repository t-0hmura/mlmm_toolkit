# mlmm/energy_diagram.py

"""
Minimal energy-diagram utility from numeric inputs.

Examples:
    pdb2reaction energy-diagram -i 0 12.5 4.3 -o energy.png
    pdb2reaction energy-diagram -i "[-205.1, -190.4, -198.7]" --label-x R TS P
"""

from __future__ import annotations

import ast
import sys
from pathlib import Path
from typing import List, Sequence

import click

from .utils import build_energy_diagram, ensure_dir

_IMG_EXTS = {".png", ".jpg", ".jpeg", ".svg", ".pdf"}


def _normalize_image_path(path: Path) -> Path:
    p = Path(path)
    if p.suffix:
        ext = p.suffix.lower()
        if ext not in _IMG_EXTS:
            raise click.BadParameter(
                f"Unsupported output extension '{p.suffix}'. "
                f"Supported: {', '.join(sorted(_IMG_EXTS))}"
            )
        return p
    return p.with_suffix(".png")


def _collect_flag_values(
    argv: Sequence[str],
    names: Sequence[str],
    stop_flags: Sequence[str],
) -> List[str]:
    vals: List[str] = []
    names_set = set(names)
    stop_set = set(stop_flags)
    i = 0
    while i < len(argv):
        tok = argv[i]
        if tok in names_set:
            j = i + 1
            while j < len(argv) and argv[j] not in stop_set:
                vals.append(argv[j])
                j += 1
            i = j
        else:
            i += 1
    return vals


def _collect_inputs(default_values: Sequence[str]) -> List[str]:
    argv = list(sys.argv[1:])
    stop_flags = (
        "-i",
        "--input",
        "-o",
        "--output",
        "--label-x",
        "--label-y",
        "-h",
        "--help",
    )
    raw = _collect_flag_values(argv, ("-i", "--input"), stop_flags)
    if raw:
        return raw
    return list(default_values)


def _collect_label_x(default_values: Sequence[str]) -> List[str]:
    argv = list(sys.argv[1:])
    stop_flags = (
        "-i",
        "--input",
        "-o",
        "--output",
        "--label-x",
        "--label-y",
        "-h",
        "--help",
    )
    raw = _collect_flag_values(argv, ("--label-x",), stop_flags)
    if raw:
        return raw
    return list(default_values)


def _flatten_numeric_token(token: str) -> List[float]:
    txt = str(token).strip()
    if not txt:
        return []

    try:
        obj = ast.literal_eval(txt)
    except Exception:
        obj = txt

    def _to_float(x) -> float:
        try:
            return float(x)
        except Exception:
            raise click.BadParameter(
                f"Invalid numeric value in -i/--input: {x!r}"
            )

    if isinstance(obj, (int, float)):
        return [float(obj)]

    if isinstance(obj, (list, tuple)):
        out: List[float] = []
        for x in obj:
            out.append(_to_float(x))
        return out

    # Support comma-separated plain text like "0, 12.5, 4.3".
    if "," in txt:
        out: List[float] = []
        for part in txt.split(","):
            part = part.strip()
            if part:
                out.append(_to_float(part))
        return out

    return [_to_float(txt)]


def _parse_numeric_inputs(tokens: Sequence[str]) -> List[float]:
    vals: List[float] = []
    for tok in tokens:
        vals.extend(_flatten_numeric_token(tok))
    if len(vals) < 2:
        raise click.BadParameter(
            "Provide at least two numeric values with -i/--input."
        )
    return vals


def _flatten_label_token(token: str) -> List[str]:
    txt = str(token).strip()
    if not txt:
        return []

    try:
        obj = ast.literal_eval(txt)
    except Exception:
        obj = txt

    if isinstance(obj, (list, tuple)):
        return [str(x) for x in obj]
    if isinstance(obj, str) and "," in obj:
        return [x.strip() for x in obj.split(",") if x.strip()]
    return [str(obj)]


def _parse_label_x(tokens: Sequence[str]) -> List[str]:
    labels: List[str] = []
    for tok in tokens:
        labels.extend(_flatten_label_token(tok))
    return [str(x) for x in labels if str(x).strip()]


@click.command(
    name="energy-diagram",
    help=(
        "Plot an energy diagram from numeric inputs only. "
        "Use -i for values (multiple numbers or list-like string)."
    ),
    context_settings={
        "help_option_names": ["-h", "--help"],
        "ignore_unknown_options": True,
        "allow_extra_args": True,
    },
)
@click.option(
    "-i",
    "--input",
    "input_values",
    type=str,
    multiple=True,
    required=True,
    help=(
        "Numeric sequence. Accepts: "
        "-i 0 12.5 4.3  or  -i \"[0, 12.5, 4.3]\"  or repeated -i."
    ),
)
@click.option(
    "-o",
    "--output",
    "output_path",
    type=click.Path(path_type=Path, dir_okay=False),
    default=Path("energy_diagram.png"),
    show_default=True,
    help="Output image path.",
)
@click.option(
    "--label-x",
    "label_x",
    type=str,
    multiple=True,
    default=(),
    help=(
        "State labels on x-axis. Accepts: "
        "--label-x R TS P  or  --label-x \"['R','TS','P']\"."
    ),
)
@click.option(
    "--label-y",
    "label_y",
    type=str,
    default="ΔE (kcal/mol)",
    show_default=True,
    help="Y-axis label.",
)
def cli(
    input_values: Sequence[str],
    output_path: Path,
    label_x: Sequence[str],
    label_y: str,
) -> None:
    input_tokens = _collect_inputs(input_values)
    energies = _parse_numeric_inputs(input_tokens)

    label_tokens = _collect_label_x(label_x)
    labels = _parse_label_x(label_tokens)
    if labels:
        if len(labels) != len(energies):
            raise click.BadParameter(
                f"--label-x count ({len(labels)}) must match value count ({len(energies)})."
            )
        labels_use = labels
    else:
        labels_use = [f"S{i + 1}" for i in range(len(energies))]

    out_img = _normalize_image_path(Path(output_path).resolve())
    ensure_dir(out_img.parent)

    fig = build_energy_diagram(
        energies=energies,
        labels=labels_use,
        ylabel=label_y,
        baseline=True,
        showgrid=False,
    )
    fig.write_image(str(out_img), scale=2)
    click.echo(f"[energy-diagram] Saved -> {out_img}")
