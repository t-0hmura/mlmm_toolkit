"""
Energy-profile utility: extract energies from XYZ trajectories and export Plotly figures / CSV.

Example:
    mlmm trj2fig -i traj.xyz -o energy.png energy.csv --unit kcal

For detailed documentation, see: docs/trj2fig.md
"""

from __future__ import annotations

import csv
import os
import time
from pathlib import Path
from typing import List, Optional, Sequence, Tuple

import click
from mlmm.core.output import emit
import plotly.graph_objs as go
from ase import Atoms
from ase.io import read
from pysisyphus.constants import AU2EV, AU2KCALPERMOL
from mlmm.io.xyz_trajectory import read_xyz_trajectory

AXIS_WIDTH = 3         # axis and tick thickness
FONT_SIZE = 18         # tick-label font size
AXIS_TITLE_SIZE = 20   # axis-title font size
LINE_WIDTH = 2         # curve width
MARKER_SIZE = 6        # marker size


def read_energies_xyz(fname: Path | str) -> List[float]:
    """Extract Hartree energies with the strict shared XYZ parser."""
    parsed = read_xyz_trajectory(fname, require_energies=True)
    return [float(value) for value in parsed["energies_ha"]]


def recompute_energies(
    traj_path: Path,
    charge: Optional[int],
    multiplicity: Optional[int],
    *,
    backend: str = "uma",
    backend_model: Optional[str] = None,
    precision: Optional[str] = None,
) -> List[float]:
    """Recalculate frame energies with the selected MLIP backend."""
    try:
        import torch
    except Exception as exc:
        raise RuntimeError("Energy recomputation requires torch.") from exc

    from mlmm.backends.mlmm_calc import _create_ml_backend
    from mlmm.core.utils import validate_charge_spin

    # Parse once with the strict shared reader before backend construction so
    # every frame has one atom identity/order.
    read_xyz_trajectory(traj_path)
    frames_obj = read(traj_path, index=":", format="xyz")
    frames = [frames_obj] if isinstance(frames_obj, Atoms) else list(frames_obj)
    if not frames:
        raise RuntimeError(f"No frames found in {traj_path}")

    q = int(charge if charge is not None else 0)
    mult = int(multiplicity if multiplicity is not None else 1)
    validate_charge_spin(
        frames[0].get_chemical_symbols(),
        q,
        mult,
        source=str(traj_path),
    )

    calc_cfg, _ = _resolved_backend_config(
        backend,
        backend_model=backend_model,
        precision=precision,
    )
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    backend_kwargs = {
        key: value
        for key, value in calc_cfg.items()
        if key
        in {
            "uma_model",
            "uma_precision",
            "orb_model",
            "orb_precision",
            "mace_model",
            "mace_dtype",
            "aimnet2_model",
        }
    }
    ml_backend = _create_ml_backend(
        str(calc_cfg["backend"]),
        model_charge=q,
        model_mult=mult,
        ml_device=device,
        **backend_kwargs,
    )

    energies_h: List[float] = []
    for atoms in frames:
        energy_ev, _, _ = ml_backend.eval(atoms, need_grad=False)
        energies_h.append(energy_ev / AU2EV)
    return energies_h


def _resolved_backend_config(
    backend: str,
    *,
    backend_model: Optional[str] = None,
    precision: Optional[str] = None,
) -> tuple[dict, dict]:
    """Resolve backend-specific model/precision and its public provenance."""
    from mlmm.backends import (
        apply_backend_model_to_calc_cfg,
        apply_precision_to_calc_cfg,
    )
    from mlmm.core.utils import calculator_provenance

    calc_cfg = {"backend": str(backend).strip().lower()}
    apply_backend_model_to_calc_cfg(calc_cfg, backend_model)
    apply_precision_to_calc_cfg(calc_cfg, precision)
    provenance = calculator_provenance(calc_cfg)
    return calc_cfg, {
        key: provenance[key]
        for key in ("mlip_backend", "mlip_model", "mlip_precision")
    }


def _parse_reference_spec(spec: str | None) -> str | int | None:
    """
    Normalize the reference specification:
      - "init" (case-insensitive) -> "init"
      - "none"/"null"             -> None
      - integer-like string       -> int
    """
    if spec is None:
        return "init"
    s = str(spec).strip()
    lower = s.lower()
    if lower in {"none", "null"}:
        return None
    if lower == "init":
        return "init"
    try:
        return int(s)
    except ValueError:
        raise ValueError(
            f'Invalid -r/--reference: {spec!r}. Use "init", "None", or an integer index.'
        )


def _resolve_reference_index(
    n_frames: int, ref_spec: str | int | None, reverse_x: bool
) -> Tuple[Optional[int], bool]:
    """
    Decide the reference index and whether to compute a ΔE series.

    Returns (reference_index or None, is_delta)
    """
    if ref_spec is None:
        return None, False  # absolute energies
    if ref_spec == "init":
        idx = 0 if not reverse_x else n_frames - 1
        return idx, True
    # integer index
    idx = int(ref_spec)
    if idx < 0 or idx >= n_frames:
        raise ValueError(f"Reference index {idx} out of range (0..{n_frames-1}).")
    return idx, True


def transform_series(
    energies_hartree: Sequence[float],
    ref_spec_raw: str | None,
    unit: str,
    reverse_x: bool,
) -> Tuple[List[float], str, bool]:
    """
    Compute the y-series and its axis label.

    Returns (values, ylabel, is_delta)
    """
    ref_spec = _parse_reference_spec(ref_spec_raw)
    ref_idx, is_delta = _resolve_reference_index(len(energies_hartree), ref_spec, reverse_x)

    scale = AU2KCALPERMOL if unit == "kcal" else 1.0
    if is_delta:
        base = energies_hartree[ref_idx]  # type: ignore[index]
        values = [float((e - base) * scale) for e in energies_hartree]
        ylabel = f"ΔE ({'kcal/mol' if unit == 'kcal' else 'hartree'})"
    else:
        values = [float(e * scale) for e in energies_hartree]
        ylabel = f"E ({'kcal/mol' if unit == 'kcal' else 'hartree'})"

    return values, ylabel, is_delta


def _axis_template() -> dict:
    return dict(
        showline=True,
        linewidth=AXIS_WIDTH,
        linecolor="#1C1C1C",
        mirror=True,
        ticks="inside",
        tickwidth=AXIS_WIDTH,
        tickcolor="#1C1C1C",
        tickfont=dict(size=FONT_SIZE, color="#1C1C1C"),
        gridcolor="lightgrey",
        gridwidth=0.5,
        zeroline=False,
    )


def build_figure(delta_or_abs: Sequence[float], ylabel: str, reverse_x: bool) -> go.Figure:
    """
    Build a Plotly figure without a title.
    """
    fig = go.Figure(
        go.Scatter(
            x=list(range(len(delta_or_abs)))),
    )
    fig.data[0].update(
        y=list(delta_or_abs),
        mode="lines+markers",
        marker=dict(size=MARKER_SIZE),
        line=dict(shape="spline", smoothing=1.0, width=LINE_WIDTH),
    )

    xaxis_conf = _axis_template() | {
        "title": dict(text="Frame", font=dict(size=AXIS_TITLE_SIZE, color="#1C1C1C"))
    }
    if reverse_x:
        xaxis_conf["autorange"] = "reversed"

    fig.update_layout(
        xaxis=xaxis_conf,
        yaxis=_axis_template() | {
            "title": dict(text=ylabel, font=dict(size=AXIS_TITLE_SIZE, color="#1C1C1C"))
        },
        plot_bgcolor="white",
        paper_bgcolor="white",
        margin=dict(l=80, r=40, t=40, b=80),
    )
    return fig


def save_outputs(
    outs: Sequence[Path],
    fig: Optional[go.Figure],
    energies: Sequence[float],
    values: Sequence[float],
    unit: str,
    is_delta: bool,
) -> None:
    """
    Write all requested outputs.
    """
    for out in outs:
        ext = out.suffix.lower()
        if ext == ".csv":
            write_csv(out, energies, values, unit, is_delta)
        elif ext == ".html":
            assert fig is not None
            fig.write_html(out)
            emit(f"[trj2fig] Saved figure -> {out}", detail=True)
        elif ext in {".png", ".jpg", ".jpeg", ".pdf", ".svg"}:
            assert fig is not None
            kw = {"engine": "kaleido"}
            if ext == ".png":
                kw["scale"] = 2  # high-resolution PNG
            fig.write_image(out, **kw)
            emit(f"[trj2fig] Saved figure -> {out}", detail=True)
        else:
            raise ValueError(f"Unsupported format: {ext}")


def write_csv(
    out: Path,
    energies_hartree: Sequence[float],
    series: Sequence[float],
    unit: str,
    is_delta: bool,
) -> None:
    """
    Save energies (hartree) and ΔE/E series to CSV.
    """
    colname = (f"delta_{unit}" if is_delta else f"energy_{unit}")
    with out.open("w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh)
        w.writerow(["frame", "energy_hartree", colname])
        for i, (eh, y) in enumerate(zip(energies_hartree, series)):
            w.writerow([i, f"{eh:.8f}", f"{y:.6f}"])
    emit(f"[trj2fig] Saved CSV -> {out}", detail=True)


def run_trj2fig(
    input_path: Path,
    outs: Sequence[Path],
    unit: str,
    reference: str,
    reverse_x: bool,
    charge: Optional[int] = None,
    multiplicity: Optional[int] = None,
    backend: str = "uma",
    backend_model: Optional[str] = None,
    precision: Optional[str] = None,
) -> dict:
    """Run trj2fig and return energies, outputs, and actual provenance."""
    traj = input_path.expanduser().resolve()
    if not traj.is_file():
        raise FileNotFoundError(traj)

    out_paths = [Path(o).expanduser().resolve() for o in outs]
    candidates = [traj, *out_paths]
    for index, left in enumerate(candidates):
        for right in candidates[index + 1 :]:
            same = left == right
            if not same and left.exists() and right.exists():
                same = os.path.samefile(left, right)
            if same:
                raise ValueError(
                    f"Input/output paths must not physically alias: {left} and {right}."
                )

    recomputed = charge is not None or multiplicity is not None
    if not recomputed:
        parsed = read_xyz_trajectory(traj, require_energies=True)
        energies = [float(value) for value in parsed["energies_ha"]]
        energy_provenance = list(parsed["energy_provenance"])
        provenance = {
            "mlip_backend": None,
            "mlip_model": None,
            "mlip_precision": None,
        }
    else:
        _, provenance = _resolved_backend_config(
            backend,
            backend_model=backend_model,
            precision=precision,
        )
        emit(
            "[trj2fig] Recomputing energies with "
            f"{provenance['mlip_backend']} "
            f"({provenance['mlip_model']}, {provenance['mlip_precision']}) ...",
            detail=True,
        )
        energies = recompute_energies(
            traj,
            charge,
            multiplicity,
            backend=backend,
            backend_model=backend_model,
            precision=precision,
        )
        energy_provenance = ["mlip-recomputed"] * len(energies)
    values, ylabel, is_delta = transform_series(energies, reference, unit, reverse_x)

    need_plot = any(Path(o).suffix.lower() != ".csv" for o in outs)
    fig = build_figure(values, ylabel, reverse_x) if need_plot else None

    save_outputs(out_paths, fig, energies, values, unit, is_delta)

    return {
        "energies_hartree": energies,
        "out_paths": out_paths,
        "energy_source": "mlip_recomputed" if recomputed else "trajectory_comment",
        "energy_provenance": energy_provenance,
        "energy_unit": "hartree",
        "charge": int(charge if charge is not None else 0) if recomputed else None,
        "multiplicity": (
            int(multiplicity if multiplicity is not None else 1)
            if recomputed
            else None
        ),
        **provenance,
    }


@click.command(
    name="trj2fig",
    help="Plot ΔE or E from an XYZ trajectory and export figure/CSV.",
    context_settings={"help_option_names": ["-h", "--help"]},
)
@click.option(
    "-i",
    "--input",
    "input_path",  # explicit internal argument name
    required=True,
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
    help="XYZ trajectory file",
)
@click.option(
    "-o",
    "--out",
    "outs",
    multiple=True,                      # allow repeating -o
    default=(),                         # default is empty (we inject the fallback later)
    type=click.Path(dir_okay=False, path_type=Path),
    help="Output file(s). You can repeat -o, and/or list extra filenames after options "
         "(.png/.jpg/.jpeg/.html/.svg/.pdf/.csv). If nothing is given, defaults to energy.png.",
)
@click.argument(
    "extra_outs",                        # also accept extra filenames provided positionally after options
    nargs=-1,
    type=click.Path(dir_okay=False, path_type=Path),
)
@click.option(
    "--unit",
    type=click.Choice(["kcal", "hartree"]),
    default="kcal", show_default=True,
    help="Energy unit.",
)
@click.option(
    "-r",
    "--reference",
    default="init", show_default=True,
    help='Reference: "init" (initial frame; last frame if --reverse-x), "None" (absolute E), or an integer index.',
)
@click.option(
    "-q",
    "--charge",
    type=int,
    default=None,
    help="Total charge. Recompute energies when supplied.",
)
@click.option(
    "-m",
    "--multiplicity",
    type=int,
    default=None,
    show_default="1",
    help="Spin multiplicity (2S+1). Recompute energies when supplied.",
)
@click.option(
    "--reverse-x/--no-reverse-x",
    "reverse_x",
    default=False,
    show_default=True,
    help="Reverse the x-axis (last frame on the left).",
)
@click.option(
    "-b",
    "--backend",
    type=click.Choice(["uma", "orb", "mace", "aimnet2"]),
    default="uma",
    show_default=True,
    help="MLIP backend used when energies are recomputed.",
)
@click.option(
    "--backend-model",
    default=None,
    show_default="the selected backend's own model",
    help="Model variant for the selected backend; defaults to its built-in model.",
)
@click.option(
    "--precision",
    type=click.Choice(["fp32", "fp64"], case_sensitive=False),
    default=None,
    show_default="per backend: uma fp32; orb, mace fp64",
    help="Backend-neutral precision used when energies are recomputed.",
)
@click.option(
    "--out-json/--no-out-json",
    default=False,
    show_default=True,
    help="Write machine-readable result.json next to the first output.",
)
def cli(
    input_path: Path,
    outs: Tuple[Path, ...],
    extra_outs: Tuple[Path, ...],
    unit: str,
    reference: str,
    charge: Optional[int],
    multiplicity: Optional[int],
    reverse_x: bool,
    backend: str,
    backend_model: Optional[str],
    precision: Optional[str],
    out_json: bool,
) -> None:
    time_start = time.perf_counter()
    # Combine outputs from -o with positional filenames that follow the options
    all_outs: List[Path] = list(outs) + list(extra_outs)
    if not all_outs:
        # Use the default when nothing is specified
        all_outs = [Path("energy.png")]
    try:
        info = run_trj2fig(
            input_path,
            all_outs,
            unit,
            reference,
            reverse_x,
            charge,
            multiplicity,
            backend=backend,
            backend_model=backend_model,
            precision=precision,
        )
    except (ValueError, RuntimeError) as exc:
        raise click.ClickException(str(exc)) from exc

    if out_json:
        from mlmm.core.utils import write_result_json

        energies = info["energies_hartree"]
        written_paths = info["out_paths"]
        out_dir = written_paths[0].parent if written_paths else Path.cwd()
        result_data = {
            "status": "ok",
            "n_frames": len(energies),
            "min_energy_hartree": float(min(energies)) if energies else None,
            "max_energy_hartree": float(max(energies)) if energies else None,
            "energy_source": info["energy_source"],
            "energy_provenance": info["energy_provenance"],
            "energy_unit": info["energy_unit"],
            "mlip_backend": info["mlip_backend"],
            "mlip_model": info["mlip_model"],
            "mlip_precision": info["mlip_precision"],
            "charge": info["charge"],
            "multiplicity": info["multiplicity"],
            # Preserve order and duplicate basenames across output
            # directories. ``files`` remains the legacy compatibility map.
            "output_files": [str(path) for path in written_paths],
            "files": {path.name: str(path) for path in written_paths},
        }
        write_result_json(out_dir, result_data, command="trj2fig")
    from mlmm.core.utils import format_elapsed

    emit(
        format_elapsed("[time] Elapsed Time for Trajectory Figure", time_start),
        narrative=True,
    )


if __name__ == "__main__":
    cli()
