#!/usr/bin/env python
# *- coding: utf-8 -*
# ---------------------------------------------------------------------
"""
trj2fig  –  Energy-profile utility for XYZ trajectories
======================================================

• Extract Hartree energies from the 2-line comment of each XYZ frame  
• Convert to ΔE relative to a reference frame (kcal mol⁻¹ or hartree)  
• Generate a aligned Plotly figure 
  (HTML / PNG / SVG / PDF, thick axes, ticks, fonts, spline curve)  
• Optionally export the ΔE table as CSV  
• Optionally write the highest internal-peak structure as .xyz or .pdb  
  (if no internal peak exists, the global maximum is used)
• `--reverse-x` flips the x-axis so the last frame appears leftmost

EXAMPLES
--------

Plot & high-res PNG **+** export peak frame (x-axis reversed):

    trj2fig -i traj.xyz -o energy.png --output-peak ts.pdb --reverse-x

Save CSV only:

    trj2fig -i traj.xyz -o energy.csv -r 5 --unit hartree

Default reference frame is 0; change it with “-r IDX”.

Author
------
Takuto Ohmura – 7 Jun 2025
"""
from __future__ import annotations

import argparse
import csv
import re
from pathlib import Path
from typing import List, Tuple

import plotly.graph_objs as go
from ase.io import read, write
from pysisyphus.constants import AU2KCALPERMOL

AXIS_WIDTH      = 3   # axis & tick thickness
FONT_SIZE       = 18  # tick-label font size
AXIS_TITLE_SIZE = 20  # axis-label font size
LINE_WIDTH      = 2   # ΔE curve width
MARKER_SIZE     = 6   # marker size


# ---------------------------------------------------------------------
#  File helpers
# ---------------------------------------------------------------------
def read_energies_xyz(fname: Path | str) -> List[float]:
    """Return a list of Hartree energies extracted from 2-line comments."""
    energies: List[float] = []
    with open(fname, encoding="utf-8") as fh:
        while (hdr := fh.readline()):
            try:
                nat = int(hdr.strip())
            except ValueError:                         # non-XYZ header reached
                break
            comment = fh.readline().strip()
            m = re.search(r"(-?\d+(?:\.\d+)?)", comment)
            if not m:
                raise RuntimeError(f"Energy not found: {comment}")
            energies.append(float(m.group(1)))
            for _ in range(nat):                       # skip coordinates
                fh.readline()
    if not energies:
        raise RuntimeError(f"No energy data in {fname}")
    return energies


def delta_e(e: List[float], ref: int, unit: str) -> Tuple[List[float], str]:
    """Convert absolute energies to ΔE and return y-axis label."""
    e0 = e[ref]
    if unit == "kcal":
        d = [(x - e0) * AU2KCALPERMOL for x in e]
        label = "ΔE (kcal/mol)"
    elif unit == "hartree":
        d = [(x - e0) for x in e]
        label = "ΔE (hartree)"
    else:
        raise ValueError(unit)
    return d, label


def write_csv(out: Path, e: List[float], d: List[float], unit: str) -> None:
    """Save energies and ΔE to CSV."""
    with out.open("w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh)
        w.writerow(["frame", "energy_hartree", f"delta_{unit}"])
        for i, (ei, di) in enumerate(zip(e, d)):
            w.writerow([i, f"{ei:.8f}", f"{di:.6f}"])
    print(f"[trj2fig] CSV → {out}")


# ---------------------------------------------------------------------
#  Plotting
# ---------------------------------------------------------------------
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


def plot_energy(
    delta: List[float],
    ylabel: str,
    out_img: Path | None,
    reverse_x: bool,
) -> None:
    """Generate a title-less Plotly figure and save it."""
    fig = go.Figure(
        go.Scatter(
            x=list(range(len(delta))),
            y=delta,
            mode="lines+markers",
            marker=dict(size=MARKER_SIZE),
            line=dict(shape="spline", smoothing=1.0, width=LINE_WIDTH),
        )
    )

    # Build axis layout
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

    if out_img:
        ext = out_img.suffix.lower()
        if ext == ".html":
            fig.write_html(out_img)
        elif ext in {".png", ".jpg", ".jpeg", ".pdf", ".svg"}:
            kw = {"engine": "kaleido"}
            if ext == ".png":
                kw["scale"] = 2                 # hi-res PNG
            fig.write_image(out_img, **kw)
        else:
            raise ValueError(f"Unsupported format: {ext}")
        print(f"[trj2fig] Figure → {out_img}")


# ---------------------------------------------------------------------
#  Peak utilities
# ---------------------------------------------------------------------
def _internal_peaks(vals: List[float]) -> List[int]:
    return [
        i
        for i in range(1, len(vals) - 1)
        if vals[i] > vals[i - 1] and vals[i] > vals[i + 1]
    ]


def find_peak_idx(energies: List[float]) -> int:
    peaks = _internal_peaks(energies)
    return (
        max(peaks, key=energies.__getitem__)
        if peaks
        else int(max(range(len(energies)), key=energies.__getitem__))
    )


def write_peak(trj: Path, idx: int, out: Path) -> None:
    if out.suffix.lower() not in {".xyz", ".pdb"}:
        raise ValueError("--output-peak must end with .xyz or .pdb")
    write(out, read(trj, index=idx, format="xyz"))
    print(f"[trj2fig] Peak frame {idx} → {out}")


# ---------------------------------------------------------------------
#  CLI
# ---------------------------------------------------------------------
def parse_cli() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        prog="trj2fig",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description="Plot ΔE from an XYZ trajectory, write a figure/CSV, and export the peak structure (no title).",
    )
    p.add_argument("-i", "--input", required=True, help="XYZ trajectory file")
    p.add_argument(
        "-o",
        "--out",
        default="energy.html",
        help="Output (.html/.png/.svg/.pdf/.csv)",
    )
    p.add_argument("--unit", choices=["kcal", "hartree"], default="kcal", help="ΔE unit")
    p.add_argument(
        "-r", "--reference", type=int, default=0, help="Reference frame index"
    )
    p.add_argument(
        "--output-peak",
        metavar="PEAK.xyz|PEAK.pdb",
        help="Write the highest peak frame",
    )
    p.add_argument(
        "--reverse-x",
        action="store_true",
        help="Reverse the x-axis (last frame on the left)",
    )
    return p.parse_args()


def main() -> None:
    args = parse_cli()
    traj = Path(args.input).expanduser().resolve()
    if not traj.is_file():
        raise FileNotFoundError(traj)
    out_path = Path(args.out).expanduser().resolve()

    energies = read_energies_xyz(traj)
    delta, ylabel = delta_e(energies, args.reference, args.unit)

    if out_path.suffix.lower() == ".csv":
        write_csv(out_path, energies, delta, args.unit)
    else:
        plot_energy(delta, ylabel, out_img=out_path, reverse_x=args.reverse_x)

    if args.output_peak:
        write_peak(
            traj, find_peak_idx(energies), Path(args.output_peak).expanduser()
        )


if __name__ == "__main__":
    main()
