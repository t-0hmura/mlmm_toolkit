# pdb2reaction/bond_summary.py
"""
bond-summary — Detect and summarize bond changes between molecular structures.

Usage:
    pdb2reaction bond-summary -i A.xyz B.xyz [C.xyz ...]

Compares consecutive pairs of structures (A→B, B→C, …) and reports
covalent bonds formed and broken. Supports XYZ, PDB, and GJF formats.
"""

from __future__ import annotations

from pathlib import Path
from typing import List

import click

from .bond_changes import compare_structures, summarize_changes


def _load_geom(path: str):
    """Load a geometry from XYZ, PDB, or GJF file."""
    p = Path(path)
    suffix = p.suffix.lower()
    if suffix == ".pdb":
        from pysisyphus.io.pdb import geom_from_pdb_str
        text = p.read_text(encoding="utf-8")
        return geom_from_pdb_str(text)
    elif suffix == ".gjf":
        from pysisyphus.io.gaussian import geom_from_gaussian_input
        return geom_from_gaussian_input(str(p))
    else:
        from pysisyphus.helpers import geom_from_xyz_file
        return geom_from_xyz_file(str(p))


@click.command("bond-summary", help="Detect bond changes between consecutive structures.")
@click.option(
    "-i", "--input", "inputs", multiple=True,
    help="Input structure files (XYZ/PDB/GJF). Repeat -i for each file.",
)
@click.argument("extra_inputs", nargs=-1, type=click.Path(exists=True))
@click.option(
    "--device", default="cpu", show_default=True,
    help="Compute device for distance calculations.",
)
@click.option(
    "--bond-factor", default=1.20, show_default=True,
    help="Scaling factor for covalent radii sum.",
)
@click.option(
    "--one-based/--zero-based", default=True, show_default=True,
    help="Use 1-based atom indices in output.",
)
def cli(inputs: tuple, extra_inputs: tuple, device: str, bond_factor: float, one_based: bool) -> None:
    """Detect and summarize bond changes between consecutive structure pairs.

    \b
    Usage:
      pdb2reaction bond-summary -i A.xyz -i B.xyz -i C.xyz
      pdb2reaction bond-summary A.xyz B.xyz C.xyz
      pdb2reaction bond-summary -i A.xyz B.xyz C.xyz
    """
    files: List[str] = list(inputs) + list(extra_inputs)
    if len(files) < 2:
        raise click.BadParameter("At least two input files are required.", param_hint="'-i'")

    geoms = []
    for f in files:
        p = Path(f)
        if not p.exists():
            raise click.FileError(f, hint="File not found.")
        geoms.append((p.name, _load_geom(f)))

    for idx in range(len(geoms) - 1):
        name1, g1 = geoms[idx]
        name2, g2 = geoms[idx + 1]

        click.echo(f"{'=' * 60}")
        click.echo(f"  {name1}  →  {name2}")
        click.echo(f"{'=' * 60}")

        try:
            result = compare_structures(g1, g2, device=device, bond_factor=bond_factor)
            summary = summarize_changes(g2, result, one_based=one_based)
            click.echo(summary)
        except AssertionError as e:
            click.echo(f"  ERROR: {e}", err=True)
        except Exception as e:
            click.echo(f"  ERROR: {e}", err=True)

        click.echo()
