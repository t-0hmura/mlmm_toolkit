# mlmm/domain/bond_summary.py
"""
bond-summary — Detect and summarize bond changes between molecular structures.

Usage:
    mlmm bond-summary -i A.xyz B.xyz [C.xyz ...]

Compares consecutive pairs of structures (A→B, B→C, …) and reports
covalent bonds formed and broken. Supports XYZ, PDB, and GJF formats.
"""

from __future__ import annotations

from pathlib import Path
from typing import List

import click

from mlmm.domain.bond_changes import compare_structures, summarize_changes


def _load_geom(path: str):
    """Load a geometry from XYZ, PDB, or GJF file."""
    p = Path(path)
    suffix = p.suffix.lower()
    if suffix == ".pdb":
        from pysisyphus.io.pdb import geom_from_pdb
        return geom_from_pdb(str(p))
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
@click.option(
    "--json/--no-json",
    "out_json",
    default=False,
    show_default=True,
    help="Write the JSON report to stdout; no file is created.",
)
def cli(inputs: tuple, extra_inputs: tuple, device: str, bond_factor: float, one_based: bool, out_json: bool) -> None:
    """Detect and summarize bond changes between consecutive structure pairs.

    \b
    Usage:
      mlmm bond-summary -i A.xyz -i B.xyz -i C.xyz
      mlmm bond-summary A.xyz B.xyz C.xyz
      mlmm bond-summary -i A.xyz B.xyz C.xyz
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

    comparisons_json: List[dict] = []
    n_failed = 0  # pairs whose comparison raised
    n_ok = 0      # pairs that compared successfully

    for idx in range(len(geoms) - 1):
        name1, g1 = geoms[idx]
        name2, g2 = geoms[idx + 1]

        if not out_json:
            click.echo(f"{'=' * 60}")
            click.echo(f"  {name1}  →  {name2}")
            click.echo(f"{'=' * 60}")

        try:
            result = compare_structures(g1, g2, device=device, bond_factor=bond_factor)
            if out_json:
                comparisons_json.append({
                    "structure_a": name1,
                    "structure_b": name2,
                    "bonds_formed": len(result.formed_covalent),
                    "bonds_broken": len(result.broken_covalent),
                })
            else:
                summary = summarize_changes(g2, result, one_based=one_based)
                click.echo(summary)
            n_ok += 1
        except Exception as e:
            # A dropped comparison must be visible: keep the per-pair stderr ERROR,
            # but also record it in the JSON envelope status + exit code below so a
            # `--json` consumer never sees status=ok with a silently missing pair.
            click.echo(f"  ERROR: {e}", err=True)
            n_failed += 1
            if out_json:
                comparisons_json.append({
                    "structure_a": name1,
                    "structure_b": name2,
                    "error": str(e),
                })

        if not out_json:
            click.echo()

    if out_json:
        import json as _json
        # Honest status: every pair ok → "ok"; some ok, some failed → "partial";
        # no pair ok → "failed". The exit code mirrors the status so non-JSON
        # callers and CI also see the failure.
        if n_failed == 0:
            status = "ok"
        elif n_ok > 0:
            status = "partial"
        else:
            status = "failed"
        # `--json` is an explicit machine-readable deliverable: it must always
        # reach stdout, even at default verbosity (which otherwise drops DETAIL
        # console output via the core.utils click.echo chokepoint). `force=True`
        # bypasses the narrative gate; see mlmm.core.utils._patch_click_echo.
        click.echo(
            _json.dumps(
                {"status": status, "comparisons": comparisons_json},
                indent=2, ensure_ascii=False,
            ),
            force=True,
        )

    if n_failed > 0:
        # Surface the failure via the process exit code (the JSON envelope above
        # already carries status=partial/failed for machine consumers).
        raise SystemExit(1)
