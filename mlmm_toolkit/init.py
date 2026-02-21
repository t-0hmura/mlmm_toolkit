"""Generate starter YAML templates for mlmm workflows."""

from __future__ import annotations

from pathlib import Path

import click


_ALL_TEMPLATE = """# Starter config for `mlmm all`
# Merge order: defaults < --config < CLI < --override-yaml

extract:
  radius: 2.6
  radius_het2het: 0.0

path_search:
  max_nodes: 10
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
  pressure: 1.0

dft:
  func_basis: wb97m-v/def2-tzvpd
  max_cycle: 100
  conv_tol: 1.0e-9
  grid_level: 3
"""


@click.command(help="Generate a starter YAML template for `mlmm all`.")
@click.option(
    "-o",
    "--out",
    "out_path",
    type=click.Path(path_type=Path, dir_okay=False),
    default=Path("mlmm_all.config.yaml"),
    show_default=True,
    help="Destination YAML file path.",
)
@click.option(
    "--force/--no-force",
    default=False,
    show_default=True,
    help="Overwrite the destination when it already exists.",
)
def cli(out_path: Path, force: bool) -> None:
    out_path = out_path.expanduser()
    if out_path.exists() and not force:
        raise click.ClickException(
            f"Refusing to overwrite existing file: {out_path}. Use --force to overwrite."
        )

    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(_ALL_TEMPLATE, encoding="utf-8")
    click.echo(f"[init] Wrote template: {out_path}")
    click.echo("[init] Next step: mlmm all --config " + str(out_path) + " --dry-run ...")

