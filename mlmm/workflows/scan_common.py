"""Shared CLI options for scan2d and scan3d commands.

scan.py is intentionally NOT routed through this factory: it has several
scan-specific quirks (the `--opt-mode` compat option, `--relax-max-cycles`
acting as an alias for `--max-cycles` with default=None, `--thresh`
default=None) that the factory cannot express without per-call-site
branching that would defeat the deduplication purpose.

The factory below is calibrated to scan2d / scan3d, where the 11 common
options share identical types, defaults, and help text. Only the
per-command help phrasing for `--dump` / `--baseline` / `--out-dir`
default differs, and that is parameterised.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any, Callable, Dict

import click

from pysisyphus.optimizers.LBFGS import LBFGS

from mlmm.core.defaults import THRESH_CHOICES


def make_scan_lbfgs(
    geom,
    lbfgs_cfg: Dict[str, Any],
    opt_cfg: Dict[str, Any],
    *,
    max_step_bohr: float,
    relax_max_cycles: int,
    out_dir: Path,
    prefix: str,
) -> LBFGS:
    # Shared LBFGS factory for scan2d / scan3d (scan.py has a different
    # closure shape and is intentionally left inline). max_step is the LBFGS
    # cap in Bohr; max_cycles is overridden per relaxation via
    # --relax-max-cycles.
    common = dict(opt_cfg)
    common["out_dir"] = str(out_dir)
    common["prefix"] = prefix
    args = {**lbfgs_cfg, **common}
    args["max_step"] = min(float(lbfgs_cfg.get("max_step", 0.30)), max_step_bohr)
    args["max_cycles"] = int(relax_max_cycles)
    return LBFGS(geom, **args)


def add_scan_common_options(
    *,
    out_dir_default: str,
    baseline_help: str,
    dump_help: str,
    max_step_help: str = "Maximum spacing between successive distance targets [Å].",
    relax_max_cycles_help: str = "Maximum LBFGS cycles per biased relaxation (also used for preopt).",
    preopt_help: str = "Run an unbiased pre-optimization.",
    thresh_default: str = "baker",
    max_step_size_default: float = 0.20,
    # Flip default to None so YAML `bias.k` is not silently clobbered by the
    # CLI default. Use-sites do `if bias_k is not None: bias_cfg["k"] = ...`,
    # so None means "fall through to YAML/BIAS_KW".
    bias_k_default: float | None = None,
    relax_max_cycles_default: int = 10000,
    one_based_help: str = "Interpret (i,j) indices in --scan-lists as 1-based (default) or 0-based.",
    include_baseline: bool = True,
    include_zmin_zmax: bool = True,
) -> Callable[[Callable], Callable]:
    """Attach the 9–12 shared scan CLI options to a Click command.

    Used by `mlmm scan2d` and `mlmm scan3d`. Each common option keeps its
    exact pre-wire surface (flag form, default value, type, help text) so
    that the only diff downstream of this refactor is the relative position
    in `--help` (= golden migrate, not a behavior change).
    """
    options = [
        click.option(
            "--one-based/--zero-based",
            "one_based",
            default=True,
            show_default=True,
            help=one_based_help,
        ),
        click.option(
            "--max-step-size",
            type=float,
            default=max_step_size_default,
            show_default=True,
            help=max_step_help,
        ),
        click.option(
            "--bias-k",
            type=float,
            default=bias_k_default,
            show_default=False,
            help=(
                "Harmonic well strength k [eV/Å^2]. "
                "Defaults to YAML bias.k (BIAS_KW['k']=300 in defaults.py) when omitted; "
                "explicit CLI value overrides YAML."
            ),
        ),
        click.option(
            "--relax-max-cycles",
            type=int,
            default=relax_max_cycles_default,
            show_default=True,
            help=relax_max_cycles_help,
        ),
        click.option(
            "--dump/--no-dump",
            "dump",
            default=False,
            show_default=True,
            help=dump_help,
        ),
        click.option(
            "-o", "--out-dir",
            type=str,
            default=out_dir_default,
            show_default=True,
            help="Base output directory.",
        ),
        click.option(
            "--thresh",
            type=click.Choice(THRESH_CHOICES, case_sensitive=False),
            default=thresh_default,
            show_default=True,
            help="Convergence preset.",
        ),
        click.option(
            "--ref-pdb",
            type=click.Path(path_type=Path, exists=True, dir_okay=False),
            default=None,
            help="Reference PDB topology to use when --input is XYZ (keeps XYZ coordinates).",
        ),
        click.option(
            "--preopt/--no-preopt",
            "preopt",
            default=False,
            show_default=True,
            help=preopt_help,
        ),
    ]
    if include_baseline:
        options.append(
            click.option(
                "--baseline",
                type=click.Choice(["min", "first"]),
                default="min",
                show_default=True,
                help=baseline_help,
            )
        )
    if include_zmin_zmax:
        options.extend(
            [
                click.option(
                    "--zmin",
                    type=float,
                    default=None,
                    show_default=False,
                    help="Lower bound of the color scale (kcal/mol).",
                ),
                click.option(
                    "--zmax",
                    type=float,
                    default=None,
                    show_default=False,
                    help="Upper bound of the color scale (kcal/mol).",
                ),
            ]
        )

    def decorator(func: Callable) -> Callable:
        for opt in reversed(options):
            func = opt(func)
        return func

    return decorator
