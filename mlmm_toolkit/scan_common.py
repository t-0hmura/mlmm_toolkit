# mlmm_toolkit/scan_common.py
"""Shared CLI options for scan, scan2d, and scan3d commands."""

from __future__ import annotations

from pathlib import Path
from typing import Callable, Dict, Tuple, Any

import click


def add_scan_common_options(
    *,
    out_dir_default: str,
    baseline_help: str,
    dump_help: str,
    max_step_help: str = "Maximum step size in either distance [Å].",
    thresh_default: str | None = "baker",
    max_step_size_default: float = 0.20,
    bias_k_default: float = 300.0,
    relax_max_cycles_default: int = 10000,
    opt_mode_default: str = "grad",
    dump_default: bool = False,
    preopt_default: bool = True,
    one_based_default: bool = True,
    include_baseline: bool = True,
    include_zmin_zmax: bool = True,
) -> Callable[[Callable], Callable]:
    """Attach the shared scan CLI options to a Click command.

    This decorator adds common options used across scan, scan2d, and scan3d commands
    to reduce code duplication.

    Parameters
    ----------
    out_dir_default : str
        Default output directory.
    baseline_help : str
        Help text for the --baseline option.
    dump_help : str
        Help text for the --dump option.
    max_step_help : str
        Help text for the --max-step-size option.
    thresh_default : str or None
        Default convergence threshold.
    max_step_size_default : float
        Default maximum step size.
    bias_k_default : float
        Default harmonic bias strength.
    relax_max_cycles_default : int
        Default maximum relaxation cycles.
    opt_mode_default : str
        Default optimization mode ("grad" or "hess").
    dump_default : bool
        Default value for --dump.
    preopt_default : bool
        Default value for --preopt.
    one_based_default : bool
        Default value for --one-based.
    include_baseline : bool
        Whether to include --baseline option.
    include_zmin_zmax : bool
        Whether to include --zmin/--zmax options.

    Returns
    -------
    Callable
        Decorator function.
    """
    thresh_note = f" Defaults to '{thresh_default}'." if thresh_default is not None else ""
    options = [
        click.option(
            "--one-based",
            "one_based",
            type=click.BOOL,
            default=one_based_default,
            show_default=True,
            help="Interpret (i,j) indices in --scan-lists as 1-based (default) or 0-based.",
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
            show_default=True,
            help="Harmonic well strength k [eV/Å^2].",
        ),
        click.option(
            "--relax-max-cycles",
            type=int,
            default=relax_max_cycles_default,
            show_default=True,
            help=(
                "Maximum optimizer cycles per grid relaxation. When explicitly provided, "
                "overrides opt.max_cycles from YAML."
            ),
        ),
        click.option(
            "--opt-mode",
            type=click.Choice(["grad", "hess"], case_sensitive=False),
            default=opt_mode_default,
            show_default=True,
            help="Relaxation mode: grad (=LBFGS) or hess (=RFO).",
        ),
        click.option(
            "--dump",
            type=click.BOOL,
            default=dump_default,
            show_default=True,
            help=dump_help,
        ),
        click.option(
            "--ref-pdb",
            type=click.Path(path_type=Path, exists=True, dir_okay=False),
            default=None,
            help="Reference PDB topology to use when the input is XYZ (keeps XYZ coordinates).",
        ),
        click.option(
            "--out-dir",
            type=str,
            default=out_dir_default,
            show_default=True,
            help="Base output directory.",
        ),
        click.option(
            "--thresh",
            type=str,
            default=thresh_default,
            show_default=False,
            help=(
                "Convergence preset (gau_loose|gau|gau_tight|gau_vtight|baker|never). "
                f"{thresh_note}"
            ),
        ),
        click.option(
            "--preopt",
            type=click.BOOL,
            default=preopt_default,
            show_default=True,
            help="Pre-optimize the initial structure without bias before the scan.",
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
                    help="Lower bound of color scale for plots (kcal/mol).",
                ),
                click.option(
                    "--zmax",
                    type=float,
                    default=None,
                    show_default=False,
                    help="Upper bound of color scale for plots (kcal/mol).",
                ),
            ]
        )

    def decorator(func: Callable) -> Callable:
        for opt in reversed(options):
            func = opt(func)
        return func

    return decorator
