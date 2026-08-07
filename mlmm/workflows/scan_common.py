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

from copy import deepcopy
import shutil
from pathlib import Path
from typing import Any, Callable, Dict, Mapping, Optional, Sequence, Tuple

import click

from pysisyphus.optimizers.LBFGS import LBFGS

from mlmm.core.defaults import LBFGS_KW, OPT_BASE_KW, THRESH_CHOICES
from mlmm.core.utils import apply_yaml_overrides


SCAN_THRESH_DEFAULT = "baker"


class OutputCollisionError(click.UsageError):
    """An output/input collision that must not create an error envelope."""


def prepare_scan_fixed_outputs(
    out_dir: Path,
    *,
    fixed_names: Sequence[str],
    protected_inputs: Sequence[Optional[Path]] = (),
) -> Path:
    """Invalidate fixed scan outputs after rejecting input collisions."""

    resolved = Path(out_dir).resolve()
    fixed = [resolved / name for name in fixed_names]
    fixed_resolved = {path.resolve(strict=False) for path in fixed}
    for protected in protected_inputs:
        if protected is None:
            continue
        if Path(protected).expanduser().resolve(strict=False) in fixed_resolved:
            raise OutputCollisionError(
                f"Input {protected} collides with a reserved scan output "
                f"under {resolved}."
            )
    resolved.mkdir(parents=True, exist_ok=True)
    for path in fixed:
        path.unlink(missing_ok=True)
    return resolved


def prepare_grid_scan_output(
    out_dir: Path,
    *,
    fixed_names: Sequence[str],
    protected_inputs: Sequence[Optional[Path]] = (),
) -> Tuple[Path, Path]:
    """Reset one grid-scan generation without touching unrelated files."""

    resolved = Path(out_dir).resolve()
    grid_dir = resolved / "grid"
    grid_resolved = grid_dir.resolve(strict=False)
    for protected in protected_inputs:
        if protected is None:
            continue
        protected_resolved = Path(protected).expanduser().resolve(strict=False)
        if (
            protected_resolved == grid_resolved
            or grid_resolved in protected_resolved.parents
        ):
            raise OutputCollisionError(
                f"Input {protected} collides with a reserved grid-scan output "
                f"under {resolved}."
            )

    resolved = prepare_scan_fixed_outputs(
        resolved,
        fixed_names=fixed_names,
        protected_inputs=protected_inputs,
    )
    if grid_dir.is_symlink() or grid_dir.is_file():
        grid_dir.unlink()
    elif grid_dir.is_dir():
        shutil.rmtree(grid_dir)
    grid_dir.mkdir()
    return resolved, grid_dir


def resolve_scan_optimizer_configs(
    yaml_cfg: Mapping[str, Any],
    *,
    opt_defaults: Mapping[str, Any] = OPT_BASE_KW,
    lbfgs_defaults: Mapping[str, Any] = LBFGS_KW,
    thresh: str,
    relax_max_cycles: int,
    is_param_explicit: Callable[[str], bool],
) -> Tuple[Dict[str, Any], Dict[str, Any]]:
    """Resolve scan optimizer settings once using Click parameter sources.

    The scan-specific ``baker`` threshold is a default layer, YAML is applied
    next, and only explicitly supplied CLI values replace the result.
    """

    opt_cfg = deepcopy(dict(opt_defaults))
    lbfgs_cfg = deepcopy(dict(lbfgs_defaults))
    opt_cfg["thresh"] = SCAN_THRESH_DEFAULT
    lbfgs_cfg["thresh"] = SCAN_THRESH_DEFAULT
    apply_yaml_overrides(
        yaml_cfg,
        [
            (opt_cfg, (("opt",),)),
            (lbfgs_cfg, (("lbfgs",), ("opt", "lbfgs"))),
        ],
    )
    if is_param_explicit("relax_max_cycles"):
        opt_cfg["max_cycles"] = int(relax_max_cycles)
        lbfgs_cfg["max_cycles"] = int(relax_max_cycles)
    if is_param_explicit("thresh"):
        opt_cfg["thresh"] = str(thresh)
        lbfgs_cfg["thresh"] = str(thresh)
    return opt_cfg, lbfgs_cfg


def build_scan_lbfgs_kwargs(
    lbfgs_cfg: Mapping[str, Any],
    opt_cfg: Mapping[str, Any],
    *,
    max_step_bohr: float,
    out_dir: Path,
    prefix: str,
) -> Dict[str, Any]:
    """Build LBFGS kwargs from the already-resolved scan configuration."""

    common = dict(opt_cfg)
    common["out_dir"] = str(out_dir)
    common["prefix"] = prefix
    args = {**dict(lbfgs_cfg), **common}
    args["max_step"] = min(float(lbfgs_cfg.get("max_step", 0.30)), max_step_bohr)
    return args


def make_scan_lbfgs(
    geom,
    lbfgs_cfg: Dict[str, Any],
    opt_cfg: Dict[str, Any],
    *,
    max_step_bohr: float,
    out_dir: Path,
    prefix: str,
) -> LBFGS:
    # Shared LBFGS factory for scan2d / scan3d (scan.py has a different
    # closure shape and is intentionally left inline). max_step is the LBFGS
    # cap in Bohr; max_cycles already comes from the effective config resolved
    # by ``resolve_scan_optimizer_configs``.
    args = build_scan_lbfgs_kwargs(
        lbfgs_cfg,
        opt_cfg,
        max_step_bohr=max_step_bohr,
        out_dir=out_dir,
        prefix=prefix,
    )
    return LBFGS(geom, **args)


def add_scan_common_options(
    *,
    out_dir_default: str,
    baseline_help: str,
    dump_help: str,
    max_step_help: str = "Maximum spacing between successive distance targets [Å].",
    relax_max_cycles_help: str = "Maximum L-BFGS cycles per biased relaxation (also used for preopt).",
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

    Used by `mlmm scan2d` and `mlmm scan3d`. Each common option has the same
    flag form, default, type, and help text in both commands.
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
            show_default="300.0",
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
