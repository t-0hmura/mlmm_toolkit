# mlmm_toolkit/cli.py

import click

from mlmm_toolkit import __version__


class DefaultGroup(click.Group):
    """Custom Click group that defaults to 'all' command when no subcommand provided."""

    def __init__(self, *args, default: str | None = None, **kwargs):
        super().__init__(*args, **kwargs)
        self._default_cmd = default

    def parse_args(self, ctx, args):
        # Show help/version without inserting default command
        if any(a in ("-h", "--help", "--version") for a in args):
            return super().parse_args(ctx, args)

        if self._default_cmd is not None:
            if not args or args[0].startswith("-"):
                args.insert(0, self._default_cmd)
        return super().parse_args(ctx, args)

    def invoke(self, ctx):
        # Add a leading blank line for subcommands (except "all") to separate CLI tool logs.
        if ctx.invoked_subcommand and ctx.invoked_subcommand != "all":
            click.echo()
        return super().invoke(ctx)


from .all import cli as all_cmd
from .mm_parm import cli as mm_parm_cmd
from .scan import cli as scan_cmd
from .opt import cli as opt_cmd
from .path_opt import cli as path_opt_cmd
from .path_search import cli as path_search_cmd
from .tsopt import cli as tsopt_cmd
from .freq import cli as freq_cmd
from .irc import cli as irc_cmd
from .trj2fig import cli as trj2fig_cmd
from .add_elem_info import cli as add_elem_info_cmd
from .dft import cli as dft_cmd
from .scan2d import cli as scan2d_cmd
from .scan3d import cli as scan3d_cmd
from .oniom_export import cli_gaussian as oniom_gaussian_cmd
from .oniom_export import cli_orca as oniom_orca_cmd
from .define_layer import cli as define_layer_cmd
from .fix_altloc import cli as fix_altloc_cmd
from .energy_diagram import cli as energy_diagram_cmd

@click.group(
    cls=DefaultGroup,
    default="all",
    help="mlmm: Execute each step by subcommands.",
    context_settings={"help_option_names": ["-h", "--help"]},
)
@click.version_option(version=__version__, prog_name="mlmm")
def cli() -> None:
    pass


@click.command(
    name="extract",
    help="Extract a binding pocket.",
    context_settings={
        "ignore_unknown_options": True,
        "allow_extra_args": True,
        "help_option_names": [],
    },
)
@click.pass_context
def extract_cmd(ctx: click.Context) -> None:
    from . import extract as _extract_mod
    args = _extract_mod.parse_args(list(ctx.args))
    _extract_mod.extract(args)


cli.add_command(all_cmd, name="all")
cli.add_command(mm_parm_cmd, name="mm-parm")
cli.add_command(scan_cmd, name="scan")
cli.add_command(opt_cmd, name="opt")
cli.add_command(path_opt_cmd, name="path-opt")
cli.add_command(path_search_cmd, name="path-search")
cli.add_command(tsopt_cmd, name="tsopt")
cli.add_command(freq_cmd, name="freq")
cli.add_command(irc_cmd, name="irc")
cli.add_command(extract_cmd, name="extract")
cli.add_command(trj2fig_cmd, name="trj2fig")
cli.add_command(add_elem_info_cmd, name="add-elem-info")
cli.add_command(dft_cmd, name="dft")
cli.add_command(scan2d_cmd, name="scan2d")
cli.add_command(scan3d_cmd, name="scan3d")
cli.add_command(oniom_gaussian_cmd, name="oniom-gaussian")
cli.add_command(oniom_orca_cmd, name="oniom-orca")
cli.add_command(define_layer_cmd, name="define-layer")
cli.add_command(fix_altloc_cmd, name="fix-altloc")
cli.add_command(energy_diagram_cmd, name="energy-diagram")

# Disable pysisyphus logging
import logging
logging.disable(logging.CRITICAL)

# Filter noisy UMA/pydmf warnings that clutter CLI output
import warnings

warnings.filterwarnings(
    "ignore",
    category=UserWarning,
    message=r"var\(\): degrees of freedom is <= 0\. Correction should be strictly less than the reduction factor.*",
    module=r"fairchem\.core\.models\.uma\.escn_moe"
)
warnings.filterwarnings(
    "ignore",
    category=UserWarning,
    message=r"t_eval update skipped due to insufficient candidates",
    module=r"dmf"
)

if __name__ == "__main__":
    cli()
