# mlmm_toolkit/cli.py

import importlib
import logging
import warnings

import click

from mlmm_toolkit import __version__


_BOOL_TRUE_LITERALS = {"1", "true", "t", "yes", "y", "on"}
_BOOL_FALSE_LITERALS = {"0", "false", "f", "no", "n", "off"}

_LAZY_SUBCOMMANDS: dict[str, tuple[str, str, str]] = {
    "all": (".all", "cli", "End-to-end workflow (extract -> MEP -> TS -> IRC -> freq -> DFT)."),
    "mm-parm": (".mm_parm", "cli", "Generate Amber parm7/rst7 topology files."),
    "scan": (".scan", "cli", "Run staged 1D scan with harmonic restraints."),
    "opt": (".opt", "cli", "Optimize one structure."),
    "path-opt": (".path_opt", "cli", "Optimize a reaction path segment."),
    "path-search": (".path_search", "cli", "Search reaction pathways recursively."),
    "tsopt": (".tsopt", "cli", "Optimize a transition-state candidate."),
    "freq": (".freq", "cli", "Run vibrational analysis and thermochemistry."),
    "irc": (".irc", "cli", "Run IRC integration from a TS geometry."),
    "trj2fig": (".trj2fig", "cli", "Plot energy profile from trajectory."),
    "add-elem-info": (".add_elem_info", "cli", "Repair/add PDB element columns."),
    "dft": (".dft", "cli", "Run single-point DFT."),
    "init": (".init", "cli", "Generate starter YAML templates."),
    "scan2d": (".scan2d", "cli", "Run 2D distance scan."),
    "scan3d": (".scan3d", "cli", "Run 3D distance scan."),
    "oniom-gaussian": (".oniom_export", "cli_gaussian", "Export Gaussian ONIOM input."),
    "oniom-orca": (".oniom_export", "cli_orca", "Export ORCA QM/MM input."),
    "define-layer": (".define_layer", "cli", "Assign ML/MM layers to a structure."),
    "fix-altloc": (".fix_altloc", "cli", "Resolve PDB alternate locations."),
    "energy-diagram": (".energy_diagram", "cli", "Draw energy diagrams from values."),
}

_COMMAND_BOOL_VALUE_OPTIONS: dict[str, frozenset[str]] = {
    "extract": frozenset({"--include-H2O", "--exclude-backbone", "--add-linkH", "--verbose"}),
    "all": frozenset(
        {
            "--include-H2O",
            "--include-h2o",
            "--exclude-backbone",
            "--add-linkH",
            "--verbose",
            "--climb",
            "--dump",
            "--pre-opt",
            "--tsopt",
            "--thermo",
            "--dft",
            "--scan-one-based",
            "--scan-preopt",
            "--scan-endopt",
        }
    ),
    "scan": frozenset({"--one-based", "--dump", "--preopt", "--endopt"}),
    "scan2d": frozenset({"--one-based", "--dump", "--preopt"}),
    "scan3d": frozenset({"--one-based", "--dump", "--preopt"}),
    "mm-parm": frozenset({"--add-ter", "--add_TER", "--add-h", "--add_H"}),
    "fix-altloc": frozenset({"--recursive", "--inplace", "--overwrite", "--force"}),
}

_COMMAND_BOOL_TOGGLE_OPTIONS: dict[str, frozenset[str]] = {
    "opt": frozenset(
        {
            "--model-indices-one-based",
            "--detect-layer",
            "--one-based",
            "--layer-opt",
            "--dump",
            "--show-config",
            "--dry-run",
        }
    ),
    "dft": frozenset(
        {
            "--model-indices-one-based",
            "--detect-layer",
            "--show-config",
            "--dry-run",
        }
    ),
    "tsopt": frozenset(
        {
            "--model-indices-one-based",
            "--detect-layer",
            "--dump",
            "--partial-hessian-flatten",
            "--show-config",
            "--dry-run",
        }
    ),
    "path-opt": frozenset(
        {
            "--model-indices-one-based",
            "--detect-layer",
            "--climb",
            "--preopt",
            "--dump",
            "--show-config",
            "--dry-run",
        }
    ),
    "path-search": frozenset(
        {
            "--model-indices-one-based",
            "--detect-layer",
            "--climb",
            "--dump",
            "--pre-opt",
            "--align",
            "--show-config",
            "--dry-run",
        }
    ),
    "freq": frozenset(
        {
            "--model-indices-one-based",
            "--detect-layer",
            "--dump",
            "--show-config",
            "--dry-run",
        }
    ),
    "irc": frozenset(
        {
            "--model-indices-one-based",
            "--detect-layer",
            "--forward",
            "--backward",
            "--show-config",
            "--dry-run",
        }
    ),
}

_COMMAND_BOOL_TOGGLE_NEGATIVE_ALIASES: dict[str, dict[str, str]] = {
    "opt": {
        "--model-indices-one-based": "--model-indices-zero-based",
        "--one-based": "--zero-based",
    },
    "dft": {
        "--model-indices-one-based": "--model-indices-zero-based",
    },
    "tsopt": {
        "--model-indices-one-based": "--model-indices-zero-based",
        "--partial-hessian-flatten": "--full-hessian-flatten",
    },
    "path-opt": {
        "--model-indices-one-based": "--model-indices-zero-based",
    },
    "path-search": {
        "--model-indices-one-based": "--model-indices-zero-based",
    },
    "freq": {
        "--model-indices-one-based": "--model-indices-zero-based",
    },
    "irc": {
        "--model-indices-one-based": "--model-indices-zero-based",
    },
}

_SUBCOMMAND_PRIMARY_HELP_OPTIONS: dict[str, frozenset[str]] = {
    "scan": frozenset(
        {
            "-i",
            "--input",
            "--real-parm7",
            "--model-pdb",
            "--detect-layer",
            "-q",
            "--charge",
            "-m",
            "--multiplicity",
            "--spec",
            "--out-dir",
            "--help-advanced",
        }
    ),
    "scan2d": frozenset(
        {
            "-i",
            "--input",
            "--real-parm7",
            "--model-pdb",
            "--detect-layer",
            "-q",
            "--charge",
            "-m",
            "--multiplicity",
            "--spec",
            "--out-dir",
            "--help-advanced",
        }
    ),
    "scan3d": frozenset(
        {
            "-i",
            "--input",
            "--real-parm7",
            "--model-pdb",
            "--detect-layer",
            "-q",
            "--charge",
            "-m",
            "--multiplicity",
            "--spec",
            "--out-dir",
            "--help-advanced",
        }
    ),
    "opt": frozenset(
        {
            "-i",
            "--input",
            "--real-parm7",
            "--model-pdb",
            "--detect-layer",
            "-q",
            "--charge",
            "-m",
            "--multiplicity",
            "--config",
            "--dry-run",
            "--out-dir",
            "--opt-mode",
            "--help-advanced",
        }
    ),
    "path-opt": frozenset(
        {
            "-i",
            "--input",
            "--real-parm7",
            "--model-pdb",
            "--detect-layer",
            "-q",
            "--charge",
            "-m",
            "--multiplicity",
            "--mep-mode",
            "--config",
            "--dry-run",
            "--max-nodes",
            "--out-dir",
            "--help-advanced",
        }
    ),
    "path-search": frozenset(
        {
            "-i",
            "--input",
            "--real-parm7",
            "--model-pdb",
            "--detect-layer",
            "-q",
            "--charge",
            "-m",
            "--multiplicity",
            "--mep-mode",
            "--config",
            "--dry-run",
            "--max-nodes",
            "--out-dir",
            "--help-advanced",
        }
    ),
    "tsopt": frozenset(
        {
            "-i",
            "--input",
            "--real-parm7",
            "--model-pdb",
            "--detect-layer",
            "-q",
            "--charge",
            "-m",
            "--multiplicity",
            "--config",
            "--dry-run",
            "--max-cycles",
            "--opt-mode",
            "--out-dir",
            "--help-advanced",
        }
    ),
    "freq": frozenset(
        {
            "-i",
            "--input",
            "--real-parm7",
            "--model-pdb",
            "--detect-layer",
            "-q",
            "--charge",
            "-m",
            "--multiplicity",
            "--temperature",
            "--pressure",
            "--out-dir",
            "--help-advanced",
        }
    ),
    "irc": frozenset(
        {
            "-i",
            "--input",
            "--real-parm7",
            "--model-pdb",
            "--detect-layer",
            "-q",
            "--charge",
            "-m",
            "--multiplicity",
            "--max-cycles",
            "--step-size",
            "--forward",
            "--backward",
            "--out-dir",
            "--help-advanced",
        }
    ),
    "dft": frozenset(
        {
            "-i",
            "--input",
            "--real-parm7",
            "--model-pdb",
            "--detect-layer",
            "-q",
            "--charge",
            "-m",
            "--multiplicity",
            "--func-basis",
            "--config",
            "--dry-run",
            "--out-dir",
            "--help-advanced",
        }
    ),
    "mm-parm": frozenset(
        {
            "-i",
            "--input",
            "--out-prefix",
            "--ligand-charge",
            "--ligand-mult",
            "--ff-set",
            "--allow-nonstandard-aa",
            "--help-advanced",
        }
    ),
    "define-layer": frozenset(
        {
            "-i",
            "--input",
            "--model-pdb",
            "--model-indices",
            "-o",
            "--output",
            "--help-advanced",
        }
    ),
    "oniom-gaussian": frozenset(
        {
            "--parm7",
            "-i",
            "--input",
            "--model-pdb",
            "-o",
            "--output",
            "--method",
            "-q",
            "--charge",
            "-m",
            "--mult",
            "--help-advanced",
        }
    ),
    "oniom-orca": frozenset(
        {
            "--parm7",
            "-i",
            "--input",
            "--model-pdb",
            "-o",
            "--output",
            "--method",
            "-q",
            "--charge",
            "-m",
            "--mult",
            "--total-charge",
            "--total-mult",
            "--help-advanced",
        }
    ),
    "add-elem-info": frozenset(
        {
            "-i",
            "--input",
            "-o",
            "--out",
            "--help-advanced",
        }
    ),
    "trj2fig": frozenset(
        {
            "-i",
            "--input",
            "-o",
            "--out",
            "--unit",
            "-r",
            "--reference",
            "--help-advanced",
        }
    ),
    "energy-diagram": frozenset(
        {
            "-i",
            "--input",
            "-o",
            "--output",
            "--help-advanced",
        }
    ),
}

_PARSER_WRAPPER_SUBCOMMANDS = frozenset({"extract", "fix-altloc"})


def _parse_bool_literal(raw: str) -> bool | None:
    token = raw.strip().lower()
    if token in _BOOL_TRUE_LITERALS:
        return True
    if token in _BOOL_FALSE_LITERALS:
        return False
    return None


def _toggle_negative_name(command: str, positive_name: str) -> str:
    aliases = _COMMAND_BOOL_TOGGLE_NEGATIVE_ALIASES.get(command)
    if aliases and positive_name in aliases:
        return aliases[positive_name]
    return f"--no-{positive_name[2:]}"


def _normalize_bool_argv(args: list[str]) -> tuple[list[str], bool]:
    if not args:
        return args, False

    command = args[0]
    bool_value_options = _COMMAND_BOOL_VALUE_OPTIONS.get(command)
    bool_toggle_options = _COMMAND_BOOL_TOGGLE_OPTIONS.get(command)
    if not bool_value_options and not bool_toggle_options:
        return args, False

    normalized: list[str] = [command]
    legacy_used = False
    i = 1
    while i < len(args):
        token = args[i]

        if token == "--":
            normalized.extend(args[i:])
            break

        if not token.startswith("--"):
            normalized.append(token)
            i += 1
            continue

        name, sep, inline_value = token.partition("=")
        if name.startswith("--no-"):
            positive_name = "--" + name[5:]
            if sep == "" and bool_value_options and positive_name in bool_value_options:
                normalized.extend([positive_name, "False"])
                i += 1
                continue
            if bool_toggle_options and positive_name in bool_toggle_options:
                if sep:
                    parsed_inline = _parse_bool_literal(inline_value)
                    if parsed_inline is not None:
                        legacy_used = True
                        normalized.append(name if parsed_inline else positive_name)
                    else:
                        normalized.append(token)
                    i += 1
                    continue
                if i + 1 < len(args):
                    parsed_next = _parse_bool_literal(args[i + 1])
                    if parsed_next is not None:
                        legacy_used = True
                        normalized.append(name if parsed_next else positive_name)
                        i += 2
                        continue
                normalized.append(name)
                i += 1
                continue
            normalized.append(token)
            i += 1
            continue

        if bool_value_options and name in bool_value_options:
            if sep:
                parsed_inline = _parse_bool_literal(inline_value)
                if parsed_inline is not None:
                    legacy_used = True
                    normalized.extend([name, "True" if parsed_inline else "False"])
                else:
                    normalized.append(token)
                i += 1
                continue

            if i + 1 < len(args):
                parsed_next = _parse_bool_literal(args[i + 1])
                if parsed_next is not None:
                    legacy_used = True
                    normalized.extend([name, "True" if parsed_next else "False"])
                    i += 2
                    continue

            normalized.extend([name, "True"])
            i += 1
            continue

        if bool_toggle_options and name in bool_toggle_options:
            if sep:
                parsed_inline = _parse_bool_literal(inline_value)
                if parsed_inline is not None:
                    legacy_used = True
                    normalized.append(name if parsed_inline else _toggle_negative_name(command, name))
                else:
                    normalized.append(token)
                i += 1
                continue

            if i + 1 < len(args):
                parsed_next = _parse_bool_literal(args[i + 1])
                if parsed_next is not None:
                    legacy_used = True
                    normalized.append(name if parsed_next else _toggle_negative_name(command, name))
                    i += 2
                    continue

            normalized.append(name)
            i += 1
            continue

        normalized.append(token)
        i += 1

    return normalized, legacy_used


def _build_unavailable_command(command_name: str, exc: ImportError) -> click.Command:
    missing = exc.name if isinstance(exc, ModuleNotFoundError) else None
    msg_lines = [
        f"Command '{command_name}' is unavailable because the module could not be imported."
    ]
    if missing:
        msg_lines.append(f"Missing dependency: {missing}")
        msg_lines.append("Install the missing dependency in your runtime environment and retry.")
    else:
        msg_lines.append(f"Import error: {exc}")

    help_text = (
        f"[Unavailable] {command_name} command.\n"
        "The command failed to import due to a missing dependency."
    )

    @click.command(name=command_name, help=help_text)
    def _unavailable() -> None:
        raise click.ClickException("\n".join(msg_lines))

    return _unavailable


def _show_advanced_subcommand_help(
    ctx: click.Context, _param: click.Parameter, value: bool
) -> None:
    """Print subcommand help with advanced options and exit."""
    if not value or ctx.resilient_parsing:
        return

    if getattr(ctx.command, "_advanced_passthrough_help", False):
        try:
            ctx.command.main(args=["--help"], standalone_mode=False)
        except SystemExit as exc:
            code = getattr(exc, "code", 1)
            if code not in (None, 0):
                raise
        ctx.exit()

    hidden = getattr(ctx.command, "_advanced_hidden_options", ())
    restored: list[click.Option] = []
    for opt in hidden:
        if opt.hidden:
            opt.hidden = False
            restored.append(opt)
    try:
        click.echo(ctx.command.get_help(ctx))
    finally:
        for opt in restored:
            opt.hidden = True
    ctx.exit()


def _ensure_help_advanced_option(command: click.Command) -> click.Command:
    """Attach --help-advanced to lazily loaded subcommands when absent."""
    passthrough_help = (
        command.context_settings.get("help_option_names") == []
        if isinstance(command.context_settings, dict)
        else False
    )
    setattr(command, "_advanced_passthrough_help", passthrough_help)

    if any(
        isinstance(param, click.Option) and "--help-advanced" in param.opts
        for param in command.params
    ):
        return command

    option = click.Option(
        ["--help-advanced"],
        is_flag=True,
        is_eager=True,
        expose_value=False,
        callback=_show_advanced_subcommand_help,
        help="Show all options (including advanced settings) and exit.",
    )
    command.params.insert(0, option)
    return command


def _configure_subcommand_help_visibility(
    command_name: str, command: click.Command
) -> click.Command:
    """Hide advanced options from default --help for selected subcommands."""
    if hasattr(command, "_advanced_hidden_options"):
        return command

    primary_options = _SUBCOMMAND_PRIMARY_HELP_OPTIONS.get(command_name)
    if not primary_options:
        return command

    hidden_options: list[click.Option] = []
    for param in command.params:
        if not isinstance(param, click.Option):
            continue
        names = set(param.opts + param.secondary_opts)
        if names & primary_options:
            continue
        if param.hidden:
            continue
        param.hidden = True
        hidden_options.append(param)

    setattr(command, "_advanced_hidden_options", tuple(hidden_options))
    return command


class DefaultGroup(click.Group):
    """Custom Click group that defaults to 'all' command when no subcommand provided."""

    def __init__(
        self,
        *args,
        default: str | None = None,
        lazy_subcommands: dict[str, tuple[str, str, str]] | None = None,
        **kwargs,
    ):
        super().__init__(*args, **kwargs)
        self._default_cmd = default
        self._lazy_subcommands = dict(lazy_subcommands or {})
        self._lazy_cache: dict[str, click.Command] = {}

    def parse_args(self, ctx, args):
        show_help_or_version = any(a in ("-h", "--help", "--version") for a in args)

        if self._default_cmd is not None and not show_help_or_version:
            if not args or args[0].startswith("-"):
                args.insert(0, self._default_cmd)

        args, used_legacy_bool = _normalize_bool_argv(args)
        if used_legacy_bool:
            click.echo(
                "[deprecation] Legacy bool syntax '--flag True/False' is supported for now; "
                "prefer '--flag/--no-flag'.",
                err=True,
            )
        return super().parse_args(ctx, args)

    def invoke(self, ctx):
        # Add a leading blank line for subcommands (except "all") to separate CLI tool logs.
        if ctx.invoked_subcommand and ctx.invoked_subcommand != "all":
            click.echo()
        return super().invoke(ctx)

    def list_commands(self, ctx):
        cmds = set(super().list_commands(ctx))
        cmds.update(self._lazy_subcommands.keys())
        return sorted(cmds)

    def get_command(self, ctx, cmd_name):
        cmd = super().get_command(ctx, cmd_name)
        if cmd is not None:
            if cmd_name not in _PARSER_WRAPPER_SUBCOMMANDS:
                cmd = _ensure_help_advanced_option(cmd)
            return _configure_subcommand_help_visibility(cmd_name, cmd)

        lazy_spec = self._lazy_subcommands.get(cmd_name)
        if lazy_spec is None:
            return None

        cached = self._lazy_cache.get(cmd_name)
        if cached is not None:
            return cached

        module_name, attr_name, _ = lazy_spec
        try:
            module = importlib.import_module(module_name, package=__package__)
            loaded_cmd = getattr(module, attr_name)
        except (ModuleNotFoundError, ImportError) as exc:
            loaded_cmd = _build_unavailable_command(cmd_name, exc)

        if cmd_name not in _PARSER_WRAPPER_SUBCOMMANDS:
            loaded_cmd = _ensure_help_advanced_option(loaded_cmd)
        loaded_cmd = _configure_subcommand_help_visibility(cmd_name, loaded_cmd)
        self._lazy_cache[cmd_name] = loaded_cmd
        return loaded_cmd

    def format_commands(self, ctx, formatter):
        rows = []
        for subcommand in self.list_commands(ctx):
            lazy_spec = self._lazy_subcommands.get(subcommand)
            if lazy_spec is not None:
                rows.append((subcommand, lazy_spec[2]))
                continue

            cmd = super().get_command(ctx, subcommand)
            if cmd is None or cmd.hidden:
                continue
            rows.append(
                (
                    subcommand,
                    cmd.get_short_help_str(formatter.width - 6 - len(subcommand)),
                )
            )

        if rows:
            with formatter.section("Commands"):
                formatter.write_dl(rows)


@click.group(
    cls=DefaultGroup,
    default="all",
    lazy_subcommands=_LAZY_SUBCOMMANDS,
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


cli.add_command(extract_cmd, name="extract")

# Disable pysisyphus logging
logging.disable(logging.CRITICAL)

# Filter noisy UMA/pydmf warnings that clutter CLI output
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
