# mlmm/cli/app.py

from __future__ import annotations

import logging
import sys
import warnings
from pathlib import Path

import click

from mlmm.cli.help_pages import (
    _configure_subcommand_help_visibility,
    _ensure_help_advanced_option,
)
from mlmm.cli.bool_compat import normalize_bool_argv
from mlmm.cli.default_group import DefaultGroup
from mlmm import __version__

_MLMM_BANNER = r"""
 __  __ _     __  __ __  __       _____           _ _    _ _
|  \/  | |   |  \/  |  \/  |     |_   _|__   ___ | | | _|_| |_
| |\/| | |   | |\/| | |\/| | _____ | |/ _ \ / _ \| | |/ / | __|
| |  | | |___| |  | | |  | | _____ | | |_| | |_| | |   <| | |_
|_|  |_|_____|_|  |_|_|  |_|       |_|\___/ \___/|_|_|\_\_|\__|
""".strip("\n")

_CONSOLE_SCRIPT_NAMES = {"mlmm"}


def _command_argv(argv: list[str]) -> list[str]:
    if not argv:
        return []
    argv0_name = Path(argv[0]).name
    if argv0_name == "__main__.py":
        return [sys.executable, "-m", "mlmm", *argv[1:]]
    if argv0_name in _CONSOLE_SCRIPT_NAMES:
        return [argv0_name, *argv[1:]]
    return list(argv)


def _display_arg(arg: str) -> str:
    """Return a readable argv token for the startup command log."""
    if arg == "":
        return '""'
    if not any(ch.isspace() or ch in {'"', "'"} for ch in arg):
        return arg
    return '"' + arg.replace('"', r'\"') + '"'


def _quoted_argv(argv: list[str]) -> str:
    """Return a readable representation of the executed argv.

    The shell removes the user's original quote characters before Python
    starts, so the literal typed string is not recoverable from ``sys.argv``.
    This keeps argv values intact while avoiding shlex's hard-to-read
    ``'"'"'`` single-quote escaping in logs.
    """
    return " ".join(_display_arg(str(arg)) for arg in _command_argv(argv))


def _has_help_or_version_request(argv: list[str]) -> bool:
    return any(arg in {"-h", "--help", "--version", "--help-advanced"} for arg in argv[1:])


def _emit_start_header(ctx: click.Context) -> None:
    from mlmm.core.utils import emit, is_child_mode, verbose_level

    if is_child_mode() or _has_help_or_version_request(sys.argv):
        return

    if verbose_level() >= 2:
        emit(f"{_MLMM_BANNER}\n\nmlmm-toolkit ver. {__version__}\n", narrative=True)
    else:
        emit(f"mlmm-toolkit ver. {__version__}\n", narrative=True)

    subcommand = (
        getattr(ctx, "invoked_subcommand", None)
        or (ctx.command.name if ctx.command is not None else None)
        or "all"
    )
    if verbose_level() >= 2:
        emit(f"[command] {_quoted_argv(sys.argv)}", narrative=True, raw_path=True)
        if subcommand != "all":
            emit(f"[mode] {subcommand}", narrative=True)
    emit("", narrative=True)


_LAZY_SUBCOMMANDS: dict[str, tuple[str, str, str]] = {
    "all": ("mlmm.workflows.all", "cli", "End-to-end workflow (extract -> MEP [-> TS -> IRC -> freq -> DFT])."),
    "mm-parm": ("mlmm.workflows.mm_parm", "cli", "Generate Amber parm7/rst7 topology files."),
    "scan": ("mlmm.workflows.scan", "cli", "Run staged 1D scan with harmonic restraints."),
    "opt": ("mlmm.workflows.opt", "cli", "Optimize one structure."),
    "path-opt": ("mlmm.workflows.path_opt", "cli", "Optimize a reaction path segment."),
    "path-search": ("mlmm.workflows.path_search", "cli", "Search reaction pathways recursively."),
    "tsopt": ("mlmm.workflows.tsopt", "cli", "Optimize a transition-state candidate."),
    "freq": ("mlmm.workflows.freq", "cli", "Run vibrational analysis and thermochemistry."),
    "irc": ("mlmm.workflows.irc", "cli", "Run IRC integration from a TS geometry."),
    "trj2fig": ("mlmm.io.trj2fig", "cli", "Plot energy profile from trajectory."),
    "add-elem-info": ("mlmm.domain.add_elem_info", "cli", "Repair/add PDB element columns."),
    "dft": ("mlmm.workflows.dft", "cli", "Run single-point DFT."),
    "sp": ("mlmm.workflows.sp", "cli", "Run single-point ML/MM ONIOM energy + forces."),
    "scan2d": ("mlmm.workflows.scan2d", "cli", "Run 2D distance scan."),
    "scan3d": ("mlmm.workflows.scan3d", "cli", "Run 3D distance scan."),
    "oniom-export": ("mlmm.workflows.oniom_export", "cli", "Export ONIOM input (Gaussian g16 or ORCA)."),
    "oniom-import": ("mlmm.workflows.oniom_import", "cli", "Import ONIOM input and reconstruct XYZ/layered PDB."),
    "define-layer": ("mlmm.workflows.define_layer", "cli", "Assign ML/MM layers to a structure."),
    "fix-altloc": ("mlmm.io.pdb_fix", "cli", "Resolve PDB alternate locations."),
    "energy-diagram": ("mlmm.io.energy_diagram", "cli", "Draw energy diagrams from values."),
    "extract": ("mlmm.workflows.extract", "cli", "Extract a binding pocket."),
    "bond-summary": ("mlmm.domain.bond_summary", "cli", "Detect bond changes between structures."),
}

# Only the ``all`` subcommand is listed here because it uses Click's
# ``type=click.BOOL`` (value-style) booleans that cannot be auto-detected
# from ``is_bool_flag``.  For all other subcommands the ``DefaultGroup``
# in ``default_group.py`` inspects the Click command's parameters at
# runtime and auto-discovers ``is_bool_flag`` / ``BoolParamType`` options,
# so they do not need to be repeated in these manual registries.
_COMMAND_BOOL_VALUE_OPTIONS: dict[str, frozenset[str]] = {
    "all": frozenset(
        {
            "--add-linkh",
            "--climb",
            "--dft",
            "--dump",
            "--exclude-backbone",
            "--include-h2o",
            "--preopt",
            "--scan-endopt",
            "--scan-one-based",
            "--scan-preopt",
            "--thermo",
            "--tsopt",
        }
    ),
}

# Manual toggle-option hints.  ``DefaultGroup._resolve_bool_options()``
# auto-detects toggle options from Click's ``is_bool_flag`` attribute,
# but entries here ensure correct normalization *before* the lazy
# subcommand is imported (needed for early argv rewriting).
_COMMAND_BOOL_TOGGLE_OPTIONS: dict[str, frozenset[str]] = {
    "add-elem-info": frozenset(
        {
            "--overwrite",
        }
    ),
    "all": frozenset(
        {
            "--auto-mm-add-ter",
            "--cmap",
            "--convert-files",
            "--detect-layer",
            "--dry-run",
            "--embedcharge",
            "--flatten",
            "--refine-path",
            "--show-config",
            "--skip-final-freq",
        }
    ),
    "bond-summary": frozenset(
        {
            "--json",
            "--one-based",
        }
    ),
    "define-layer": frozenset(
        {
            "--one-based",
        }
    ),
    "dft": frozenset(
        {
            "--cmap",
            "--convert-files",
            "--detect-layer",
            "--dry-run",
            "--embedcharge",
            "--lowmem",
            "--model-indices-one-based",
            "--out-json",
            "--show-config",
        }
    ),
    "energy-diagram": frozenset(
        {
            "--out-json",
        }
    ),
    "extract": frozenset(
        {
            "--add-linkh",
            "--exclude-backbone",
            "--include-h2o",
            "--out-json",
        }
    ),
    "fix-altloc": frozenset(
        {
            "--force",
            "--inplace",
            "--overwrite",
            "--recursive",
        }
    ),
    "freq": frozenset(
        {
            "--cmap",
            "--convert-files",
            "--detect-layer",
            "--dry-run",
            "--dump",
            "--embedcharge",
            "--model-indices-one-based",
            "--out-json",
            "--show-config",
        }
    ),
    "irc": frozenset(
        {
            "--backward",
            "--cmap",
            "--convert-files",
            "--detect-layer",
            "--dry-run",
            "--embedcharge",
            "--forward",
            "--model-indices-one-based",
            "--out-json",
            "--show-config",
        }
    ),
    "mm-parm": frozenset(
        {
            "--add-h",
            "--add-ter",
            "--keep-temp",
        }
    ),
    "oniom-export": frozenset(
        {
            "--convert-orcaff",
            "--element-check",
        }
    ),
    "opt": frozenset(
        {
            "--cmap",
            "--convert-files",
            "--detect-layer",
            "--dry-run",
            "--dump",
            "--embedcharge",
            "--flatten",
            "--microiter",
            "--mm-only",
            "--model-indices-one-based",
            "--one-based",
            "--out-json",
            "--show-config",
        }
    ),
    "path-opt": frozenset(
        {
            "--climb",
            "--cmap",
            "--convert-files",
            "--detect-layer",
            "--dry-run",
            "--dump",
            "--embedcharge",
            "--fix-ends",
            "--model-indices-one-based",
            "--out-json",
            "--preopt",
            "--show-config",
        }
    ),
    "path-search": frozenset(
        {
            "--align",
            "--climb",
            "--cmap",
            "--convert-files",
            "--detect-layer",
            "--dry-run",
            "--dump",
            "--embedcharge",
            "--model-indices-one-based",
            "--preopt",
            "--show-config",
        }
    ),
    "scan": frozenset(
        {
            "--cmap",
            "--convert-files",
            "--detect-layer",
            "--dry-run",
            "--dump",
            "--embedcharge",
            "--endopt",
            "--model-indices-one-based",
            "--one-based",
            "--out-json",
            "--preopt",
            "--print-parsed",
        }
    ),
    "scan2d": frozenset(
        {
            "--cmap",
            "--convert-files",
            "--detect-layer",
            "--dump",
            "--embedcharge",
            "--model-indices-one-based",
            "--one-based",
            "--out-json",
            "--preopt",
            "--print-parsed",
        }
    ),
    "scan3d": frozenset(
        {
            "--cmap",
            "--convert-files",
            "--detect-layer",
            "--dump",
            "--embedcharge",
            "--model-indices-one-based",
            "--one-based",
            "--out-json",
            "--preopt",
            "--print-parsed",
        }
    ),
    "trj2fig": frozenset(
        {
            "--reverse-x",
        }
    ),
    "tsopt": frozenset(
        {
            "--cmap",
            "--convert-files",
            "--detect-layer",
            "--dry-run",
            "--dump",
            "--embedcharge",
            "--flatten",
            "--microiter",
            "--ml-only-hessian-dimer",
            "--model-indices-one-based",
            "--out-json",
            "--partial-hessian-flatten",
            "--show-config",
            "--skip-final-freq",
        }
    ),
}

_COMMAND_BOOL_SINGLE_FLAG_OPTIONS: dict[str, frozenset[str]] = {
    "all": frozenset(
        {
            "--auto-mm-keep-temp",
        }
    ),
}

_COMMAND_BOOL_TOGGLE_NEGATIVE_ALIASES: dict[str, dict[str, str]] = {
    "scan": {
        "--one-based": "--zero-based",
    },
    "scan2d": {
        "--one-based": "--zero-based",
    },
    "scan3d": {
        "--one-based": "--zero-based",
    },
    "opt": {
        "--model-indices-one-based": "--model-indices-zero-based",
        "--one-based": "--zero-based",
        "--microiter": "--no-microiter",
    },
    "dft": {
        "--model-indices-one-based": "--model-indices-zero-based",
    },
    "tsopt": {
        "--model-indices-one-based": "--model-indices-zero-based",
        "--partial-hessian-flatten": "--full-hessian-flatten",
        "--ml-only-hessian-dimer": "--no-ml-only-hessian-dimer",
        "--microiter": "--no-microiter",
        "--convert-files": "--no-convert-files",
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
            "--parm",
            "--model-pdb",
            "--detect-layer",
            "-q",
            "--charge",
            "-l",
            "--ligand-charge",
            "-m",
            "--multiplicity",
            "-b", "--backend",
            "--embedcharge",
            "-s", "--scan-lists",
            "--config",
            "-o", "--out-dir",
            "--help-advanced",
        }
    ),
    "scan2d": frozenset(
        {
            "-i",
            "--input",
            "--parm",
            "--model-pdb",
            "--detect-layer",
            "-q",
            "--charge",
            "-l",
            "--ligand-charge",
            "-m",
            "--multiplicity",
            "-b", "--backend",
            "--embedcharge",
            "-s", "--scan-lists",
            "--config",
            "-o", "--out-dir",
            "--help-advanced",
        }
    ),
    "scan3d": frozenset(
        {
            "-i",
            "--input",
            "--parm",
            "--model-pdb",
            "--detect-layer",
            "-q",
            "--charge",
            "-l",
            "--ligand-charge",
            "-m",
            "--multiplicity",
            "-b", "--backend",
            "--embedcharge",
            "-s", "--scan-lists",
            "--csv",
            "--config",
            "-o", "--out-dir",
            "--help-advanced",
        }
    ),
    "opt": frozenset(
        {
            "-i",
            "--input",
            "--parm",
            "--model-pdb",
            "--detect-layer",
            "-q",
            "--charge",
            "-l",
            "--ligand-charge",
            "-m",
            "--multiplicity",
            "-b", "--backend",
            "--embedcharge",
            "--config",
            "-o", "--out-dir",
            "--opt-mode",
            "--max-cycles",
            "--help-advanced",
        }
    ),
    "path-opt": frozenset(
        {
            "-i",
            "--input",
            "--parm",
            "--model-pdb",
            "--detect-layer",
            "-q",
            "--charge",
            "-l",
            "--ligand-charge",
            "-m",
            "--multiplicity",
            "-b", "--backend",
            "--embedcharge",
            "--mep-mode",
            "--fix-ends",
            "--config",
            "--max-nodes",
            "-o", "--out-dir",
            "--max-cycles",
            "--help-advanced",
        }
    ),
    "path-search": frozenset(
        {
            "-i",
            "--input",
            "--parm",
            "--model-pdb",
            "--detect-layer",
            "-q",
            "--charge",
            "-l",
            "--ligand-charge",
            "-m",
            "--multiplicity",
            "-b", "--backend",
            "--embedcharge",
            "--mep-mode",
            "--refine-mode",
            "--config",
            "--max-nodes",
            "-o", "--out-dir",
            "--max-cycles",
            "--help-advanced",
        }
    ),
    "tsopt": frozenset(
        {
            "-i",
            "--input",
            "--parm",
            "--model-pdb",
            "--detect-layer",
            "-q",
            "--charge",
            "-l",
            "--ligand-charge",
            "-m",
            "--multiplicity",
            "-b", "--backend",
            "--embedcharge",
            "--config",
            "--max-cycles",
            "--opt-mode",
            "-o", "--out-dir",
            "--hessian-calc-mode",
            "--help-advanced",
        }
    ),
    "freq": frozenset(
        {
            "-i",
            "--input",
            "--parm",
            "--model-pdb",
            "--detect-layer",
            "-q",
            "--charge",
            "-l",
            "--ligand-charge",
            "-m",
            "--multiplicity",
            "-b", "--backend",
            "--embedcharge",
            "--temperature",
            "--pressure",
            "--config",
            "-o", "--out-dir",
            "--hessian-calc-mode",
            "--help-advanced",
        }
    ),
    "irc": frozenset(
        {
            "-i",
            "--input",
            "--parm",
            "--model-pdb",
            "--detect-layer",
            "-q",
            "--charge",
            "-l",
            "--ligand-charge",
            "-m",
            "--multiplicity",
            "-b", "--backend",
            "--embedcharge",
            "--max-cycles",
            "--step-size",
            "--forward",
            "--backward",
            "--config",
            "-o", "--out-dir",
            "--hessian-calc-mode",
            "--help-advanced",
        }
    ),
    "dft": frozenset(
        {
            "-i",
            "--input",
            "--parm",
            "--model-pdb",
            "--detect-layer",
            "-q",
            "--charge",
            "-l",
            "--ligand-charge",
            "-m",
            "--multiplicity",
            "-b", "--backend",
            "--embedcharge",
            "--func-basis",
            "--config",
            "-o", "--out-dir",
            "--help-advanced",
        }
    ),
    "sp": frozenset(
        {
            "-i",
            "--input",
            "--parm",
            "--model-pdb",
            "--detect-layer",
            "-q",
            "--charge",
            "-l",
            "--ligand-charge",
            "-m",
            "--multiplicity",
            "-b", "--backend",
            "--embedcharge",
            "--hess",
            "--hessian-calc-mode",
            "--config",
            "-o", "--out-dir",
            "--help-advanced",
        }
    ),
    "mm-parm": frozenset(
        {
            "-i",
            "--input",
            "--out-prefix",
            "-l",
            "--ligand-charge",
            "--ligand-mult",
            "--ff-set",
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
    "oniom-export": frozenset(
        {
            "--parm",
            "-i",
            "--input",
            "--model-pdb",
            "-o",
            "--output",
            "--mode",
            "--method",
            "-q",
            "--charge",
            "-m",
            "--multiplicity",
            "--help-advanced",
        }
    ),
    "oniom-import": frozenset(
        {
            "-i",
            "--input",
            "--mode",
            "-o",
            "--out-prefix",
            "--ref-pdb",
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
            "-q",
            "--charge",
            "-m",
            "--multiplicity",
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
    "extract": frozenset(
        {
            "-i",
            "--input",
            "-c",
            "--center",
            "-o",
            "--output",
            "-r",
            "--radius",
            "-l",
            "--ligand-charge",
            "-v",
            "--verbose",
            "--help-advanced",
        }
    ),
    "fix-altloc": frozenset(
        {
            "-i",
            "--input",
            "-o",
            "--out",
            "--recursive",
            "--no-recursive",
            "--help-advanced",
        }
    ),
    "bond-summary": frozenset(
        {
            "-i",
            "--input",
            "--device",
            "--bond-factor",
            "--help-advanced",
        }
    ),
}

_PARSER_WRAPPER_SUBCOMMANDS = frozenset()


_PARSER_WRAPPER_BOOL_OPTION_PROVIDERS: dict[str, object] = {}

_COMMAND_GROUPS: dict[str, tuple[str, ...]] = {
    "Pipelines": ("all",),
    "Pipeline stages": (
        "opt", "sp", "tsopt", "freq", "irc", "dft",
        "scan", "scan2d", "scan3d", "path-opt", "path-search",
    ),
    "Inputs & topology": (
        "extract", "define-layer", "mm-parm",
        "oniom-export", "oniom-import", "fix-altloc",
        "add-elem-info",
    ),
    "Analysis": ("energy-diagram", "bond-summary", "trj2fig"),
}

_DEFAULT_GROUP_KWARGS = {
    "command_bool_value_options": _COMMAND_BOOL_VALUE_OPTIONS,
    "command_bool_toggle_options": _COMMAND_BOOL_TOGGLE_OPTIONS,
    "command_bool_toggle_negative_aliases": _COMMAND_BOOL_TOGGLE_NEGATIVE_ALIASES,
    "command_bool_single_flag_options": _COMMAND_BOOL_SINGLE_FLAG_OPTIONS,
    "parser_wrapper_subcommands": _PARSER_WRAPPER_SUBCOMMANDS,
    "parser_wrapper_bool_option_providers": _PARSER_WRAPPER_BOOL_OPTION_PROVIDERS,
    "subcommand_primary_help_options": _SUBCOMMAND_PRIMARY_HELP_OPTIONS,
    "normalize_bool_argv": normalize_bool_argv,
    "ensure_help_advanced_option": _ensure_help_advanced_option,
    "configure_subcommand_help_visibility": _configure_subcommand_help_visibility,
    "command_groups": _COMMAND_GROUPS,
}


def _verbose_callback(ctx: click.Context, param: click.Parameter, value: int) -> int:
    # `-v/--verbose LEVEL` (0-3) is the single, unified verbosity control:
    #   0 silent · 1 milestones · 2 (default) +cycle tables/timing/VRAM/paths
    #   · 3 everything (config dumps, per-file paths, DEBUG logging).
    # Injected into every subcommand (see `_ensure_verbose_option`), so this
    # eager callback is where each invocation sets its level, enables the
    # console gate, and emits the start header once the level is known.
    from mlmm.core.utils import set_verbose_level, set_console_gating
    set_verbose_level(value)
    set_console_gating(True)
    if value >= 3:
        from mlmm.core.logging import setup_logging
        setup_logging(value)
    _emit_start_header(ctx)
    return value


# Subcommands that run an iterative optimizer, so their -v 2 genuinely adds
# per-cycle optimizer tables and GPU VRAM accounting. Every other subcommand
# is a single-shot compute (sp/freq/dft) or a lightweight IO/analysis step:
# its -v 2 adds detailed step logging + deliverable paths but no optimizer
# cycle tables or VRAM lines, so claiming them in the help would mislead.
_OPTIMIZER_SUBCOMMANDS = frozenset({
    "all", "opt", "tsopt", "irc",
    "scan", "scan2d", "scan3d",
    "path-opt", "path-search",
})


def _verbose_help(name: str | None) -> str:
    """`-v/--verbose` help text, specialised per subcommand so the level-2
    description never claims optimizer cycle tables / VRAM for the non-optimizer
    commands (extract, sp, freq, dft, file IO, analysis)."""
    if name in _OPTIMIZER_SUBCOMMANDS:
        level2 = "optimizer cycle tables, per-stage timing, VRAM, deliverable paths"
    else:
        level2 = "detailed step logging and deliverable paths"
    return (
        "Console verbosity 0-3 (default 2). 0=silent; 1=milestones only; "
        f"2=+{level2}; "
        "3=everything (full config blocks, per-file paths, DEBUG logging)."
    )


def _ensure_verbose_option(
    command: click.Command, cmd_name: str | None = None
) -> click.Command:
    """Attach the unified `-v/--verbose LEVEL` control to a subcommand when absent.

    Verbosity is a per-subcommand option (not a root-group option) so it can be
    written in the natural position, e.g. ``mlmm opt -v 0 ...``.  The eager
    callback sets the level, enables the console gate, and emits the start header.
    ``cmd_name`` is the registry key (``command.name`` is unreliable: lazy
    subcommands share the underlying ``cli`` function name) and selects the
    per-subcommand help text.
    """
    if any(
        isinstance(param, click.Option)
        and ("--verbose" in param.opts or "-v" in param.opts)
        for param in command.params
    ):
        return command

    option = click.Option(
        ["-v", "--verbose"],
        type=click.IntRange(0, 3),
        default=2,
        metavar="LEVEL",
        is_eager=True,
        expose_value=False,
        callback=_verbose_callback,
        help=_verbose_help(cmd_name),
    )
    command.params.insert(0, option)
    return command


@click.group(
    cls=DefaultGroup,
    default="all",
    lazy_subcommands=_LAZY_SUBCOMMANDS,
    **_DEFAULT_GROUP_KWARGS,
    ensure_verbose_option=_ensure_verbose_option,
    help="mlmm: Run workflow steps via subcommands.",
    context_settings={"help_option_names": ["-h", "--help"]},
)
@click.version_option(version=__version__, prog_name="mlmm")
def cli() -> None:
    # `-v/--verbose` is a per-subcommand option (injected by
    # `_ensure_verbose_option`); the start header is emitted by that option's
    # eager callback once the level is known, not from this group body.
    from mlmm.core.utils import set_base_dir
    set_base_dir(Path.cwd())


# Pysisyphus log suppression is handled by DefaultGroup._silence_pysisyphus_loggers()
# which runs after lazy subcommand import (when pysisyphus __init__ handlers are created).

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
warnings.filterwarnings(
    "ignore",
    category=UserWarning,
    message=r"Setting global torch default dtype to torch\.float32\.",
)
# Suppress fairchem dataset_list deprecation warning (baked into UMA checkpoint
# config; cannot be fixed caller-side until fairchem removes dataset_list from
# checkpoints). The message is emitted from `escn_md.resolve_dataset_mapping`
# via the bare `logging.warning(...)` call, which routes through the root
# logger (visible as `WARNING:root:...`), not a `fairchem.*` named logger.
# Filter the root logger but require the verbatim opening clause so unrelated
# user log messages that happen to contain "dataset_list" survive.
_DATASET_LIST_MSG = "If 'dataset_list' is provided in the config"
logging.getLogger().addFilter(
    lambda record: _DATASET_LIST_MSG not in record.getMessage()
)

if __name__ == "__main__":
    cli()
