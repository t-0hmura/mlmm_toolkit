# mlmm_toolkit/cli.py

import logging
import warnings

import click

from .advanced_help import (
    _configure_subcommand_help_visibility,
    _ensure_help_advanced_option,
)
from .bool_compat import normalize_bool_argv
from .default_group import DefaultGroup
from mlmm_toolkit import __version__

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
    "scan2d": (".scan2d", "cli", "Run 2D distance scan."),
    "scan3d": (".scan3d", "cli", "Run 3D distance scan."),
    "oniom-export": (".oniom_export", "cli", "Export ONIOM input (Gaussian g16 or ORCA)."),
    "oniom-import": (".oniom_import", "cli", "Import ONIOM input and reconstruct XYZ/layered PDB."),
    "define-layer": (".define_layer", "cli", "Assign ML/MM layers to a structure."),
    "fix-altloc": (".fix_altloc", "cli", "Resolve PDB alternate locations."),
    "energy-diagram": (".energy_diagram", "cli", "Draw energy diagrams from values."),
    "extract": (".extract", "cli", "Extract a binding pocket."),
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
            "--include-H2O",
            "--include-h2o",
            "--exclude-backbone",
            "--add-linkH",
            "--verbose",
            "--climb",
            "--dump",
            "--preopt",
            "--tsopt",
            "--thermo",
            "--dft",
            "--scan-one-based",
            "--scan-preopt",
            "--scan-endopt",
        }
    ),
}

# Manual toggle-option hints.  ``DefaultGroup._resolve_bool_options()``
# auto-detects toggle options from Click's ``is_bool_flag`` attribute,
# but entries here ensure correct normalization *before* the lazy
# subcommand is imported (needed for early argv rewriting).
_COMMAND_BOOL_TOGGLE_OPTIONS: dict[str, frozenset[str]] = {
    "all": frozenset(
        {
            "--show-config",
            "--dry-run",
            "--detect-layer",
            "--auto-mm-add-ter",
            "--convert-files",
            "--embedcharge",
        }
    ),
    "mm-parm": frozenset({"--add-ter", "--add-h", "--allow-nonstandard-aa", "--keep-temp"}),
    "scan": frozenset({"--one-based", "--dump", "--preopt", "--endopt", "--convert-files", "--embedcharge", "--detect-layer", "--model-indices-one-based", "--print-parsed"}),
    "scan2d": frozenset({"--one-based", "--dump", "--preopt", "--convert-files", "--embedcharge", "--detect-layer", "--model-indices-one-based", "--print-parsed"}),
    "scan3d": frozenset({"--one-based", "--dump", "--preopt", "--convert-files", "--embedcharge", "--detect-layer", "--model-indices-one-based", "--print-parsed"}),
    "opt": frozenset(
        {
            "--model-indices-one-based",
            "--detect-layer",
            "--one-based",
            "--dump",
            "--microiter",
            "--flatten",
            "--convert-files",
            "--show-config",
            "--dry-run",
            "--embedcharge",
        }
    ),
    "dft": frozenset(
        {
            "--model-indices-one-based",
            "--detect-layer",
            "--convert-files",
            "--show-config",
            "--dry-run",
            "--embedcharge",
        }
    ),
    "tsopt": frozenset(
        {
            "--model-indices-one-based",
            "--detect-layer",
            "--dump",
            "--microiter",
            "--partial-hessian-flatten",
            "--ml-only-hessian-dimer",
            "--show-config",
            "--dry-run",
            "--convert-files",
            "--embedcharge",
        }
    ),
    "path-opt": frozenset(
        {
            "--model-indices-one-based",
            "--detect-layer",
            "--climb",
            "--fix-ends",
            "--preopt",
            "--dump",
            "--convert-files",
            "--show-config",
            "--dry-run",
            "--embedcharge",
        }
    ),
    "path-search": frozenset(
        {
            "--model-indices-one-based",
            "--detect-layer",
            "--climb",
            "--dump",
            "--preopt",
            "--align",
            "--convert-files",
            "--show-config",
            "--dry-run",
            "--embedcharge",
        }
    ),
    "freq": frozenset(
        {
            "--model-indices-one-based",
            "--detect-layer",
            "--dump",
            "--convert-files",
            "--show-config",
            "--dry-run",
            "--embedcharge",
        }
    ),
    "irc": frozenset(
        {
            "--model-indices-one-based",
            "--detect-layer",
            "--forward",
            "--backward",
            "--convert-files",
            "--show-config",
            "--dry-run",
            "--embedcharge",
        }
    ),
    "oniom-export": frozenset({"--element-check", "--convert-orcaff"}),
    "trj2fig": frozenset({"--reverse-x"}),
    "add-elem-info": frozenset({"--overwrite"}),
    "extract": frozenset({"--include-H2O", "--exclude-backbone", "--add-linkH", "--verbose"}),
    "fix-altloc": frozenset({"--recursive", "--inplace", "--overwrite", "--force"}),
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
            "-m",
            "--multiplicity",
            "--backend",
            "--spec",
            "--out-dir",
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
            "-m",
            "--multiplicity",
            "--backend",
            "--spec",
            "--out-dir",
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
            "-m",
            "--multiplicity",
            "--backend",
            "--csv",
            "--spec",
            "--out-dir",
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
            "-m",
            "--multiplicity",
            "--backend",
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
            "--parm",
            "--model-pdb",
            "--detect-layer",
            "-q",
            "--charge",
            "-m",
            "--multiplicity",
            "--backend",
            "--mep-mode",
            "--fix-ends",
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
            "--parm",
            "--model-pdb",
            "--detect-layer",
            "-q",
            "--charge",
            "-m",
            "--multiplicity",
            "--backend",
            "--mep-mode",
            "--refine-mode",
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
            "--parm",
            "--model-pdb",
            "--detect-layer",
            "-q",
            "--charge",
            "-m",
            "--multiplicity",
            "--backend",
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
            "--parm",
            "--model-pdb",
            "--detect-layer",
            "-q",
            "--charge",
            "-m",
            "--multiplicity",
            "--backend",
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
            "--parm",
            "--model-pdb",
            "--detect-layer",
            "-q",
            "--charge",
            "-m",
            "--multiplicity",
            "--backend",
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
            "--parm",
            "--model-pdb",
            "--detect-layer",
            "-q",
            "--charge",
            "-m",
            "--multiplicity",
            "--backend",
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
            "--ligand-charge",
            "--help-advanced",
        }
    ),
}

_PARSER_WRAPPER_SUBCOMMANDS = frozenset()


_PARSER_WRAPPER_BOOL_OPTION_PROVIDERS: dict[str, object] = {}

_DEFAULT_GROUP_KWARGS = {
    "command_bool_value_options": _COMMAND_BOOL_VALUE_OPTIONS,
    "command_bool_toggle_options": _COMMAND_BOOL_TOGGLE_OPTIONS,
    "command_bool_toggle_negative_aliases": _COMMAND_BOOL_TOGGLE_NEGATIVE_ALIASES,
    "parser_wrapper_subcommands": _PARSER_WRAPPER_SUBCOMMANDS,
    "parser_wrapper_bool_option_providers": _PARSER_WRAPPER_BOOL_OPTION_PROVIDERS,
    "subcommand_primary_help_options": _SUBCOMMAND_PRIMARY_HELP_OPTIONS,
    "normalize_bool_argv": normalize_bool_argv,
    "ensure_help_advanced_option": _ensure_help_advanced_option,
    "configure_subcommand_help_visibility": _configure_subcommand_help_visibility,
}


@click.group(
    cls=DefaultGroup,
    default="all",
    lazy_subcommands=_LAZY_SUBCOMMANDS,
    **_DEFAULT_GROUP_KWARGS,
    help="mlmm: Run workflow steps via subcommands.",
    context_settings={"help_option_names": ["-h", "--help"]},
)
@click.version_option(version=__version__, prog_name="mlmm")
def cli() -> None:
    click.echo(f"mlmm ver. {__version__}\n")


# Silence pysisyphus logger without muting application/global logging.
logging.getLogger("pysisyphus").setLevel(logging.CRITICAL)
logging.getLogger("pysisyphus").propagate = False

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
