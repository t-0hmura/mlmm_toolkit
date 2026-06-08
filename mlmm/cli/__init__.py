"""L1 Interface — Click CLI entry point + helpers.

Modules:
- ``app`` — Click group definition, ``_LAZY_SUBCOMMANDS`` registry, bool registries.
- ``bool_compat`` — legacy ``--flag True/False`` syntax normalization.
- ``common_options`` — shared Click option factories (ML region + charge/spin triple,
  MLIP backend precision, etc.) wired across multiple subcommands.
- ``decorators`` — ``_write_error_json``, ``render_cli_exception``, YAML config loaders.
- ``default_group`` — ``DefaultGroup`` (Click MultiCommand with a default subcommand).
- ``help_pages`` — ``--help`` / ``--help-advanced`` progressive disclosure.
- ``preflight`` — shared preflight (file existence, required commands, required Python modules).

The ``cli`` callable is re-exported from ``mlmm.cli.app`` so that
``from mlmm.cli import cli`` (used by ``__main__.py`` and the ``console_scripts``
entry point in ``pyproject.toml``) keeps working.
"""
from mlmm.cli.app import cli  # noqa: F401
