"""Shared Click group helpers for lazy subcommand loading and bool normalization."""

from __future__ import annotations

import importlib
from collections.abc import Callable, Mapping

import click

from .bool_compat import normalize_argv_option_names as _normalize_argv_option_names

LazySubcommands = Mapping[str, tuple[str, str, str]]
BoolOptionsByCommand = Mapping[str, frozenset[str]]
BoolNegativeAliasesByCommand = Mapping[str, Mapping[str, str]]
BoolSingleFlagOptionsByCommand = Mapping[str, frozenset[str]]
PrimaryHelpOptionsByCommand = Mapping[str, frozenset[str]]
ParserWrapperBoolOptionProviders = Mapping[str, Callable[[], frozenset[str]]]

NormalizeBoolArgvFunc = Callable[
    [
        list[str],
        BoolOptionsByCommand,
        BoolOptionsByCommand,
        BoolNegativeAliasesByCommand,
        BoolSingleFlagOptionsByCommand,
    ],
    tuple[list[str], bool],
]
EnsureHelpAdvancedOptionFunc = Callable[[click.Command], click.Command]
ConfigureSubcommandHelpVisibilityFunc = Callable[
    [str, click.Command, PrimaryHelpOptionsByCommand], click.Command
]
BuildUnavailableCommandFunc = Callable[[str, ImportError], click.Command]


def build_unavailable_command(command_name: str, exc: ImportError) -> click.Command:
    """Return a placeholder command that reports import failure details at runtime."""
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


class DefaultGroup(click.Group):
    """Click group with default subcommand + lazy loading + bool compatibility normalization."""

    def __init__(
        self,
        *args,
        default: str | None = None,
        lazy_subcommands: LazySubcommands | None = None,
        command_bool_value_options: BoolOptionsByCommand | None = None,
        command_bool_toggle_options: BoolOptionsByCommand | None = None,
        command_bool_toggle_negative_aliases: BoolNegativeAliasesByCommand | None = None,
        command_bool_single_flag_options: BoolSingleFlagOptionsByCommand | None = None,
        parser_wrapper_subcommands: frozenset[str] | None = None,
        parser_wrapper_bool_option_providers: ParserWrapperBoolOptionProviders | None = None,
        subcommand_primary_help_options: PrimaryHelpOptionsByCommand | None = None,
        normalize_bool_argv: NormalizeBoolArgvFunc | None = None,
        ensure_help_advanced_option: EnsureHelpAdvancedOptionFunc | None = None,
        configure_subcommand_help_visibility: ConfigureSubcommandHelpVisibilityFunc | None = None,
        build_unavailable_command: BuildUnavailableCommandFunc = build_unavailable_command,
        **kwargs,
    ):
        super().__init__(*args, **kwargs)
        if normalize_bool_argv is None:
            raise ValueError("normalize_bool_argv is required")
        if ensure_help_advanced_option is None:
            raise ValueError("ensure_help_advanced_option is required")
        if configure_subcommand_help_visibility is None:
            raise ValueError("configure_subcommand_help_visibility is required")

        self._default_cmd = default
        self._lazy_subcommands = dict(lazy_subcommands or {})
        self._lazy_cache: dict[str, click.Command] = {}
        self._command_bool_value_options = dict(command_bool_value_options or {})
        self._command_bool_toggle_options = dict(command_bool_toggle_options or {})
        self._command_bool_toggle_negative_aliases = dict(
            command_bool_toggle_negative_aliases or {}
        )
        self._command_bool_single_flag_options = dict(command_bool_single_flag_options or {})
        self._resolved_bool_options_by_command: dict[
            str, tuple[frozenset[str], frozenset[str], dict[str, str], frozenset[str]]
        ] = {}
        self._parser_wrapper_subcommands = set(parser_wrapper_subcommands or set())
        self._parser_wrapper_bool_option_providers = dict(
            parser_wrapper_bool_option_providers or {}
        )
        self._resolved_parser_wrapper_bool_options: dict[str, frozenset[str]] = {}
        self._subcommand_primary_help_options = dict(subcommand_primary_help_options or {})
        self._normalize_bool_argv = normalize_bool_argv
        self._ensure_help_advanced_option = ensure_help_advanced_option
        self._configure_subcommand_help_visibility = configure_subcommand_help_visibility
        self._build_unavailable_command = build_unavailable_command

    def _resolve_parser_wrapper_bool_options(self, command_name: str) -> frozenset[str]:
        cached = self._resolved_parser_wrapper_bool_options.get(command_name)
        if cached is not None:
            return cached

        provider = self._parser_wrapper_bool_option_providers.get(command_name)
        if provider is None:
            resolved = frozenset()
            self._resolved_parser_wrapper_bool_options[command_name] = resolved
            return resolved

        try:
            provided = provider()
        except Exception:
            resolved = frozenset()
            self._resolved_parser_wrapper_bool_options[command_name] = resolved
            return resolved

        normalized: set[str] = set()
        for name in provided:
            if not name.startswith("--"):
                continue
            if name.startswith("--no-"):
                normalized.add(f"--{name[5:]}")
            else:
                normalized.add(name)

        resolved = frozenset(normalized)
        self._resolved_parser_wrapper_bool_options[command_name] = resolved
        return resolved

    @staticmethod
    def _long_option_names(option_names: list[str]) -> tuple[str, ...]:
        return tuple(name for name in option_names if name.startswith("--"))

    def _resolve_bool_options(
        self, ctx: click.Context, command_name: str
    ) -> tuple[frozenset[str], frozenset[str], dict[str, str], frozenset[str]]:
        cached = self._resolved_bool_options_by_command.get(command_name)
        if cached is not None:
            return cached

        bool_value_options = set(
            self._command_bool_value_options.get(command_name, frozenset())
        )
        bool_toggle_options = set(
            self._command_bool_toggle_options.get(command_name, frozenset())
        )
        bool_toggle_negative_aliases = dict(
            self._command_bool_toggle_negative_aliases.get(command_name, {})
        )
        bool_single_flag_options = set(
            self._command_bool_single_flag_options.get(command_name, frozenset())
        )

        bool_toggle_options.update(
            self._resolve_parser_wrapper_bool_options(command_name)
        )

        command = self.get_command(ctx, command_name)
        if command is not None:
            for param in command.params:
                if not isinstance(param, click.Option):
                    continue

                positive_long_names = self._long_option_names(param.opts)
                if not positive_long_names:
                    continue

                negative_long_names = self._long_option_names(param.secondary_opts)
                if param.is_bool_flag:
                    if negative_long_names:
                        bool_toggle_options.update(positive_long_names)
                        first_negative = negative_long_names[0]
                        for index, positive_name in enumerate(positive_long_names):
                            negative_name = (
                                negative_long_names[index]
                                if index < len(negative_long_names)
                                else first_negative
                            )
                            bool_toggle_negative_aliases.setdefault(
                                positive_name, negative_name
                            )
                    else:
                        bool_single_flag_options.update(positive_long_names)
                    continue

                if isinstance(param.type, click.types.BoolParamType):
                    bool_value_options.update(positive_long_names)

        resolved = (
            frozenset(bool_value_options),
            frozenset(bool_toggle_options),
            bool_toggle_negative_aliases,
            frozenset(bool_single_flag_options),
        )
        self._resolved_bool_options_by_command[command_name] = resolved
        return resolved

    def parse_args(self, ctx, args):
        # Normalize long option names to lowercase before any other processing
        # so that --add-linkH, --include-H2O etc. are accepted case-insensitively.
        args = _normalize_argv_option_names(args)

        show_help_or_version = any(a in ("-h", "--help", "--version") for a in args)

        if self._default_cmd is not None and not show_help_or_version:
            if not args or args[0].startswith("-"):
                args.insert(0, self._default_cmd)

        bool_value_options = self._command_bool_value_options
        bool_toggle_options = self._command_bool_toggle_options
        bool_toggle_negative_aliases = self._command_bool_toggle_negative_aliases
        bool_single_flag_options = self._command_bool_single_flag_options
        if args and not args[0].startswith("-"):
            command_name = args[0]
            (
                command_bool_value_options,
                command_bool_toggle_options,
                command_toggle_negative_aliases,
                command_bool_single_flag_options,
            ) = self._resolve_bool_options(ctx, command_name)
            bool_value_options = {command_name: command_bool_value_options}
            bool_toggle_options = {command_name: command_bool_toggle_options}
            bool_toggle_negative_aliases = {
                command_name: command_toggle_negative_aliases
            }
            bool_single_flag_options = {command_name: command_bool_single_flag_options}

        args, _ = self._normalize_bool_argv(
            args,
            bool_value_options,
            bool_toggle_options,
            bool_toggle_negative_aliases,
            bool_single_flag_options,
        )
        return super().parse_args(ctx, args)

    @staticmethod
    def _silence_pysisyphus_loggers():
        """Remove pysisyphus file handlers and suppress log output.

        pysisyphus creates FileHandler + StreamHandler in its various
        ``__init__.py`` files at import time, overriding any prior
        level settings.  This method must run *after* the lazy import
        has loaded the subcommand module (and thus pysisyphus).
        """
        import logging as _logging

        _PYSIS_LOGGERS = (
            "pysisyphus", "calculator", "cos", "dimer", "dynamics",
            "gdiis", "internal_coords", "irc", "optimizer",
            "tsoptimizer", "mwfn", "stocastic", "wfoverlap",
        )
        for name in _PYSIS_LOGGERS:
            lg = _logging.getLogger(name)
            lg.setLevel(_logging.CRITICAL + 1)
            lg.propagate = False
            for h in lg.handlers[:]:
                lg.removeHandler(h)
                try:
                    h.close()
                except Exception:
                    pass

    def invoke(self, ctx):
        # Add a leading blank line for subcommands (except "all") to separate CLI tool logs.
        if ctx.invoked_subcommand and ctx.invoked_subcommand != "all":
            click.echo()
        # Suppress pysisyphus loggers AFTER lazy import has loaded the
        # subcommand module (which triggers pysisyphus __init__ file handlers).
        self._silence_pysisyphus_loggers()
        return super().invoke(ctx)

    def list_commands(self, ctx):
        cmds = set(super().list_commands(ctx))
        cmds.update(self._lazy_subcommands.keys())
        return sorted(cmds)

    def get_command(self, ctx, cmd_name):
        cmd = super().get_command(ctx, cmd_name)
        if cmd is not None:
            if cmd_name not in self._parser_wrapper_subcommands:
                cmd = self._ensure_help_advanced_option(cmd)
            return self._configure_subcommand_help_visibility(
                cmd_name, cmd, self._subcommand_primary_help_options
            )

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
            loaded_cmd = self._build_unavailable_command(cmd_name, exc)

        if cmd_name not in self._parser_wrapper_subcommands:
            loaded_cmd = self._ensure_help_advanced_option(loaded_cmd)
        loaded_cmd = self._configure_subcommand_help_visibility(
            cmd_name, loaded_cmd, self._subcommand_primary_help_options
        )
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
