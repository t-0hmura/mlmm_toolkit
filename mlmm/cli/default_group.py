"""Shared Click group helpers for lazy subcommand loading and bool normalization."""

from __future__ import annotations

import contextlib
import importlib
import os
import select
import shlex
import signal
import subprocess
import sys
import threading
from collections.abc import Callable, Mapping
from pathlib import Path

import click

from mlmm.cli.bool_compat import normalize_argv_option_names as _normalize_argv_option_names
from mlmm.core.prerequisites import exit_if_missing_prerequisite

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
# Verbose-option injector also receives the registry command name so the help
# text can be specialised per subcommand (lazy commands share cmd.name="cli").
EnsureVerboseOptionFunc = Callable[[click.Command, "str | None"], click.Command]
ConfigureSubcommandHelpVisibilityFunc = Callable[
    [str, click.Command, PrimaryHelpOptionsByCommand], click.Command
]
BuildUnavailableCommandFunc = Callable[[str, ImportError], click.Command]

_INTERNAL_IMPORT_ROOTS = frozenset(
    {"mlmm", "pysisyphus", "thermoanalysis", "hessian_ff"}
)


_C_FLUSH_INTERVAL_SECONDS = 5.0


def _write_all(write, data) -> None:
    view = memoryview(data)
    while view:
        view = view[write(view):]


class _RunLogSink:
    """Append the command's output to ``<out-dir>/run.log`` as it is produced.

    Output that arrives before the command creates its output directory is
    held in memory and written once the directory exists. A rerun into the
    same directory appends after the previous run's log.
    """

    def __init__(self, directory: Path, header: str):
        self._path = directory / "run.log"
        self._pending = bytearray(header.encode("utf-8"))
        self._handle = None
        self._closed = False
        self._lock = threading.Lock()
        self.error: OSError | None = None

    def _ready(self) -> bool:
        if self._handle is not None:
            return True
        if not self._path.parent.is_dir():
            return False
        handle = open(self._path, "ab", buffering=0)
        try:
            if handle.tell() > 0:
                _write_all(handle.write, b"\n")
            _write_all(handle.write, bytes(self._pending))
        except OSError:
            handle.close()
            raise
        self._pending.clear()
        self._handle = handle
        return True

    def write(self, data: bytes) -> None:
        with self._lock:
            if self._closed or self.error is not None:
                return
            try:
                if self._ready():
                    _write_all(self._handle.write, data)
                else:
                    self._pending += data
            except OSError as exc:
                self.error = exc

    def hand_over(self) -> None:
        """Leave the log to the run-log helper process; write nothing more."""
        with self._lock:
            self._closed = True

    def close(self) -> None:
        with self._lock:
            if self._closed:
                return
            self._closed = True
            try:
                if self.error is None:
                    self._ready()
            except OSError as exc:
                self.error = exc
            if self._handle is not None:
                self._handle.close()
                self._handle = None


class _TeeText:
    """Mirror text to a replaced sys.stdout/sys.stderr and to the run log."""

    def __init__(self, stream, sink: _RunLogSink):
        self._stream = stream
        self._sink = sink

    def write(self, value):
        self._sink.write(value.encode("utf-8", errors="replace"))
        return self._stream.write(value)

    def flush(self):
        self._stream.flush()

    def __getattr__(self, name):
        return getattr(self._stream, name)


class _FdTee:
    """Copy everything written to file descriptors 1 and 2 (by Python, C
    libraries or child processes) to their original destinations and to the
    run log.

    A helper process does the copying, so the output written just before the
    command dies by a signal (an abort or a segfault in C code) still reaches
    both.
    """

    def __init__(self, sink: _RunLogSink):
        pipes = [os.pipe() for _ in range(4)]
        (out_read, out_write), (err_read, err_write), (control_read, control_write), \
            (ready_read, ready_write) = pipes
        child_ends = (out_read, err_read, control_read, ready_write)
        try:
            self._process = subprocess.Popen(
                [
                    sys.executable, "-c", _PUMP_SCRIPT,
                    str(Path(__file__).resolve().parents[2]), __name__,
                    *(str(fd) for fd in child_ends),
                    str(sink._path.parent), bytes(sink._pending).decode("utf-8"),
                ],
                stdin=subprocess.DEVNULL,
                pass_fds=child_ends,
            )
        except (OSError, subprocess.SubprocessError, ValueError) as exc:
            for fd in (fd for pipe in pipes for fd in pipe):
                os.close(fd)
            raise OSError(f"cannot start the run-log helper: {exc}") from exc
        for fd in child_ends:
            os.close(fd)
        started = bool(select.select([ready_read], [], [], 30.0)[0]) and os.read(ready_read, 1) == b"R"
        os.close(ready_read)
        if not started:
            self._process.kill()
            self._process.wait()
            for fd in (out_write, err_write, control_write):
                os.close(fd)
            raise OSError("the run-log helper did not start")
        sink.hand_over()
        self._saved = [os.dup(1), os.dup(2)]
        os.dup2(out_write, 1)
        os.dup2(err_write, 2)
        os.close(out_write)
        os.close(err_write)
        self._control = control_write

    def restore(self) -> None:
        os.dup2(self._saved[0], 1)
        os.dup2(self._saved[1], 2)

    def finish(self) -> None:
        # Closing the control pipe makes the helper copy what is left in the
        # pipes and exit, even when a process that outlives the command (such
        # as multiprocessing's resource tracker) still holds them open.
        os.close(self._control)
        with contextlib.suppress(subprocess.TimeoutExpired):
            self._process.wait(timeout=5.0)
        for fd in self._saved:
            os.close(fd)


_PUMP_SCRIPT = (
    "import importlib, sys; sys.path.insert(0, sys.argv[1]); "
    "importlib.import_module(sys.argv[2])._pump_run_log(sys.argv[3:])"
)


def _pump_run_log(argv: list[str]) -> None:
    """Run the helper process of :class:`_FdTee`."""
    out_read, err_read, control_read, ready_write = (int(value) for value in argv[:4])
    # Stop when the command closes the control pipe or dies, not on a signal
    # sent to its whole process group (Ctrl-C, a batch system's SIGTERM).
    for signum in (signal.SIGINT, signal.SIGTERM, signal.SIGHUP):
        signal.signal(signum, signal.SIG_IGN)
    sink = _RunLogSink(Path(argv[4]), argv[5])
    os.write(ready_write, b"R")
    os.close(ready_write)
    destinations = {out_read: 1, err_read: 2}

    def copy(fd: int) -> bool:
        try:
            data = os.read(fd, 65536)
        except BlockingIOError:
            return False
        except OSError:
            data = b""
        if not data:
            os.close(fd)
            del destinations[fd]
            return False
        try:
            _write_all(lambda chunk: os.write(destinations[fd], chunk), data)
        except OSError:
            pass
        sink.write(data)
        return True

    while True:
        readable = select.select([*destinations, control_read], [], [])[0]
        for fd in readable:
            if fd in destinations:
                copy(fd)
        if control_read in readable:
            break
    # The command has finished or died: copy what is left in the pipes.
    for fd in list(destinations):
        os.set_blocking(fd, False)
        while copy(fd):
            pass
    sink.close()
    if sink.error is not None:
        message = f"[run.log] WARNING: Could not save run.log: {sink.error}\n"
        with contextlib.suppress(OSError):
            _write_all(lambda chunk: os.write(2, chunk), message.encode("utf-8"))


def _libc_fflush():
    """Return a callable that flushes the C stdio buffers, or None."""
    try:
        import ctypes

        fflush = ctypes.CDLL(None).fflush
    except (OSError, AttributeError):
        return None
    fflush.argtypes = [ctypes.c_void_p]
    return lambda: fflush(None)


def _start_fd_tees(sink: _RunLogSink, fflush) -> list[_FdTee]:
    try:
        if sys.stdout.fileno() != 1 or sys.stderr.fileno() != 2:
            return []
    except (AttributeError, OSError, ValueError):
        return []
    for stream in (sys.stdout, sys.stderr):
        stream.flush()
    if fflush is not None:
        fflush()
    try:
        return [_FdTee(sink)]
    except OSError:
        return []


@contextlib.contextmanager
def _capture_run_log(sink: _RunLogSink):
    """Send stdout and stderr to the run log while the command runs."""
    fflush = _libc_fflush()
    tees = _start_fd_tees(sink, fflush)
    if not tees:
        # In-process callers (tests, notebooks) replace sys.stdout/sys.stderr.
        with contextlib.redirect_stdout(_TeeText(sys.stdout, sink)), \
                contextlib.redirect_stderr(_TeeText(sys.stderr, sink)):
            yield
        return

    # Python output reaches run.log line by line; C stdio buffers (for
    # example Ipopt in DMF) are flushed every few seconds.
    line_buffering = getattr(sys.stdout, "line_buffering", None)
    with contextlib.suppress(AttributeError, ValueError):
        sys.stdout.reconfigure(line_buffering=True)
    stop = threading.Event()
    flusher = None
    if fflush is not None:
        def _flush_c_streams() -> None:
            while not stop.wait(_C_FLUSH_INTERVAL_SECONDS):
                fflush()

        flusher = threading.Thread(target=_flush_c_streams, name="run-log-flush", daemon=True)
        flusher.start()
    try:
        yield
    finally:
        stop.set()
        if flusher is not None:
            flusher.join(timeout=1.0)
        for stream in (sys.stdout, sys.stderr):
            with contextlib.suppress(Exception):
                stream.flush()
        if fflush is not None:
            fflush()
        for tee in tees:
            tee.restore()
        for tee in tees:
            tee.finish()
        if line_buffering is not None:
            with contextlib.suppress(AttributeError, ValueError):
                sys.stdout.reconfigure(line_buffering=line_buffering)


def _is_external_missing_dependency(exc: ModuleNotFoundError) -> bool:
    """Return whether *exc* names a dependency outside this distribution.

    A missing ``mlmm`` module is an internal packaging/import defect and must
    retain its traceback.  Only a genuinely external missing dependency is
    eligible for the friendly lazy-command placeholder.
    """

    missing = str(exc.name or "")
    missing_root = missing.split(".", 1)[0]
    return bool(missing_root) and missing_root not in _INTERNAL_IMPORT_ROOTS


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
        ensure_verbose_option: EnsureVerboseOptionFunc | None = None,
        configure_subcommand_help_visibility: ConfigureSubcommandHelpVisibilityFunc | None = None,
        build_unavailable_command: BuildUnavailableCommandFunc = build_unavailable_command,
        command_groups: "dict[str, tuple[str, ...]] | None" = None,
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
        self._ensure_verbose_option = ensure_verbose_option
        self._configure_subcommand_help_visibility = configure_subcommand_help_visibility
        self._build_unavailable_command = build_unavailable_command
        # Optional command grouping for `--help`. Keys are display section
        # titles ("Pipelines" / "Pipeline stages" / ...) in render order;
        # values are subcommand names to bucket under each title. Subcommands
        # absent from every bucket fall through to a final "Other" section.
        self._command_groups: dict[str, tuple[str, ...]] = {
            title: tuple(members) for title, members in (command_groups or {}).items()
        }

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
                    elif param.default is not True:
                        # A positive-only, default-on flag is an explicit marker,
                        # not a toggle. Do not synthesize a --no-* spelling that
                        # the command deliberately does not expose.
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

    def _reject_conflicting_aliases(
        self, ctx: click.Context, args: list[str]
    ) -> None:
        """Reject different values supplied through aliases of one option."""

        if not args or args[0].startswith("-"):
            return
        command = self.get_command(ctx, args[0])
        if command is None:
            return
        for option in command.params:
            if not isinstance(option, click.Option) or option.multiple or option.nargs != 1:
                continue
            positives = [name for name in option.opts if name.startswith("--")]
            negatives = [name for name in option.secondary_opts if name.startswith("--")]
            if len(positives) < 2 and len(negatives) < 2:
                continue
            accepted = set(positives) | set(negatives)
            seen: list[tuple[str, object]] = []
            index = 1
            while index < len(args):
                token = args[index]
                name, separator, inline = token.partition("=")
                if name not in accepted:
                    index += 1
                    continue
                if option.is_bool_flag:
                    value: object = name in positives
                    index += 1
                else:
                    if separator:
                        raw = inline
                        index += 1
                    elif index + 1 < len(args):
                        raw = args[index + 1]
                        index += 2
                    else:
                        index += 1
                        continue
                    try:
                        value = option.type.convert(raw, option, ctx)
                    except Exception:
                        value = raw
                seen.append((name, value))
            if len({name for name, _ in seen}) > 1 and any(
                value != seen[0][1] for _, value in seen[1:]
            ):
                aliases = ", ".join(sorted({name for name, _ in seen}))
                raise click.UsageError(
                    f"Conflicting values were supplied through aliases: {aliases}."
                )

    def parse_args(self, ctx, args):
        # Reset the console gate and pipeline mode at the start of every
        # top-level parse; the eager -v callback (only on the real CLI group)
        # re-enables the gate inside super().parse_args(), and the `all` command
        # re-sets pipeline mode. This scopes both to one invocation rather than
        # letting them leak between repeated runs / CliRunner test calls.
        from mlmm.core.utils import (
            set_allow_charge_mult_mismatch,
            set_console_gating,
            set_pipeline_mode,
        )
        set_allow_charge_mult_mismatch(False)
        set_console_gating(False)
        set_pipeline_mode(False)
        # Normalize long option names to lowercase before any other processing
        # so that --add-linkH, --include-H2O etc. are accepted case-insensitively.
        args = _normalize_argv_option_names(args)

        root_help_or_version = bool(
            args and args[0] in ("-h", "--help", "--version")
        )

        if self._default_cmd is not None and not root_help_or_version:
            # Without this guard, `mlmm -v tsopt ...` would prepend the
            # default command and swallow `tsopt` as an argument to it.
            # Collect every spelling of every top-level option declared on
            # this group so that a leading group-level flag (e.g. `-v` /
            # does not trigger the default-command fallback.
            top_level_opts = set()
            for p in self.params:
                top_level_opts.update(getattr(p, "opts", ()) or ())
                top_level_opts.update(getattr(p, "secondary_opts", ()) or ())
            first = args[0] if args else ""
            if first.startswith("--"):
                first_opt = first.split("=", 1)[0]
                is_top_level = first_opt in top_level_opts
            elif first.startswith("-") and len(first) >= 2:
                # For a short-option cluster, the first two characters identify
                # the option Click examines first.
                is_top_level = first[:2] in top_level_opts
            else:
                is_top_level = False
            first_opt = first.split("=", 1)[0]
            known_commands = set(self.commands) | set(self._lazy_subcommands)
            command_after_verbose = (
                (first_opt == "--verbose" or first[:2] == "-v")
                and any(arg in known_commands for arg in args[1:])
            )
            if not args or (
                args[0].startswith("-")
                and not is_top_level
                and not command_after_verbose
            ):
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
        self._reject_conflicting_aliases(ctx, args)
        # Preserve the normalized token stream for commands with historical
        # variadic options (for example ``-i A B``).  Nested Click contexts
        # share ``meta``, unlike ``sys.argv`` under in-process invocation.
        ctx.meta["mlmm.cli.raw_args"] = tuple(args)
        result = super().parse_args(ctx, args)
        # Any help/version request must print in full: the subcommand's eager
        # --help-advanced renderer fires later (during invoke) and routes through
        # the patched echo, so re-disable the gate that the -v callback just
        # turned on. (Plain -h/--help/--version exit during the super() call
        # above and are unaffected, but --help-advanced is not a built-in flag.)
        if any(a in ("-h", "--help", "--version", "--help-advanced") for a in args):
            set_console_gating(False)
        return result

    @staticmethod
    def _silence_pysisyphus_loggers():
        """Remove pysisyphus file handlers and suppress log output.

        pysisyphus creates FileHandler + StreamHandler in its various
        ``__init__.py`` files at import time, overriding any prior
        level settings.  This method must run *after* the lazy import
        has loaded the subcommand module (and thus pysisyphus).
        """
        import logging as _logging

        # pysisyphus has no public aggregate logger, so enumerate the loggers
        # whose handlers are installed during its module imports.
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

    def _invoke_command(self, ctx):
        """Run the selected command; an error that escapes it and names a missing system
        prerequisite ends with the install guidance instead of a traceback."""
        try:
            return super().invoke(ctx)
        except Exception as exc:
            exit_if_missing_prerequisite(exc)
            raise

    def invoke(self, ctx):
        # Add a leading blank line for subcommands (except "all") to separate CLI tool logs.
        if ctx.invoked_subcommand and ctx.invoked_subcommand != "all":
            click.echo()
        # Suppress pysisyphus loggers AFTER lazy import has loaded the
        # subcommand module (which triggers pysisyphus __init__ file handlers).
        self._silence_pysisyphus_loggers()
        run_log_dir = self._run_log_directory(ctx)
        if run_log_dir is None:
            return self._invoke_command(ctx)

        raw = list(ctx.meta.get("mlmm.cli.raw_args") or ())
        sink = _RunLogSink(run_log_dir, f"[command] {shlex.join(['mlmm', *raw])}\n")
        try:
            with _capture_run_log(sink):
                return self._invoke_command(ctx)
        finally:
            sink.close()
            if sink.error is not None:
                click.echo(f"[run.log] WARNING: Could not save run.log: {sink.error}", err=True)

    def _run_log_directory(self, ctx) -> Path | None:
        """Resolve the selected command's Click-owned output directory."""

        from mlmm.core.utils import is_child_mode

        raw = list(ctx.meta.get("mlmm.cli.raw_args") or ())
        if (
            is_child_mode()
            or not raw
            or any(token in ("-h", "--help", "--help-advanced", "--version", "--dry-run")
                   for token in raw)
        ):
            return None
        command_name = raw[0]
        command = self.get_command(ctx, command_name)
        if command is None:
            return None
        output_option = next(
            (
                param
                for param in command.params
                if isinstance(param, click.Option) and param.name == "out_dir"
            ),
            None,
        )
        if output_option is None:
            return None
        try:
            sub_ctx = click.Context(
                command,
                info_name=command_name,
                parent=ctx,
                resilient_parsing=True,
            )
            parsed, _remaining, _order = command.make_parser(sub_ctx).parse_args(
                args=raw[1:]
            )
            value = parsed.get(output_option.name, output_option.default)
        except Exception:
            return None
        output_spellings = set(output_option.opts) | set(output_option.secondary_opts)
        output_explicit = any(
            token in output_spellings
            or any(token.startswith(f"{spelling}=") for spelling in output_spellings)
            or (token.startswith("-o") and "-o" in output_spellings and token != "-o")
            for token in raw[1:]
        )
        if not output_explicit and parsed.get("config_yaml"):
            try:
                from mlmm.core.utils import load_yaml_dict

                config = load_yaml_dict(Path(parsed["config_yaml"]))
                section_names = (
                    command_name,
                    command_name.replace("-", "_"),
                    "search" if command_name == "path-search" else command_name,
                )
                for section_name in section_names:
                    section = config.get(section_name)
                    if isinstance(section, Mapping) and section.get("out_dir") is not None:
                        value = section["out_dir"]
                        break
            except Exception:
                pass
        if value is None:
            return None
        return Path(os.path.abspath(os.path.expanduser(os.fspath(value))))

    def list_commands(self, ctx):
        # Preserve the semantic order declared in `_LAZY_SUBCOMMANDS` so the
        # canonical pipeline (`all`, `extract`, `mm-parm`, `opt`, `tsopt`, ...)
        # shows up first in `--help` instead of being scattered alphabetically.
        # Any externally-registered subcommands appear at the end, sorted.
        lazy_order = list(self._lazy_subcommands.keys())
        all_cmds = set(super().list_commands(ctx))
        all_cmds.update(lazy_order)
        extras = sorted(all_cmds - set(lazy_order))
        return lazy_order + extras

    def get_command(self, ctx, cmd_name):
        cmd = super().get_command(ctx, cmd_name)
        if cmd is not None:
            if cmd_name not in self._parser_wrapper_subcommands:
                cmd = self._ensure_help_advanced_option(cmd)
            cmd = self._configure_subcommand_help_visibility(
                cmd_name, cmd, self._subcommand_primary_help_options
            )
            # Inject `-v` AFTER the advanced-hide pass so verbosity stays
            # visible in the primary --help (it is a core, commonly-used option).
            if (
                cmd_name not in self._parser_wrapper_subcommands
                and self._ensure_verbose_option is not None
            ):
                cmd = self._ensure_verbose_option(cmd, cmd_name)
            return cmd

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
        except ModuleNotFoundError as exc:
            if not _is_external_missing_dependency(exc):
                raise
            loaded_cmd = self._build_unavailable_command(cmd_name, exc)

        if cmd_name not in self._parser_wrapper_subcommands:
            loaded_cmd = self._ensure_help_advanced_option(loaded_cmd)
        loaded_cmd = self._configure_subcommand_help_visibility(
            cmd_name, loaded_cmd, self._subcommand_primary_help_options
        )
        # Inject `-v` AFTER the advanced-hide pass so verbosity stays visible
        # in the primary --help (it is a core, commonly-used option).
        if (
            cmd_name not in self._parser_wrapper_subcommands
            and self._ensure_verbose_option is not None
        ):
            loaded_cmd = self._ensure_verbose_option(loaded_cmd, cmd_name)
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

        if not rows:
            return

        groups = self._command_groups
        if not groups:
            with formatter.section("Commands"):
                formatter.write_dl(rows)
            return

        # Render commands bucketed by `_command_groups`; any subcommand not
        # named in a bucket falls through to a final "Other" section so we
        # never silently hide subcommands.
        row_by_name = dict(rows)
        seen: set[str] = set()
        for group_title, members in groups.items():
            group_rows = [
                (m, row_by_name[m]) for m in members if m in row_by_name
            ]
            if not group_rows:
                continue
            seen.update(name for name, _ in group_rows)
            with formatter.section(group_title):
                formatter.write_dl(group_rows)

        leftover = [(name, help_) for name, help_ in rows if name not in seen]
        if leftover:
            with formatter.section("Other"):
                formatter.write_dl(leftover)
