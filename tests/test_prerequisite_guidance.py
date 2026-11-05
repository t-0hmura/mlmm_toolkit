"""A missing compiler or Chrome library ends with a message that says what to install."""

from __future__ import annotations

import re
from pathlib import Path

import click
import pytest
from click.testing import CliRunner

from mlmm.core import prerequisites as pre
from mlmm.io import plotly_image

# The patterns a gate or a user can search for; one per kind.
HEADERS = {
    pre.KIND_C_COMPILER: r"Missing system prerequisite: C compiler \(gcc\)\.",
    pre.KIND_CXX_COMPILER: r"Missing system prerequisite: C\+\+ compiler \(g\+\+\)\.",
    pre.KIND_CUDA: r"Missing system prerequisite: CUDA build tools\.",
    pre.KIND_CHROME: r"Missing system prerequisite: Chrome for PNG export\.",
}

# Errors exactly as the container run printed them (torch 2.13.0, Triton 3.7.1, Kaleido 1.5.0).
TRITON_C = (
    "Failed to find C compiler. Please specify via CC environment variable "
    "or set triton.knobs.build.impl."
)
INDUCTOR_C = f"RuntimeError: {TRITON_C}"
INDUCTOR_CXX = (
    "InvalidCxxCompiler: No working C++ compiler found in "
    "torch._inductor.config.cpp.cxx: (None, 'g++')"
)
BROWSER = (
    "BrowserFailedError: ('The browser seemed to close immediately after starting.', "
    "'You can set the `logging.Logger` level lower to see more output.')"
)
# Plotly turns Kaleido's ChromeNotFoundError into a RuntimeError with this text (plotly.io._kaleido, Plotly 7.1.0).
PLOTLY_NO_CHROME = (
    "\n\nKaleido requires Google Chrome to be installed.\n\n"
    "Either download and install Chrome yourself following Google's instructions for your operating system,\n"
    "or install it from your terminal by running:\n\n    $ plotly_get_chrome\n\n"
)


class InductorError(Exception):
    """Stand-in with the name and message shape of torch._inductor.exc.InductorError."""


def _wrapped(inner: str) -> Exception:
    return InductorError(inner)


@pytest.fixture(autouse=True)
def _default_console_state(monkeypatch):
    """Run every test at the default verbosity with the console gate off.

    The CLI keeps its verbosity and console gate in module globals and never resets them, so other tests leave them
    behind: test_utils.py-style tests end at -v 0, which silences even errors once any CLI test has installed the patched
    click.echo. These tests read the messages from stderr, so they must not depend on what ran before; monkeypatch puts
    the previous values back afterwards.
    """
    from mlmm.core import utils

    for module, name, value in (
        (utils, "_VERBOSE_LEVEL", 2),
        (utils, "_GATE_ACTIVE", False),
        (utils, "_PIPELINE_MODE", False),
        (utils, "_CHILD_MODE", False),
    ):
        monkeypatch.setattr(module, name, value)


# ---------------------------------------------------------------------------
# The message
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("kind", sorted(HEADERS))
def test_message_names_what_is_missing_and_where_to_read_more(kind):
    text = pre.guidance(kind, "RuntimeError: original failure")
    assert re.match(HEADERS[kind], text)
    assert [other for other, pattern in HEADERS.items() if re.search(pattern, text)] == [kind]
    assert f"  Details: {pre.DOCS_URL}" in text
    assert text.endswith("  Original error: RuntimeError: original failure")


def test_message_links_to_a_docs_section_that_exists():
    page = Path(__file__).parents[1] / "docs" / "installation.md"
    assert pre.DOCS_URL.endswith("/installation.html#system-prerequisites")
    assert "(system-prerequisites)=\n" in page.read_text(encoding="utf-8")


def test_compiler_messages_give_the_install_commands():
    gpu = pre.guidance(pre.KIND_C_COMPILER)
    assert "UMA on a GPU, and ORB with --compile, compile code at run time for speed" in gpu
    cpu = pre.guidance(pre.KIND_CXX_COMPILER)
    assert "With --compile, UMA and ORB compile CPU kernels at run time for speed" in cpu
    for kind in (pre.KIND_C_COMPILER, pre.KIND_CXX_COMPILER, pre.KIND_CUDA):
        text = pre.guidance(kind)
        assert "conda install -c conda-forge cuda-nvcc" in text
        assert "reactivate the environment" in text


def test_chrome_message_gives_both_steps_and_both_package_lists():
    text = pre.guidance(pre.KIND_CHROME)
    assert "plotly_get_chrome -y" in text
    assert text.index("sudo apt-get update") < text.index("sudo apt-get install -y ca-certificates")
    assert "libnss3" in text and "libgbm1" in text and "xdg-utils" in text
    assert "sudo dnf install -y alsa-lib" in text
    assert "nss" in text and "mesa-libgbm" in text and "vulkan-loader" in text
    # RHEL 9 images carry libcurl-minimal, which provides libcurl.so.4 and conflicts with the libcurl package.
    assert "libcurl" not in text.split("Rocky/RHEL:", 1)[1].split()


def test_original_error_is_one_bounded_line():
    long = "RuntimeError: " + "x" * 500 + "\nsecond line"
    line = pre.guidance(pre.KIND_C_COMPILER, long).splitlines()[-1]
    assert line.startswith("  Original error: RuntimeError: xxx")
    assert line.endswith("...") and len(line) < 330


# ---------------------------------------------------------------------------
# Recognizing the errors
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    ("error", "kind"),
    [
        (RuntimeError(TRITON_C), pre.KIND_C_COMPILER),
        (_wrapped(INDUCTOR_C), pre.KIND_C_COMPILER),
        (_wrapped(INDUCTOR_CXX), pre.KIND_CXX_COMPILER),
        (_wrapped("fatal error: cuda.h: No such file or directory"), pre.KIND_CUDA),
        (_wrapped("Cannot find ptxas"), pre.KIND_CUDA),
        (RuntimeError("ptxas: command not found"), pre.KIND_CUDA),
        # Triton's launcher build without cuda.h (the conda env of 2026-10-25 without cuda-nvcc).
        (_wrapped("CalledProcessError: Command '['/env/bin/gcc', '/tmp/tmp9c2eap4f/cuda_utils.c', '-O3', "
                  "'-shared']' returned non-zero exit status 1."), pre.KIND_CUDA),
        (RuntimeError(BROWSER), pre.KIND_CHROME),
        (pre.MissingPrerequisiteError(pre.KIND_CHROME, "Plotly image export failed"), pre.KIND_CHROME),
    ],
)
def test_known_failures_get_their_message(error, kind):
    text = pre.prerequisite_message(error)
    assert text is not None and re.search(HEADERS[kind], text)


def test_chrome_error_class_names_are_recognized():
    class BrowserFailedError(RuntimeError): ...

    class ChromeNotFoundError(RuntimeError): ...

    assert pre.is_chrome_failure(BROWSER)
    assert pre.prerequisite_message(BrowserFailedError("closed")) is not None
    assert pre.prerequisite_message(ChromeNotFoundError("Browser not found")) is not None


@pytest.mark.parametrize(
    "error",
    [
        ValueError("bad value"),
        RuntimeError("Plotly image export failed: TimeoutError: browser configuration lock is held"),
        RuntimeError("compiler flags are invalid"),
        RuntimeError("CUDA out of memory"),
        RuntimeError("ptxas fatal: invalid value for --gpu-name"),
        click.UsageError("missing option"),
        click.ClickException("Path pre-alignment failed: shapes differ"),
        plotly_image.PlotlyImageTimeoutError("Plotly image export timed out after 600 s"),
    ],
)
def test_other_errors_are_left_alone(error):
    assert pre.prerequisite_message(error) is None


def test_cause_and_context_are_searched():
    try:
        try:
            raise _wrapped(INDUCTOR_C)
        except InductorError:
            raise RuntimeError("[dmf] optimization failed")
    except RuntimeError as implicit:
        assert re.search(HEADERS[pre.KIND_C_COMPILER], pre.prerequisite_message(implicit))

    explicit = click.ClickException("Path pre-alignment failed")
    explicit.__cause__ = _wrapped(INDUCTOR_CXX)
    assert re.search(HEADERS[pre.KIND_CXX_COMPILER], pre.prerequisite_message(explicit))

    loop = RuntimeError("loop")
    loop.__cause__ = loop
    assert pre.prerequisite_message(loop) is None


def test_the_original_error_line_quotes_the_innermost_match():
    try:
        try:
            raise _wrapped(INDUCTOR_C)
        except InductorError as inner:
            raise click.ClickException(f"Path pre-alignment failed: {inner}") from inner
    except click.ClickException as wrapper:
        last = pre.prerequisite_message(wrapper).splitlines()[-1]
    assert last.startswith("  Original error: InductorError: RuntimeError: Failed to find C compiler")

    figure = pre.MissingPrerequisiteError(pre.KIND_CHROME, "Plotly image export failed: BrowserFailedError: x")
    wrapped = click.ClickException(str(figure))
    wrapped.__cause__ = figure
    last = pre.prerequisite_message(wrapped).splitlines()[-1]
    assert last == "  Original error: Plotly image export failed: BrowserFailedError: x"


def test_a_suppressed_context_is_not_searched():
    try:
        try:
            raise _wrapped(INDUCTOR_C)
        except InductorError:
            raise RuntimeError("something else") from None
    except RuntimeError as unrelated:
        assert pre.prerequisite_message(unrelated) is None


def test_the_prerequisite_error_survives_pickling():
    import pickle

    original = pre.MissingPrerequisiteError(pre.KIND_CHROME, "Plotly image export failed")
    restored = pickle.loads(pickle.dumps(original))
    assert restored.kind == pre.KIND_CHROME and str(restored) == str(original)


def test_exit_if_missing_prerequisite_exits_one_and_prints_to_stderr(capsys):
    pre.exit_if_missing_prerequisite(ValueError("not ours"))
    assert capsys.readouterr().err == ""
    with pytest.raises(SystemExit) as stop:
        pre.exit_if_missing_prerequisite(_wrapped(INDUCTOR_C))
    assert stop.value.code == 1
    assert re.search(HEADERS[pre.KIND_C_COMPILER], capsys.readouterr().err)


# The patterns above are copied from the container run; these use the real classes.


def test_triton_without_a_compiler_raises_the_matched_error(monkeypatch, tmp_path):
    build = pytest.importorskip("triton.runtime.build")
    monkeypatch.delenv("CC", raising=False)
    monkeypatch.setenv("PATH", str(tmp_path))
    with pytest.raises(RuntimeError) as caught:
        build.compile_module_from_src(f"/* {tmp_path.name} */", "probe")
    assert re.search(HEADERS[pre.KIND_C_COMPILER], pre.prerequisite_message(caught.value))


def test_inductor_without_a_cxx_compiler_raises_the_matched_error():
    inductor = pytest.importorskip("torch._inductor.exc")
    error = inductor.InductorError(inductor.InvalidCxxCompiler(), None)
    assert re.search(HEADERS[pre.KIND_CXX_COMPILER], pre.prerequisite_message(error))


def test_kaleido_errors_are_the_matched_errors():
    errors = pytest.importorskip("choreographer.errors")
    for error in (
        errors.BrowserFailedError("The browser seemed to close immediately after starting."),
        errors.ChromeNotFoundError("Browser not found. You can use get_chrome()."),
    ):
        assert re.search(HEADERS[pre.KIND_CHROME], pre.prerequisite_message(error))


# ---------------------------------------------------------------------------
# PNG export (the one place that writes figures)
# ---------------------------------------------------------------------------


class _Figure:
    def to_json(self) -> str:
        return "{}"


class _Connection:
    def __init__(self, result=None) -> None:
        self.result = result

    def poll(self) -> bool:
        return self.result is not None

    def recv(self):
        return self.result

    def close(self) -> None:
        pass


class _Process:
    exitcode = 0

    def __init__(self, **kwargs) -> None:
        pass

    def start(self) -> None:
        pass

    def join(self, timeout=None) -> None:
        pass

    def is_alive(self) -> bool:
        return False


def _render_returns(monkeypatch, result) -> None:
    class Context:
        def Pipe(self, duplex=False):
            return _Connection(result), _Connection()

        Process = _Process

    monkeypatch.setattr(plotly_image.mp, "get_context", lambda method: Context())


def test_png_export_that_cannot_start_chrome_raises_the_prerequisite_error(monkeypatch, tmp_path):
    _render_returns(monkeypatch, ("error", BROWSER))
    with pytest.raises(pre.MissingPrerequisiteError) as caught:
        plotly_image.write_plotly_image(_Figure(), tmp_path / "plot.png")
    assert isinstance(caught.value, RuntimeError)
    assert caught.value.kind == pre.KIND_CHROME
    text = str(caught.value)  # one line: all and path-search print a skipped figure's error inline
    assert "\n" not in text and pre.DOCS_URL in text
    assert re.search(HEADERS[pre.KIND_CHROME], text)
    assert "Plotly image export failed: BrowserFailedError" in text
    assert not (tmp_path / "plot.png").exists()


def test_png_export_without_chrome_raises_the_prerequisite_error(monkeypatch, tmp_path):
    error = RuntimeError(PLOTLY_NO_CHROME)
    _render_returns(monkeypatch, ("error", f"{type(error).__name__}: {error}"))  # what the worker sends
    with pytest.raises(pre.MissingPrerequisiteError) as caught:
        plotly_image.write_plotly_image(_Figure(), tmp_path / "plot.png")
    assert caught.value.kind == pre.KIND_CHROME
    assert "Kaleido requires Google Chrome to be installed" in str(caught.value)


def test_png_export_failing_for_another_reason_keeps_its_plain_error(monkeypatch, tmp_path):
    _render_returns(monkeypatch, ("error", "ValueError: render failed"))
    with pytest.raises(RuntimeError, match="Plotly image export failed: ValueError: render failed") as caught:
        plotly_image.write_plotly_image(_Figure(), tmp_path / "plot.png")
    assert not isinstance(caught.value, pre.MissingPrerequisiteError)


# ---------------------------------------------------------------------------
# What the commands print: the guidance, exit status 1, no traceback
# ---------------------------------------------------------------------------


def _assert_guidance_only(output: str, kind: str) -> None:
    assert re.search(HEADERS[kind], output), output
    assert "Traceback" not in output and "Unhandled" not in output, output


def test_energy_diagram_without_chrome_libraries_prints_the_guidance(monkeypatch, tmp_path):
    from mlmm.cli import cli as root_cli
    from mlmm.io import energy_diagram

    def no_chrome(fig, path, **kwargs):
        raise pre.MissingPrerequisiteError(pre.KIND_CHROME, f"Plotly image export failed: {BROWSER}")

    monkeypatch.setattr(energy_diagram, "write_plotly_image", no_chrome)
    out = tmp_path / "energy.png"
    result = CliRunner().invoke(root_cli, ["energy-diagram", "-i", "0", "-i", "5", "-o", str(out)])
    assert result.exit_code == 1
    # CliRunner turns an escaping exception into exit 1 without a traceback.
    assert isinstance(result.exception, SystemExit)
    _assert_guidance_only(result.output, pre.KIND_CHROME)
    assert "sudo apt-get install -y" in result.output and "Original error: Plotly image export failed" in result.output
    assert not out.exists()


def test_an_error_that_escapes_a_command_reaches_the_guidance(monkeypatch, tmp_path):
    """The all pipeline has no handler of its own: a wrapped error is found through its cause."""
    from mlmm.cli import cli as root_cli
    from mlmm.io import energy_diagram

    def wrapped(**kwargs):
        try:
            raise _wrapped(INDUCTOR_CXX)
        except InductorError as exc:
            raise click.ClickException(f"Path pre-alignment failed: {exc}") from exc

    monkeypatch.setattr(energy_diagram, "build_energy_diagram", wrapped)
    result = CliRunner().invoke(
        root_cli, ["energy-diagram", "-i", "0", "-i", "5", "-o", str(tmp_path / "e.png")],
    )
    assert result.exit_code == 1
    assert isinstance(result.exception, SystemExit)
    _assert_guidance_only(result.output, pre.KIND_CXX_COMPILER)
    assert "\nError:" not in result.output


def test_a_usage_error_keeps_its_exit_status_and_message(tmp_path):
    from mlmm.cli import cli as root_cli

    result = CliRunner().invoke(
        root_cli, ["energy-diagram", "-i", "0", "-i", "5", "-o", str(tmp_path / "e.xyz")],
    )
    assert result.exit_code == 2
    assert "Unsupported output extension" in result.output
    assert "Missing system prerequisite" not in result.output


# ---------------------------------------------------------------------------
# The shared renderer of the commands (mlmm/cli/decorators.py)
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    ("inner", "kind"),
    [(INDUCTOR_C, pre.KIND_C_COMPILER), (INDUCTOR_CXX, pre.KIND_CXX_COMPILER)],
)
def test_render_cli_exception_prints_the_guidance(tmp_path, capsys, inner, kind):
    import json

    from mlmm.cli.decorators import render_cli_exception

    with pytest.raises(SystemExit) as stop:
        render_cli_exception(
            _wrapped(inner), label="TS optimization", out_dir=tmp_path, command="tsopt",
        )
    assert stop.value.code == 1
    _assert_guidance_only(capsys.readouterr().err, kind)
    assert json.loads((tmp_path / "result.json").read_text())["error_type"] == "InductorError"


def test_render_cli_exception_still_shows_the_traceback_of_other_errors(capsys):
    from mlmm.cli.decorators import render_cli_exception

    with pytest.raises(SystemExit) as stop:
        render_cli_exception(RuntimeError("an internal bug"), label="optimization")
    assert stop.value.code == 1
    err = capsys.readouterr().err
    assert "Unhandled error during optimization" in err and "an internal bug" in err
    assert "Missing system prerequisite" not in err


def test_a_stage_wrapped_by_all_prints_the_guidance_once(monkeypatch, tmp_path):
    """A stage prints the guidance and exits; all turns that exit into its own error line."""
    from mlmm.cli import cli as root_cli
    from mlmm.cli.decorators import render_cli_exception
    from mlmm.io import energy_diagram

    def stage(**kwargs):
        try:
            try:
                raise _wrapped(INDUCTOR_C)
            except Exception as exc:
                render_cli_exception(exc, label="path optimization")
        except SystemExit as stop:
            raise click.ClickException(f"[all] path-opt exit code {stop.code}.")

    monkeypatch.setattr(energy_diagram, "build_energy_diagram", stage)
    result = CliRunner().invoke(
        root_cli, ["energy-diagram", "-i", "0", "-i", "5", "-o", str(tmp_path / "e.png")],
    )
    assert result.exit_code == 1
    assert result.output.count("Missing system prerequisite") == 1, result.output
    assert "Error: [all] path-opt exit code 1." in result.output


def test_path_opt_dmf_without_a_compiler_keeps_the_error_record(monkeypatch, tmp_path):
    """The guidance exit still leaves the DMF failure in result.json."""
    import json
    from types import SimpleNamespace

    from mlmm.cli import cli as root_cli
    from mlmm.workflows import path_opt

    toy = Path(__file__).resolve().parents[1] / "examples" / "toy_system"
    if not (toy / "p_toy.parm7").exists():
        pytest.skip("toy_system example inputs not present")

    def no_compiler(*args, **kwargs):
        raise RuntimeError(TRITON_C)

    monkeypatch.setattr(path_opt, "mlmm", lambda **kwargs: SimpleNamespace(core=None))
    monkeypatch.setattr(path_opt, "_run_dmf_mep", no_compiler)
    out_dir = tmp_path / "out"
    result = CliRunner().invoke(root_cli, [
        "path-opt", "-i", str(toy / "r_toy.pdb"), str(toy / "p_toy.pdb"),
        "--parm", str(toy / "p_toy.parm7"), "-q", "-1", "-m", "1",
        "--mep-mode", "dmf", "--no-preopt", "--out-dir", str(out_dir),
    ])
    assert result.exit_code == 1
    _assert_guidance_only(result.output, pre.KIND_C_COMPILER)
    record = json.loads((out_dir / "result.json").read_text(encoding="utf-8"))
    assert (record["error_type"], record["error_label"]) == ("RuntimeError", "DMFError")
