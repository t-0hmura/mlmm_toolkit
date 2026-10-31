"""Guidance for system prerequisites that pip and conda do not install.

UMA builds code at run time (a C compiler on a GPU, a C++ compiler on a CPU) and PNG
export starts Chrome (which loads operating-system libraries). When one is missing,
torch, Triton, or Kaleido stop with an error that does not say what to install. This
module recognizes those errors and returns one fixed message: what is missing, how to
install it, and the docs page. It only describes the failure; it never retries, falls
back, or changes an exit status.
"""

from __future__ import annotations

import re
import textwrap
from typing import Iterator, Optional

import click

KIND_C_COMPILER = "c-compiler"
KIND_CXX_COMPILER = "cxx-compiler"
KIND_CHROME = "chrome"

DOCS_URL = "https://t-0hmura.github.io/mlmm_toolkit/installation.html#system-prerequisites"

# Messages of the third-party errors: Triton (GPU), Inductor (CPU), Kaleido (Chrome).
_CHROME_ERROR = re.compile(r"BrowserFailedError|BrowserDepsError|ChromeNotFoundError")
_PATTERNS = (
    (KIND_C_COMPILER, re.compile(r"Failed to find C compiler")),
    (KIND_CXX_COMPILER, re.compile(r"No working C\+\+ compiler found|InvalidCxxCompiler")),
    (KIND_CHROME, _CHROME_ERROR),
)

_ORIGINAL_ERROR_MAX = 300

# The operating-system packages of Chrome for Testing (its .deb and .rpm dependencies). libcurl is left out of the
# dnf list: RHEL 9 images carry libcurl-minimal, which provides libcurl.so.4 (the only curl library Chrome needs)
# and conflicts with libcurl.
_CHROME_APT = (
    "ca-certificates fonts-liberation libasound2t64 libatk-bridge2.0-0 libatk1.0-0 libatspi2.0-0 "
    "libcairo2 libcups2 libcurl4 libdbus-1-3 libexpat1 libgbm1 libglib2.0-0 libgtk-3-0 "
    "libnspr4 libnss3 libpango-1.0-0 libudev1 libvulkan1 libx11-6 libxcb1 libxcomposite1 "
    "libxdamage1 libxext6 libxfixes3 libxkbcommon0 libxrandr2 wget xdg-utils"
)
_CHROME_DNF = (
    "alsa-lib at-spi2-atk at-spi2-core atk ca-certificates cairo cups-libs dbus-libs expat "
    "glib2 gtk3 libX11 libXcomposite libXdamage libXext libXfixes libXrandr libgcc "
    "liberation-fonts libxcb libxkbcommon mesa-libgbm nspr nss nss-util pango systemd-libs "
    "vulkan-loader wget xdg-utils"
)


def _command(words: str, indent: int) -> str:
    """One shell command wrapped with line continuations, every line indented by ``indent``."""
    lines = textwrap.wrap(words, width=92 - indent, break_long_words=False, break_on_hyphens=False)
    return " \\\n".join(" " * indent + ("  " if n else "") + line for n, line in enumerate(lines))


_MESSAGES = {
    KIND_C_COMPILER: (
        "Missing system prerequisite: C compiler (gcc).\n"
        "  UMA on a GPU builds GPU code at run time through PyTorch (Triton), and no C compiler was found.\n"
        "  Install one, then run the command again:\n"
        "    Ubuntu/Debian:  sudo apt-get update && sudo apt-get install -y gcc\n"
        "    Rocky/RHEL:     sudo dnf install -y gcc\n"
        "    No sudo (HPC):  module load gcc   or   conda install -c conda-forge c-compiler"
    ),
    KIND_CXX_COMPILER: (
        "Missing system prerequisite: C++ compiler (g++).\n"
        "  UMA on a CPU builds its kernels at run time through PyTorch (Inductor), and no C++ compiler was found.\n"
        "  Install one, then run the command again:\n"
        "    Ubuntu/Debian:  sudo apt-get update && sudo apt-get install -y g++\n"
        "    Rocky/RHEL:     sudo dnf install -y gcc-c++\n"
        "    No sudo (HPC):  module load gcc   or   conda install -c conda-forge cxx-compiler"
    ),
    KIND_CHROME: "\n".join((
        "Missing system prerequisite: Chrome for PNG export.",
        "  PNG figures are drawn by headless Chrome (Plotly/Kaleido). Chrome is not installed, or it",
        "  cannot start because an operating-system library is missing.",
        "  1. Install Chrome (needs internet access):  plotly_get_chrome -y",
        "  2. Install the libraries Chrome loads:",
        "     Ubuntu/Debian (Ubuntu 22.04 and Debian 12: libasound2 instead of libasound2t64):",
        "       sudo apt-get update",
        _command("sudo apt-get install -y " + _CHROME_APT, 7),
        "     Rocky/RHEL:",
        _command("sudo dnf install -y " + _CHROME_DNF, 7),
        "     No sudo (HPC):  ask the administrator to install them.",
    )),
}


_BRIEF = {
    KIND_C_COMPILER: "Install gcc",
    KIND_CXX_COMPILER: "Install g++",
    KIND_CHROME: "Run plotly_get_chrome -y, then install the operating-system libraries of Chrome",
}


class MissingPrerequisiteError(RuntimeError):
    """A system prerequisite is missing. Its text is the one-line form, because commands such as ``all``
    print a skipped figure's error inline; the command that stops on it prints :func:`guidance`."""

    def __init__(self, kind: str, original: str = "") -> None:
        self.kind = kind
        self.original = original
        super().__init__(brief(kind, original))

    def __reduce__(self):
        return (type(self), (self.kind, self.original))


def _one_line(text: str) -> str:
    text = " ".join(text.split())
    return text if len(text) <= _ORIGINAL_ERROR_MAX else text[: _ORIGINAL_ERROR_MAX - 3] + "..."


def guidance(kind: str, original: str = "") -> str:
    """The fixed message for ``kind``, then the docs page and the original error line."""
    lines = [_MESSAGES[kind], f"  Details: {DOCS_URL}"]
    if original.strip():
        lines.append(f"  Original error: {_one_line(original)}")
    return "\n".join(lines)


def brief(kind: str, original: str = "") -> str:
    """The message of ``kind`` on one line: what is missing, the docs page, and the original error."""
    text = f"{_MESSAGES[kind].split(chr(10), 1)[0]} {_BRIEF[kind]}: {DOCS_URL}"
    return f"{text} ({_one_line(original)})" if original.strip() else text


def _chain(exc: BaseException) -> Iterator[BaseException]:
    seen = set()
    while exc is not None and id(exc) not in seen:
        seen.add(id(exc))
        yield exc
        following = exc.__cause__
        if following is None and not exc.__suppress_context__:
            following = exc.__context__
        exc = following


def is_chrome_failure(detail: str) -> bool:
    """Whether the text of a Plotly/Kaleido export error says Chrome is missing or cannot start."""
    return _CHROME_ERROR.search(detail) is not None


def prerequisite_message(exc: BaseException) -> Optional[str]:
    """The guidance for a known missing-prerequisite failure, or ``None`` for any other error
    (and for one whose guidance an inner handler has already printed). The original error line
    quotes the innermost matching error, which is the one torch, Triton, or Kaleido raised."""
    found = None
    for err in _chain(exc):
        if getattr(err, "prerequisite_reported", False):
            return None
        if isinstance(err, MissingPrerequisiteError):
            return guidance(err.kind, err.original)
        text = f"{type(err).__name__}: {err}"
        for kind, pattern in _PATTERNS:
            if pattern.search(text):
                found = guidance(kind, text.splitlines()[0])
                break
    return found


def exit_if_missing_prerequisite(exc: BaseException) -> None:
    """Print the guidance and exit with status 1 for a known missing prerequisite; else return."""
    message = prerequisite_message(exc)
    if message is None:
        return
    click.echo(message, err=True)
    stop = SystemExit(1)
    stop.prerequisite_reported = True  # a handler that wraps this exit must not print the guidance again
    raise stop
