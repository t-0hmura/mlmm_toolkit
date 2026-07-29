"""Shared configuration, structure-I/O, reporting, and plotting utilities."""

import ast
import logging
import math
import os
import re
import sys
import time
import tempfile
from collections import Counter
from collections.abc import Iterable as _Iterable, Mapping, Sequence as _Sequence
from dataclasses import dataclass, field
from numbers import Real, Integral
from pathlib import Path
from typing import Any, Callable, Dict, Optional, Sequence, List, Tuple

import click
import numpy as np
import yaml
import plotly.graph_objs as go

from pysisyphus.helpers import geom_loader
from pysisyphus.constants import ANG2BOHR

from mlmm.domain.add_elem_info import guess_element
from mlmm.core.output import _TAG_AWARE_MARKER, emit
from mlmm.core.result_commit import commit_payloads
from mlmm.io.structure_formats import (
    CIF_SUFFIXES,
    CoordinateTemplate,
    cleanup_normalized_structure,
    coordinate_template_for,
    is_cif_path,
    normalize_structure_to_pdb,
    pdb_requires_normalization,
    register_coordinate_template,
    render_pdb_coordinate_frames,
    render_mmcif_frames,
    unregister_coordinate_template,
    validate_coordinate_template_symbols,
    write_pdb_as_mmcif,
)

logger = logging.getLogger(__name__)


# CLI verbosity state (set by the top-level --verbose callback in cli/app.py).
# `pretty_block` and any other config-echo helpers consult `is_verbose()` so
# that running without `-v` keeps stdout focused on milestones / errors /
# user-set parameters, and `-v` (any positive count) restores the full
# config dump for debugging.
#
# The console output flows through two levels, decided in the single
# `_patch_click_echo` chokepoint:
#   default        : NARRATIVE only — stage banners, per-stage one-line status
#                    and the final summary. The bulk of the per-stage chatter
#                    (detail) is dropped. Warnings/errors (err=True) and
#                    explicit `force=True` deliverables always pass through.
#   -v/--verbose   : every line, plus pretty_block dumps and DEBUG logging.
# Narrative lines are flagged with `narrative=True` on `click.echo`; the patched
# echo strips that flag before delegating. `emit()` below is the helper used to
# tag those lines.
_VERBOSE_LEVEL: int = 0

# The narrative/detail gate stays disabled until the real CLI entry point calls
# `set_console_gating(True)` (from the -v group callback). Until then every
# echo prints (the patch only rewrites paths) so importing this module and
# calling `click.echo` from library code or tests is not silently muted. Once
# enabled, -v applies globally; default-verbosity narrative suppression is
# scoped to the `all` pipeline (see _PIPELINE_MODE / child mode).
_GATE_ACTIVE: bool = False

# True while the `all` command runs. Default-verbosity suppression of DETAIL
# lines is limited to the pipeline (the `all` parent and its in-proc child
# stages); a standalone leaf/report command (whose stdout is the deliverable,
# e.g. bond-summary / energy-diagram) keeps full output at default. -v is
# unaffected and applies everywhere.
_PIPELINE_MODE: bool = False

# Child-invocation state: True when this process is dispatching one of its
# subcommands through `_run_cli_main` (in-proc), so the child's group-callback
# banner / `[calc] Resolved device:` echo can be skipped to avoid repeating
# the same line for each in-process child stage.
_CHILD_MODE: bool = False


def set_verbose_level(level: int) -> None:
    """Record the CLI --verbose level (0=silent .. 3=full) for gating."""
    global _VERBOSE_LEVEL
    _VERBOSE_LEVEL = max(0, min(3, int(level)))


def verbose_level() -> int:
    """Current console verbosity: 0=silent, 1=milestones, 2=default(+detail), 3=full."""
    return _VERBOSE_LEVEL


def is_verbose() -> bool:
    """True iff verbosity reaches the detail tier (level >= 2). Back-compat shim
    for call sites that gate optimizer/SCF detail (cycle tables, PySCF logger)."""
    return _VERBOSE_LEVEL >= 2


def set_console_gating(value: bool) -> None:
    """Turn the narrative/detail console gate on or off.

    Invoked once per run from the real CLI group callbacks, so that only a
    genuine ``mlmm`` invocation suppresses detail output; library/test callers
    that import this module keep their echo output intact.
    """
    global _GATE_ACTIVE
    _GATE_ACTIVE = bool(value)


def is_console_gating() -> bool:
    """True iff the narrative/detail console gate is currently engaged."""
    return _GATE_ACTIVE


def set_pipeline_mode(value: bool) -> None:
    """Mark the `all` pipeline as running (scopes default-suppression)."""
    global _PIPELINE_MODE
    _PIPELINE_MODE = bool(value)


def is_pipeline_mode() -> bool:
    """True iff the `all` pipeline (parent or child stage) is running."""
    return _PIPELINE_MODE or _CHILD_MODE


def set_child_mode(value: bool) -> None:
    """Toggle child-invocation mode for in-proc subcommand dispatch."""
    global _CHILD_MODE
    _CHILD_MODE = bool(value)


def is_child_mode() -> bool:
    """True iff we are dispatching one of our own subcommands in-process."""
    return _CHILD_MODE


def echo_run_summary(items: Dict[str, Any]) -> None:
    """Echo a compact `[key] value` run summary, then a blank line.

    Skipped in child mode (`all` already printed its own summary) and when
    items is empty. Each subcommand calls this at its entry point so a
    default-verbosity run still surfaces the input file, backend, opt mode,
    and output dir without dumping the full per-stage config block (which
    only fires under `-v`).
    """
    if is_child_mode() or not items:
        return
    for key, value in items.items():
        if value is None or value == "":
            continue
        emit(f"[{key}] {value}", narrative=True)
    emit("", narrative=False)


def ensure_dir(path: Path) -> None:
    """Create a directory (parents ok); noop if it already exists."""
    path.mkdir(parents=True, exist_ok=True)


def echo_resolved_device() -> None:
    """Print the resolved torch device (cuda/cpu) as a cosmetic CLI breadcrumb.

    Best-effort: silently no-ops on torch-import or CUDA-probe failure
    (the echo is informational only; a missing torch is a separate
    error path that fires elsewhere when the calculator is built).
    """
    # Skip in child mode: the parent `mlmm all` already printed this once;
    # reprinting at every stage entry adds 4-8 identical lines per run.
    if is_child_mode():
        return
    try:
        import torch as _torch
        _resolved_dev = "cuda" if _torch.cuda.is_available() else "cpu"
    except (ImportError, AttributeError, RuntimeError):
        return
    import click as _click
    # Device breadcrumb is a level-3 (debug) line: untagged so it only shows at -v 3.
    _click.echo(f"[calc] Resolved device: {_resolved_dev}")


def optimizer_cycle_count(optimizer: Any) -> Optional[int]:
    """Return the number of optimizer cycles spent, if the optimizer exposes it."""
    cur_cycle = getattr(optimizer, "cur_cycle", None)
    if cur_cycle is None:
        return None
    try:
        return max(int(cur_cycle) + 1, 0)
    except (TypeError, ValueError):
        return None


def optimizer_terminal_status(optimizer: Any) -> str:
    """Map a pysisyphus optimizer (or a product-local runner) terminal state to
    the public status vocabulary.

    Returns ``"stalled"`` for an energy-plateau outcome,
    ``"converged"`` for a genuine stationary point, and ``"not_converged"``
    otherwise.  ``stalled`` takes precedence so a plateau is never reported as
    converged; legacy callers that only read ``is_converged`` still see
    ``False`` for a stall.
    """
    status = getattr(optimizer, "termination_status", None)
    if status in ("stalled", "converged", "not_converged"):
        return status
    if getattr(optimizer, "is_stalled", False):
        return "stalled"
    return "converged" if getattr(optimizer, "is_converged", False) else "not_converged"


def finalize_microiter_macro_convergence(
    macro_optimizer: Any,
    *,
    macro_converged: bool,
    latest_micro_stalled: bool,
    latest_micro_stop_reason: str = "",
) -> bool:
    """Fold a stalled latest micro (MM) relaxation into the macro terminal state.

    A stalled final micro (MM) relaxation is a real energy plateau, so it must
    never read as clean macro convergence. Surface it as an
    energy-plateau stall on ``macro_optimizer`` (carrying its reason) whether the
    macro would otherwise converge (a demotion) OR merely ran out of macro
    cycles -- otherwise a stalled final MM relaxation on a non-converged macro is
    lost as a reasonless ``not_converged``.  A macro that already stalled or
    requested its own (more specific) stop keeps that reason.  Returns the
    terminal macro-convergence flag, which is never ``True`` once a stall is
    surfaced.
    """
    macro_conv = bool(
        getattr(macro_optimizer, "is_converged", False) or macro_converged
    )
    # Write the terminal verdict back onto the optimizer.  The microiteration
    # driver calls ``check_convergence()`` directly instead of ``run()``, and
    # ``check_convergence`` only *returns* its verdict -- it never assigns
    # ``self.is_converged`` (that happens inside ``run``).  So a macro loop that
    # converged still carries ``is_converged=False``, and every consumer reads
    # the attribute, not this return value: ``optimizer_terminal_status`` and
    # ``OptimizerOutcome.from_optimizer`` both do ``getattr(optimizer,
    # "is_converged", False)``.  Without this write-back a converged TS with an
    # exact-PHVA-validated n_imag=1 saddle is reported ``not_converged`` and the
    # ``all`` pipeline refuses to start its IRC.
    if not latest_micro_stalled:
        macro_optimizer.is_converged = macro_conv
        return macro_conv
    # ``request_stall`` also sets ``stop_requested``; a macro that already
    # stalled or made a clean ``request_stop`` therefore short-circuits here and
    # keeps its own reason rather than being overwritten by the micro reason.
    if getattr(macro_optimizer, "stop_requested", False):
        macro_optimizer.is_converged = macro_conv
        return macro_conv
    # ``request_stall`` sets ``is_converged = False`` itself.
    macro_optimizer.request_stall(
        latest_micro_stop_reason
        or "energy plateau in the latest micro (MM) relaxation"
    )
    return False


def emit_optimizer_terminal_status(
    label: str,
    *,
    converged: Optional[bool],
    cycles: Optional[int],
    max_cycles: Optional[int],
    stalled: bool = False,
    stop_reason: Optional[str] = None,
) -> None:
    """Emit a consistent optimizer terminal status at detail verbosity.

    ``stalled`` renders the energy-plateau outcome and takes
    precedence over the convergence/max-cycle branches so a stalled run is
    never printed as ``Converged!``.
    """
    prefix = f"[{label}]"
    if stalled:
        if stop_reason:
            emit(f"{prefix} Stalled (energy plateau; not converged): {stop_reason}", detail=True)
        else:
            emit(f"{prefix} Stalled (energy plateau; not converged).", detail=True)
    elif converged is True:
        emit(f"{prefix} Converged!", detail=True)
    elif cycles is not None and max_cycles is not None and cycles >= max_cycles:
        emit(f"{prefix} Reached max cycles ({cycles}/{max_cycles}).", detail=True)
    elif converged is False:
        if cycles is None:
            emit(f"{prefix} Stopped without convergence.", detail=True)
        else:
            emit(f"{prefix} Stopped without convergence (cycles={cycles}).", detail=True)
    elif cycles is not None:
        emit(f"{prefix} Finished (cycles={cycles}).", detail=True)

    if cycles is not None:
        emit(f"{prefix} Total cycles: {cycles}", detail=True)


def _parse_freeze_atoms(arg: Optional[str]) -> List[int]:
    """Parse comma-separated 1-based indices (e.g., ``"1,3,5"``) into a sorted 0-based list.

    Canonical home for the ``--freeze-atoms`` CLI parser. ``workflows/opt.py``
    re-exports this symbol for backward compatibility so every subcommand
    (``sp``, ``opt``, ``tsopt``, ``path-opt``, ``freq``, ``irc``, ``scan*``,
    ``path-search``, ``dft``) can import it from a single place.
    """
    if arg is None:
        return []

    items = [chunk.strip() for chunk in str(arg).split(",")]
    indices: List[int] = []
    for idx, chunk in enumerate(items, start=1):
        if not chunk:
            continue
        try:
            value = int(chunk)
        except ValueError as exc:
            raise click.BadParameter(
                f"Invalid integer in --freeze-atoms entry #{idx}: '{chunk}'"
            ) from exc
        if value <= 0:
            raise click.BadParameter(
                f"--freeze-atoms expects 1-based positive indices; got {value}"
            )
        indices.append(value - 1)
    return sorted(set(indices))


def read_xyz_as_blocks(path: Path, *, strict: bool = False) -> List[List[str]]:
    """Read an XYZ-style trajectory into blocks of lines.

    When *strict* is True, malformed headers or truncated frames raise a ClickException.
    """
    try:
        lines = path.read_text(encoding="utf-8").splitlines()
    except Exception as e:
        import click
        raise click.ClickException(f"Failed to read {path}: {e}")

    blocks: List[List[str]] = []
    i = 0
    while i < len(lines):
        if not lines[i].strip():
            i += 1
            continue
        try:
            n_atoms = int(lines[i].strip().split()[0])
        except Exception:
            if strict:
                import click
                raise click.ClickException(f"[xyz] Malformed XYZ/TRJ header at line {i+1} of {path}")
            break
        end = i + n_atoms + 2
        if end > len(lines):
            if strict:
                import click
                raise click.ClickException(f"[xyz] Incomplete XYZ frame at line {i+1} of {path}")
            break
        blocks.append(lines[i:end])
        i = end
    return blocks


def parse_xyz_block(
    block: Sequence[str],
    *,
    path: Path,
    frame_idx: int,
) -> Tuple[List[str], "np.ndarray"]:
    """Parse a single XYZ frame block into (elements, coords_angstrom)."""
    import click

    if not block:
        raise click.ClickException(f"[xyz] Empty XYZ frame in {path}")
    try:
        nat = int(block[0].strip().split()[0])
    except Exception:
        raise click.ClickException(
            f"[xyz] Malformed XYZ/TRJ header in frame {frame_idx} of {path}"
        )
    if len(block) < 2 + nat:
        raise click.ClickException(
            f"[xyz] Incomplete XYZ frame {frame_idx} in {path} (expected {nat} atoms)."
        )
    elems: List[str] = []
    coords: List[List[float]] = []
    for k in range(nat):
        parts = block[2 + k].split()
        if len(parts) < 4:
            raise click.ClickException(
                f"[xyz] Malformed atom line in frame {frame_idx} of {path}"
            )
        elems.append(parts[0])
        coords.append([float(parts[1]), float(parts[2]), float(parts[3])])
    return elems, np.array(coords, dtype=float)


def xyz_blocks_first_last(
    blocks: Sequence[Sequence[str]],
    *,
    path: Path,
) -> Tuple[List[str], "np.ndarray", "np.ndarray"]:
    """Return (elements, first_coords_ang, last_coords_ang) from pre-parsed XYZ blocks."""
    import click

    if not blocks:
        raise click.ClickException(f"[xyz] No frames found in {path}")
    first_elems, first_coords = parse_xyz_block(blocks[0], path=path, frame_idx=1)
    last_elems, last_coords = parse_xyz_block(blocks[-1], path=path, frame_idx=len(blocks))
    if first_elems != last_elems:
        raise click.ClickException(f"[xyz] Element list changed across frames in {path}")
    return first_elems, first_coords, last_coords


def read_xyz_first_last(trj_path: Path) -> Tuple[List[str], "np.ndarray", "np.ndarray"]:
    """Lightweight XYZ trajectory reader: return (elements, first_coords[Å], last_coords[Å])."""
    blocks = read_xyz_as_blocks(trj_path, strict=True)
    return xyz_blocks_first_last(blocks, path=trj_path)


def close_matplotlib_figures() -> None:
    """Best-effort cleanup for matplotlib figures to avoid open-figure warnings."""
    try:
        import matplotlib.pyplot as plt
        plt.close("all")
    except Exception:
        pass


def distance_A_from_coords(coords_bohr: "np.ndarray", i: int, j: int) -> float:
    """Return interatomic distance in Å given coords in Bohr."""
    diff = coords_bohr[i] - coords_bohr[j]
    return float(np.linalg.norm(diff) / ANG2BOHR)


def distance_tag(value_A: float, *, digits: int = 2, pad: int = 3) -> str:
    """Format a distance in Å as a zero-padded integer tag (default: ×10^2)."""
    scale = 10 ** digits
    return f"{int(round(value_A * scale)):0{pad}d}"


def values_from_bounds(low: float, high: float, h: float) -> "np.ndarray":
    """Return evenly spaced values from low→high with step cap h (inclusive)."""
    if h <= 0.0:
        raise click.BadParameter("--max-step-size must be > 0.")
    delta = abs(high - low)
    if delta < 1e-12:
        return np.array([low], dtype=float)
    N = int(math.ceil(delta / h))
    return np.linspace(low, high, N + 1, dtype=float)


def geom_from_xyz_string(
    xyz_text: str,
    *,
    coord_type: str,
    freeze_atoms: Optional[Sequence[int]] = None,
) -> Any:
    """Load a pysisyphus Geometry from an XYZ text string (tempfile-backed)."""
    s = xyz_text if xyz_text.endswith("\n") else (xyz_text + "\n")
    freeze_atoms = list(freeze_atoms) if freeze_atoms is not None else []
    tmp = tempfile.NamedTemporaryFile("w+", suffix=".xyz", delete=False)
    try:
        tmp.write(s)
        tmp.flush()
        tmp.close()

        g = geom_loader(
            Path(tmp.name),
            coord_type=coord_type,
            freeze_atoms=freeze_atoms,
        )
        try:
            g.freeze_atoms = np.array(sorted(set(map(int, freeze_atoms))), dtype=int)
        except Exception:
            click.echo(
                "[geom] WARNING: Failed to attach freeze_atoms to geometry.",
                err=True,
            )
        return g
    finally:
        try:
            os.unlink(tmp.name)
        except Exception:
            logger.debug("Failed to unlink temp file %s", tmp.name, exc_info=True)


def append_xyz_trajectory(dst_path: Path, src_path: Path, *, reset: bool = False) -> bool:
    """Append an XYZ trajectory segment to a concatenated trajectory file."""
    if not src_path.exists():
        return False
    mode = "w" if reset else "a"
    with src_path.open("r", encoding="utf-8") as src, dst_path.open(mode, encoding="utf-8") as dst:
        while True:
            chunk = src.read(1024 * 1024)
            if not chunk:
                break
            dst.write(chunk)
    return True


def snapshot_geometry(geom: Any, *, coord_type_default: str) -> Any:
    """Create an independent pysisyphus Geometry snapshot from the given Geometry."""
    s = geom.as_xyz()
    return geom_from_xyz_string(
        s,
        coord_type=getattr(geom, "coord_type", coord_type_default),
        freeze_atoms=getattr(geom, "freeze_atoms", []),
    )


def unbiased_energy_hartree(geom, base_calc) -> float:
    """Evaluate the underlying ML/MM energy (Hartree) without harmonic bias."""
    # ``geom.coords`` is an internal-coordinate vector for redund/dlc/tric;
    # calculator calls always require Cartesian coordinates in bohr.
    coords_bohr = np.asarray(geom.coords3d)
    elems = getattr(geom, "atoms", None)
    if elems is None:
        return float("nan")
    try:
        return float(base_calc.get_energy(elems, coords_bohr)["energy"])
    except Exception as exc:
        click.echo(
            f"[energy] WARNING: bare ML/MM energy evaluation failed: {exc}",
            err=True,
        )
        return float("nan")


def calculator_provenance(calc_cfg: Mapping[str, Any]) -> Dict[str, Any]:
    """Return backend-neutral ML/MM calculator provenance for JSON outputs."""
    from mlmm.core.defaults import MLMM_CALC_KW

    backend = str(calc_cfg.get("backend") or MLMM_CALC_KW["backend"]).lower()
    model_keys = {
        "uma": "uma_model",
        "orb": "orb_model",
        "mace": "mace_model",
        "aimnet2": "aimnet2_model",
    }
    if backend == "custom":
        calc_file = calc_cfg.get("calc_file")
        factory = calc_cfg.get("calc_factory") or "get_calculator"
        model = f"{Path(calc_file).name}:{factory}" if calc_file else str(factory)
        precision = None
    else:
        key = model_keys.get(backend)
        model = calc_cfg.get(key) if key is not None else None
        if model is None and key is not None:
            model = MLMM_CALC_KW.get(key)
        precision_keys = {
            "uma": "uma_precision",
            "orb": "orb_precision",
            "mace": "mace_dtype",
        }
        precision_key = precision_keys.get(backend)
        if precision_key is None:
            precision = "fp32" if backend == "aimnet2" else None
        else:
            precision = calc_cfg.get(precision_key)
            if precision is None:
                precision = MLMM_CALC_KW.get(precision_key)

        token = "" if precision is None else str(precision).strip().lower()
        if token in {"fp64", "float64", "double", "highest"}:
            precision = "fp64"
        elif token in {"fp32", "float32", "float32-high", "float32-highest", "single"}:
            precision = "fp32"
        else:
            precision = token or None

    return {
        "mlip_backend": backend,
        "mlip_model": None if model is None else str(model),
        "mlip_precision": None if precision is None else str(precision),
        "mm_backend": str(calc_cfg.get("mm_backend") or MLMM_CALC_KW["mm_backend"]),
        "link_atom_method": str(
            calc_cfg.get("link_atom_method") or MLMM_CALC_KW["link_atom_method"]
        ),
        "use_cmap": bool(calc_cfg.get("use_cmap", MLMM_CALC_KW["use_cmap"])),
    }


def calculator_run_label(calc_cfg: Mapping[str, Any]) -> str:
    """Format backend, model, and effective precision for concise run headers."""
    provenance = calculator_provenance(calc_cfg)
    backend = provenance["mlip_backend"]
    model = provenance["mlip_model"]
    precision = provenance["mlip_precision"]
    details = [str(value) for value in (model, precision) if value not in (None, "")]
    return f"{backend} ({', '.join(details)})" if details else str(backend)


def pretty_block(title: str, content: Dict[str, Any], *, force: bool = False) -> str:
    """Return a YAML-formatted block with an underlined title.

    Returns an empty string below verbosity level 3 so that the default
    CLI output stays focused on milestones and user-set parameters; the
    full config dump is restored under `-v 3` for debugging.

    ``force=True`` bypasses that gate. Use it for output the user asked for
    explicitly (``--print-parsed``, ``--show-config``): a flag whose whole purpose is
    to print something must not render nothing at the default verbosity.
    """
    if not force and verbose_level() < 3:
        return ""
    if not content:
        return ""  # suppress empty blocks entirely
    if _base_dir is not None:
        content = _shorten_paths(content)
    body = yaml.safe_dump(_to_yaml_safe(content), sort_keys=False, allow_unicode=True).strip()
    return f"\n{title}\n" + "-" * len(title) + "\n" + body


# Module-level base directory for relative path display.
_base_dir: Path | None = None
_original_click_echo = None

# Raw-stdout verbosity tap for bundled optimizers whose per-cycle tables and
# summaries are written via sys.stdout (bypassing the click.echo gate). At -v 1
# only convergence verdicts are kept. At -v 2 optimizer tables are kept but
# high-volume DLC/trust-radius chatter and IPOPT metric rows that grep as
# "error" are hidden. At -v 3 raw optimizer stdout is passed through unchanged.
_PYSIS_L1_ALLOW = re.compile(
    r"^(?:Converged!|Final summary:"
    r"|max\(forces,\s*\w+\):|rms\(forces,\s*\w+\):|energy:\s"
    r"|Path with \d+ moving images\.|Number of cycles exceeded!"
    r"|Operator indicated convergence!|Insignificant coordinate change"
    r"|Energy plateau detected|Wrote final geometr)"
)
_PYSIS_L2_DENY = re.compile(
    r"^(?:\d+\s+(?:[-\d.]|nan)"          # compact table data row (cycle 0 col = nan*)
    r"|cycle\s+\S.*energy"               # compact table header
    r"|[-=]{5,}\s*$"                      # separator rule
    r"|If not specified otherwise, all quantities"
    r"|Spent\s+[\d.]+\s+s\s+preparing"
    r"|Convergence thresholds|'Superscript"
    r"|max\(\|force\|\)\s*<=|rms\(force\)\s*<="
    r"|max\(\|step\|\)\s*<=|rms\(step\)\s*<="
    r"|Rebuilt internal coordinates|Interfragment distances increased"
    r"|Dumped latest coordinates|String=|Overall NLP error\.+:)"
)
_PYSIS_V2_DENY = re.compile(
    r"^(?:Rebuilt internal coordinates|Interfragment distances increased"
    r"|Dumped latest coordinates"
    r"|Overall NLP error\.+:"
    r"|Unexpected energy increase"
    r"|Current trust radius:|Decreasing trust radius\.|Increasing trust radius\."
    r"|Keeping current trust radius|Updated trust radius:)"
)


def _pysis_stdout_visible(stripped: str) -> bool:
    """Whether a raw-stdout optimizer line is visible at the current level."""
    if _VERBOSE_LEVEL <= 0:
        return False                       # -v 0: silent
    if _VERBOSE_LEVEL >= 3:
        return True                        # -v 3: full raw optimizer output
    if _VERBOSE_LEVEL >= 2:
        return not _PYSIS_V2_DENY.match(stripped)
    if _PYSIS_L1_ALLOW.match(stripped):    # -v 1: keep the convergence verdict
        return True
    return not _PYSIS_L2_DENY.match(stripped)  # drop only the table noise


def _patch_click_echo() -> None:
    """Monkey-patch click.echo to shorten absolute paths in output."""
    import click as _click
    global _original_click_echo
    if _original_click_echo is not None:
        return  # already patched
    _original_click_echo = _click.echo

    _last_was_blank = [False]
    _last_visible_line = [""]
    _raw_path_echo_depth = [0]

    def _is_hessian_status_line(line: str) -> bool:
        return (
            line.startswith("[hessian]")
            or line.startswith("[HessianTiming]")
            or line.startswith("[HessianVRAM]")
        )

    def _wants_blank_before(first_visible: str, prev_visible: str) -> bool:
        if not first_visible:
            return False
        starts_hessian_status = _is_hessian_status_line(first_visible)
        return (
            first_visible.startswith("======")
            or first_visible.startswith("[time]")
            or first_visible.startswith("[stage]")
            or (
                starts_hessian_status
                and not _is_hessian_status_line(prev_visible)
            )
            or (
                first_visible.startswith("[Imaginary modes]")
                and not prev_visible.startswith("[Imaginary modes]")
            )
            or first_visible.startswith("Convergence thresholds (non mass-weighted gradient)")
            or first_visible.startswith("IRC steps exceeded.")
            or first_visible.startswith("Transition vector is mode")
            or first_visible.startswith("Wrote final geometry")
            or (
                _is_hessian_status_line(prev_visible)
                and not _is_hessian_status_line(first_visible)
            )
        )

    def _wants_blank_after(first_visible: str) -> bool:
        return first_visible.startswith("======") or first_visible.startswith("[stage]")

    def _patched_echo(message=None, **kwargs):
        # Pull the verbosity tags before handing off (click.echo would choke on
        # unknown kwargs). A line's required level: narrative -> 1 (milestone),
        # detail -> 2 (cycle tables / timing / deliverable paths), untagged -> 3
        # (config dumps / per-file path bullets / device breadcrumbs). Untagged
        # lines fall through to level 3, never to "dropped", so `-v 3`
        # reproduces every original line.
        narrative = bool(kwargs.pop("narrative", False))
        detail = bool(kwargs.pop("detail", False))
        # `force=True` bypasses the level gate entirely: it marks an explicitly
        # requested machine-readable deliverable (e.g. `--json` output) that must
        # always reach stdout (except at -v 0). Path-shortening still applies.
        force = bool(kwargs.pop("force", False))
        raw_path = bool(kwargs.pop("raw_path", False))
        is_err = bool(kwargs.get("err", False))
        is_blank = (message is None or (isinstance(message, str) and message.strip() == ""))
        # Level gate, active only after the real CLI calls set_console_gating():
        #   level 0 (-v 0): completely silent (even errors/forced output).
        #   err=True / force=True: shown at level >= 1.
        #   otherwise, inside the `all` pipeline a line prints iff the current
        #   level >= its required level (narrative 1 / detail 2 / untagged 3).
        #   Blank lines pass through (spacing); a standalone leaf/report command
        #   keeps full output (its stdout is the deliverable).
        if _GATE_ACTIVE:
            if _VERBOSE_LEVEL <= 0:
                return  # -v 0: completely silent
            if (
                not is_err
                and not force
                and not is_blank
                and is_pipeline_mode()
            ):
                required = 1 if narrative else (2 if detail else 3)
                if _VERBOSE_LEVEL < required:
                    return
        if not raw_path and message is not None and _base_dir is not None and isinstance(message, str):
            bd = str(_base_dir)
            if bd in message:
                message = message.replace(bd + "/", "").replace(bd, ".")
        if isinstance(message, str):
            first_visible = next((line.strip() for line in message.splitlines() if line.strip()), "")
            prev_visible = _last_visible_line[0]
            if not _last_was_blank[0] and not message.startswith("\n") and _wants_blank_before(first_visible, prev_visible):
                message = "\n" + message
            if first_visible and _wants_blank_after(first_visible) and not message.endswith("\n"):
                message += "\n"
        # Suppress consecutive blank lines
        if isinstance(message, str) and _last_was_blank[0] and message.startswith("\n"):
            message = message.lstrip("\n")
        if is_blank and _last_was_blank[0]:
            return
        ends_with_nl = isinstance(message, str) and message.endswith("\n")
        _last_was_blank[0] = is_blank or ends_with_nl
        if raw_path:
            _raw_path_echo_depth[0] += 1
        try:
            _original_click_echo(message, **kwargs)
            if isinstance(message, str):
                for line in reversed(message.splitlines()):
                    if line.strip():
                        _last_visible_line[0] = line.strip()
                        break
        finally:
            if raw_path:
                _raw_path_echo_depth[0] -= 1

    setattr(_patched_echo, _TAG_AWARE_MARKER, True)
    _click.echo = _patched_echo

    # Wrap sys.stdout so the bundled pysisyphus optimizer's raw output obeys
    # the verbosity level: silent at -v 0, verdict-only at -v 1, table with
    # high-volume chatter suppressed at -v 2, full raw output at -v 3.
    import sys as _sys
    _real_stdout = _sys.stdout

    class _FilteredStdout:
        def __init__(self, stream):
            self._stream = stream
            self._last_was_blank = False
            self._current_line_has_content = False
            self._has_visible_output = False
            self._suppress_next_nl = False

        @staticmethod
        def _starts_raw_section(stripped: str) -> bool:
            return (
                stripped.startswith("Spent ")
                or stripped.startswith("Path with ")
                or stripped.startswith("Final summary:")
                or stripped.startswith("Convergence thresholds")
                or stripped.startswith("IRC steps exceeded.")
                or stripped.startswith("Transition vector is mode")
                or stripped.startswith("Wrote final geometry")
                or stripped.startswith("======")
            )

        def write(self, s):
            original_len = len(s) if isinstance(s, str) else None
            if _GATE_ACTIVE and isinstance(s, str):
                stripped = s.strip()
                if not stripped:
                    # blank/whitespace: swallow entirely at -v 0, or when it
                    # trails a line we just suppressed (avoid a blank gap).
                    if _VERBOSE_LEVEL <= 0 or self._suppress_next_nl:
                        self._suppress_next_nl = False
                        return len(s)
                elif not _pysis_stdout_visible(stripped):
                    self._suppress_next_nl = True
                    return len(s)
                else:
                    self._suppress_next_nl = False
                    prev_visible = _last_visible_line[0]
                    visible_before = self._has_visible_output or bool(prev_visible)
                    needs_blank = (
                        self._starts_raw_section(stripped)
                        or _wants_blank_before(stripped, prev_visible)
                    )
                    if (
                        visible_before
                        and needs_blank
                        and not self._last_was_blank
                        and not _last_was_blank[0]
                    ):
                        self._stream.write("\n")
                        self._last_was_blank = True
                        _last_was_blank[0] = True
                        self._current_line_has_content = False
            if _raw_path_echo_depth[0] <= 0 and _base_dir is not None and isinstance(s, str):
                bd = str(_base_dir)
                if bd in s:
                    s = s.replace(bd + "/", "").replace(bd, ".")
            if isinstance(s, str):
                if self._last_was_blank and s.startswith("\n"):
                    s = s.lstrip("\n")
                    if not s:
                        return original_len
                s = re.sub(r"\n{3,}", "\n\n", s)
            if s == "\n":
                if self._current_line_has_content:
                    self._current_line_has_content = False
                    self._last_was_blank = False
                    self._has_visible_output = True
                    _last_was_blank[0] = False
                    return self._stream.write(s)
                if self._last_was_blank:
                    _last_was_blank[0] = True
                    return len(s)
                self._last_was_blank = True
            elif s.strip() == "":
                pass
            else:
                self._current_line_has_content = True
                self._last_was_blank = False
                self._has_visible_output = True
                if s.endswith("\n\n"):
                    self._current_line_has_content = False
                    self._last_was_blank = True
                elif s.endswith("\n"):
                    self._current_line_has_content = False
            _last_was_blank[0] = self._last_was_blank
            if isinstance(s, str):
                for line in reversed(s.splitlines()):
                    if line.strip():
                        _last_visible_line[0] = line.strip()
                        break
            return self._stream.write(s)

        def flush(self):
            self._stream.flush()

        def __getattr__(self, name):
            return getattr(self._stream, name)

    _sys.stdout = _FilteredStdout(_real_stdout)


def set_base_dir(path: Path | str | None) -> None:
    """Set the base directory for relative path display.

    Also monkey-patches ``click.echo`` so that any absolute path under
    *base_dir* is automatically shortened to a relative path in all
    CLI output.
    """
    global _base_dir
    _base_dir = Path(path).resolve() if path else None
    _patch_click_echo()


def rel_display(path: Path | str) -> str:
    """Return a display string for *path*, relative to the base dir when possible."""
    p = Path(path)
    if _base_dir is not None:
        try:
            return str(p.resolve().relative_to(_base_dir))
        except ValueError:
            pass
    return str(p)


def _shorten_paths(content: Dict[str, Any]) -> Dict[str, Any]:
    """Replace absolute path strings with relative paths in a config dict."""
    out: Dict[str, Any] = {}
    for k, v in content.items():
        if isinstance(v, str) and v.startswith("/") and ("/" in v[1:]):
            out[k] = rel_display(v)
        else:
            out[k] = v
    return out


def _to_yaml_safe(value: Any) -> Any:
    """Recursively convert NumPy values/containers into YAML-safe builtins."""
    if isinstance(value, np.generic):
        return value.item()
    if isinstance(value, np.ndarray):
        return [_to_yaml_safe(v) for v in value.tolist()]
    if isinstance(value, Mapping):
        out: Dict[Any, Any] = {}
        for k, v in value.items():
            nk = _to_yaml_safe(k)
            if isinstance(nk, (list, tuple, set, dict)):
                nk = str(nk)
            out[nk] = _to_yaml_safe(v)
        return out
    if isinstance(value, tuple):
        return [_to_yaml_safe(v) for v in value]
    if isinstance(value, list):
        return [_to_yaml_safe(v) for v in value]
    if isinstance(value, set):
        return [_to_yaml_safe(v) for v in sorted(value, key=lambda x: str(x))]
    return value


# Backend-specific key prefixes in MLMM_CALC_KW.
# Keys with these prefixes are only relevant when the corresponding backend is active.
_BACKEND_KEY_PREFIXES: Dict[str, tuple] = {
    "uma": ("uma_model", "uma_task_name", "uma_precision"),
    "orb": ("orb_model", "orb_precision"),
    "mace": ("mace_model", "mace_dtype"),
    "aimnet2": ("aimnet2_model",),
}


def filter_calc_for_echo(calc_cfg: Dict[str, Any]) -> Dict[str, Any]:
    """Remove backend-specific keys that are irrelevant for the active backend.

    Also hides xTB/embedcharge keys when embedcharge is disabled,
    and freeze_atoms (already shown in the geom block).
    """
    cfg = dict(calc_cfg)
    cfg.pop("freeze_atoms", None)
    active = cfg.get("backend", "uma")

    # Remove keys belonging to inactive ML backends
    for backend, keys in _BACKEND_KEY_PREFIXES.items():
        if backend != active:
            for k in keys:
                cfg.pop(k, None)

    # Hide xTB-specific keys when embedcharge is disabled
    if not cfg.get("embedcharge"):
        for k in list(cfg):
            if k.startswith("xtb_"):
                cfg.pop(k)
        cfg.pop("embedcharge_step", None)
        cfg.pop("embedcharge_cutoff", None)

    return cfg


def strip_inherited_keys(
    child_cfg: Dict[str, Any],
    base_cfg: Dict[str, Any],
    *,
    mode: str = "present",
) -> Dict[str, Any]:
    """Return child_cfg without inherited keys (for concise logs).

    Parameters
    ----------
    child_cfg : Dict[str, Any]
        The child configuration dictionary to trim.
    base_cfg : Dict[str, Any]
        The base configuration dictionary to compare against.
    mode : str
        - "present": Remove keys that exist in base_cfg regardless of value.
        - "same": Remove keys only when the value matches base_cfg.

    Returns
    -------
    Dict[str, Any]
        A new dictionary with inherited keys removed.
    """
    if mode not in {"present", "same"}:
        raise ValueError(f"Unknown strip_inherited_keys mode: {mode}")
    trimmed: Dict[str, Any] = {}
    for key, value in child_cfg.items():
        if key in base_cfg:
            if mode == "present":
                continue
            if base_cfg.get(key) == value:
                continue
        trimmed[key] = value
    return trimmed


def _summarize_atom_indices(items: Sequence[Any]) -> str:
    """Return a compact single-line summary for atom indices."""
    if not items:
        return ""

    count = len(items)
    if count <= 64:
        return f"{count} atoms [{','.join(map(str, items))}]"

    head = ",".join(map(str, items[:5]))
    tail = ",".join(map(str, items[-5:]))
    return f"{count} atoms [{head},...,{tail}]"


def format_freeze_atoms_for_echo(
    cfg: Dict[str, Any],
    *,
    key: str = "freeze_atoms",
) -> Dict[str, Any]:
    """
    Normalize freeze-atoms fields for concise CLI echo output.
    """
    g = dict(cfg)
    freeze_atoms = g.get(key)
    if freeze_atoms is None:
        return g

    if isinstance(freeze_atoms, str):
        return g

    try:
        items = list(freeze_atoms)
    except TypeError:
        return g

    # Display as 1-based (internal is 0-based)
    items_1based = [i + 1 for i in items]
    g[key] = _summarize_atom_indices(items_1based)
    return g


def format_elapsed(prefix: str, start_time: float, end_time: Optional[float] = None) -> str:
    """Return a formatted elapsed-time string with the provided ``prefix`` label."""
    finish = end_time if end_time is not None else time.perf_counter()
    elapsed = max(0.0, finish - start_time)
    hours, rem = divmod(elapsed, 3600)
    minutes, seconds = divmod(rem, 60)
    return f"{prefix}: {int(hours):02d}:{int(minutes):02d}:{seconds:06.3f}"


def normalize_freeze_atoms(raw: Any) -> List[int]:
    """Normalize freeze_atoms values (string/list/iterable) into a list of integers.

    Parameters
    ----------
    raw : Any
        Input value that can be a string (e.g., "1,2,3" or "1 2 3"),
        a list of integers, or any iterable of numeric values.

    Returns
    -------
    List[int]
        List of integer indices.

    Examples
    --------
    >>> normalize_freeze_atoms("1, 2, 3")
    [1, 2, 3]
    >>> normalize_freeze_atoms([1, 2, 3])
    [1, 2, 3]
    >>> normalize_freeze_atoms(None)
    []
    """
    import re

    if raw is None:
        return []
    if isinstance(raw, str):
        tokens = re.findall(r"-?\d+", raw)
        return [int(tok) for tok in tokens]
    try:
        items = list(raw)
    except TypeError:
        return []
    out: List[int] = []
    for item in items:
        try:
            out.append(int(item))
        except (TypeError, ValueError) as exc:
            # Never fall back to an empty list for an unparsable entry: that turned
            # `geom.freeze_atoms: [1, 2, three]` into a run with NOTHING frozen, silently
            # making a hard-freeze / PHVA result an unconstrained one.
            raise ValueError(
                f"freeze_atoms: cannot interpret {item!r} as an atom index"
            ) from exc
    return out


def merge_freeze_atom_indices(
    geom_cfg: Dict[str, Any],
    *indices: _Iterable[int],
) -> List[int]:
    """Merge one or more iterables of indices into ``geom_cfg['freeze_atoms']``.

    Existing entries are preserved, duplicates removed, and the result sorted.
    The updated list is returned.
    """
    merged: set[int] = set()
    base = geom_cfg.get("freeze_atoms", None)
    merged.update(normalize_freeze_atoms(base))
    for seq in indices:
        merged.update(normalize_freeze_atoms(seq))
    result = sorted(merged)
    geom_cfg["freeze_atoms"] = result
    return result


def parse_pdb_coords(pdb_path):
    """Parse ATOM/HETATM records from *pdb_path* and separate link hydrogen (HL) atoms.

    Returns:
        A tuple (others, lkhs) where:
            - others: list of tuples (index, x, y, z, line) for all atoms except the
              'HL' atom of residue 'LKH'. ``index`` is the 0-based position in the
              atom sequence as loaded from the *first* MODEL (or the full file if no
              MODEL records are present).
            - lkhs: list of tuples (x, y, z, line) for atoms where residue name is
              'LKH' and atom name is 'HL' in the same MODEL selection.

    Notes
    -----
        - Coordinates are read from standard PDB columns:
          X: columns 31-38, Y: 39-46, Z: 47-54 (1-based indexing).
        - If multiple MODEL blocks are present, only the first model is considered,
          matching typical geom_loader behavior.
    """
    with open(pdb_path, "r") as f:
        lines = f.readlines()

    others = []
    lkhs = []
    model_seen = False
    in_first_model = True
    atom_index = 0
    for line in lines:
        if line.startswith("MODEL"):
            if not model_seen:
                model_seen = True
                in_first_model = True
            else:
                in_first_model = False
            continue
        if line.startswith("ENDMDL"):
            if model_seen and in_first_model:
                break
            continue
        if model_seen and not in_first_model:
            continue
        if not (line.startswith("ATOM") or line.startswith("HETATM")):
            continue

        current_index = atom_index
        atom_index += 1

        name = line[12:16].strip()
        resname = line[17:20].strip()
        try:
            x = float(line[30:38])
            y = float(line[38:46])
            z = float(line[46:54])
        except ValueError:
            continue

        if resname == "LKH" and name == "HL":
            lkhs.append((x, y, z, line))
        else:
            others.append((current_index, x, y, z, line))
    return others, lkhs


def nearest_index(point, pool):
    """Find the nearest point in *pool* to *point* using Euclidean distance.

    Args:
        point: Tuple (x, y, z) representing the query coordinate.
        pool: Iterable of tuples (index, x, y, z, line) to search.

    Returns:
        A tuple (index, distance) where:
            - index is the 0-based index of the nearest entry in *pool* (or -1 if *pool* is empty).
            - distance is the Euclidean distance to that entry (``inf`` if *pool* is empty).
    """
    x, y, z = point
    best_i = -1
    best_d2 = float("inf")
    for atom_index, a, b, c, _ in pool:
        d2 = (a - x) ** 2 + (b - y) ** 2 + (c - z) ** 2
        if d2 < best_d2:
            best_d2 = d2
            best_i = atom_index
    return best_i, math.sqrt(best_d2)




def apply_layer_freeze_constraints(
    geom_cfg: Dict[str, Any],
    calc_cfg: Dict[str, Any],
    layer_info: Optional[Dict[str, Sequence[int]]],
    *,
    echo_fn: Optional[Callable[[str], None]] = None,
) -> List[int]:
    """Merge frozen-layer atoms into geometry/calculator freeze lists."""
    if echo_fn is not None:
        echo_fn("")  # blank line after layer detection summary
    base_freeze = normalize_freeze_atoms(geom_cfg.get("freeze_atoms"))
    frozen_from_layer = normalize_freeze_atoms((layer_info or {}).get("frozen_indices", []))

    if frozen_from_layer:
        before = set(base_freeze)
        merged = sorted(before | set(frozen_from_layer))
        added = len(set(merged) - before)
        if echo_fn is not None:
            echo_fn(
                f"[layer] Applied freeze constraints from frozen layer: "
                f"total={len(merged)} (added_from_layer={added})"
            )
    else:
        merged = sorted(set(base_freeze))

    geom_cfg["freeze_atoms"] = merged
    calc_cfg["freeze_atoms"] = merged
    return merged


def deep_update(dst: Dict[str, Any], src: Optional[Dict[str, Any]]) -> Dict[str, Any]:
    """
    Recursively update mapping *dst* with *src*, returning *dst*.
    """
    for k, v in (src or {}).items():
        if isinstance(v, dict) and isinstance(dst.get(k), dict):
            deep_update(dst[k], v)
        else:
            dst[k] = v
    return dst


def collect_option_values(
    argv: _Sequence[str],
    names: _Sequence[str],
) -> List[str]:
    """Collect variadic option values in raw occurrence order.

    In addition to grouped/repeated options, accept Click's
    ``--long=value`` and attached-short ``-ivalue`` spellings.
    """
    vals: List[str] = []
    names_set = set(names)
    long_names = tuple(name for name in names_set if name.startswith("--"))
    short_names = tuple(
        name for name in names_set if name.startswith("-") and not name.startswith("--")
    )
    i = 0
    while i < len(argv):
        tok = argv[i]
        matched = tok in names_set
        inline: Optional[str] = None
        if not matched:
            for name in long_names:
                if tok.startswith(name + "="):
                    matched = True
                    inline = tok.split("=", 1)[1]
                    break
        if not matched:
            for name in short_names:
                if tok.startswith(name) and tok != name:
                    matched = True
                    inline = tok[len(name):]
                    break
        if not matched:
            i += 1
            continue
        if inline is not None:
            vals.append(inline)
        i += 1
        while i < len(argv) and not argv[i].startswith("-"):
            vals.append(argv[i])
            i += 1
    return vals


def current_cli_args(ctx: Optional[click.Context] = None) -> List[str]:
    """Return normalized arguments for the current top-level CLI invocation.

    ``CliRunner`` and other in-process callers do not rewrite ``sys.argv``.
    ``DefaultGroup`` therefore records the normalized token stream in Click's
    shared context metadata; direct command invocation falls back to the real
    process arguments for backward compatibility.
    """
    if ctx is None:
        ctx = click.get_current_context(silent=True)
    if ctx is not None:
        recorded = ctx.meta.get("mlmm.cli.raw_args")
        if recorded is not None:
            return [str(value) for value in recorded]
    return list(sys.argv[1:])


def reject_option_like_extra_args(
    extra_args: _Sequence[str],
    *,
    allowed_options: _Sequence[str] = (),
    allowed_values: _Sequence[str] = (),
    consumed_values: _Sequence[Any] = (),
) -> None:
    """Reject residual tokens not claimed by a legacy variadic option.

    Some commands accept ``-i A B`` or staged ``-s VALUE1 VALUE2`` syntax
    through Click's ``allow_extra_args`` compatibility mode.  Only values
    recovered from those declared variadic options may remain unparsed.
    """
    allowed = frozenset(str(value) for value in allowed_options)
    remaining = Counter(str(value) for value in allowed_values)
    for raw in consumed_values:
        value = str(raw)
        if remaining[value] > 0:
            remaining[value] -= 1
    for raw in extra_args:
        value = str(raw)
        if remaining[value] > 0:
            remaining[value] -= 1
            continue
        if value.startswith("-") and value not in allowed:
            raise click.UsageError(f"No such option: {value}")
        if value not in allowed:
            raise click.UsageError(f"Unexpected extra argument: {value}")


def collect_single_option_values(
    argv: _Sequence[str],
    names: _Sequence[str],
    label: str,
) -> List[str]:
    """Collect values following a flag that must appear at most once."""
    vals: List[str] = []
    seen = 0
    i = 0
    while i < len(argv):
        tok = argv[i]
        if tok in names:
            seen += 1
            j = i + 1
            while j < len(argv) and not argv[j].startswith("-"):
                vals.append(argv[j])
                j += 1
            i = j
        else:
            i += 1
    if seen > 1:
        raise click.BadParameter(
            f"Use a single {label} followed by multiple values; repeated flags are not accepted."
        )
    return vals


def load_pdb_atom_metadata(pdb_path: Path) -> List[Dict[str, Any]]:
    """Return atom metadata in file order, restoring original CIF identifiers."""
    atoms: List[Dict[str, Any]] = []
    with open(pdb_path, "r") as f:
        for line in f:
            if not (line.startswith("ATOM") or line.startswith("HETATM")):
                continue

            serial_txt = line[6:11].strip()
            resseq_txt = line[22:26].strip()
            atom_name = line[12:16].strip()
            altloc = line[16:17].strip()
            res_name = line[17:20].strip()
            chain_id = line[21:22].strip()
            icode = line[26:27].strip()
            element_txt = line[76:78].strip()
            is_hetatm = line.startswith("HETATM")

            try:
                serial = int(serial_txt) if serial_txt else None
            except ValueError:
                serial = None
            try:
                resseq = int(resseq_txt) if resseq_txt else None
            except ValueError:
                resseq = None

            if not element_txt:
                inferred = guess_element(atom_name, res_name, is_hetatm)
                element_txt = inferred or ""

            atoms.append(
                {
                    "serial": serial,
                    "name": atom_name,
                    "altloc": altloc,
                    "resname": res_name,
                    "resseq": resseq,
                    "chain": chain_id,
                    "icode": icode,
                    "element": element_txt,
                    "is_hetatm": is_hetatm,
                }
            )
    template = coordinate_template_for(pdb_path)
    if template is not None:
        if len(atoms) != template.natoms:
            raise ValueError(
                f"PDB metadata atom count ({len(atoms)}) does not match retained "
                f"template ({template.natoms}) for {pdb_path}."
            )
        for meta, record in zip(atoms, template.records):
            meta.update(
                {
                    "name": record.atom_name,
                    "altloc": record.altloc,
                    "resname": record.resname,
                    "resseq": (
                        int(record.resseq)
                        if re.fullmatch(r"[-+]?\d+", record.resseq)
                        else record.resseq
                    ),
                    "chain": record.chain_id,
                    "icode": record.icode,
                    "element": record.element,
                    "is_hetatm": record.group_pdb.upper() == "HETATM",
                }
            )
    return atoms


def _split_atom_spec_tokens(spec: str) -> List[str]:
    return [
        token
        for token in re.split(r"[\s/:`,\\]+", spec.strip().replace(" ", ","))
        if token
    ]


def resolve_atom_spec_index(spec: str, atom_meta: _Sequence[Dict[str, Any]]) -> int:
    """Resolve 3-field or ``CHAIN:RESNAME:RESSEQ[ICODE]:ATOM`` selectors."""
    tokens = _split_atom_spec_tokens(spec)
    if len(tokens) not in {3, 4}:
        raise ValueError(
            f"Atom spec '{spec}' must have 3 fields (resname, resseq, atomname) "
            "or 4 fields including chain ID."
        )

    tokens_upper = [t.upper() for t in tokens]
    canonical_four = spec.count(":") == 3 and all(
        part.strip() for part in spec.split(":")
    )
    canonical_parts = [part.strip() for part in spec.split(":")] if canonical_four else []
    matches: List[int] = []
    for idx, meta in enumerate(atom_meta):
        resname = (meta.get("resname") or "").strip().upper()
        resseq = meta.get("resseq")
        atom = (meta.get("name") or "").strip().upper()
        chain_text = (meta.get("chain") or "").strip()
        chain = chain_text.upper()
        if resseq is None:
            continue
        resseq_text = str(resseq)
        if canonical_four:
            chain_token, resname_token, resseq_token, atom_token = canonical_parts
            numbered = re.fullmatch(
                r"(?P<number>[-+]?\d+)(?P<icode>[A-Za-z]?)", resseq_token
            )
            if numbered is not None:
                try:
                    same_resseq = int(numbered.group("number")) == int(resseq_text)
                except ValueError:
                    same_resseq = numbered.group("number") == resseq_text
                requested_icode = numbered.group("icode").upper()
                if requested_icode:
                    same_resseq = same_resseq and requested_icode == str(
                        meta.get("icode") or ""
                    ).upper()
            else:
                same_resseq = resseq_token.upper() == resseq_text.upper()
            is_match = (
                chain_token == chain_text
                and resname_token.upper() == resname
                and same_resseq
                and atom_token.upper() == atom
            )
        else:
            normalized_tokens = [
                str(int(token)) if re.fullmatch(r"[-+]?\d+", token) else token
                for token in tokens_upper
            ]
            expected = [resname, resseq_text, atom]
            if len(tokens) == 4:
                expected.append(chain)
            is_match = Counter(normalized_tokens) == Counter(expected)
        if is_match:
            matches.append(idx)

    if len(matches) == 1:
        return matches[0]
    if len(matches) > 1:
        raise ValueError(
            f"Atom spec '{spec}' matches {len(matches)} atoms; add chain ID as "
            "CHAIN:RESNAME:RESSEQ[ICODE]:ATOM or use an explicit atom index."
        )
    if len(tokens) == 4 and not canonical_four:
        raise ValueError(
            f"Atom spec '{spec}' did not match any atom. Use the positional "
            "CHAIN:RESNAME:RESSEQ[ICODE]:ATOM form for chain-qualified selectors."
        )
    raise ValueError(f"Atom spec '{spec}' did not match any atom.")


def atom_label_from_meta(atom_meta: _Sequence[Dict[str, Any]], index: int) -> str:
    if index < 0 or index >= len(atom_meta):
        return f"idx{index}"
    meta = atom_meta[index]
    resname = (meta.get("resname") or "?").strip() or "?"
    resseq = meta.get("resseq")
    resseq_txt = "?" if resseq is None else str(resseq)
    icode = (meta.get("icode") or "").strip()
    if icode:
        resseq_txt += icode
    atom = (meta.get("name") or "?").strip() or "?"
    chain = (meta.get("chain") or "").strip()
    prefix = f"{chain}:" if chain else ""
    return f"{prefix}{resname}:{resseq_txt}:{atom}"


def axis_label_csv(
    axis_name: str,
    i_idx: int,
    j_idx: int,
    one_based: bool,
    atom_meta: Optional[_Sequence[Dict[str, Any]]] = None,
    pair_raw: Optional[Tuple[Any, Any, float, float]] = None,
) -> str:
    if pair_raw and (isinstance(pair_raw[0], str) or isinstance(pair_raw[1], str)) and atom_meta:
        i_label = atom_label_from_meta(atom_meta, i_idx)
        j_label = atom_label_from_meta(atom_meta, j_idx)
        return f"{axis_name}_{i_label}_{j_label}_A"
    i_disp = i_idx + 1 if one_based else i_idx
    j_disp = j_idx + 1 if one_based else j_idx
    return f"{axis_name}_{i_disp}_{j_disp}_A"


def axis_label_html(label: str) -> str:
    parts = label.split("_")
    if len(parts) >= 4 and parts[-1] == "A":
        axis = parts[0]
        i_disp = parts[1]
        j_disp = parts[2]
        return f"{axis} ({i_disp},{j_disp}) (Å)"
    return label


def resolve_scan_index(
    value: Any,
    *,
    one_based: bool,
    atom_meta: Optional[_Sequence[Dict[str, Any]]],
    context: str,
) -> int:
    """Resolve an index or atom-spec string for scan lists with consistent errors."""
    if isinstance(value, Integral):
        idx_val = int(value)
        if one_based:
            idx_val -= 1
        if idx_val < 0:
            raise click.BadParameter(
                f"Negative atom index after base conversion in {context}: {idx_val} (0-based expected)."
            )
        # Out-of-range upper bound: when atom metadata is available we know
        # the atom count, so reject an index past the end here (clean error,
        # consistent for bracketed and unbracketed specs) instead of letting
        # --dry-run pass and only crashing in the real run.
        if atom_meta and idx_val >= len(atom_meta):
            raise click.BadParameter(
                f"{context}: atom index {idx_val} (0-based) is out of range "
                f"for a {len(atom_meta)}-atom structure."
            )
        return idx_val
    if isinstance(value, str):
        if not atom_meta:
            raise click.BadParameter(
                f"{context} uses a string atom spec, but no PDB metadata is available."
            )
        try:
            return resolve_atom_spec_index(value, atom_meta)
        except ValueError as exc:
            raise click.BadParameter(f"{context} {exc}")
    raise click.BadParameter(f"{context} must be an int index or atom spec string.")


def parse_scan_list_triples(
    raw: str,
    *,
    one_based: bool,
    atom_meta: Optional[_Sequence[Dict[str, Any]]],
    option_name: str,
    return_one_based: bool = False,
) -> Tuple[List[Tuple[int, int, float]], List[Tuple[Any, Any, float]]]:
    """Parse --scan-lists entries into indices (0-based by default).

    Accepts both 3-tuples ``(i, j, target)`` and 4-tuples
    ``(i, j, start, end)`` for bidirectional scans.  4-tuples are
    expanded into two 3-tuple stages (initial→start, then initial→end)
    by the caller in scan.py.

    The returned *parsed* list contains tuples of length 3 **or** 4:
    ``(i, j, target)`` or ``(i, j, start, end)``.
    """
    try:
        obj = ast.literal_eval(raw)
    except Exception as e:
        raise click.BadParameter(f"Invalid literal for {option_name}: {e}")

    if not isinstance(obj, (list, tuple)):
        raise click.BadParameter(f"{option_name} must be a list/tuple of (i,j,target) or (i,j,start,end).")
    if len(obj) == 0:
        raise click.BadParameter(f"{option_name} must contain at least one atom pair.")

    parsed: list = []
    seen_pairs: set[tuple[int, int]] = set()
    for entry_idx, t in enumerate(obj, start=1):
        is_3 = (
            isinstance(t, (list, tuple))
            and len(t) == 3
            and isinstance(t[2], Real)
        )
        is_4 = (
            isinstance(t, (list, tuple))
            and len(t) == 4
            and isinstance(t[2], Real)
            and isinstance(t[3], Real)
        )
        if not (is_3 or is_4):
            raise click.BadParameter(
                f"{option_name} entry {entry_idx} must be (i,j,target) or (i,j,start,end): got {t}"
            )

        i = resolve_scan_index(
            t[0],
            one_based=one_based,
            atom_meta=atom_meta,
            context=f"{option_name} entry {entry_idx} (i)",
        )
        j = resolve_scan_index(
            t[1],
            one_based=one_based,
            atom_meta=atom_meta,
            context=f"{option_name} entry {entry_idx} (j)",
        )
        if i == j:
            raise click.BadParameter(
                f"{option_name} entry {entry_idx} selects the same atom twice."
            )
        pair_key = tuple(sorted((i, j)))
        if pair_key in seen_pairs:
            raise click.BadParameter(
                f"{option_name} entry {entry_idx} repeats atom pair "
                f"{pair_key}; each simultaneous scan axis must be unique."
            )
        seen_pairs.add(pair_key)
        if return_one_based:
            i += 1
            j += 1
        if is_4:
            parsed.append((i, j, float(t[2]), float(t[3])))
        else:
            parsed.append((i, j, float(t[2])))

    return parsed, list(obj)


def parse_dist_freeze_list(
    raw: str,
    *,
    one_based: bool,
    atom_meta: Optional[_Sequence[Dict[str, Any]]],
    option_name: str = "--dist-freeze",
) -> List[Tuple[int, int, Optional[float]]]:
    """Parse ``--dist-freeze`` entries: ``(i,j)`` or ``(i,j,target_A)``.

    Uses the same :func:`resolve_scan_index` as ``--scan-lists``, so string
    atom specs (e.g. ``'A:SER123:OG'``) are supported when PDB metadata is
    available.
    """
    try:
        obj = ast.literal_eval(raw)
    except Exception as e:
        raise click.BadParameter(f"Invalid literal for {option_name}: {e}")

    if not isinstance(obj, (list, tuple)):
        raise click.BadParameter(f"{option_name} must be a list/tuple of (i,j) or (i,j,target).")

    # Single tuple → wrap in list
    if obj and not isinstance(obj[0], (list, tuple)):
        obj = [obj]

    parsed: List[Tuple[int, int, Optional[float]]] = []
    seen_pairs: set[tuple[int, int]] = set()
    for entry_idx, t in enumerate(obj, start=1):
        if not (isinstance(t, (list, tuple)) and len(t) in (2, 3)):
            raise click.BadParameter(
                f"{option_name} entry {entry_idx} must be (i,j) or (i,j,target): got {t}"
            )
        i = resolve_scan_index(
            t[0], one_based=one_based, atom_meta=atom_meta,
            context=f"{option_name} entry {entry_idx} (i)",
        )
        j = resolve_scan_index(
            t[1], one_based=one_based, atom_meta=atom_meta,
            context=f"{option_name} entry {entry_idx} (j)",
        )
        if i == j:
            raise click.BadParameter(
                f"{option_name} entry {entry_idx} selects the same atom twice."
            )
        pair_key = tuple(sorted((i, j)))
        if pair_key in seen_pairs:
            raise click.BadParameter(
                f"{option_name} entry {entry_idx} repeats atom pair {pair_key}."
            )
        seen_pairs.add(pair_key)
        target: Optional[float] = None
        if len(t) == 3:
            if not isinstance(t[2], Real):
                raise click.BadParameter(
                    f"Target distance must be numeric in {option_name} entry {entry_idx}: {t}"
                )
            target = float(t[2])
            if target <= 0.0:
                raise click.BadParameter(
                    f"Target distance must be > 0 in {option_name} entry {entry_idx}: {t}"
                )
        parsed.append((i, j, target))
    return parsed


def parse_dist_freeze_spec(
    spec_path: Path,
    *,
    one_based_default: bool,
    atom_meta: Optional[_Sequence[Dict[str, Any]]],
    option_name: str = "--dist-freeze",
) -> List[Tuple[int, int, Optional[float]]]:
    """Parse a YAML/JSON dist-freeze spec file.

    Expected format::

        constraints:       # or "pairs" / "stages"
          - [1, 5, 1.4]   # (i, j, target_A) — target optional
          - [2, 6]         # freeze at current distance
        one_based: true    # optional, defaults to CLI value
    """
    spec_cfg = _load_scan_spec_root(spec_path, option_name=option_name)
    key, raw_list = _first_spec_field(spec_cfg, ("constraints", "pairs", "stages"))
    if key is None:
        raise click.BadParameter(
            f"{option_name} spec must define 'constraints', 'pairs', or 'stages'."
        )
    if not isinstance(raw_list, (list, tuple)) or len(raw_list) == 0:
        raise click.BadParameter(
            f"{option_name} field '{key}' must be a non-empty list."
        )

    one_based = _spec_one_based(
        spec_cfg.get("one_based"), default=one_based_default, option_name=option_name,
    )
    return parse_dist_freeze_list(
        repr(list(raw_list)),
        one_based=one_based,
        atom_meta=atom_meta,
        option_name=f"{option_name} {key}",
    )


def parse_scan_list_quads(
    raw: str,
    *,
    expected_len: int,
    one_based: bool,
    atom_meta: Optional[_Sequence[Dict[str, Any]]],
    option_name: str,
) -> Tuple[List[Tuple[int, int, float, float]], List[Tuple[Any, Any, float, float]]]:
    """Parse --scan-lists quadruples into 0-based indices."""
    try:
        obj = ast.literal_eval(raw)
    except Exception as e:
        raise click.BadParameter(f"Invalid literal for {option_name}: {e}")

    if not (isinstance(obj, (list, tuple)) and len(obj) == expected_len):
        quads = ",".join([f"(i{n},j{n},low{n},high{n})" for n in range(1, expected_len + 1)])
        raise click.BadParameter(
            f"{option_name} must contain exactly {expected_len} quadruples: [{quads}]"
        )

    parsed: List[Tuple[int, int, float, float]] = []
    seen_pairs: set[tuple[int, int]] = set()
    for entry_idx, q in enumerate(obj, start=1):
        if not (
            isinstance(q, (list, tuple))
            and len(q) == 4
            and isinstance(q[2], Real)
            and isinstance(q[3], Real)
        ):
            raise click.BadParameter(f"{option_name} entry must be (i,j,low,high): got {q}")

        i = resolve_scan_index(
            q[0],
            one_based=one_based,
            atom_meta=atom_meta,
            context=f"{option_name} entry {entry_idx} (i)",
        )
        j = resolve_scan_index(
            q[1],
            one_based=one_based,
            atom_meta=atom_meta,
            context=f"{option_name} entry {entry_idx} (j)",
        )
        if i == j:
            raise click.BadParameter(
                f"{option_name} entry {entry_idx} selects the same atom twice."
            )
        pair_key = tuple(sorted((i, j)))
        if pair_key in seen_pairs:
            raise click.BadParameter(
                f"{option_name} entry {entry_idx} repeats atom pair "
                f"{pair_key}; each scan axis must be unique."
            )
        seen_pairs.add(pair_key)
        parsed.append((i, j, float(q[2]), float(q[3])))

    for i, j, low, high in parsed:
        if low <= 0.0 or high <= 0.0:
            raise click.BadParameter(f"Distances must be positive: {(i, j, low, high)}")

    return parsed, list(obj)


def _load_scan_spec_root(
    spec_path: Path,
    *,
    option_name: str = "--scan-lists",
) -> Mapping[str, Any]:
    """Load a scan spec (YAML/JSON) and ensure mapping root."""
    try:
        with open(spec_path, "r", encoding="utf-8") as handle:
            data = yaml.safe_load(handle)
    except Exception as exc:
        raise click.BadParameter(
            f"Failed to parse {option_name} file '{spec_path}': {exc}"
        )

    if data is None:
        raise click.BadParameter(f"{option_name} file '{spec_path}' is empty.")
    if not isinstance(data, Mapping):
        raise click.BadParameter(
            f"{option_name} file '{spec_path}' must have a mapping at the YAML/JSON root."
        )
    # A misspelled root key used to be dropped in silence, so the whole spec file became a
    # no-op and the run continued with no scan/freeze constraints at all.
    unknown = sorted(set(data) - {"stages", "pairs", "constraints", "one_based"})
    if unknown:
        raise click.BadParameter(
            f"{option_name} file '{spec_path}' has unrecognized root key(s): "
            f"{', '.join(unknown)}. Expected any of: stages, pairs, constraints, one_based."
        )
    return data


def _spec_one_based(
    value: Any,
    *,
    default: bool,
    option_name: str = "--scan-lists",
) -> bool:
    """Resolve one_based value from spec with CLI fallback."""
    if value is None:
        return bool(default)
    if isinstance(value, bool):
        return value
    if isinstance(value, str):
        key = value.strip().lower()
        if key in {"1", "true", "yes", "y", "on"}:
            return True
        if key in {"0", "false", "no", "n", "off"}:
            return False
    raise click.BadParameter(
        f"{option_name} field 'one_based' must be a boolean (true/false)."
    )


def _first_spec_field(
    spec_cfg: Mapping[str, Any],
    candidates: _Sequence[str],
) -> Tuple[Optional[str], Any]:
    for key in candidates:
        if key in spec_cfg:
            return key, spec_cfg[key]
    return None, None


def is_scan_spec_file(value: str) -> bool:
    """Return True if *value* looks like an existing YAML/JSON scan spec file."""
    p = Path(value)
    return p.is_file() and p.suffix.lower() in {".yaml", ".yml", ".json"}


def parse_scan_spec_stages(
    spec_path: Path,
    *,
    one_based_default: bool,
    atom_meta: Optional[_Sequence[Dict[str, Any]]],
    option_name: str = "--scan-lists",
    return_bidirectional_markers: bool = False,
) -> Any:
    """Parse staged 1D scan spec into 0-based executable stages.

    A stage containing only ``(i, j, target)`` entries stays simultaneous.
    The legacy ``(i, j, start, end)`` form expands exactly like the inline
    CLI form: snapshot before the first leg, restore before the second.
    """
    spec_cfg = _load_scan_spec_root(spec_path, option_name=option_name)
    stages_key, stages_raw = _first_spec_field(spec_cfg, ("stages",))
    if stages_key is None:
        raise click.BadParameter(f"{option_name} must define 'stages'.")
    if not isinstance(stages_raw, (list, tuple)) or len(stages_raw) == 0:
        raise click.BadParameter(f"{option_name} field '{stages_key}' must be a non-empty list.")

    one_based = _spec_one_based(
        spec_cfg.get("one_based"), default=one_based_default, option_name=option_name
    )
    stages: List[List[Tuple[int, int, float]]] = []
    reset_before: set[int] = set()
    snapshot_before: set[int] = set()
    for stage_idx, stage_raw in enumerate(stages_raw, start=1):
        if not isinstance(stage_raw, (list, tuple)):
            raise click.BadParameter(
                f"{option_name} {stages_key}[{stage_idx}] must be a list of (i,j,target) entries."
            )
        parsed, _ = parse_scan_list_triples(
            repr(list(stage_raw)),
            one_based=one_based,
            atom_meta=atom_meta,
            option_name=f"{option_name} {stages_key}[{stage_idx}]",
        )
        if not parsed:
            raise click.BadParameter(
                f"{option_name} {stages_key}[{stage_idx}] must contain at least one (i,j,target) triple."
            )
        for entry in parsed:
            if any(float(distance) <= 0.0 for distance in entry[2:]):
                raise click.BadParameter(
                    f"Non-positive target distance in {option_name} "
                    f"{stages_key}[{stage_idx}]: {entry}."
                )
        if any(len(entry) == 4 for entry in parsed):
            for entry in parsed:
                if len(entry) == 4:
                    i, j, start, end = entry
                    first_leg = len(stages)
                    stages.append([(i, j, start)])
                    snapshot_before.add(first_leg)
                    reset_before.add(first_leg + 1)
                    stages.append([(i, j, end)])
                else:
                    stages.append([entry])
        else:
            stages.append(parsed)
    if return_bidirectional_markers:
        return (
            stages,
            one_based,
            frozenset(snapshot_before),
            frozenset(reset_before),
        )
    return stages, one_based


def parse_scan_spec_quads(
    spec_path: Path,
    *,
    expected_len: int,
    one_based_default: bool,
    atom_meta: Optional[_Sequence[Dict[str, Any]]],
    option_name: str = "--scan-lists",
) -> Tuple[List[Tuple[int, int, float, float]], List[Tuple[Any, Any, float, float]], bool]:
    """Parse 2D/3D scan spec into 0-based quad tuples."""
    spec_cfg = _load_scan_spec_root(spec_path, option_name=option_name)
    pairs_key, pairs_raw = _first_spec_field(spec_cfg, ("pairs",))
    if pairs_key is None:
        raise click.BadParameter(f"{option_name} must define 'pairs'.")
    if not isinstance(pairs_raw, (list, tuple)):
        raise click.BadParameter(f"{option_name} field '{pairs_key}' must be a list.")

    one_based = _spec_one_based(
        spec_cfg.get("one_based"), default=one_based_default, option_name=option_name
    )
    parsed, raw_pairs = parse_scan_list_quads(
        repr(list(pairs_raw)),
        expected_len=expected_len,
        one_based=one_based,
        atom_meta=atom_meta,
        option_name=f"{option_name} {pairs_key}",
    )
    return parsed, raw_pairs, one_based


PDB_ATOM_META_HEADER = f"{'id':>5} {'atom':<4} {'res':<4} {'resid':>4} {'el':<2}"


def format_pdb_atom_metadata(atom_meta: _Sequence[Dict[str, Any]], index: int) -> str:
    """Format metadata for atom *index* as aligned text: serial name resname resseq element."""
    fallback_serial = index + 1
    if index < 0 or index >= len(atom_meta):
        return f"{fallback_serial:>5} {'?':<4} {'?':<4} {'?':>4} {'?':<2}"

    meta = atom_meta[index]
    serial = meta.get("serial") or fallback_serial
    name = meta.get("name") or "?"
    resname = meta.get("resname") or "?"
    resseq = meta.get("resseq")
    resseq_str = "?" if resseq is None else str(resseq)
    element = (meta.get("element") or "?").strip() or "?"

    return f"{serial:>5} {name:<4} {resname:<4} {resseq_str:>4} {element:<2}"


def normalize_choice(
    value: str,
    *,
    param: str,
    alias_groups: Sequence[Tuple[Sequence[str], str]],
    allowed_hint: str,
) -> str:
    """Normalize a mode choice using alias groups and raise error on failure.

    Parameters
    ----------
    value : str
        The value to normalize.
    param : str
        Parameter name for error messages.
    alias_groups : Sequence[Tuple[Sequence[str], str]]
        Sequence of (aliases, canonical) pairs where aliases is a sequence of strings.
    allowed_hint : str
        Description of allowed values for error messages.

    Returns
    -------
    str
        The canonical value corresponding to the matched alias.

    Raises
    ------
    click.BadParameter
        If the value does not match any alias.
    """
    key = (value or "").strip().lower()
    for aliases, canonical in alias_groups:
        if any(key == alias.lower() for alias in aliases):
            return canonical

    hint = allowed_hint.strip()
    detail = f" Allowed: {hint}." if hint else ""
    raise click.BadParameter(f"Unknown value for {param} '{value}'.{detail}")


def _get_mapping_section(cfg: Mapping[str, Any], path: _Sequence[str]) -> Optional[Dict[str, Any]]:
    cur: Any = cfg
    for key in path:
        if not isinstance(cur, Mapping):
            return None
        cur = cur.get(key)
        if cur is None:
            return None
    return cur if isinstance(cur, dict) else None


def apply_yaml_overrides(
    yaml_cfg: Mapping[str, Any],
    overrides: _Sequence[Tuple[Dict[str, Any], _Sequence[_Sequence[str]]]],
) -> None:
    """Apply YAML overrides to multiple target dictionaries.

    Parameters
    ----------
    yaml_cfg : Mapping[str, Any]
        Parsed YAML configuration (root-level mapping).
    overrides : Sequence[Tuple[Dict[str, Any], Sequence[Sequence[str]]]]
        Each entry consists of the target dictionary to update followed by one or
        more candidate key paths. The first existing path is used. For example::

            apply_yaml_overrides(
                yaml_cfg,
                [
                    (geom_cfg, (("geom",),)),
                    (lbfgs_cfg, (("stopt", "lbfgs"), ("lbfgs",))),
                ],
            )

        Candidate paths are checked in order and the first mapping is applied.
    """
    for target, paths in overrides:
        for path in paths:
            norm_path = tuple(path)
            section = _get_mapping_section(yaml_cfg, norm_path)
            if section is not None:
                deep_update(target, section)
                break
            # A present-but-unusable section is a silent no-op otherwise: the user wrote
            # `geom:` as a list/scalar/empty and the whole block is dropped.
            if len(norm_path) == 1 and norm_path[0] in yaml_cfg:
                click.echo(
                    f"[config] WARNING: YAML section '{norm_path[0]}' is not a mapping; ignored.",
                    err=True,
                )


def yaml_section_has_key(
    yaml_cfg: Mapping[str, Any],
    paths: _Sequence[_Sequence[str]],
    key: str,
) -> bool:
    """Return True when any candidate YAML section explicitly defines ``key``."""
    for path in paths:
        section = _get_mapping_section(yaml_cfg, tuple(path))
        if isinstance(section, Mapping) and (key in section):
            return True
    return False


def load_yaml_dict(path: Optional[Path]) -> Dict[str, Any]:
    """
    Load a YAML file whose root must be a mapping. Return an empty dict if *path* is None.
    """
    if not path:
        return {}

    try:
        with open(path, "r") as f:
            data = yaml.safe_load(f) or {}
    except yaml.YAMLError as e:
        # Surface a malformed --config/--override YAML as a clean Click
        # error (one line, exit 2) instead of a raw parser traceback. This
        # is parsed during Click option handling, above the subcommand's
        # own exception wrapper, so the conversion must happen here.
        raise click.BadParameter(f"invalid YAML in '{path}': {e}")

    if not isinstance(data, dict):
        # ValueError (not click.BadParameter): public API contract /
        # test_load_yaml_dict_rejects_non_mapping_root.
        raise ValueError(f"YAML root must be a mapping, got: {type(data)}")

    return data


# Plotly: Energy diagram builder
def build_energy_diagram(
    energies: Sequence[float],
    labels: Sequence[str],
    ylabel: str = "ΔE",
    baseline: bool = False,
    showgrid: bool = False,
) -> go.Figure:
    """
    Plot an energy diagram using Plotly.

    Parameters
    ----------
    energies : Sequence[float]
        Energies for each state (same unit). Values are plotted without conversion.
    labels : Sequence[str]
        Labels corresponding to each state (for example, ["R", "TS1", "IM1", "TS2", "P"]).
        Must be the same length as ``energies``.
    ylabel : str, optional
        Y-axis label (for example, "ΔE" or "ΔG"). Defaults to ``"ΔE"``.
    baseline : bool, optional
        If ``True``, draw a dotted baseline at the energy of the first state across the plot.
    showgrid : bool, optional
        If ``True``, show grid lines on both axes. Defaults to ``False``.

    Returns
    -------
    plotly.graph_objs.Figure
        Figure containing the energy diagram.

    Notes
    -----
    - Each state is rendered as a thick horizontal segment (width ``HLINE_WIDTH``).
    - Adjacent states are connected by dotted diagonal segments from the right end of
      the left state to the left end of the right state.
    - Segment length automatically shrinks with additional states so that gaps remain
      between neighbors.
    - X-axis ticks are centered on each state and labeled using ``labels``.
    """
    if len(energies) == 0:
        raise ValueError("`energies` must contain at least one value.")
    if len(energies) != len(labels):
        raise ValueError("`energies` and `labels` must have the same length.")

    n = len(energies)
    energies = [float(e) for e in energies]

    AXIS_WIDTH = 3
    FONT_SIZE = 18
    AXIS_TITLE_SIZE = 20
    HLINE_WIDTH = 6           # Width of the horizontal state segments
    CONNECTOR_WIDTH = 2       # Width of the dotted connectors
    LINE_COLOR = "#1C1C1C"
    GRID_COLOR = "lightgrey"

    # Geometry along the X axis (centers and segment lengths)
    # Place segment centers at 0.5, 1.5, 2.5, ... (equally spaced)
    centers = [i + 0.5 for i in range(n)]

    # Shorten the segment as n grows (min 0.35, max 0.85)
    # Examples: n=5 -> 0.7, n=10 -> 0.5, n>=20 -> 0.35
    seg_width = min(0.85, max(0.35, 0.90 - 0.04 * n))
    half = seg_width / 2.0

    lefts = [c - half for c in centers]
    rights = [c + half for c in centers]

    fig = go.Figure()

    # Baseline (dotted line at the first energy level)
    if baseline:
        fig.add_trace(
            go.Scatter(
                x=[lefts[0], rights[-1]],
                y=[energies[0], energies[0]],
                mode="lines",
                line=dict(color=GRID_COLOR, dash="dot", width=2),
                hoverinfo="skip",
                showlegend=False,
            )
        )

    # Horizontal segments for each state
    for i, (e, lab) in enumerate(zip(energies, labels)):
        fig.add_trace(
            go.Scatter(
                x=[lefts[i], rights[i]],
                y=[e, e],
                mode="lines",
                line=dict(color=LINE_COLOR, width=HLINE_WIDTH),
                hovertemplate=f"{lab}: %{{y:.6f}}<extra></extra>",
                showlegend=False,
            )
        )

    # Dotted diagonals between adjacent states (right end -> left end)
    for i in range(n - 1):
        fig.add_trace(
            go.Scatter(
                x=[rights[i], lefts[i + 1]],
                y=[energies[i], energies[i + 1]],
                mode="lines",
                line=dict(color=LINE_COLOR, width=CONNECTOR_WIDTH, dash="dot"),
                hoverinfo="skip",
                showlegend=False,
            )
        )

    # Add a small margin beyond the first/last segments on X
    xpad = max(0.08, 0.15 * (1.0 - seg_width))
    x_min = lefts[0] - xpad
    x_max = rights[-1] + xpad

    # Add vertical padding above and below
    y_min = min(energies)
    y_max = max(energies)
    span = max(1e-6, y_max - y_min)  # Avoid zero span even if all values match
    ypad_low = 0.10 * span
    ypad_high = 0.20 * span
    y_range = [y_min - ypad_low, y_max + ypad_high]

    xaxis_config = dict(
        range=[x_min, x_max],
        showline=True,
        linewidth=AXIS_WIDTH,
        linecolor=LINE_COLOR,
        mirror=True,
        ticks="inside",
        tickwidth=AXIS_WIDTH,
        tickcolor=LINE_COLOR,
        tickfont=dict(size=FONT_SIZE, color=LINE_COLOR),
        showgrid=showgrid,
        gridcolor=GRID_COLOR,
        gridwidth=0.5,
        zeroline=False,
        tickmode="array",
        tickvals=centers,
        ticktext=list(labels),
        title=dict(text="", font=dict(size=AXIS_TITLE_SIZE, color=LINE_COLOR)),
    )

    yaxis_config = dict(
        range=y_range,
        showline=True,
        linewidth=AXIS_WIDTH,
        linecolor=LINE_COLOR,
        mirror=True,
        ticks="inside",
        tickwidth=AXIS_WIDTH,
        tickcolor=LINE_COLOR,
        tickfont=dict(size=FONT_SIZE, color=LINE_COLOR),
        showgrid=showgrid,
        gridcolor=GRID_COLOR,
        gridwidth=0.5,
        zeroline=False,
        title=dict(text=ylabel, font=dict(size=AXIS_TITLE_SIZE, color=LINE_COLOR)),
    )

    fig.update_layout(
        xaxis=xaxis_config,
        yaxis=yaxis_config,
        plot_bgcolor="white",
        paper_bgcolor="white",
        margin=dict(l=80, r=40, t=40, b=80),
    )

    return fig


def convert_xyz_to_pdb(
    xyz_path: Path,
    ref_pdb_path: Path,
    out_pdb_path: Path,
) -> None:
    """Overlay coordinates from *xyz_path* onto the topology of *ref_pdb_path* and write to *out_pdb_path*.

    The reference PDB is used as a text template: only the coordinate columns
    (31–54) of ATOM/HETATM records are replaced with coordinates from the XYZ
    frames.  All other PDB metadata (atom names, residue info, element columns,
    chain IDs, B-factors, etc.) are preserved verbatim, avoiding element
    misidentification bugs in external PDB parsers (e.g., ASE reading ``ZN``
    atom names as nitrogen).

    Every frame is validated before the destination, a CIF companion, or the
    coordinate-template registry is changed.  Publication replaces the exact
    destination path atomically.
    """
    from ase.io import read as ase_read
    traj = ase_read(str(xyz_path), index=":", format="xyz")
    if not traj:
        raise ValueError(f"No frames found in {xyz_path}.")
    symbols = [frame.get_chemical_symbols() for frame in traj]
    positions_by_frame = [
        np.asarray(frame.get_positions(), dtype=float) for frame in traj
    ]
    template = coordinate_template_for(ref_pdb_path)
    if template is not None:
        validate_coordinate_template_symbols(symbols, template)
    pdb_content = render_pdb_coordinate_frames(
        ref_pdb_path,
        symbols,
        positions_by_frame,
    ).encode("utf-8")

    cif_content: Optional[bytes] = None
    if template is not None:
        cif_content = render_mmcif_frames(
            positions_by_frame,
            template,
        ).encode("utf-8")

    payloads = {Path(out_pdb_path): pdb_content}
    if cif_content is not None:
        payloads[Path(out_pdb_path).with_suffix(".cif")] = cif_content
    commit_payloads(Path(out_pdb_path), payloads)

    # Propagate original identifiers from a CIF/oversized-PDB bridge and emit
    # a public mmCIF companion directly from the unrounded XYZ coordinates.
    if template is not None:
        register_coordinate_template(out_pdb_path, template)
    else:
        unregister_coordinate_template(out_pdb_path)


_CONVERT_FILES_ENABLED: bool = True


def set_convert_file_enabled(enabled: bool) -> None:
    """Globally enable or disable XYZ/TRJ conversions to PDB outputs."""
    global _CONVERT_FILES_ENABLED
    _CONVERT_FILES_ENABLED = bool(enabled)


def is_convert_file_enabled() -> bool:
    """Check if convert-files is globally enabled."""
    return _CONVERT_FILES_ENABLED


def pdb_keys_from_line(line: str) -> Tuple[Tuple, Tuple]:
    """Extract robust keys from a PDB ATOM/HETATM record.

    Returns:
        key_full: (chain, resseq, icode, resname, atomname, altloc)
        key_simple: (chain, resseq, icode, atomname)
    """
    atomname = line[12:16].strip()
    altloc = line[16:17].strip()
    resname = line[17:20].strip()
    chain = line[21:22].strip()
    resseq_str = line[22:26].strip()
    try:
        resseq = int(resseq_str)
    except ValueError:
        resseq = -10**9  # unlikely sentinel when missing
    icode = line[26:27].strip()
    key_full = (chain, resseq, icode, resname, atomname, altloc)
    key_simple = (chain, resseq, icode, atomname)
    return key_full, key_simple


def collect_ml_atom_keys(model_pdb: Path) -> Tuple[set, set]:
    """Collect ML-region atom keys from model_pdb.

    Returns:
        keys_full: Set of (chain, resseq, icode, resname, atomname, altloc)
        keys_simple: Set of (chain, resseq, icode, atomname)
    """
    from typing import Set as SetType
    keys_full: SetType[Tuple] = set()
    keys_simple: SetType[Tuple] = set()
    try:
        with model_pdb.open("r") as fh:
            for line in fh:
                if line.startswith("ATOM") or line.startswith("HETATM"):
                    kf, ks = pdb_keys_from_line(line)
                    keys_full.add(kf)
                    keys_simple.add(ks)
    except OSError as exc:
        # Do not degrade to empty sets: the caller then writes B=10.00 for every atom, i.e.
        # publishes a PDB whose ML layer has silently vanished.
        raise OSError(
            f"Failed to read model PDB for layer annotation '{model_pdb}': {exc}"
        ) from exc
    return keys_full, keys_simple


def format_pdb_with_bfactor(line: str, b: float) -> str:
    """Return PDB line with B-factor field (cols 61-66) set to b (6.2f)."""
    if len(line) < 66:
        line = line.rstrip("\n")
        line = line + " " * max(0, 66 - len(line))
        line = line + "\n"
    bf_str = f"{b:6.2f}"
    # Preserve occupancy (cols 55-60), overwrite tempFactor (61-66).
    new_line = line[:60] + bf_str + line[66:]
    return new_line


def annotate_pdb_bfactors_inplace(
    pdb_path: Path,
    model_pdb: Path,
    freeze_indices_0based: Sequence[int],
    beta_ml: float = 0.0,
    beta_frz: float = 20.0,
    beta_both: float = 0.0,
    *,
    _emit_cif: bool = True,
) -> None:
    """Overwrite B-factors in-place using 3-layer encoding (ML=0, MovableMM=10, FrozenMM=20).

    - ML-region atoms: beta_ml (default 0.00)
    - frozen atoms: beta_frz (default 20.00)
    - ML ∩ frozen: beta_both (default 0.00, ML takes precedence)

    Indexing for 'frozen' is 0-based and resets at each MODEL.
    """
    ml_full, ml_simple = collect_ml_atom_keys(model_pdb)
    frozen_set = set(int(i) for i in (freeze_indices_0based or []))

    try:
        lines = pdb_path.read_text().splitlines(keepends=True)
    except OSError as exc:
        raise OSError(f"Failed to read PDB for B-factor annotation '{pdb_path}': {exc}") from exc

    out_lines: List[str] = []
    atom_idx = 0  # resets per MODEL

    for line in lines:
        rec = line[:6]
        if rec.startswith("MODEL"):
            # reset atom counter for each model
            atom_idx = 0
            out_lines.append(line)
            continue
        if rec.startswith("ATOM  ") or rec.startswith("HETATM"):
            kf, ks = pdb_keys_from_line(line)
            is_ml = (kf in ml_full) or (ks in ml_simple)
            is_frz = (atom_idx in frozen_set)
            if is_ml and is_frz:
                out_lines.append(format_pdb_with_bfactor(line, beta_both))
            elif is_ml:
                out_lines.append(format_pdb_with_bfactor(line, beta_ml))
            elif is_frz:
                out_lines.append(format_pdb_with_bfactor(line, beta_frz))
            else:
                out_lines.append(format_pdb_with_bfactor(line, 10.0))
            atom_idx += 1
        else:
            out_lines.append(line)

    try:
        pdb_path.write_text("".join(out_lines))
    except Exception as exc:
        raise OSError(f"Failed to annotate B factors in '{pdb_path}': {exc}") from exc

    # Keep a bridged CIF companion synchronized with the layer B factors. This
    # generic path uses the annotated PDB coordinates; callers that still hold
    # the original XYZ overwrite it below with the unrounded coordinates.
    template = coordinate_template_for(pdb_path)
    if template is not None and _emit_cif:
        write_pdb_as_mmcif(pdb_path, template, pdb_path.with_suffix(".cif"))


def convert_and_annotate_xyz_to_pdb(
    src_xyz_or_trj: Path,
    ref_pdb: Path,
    dst_pdb: Path,
    model_pdb: Path,
    freeze_indices_0based: Sequence[int],
) -> None:
    """Convert an XYZ/TRJ file to PDB and annotate B-factors with the 3-layer encoding.

    The complete trajectory is rendered and validated first. Annotation runs
    on a private PDB, then the annotated PDB and retained-metadata CIF are
    staged together and published with the PDB authoritative. Layer values
    match the `opt` workflow's PDB output:

      - ML-region atoms: 0.00
      - movable MM atoms: 10.00
      - frozen MM atoms: 20.00
      - ML ∩ frozen: 0.00 (ML takes precedence)
    """
    from ase.io import read as ase_read
    from mlmm.io.structure_formats import _pdb_frame_data

    trajectory = ase_read(str(src_xyz_or_trj), index=":", format="xyz")
    if not trajectory:
        raise ValueError(f"No frames found in {src_xyz_or_trj}.")
    symbols = [frame.get_chemical_symbols() for frame in trajectory]
    frames = [np.asarray(frame.get_positions(), dtype=float) for frame in trajectory]
    template = coordinate_template_for(ref_pdb)
    if template is not None:
        validate_coordinate_template_symbols(symbols, template)
    unannotated = render_pdb_coordinate_frames(ref_pdb, symbols, frames)

    # B-factor annotation is performed only on a private file.  No public PDB,
    # CIF companion, or registry entry changes until every frame, annotation,
    # and companion serialization has succeeded.
    with tempfile.TemporaryDirectory(prefix="mlmm_annotated_pdb_") as tmp_dir:
        private_pdb = Path(tmp_dir) / "annotated.pdb"
        private_pdb.write_text(unannotated, encoding="utf-8")
        annotate_pdb_bfactors_inplace(
            private_pdb,
            model_pdb=model_pdb,
            freeze_indices_0based=freeze_indices_0based,
            _emit_cif=False,
        )
        annotated_payload = private_pdb.read_bytes()
        _, occupancies, bfactors = _pdb_frame_data(private_pdb)

    payloads = {Path(dst_pdb): annotated_payload}
    if template is not None:
        cif_payload = render_mmcif_frames(
            frames,
            template,
            occupancy_frames=occupancies,
            bfactor_frames=bfactors,
        ).encode("utf-8")
        payloads[Path(dst_pdb).with_suffix(".cif")] = cif_payload

    commit_payloads(Path(dst_pdb), payloads)
    if template is not None:
        register_coordinate_template(dst_pdb, template)
    else:
        unregister_coordinate_template(dst_pdb)




@dataclass
class PreparedInputStructure:
    source_path: Path
    geom_path: Path
    original_path: Optional[Path] = None
    structure_template: Optional[CoordinateTemplate] = None
    _normalized_structures: List[Tuple[Path, Path]] = field(default_factory=list)

    @property
    def is_cif(self) -> bool:
        return bool(
            self.original_path is not None
            and self.original_path.suffix.lower() in CIF_SUFFIXES
        )

    @property
    def display_path(self) -> Path:
        return self.original_path or self.source_path

    def cleanup(self) -> None:
        for internal_path, tmp_dir in self._normalized_structures:
            cleanup_normalized_structure(internal_path, tmp_dir)
        self._normalized_structures.clear()

    def __enter__(self) -> "PreparedInputStructure":
        return self

    def __exit__(self, exc_type, exc, tb) -> None:
        self.cleanup()

    def __del__(self) -> None:
        try:
            self.cleanup()
        except Exception:
            pass


def prepare_input_structure(path: Path) -> PreparedInputStructure:
    """Normalize mmCIF and PDB-overflow inputs to an internal safe PDB."""
    path = Path(path)
    registered_template = coordinate_template_for(path)
    if registered_template is not None:
        return PreparedInputStructure(
            source_path=path,
            geom_path=path,
            original_path=registered_template.source_path,
            structure_template=registered_template,
        )
    if is_cif_path(path) or (
        path.suffix.lower() == ".pdb" and pdb_requires_normalization(path)
    ):
        internal, structure_template, tmp_dir = normalize_structure_to_pdb(path)
        return PreparedInputStructure(
            source_path=internal,
            geom_path=internal,
            original_path=path.resolve(),
            structure_template=structure_template,
            _normalized_structures=[(internal, tmp_dir)],
        )
    return PreparedInputStructure(
        source_path=path,
        geom_path=path,
        original_path=path,
    )


def _count_atoms_in_file(path: Path) -> int:
    """Count atoms in a structure file (PDB or XYZ)."""
    suffix = path.suffix.lower()
    if suffix == ".pdb":
        count = 0
        with open(path, "r") as f:
            for line in f:
                if line.startswith(("ATOM  ", "HETATM")):
                    count += 1
        return count
    elif suffix == ".xyz":
        # XYZ format: first line is atom count
        with open(path, "r") as f:
            first_line = f.readline().strip()
            try:
                return int(first_line)
            except ValueError:
                return 0
    return 0


def _ordered_elements_in_file(path: Path) -> Optional[List[str]]:
    """Read the first structure's ordered element symbols when supported."""
    try:
        from ase.io import read as ase_read

        atoms = ase_read(str(path), index=0)
        return [str(symbol).capitalize() for symbol in atoms.get_chemical_symbols()]
    except Exception:
        logger.debug("Failed to read ordered elements from %s", path, exc_info=True)
        return None


def apply_ref_pdb_override(
    prepared_input: PreparedInputStructure,
    ref_pdb: Optional[Path],
) -> Optional[Path]:
    """Use a reference PDB topology while keeping XYZ coordinates for geometry loading.

    When --ref-pdb is provided:
    - geom_path remains the original input (xyz) for high-precision coordinates
    - source_path is updated to ref_pdb for topology/residue information
    """
    import click
    if ref_pdb is None:
        return None
    ref_pdb = Path(ref_pdb).resolve()
    if ref_pdb.suffix.lower() not in ({".pdb"} | set(CIF_SUFFIXES)):
        raise click.BadParameter("--ref-pdb must be a .pdb, .cif, or .mmcif file.")
    prepared_ref = prepare_input_structure(ref_pdb)
    geom_count = _count_atoms_in_file(prepared_input.geom_path)
    ref_count = _count_atoms_in_file(prepared_ref.geom_path)
    if geom_count != ref_count:
        prepared_ref.cleanup()
        raise click.BadParameter(
            f"atom count mismatch: {prepared_input.geom_path.name} has {geom_count} atoms, "
            f"but --ref-pdb {ref_pdb.name} has {ref_count} atoms."
        )
    geom_elements = _ordered_elements_in_file(prepared_input.geom_path)
    ref_elements = _ordered_elements_in_file(prepared_ref.geom_path)
    if (
        geom_elements is not None
        and ref_elements is not None
        and geom_elements != ref_elements
    ):
        mismatch = next(
            (
                index
                for index, (geom_element, ref_element) in enumerate(
                    zip(geom_elements, ref_elements)
                )
                if geom_element != ref_element
            ),
            0,
        )
        prepared_ref.cleanup()
        raise click.BadParameter(
            "atom-order element mismatch at 1-based atom "
            f"{mismatch + 1}: geometry={geom_elements[mismatch]}, "
            f"--ref-pdb={ref_elements[mismatch]}. The reference topology must "
            "use the identical atom order."
        )
    prepared_input.source_path = prepared_ref.source_path
    prepared_input.structure_template = prepared_ref.structure_template
    if prepared_ref.structure_template is not None:
        prepared_input.original_path = ref_pdb
        prepared_input._normalized_structures.extend(prepared_ref._normalized_structures)
        prepared_ref._normalized_structures.clear()
    return prepared_input.source_path


def validate_endpoint_atom_identities(
    inputs: Sequence[PreparedInputStructure],
) -> None:
    """Require identical ordered topology identities for path endpoints."""

    import click

    fields = (
        "is_hetatm",
        "chain",
        "resname",
        "resseq",
        "icode",
        "name",
        "altloc",
        "element",
    )
    reference = load_pdb_atom_metadata(inputs[0].source_path)
    reference_identity = [
        tuple(atom.get(field) for field in fields) for atom in reference
    ]
    for input_index, prepared in enumerate(inputs[1:], start=2):
        atoms = load_pdb_atom_metadata(prepared.source_path)
        identity = [tuple(atom.get(field) for field in fields) for atom in atoms]
        if len(identity) != len(reference_identity):
            raise click.BadParameter(
                "Path endpoint atom count mismatch: "
                f"input #1 has {len(reference_identity)}, input #{input_index} "
                f"has {len(identity)}."
            )
        if identity != reference_identity:
            mismatch = next(
                index
                for index, (left, right) in enumerate(
                    zip(reference_identity, identity), start=1
                )
                if left != right
            )
            raise click.BadParameter(
                "Path endpoint ordered atom identity mismatch at atom "
                f"{mismatch} between input #1 and input #{input_index}."
            )


# Charge/spin preparation consumes the workflow-level charge-summary service.
# Workflow subcommands import it from ``mlmm.workflows.charge_prep`` to keep
# ``core`` independent of ``extract``.


def read_bfactors_from_pdb(pdb_path: Path) -> List[float]:
    """
    Read B-factor (temperature factor) values from a PDB file.

    Returns a list of B-factors in atom order (0-indexed).
    Only ATOM and HETATM records are processed.
    """
    bfactors: List[float] = []
    with open(pdb_path, "r") as f:
        for line in f:
            if line.startswith(("ATOM  ", "HETATM")):
                # B-factor is at columns 61-66 (1-indexed), i.e., [60:66]
                try:
                    bfac = float(line[60:66].strip())
                except (ValueError, IndexError):
                    bfac = 0.0
                bfactors.append(bfac)
    return bfactors


def parse_layer_indices_from_bfactors(
    bfactors: List[float],
    tolerance: float = 1.0,
) -> Dict[str, List[int]]:
    """
    Parse B-factor values into layer indices for 3-layer ML/MM system.

    B-factor encoding:
        0.0 (±tolerance): ML atoms
        10.0 (±tolerance): Movable MM atoms
        20.0 (±tolerance): Frozen MM atoms

    Parameters
    ----------
    bfactors : List[float]
        B-factor values for each atom (0-indexed).
    tolerance : float
        Tolerance for B-factor matching (default: 1.0).

    Returns
    -------
    Dict[str, List[int]]
        Dictionary with keys:
        - "ml_indices": ML region atoms
        - "hess_mm_indices": Compatibility key (empty in 3-layer encoding)
        - "movable_mm_indices": Movable MM atoms
        - "frozen_indices": Frozen atoms
        - "unassigned_indices": Atoms with B-factors not matching any layer
    """
    from mlmm.core.defaults import BFACTOR_ML, BFACTOR_HESS_MM, BFACTOR_MOVABLE_MM, BFACTOR_FROZEN

    ml_indices: List[int] = []
    hess_mm_indices: List[int] = []
    movable_mm_indices: List[int] = []
    frozen_indices: List[int] = []
    unassigned_indices: List[int] = []

    for i, bfac in enumerate(bfactors):
        if abs(bfac - BFACTOR_ML) <= tolerance:
            ml_indices.append(i)
        elif abs(bfac - BFACTOR_FROZEN) <= tolerance:
            frozen_indices.append(i)
        elif abs(bfac - BFACTOR_MOVABLE_MM) <= tolerance:
            movable_mm_indices.append(i)
        elif (
            BFACTOR_HESS_MM != BFACTOR_MOVABLE_MM
            and abs(bfac - BFACTOR_HESS_MM) <= tolerance
        ):
            hess_mm_indices.append(i)
        else:
            unassigned_indices.append(i)

    return {
        "ml_indices": ml_indices,
        "hess_mm_indices": hess_mm_indices,
        "movable_mm_indices": movable_mm_indices,
        "frozen_indices": frozen_indices,
        "unassigned_indices": unassigned_indices,
    }


def has_valid_layer_bfactors(bfactors: List[float], tolerance: float = 1.0) -> bool:
    """
    Check if PDB B-factors contain valid 3-layer encoding.

    Returns True if at least one atom has ML B-factor and the B-factors are
    predominantly in the expected range (0, 10, 20).
    """
    from mlmm.core.defaults import BFACTOR_ML, BFACTOR_HESS_MM, BFACTOR_MOVABLE_MM, BFACTOR_FROZEN

    valid_bfactors = {BFACTOR_ML, BFACTOR_MOVABLE_MM, BFACTOR_FROZEN, BFACTOR_HESS_MM}
    has_ml = False
    valid_count = 0

    for bfac in bfactors:
        for valid in valid_bfactors:
            if abs(bfac - valid) <= tolerance:
                valid_count += 1
                if abs(bfac - BFACTOR_ML) <= tolerance:
                    has_ml = True
                break

    # Consider valid if:
    # 1. Has at least one ML atom
    # 2. At least 80% of atoms have valid B-factors
    return has_ml and (valid_count / max(len(bfactors), 1) >= 0.8)


def parse_indices_string(indices_str: str, one_based: bool = True) -> List[int]:
    """
    Parse a comma-separated index string into a sorted list of 0-based ints.

    Supports ranges like "1-5" (inclusive). By default, inputs are 1-based.
    """
    import click
    if indices_str is None:
        return []
    tokens = [tok.strip() for tok in str(indices_str).replace(" ", ",").split(",") if tok.strip()]
    indices: List[int] = []
    for token in tokens:
        if "-" in token and not token.startswith("-"):
            parts = token.split("-")
            if len(parts) == 2 and parts[0] and parts[1]:
                try:
                    start = int(parts[0])
                    end = int(parts[1])
                except ValueError as exc:
                    raise click.BadParameter(f"Invalid range token in --model-indices: '{token}'") from exc
                if one_based:
                    start -= 1
                    end -= 1
                if start < 0 or end < 0 or start > end:
                    raise click.BadParameter(f"Invalid range in --model-indices: '{token}'")
                indices.extend(range(start, end + 1))
                continue
        try:
            value = int(token)
        except ValueError as exc:
            raise click.BadParameter(f"Invalid index in --model-indices: '{token}'") from exc
        if one_based:
            value -= 1
        if value < 0:
            raise click.BadParameter(f"--model-indices expects positive indices; got {value + (1 if one_based else 0)}")
        indices.append(value)
    return sorted(set(indices))


def write_model_pdb_from_indices(
    input_pdb_path: Path,
    output_pdb_path: Path,
    indices: Sequence[int],
) -> None:
    """
    Write a model PDB containing only atoms at the specified 0-based indices.
    """
    import click
    if not indices:
        raise ValueError("No indices provided to build model PDB.")
    n_atoms = _count_atoms_in_file(input_pdb_path)
    if n_atoms <= 0:
        raise ValueError(f"No atoms found in input PDB: {input_pdb_path}")
    for idx in indices:
        if idx < 0 or idx >= n_atoms:
            raise click.BadParameter(
                f"model index out of range: {idx} (valid: 0 <= idx < {n_atoms})"
            )

    keep = set(int(i) for i in indices)
    lines_out: List[str] = []
    atom_idx = 0
    with open(input_pdb_path, "r") as f:
        for line in f:
            if line.startswith(("ATOM  ", "HETATM")):
                if atom_idx in keep:
                    # Auto-fill element column (77-78) if missing
                    raw = line.rstrip("\n")
                    elem_field = raw[76:78].strip() if len(raw) >= 78 else ""
                    if not elem_field:
                        atom_name = raw[12:16].strip()
                        res_name = raw[17:20].strip()
                        is_hetatm = raw.startswith("HETATM")
                        elem = guess_element(atom_name, res_name, is_hetatm)
                        if elem:
                            padded = raw.ljust(76) + f"{elem:>2}" + "\n"
                            lines_out.append(padded)
                        else:
                            lines_out.append(line)
                    else:
                        lines_out.append(line)
                atom_idx += 1
    if not lines_out:
        raise ValueError("Model PDB would be empty; check indices and input PDB.")
    if not lines_out[-1].endswith("\n"):
        lines_out[-1] = lines_out[-1] + "\n"
    lines_out.append("END\n")
    with open(output_pdb_path, "w") as f:
        f.writelines(lines_out)


def build_model_pdb_from_indices(
    input_pdb_path: Path,
    out_dir: Path,
    indices: Sequence[int],
    *,
    label: str = "model_from_indices",
) -> Path:
    """
    Create a temporary model PDB under out_dir using explicit indices.
    """
    out_dir.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile(
        mode="w",
        suffix=".pdb",
        prefix=f"{label}_",
        dir=out_dir,
        delete=False,
    ) as tmp:
        tmp_path = Path(tmp.name)
    write_model_pdb_from_indices(input_pdb_path, tmp_path, indices)
    return tmp_path


def build_model_pdb_from_bfactors(
    input_pdb_path: Path,
    out_dir: Path,
    *,
    tolerance: Optional[float] = None,
    label: str = "model_from_bfactor",
) -> Tuple[Path, Dict[str, List[int]]]:
    """
    Create a model PDB using ML indices derived from B-factors.

    Returns (model_pdb_path, layer_info).
    """
    from mlmm.core.defaults import BFACTOR_TOLERANCE
    tol = BFACTOR_TOLERANCE if tolerance is None else float(tolerance)
    bfactors = read_bfactors_from_pdb(input_pdb_path)
    if not bfactors:
        raise ValueError(f"No ATOM/HETATM records found in {input_pdb_path}.")
    if not has_valid_layer_bfactors(bfactors, tolerance=tol):
        raise ValueError(
            "Invalid or missing layer B-factors (expected ~0/10/20). "
            "Provide --no-detect-layer with --model-pdb/--model-indices."
        )
    layer_info = parse_layer_indices_from_bfactors(bfactors, tolerance=tol)
    ml_indices = layer_info.get("ml_indices") or []
    if not ml_indices:
        raise ValueError("No ML atoms detected from B-factors (value ~0).")
    out_dir.mkdir(parents=True, exist_ok=True)
    tmp_path = out_dir / f"{label}.pdb"
    write_model_pdb_from_indices(input_pdb_path, tmp_path, ml_indices)
    return tmp_path, layer_info


def _resolve_model_indices_for_layer_filter(
    input_pdb_path: Path,
    model_pdb_path: Path,
) -> List[int]:
    """Map an explicit model PDB to zero-based input ordinals for layer filtering."""
    full_atoms = load_pdb_atom_metadata(input_pdb_path)
    model_atoms = load_pdb_atom_metadata(model_pdb_path)
    if not model_atoms:
        raise ValueError(f"No atoms found in model PDB: {model_pdb_path}")

    identity_fields = ("chain", "resname", "resseq", "icode", "name")
    by_identity: Dict[Tuple[Any, ...], List[int]] = {}
    for index, atom in enumerate(full_atoms):
        key = tuple(atom.get(field) for field in identity_fields)
        by_identity.setdefault(key, []).append(index)

    resolved: List[int] = []
    seen: set[Tuple[Any, ...]] = set()
    for atom in model_atoms:
        key = tuple(atom.get(field) for field in identity_fields)
        label = (
            f"{atom.get('chain') or '-'}:{atom.get('resname')}:"
            f"{atom.get('resseq')}{atom.get('icode') or ''}:"
            f"{atom.get('name')}"
        )
        if key in seen:
            raise ValueError(f"model_pdb contains duplicate atom identity {label}.")
        seen.add(key)
        matches = by_identity.get(key, [])
        if len(matches) != 1:
            if not matches:
                raise ValueError(
                    f"model_pdb atom {label} is absent from the input PDB."
                )
            raise ValueError(
                f"model_pdb atom {label} matches {len(matches)} input atoms; "
                "retain unique chain and insertion-code identifiers."
            )
        full_index = matches[0]
        full_element = str(full_atoms[full_index].get("element") or "").upper()
        model_element = str(atom.get("element") or "").upper()
        if full_element != model_element:
            raise ValueError(
                f"model_pdb atom {label} has element {model_element}, but the "
                f"matching input atom has element {full_element}."
            )
        resolved.append(full_index)

    if any(right <= left for left, right in zip(resolved, resolved[1:])):
        raise ValueError(
            "model_pdb atom order differs from the input PDB; preserve full-system "
            "file order when creating the ML-region PDB."
        )
    return resolved


def resolve_ml_layer_assignment(
    *,
    source_path: Path,
    out_dir_path: Path,
    model_pdb: Optional[Any],
    model_indices: Optional[Sequence[int]],
    detect_layer: bool,
    hess_cutoff: Optional[float],
    movable_cutoff: Optional[float],
    calc_cfg: Dict[str, Any],
    echo_fn=None,
) -> Tuple[Path, Optional[Dict[str, List[int]]]]:
    """Resolve the ML-region model PDB path + layer_info dict.

    Shared by ML/MM per-stage workflows. Explicit ``model_pdb`` or
    ``model_indices`` defines ML membership; when layer detection remains
    enabled, valid input B-factors still define the movable/frozen MM layers.
    With no explicit membership, B-factors define all layers. Mutates
    ``calc_cfg`` in place to set ``use_bfactor_layers``, ``hess_cutoff``,
    ``movable_cutoff``, and ``model_pdb``. Raises ``click.ClickException`` on
    user-input failure so each caller can map to the same exit-code semantics.

    Returns:
        (model_pdb_path, layer_info) — layer_info is present when valid
        B-factor movable/frozen layers were read, including with explicit ML
        membership.
    """
    import click as _click  # local import to keep core.utils click-light

    echo = echo_fn if echo_fn is not None else _click.echo
    detect_layer_eff = detect_layer
    if model_pdb is None and calc_cfg.get("model_pdb"):
        model_pdb = calc_cfg.get("model_pdb")

    # movable_cutoff implies full distance-based layer assignment.
    # hess_cutoff alone can be combined with --detect-layer.
    if movable_cutoff is not None:
        if detect_layer_eff:
            echo("[layer] --movable-cutoff provided; disabling --detect-layer.", err=True)
        detect_layer_eff = False

    layer_source_pdb = source_path
    if detect_layer_eff and layer_source_pdb.suffix.lower() != ".pdb":
        raise _click.ClickException("--detect-layer requires a PDB input (or --ref-pdb).")

    model_pdb_path: Optional[Path] = None
    layer_info: Optional[Dict[str, List[int]]] = None
    explicit_region = model_pdb is not None or bool(model_indices)

    if explicit_region:
        explicit_model_indices: Optional[List[int]] = None
        if model_pdb is not None:
            model_pdb_path = Path(model_pdb)
            try:
                explicit_model_indices = _resolve_model_indices_for_layer_filter(
                    layer_source_pdb,
                    model_pdb_path,
                )
            except Exception as e:
                raise _click.ClickException(str(e)) from e
        else:
            if layer_source_pdb.suffix.lower() != ".pdb":
                raise _click.ClickException(
                    "--model-indices requires a PDB input (or --ref-pdb)."
                )
            try:
                model_pdb_path = build_model_pdb_from_indices(
                    layer_source_pdb, out_dir_path, list(model_indices or [])
                )
            except Exception as e:
                raise _click.ClickException(str(e)) from e

        if detect_layer_eff:
            try:
                from mlmm.core.defaults import BFACTOR_TOLERANCE

                bfactors = read_bfactors_from_pdb(layer_source_pdb)
                if not has_valid_layer_bfactors(
                    bfactors, tolerance=BFACTOR_TOLERANCE
                ):
                    raise ValueError(
                        "Invalid or missing layer B-factors (expected ~0/10/20)."
                    )
                layer_info = parse_layer_indices_from_bfactors(
                    bfactors, tolerance=BFACTOR_TOLERANCE
                )
                if model_indices:
                    explicit_ml_indices = {int(index) for index in model_indices}
                else:
                    explicit_ml_indices = set(explicit_model_indices or [])

                # The explicit region owns ML membership. B-factors only assign
                # the remaining atoms to MM sublayers, so an explicit ML atom
                # can never also become a geometry-level frozen atom.
                mm_layer_keys = (
                    "hess_mm_indices",
                    "movable_mm_indices",
                    "frozen_indices",
                )
                normalized_layer_info = {
                    key: sorted(
                        set(int(index) for index in layer_info.get(key, []))
                        - explicit_ml_indices
                    )
                    for key in mm_layer_keys
                }
                assigned = explicit_ml_indices | {
                    index
                    for key in mm_layer_keys
                    for index in normalized_layer_info[key]
                }
                normalized_layer_info["ml_indices"] = sorted(explicit_ml_indices)
                normalized_layer_info["unassigned_indices"] = sorted(
                    set(range(len(bfactors))) - assigned
                )
                layer_info = normalized_layer_info
                calc_cfg["use_bfactor_layers"] = True
                echo(
                    "[layer] Using explicit ML membership with B-factor "
                    "movable/frozen MM layers."
                )
            except Exception as e:
                echo(
                    f"[layer] WARNING: {e} Explicit ML membership remains active.",
                    err=True,
                )
                calc_cfg["use_bfactor_layers"] = False
        else:
            calc_cfg["use_bfactor_layers"] = False
    elif detect_layer_eff:
        try:
            model_pdb_path, layer_info = build_model_pdb_from_bfactors(layer_source_pdb, out_dir_path)
            calc_cfg["use_bfactor_layers"] = True
            echo(
                f"[layer] Detected B-factor layers: ML={len(layer_info.get('ml_indices', []))}, "
                f"MovableMM={len(layer_info.get('movable_mm_indices', []))}, "
                f"FrozenMM={len(layer_info.get('frozen_indices', []))}"
            )
        except Exception as e:
            if model_pdb is None and not model_indices:
                raise _click.ClickException(str(e)) from e
            echo(f"[layer] WARNING: {e} Falling back to explicit ML region.", err=True)
            detect_layer_eff = False

    if not explicit_region and not detect_layer_eff:
        if model_pdb is None and not model_indices:
            raise _click.ClickException(
                "Provide --model-pdb or --model-indices when --no-detect-layer."
            )
        if model_pdb is not None:
            model_pdb_path = Path(model_pdb)
        else:
            if layer_source_pdb.suffix.lower() != ".pdb":
                raise _click.ClickException("--model-indices requires a PDB input (or --ref-pdb).")
            try:
                model_pdb_path = build_model_pdb_from_indices(
                    layer_source_pdb, out_dir_path, list(model_indices or [])
                )
            except Exception as e:
                raise _click.ClickException(str(e)) from e
        calc_cfg["use_bfactor_layers"] = False

    if model_pdb_path is None:
        raise _click.ClickException("Failed to resolve model PDB for the ML region.")

    calc_cfg["model_pdb"] = str(model_pdb_path)
    if hess_cutoff is not None:
        calc_cfg["hess_cutoff"] = hess_cutoff
    if movable_cutoff is not None:
        calc_cfg["movable_cutoff"] = movable_cutoff
        calc_cfg["use_bfactor_layers"] = False

    return model_pdb_path, layer_info


def write_layer_bfactors_to_pdb(
    input_pdb_path: Path,
    output_pdb_path: Path,
    ml_indices: List[int],
    hess_mm_indices: Optional[List[int]] = None,
    movable_mm_indices: Optional[List[int]] = None,
    frozen_indices: Optional[List[int]] = None,
) -> None:
    """
    Write a PDB file with B-factors set according to 3-layer assignments.

    B-factor encoding:
        ML atoms: 0.0
        Movable MM atoms: 10.0
        Frozen MM atoms: 20.0
        Hessian MM atoms: encoded with the same B-factor as movable MM

    Parameters
    ----------
    input_pdb_path : Path
        Source PDB file to read atom records from.
    output_pdb_path : Path
        Output PDB file path.
    ml_indices : List[int]
        0-based indices of ML region atoms.
    hess_mm_indices : Optional[List[int]]
        0-based indices of MM atoms with Hessian (written as movable B-factor).
    movable_mm_indices : Optional[List[int]]
        0-based indices of movable MM atoms without Hessian.
    frozen_indices : Optional[List[int]]
        0-based indices of frozen atoms.

    Notes
    -----
    Supports multi-MODEL PDB files (e.g., trajectories): atom index resets
    at each MODEL record.
    """
    from mlmm.core.defaults import BFACTOR_ML, BFACTOR_HESS_MM, BFACTOR_MOVABLE_MM, BFACTOR_FROZEN

    ml_set = set(ml_indices or [])
    hess_mm_set = set(hess_mm_indices or [])
    movable_mm_set = set(movable_mm_indices or [])
    frozen_set = set(frozen_indices or [])

    lines_out: List[str] = []
    atom_idx = 0

    with open(input_pdb_path, "r") as f:
        for line in f:
            rec = line[:6]
            # Reset atom counter at each MODEL record (for trajectory files)
            if rec.startswith("MODEL"):
                atom_idx = 0
                lines_out.append(line)
                continue

            if line.startswith(("ATOM  ", "HETATM")):
                # Determine B-factor for this atom
                if atom_idx in ml_set:
                    bfac = BFACTOR_ML
                elif atom_idx in hess_mm_set:
                    bfac = BFACTOR_HESS_MM
                elif atom_idx in movable_mm_set:
                    bfac = BFACTOR_MOVABLE_MM
                elif atom_idx in frozen_set:
                    bfac = BFACTOR_FROZEN
                else:
                    # Default: treat as movable MM (layer 3)
                    bfac = BFACTOR_MOVABLE_MM

                # Replace B-factor (columns 61-66, 1-indexed)
                # PDB format: columns 61-66 are B-factor with format %6.2f
                # Ensure line is long enough before modifying
                if len(line) >= 66:
                    new_line = line[:60] + f"{bfac:6.2f}" + line[66:]
                else:
                    # Pad line if too short
                    padded = line.rstrip("\n").ljust(66)
                    new_line = padded[:60] + f"{bfac:6.2f}" + "\n"
                lines_out.append(new_line)
                atom_idx += 1
            else:
                lines_out.append(line)

    with open(output_pdb_path, "w") as f:
        f.writelines(lines_out)


def update_pdb_bfactors_from_layers(
    pdb_path: Path,
    ml_indices: List[int],
    hess_mm_indices: Optional[List[int]] = None,
    movable_mm_indices: Optional[List[int]] = None,
    frozen_indices: Optional[List[int]] = None,
) -> None:
    """
    Update B-factors in a PDB file in-place based on layer assignments.

    This is a convenience wrapper that reads and writes to the same file.
    """
    import tempfile
    import shutil

    with tempfile.NamedTemporaryFile(mode="w", suffix=".pdb", delete=False) as tmp:
        tmp_path = Path(tmp.name)

    try:
        write_layer_bfactors_to_pdb(
            pdb_path,
            tmp_path,
            ml_indices,
            hess_mm_indices,
            movable_mm_indices,
            frozen_indices,
        )
        shutil.move(str(tmp_path), str(pdb_path))
    finally:
        if tmp_path.exists():
            tmp_path.unlink()



def _collect_environment_info() -> dict:
    """Collect compute environment info (CPU, RAM, GPU, VRAM, resolved device)."""
    import platform
    env: dict = {}
    try:
        import torch
        cuda_ok = torch.cuda.is_available()
        env["device"] = "cuda" if cuda_ok else "cpu"
        if cuda_ok:
            try:
                env["gpu_name"] = torch.cuda.get_device_name(0)
                props = torch.cuda.get_device_properties(0)
                vram = getattr(props, "total_memory", None) or getattr(props, "total_mem", None)
                if vram:
                    env["gpu_vram_gb"] = round(vram / 1e9, 1)
            except Exception:
                pass
            env["cuda_version"] = getattr(torch.version, "cuda", None) or "unknown"
    except Exception:
        env["device"] = "cpu"
    try:
        import os
        cpu_info = platform.processor()
        if not cpu_info or cpu_info == "x86_64":
            try:
                with open("/proc/cpuinfo") as f:
                    for line in f:
                        if "model name" in line:
                            cpu_info = line.split(":")[1].strip()
                            break
            except Exception:
                pass
        env["cpu"] = cpu_info or "unknown"
        env["n_cpus"] = os.cpu_count()
        try:
            import psutil
            env["ram_gb"] = round(psutil.virtual_memory().total / 1e9, 1)
        except ImportError:
            pass
    except Exception:
        pass
    return env


# Schema version for result/summary JSON. Version 2.0 removes the UMA-specific
# all-workflow energy keys in favor of backend-neutral MLIP keys.
RESULT_JSON_SCHEMA_VERSION = "2.0"

# Union of public command-specific values for ``status``. Each command exposes
# a narrower enum documented in docs/json-output.md.
RESULT_JSON_STATUS_VALUES = (
    "completed",
    "converged",
    "error",
    "failed",
    "not_converged",
    "ok",
    "partial",
    "stalled",
    "success",
    "unknown",
    "unverified",
)


def write_result_json(
    out_dir: Path,
    data: dict,
    *,
    command: str,
    elapsed_seconds: Optional[float] = None,
    filename: str = "result.json",
    also_write_summary_json: bool = True,
) -> Path:
    """Write a machine-readable result.json for a subcommand.

    The ``data`` dict is augmented with common envelope fields
    (``command``, ``mlmm_version``, ``schema_version``, ``status``,
    ``elapsed_seconds``, ``files``, ``environment``) and serialized as
    indented JSON.

    When ``also_write_summary_json`` is True (default) the same payload
    is mirrored to ``summary.json`` alongside the legacy ``result.json``
    so downstream consumers can use a single filename across every
    subcommand (the ``all`` and ``path-search`` runners already write
    ``summary.json``; the per-stage subcommands continue to write
    ``result.json`` for backward compatibility).

    Returns the primary path.  Any serialization, staging, or publication
    failure raises ``ResultCommitError`` (an ``OSError`` subclass).
    """
    try:
        from mlmm._version import __version__
    except ImportError:
        __version__ = "unknown"

    data = dict(data)
    data.setdefault("command", command)
    data.setdefault("mlmm_version", __version__)
    data.setdefault("schema_version", RESULT_JSON_SCHEMA_VERSION)
    if elapsed_seconds is not None:
        data["elapsed_seconds"] = round(elapsed_seconds, 3)
    data.setdefault("environment", _collect_environment_info())

    # Convert non-serializable objects for json.dump
    def _to_json(obj):
        if isinstance(obj, Path):
            return str(obj)
        if isinstance(obj, dict):
            return {k: _to_json(v) for k, v in obj.items()}
        if isinstance(obj, (list, tuple)):
            return [_to_json(i) for i in obj]
        if isinstance(obj, float) and not math.isfinite(obj):
            return None
        try:
            import numpy as _np
            if isinstance(obj, _np.generic):
                return _to_json(obj.item())
            if isinstance(obj, _np.ndarray):
                return _to_json(obj.tolist())
        except ImportError:
            pass
        try:
            import torch as _th
            if isinstance(obj, _th.Tensor):
                return _to_json(obj.detach().cpu().tolist())
        except ImportError:
            pass
        return obj

    from mlmm.core.result_commit import commit_json_exact, with_current_run_id

    payload = with_current_run_id(_to_json(data))
    dest = Path(out_dir) / filename
    mirrors = ()
    if also_write_summary_json and Path(filename).name != "summary.json":
        mirrors = (Path(out_dir) / "summary.json",)
    return commit_json_exact(dest, payload, mirrors=mirrors)


_ALLOW_CHARGE_MULT_MISMATCH = False


def set_allow_charge_mult_mismatch(value: bool = True) -> None:
    """Process-global toggle for ``--allow-charge-mult-mismatch`` (set by the CLI eager callback)."""
    global _ALLOW_CHARGE_MULT_MISMATCH
    _ALLOW_CHARGE_MULT_MISMATCH = bool(value)


def validate_charge_spin(elements, charge, multiplicity, source: Optional[str] = None):
    """Raise ValueError if sum_Z(elements) - charge has the wrong parity for multiplicity,
    unless ``--allow-charge-mult-mismatch`` was set (then log a warning and skip).

    ``source`` is an optional label (typically the ML-region/model PDB path) whose atoms
    were counted; it is appended to both the raise message and the skip warning so a
    full-system count (thousands of atoms) is distinguishable from an ML-only count.
    """
    from pysisyphus.elem_data import ATOMIC_NUMBERS

    multiplicity = int(multiplicity)
    if multiplicity < 1:
        raise ValueError(
            f"Spin multiplicity must be an integer >= 1, got {multiplicity}."
        )
    sum_z = sum(ATOMIC_NUMBERS[str(e).lower()] for e in elements)
    total = sum_z - int(charge)
    unpaired = multiplicity - 1
    counted_atoms = len(elements)
    source_suffix = f", source={source}" if source is not None else ""
    if total < unpaired or (total - unpaired) % 2:
        if _ALLOW_CHARGE_MULT_MISMATCH:
            import logging
            logging.getLogger(__name__).warning(
                "ML-region electron-parity check SKIPPED (--allow-charge-mult-mismatch): "
                "sum_Z=%d, charge=%d, total_electrons=%d, multiplicity=%d, counted_atoms=%d%s "
                "-- proceeding; make sure this charge/multiplicity is intentional.",
                sum_z, charge, total, multiplicity, counted_atoms, source_suffix,
            )
            return
        raise ValueError(
            f"ML region electron count inconsistent: sum_Z={sum_z}, charge={charge}, "
            f"total_electrons={total}, multiplicity={multiplicity}, counted_atoms={counted_atoms}"
            f"{source_suffix}. Adjust charge (-q) or "
            f"multiplicity (-m) so the electron count matches the spin state (e.g. -m 2 for an odd "
            f"electron count). Common cause: a covalently-modified residue whose ML/MM cut was not "
            f"capped -- include the bonded partner in the ML region. A very large counted_atoms "
            f"(e.g. the whole system) usually means an ML-region PDB was not restricted to the "
            f"B-factor=0 layer. If this charge/multiplicity is "
            f"intentional, pass --allow-charge-mult-mismatch to skip this check."
        )


# Compatibility re-export: the bounded-peak Hessian symmetrizer lives in the
# bundled-engine layer (``pysisyphus.normal_modes``) so the pure
# normal-mode kernel there stays free of any upward ``mlmm`` import. It is
# re-exported here so existing callers of
# ``mlmm.core.utils.symmetrize_inplace`` (backends/mlmm_calc, workflows/tsopt,
# tests) keep resolving to the SAME function object.
from pysisyphus.normal_modes import symmetrize_inplace  # noqa: F401,E402


# ---------------------------------------------------------------------------
# XYZTrajectoryWriter — per-frame streaming XYZ writer (tail-able)
# ---------------------------------------------------------------------------
# Uses the same open / write+flush / close convention as
# `pysisyphus.optimizers.Optimizer.out_trj_handle`, so live monitoring with
# `tail -f scan_trj.xyz` remains available.


class XYZTrajectoryWriter:
    """Per-frame streaming XYZ writer (open / write+flush / close).

    Usage::

        with XYZTrajectoryWriter(path, mode="w") as w:
            for geom in geoms:
                w.write_raw(_coords3d_to_xyz_string(geom))

    Parameters
    ----------
    path : str | os.PathLike
        Output trajectory path. Written in real time (tail-able).
    mode : {"w", "a"}
        ``"w"`` truncates on open; ``"a"`` appends. Matches ``open()`` semantics.

    Attributes
    ----------
    frames_written : int
        Cumulative count of ``write_*`` calls since ``__enter__``.
    """

    __slots__ = ("path", "_mode", "_handle", "frames_written")

    def __init__(self, path, *, mode: str = "w") -> None:
        if mode not in ("w", "a"):
            raise ValueError(
                f"XYZTrajectoryWriter mode must be 'w' or 'a', got {mode!r}"
            )
        self.path = Path(path)
        self._mode = mode
        self._handle = None
        self.frames_written = 0

    def __enter__(self):
        self.path.parent.mkdir(parents=True, exist_ok=True)
        self._handle = self.path.open(self._mode, encoding="utf-8")
        return self

    def __exit__(self, exc_type, exc, tb) -> None:
        if self._handle is not None:
            try:
                self._handle.flush()
            finally:
                self._handle.close()
            self._handle = None

    def write_block(self, frame_xyz: str, energy=None) -> None:
        """Write a pre-formatted XYZ frame; optionally rewrite line 2 with energy."""
        if energy is not None:
            lines = frame_xyz.splitlines()
            if len(lines) >= 2 and lines[0].strip().isdigit():
                lines[1] = f"{float(energy):.12f}"
            frame_xyz = "\n".join(lines)
        self.write_raw(frame_xyz)

    def write_raw(self, block: str, energy=None) -> None:
        """Write a pre-formatted XYZ block as-is; flush after every frame."""
        if self._handle is None:
            raise RuntimeError(
                "XYZTrajectoryWriter must be used as a context manager."
            )
        if energy is not None:
            self.write_block(block, energy=energy)
            return
        if not block.endswith("\n"):
            block = block + "\n"
        self._handle.write(block)
        self._handle.flush()
        self.frames_written += 1
