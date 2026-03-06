# mlmm_toolkit/summary_log.py

"""
User-friendly summary log writer used by ``path_search`` and ``all``.

The goal is to provide a compact, readable ``summary.log`` alongside the
``summary.yaml``. The log aggregates MEP details, segment barriers,
post-processing energies, 3-layer system info, and key output paths.
"""

from __future__ import annotations

import textwrap
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional, Sequence

from pysisyphus.constants import AU2KCALPERMOL
from . import __version__
from .defaults import MLMM_CALC_KW, BFACTOR_ML, BFACTOR_MOVABLE_MM, BFACTOR_FROZEN

REQUIRED_SUMMARY_PAYLOAD_KEYS: tuple[str, ...] = (
    "root_out_dir",
    "path_module_dir",
    "pipeline_mode",
    "segments",
    "energy_diagrams",
)


def normalize_summary_payload(payload: Dict[str, Any] | None) -> Dict[str, Any]:
    """Return a defensive payload with stable defaults for summary rendering."""
    raw = payload if isinstance(payload, dict) else {}
    out: Dict[str, Any] = dict(raw)
    out.setdefault("root_out_dir", "-")
    out.setdefault("path_module_dir", "-")
    out.setdefault("pipeline_mode", "-")
    out.setdefault("segments", [])
    out.setdefault("post_segments", [])
    out.setdefault("energy_diagrams", [])
    out.setdefault("mep", {})
    out.setdefault("key_files", {})
    return out


def _fmt_bool(val: Optional[Any]) -> str:
    if val is None:
        return "-"
    return "True" if bool(val) else "False"


def _shorten_path(path: Optional[Path], root_out: Optional[Path]) -> str:
    """Return a path string, preferring a relative form to ``root_out`` or its parent."""
    if not path:
        return "(not available)"

    path_obj = Path(path)

    if root_out:
        for base in (root_out, root_out.parent):
            try:
                return str(path_obj.relative_to(base))
            except ValueError:
                continue

    return str(path_obj)


def _format_energy_rows(
    labels: Sequence[str],
    energies_au: Optional[Sequence[Optional[float]]],
    energies_kcal: Optional[Sequence[Optional[float]]],
) -> List[str]:
    rows: List[str] = []
    try:
        energies_au_list = list(energies_au) if energies_au is not None else []
    except Exception:
        energies_au_list = []
    try:
        energies_kcal_list = list(energies_kcal) if energies_kcal is not None else []
    except Exception:
        energies_kcal_list = []
    base_e = energies_au_list[0] if energies_au_list else None

    for i, lab in enumerate(labels):
        abs_e = energies_au_list[i] if i < len(energies_au_list) else None
        rel_e = energies_kcal_list[i] if i < len(energies_kcal_list) else None
        if rel_e is None and abs_e is not None and base_e is not None:
            rel_e = (abs_e - base_e) * AU2KCALPERMOL

        abs_txt = f"{abs_e:14.6f}" if abs_e is not None else f"{'n/a':>14}"
        rel_txt = f"{rel_e:14.4f}" if rel_e is not None else f"{'n/a':>14}"
        rows.append(f"        {lab:<8}{abs_txt}    {rel_txt}")
    return rows


def _format_bond_changes(text: str, indent: int = 6) -> List[str]:
    if not text:
        return ["".rjust(indent) + "(no covalent changes detected)"]
    blocks = [ln.rstrip() for ln in textwrap.dedent(text).splitlines() if ln.strip()]
    return ["".rjust(indent) + ln for ln in blocks]


def _format_ts_imag_info(ts_info: Any) -> List[str]:
    if ts_info is None:
        return []

    lines: List[str] = ["    TS imaginary freq:"]
    n_imag: Optional[int] = None
    nu_imag: Optional[float] = None
    min_abs: Optional[float] = None

    if isinstance(ts_info, dict):
        n_imag = ts_info.get("n_imag")
        nu_imag = ts_info.get("nu_imag_max_cm") or ts_info.get("nu_imag_cm")
        min_abs = ts_info.get("min_abs_imag_cm")
        if nu_imag is None and ts_info.get("ts_imag_freq_cm"):
            nu_imag = ts_info.get("ts_imag_freq_cm")
    else:
        try:
            nu_imag = float(ts_info)
            n_imag = 1 if nu_imag is not None else None
        except Exception:
            nu_imag = None

    n_imag_txt = str(n_imag) if n_imag is not None else "-"
    lines.append(f"      n_imag       : {n_imag_txt}")

    nu_label = "\u03bd_imag (max)"
    if nu_imag is not None:
        lines.append(f"      {nu_label} : {nu_imag:.1f} cm^-1")
    else:
        lines.append(f"      {nu_label} : -")

    magnitude = min_abs if min_abs is not None else (abs(nu_imag) if nu_imag is not None else None)
    note: Optional[str] = None
    if n_imag is not None:
        if n_imag == 1:
            if magnitude is not None and magnitude < 100.0:
                note = "WARNING      : Imaginary frequency magnitude is small; TS may be poorly optimized."
            else:
                note = "NOTE         : OK (single imaginary mode)"
        elif n_imag == 0:
            note = "WARNING      : No imaginary frequency; structure may not be a TS."
        else:
            note = "WARNING      : Multiple imaginary frequencies; TS may be poorly optimized."
    elif nu_imag is not None:
        if magnitude is not None and magnitude < 100.0:
            note = "WARNING      : Imaginary frequency magnitude is small; TS may be poorly optimized."
        else:
            note = "NOTE         : Single imaginary frequency (count unavailable)"

    if note:
        lines.append(f"      {note}")

    return lines


def _format_layer_info(payload: Dict[str, Any]) -> List[str]:
    """Format 3-layer ML/MM system information with Hessian-target MM subset info."""
    lines: List[str] = []

    ml_atoms = payload.get("ml_atoms")
    hess_mm_atoms = payload.get("hess_mm_atoms")
    movable_mm_atoms = payload.get("movable_mm_atoms")
    frozen_atoms = payload.get("frozen_atoms")

    n_ml = len(ml_atoms) if ml_atoms else 0
    n_hess = len(hess_mm_atoms) if hess_mm_atoms else 0
    n_movable = len(movable_mm_atoms) if movable_mm_atoms else 0
    n_frozen = len(frozen_atoms) if frozen_atoms else 0

    n_movable_total = n_hess + n_movable
    lines.append("  3-Layer ML/MM System:")
    lines.append(f"    ML (B={BFACTOR_ML:.0f})             : {n_ml:6d} atoms")
    lines.append(f"    Movable MM (B={BFACTOR_MOVABLE_MM:.0f})    : {n_movable_total:6d} atoms")
    lines.append(f"      Hessian-target subset             : {n_hess:6d} atoms")
    lines.append(f"    Frozen MM (B={BFACTOR_FROZEN:.0f})     : {n_frozen:6d} atoms")
    lines.append(f"    Total                               : {n_ml + n_movable_total + n_frozen:6d} atoms")

    hess_cutoff = payload.get("hess_cutoff")
    movable_cutoff = payload.get("movable_cutoff")
    use_bfactor = payload.get("use_bfactor_layers")

    if hess_cutoff is not None or movable_cutoff is not None:
        lines.append(
            f"    hess_cutoff    : {hess_cutoff} A"
            if hess_cutoff is not None
            else "    hess_cutoff    : -"
        )
        lines.append(
            f"    movable_cutoff : {movable_cutoff} A"
            if movable_cutoff is not None
            else "    movable_cutoff : -"
        )
    elif use_bfactor:
        lines.append("    Layer source   : B-factor based")

    return lines


def _emit_energy_block(
    lines: List[str],
    title: str,
    payload: Optional[Dict[str, Any]],
    root_out: Optional[Path],
) -> None:
    if not payload:
        return
    labels: Sequence[str] = payload.get("labels") or ["R", "TS", "P"]
    energies_au = payload.get("energies_au")
    energies_kcal = payload.get("energies_kcal")
    lines.append(f"    -- {title} --")
    lines.append("       State   Abs [Eh]          Rel [kcal/mol]")
    lines.extend(_format_energy_rows(labels, energies_au, energies_kcal))

    diagram = payload.get("diagram") or payload.get("image")
    if diagram:
        lines.append(f"       Diagram  : {_shorten_path(diagram, root_out)}")
    structs: Dict[str, Any] = payload.get("structures", {})
    if structs:
        lines.append("       Structures:")
        for key in ("R", "TS", "P"):
            if key in structs:
                lines.append(f"         {key}: {_shorten_path(structs.get(key), root_out)}")


def _tree_rel_path(root: Path, p: Path) -> str:
    try:
        return p.relative_to(root).as_posix()
    except ValueError:
        return p.name


def _tree_annotate(annotations: Dict[str, str], rel: str) -> str:
    note = annotations.get(rel)
    return f"  # {note}" if note else ""


def _tree_leaf_files(dir_path: Path) -> Optional[List[str]]:
    try:
        inner_children = sorted(dir_path.iterdir(), key=lambda p: p.name.lower())
    except Exception:
        return None

    if any(p.is_dir() for p in inner_children):
        return None
    return [p.name for p in inner_children if p.is_file()]


def _walk_directory_tree(
    dir_path: Path,
    prefix: str,
    depth: int,
    *,
    root: Path,
    annotations: Dict[str, str],
    max_depth: int,
    max_entries: int,
    lines: List[str],
    entries_seen_ref: List[int],
) -> bool:
    try:
        children = sorted(
            dir_path.iterdir(),
            key=lambda p: (p.is_file(), p.name.lower()),
        )
    except Exception:
        return False

    for idx, child in enumerate(children):
        connector = "\u2514\u2500" if idx == len(children) - 1 else "\u251c\u2500"
        rel = _tree_rel_path(root, child)
        if child.is_dir():
            leaf_names = _tree_leaf_files(child) if depth < max_depth else None
            if leaf_names is not None:
                lines.append(f"{prefix}{connector} {child.name}/{_tree_annotate(annotations, rel)}")
                entries_seen_ref[0] += 1
                if entries_seen_ref[0] >= max_entries:
                    lines.append(
                        f"{prefix}   ... (truncated after {max_entries} entries)"
                    )
                    return True

                next_prefix = prefix + ("   " if idx == len(children) - 1 else "\u2502  ")
                grouped = ",".join(leaf_names)
                lines.append(f"{next_prefix}\u2514\u2500 {{{grouped}}}")
                entries_seen_ref[0] += 1
                if entries_seen_ref[0] >= max_entries:
                    lines.append(
                        f"{next_prefix}   ... (truncated after {max_entries} entries)"
                    )
                    return True
                continue

        name = child.name + ("/" if child.is_dir() else "")
        lines.append(f"{prefix}{connector} {name}{_tree_annotate(annotations, rel)}")
        entries_seen_ref[0] += 1
        if entries_seen_ref[0] >= max_entries:
            lines.append(f"{prefix}   ... (truncated after {max_entries} entries)")
            return True
        if child.is_dir() and depth < max_depth:
            next_prefix = prefix + ("   " if idx == len(children) - 1 else "\u2502  ")
            if _walk_directory_tree(
                child,
                next_prefix,
                depth + 1,
                root=root,
                annotations=annotations,
                max_depth=max_depth,
                max_entries=max_entries,
                lines=lines,
                entries_seen_ref=entries_seen_ref,
            ):
                return True
    return False


def _format_directory_tree(
    root: Path,
    annotations: Dict[str, str],
    max_depth: int = 4,
    max_entries: int = 200,
) -> List[str]:
    """Render a compact directory tree rooted at ``root``."""

    lines: List[str] = []
    entries_seen_ref = [0]
    lines.append(f"  {root.name}/" + _tree_annotate(annotations, "."))
    _walk_directory_tree(
        root,
        "  ",
        1,
        root=root,
        annotations=annotations,
        max_depth=max_depth,
        max_entries=max_entries,
        lines=lines,
        entries_seen_ref=entries_seen_ref,
    )
    return lines


def _segment_table_value(entry: Dict[str, Any], key: str, col_width: int) -> str:
    if entry.get("kind") == "bridge" and not key.startswith("mep_"):
        return "---".rjust(col_width)
    val = entry.get(key)
    if val is None:
        return "---".rjust(col_width)
    return f"{val:>{col_width}.2f}"


def _classify_diagram_method(diag: Dict[str, Any]) -> str:
    name = str(diag.get("name", "")).lower()
    ylabel_txt = str(diag.get("ylabel", "")).lower()

    if "g_dft" in name or "gibbs_dft" in name or ("gibbs" in ylabel_txt and "dft" in name):
        return "gibbs_dft_uma"
    if "dft" in name:
        return "dft"
    if "g_uma" in name or "gibbs" in name or "gibbs" in ylabel_txt:
        return "gibbs_uma"
    if "uma" in name:
        return "uma"
    return "mep"


def _format_diag_row(
    diag: Optional[Dict[str, Any]],
    label: str,
    col_width: int,
    states: Sequence[str],
    label_width: int,
) -> str:
    if not diag:
        values = " ".join("---".rjust(col_width) for _ in states)
        return f"    {label:<{label_width}} {values}"

    try:
        labels_iter = list(diag.get("labels", []))
    except Exception:
        labels_iter = []
    labels_map = {lab: i for i, lab in enumerate(labels_iter)}
    energies_raw = diag.get("energies_kcal", [])
    try:
        energies = list(energies_raw) if energies_raw is not None else []
    except Exception:
        energies = []
    row_vals: List[str] = []
    for st in states:
        idx = labels_map.get(st)
        val = energies[idx] if idx is not None and idx < len(energies) else None
        row_vals.append(f"{val:>{col_width}.2f}" if val is not None else "---".rjust(col_width))
    return f"    {label:<{label_width}} {' '.join(row_vals)}"


def write_summary_log(dest: Path, payload: Dict[str, Any]) -> None:
    """Write a user-friendly summary.log at ``dest`` from a pre-collected payload."""
    missing_required = [
        key for key in REQUIRED_SUMMARY_PAYLOAD_KEYS if not isinstance(payload, dict) or key not in payload
    ]
    payload = normalize_summary_payload(payload)

    root_out = payload.get("root_out_dir") or "-"
    root_out_path = Path(root_out) if root_out not in (None, "-") else None
    path_module = payload.get("path_module_dir") or "-"
    pipeline_mode = payload.get("pipeline_mode") or "-"
    charge = payload.get("charge")
    spin = payload.get("spin")
    command = payload.get("command") or payload.get("cli_command")

    lines: List[str] = []
    lines.append("========================================================================")
    lines.append("mlmm_toolkit summary.log")
    lines.append("========================================================================")
    if command:
        lines.append(f"Input              : {command}")
    if missing_required:
        lines.append(
            "Preflight note     : missing payload keys replaced with defaults -> "
            + ", ".join(missing_required)
        )
    lines.append(f"Root out_dir       : {root_out}")
    path_module_disp = (
        _shorten_path(path_module, root_out_path)
        if path_module not in (None, "-")
        else path_module
    )
    lines.append(f"Path module dir    : {path_module_disp}")
    lines.append(f"Pipeline mode      : {pipeline_mode}")
    lines.append(f"refine-path        : {_fmt_bool(payload.get('refine_path'))}")
    lines.append(f"TSOPT/IRC          : {_fmt_bool(payload.get('tsopt'))}")
    lines.append(f"Thermochemistry    : {_fmt_bool(payload.get('thermo'))}")
    lines.append(f"DFT single-point   : {_fmt_bool(payload.get('dft'))}")
    opt_mode_disp = payload.get("opt_mode") or "-"
    lines.append(
        f"Opt mode           : {opt_mode_disp}  (grad: lbfgs/dimer; hess: rfo/rsirfo)"
    )
    lines.append(f"MEP mode           : {payload.get('mep_mode') or '-'}")

    version_base = payload.get("code_version") or __version__
    version_txt = f"mlmm_toolkit {version_base}"
    lines.append(f"Code version       : {version_txt}")
    uma_model = payload.get("uma_model") or MLMM_CALC_KW.get("uma_model") or "-"
    lines.append(f"UMA model          : {uma_model}")
    lines.append(f"Total charge (ML)  : {charge if charge is not None else '-'}")
    lines.append(f"Multiplicity (2S+1): {spin if spin is not None else '-'}")

    freeze_atoms_raw = payload.get("freeze_atoms")
    if freeze_atoms_raw is None:
        freeze_atoms_iter: List[Any] = []
    else:
        try:
            freeze_atoms_iter = list(freeze_atoms_raw)
        except Exception:
            freeze_atoms_iter = []
    try:
        freeze_atoms_list = sorted({int(i) for i in freeze_atoms_iter})
    except Exception:
        freeze_atoms_list = []
    if freeze_atoms_list:
        lines.append(
            "Freeze atoms (0-based): " + ",".join(map(str, freeze_atoms_list))
        )
    lines.append("")

    # 3-layer ML/MM system info
    if any(payload.get(k) for k in ["ml_atoms", "hess_mm_atoms", "movable_mm_atoms", "frozen_atoms"]):
        lines.extend(_format_layer_info(payload))
        lines.append("")

    delta = "\u0394"
    dagger = "\u2021"

    mep = payload.get("mep", {}) or {}
    diag = mep.get("diagram") or {}
    lines.append("[1] Global MEP overview")
    lines.append(f"  Number of MEP images : {mep.get('n_images', '-')}")
    lines.append(f"  Number of segments   : {mep.get('n_segments', '-')}")
    if mep.get("traj_pdb"):
        lines.append(
            f"  MEP trajectory (PDB) : {_shorten_path(mep.get('traj_pdb'), root_out_path)}"
        )
    if mep.get("mep_plot"):
        lines.append(
            f"  MEP energy plot      : {_shorten_path(mep.get('mep_plot'), root_out_path)}"
        )
    lines.append("")
    lines.append(f"  MEP energy diagram ({delta}E, kcal/mol)")
    if diag:
        if diag.get("image"):
            lines.append(
                f"    Image : {_shorten_path(diag.get('image'), root_out_path)}"
            )
        lines.append(f"    State    {delta}E [kcal/mol]")
        labels = diag.get("labels", [])
        energies = diag.get("energies_kcal", [])
        for i, lab in enumerate(labels):
            rel = energies[i] if i < len(energies) else None
            rel_txt = f"{rel:9.4f}" if rel is not None else "   n/a"
            lines.append(f"        {lab:<8}{rel_txt}")
    else:
        lines.append("    (no diagram available)")

    segments: Iterable[Dict[str, Any]] = payload.get("segments", []) or []
    lines.append("")
    lines.append("[2] Segment-level MEP summary (ML/MM path)")
    if segments:
        for seg in segments:
            idx = int(seg.get("index", 0) or 0)
            tag = seg.get("tag", f"seg_{idx:03d}")
            kind = seg.get("kind", "seg")
            lines.append(f"  - Segment {idx:02d} [{kind}]  tag={tag}")
            barrier = seg.get("barrier_kcal")
            delta_e = seg.get("delta_kcal")
            b_txt = f"{barrier:7.2f}" if barrier is not None else "   n/a"
            d_txt = f"{delta_e:7.2f}" if delta_e is not None else "   n/a"
            lines.append(f"      {delta}E{dagger} = {b_txt} kcal/mol,  {delta}E = {d_txt} kcal/mol")
            lines.append("      Bond changes:")
            lines.extend(_format_bond_changes(str(seg.get("bond_changes", ""))))
    else:
        lines.append("  (no segment reports)")

    post_segments: Iterable[Dict[str, Any]] = payload.get("post_segments", []) or []
    segment_entries: Dict[int, Dict[str, Any]] = {}
    for seg in segments:
        idx = int(seg.get("index", 0) or 0)
        tag = seg.get("tag", f"seg_{idx:03d}")
        kind = seg.get("kind", "seg")
        entry = segment_entries.setdefault(
            idx, {"index": idx, "tag": tag, "kind": kind}
        )
        entry.setdefault("tag", tag)
        entry.setdefault("kind", kind)
        if seg.get("barrier_kcal") is not None:
            entry["mep_barrier"] = seg.get("barrier_kcal")
        if seg.get("delta_kcal") is not None:
            entry["mep_delta"] = seg.get("delta_kcal")
    lines.append("")
    lines.append("[3] Per-segment post-processing (TSOPT / Thermo / DFT)")
    if post_segments:
        for seg in post_segments:
            idx = int(seg.get("index", 0) or 0)
            tag = seg.get("tag", f"seg_{idx:02d}")
            kind = seg.get("kind", "seg")
            lines.append(f"  === Segment {idx:02d} ({kind}) tag={tag} ===")
            if seg.get("post_dir"):
                lines.append(
                    f"    Post-process dir : {_shorten_path(seg.get('post_dir'), root_out_path)}"
                )
            ts_imag = seg.get("ts_imag") or seg.get("ts_imag_freq_cm")
            lines.extend(_format_ts_imag_info(ts_imag))
            if seg.get("irc_plot"):
                lines.append(
                    f"    IRC plot         : {_shorten_path(seg.get('irc_plot'), root_out_path)}"
                )
            if seg.get("irc_traj"):
                lines.append(
                    f"    IRC trajectory   : {_shorten_path(seg.get('irc_traj'), root_out_path)}"
                )
            _emit_energy_block(
                lines, "ML/MM energies (TSOPT+IRC)", seg.get("uma"), root_out_path
            )
            _emit_energy_block(lines, "ML/MM Gibbs (thermo)", seg.get("gibbs_uma"), root_out_path)
            _emit_energy_block(lines, "DFT single-point", seg.get("dft"), root_out_path)
            _emit_energy_block(
                lines, "DFT//ML/MM Gibbs", seg.get("gibbs_dft_uma"), root_out_path
            )

            entry = segment_entries.setdefault(
                idx, {"index": idx, "tag": tag, "kind": kind}
            )
            entry.setdefault("tag", tag)
            entry.setdefault("kind", kind)
            if seg.get("mep_barrier_kcal") is not None:
                entry["mep_barrier"] = seg.get("mep_barrier_kcal")
            if seg.get("mep_delta_kcal") is not None:
                entry["mep_delta"] = seg.get("mep_delta_kcal")
            if seg.get("uma"):
                uma_payload = seg.get("uma") or {}
                if uma_payload.get("barrier_kcal") is not None:
                    entry["uma_barrier"] = uma_payload.get("barrier_kcal")
                if uma_payload.get("delta_kcal") is not None:
                    entry["uma_delta"] = uma_payload.get("delta_kcal")
            if seg.get("gibbs_uma"):
                g_payload = seg.get("gibbs_uma") or {}
                if g_payload.get("barrier_kcal") is not None:
                    entry["gibbs_uma_barrier"] = g_payload.get("barrier_kcal")
                if g_payload.get("delta_kcal") is not None:
                    entry["gibbs_uma_delta"] = g_payload.get("delta_kcal")
            if seg.get("dft"):
                dft_payload = seg.get("dft") or {}
                if dft_payload.get("barrier_kcal") is not None:
                    entry["dft_barrier"] = dft_payload.get("barrier_kcal")
                if dft_payload.get("delta_kcal") is not None:
                    entry["dft_delta"] = dft_payload.get("delta_kcal")
            if seg.get("gibbs_dft_uma"):
                gd_payload = seg.get("gibbs_dft_uma") or {}
                if gd_payload.get("barrier_kcal") is not None:
                    entry["gibbs_dft_uma_barrier"] = gd_payload.get("barrier_kcal")
                if gd_payload.get("delta_kcal") is not None:
                    entry["gibbs_dft_uma_delta"] = gd_payload.get("delta_kcal")
    else:
        lines.append("  (no post-processing results)")

    if segment_entries:
        table_rows = [
            (f"MEP {delta}E{dagger} [kcal/mol]", "mep_barrier"),
            (f"MEP {delta}E  [kcal/mol]", "mep_delta"),
            (f"UMA {delta}E{dagger} [kcal/mol]", "uma_barrier"),
            (f"UMA {delta}E  [kcal/mol]", "uma_delta"),
            (f"UMA {delta}G{dagger} [kcal/mol]", "gibbs_uma_barrier"),
            (f"UMA {delta}G  [kcal/mol]", "gibbs_uma_delta"),
            (f"DFT//UMA {delta}E{dagger} [kcal/mol]", "dft_barrier"),
            (f"DFT//UMA {delta}E  [kcal/mol]", "dft_delta"),
            (f"DFT//UMA {delta}G{dagger} [kcal/mol]", "gibbs_dft_uma_barrier"),
            (f"DFT//UMA {delta}G  [kcal/mol]", "gibbs_dft_uma_delta"),
        ]
        sorted_entries = [segment_entries[k] for k in sorted(segment_entries.keys())]
        headers = [f"{int(e.get('index', 0)):d}({e.get('tag', '-')})" for e in sorted_entries]
        label_width = max(len(label) for label, _ in table_rows) + 2
        col_width = max(max(len(h) for h in headers), 8)

        lines.append("")
        lines.append("  Segment overview table")
        lines.append(
            "    "
            + f"{'Seg':<{label_width}} "
            + " ".join(f"{h:>{col_width}}" for h in headers)
        )
        for label, key in table_rows:
            values = " ".join(_segment_table_value(entry, key, col_width) for entry in sorted_entries)
            lines.append(f"    {label:<{label_width}} {values}")

    lines.append("")
    lines.append("[4] Energy diagrams (overview)")
    diagrams: Iterable[Dict[str, Any]] = payload.get("energy_diagrams", []) or []
    diag_by_method: Dict[str, Dict[str, Any]] = {}
    state_order: List[str] = []

    if diagrams:
        for diag_payload in diagrams:
            image_path = diag_payload.get("image") or diag_payload.get("diagram")
            if image_path and ("post_seg" in str(image_path) or "tsopt_seg_" in str(image_path)):
                continue

            name = diag_payload.get("name", "diagram")
            ylabel = diag_payload.get("ylabel", f"{delta}E (kcal/mol)")
            lines.append(f"  {name}  (ylabel: {ylabel})")
            labels = diag_payload.get("labels", [])
            energies = diag_payload.get("energies_kcal", [])
            energy_label = f"{delta}G [kcal/mol]" if f"{delta}G" in str(ylabel) else f"{delta}E [kcal/mol]"
            lines.append(f"    State   {energy_label}")
            for i, lab in enumerate(labels):
                rel = energies[i] if i < len(energies) else None
                rel_txt = f"{rel:7.3f}" if rel is not None else "   n/a"
                lines.append(f"        {lab:<8}{rel_txt}")
            if diag_payload.get("image"):
                lines.append(
                    f"    Image : {_shorten_path(diag_payload.get('image'), root_out_path)}"
                )

            method_key = _classify_diagram_method(diag_payload)
            diag_by_method.setdefault(method_key, diag_payload)
            if not state_order and labels:
                state_order = list(labels)
    else:
        lines.append("  (no energy diagrams recorded)")

    if state_order and diag_by_method:
        lines.append("")
        lines.append("  Energy diagram overview table")

        table_rows = [
            (f"MEP {delta}E  [kcal/mol]", "mep"),
            (f"UMA {delta}E  [kcal/mol]", "uma"),
            (f"UMA {delta}G  [kcal/mol]", "gibbs_uma"),
            (f"DFT//UMA {delta}E  [kcal/mol]", "dft"),
            (f"DFT//UMA {delta}G  [kcal/mol]", "gibbs_dft_uma"),
        ]

        label_width = max(len(label) for label, _ in table_rows) + 2
        col_width = max(max(len(st) for st in state_order), 7)

        lines.append(
            "    "
            + f"{'State':<{label_width}} "
            + " ".join(f"{st:>{col_width}}" for st in state_order)
        )

        for label, method in table_rows:
            diag_payload = diag_by_method.get(method)
            lines.append(_format_diag_row(diag_payload, label, col_width, state_order, label_width))

    lines.append("")
    lines.append("[5] Output directory structure")

    key_files = payload.get("key_files") or {}
    annotations: Dict[str, str] = {Path(k).as_posix(): v for k, v in key_files.items()}

    default_notes = {
        "pockets": "Extracted pocket PDBs",
        "scan": "Staged scan outputs",
        "path_search": "Recursive GSM outputs",
        "path_opt": "Single-pass GSM outputs",
        "tsopt_single": "Single-structure TSOPT-only outputs",
        "mep_plot.png": "ML/MM MEP energy plot",
        "energy_diagram_MEP.png": "Compressed MEP diagram",
        "energy_diagram_UMA_all.png": "UMA R-TS-P energies (all segments)",
        "energy_diagram_G_UMA_all.png": "UMA Gibbs R-TS-P (all segments)",
        "energy_diagram_DFT_all.png": "DFT R-TS-P (all segments)",
        "energy_diagram_G_DFT_plus_UMA_all.png": "DFT//UMA Gibbs R-TS-P (all segments)",
        "irc_plot_all.png": "Aggregated IRC plot",
    }

    if root_out_path:
        path_dir = payload.get("path_dir")
        if path_dir:
            try:
                rel = Path(path_dir).relative_to(root_out_path).as_posix()
                annotations.setdefault(rel, "Primary path module outputs")
            except ValueError:
                pass

        for rel, desc in default_notes.items():
            if (root_out_path / rel).exists():
                annotations.setdefault(rel, desc)

        if root_out_path.exists():
            lines.extend(_format_directory_tree(root_out_path, annotations))
        else:
            lines.append("  (root output directory not found on disk)")
    else:
        lines.append("  (root output directory unknown)")

    dest.parent.mkdir(parents=True, exist_ok=True)
    dest.write_text("\n".join(lines) + "\n", encoding="utf-8")
