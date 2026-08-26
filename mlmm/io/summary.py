"""
User-friendly summary log writer used by ``path_search`` and ``all``.

The goal is to provide a compact, readable ``summary.log`` alongside the
``summary.json``. The log aggregates MEP details, segment barriers,
post-processing energies, 3-layer system info, and key output paths.
"""

from __future__ import annotations

import re
import textwrap
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional, Sequence

from pysisyphus.constants import AU2KCALPERMOL
from mlmm import __version__
from mlmm.core.defaults import (
    BFACTOR_ML,
    BFACTOR_MOVABLE_MM,
    BFACTOR_FROZEN,
    SEGMENTS_DIRNAME,
    WORK_DIRNAME,
    TS_IMAG_SOFT_WARN_CM,
)
from mlmm.core.output import mlip_model_label


def format_result_warning(
    reason: Any,
    *,
    refine_path: bool = False,
    flatten: bool = False,
) -> str:
    """Translate an internal scientific-status reason into user-facing English."""
    raw = str(reason or "").strip()
    lowered = raw.lower()
    segment_match = re.search(r"(?:^|:)segment_(\d+)(?::|$)", lowered)
    if segment_match is None:
        segment_match = re.fullmatch(r"missing:segment_(\d+)", lowered)
    segment = segment_match.group(1) if segment_match else None
    human_segment_match = re.match(r"segment\s+(\d+):\s*(.+)$", lowered)
    if segment is None and human_segment_match is not None:
        segment = human_segment_match.group(1)
    human_detail = (
        human_segment_match.group(2) if human_segment_match is not None else lowered
    )

    def scoped(message: str) -> str:
        if segment is None:
            return message[0].upper() + message[1:]
        return f"Segment {segment}: {message}"

    direction_matches = re.findall(
        r"(?:^|[;:])irc:(?:irc:)?(forward|backward):([^;]+)", lowered
    )
    if direction_matches:
        messages: List[str] = []
        direction_text = {
            "not_converged": "IRC did not converge. Review its trajectory and IRC log.",
            "convergence_unknown": (
                "IRC convergence could not be confirmed. Review its trajectory and IRC log."
            ),
            "energy_invalid": (
                "IRC did not produce a valid energy profile. Review its trajectory and IRC log."
            ),
        }
        for direction, direction_code in direction_matches:
            detail = direction_text.get(direction_code.strip())
            if detail is not None:
                messages.append(f"{direction.capitalize()} {detail}")
        if messages:
            return scoped(" ".join(dict.fromkeys(messages)))

    code = lowered.rsplit(":", 1)[-1].strip()
    if code.startswith("ts_optimization_"):
        status = code.removeprefix("ts_optimization_")
        if status == "not_converged":
            return scoped("TS optimization did not converge. Review the TS trajectory.")
        return scoped(
            f"TS optimization status is {status.replace('_', ' ')}. Review the TS trajectory."
        )
    if code.startswith("terminal_hessian_"):
        status = code.removeprefix("terminal_hessian_").replace("_", " ")
        return scoped(
            f"terminal Hessian status is {status}. Review the TSOPT result."
        )
    if code == "imaginary_mode_count_unavailable":
        return scoped("imaginary-mode validation was unavailable. Review the TSOPT result.")
    if code == "no_imaginary_reaction_mode":
        message = "no imaginary mode was detected."
        if not refine_path:
            message += " Consider --refine-path."
        return scoped(message)
    if code in {"tsopt_status_unknown", "status_unknown"}:
        return scoped("TSOPT status could not be confirmed. Review the TSOPT result.")
    aggregate_messages = {
        "no usable path segments or energy diagrams were produced": (
            "No usable path result was produced. Review the MEP output."
        ),
        "no usable energy diagram was produced": (
            "No usable energy diagram was produced. Review the MEP and TSOPT results."
        ),
        "requested post-processing produced no segment records": (
            "Post-processing produced no segment records. Review the TSOPT/IRC output."
        ),
        "requested post-processing record is missing": (
            "the post-processing record is missing. Review the segment output."
        ),
        "tsopt/irc refined mlip energies are missing": (
            "TSOPT/IRC energies are missing. Review the TSOPT and IRC outputs."
        ),
        "irc trajectory is missing": (
            "the IRC trajectory is missing. Review the IRC output."
        ),
        "mlip thermochemistry result is missing": (
            "the ML/MM thermochemistry result is missing. Review the frequency output."
        ),
        "ts imaginary-mode validation is missing": (
            "TS imaginary-mode validation is missing. Review the TSOPT result."
        ),
        "dft result is missing": (
            "the DFT result is missing. Review the DFT output."
        ),
        "dft//mlip/mm thermochemistry result is missing": (
            "the DFT//MLIP/MM thermochemistry result is missing. Review the DFT and frequency outputs."
        ),
        "dft failed for one or more ts-only states": (
            "DFT failed for one or more TS-only states. Review the DFT output."
        ),
    }
    if human_detail in aggregate_messages:
        return scoped(aggregate_messages[human_detail])
    imag_count_match = re.match(
        r"ts imaginary-mode validation found n_imag=(\d+), expected 1$",
        human_detail,
    )
    if imag_count_match is not None:
        n_imag = int(imag_count_match.group(1))
        recovery: List[str] = []
        if n_imag > 1 and not flatten:
            recovery.append("--flatten")
        if not refine_path:
            recovery.append("--refine-path")
        message = f"TS imaginary-mode validation found n_imag={n_imag}."
        if recovery:
            message += f" Consider {' '.join(recovery)}."
        return scoped(message)
    if code == "irc_endpoint_connectivity_unvalidated":
        return (
            "Bond-topology matching between the two IRC endpoints and the two "
            "input endpoint structures could not be validated. Review both IRC "
            "endpoint structures before using this result."
        )
    priority_messages = {
        "mep_not_converged": (
            "MEP optimization did not converge. Review the MEP trajectory and convergence log."
        ),
        "mep_convergence_unknown": (
            "MEP convergence could not be confirmed. Review the MEP result and log."
        ),
        "irc_missing": "IRC results are missing. Confirm that TS optimization and IRC completed.",
        "post_missing": "requested post-processing is incomplete. Review the post-processing log.",
        "not_converged": "the calculation did not converge. Review the trajectory and convergence log.",
        "convergence_unknown": "convergence could not be confirmed. Review the result and log.",
        "endpoint_hei": (
            "no reactive segment was identified; only an endpoint/HEI path is available. "
            "Review the path before using its barrier."
        ),
        "engine_nonconverged": "the path-search engine did not converge. Review the path-search log.",
        "irc_result_missing": "IRC result metadata are missing. Confirm that IRC completed and wrote result.json.",
        "irc_result_unreadable": "IRC result metadata could not be read. Review result.json and the IRC log.",
        "irc_status_unknown": "IRC completion status could not be confirmed. Review result.json and the IRC log.",
        "irc_partial": "IRC completed only partially. Review both directional trajectories and the IRC log.",
        "irc_failed": "IRC failed. Review the IRC log and generated trajectories.",
    }
    if "endpoint_hei" in lowered and "engine_nonconverged" in lowered:
        return (
            "No reactive segment was identified and the path-search engine did not converge. "
            "Review the path and path-search log."
        )
    endpoint_match = re.search(r":endpoint_opt:([a-z0-9_]+)$", lowered)
    if endpoint_match:
        label = endpoint_match.group(1).removesuffix("_converged").replace("_", " ")
        subject = (
            f"{label} optimization"
            if label.startswith("endpoint ")
            else f"{label} endpoint optimization"
        )
        return scoped(
            f"the {subject} did not converge or could not be confirmed. "
            "Review the endpoint structure and optimizer log."
        )
    if lowered.startswith("missing:segment_"):
        return scoped("the expected segment result is missing. Review the workflow outputs.")
    if code in priority_messages:
        return scoped(priority_messages[code])
    if not raw:
        return "Result validation did not complete. Review the run details."
    if "_" in raw and " " not in raw:
        raw = raw.replace("_", " ")
    message = raw[0].upper() + raw[1:]
    if message[-1] not in ".!?":
        message += "."
    return message

_CITATION_RECORDS: Dict[str, tuple[str, str]] = {
    "software": (
        "mlmm-toolkit",
        "Ohmura, T.; Inoue, S.; Terada, T. ML/MM toolkit — Towards Accelerated "
        "Mechanistic Investigation of Enzymatic Reactions. ChemRxiv (2025). "
        "https://doi.org/10.26434/chemrxiv-2025-jft1k",
    ),
    "pysisyphus": (
        "pysisyphus engine",
        "Steinmetzer, J.; Kupfer, S.; Gräfe, S. pysisyphus: Exploring potential "
        "energy surfaces in ground and excited states. Int. J. Quantum Chem. 121, "
        "e26390 (2021). https://doi.org/10.1002/qua.26390",
    ),
    "uma": (
        "UMA",
        "Wood, B. M.; Dzamba, M.; Fu, X.; Gao, M.; Shuaibi, M.; "
        "Barroso-Luque, L.; Abdelmaqsoud, K.; Gharakhanyan, V.; Kitchin, "
        "J. R.; Levine, D. S.; Michel, K.; Sriram, A.; Cohen, T. S.; Das, A.; "
        "Sahoo, S. J.; Rizvi, A.; Ulissi, Z. W.; Zitnick, C. L. UMA: A Family "
        "of Universal Models for Atoms. Advances in Neural Information "
        "Processing Systems 38, "
        "129391-129427 (2025). "
        "https://doi.org/10.52202/085713-4310",
    ),
    "orb_v3": (
        "Orb-v3",
        "Rhodes, B.; Vandenhaute, S.; Šimkus, V.; Gin, J.; Godwin, J.; "
        "Duignan, T.; Neumann, M. Orb-v3: Atomistic Simulation at Scale. "
        "arXiv 2025, arXiv:2504.06231. "
        "https://doi.org/10.48550/arXiv.2504.06231",
    ),
    "mace": (
        "MACE",
        "Batatia, I.; Kovács, D. P.; Simm, G. N. C.; Ortner, C.; Csányi, G. "
        "MACE: Higher Order Equivariant Message Passing Neural Networks for Fast "
        "and Accurate Force Fields. Advances in Neural Information Processing "
        "Systems 35, 11423-11436 (2022). "
        "https://doi.org/10.52202/068431-0830",
    ),
    "mace_design": (
        "MACE",
        "Batatia, I.; Batzner, S.; Kovács, D. P.; Musaelian, A.; Simm, G. N. "
        "C.; Drautz, R.; Ortner, C.; Kozinsky, B.; Csányi, G. The Design "
        "Space of E(3)-Equivariant Atom-Centered Interatomic Potentials. "
        "arXiv 2022, arXiv:2205.06643. "
        "https://doi.org/10.48550/arXiv.2205.06643",
    ),
    "omol25": (
        "OMol25",
        "Levine, D. S.; Shuaibi, M.; Spotte-Smith, E. W. C.; Taylor, M. G.; "
        "Hasyim, M. R.; Michel, K.; Batatia, I.; Csányi, G.; Dzamba, M.; "
        "Eastman, P.; Frey, N. C.; Fu, X.; Gharakhanyan, V.; Krishnapriyan, A. "
        "S.; Rackers, J. A.; Raja, S.; Rizvi, A.; Rosen, A. S.; Ulissi, Z.; "
        "Vargas, S.; Zitnick, C. L.; Blau, S. M.; Wood, B. M. The Open "
        "Molecules 2025 (OMol25) Dataset, Evaluations, and Models. arXiv 2025, "
        "arXiv:2505.08762. "
        "https://doi.org/10.48550/arXiv.2505.08762",
    ),
    "gsm_peters": (
        "Growing String Method (GSM)",
        "Peters, B.; Heyden, A.; Bell, A. T.; Chakraborty, A. A growing string method "
        "for determining transition states: Comparison to the nudged elastic band "
        "and string methods. J. Chem. Phys. 120, 7877-7886 (2004). "
        "https://doi.org/10.1063/1.1691018",
    ),
    "gsm_zimmerman": (
        "Growing String Method (GSM)",
        "Zimmerman, P. M. Reliable Transition State Searches Integrated with "
        "the Growing String Method. J. Chem. Theory Comput. 9, 3043-3050 "
        "(2013). https://doi.org/10.1021/ct400319w",
    ),
    "dmf": (
        "Direct Max Flux (DMF)",
        "Koda, S.-i.; Saito, S. Locating Transition States by Variational Reaction "
        "Path Optimization with an Energy-Derivative-Free Objective Function. "
        "J. Chem. Theory Comput. 20, 2798-2811 (2024). "
        "https://doi.org/10.1021/acs.jctc.3c01246",
    ),
    "fbenm": (
        "FB-ENM initialization for DMF",
        "Koda, S.-i.; Saito, S. Flat-Bottom Elastic Network Model for Generating "
        "Improved Plausible Reaction Paths. J. Chem. Theory Comput. 20, 7176-7187 "
        "(2024). https://doi.org/10.1021/acs.jctc.4c00792",
    ),
    "cfbenm": (
        "Correlated FB-ENM (CFB-ENM) initialization for DMF",
        "Koda, S.-i.; Saito, S. Correlated Flat-Bottom Elastic Network Model for "
        "Improved Bond Rearrangement in Reaction Paths. J. Chem. Theory Comput. "
        "21, 3513-3522 (2025). https://doi.org/10.1021/acs.jctc.4c01549",
    ),
    "rfo": (
        "RFO / P-RFO",
        "Banerjee, A.; Adams, N.; Simons, J.; Shepard, R. Search for stationary "
        "points on surfaces. J. Phys. Chem. 89, 52-57 (1985). "
        "https://doi.org/10.1021/j100247a015",
    ),
    "lbfgs": (
        "Limited-memory BFGS (L-BFGS)",
        "Liu, D. C.; Nocedal, J. On the limited memory BFGS method for large "
        "scale optimization. Math. Program. 45, 503-528 (1989). "
        "https://doi.org/10.1007/BF01589116",
    ),
    "rsprfo": (
        "RS-P-RFO",
        "Besalú, E.; Bofill, J. M. On the automatic restricted-step rational-"
        "function-optimization method. Theor. Chem. Acc. 100, 265-274 (1998). "
        "https://doi.org/10.1007/s002140050387",
    ),
    "rsirfo": (
        "RS-I-RFO",
        "Besalú, E.; Bofill, J. M. On the automatic restricted-step rational-"
        "function-optimization method. Theor. Chem. Acc. 100, 265-274 (1998). "
        "https://doi.org/10.1007/s002140050387",
    ),
    "trim": (
        "Trust-Region Image Minimization (TRIM)",
        "Helgaker, T. Transition-state optimizations by trust-region image "
        "minimization. Chem. Phys. Lett. 182, 503-510 (1991). "
        "https://doi.org/10.1016/0009-2614(91)90115-P",
    ),
    "dimer": (
        "Dimer transition-state search",
        "Henkelman, G.; Jónsson, H. A dimer method for finding saddle points on "
        "high dimensional potential surfaces using only first derivatives. "
        "J. Chem. Phys. 111, 7010-7022 (1999). "
        "https://doi.org/10.1063/1.480097",
    ),
    "eulerpc": (
        "Euler predictor-corrector IRC (EulerPC)",
        "Hratchian, H. P.; Frisch, M. J.; Schlegel, H. B. Steepest Descent "
        "Reaction Path Integration Using a First-Order Predictor-Corrector "
        "Method. J. Chem. Phys. 133, 224101 (2010). "
        "https://doi.org/10.1063/1.3514202",
    ),
    "eulerpc_hessian": (
        "Euler predictor-corrector IRC (EulerPC)",
        "Hratchian, H. P.; Schlegel, H. B. Accurate Reaction Paths Using a "
        "Hessian-Based Predictor-Corrector Integrator. J. Chem. Phys. 120, "
        "9918-9924 (2004). https://doi.org/10.1063/1.1724823",
    ),
    "qrrho": (
        "quasi-RRHO thermochemistry",
        "Grimme, S. Supramolecular Binding Thermodynamics by Dispersion-Corrected "
        "Density Functional Theory. Chem. Eur. J. 18, 9955-9964 (2012). "
        "https://doi.org/10.1002/chem.201200497",
    ),
}


def _method_citation_record_keys(payload: Dict[str, Any]) -> List[str]:
    """Resolve citations from methods actually selected in this run."""

    keys = ["software", "pysisyphus"]
    mlip_backend = str(payload.get("mlip_backend") or "").strip().lower()
    mlip_model = str(payload.get("mlip_model") or "").strip().lower()
    if mlip_backend == "uma":
        keys.extend(("uma", "omol25"))
    elif mlip_backend == "orb":
        if "orb_v3" in mlip_model or "orb-v3" in mlip_model:
            keys.append("orb_v3")
        if "omol" in mlip_model or "orbmol" in mlip_model:
            keys.append("omol25")
    elif mlip_backend == "mace":
        keys.extend(("mace", "mace_design"))
        if "omol" in mlip_model:
            keys.append("omol25")
    pipeline_mode = str(payload.get("pipeline_mode") or "").strip().lower()
    if pipeline_mode != "tsopt-only":
        mep_mode = str(payload.get("mep_mode") or "").strip().lower()
        if mep_mode == "gsm":
            keys.extend(("gsm_peters", "gsm_zimmerman"))
        elif mep_mode == "dmf":
            keys.extend(("dmf", "fbenm"))
            if bool(payload.get("dmf_correlated")):
                keys.append("cfbenm")

        path_opt_mode = str(
            payload.get("path_opt_mode") or payload.get("opt_mode") or ""
        ).strip().lower()
        if path_opt_mode in {"grad", "lbfgs"}:
            keys.append("lbfgs")
        elif path_opt_mode in {"hess", "rfo", "rsprfo", "rsirfo"}:
            keys.append("rfo")

    post_segments = payload.get("post_segments") or []
    tsopt_used = bool(payload.get("tsopt_executed")) or any(
        isinstance(segment, dict) and "endpoint_opt" in segment
        for segment in post_segments
    )
    thermo_used = bool(payload.get("thermo_executed")) or any(
        isinstance(segment, dict)
        and (
            "thermo_symmetry" in segment
            or "gibbs_mlip" in segment
            or "gibbs_dft_mlip" in segment
        )
        for segment in post_segments
    )
    if thermo_used:
        keys.append("qrrho")

    if tsopt_used:
        keys.extend(("eulerpc", "eulerpc_hessian"))
        legacy_post_mode = str(
            payload.get("post_opt_mode")
            or payload.get("opt_mode_post")
            or payload.get("opt_mode")
            or ""
        ).strip().lower()
        ts_opt_mode = str(
            payload.get("ts_opt_mode") or legacy_post_mode
        ).strip().lower()
        endpoint_opt_mode = str(
            payload.get("endpoint_opt_mode") or legacy_post_mode
        ).strip().lower()

        if ts_opt_mode in {"grad", "lbfgs", "dimer"}:
            keys.extend(("lbfgs", "dimer"))
        elif ts_opt_mode in {"hess", "rfo", "rsprfo"}:
            keys.extend(("rfo", "rsprfo"))
        elif ts_opt_mode == "rsirfo":
            keys.extend(("rfo", "rsirfo"))
        elif ts_opt_mode == "trim":
            keys.append("trim")

        if endpoint_opt_mode in {"grad", "lbfgs", "dimer"}:
            keys.append("lbfgs")
        elif endpoint_opt_mode in {
            "hess",
            "rfo",
            "rsprfo",
            "rsirfo",
            "trim",
        }:
            keys.append("rfo")

    return list(dict.fromkeys(keys))


def format_method_citations(
    payload: Dict[str, Any], *, header: str = "[6] Methods and citations"
) -> List[str]:
    """Return the citation block, headed for its destination.

    ``summary.log`` numbers its sections, so the default header is section 6.
    Standard output heads its blocks with ``====== ... ======`` instead; passing
    the default section header there would leak the log section number.
    """

    lines = [header, "Please cite the software and methods used:"]
    previous_method = None
    for index, reference in enumerate(method_references(payload), start=1):
        if reference["method"] != previous_method:
            lines.append(f"- {reference['method']}:")
            previous_method = reference["method"]
        lines.append(f"[{index}] {reference['citation']}")
    return lines


def method_references(payload: Dict[str, Any]) -> List[Dict[str, str]]:
    """Return machine-readable references for ``summary.json``."""

    references: List[Dict[str, str]] = []
    for key in _method_citation_record_keys(payload):
        label, citation = _CITATION_RECORDS[key]
        doi_url = citation.rsplit("https://doi.org/", 1)[-1]
        references.append(
            {
                "method": label,
                "citation": citation,
                "doi": doi_url,
            }
        )
    return references


def emit_method_citations(payload: Dict[str, Any]) -> None:
    """Echo the citation block on stdout, headed like every other section.

    Deliberately a bare ``print``: the citation block is a required release
    surface, and the console emitter is gated on a verbosity level that defaults
    to 0 (silent) outside the CLI entry point, which would drop it. Only the
    header differs from the ``summary.log`` copy.
    """

    print("\n".join(format_method_citations(
        payload, header="====== Citations & References ======"
    )))


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
    zero_cutoff: Optional[float] = None

    if isinstance(ts_info, dict):
        n_imag = ts_info.get("n_imag")
        nu_imag = ts_info.get("nu_imag_max_cm") or ts_info.get("nu_imag_cm")
        min_abs = ts_info.get("min_abs_imag_cm")
        zero_cutoff = ts_info.get("frequency_zero_cutoff_cm")
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
    if zero_cutoff is not None:
        lines.append(f"      zero cutoff  : {float(zero_cutoff):.2f} cm^-1")

    nu_label = "\u03bd_imag (max)"
    if nu_imag is not None:
        lines.append(f"      {nu_label} : {nu_imag:.1f} cm^-1")
    else:
        lines.append(f"      {nu_label} : -")

    magnitude = min_abs if min_abs is not None else (abs(nu_imag) if nu_imag is not None else None)
    note: Optional[str] = None
    if n_imag is not None:
        if n_imag == 1:
            if magnitude is not None and magnitude < TS_IMAG_SOFT_WARN_CM:
                note = "WARNING      : Imaginary frequency magnitude is small; TS may be poorly optimized."
            else:
                note = "NOTE         : OK (single imaginary mode)"
        elif n_imag == 0:
            note = "WARNING      : No imaginary frequency; structure may not be a TS."
        else:
            note = "WARNING      : Multiple imaginary frequencies; TS may be poorly optimized."
    elif nu_imag is not None:
        if magnitude is not None and magnitude < TS_IMAG_SOFT_WARN_CM:
            note = "WARNING      : Imaginary frequency magnitude is small; TS may be poorly optimized."
        else:
            note = "NOTE         : Single imaginary frequency (count unavailable)"

    if note:
        lines.append(f"      {note}")

    return lines


def _format_thermo_symmetry(provenance: Any) -> List[str]:
    if not isinstance(provenance, dict):
        return []
    entries: List[str] = []
    for label in ("R", "TS", "P", "E1", "E2"):
        state = provenance.get(label)
        if not isinstance(state, dict):
            continue
        number = state.get("symmetry_number")
        source = state.get("symmetry_number_source")
        if number is not None:
            entries.append(
                f"{label}={number}" + (f" ({source})" if source else "")
            )
    return ["    Thermo symmetry   : " + ", ".join(entries)] if entries else []


def _format_layer_info(payload: Dict[str, Any]) -> List[str]:
    """Format compact 3-layer ML/MM system information."""
    lines: List[str] = []

    counts = payload.get("layer_counts")
    ml_atoms = payload.get("ml_atoms")
    hess_mm_atoms = payload.get("hess_mm_atoms")
    movable_mm_atoms = payload.get("movable_mm_atoms")
    frozen_atoms = payload.get("frozen_atoms")

    if isinstance(counts, dict):
        n_ml = int(counts.get("ml", 0) or 0)
        n_hess = 0
        n_movable = int(counts.get("movable", 0) or 0)
        n_frozen = int(counts.get("frozen", 0) or 0)
    else:
        n_ml = len(ml_atoms) if ml_atoms else 0
        n_hess = len(hess_mm_atoms) if hess_mm_atoms else 0
        n_movable = len(movable_mm_atoms) if movable_mm_atoms else 0
        n_frozen = len(frozen_atoms) if frozen_atoms else 0

    n_movable_total = n_hess + n_movable
    lines.append("  3-Layer ML/MM System:")
    lines.append(f"    ML (B={BFACTOR_ML:.0f})             : {n_ml:6d} atoms")
    lines.append(f"    Movable MM (B={BFACTOR_MOVABLE_MM:.0f})    : {n_movable_total:6d} atoms")
    if not isinstance(counts, dict):
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
    for endpoint in (1, 2):
        value = payload.get(f"barrier_from_endpoint_{endpoint}_kcal")
        if value is not None:
            lines.append(
                f"       Barrier E{endpoint}->TS: {float(value):.4f} kcal/mol"
            )

    diagram = payload.get("diagram") or payload.get("image")
    if diagram:
        lines.append(f"       Diagram  : {_shorten_path(diagram, root_out)}")
    structs: Dict[str, Any] = payload.get("structures", {})
    if structs:
        lines.append("       Structures:")
        for key in labels:
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
    if dir_path.is_symlink():
        return None
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
        is_recursive_dir = child.is_dir() and not child.is_symlink()
        if is_recursive_dir:
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

        name = child.name + ("@" if child.is_symlink() else ("/" if is_recursive_dir else ""))
        lines.append(f"{prefix}{connector} {name}{_tree_annotate(annotations, rel)}")
        entries_seen_ref[0] += 1
        if entries_seen_ref[0] >= max_entries:
            lines.append(f"{prefix}   ... (truncated after {max_entries} entries)")
            return True
        if is_recursive_dir and depth < max_depth:
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


def _format_selected_directory_tree(
    root: Path,
    relative_paths: Sequence[str],
    annotations: Dict[str, str],
) -> List[str]:
    """Render only caller-selected current-run files and their ancestors."""

    tree: Dict[str, Any] = {}
    for raw in relative_paths:
        relative = Path(str(raw))
        if relative.is_absolute() or ".." in relative.parts or not relative.parts:
            continue
        node = tree
        for part in relative.parts[:-1]:
            child = node.setdefault(part, {})
            if not isinstance(child, dict):
                break
            node = child
        else:
            node.setdefault(relative.parts[-1], None)

    lines = [f"  {root.name}/" + _tree_annotate(annotations, ".")]

    def render(node: Dict[str, Any], prefix: str, parts: tuple[str, ...]) -> None:
        entries = sorted(
            node.items(), key=lambda item: (item[1] is None, item[0].lower())
        )
        for index, (name, child) in enumerate(entries):
            last = index == len(entries) - 1
            connector = "└─" if last else "├─"
            rel_parts = parts + (name,)
            rel = Path(*rel_parts).as_posix()
            is_dir = isinstance(child, dict)
            shown = name + ("/" if is_dir else "")
            lines.append(
                f"{prefix}{connector} {shown}{_tree_annotate(annotations, rel)}"
            )
            if is_dir:
                render(child, prefix + ("   " if last else "│  "), rel_parts)

    render(tree, "  ", ())
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
        return "gibbs_dft_mlip"
    if "dft" in name:
        return "dft"
    if "g_mlip" in name or "gibbs" in name or "gibbs" in ylabel_txt:
        return "gibbs_mlip"
    if "mlip" in name:
        return "mlip"
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
    ts_only = pipeline_mode == "tsopt-only"
    charge = payload.get("charge")
    spin = payload.get("spin")
    command = payload.get("command") or payload.get("cli_command")

    lines: List[str] = []
    lines.append("========================================================================")
    lines.append("mlmm summary.log")
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
    dft_enabled = payload.get("dft")
    dft_status_str = _fmt_bool(dft_enabled)
    if dft_enabled:
        dft_result = payload.get("dft_status")
        if dft_result == "failed":
            dft_status_str = "True (Failed)"
        elif dft_result == "converged":
            dft_status_str = "True (Converged)"
    lines.append(f"DFT single-point   : {dft_status_str}")
    dft_func_basis = payload.get("dft_func_basis")
    if dft_func_basis:
        lines.append(f"DFT functional/basis: {dft_func_basis}")
    opt_mode_disp = payload.get("opt_mode") or "-"
    lines.append(f"Opt mode           : {opt_mode_disp}")
    lines.append(f"MEP mode           : {payload.get('mep_mode') or '-'}")

    version_base = payload.get("code_version") or __version__
    version_txt = f"mlmm {version_base}"
    lines.append(f"Code version       : {version_txt}")
    mlip_backend = payload.get("mlip_backend") or "-"
    mlip_model = payload.get("mlip_model") or "-"
    mlip_precision = payload.get("mlip_precision") or "-"
    backend_label = {
        "uma": "UMA", "orb": "ORB", "mace": "MACE", "aimnet2": "AIMNet2",
    }.get(str(mlip_backend).lower(), mlip_backend)
    model_label = payload.get("mlip_model_label") or mlip_model_label(
        mlip_backend, mlip_model, payload.get("mlip_task")
    )
    lines.append(f"MLIP backend        : {backend_label}")
    lines.append(f"MLIP model          : {model_label}")
    lines.append(f"MLIP precision      : {mlip_precision}")
    execution_status = payload.get("execution_status")
    scientific_status = payload.get("scientific_status") or payload.get("status")
    if execution_status is not None:
        lines.append(f"Execution status    : {execution_status}")
    if scientific_status is not None:
        lines.append(f"Scientific status   : {scientific_status}")
    status_reasons = (
        payload.get("scientific_status_reasons")
        or payload.get("status_reasons")
        or []
    )
    if scientific_status not in (None, "success"):
        reasons = list(status_reasons) or [None]
        for reason in reasons:
            lines.append(
                "RESULT WARNING      : "
                + format_result_warning(
                    reason,
                    refine_path=bool(payload.get("refine_path")),
                    flatten=bool(payload.get("flatten")),
                )
            )
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
            "Freeze atoms (1-based): " + ",".join(str(i + 1) for i in freeze_atoms_list)
        )
    lines.append("")

    # 3-layer ML/MM system info
    if payload.get("layer_counts") or any(
        payload.get(k)
        for k in ["ml_atoms", "hess_mm_atoms", "movable_mm_atoms", "frozen_atoms"]
    ):
        lines.extend(_format_layer_info(payload))
        lines.append("")

    delta = "\u0394"
    dagger = "\u2021"

    mep = payload.get("mep", {}) or {}
    diag = mep.get("diagram") or {}
    lines.append(
        "[1] Refined TS/IRC overview" if ts_only else "[1] Global MEP overview"
    )
    if ts_only:
        lines.append(
            f"  Number of IRC frames : "
            f"{payload.get('n_images', mep.get('n_images', '-'))}"
        )
    else:
        lines.append(f"  Number of MEP images : {mep.get('n_images', '-')}")
    lines.append(
        f"  Number of segments   : "
        f"{payload.get('n_segments', mep.get('n_segments', '-')) if ts_only else mep.get('n_segments', '-')}"
    )
    if mep.get("traj_pdb"):
        lines.append(
            f"  MEP trajectory (PDB) : {_shorten_path(mep.get('traj_pdb'), root_out_path)}"
        )
    if mep.get("mep_plot"):
        lines.append(
            f"  MEP energy plot      : {_shorten_path(mep.get('mep_plot'), root_out_path)}"
        )
    lines.append("")
    lines.append(
        f"  Refined TS/endpoint energy diagram ({delta}E, kcal/mol)"
        if ts_only
        else f"  MEP energy diagram ({delta}E, kcal/mol)"
    )
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
    lines.append(
        "[2] Refined TS/endpoint summary (ML/MM)"
        if ts_only
        else "[2] Segment-level MEP summary (ML/MM path)"
    )
    if segments:
        for seg in segments:
            idx = int(seg.get("index", 0) or 0)
            tag = seg.get("tag", f"seg_{idx:03d}")
            kind = seg.get("kind", "seg")
            lines.append(f"  - Segment {idx:02d} [{kind}]  tag={tag}")
            barrier = seg.get("barrier_kcal")
            delta_e = seg.get("delta_kcal")
            if kind == "tsopt" and (
                seg.get("barrier_from_endpoint_1_kcal") is not None
                or seg.get("barrier_from_endpoint_2_kcal") is not None
            ):
                for endpoint in (1, 2):
                    value = seg.get(f"barrier_from_endpoint_{endpoint}_kcal")
                    value_text = (
                        f"{float(value):7.2f}" if value is not None else "   n/a"
                    )
                    lines.append(
                        f"      {delta}E{dagger}(E{endpoint}->TS) = "
                        f"{value_text} kcal/mol  [chemically unassigned endpoint]"
                    )
            else:
                b_txt = f"{barrier:7.2f}" if barrier is not None else "   n/a"
                d_txt = f"{delta_e:7.2f}" if delta_e is not None else "   n/a"
                lines.append(
                    f"      {delta}E{dagger} = {b_txt} kcal/mol,  "
                    f"{delta}E = {d_txt} kcal/mol  [MEP]"
                )
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
        prefix = "mlip" if kind == "tsopt" else "mep"
        if seg.get("barrier_kcal") is not None:
            entry[f"{prefix}_barrier"] = seg.get("barrier_kcal")
        if seg.get("delta_kcal") is not None:
            entry[f"{prefix}_delta"] = seg.get("delta_kcal")
        for endpoint in (1, 2):
            value = seg.get(f"barrier_from_endpoint_{endpoint}_kcal")
            if value is not None:
                entry[f"{prefix}_barrier_e{endpoint}"] = value
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
            lines.extend(_format_thermo_symmetry(seg.get("thermo_symmetry")))
            if seg.get("irc_plot"):
                lines.append(
                    f"    IRC plot         : {_shorten_path(seg.get('irc_plot'), root_out_path)}"
                )
            if seg.get("irc_traj"):
                lines.append(
                    f"    IRC trajectory   : {_shorten_path(seg.get('irc_traj'), root_out_path)}"
                )
            _emit_energy_block(
                lines, "ML/MM energies (TSOPT+IRC)", seg.get("mlip"), root_out_path
            )
            _emit_energy_block(lines, "ML/MM Gibbs (thermo)", seg.get("gibbs_mlip"), root_out_path)
            _emit_energy_block(
                lines,
                "model-region DFT single-point",
                seg.get("dft"),
                root_out_path,
            )
            _emit_energy_block(
                lines, "DFT//MLIP/MM Gibbs", seg.get("gibbs_dft_mlip"), root_out_path
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
            if seg.get("mlip"):
                mlip_payload = seg.get("mlip") or {}
                if mlip_payload.get("barrier_kcal") is not None:
                    entry["mlip_barrier"] = mlip_payload.get("barrier_kcal")
                if mlip_payload.get("delta_kcal") is not None:
                    entry["mlip_delta"] = mlip_payload.get("delta_kcal")
                for endpoint in (1, 2):
                    value = mlip_payload.get(
                        f"barrier_from_endpoint_{endpoint}_kcal"
                    )
                    if value is not None:
                        entry[f"mlip_barrier_e{endpoint}"] = value
            if seg.get("gibbs_mlip"):
                g_payload = seg.get("gibbs_mlip") or {}
                if g_payload.get("barrier_kcal") is not None:
                    entry["gibbs_mlip_barrier"] = g_payload.get("barrier_kcal")
                if g_payload.get("delta_kcal") is not None:
                    entry["gibbs_mlip_delta"] = g_payload.get("delta_kcal")
                for endpoint in (1, 2):
                    value = g_payload.get(
                        f"barrier_from_endpoint_{endpoint}_kcal"
                    )
                    if value is not None:
                        entry[f"gibbs_mlip_barrier_e{endpoint}"] = value
            if seg.get("dft"):
                dft_payload = seg.get("dft") or {}
                if dft_payload.get("barrier_kcal") is not None:
                    entry["dft_barrier"] = dft_payload.get("barrier_kcal")
                if dft_payload.get("delta_kcal") is not None:
                    entry["dft_delta"] = dft_payload.get("delta_kcal")
                for endpoint in (1, 2):
                    value = dft_payload.get(
                        f"barrier_from_endpoint_{endpoint}_kcal"
                    )
                    if value is not None:
                        entry[f"dft_barrier_e{endpoint}"] = value
            if seg.get("gibbs_dft_mlip"):
                gd_payload = seg.get("gibbs_dft_mlip") or {}
                if gd_payload.get("barrier_kcal") is not None:
                    entry["gibbs_dft_mlip_barrier"] = gd_payload.get("barrier_kcal")
                if gd_payload.get("delta_kcal") is not None:
                    entry["gibbs_dft_mlip_delta"] = gd_payload.get("delta_kcal")
                for endpoint in (1, 2):
                    value = gd_payload.get(
                        f"barrier_from_endpoint_{endpoint}_kcal"
                    )
                    if value is not None:
                        entry[f"gibbs_dft_mlip_barrier_e{endpoint}"] = value
    else:
        lines.append("  (no post-processing results)")

    if segment_entries:
        table_rows = [
            (f"MEP {delta}E{dagger} [kcal/mol]", "mep_barrier"),
            (f"MEP {delta}E  [kcal/mol]", "mep_delta"),
            (f"ML/MM {delta}E{dagger} [kcal/mol]", "mlip_barrier"),
            (f"ML/MM {delta}E  [kcal/mol]", "mlip_delta"),
            (f"ML/MM {delta}G{dagger} [kcal/mol]", "gibbs_mlip_barrier"),
            (f"ML/MM {delta}G  [kcal/mol]", "gibbs_mlip_delta"),
            (f"model-region DFT {delta}E{dagger} [kcal/mol]", "dft_barrier"),
            (f"model-region DFT {delta}E  [kcal/mol]", "dft_delta"),
            (f"DFT//MLIP/MM {delta}G{dagger} [kcal/mol]", "gibbs_dft_mlip_barrier"),
            (f"DFT//MLIP/MM {delta}G  [kcal/mol]", "gibbs_dft_mlip_delta"),
        ]
        if ts_only:
            table_rows = [
                (f"ML/MM {delta}E{dagger} E1->TS [kcal/mol]", "mlip_barrier_e1"),
                (f"ML/MM {delta}E{dagger} E2->TS [kcal/mol]", "mlip_barrier_e2"),
                (f"ML/MM {delta}G{dagger} E1->TS [kcal/mol]", "gibbs_mlip_barrier_e1"),
                (f"ML/MM {delta}G{dagger} E2->TS [kcal/mol]", "gibbs_mlip_barrier_e2"),
                (f"model-region DFT {delta}E{dagger} E1->TS [kcal/mol]", "dft_barrier_e1"),
                (f"model-region DFT {delta}E{dagger} E2->TS [kcal/mol]", "dft_barrier_e2"),
                (f"DFT//MLIP/MM {delta}G{dagger} E1->TS [kcal/mol]", "gibbs_dft_mlip_barrier_e1"),
                (f"DFT//MLIP/MM {delta}G{dagger} E2->TS [kcal/mol]", "gibbs_dft_mlip_barrier_e2"),
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
            (f"ML/MM {delta}E  [kcal/mol]", "mlip"),
            (f"ML/MM {delta}G  [kcal/mol]", "gibbs_mlip"),
            (f"model-region DFT {delta}E  [kcal/mol]", "dft"),
            (f"DFT//MLIP/MM {delta}G  [kcal/mol]", "gibbs_dft_mlip"),
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

    # Annotations follow the systematized `all` layout: user-facing deliverables
    # live at the output root (and under segments/seg_NN/), while all pipeline
    # scratch is confined to _work/ (safe to rm -rf).
    state_triplet = "E1-TS-E2" if ts_only else "R-TS-P"
    default_notes = {
        # Root deliverables — directories
        SEGMENTS_DIRNAME: f"Per-segment deliverables ({state_triplet}, IRC, freq, DFT)",
        "mm_parm": "AMBER MM topology (parm7/rst7) — reuse via --parm",
        "layered": "B-factor-layered PDBs (ML/MM region markup) for inspection/reuse",
        WORK_DIRNAME: "Pipeline working files (scratch; safe to rm -rf)",
        # Root deliverables — files
        "summary.json": "Machine-readable results (JSON)",
        "summary.log": "Human-readable results summary",
        "ml_region.pdb": "ML-region definition — reuse via --model-pdb",
        "ml_region_without_linkH.xyz": "ML-region coordinates without link H",
        "ml_region_with_linkH.xyz": "ML-region coordinates with parm7-derived link H",
        "mep.pdb": "Full MEP as single PDB (all segments)",
        "mep_trj.xyz": "Full MEP as XYZ trajectory",
        "mep_plot.png": "ML/MM MEP energy plot",
        "energy_diagram_MEP.png": "Compressed MEP diagram",
        "energy_diagram_MLIP_all.png": f"ML/MM {state_triplet} energies (all segments)",
        "energy_diagram_G_MLIP_all.png": f"ML/MM Gibbs {state_triplet} (all segments)",
        "energy_diagram_DFT_all.png": f"DFT {state_triplet} (all segments)",
        "energy_diagram_G_DFT_plus_MLIP_all.png": f"DFT//MLIP/MM Gibbs {state_triplet} (all segments)",
        "irc_plot_all.png": "Aggregated IRC plot",
        # _work/ scratch subdirectories
        f"{WORK_DIRNAME}/pockets": "Extracted pocket PDBs",
        f"{WORK_DIRNAME}/add_elem_info": "Inputs with element fields repaired (preflight)",
        f"{WORK_DIRNAME}/scan": "Staged scan outputs",
    }

    selected_current = payload.get("current_output_paths")
    if root_out_path and isinstance(selected_current, (list, tuple)):
        selected_relative = [Path(str(path)).as_posix() for path in selected_current]
        for rel in selected_relative:
            note = default_notes.get(rel)
            if note:
                annotations.setdefault(rel, note)
        lines.extend(
            _format_selected_directory_tree(
                root_out_path,
                selected_relative,
                annotations,
            )
        )
    elif root_out_path:
        path_dir = payload.get("path_dir")
        if path_dir:
            try:
                rel = Path(path_dir).relative_to(root_out_path).as_posix()
                annotations.setdefault(rel, "Primary path module outputs (scratch)")
            except ValueError:
                pass

        for rel, desc in default_notes.items():
            if (root_out_path / rel).exists():
                annotations.setdefault(rel, desc)

        import re as _re

        # Per-segment deliverables under segments/seg_NN/
        seg_parent = root_out_path / SEGMENTS_DIRNAME
        if seg_parent.exists():
            _seg_subdir_notes = {
                "structures": "State structures and IRC endpoints",
                "irc": "IRC trajectories (forward/backward/finished)",
                "ts": "TS optimization output",
                "ts/vib": "TS vibrational analysis",
                "freq": "Frequency and thermochemistry",
                "freq/R": "Reactant freq/thermo",
                "freq/TS": "TS freq/thermo",
                "freq/P": "Product freq/thermo",
                "freq/E1": "Endpoint 1 freq/thermo",
                "freq/E2": "Endpoint 2 freq/thermo",
                "dft": "Single-point DFT refinement",
                "dft/R": "Reactant DFT single point",
                "dft/TS": "TS DFT single point",
                "dft/P": "Product DFT single point",
                "dft/E1": "Endpoint 1 DFT single point",
                "dft/E2": "Endpoint 2 DFT single point",
            }
            for seg_child in sorted(seg_parent.iterdir()):
                if not (seg_child.is_dir() and _re.match(r"seg_\d+$", seg_child.name)):
                    continue
                seg_num = seg_child.name.replace("seg_", "")
                seg_rel = f"{SEGMENTS_DIRNAME}/{seg_child.name}"
                annotations.setdefault(
                    seg_rel,
                    f"Refined TS and optimized IRC endpoints of segment {seg_num}",
                )
                for subdir_name, desc in _seg_subdir_notes.items():
                    if (seg_child / subdir_name).exists():
                        annotations.setdefault(f"{seg_rel}/{subdir_name}", desc)

        # Dynamic annotations for path module internal directories
        if path_dir:
            try:
                path_dir_path = Path(path_dir)
                if path_dir_path.exists():
                    path_dir_rel = path_dir_path.relative_to(root_out_path).as_posix()
                    for child in sorted(path_dir_path.iterdir()):
                        if not child.is_dir():
                            # Annotate key files inside path module dir
                            crel = f"{path_dir_rel}/{child.name}"
                            if _re.match(r"hei_seg_\d+\.", child.name):
                                annotations.setdefault(crel, "Highest-energy image (approx. TS)")
                            elif _re.match(r"hei_w_ref_seg_\d+\.", child.name):
                                annotations.setdefault(crel, "HEI with protein reference frame")
                            elif _re.match(r"mep_seg_\d+_trj\.", child.name):
                                annotations.setdefault(crel, "Per-segment MEP trajectory")
                            elif _re.match(r"mep_w_ref_seg_\d+\.", child.name):
                                annotations.setdefault(crel, "Per-segment MEP with protein reference")
                            continue

                        crel = f"{path_dir_rel}/{child.name}"

                        # post_seg_XX/ directories
                        if child.name.startswith("post_seg_"):
                            annotations.setdefault(
                                crel,
                                f"Post-processing: TSOPT, IRC, freq for {child.name}",
                            )
                            _subdir_notes = {
                                "structures": "Optimized state structures (IRC endpoints)",
                                "irc": "IRC trajectories and plots",
                                "ts": "TS optimization output",
                                "ts/vib": "TS vibrational analysis",
                                "freq": "Frequency and thermochemistry",
                                "freq/R": "Reactant freq/thermo",
                                "freq/TS": "TS freq/thermo",
                                "freq/P": "Product freq/thermo",
                                "freq/E1": "Endpoint 1 freq/thermo",
                                "freq/E2": "Endpoint 2 freq/thermo",
                            }
                            for subdir_name, desc in _subdir_notes.items():
                                sub = child / subdir_name
                                if sub.exists():
                                    annotations.setdefault(f"{crel}/{subdir_name}", desc)
                        # init optimization dirs
                        elif _re.match(r"init\d+_lbfgs_opt$", child.name):
                            idx = _re.search(r"init(\d+)", child.name).group(1)
                            annotations.setdefault(crel, f"Initial optimization of endpoint {idx}")
                        # GSM/NEB path dirs
                        elif _re.match(r"seg_\d+_mep$", child.name):
                            annotations.setdefault(crel, "Initial GSM/NEB path")
                        elif _re.match(r"seg_\d+_refine_mep$", child.name):
                            annotations.setdefault(crel, "Refined GSM/NEB path")
                        elif _re.match(r"seg_\d+_left_lbfgs_opt$", child.name):
                            annotations.setdefault(crel, "Optimized left (R) endpoint")
                        elif _re.match(r"seg_\d+_right_lbfgs_opt$", child.name):
                            annotations.setdefault(crel, "Optimized right (P) endpoint")
                        elif _re.search(r"_bridge_mep$", child.name):
                            annotations.setdefault(crel, "Bridge MEP (non-reactive conformational change)")
                        elif child.name == "align_refine":
                            annotations.setdefault(crel, "Pair-wise endpoint alignment")
            except (ValueError, OSError):
                pass

        if root_out_path.exists():
            lines.extend(_format_directory_tree(root_out_path, annotations))
        else:
            lines.append("  (root output directory not found on disk)")
    else:
        lines.append("  (root output directory unknown)")

    lines.append("")
    lines.extend(format_method_citations(payload))

    from mlmm.core.result_commit import commit_exact_bytes

    commit_exact_bytes(dest, ("\n".join(lines) + "\n").encode("utf-8"))
