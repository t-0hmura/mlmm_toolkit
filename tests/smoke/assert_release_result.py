#!/usr/bin/env python3
"""Strict semantic checks for mlmm release-smoke outputs."""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import numpy as np

from mlmm.core.utils import mlip_model_label


def require_finite(value: object, path: str = "root") -> None:
    if isinstance(value, bool) or value is None or isinstance(value, str):
        return
    if isinstance(value, (int, float)):
        if not math.isfinite(float(value)):
            raise SystemExit(f"non-finite value at {path}")
        return
    if isinstance(value, dict):
        for key, child in value.items():
            require_finite(child, f"{path}.{key}")
    elif isinstance(value, list):
        for index, child in enumerate(value):
            require_finite(child, f"{path}[{index}]")


def count_xyz_frames(path: Path) -> int:
    lines = path.read_text(encoding="utf-8").splitlines()
    index = 0
    frames = 0
    while index < len(lines):
        try:
            atoms = int(lines[index].strip())
        except ValueError as exc:
            raise SystemExit(f"invalid XYZ frame in {path}") from exc
        index += atoms + 2
        frames += 1
    if index != len(lines):
        raise SystemExit(f"truncated XYZ trajectory: {path}")
    return frames


def check_dft_states(segment_root: Path, segment: dict) -> None:
    """Verify the requested CPU HF/STO-3G SCF results, not artifact names."""
    import yaml

    dft = segment.get("dft") or {}
    energies = dft.get("energies_au") or []
    if dft.get("status") == "failed" or len(energies) != 3 or not all(
        isinstance(value, (int, float)) and not isinstance(value, bool)
        and math.isfinite(float(value)) for value in energies
    ):
        raise SystemExit(f"DFT R/TS/P energies are incomplete: {dft!r}")
    states = ("E1", "TS", "E2") if segment.get("kind") == "tsopt" else ("R", "TS", "P")
    for state in states:
        path = segment_root / "dft" / state / "result.yaml"
        data = yaml.safe_load(path.read_text(encoding="utf-8"))
        energy = data.get("energy") or {}
        value = energy.get("hartree")
        if energy.get("converged") is not True or isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value):
            raise SystemExit(f"DFT {state} has no finite converged SCF result: {energy!r}")
        request = data.get("input") or {}
        if str(request.get("xc")).lower() != "hf" or str(request.get("basis")).lower() != "sto-3g":
            raise SystemExit(f"DFT {state} did not use requested HF/STO-3G: {request!r}")
        if request.get("max_cycle") != 40 or request.get("grid_level") != 0 or request.get("conv_tol") != 1e-5:
            raise SystemExit(f"DFT {state} did not preserve requested SCF settings: {request!r}")
        if energy.get("used_gpu") is not False:
            raise SystemExit(f"DFT {state} did not use the requested CPU calculation")


def check_all(root: Path, require_thermo: bool, require_dft: bool) -> None:
    summary = json.loads((root / "summary.json").read_text(encoding="utf-8"))
    require_finite(summary)

    # The pipeline must have run to the end.  This is the part a smoke can
    # guarantee, and the part a broken pipeline fails.
    execution = summary.get("execution_status")
    if execution != "completed":
        raise SystemExit(f"all pipeline did not complete: execution_status={execution!r}")

    # The test73 lane runs with strict convergence criteria (`gau` optimizer
    # and `baker` TS verification) and requires the reported outcome to be TRUE:
    # either a real `success`, or a `partial` that clearly states what is missing.
    # A silent degradation, a `failed`, or a `partial` with no stated reason will
    # fail the lane -- this ensures the reported results accurately reflect the
    # actual computed state.
    status = summary.get("status")
    scientific = summary.get("scientific_status")
    if status not in ("success", "partial") or scientific not in ("success", "partial"):
        raise SystemExit(
            f"all summary is neither success nor an explained partial: "
            f"status={status!r} scientific_status={scientific!r}"
        )
    reasons = summary.get("scientific_status_reasons") or []
    if scientific == "partial" and not reasons:
        raise SystemExit("all summary is partial but states no reason")
    if scientific == "success" and reasons:
        raise SystemExit(f"all summary claims success yet states reasons: {reasons}")
    if any("irc:" in str(reason) and "not_converged" in str(reason) for reason in reasons):
        raise SystemExit(f"normal IRC stopping was misreported as nonconvergence: {reasons}")

    segments = summary.get("post_segments") or []
    if not segments:
        raise SystemExit("all summary has no post-TS segments")
    for segment in segments:
        if (segment.get("tsopt") or {}).get("optimization_status") != "converged":
            raise SystemExit("TS numerical optimization did not converge")
        tag = str(segment["tag"])
        # Use the directory the producer published for this segment. ``tag`` is
        # the path-search segment id (``seg_%03d``) and is NOT a directory name;
        # the per-segment deliverable dir is indexed separately (``seg_%02d``).
        # Rebuilding the path from ``tag`` asserts an invariant the product has
        # never promised, and fails on the very first segment.
        published = segment.get("post_dir")
        if not published:
            raise SystemExit(f"segment {tag!r} published no post_dir")
        segment_root = Path(published)
        if not segment_root.is_dir():
            raise SystemExit(f"published post_dir is not a directory: {segment_root}")
        for direction in ("forward", "backward"):
            trajectory = segment_root / "irc" / f"{direction}_irc_trj.xyz"
            if not trajectory.is_file() or count_xyz_frames(trajectory) < 2:
                raise SystemExit(f"missing/nontrivial IRC branch: {trajectory}")
        irc = segment.get("irc") or {}
        if {"usable", "scientific_status", "forward_status", "backward_status"} & irc.keys():
            raise SystemExit(f"all retained an independent IRC acceptance field: {irc!r}")
        endpoint_opt = segment.get("endpoint_opt") or {}
        endpoint_keys = (
            ("endpoint_1_converged", "endpoint_2_converged")
            if segment.get("kind") == "tsopt"
            else ("reactant_converged", "product_converged")
        )
        if any(endpoint_opt.get(key) is not True for key in endpoint_keys):
            raise SystemExit(f"optimized endpoints did not converge: {endpoint_opt!r}")
        if require_thermo:
            missing = [
                state
                for state in ("R", "TS", "P")
                if not (segment_root / "freq" / state / "thermoanalysis.yaml").is_file()
                or (segment_root / "freq" / state / "thermoanalysis.yaml").stat().st_size == 0
            ]
            if missing:
                raise SystemExit(f"thermochemistry missing for {missing}")
            if not isinstance(segment.get("gibbs_mlip"), dict):
                raise SystemExit(f"MLIP/MM thermochemistry is missing for {tag}")
        if require_dft:
            check_dft_states(segment_root, segment)
            if require_thermo and not isinstance(segment.get("gibbs_dft_mlip"), dict):
                raise SystemExit(f"DFT//MLIP thermochemistry is missing for {tag}")


def check_opt_config(
    root: Path,
    *,
    expected_model: str,
    expected_max_cycles: int,
    expected_precision: str,
    expected_link_atom_method: str,
    expected_thresh: str,
) -> None:
    payload = json.loads((root / "result.json").read_text(encoding="utf-8"))
    require_finite(payload)
    if payload.get("command") != "opt":
        raise SystemExit(f"unexpected result command: {payload.get('command')!r}")
    if payload.get("mlip_model") != expected_model:
        raise SystemExit(f"runtime model did not come from YAML: {payload.get('mlip_model')!r}")
    if int(payload.get("max_cycles", -1)) != expected_max_cycles:
        raise SystemExit(f"CLI max-cycles did not override YAML: {payload.get('max_cycles')!r}")
    if payload.get("mlip_precision") != expected_precision:
        raise SystemExit(f"runtime precision did not come from YAML: {payload.get('mlip_precision')!r}")
    if payload.get("link_atom_method") != expected_link_atom_method:
        raise SystemExit(
            "runtime link-atom method did not come from YAML: "
            f"{payload.get('link_atom_method')!r}"
        )
    if payload.get("thresh") != expected_thresh:
        raise SystemExit(f"runtime threshold did not come from YAML: {payload.get('thresh')!r}")


def check_sp_hessian(root: Path) -> None:
    payload = json.loads((root / "result.json").read_text(encoding="utf-8"))
    require_finite(payload)
    if payload.get("status") != "ok" or payload.get("mlip_backend") != "custom":
        raise SystemExit("custom calculator SP did not complete successfully")
    if payload.get("mlip_model") != "harmonic_calc.py:get_calculator":
        raise SystemExit(f"custom calculator provenance is wrong: {payload.get('mlip_model')!r}")
    if payload.get("mlip_precision") is not None:
        raise SystemExit("custom calculator incorrectly reports a built-in MLIP precision")
    hessian = np.load(root / "hessian.npy", allow_pickle=False)
    forces = np.load(root / "forces.npy", allow_pickle=False)
    if forces.ndim != 2 or forces.shape[1] != 3 or not np.isfinite(forces).all():
        raise SystemExit(f"custom calculator forces are invalid: {forces.shape}")
    if Path(str(payload.get("forces_path", ""))).resolve() != (root / "forces.npy").resolve():
        raise SystemExit("custom calculator result does not identify forces.npy")
    if Path(str(payload.get("hessian_path", ""))).resolve() != (root / "hessian.npy").resolve():
        raise SystemExit("custom calculator result does not identify hessian.npy")
    full_size = 3 * int(np.atleast_2d(forces).shape[0])
    if hessian.ndim != 2 or hessian.shape[0] != hessian.shape[1]:
        raise SystemExit(f"Hessian is not square: {hessian.shape}")
    if not 0 < hessian.shape[0] <= full_size:
        raise SystemExit(f"Hessian size {hessian.shape[0]} is incompatible with {full_size} Cartesian DOFs")
    if not np.isfinite(hessian).all():
        raise SystemExit("custom finite-difference Hessian contains non-finite values")
    if not np.allclose(hessian, hessian.T, atol=1.0e-8, rtol=0.0):
        raise SystemExit("custom finite-difference Hessian is not symmetric")
    if not np.any(np.abs(hessian) > 1.0e-8):
        raise SystemExit("custom finite-difference Hessian is identically zero")


def read_xyz_coords(path: Path) -> "list[np.ndarray]":
    """Return one (n_atoms, 3) array of Angstrom coordinates per XYZ frame."""
    lines = path.read_text(encoding="utf-8").splitlines()
    frames: list = []
    index = 0
    while index < len(lines):
        try:
            atoms = int(lines[index].strip())
        except ValueError as exc:
            raise SystemExit(f"invalid XYZ frame in {path}") from exc
        rows = []
        for offset in range(atoms):
            parts = lines[index + 2 + offset].split()
            if len(parts) < 4:
                raise SystemExit(f"invalid XYZ atom line in {path}")
            rows.append([float(value) for value in parts[1:4]])
        frames.append(np.asarray(rows, dtype=float))
        index += atoms + 2
    if index != len(lines):
        raise SystemExit(f"truncated XYZ trajectory: {path}")
    return frames


def check_dmf_frozen_atoms(root: Path, frozen_1based: str) -> None:
    """Assert the DMF harmonic restraint pinned the requested atoms.

    The rest of the release smoke never enters the DMF ``fix_atoms`` branch, so
    this lane is its only coverage. The check compares the optimized path against
    the FB-ENM interpolation the per-image restraint references, which keeps it
    independent of the MLIP energies and of whether IPOPT converged.
    """
    payload = json.loads((root / "result.json").read_text(encoding="utf-8"))
    require_finite(payload)
    energies = payload.get("image_energies_hartree")
    if not isinstance(energies, list) or not energies:
        raise SystemExit("DMF run published no image_energies_hartree")

    initial = read_xyz_coords(root / "dmf_initial_trj.xyz")
    final = read_xyz_coords(root / "final_geometries_trj.xyz")
    if not final or len(initial) != len(final):
        raise SystemExit(
            f"DMF frame counts disagree: {len(initial)} interpolated vs {len(final)} final"
        )

    frozen = sorted({int(token) - 1 for token in frozen_1based.split(",") if token.strip()})
    if not frozen:
        raise SystemExit("dmf-freeze requires at least one frozen atom index")
    n_atoms = final[0].shape[0]
    if any(index < 0 or index >= n_atoms for index in frozen):
        raise SystemExit(f"frozen indices outside a {n_atoms}-atom geometry: {frozen}")
    free = [index for index in range(n_atoms) if index not in set(frozen)]
    if not free:
        raise SystemExit("dmf-freeze needs at least one unfrozen atom")

    frozen_drift = 0.0
    free_drift = 0.0
    for before, after in zip(initial, final):
        if before.shape != after.shape:
            raise SystemExit("DMF interpolated and final frames differ in atom count")
        displacement = np.linalg.norm(after - before, axis=1)
        frozen_drift = max(frozen_drift, float(displacement[frozen].max()))
        free_drift = max(free_drift, float(displacement[free].max()))

    if frozen_drift > 0.05:
        raise SystemExit(
            f"frozen atoms drifted {frozen_drift:.4f} A from their DMF reference"
        )
    # Positive control: a run in which nothing moved would satisfy the bound above
    # without exercising the restraint at all.
    if free_drift <= frozen_drift:
        raise SystemExit(
            f"no unfrozen atom moved more than the frozen ones "
            f"(free {free_drift:.4f} A, frozen {frozen_drift:.4f} A); the lane proves nothing"
        )


def check_irc_direction_status_contract(payload: dict) -> None:
    """IRC retains stop diagnostics and candidates, not a scientific verdict."""
    for removed in ("forward_converged", "backward_converged",
                    "forward_endpoint_stationary", "backward_endpoint_stationary",
                    "forward_status", "backward_status", "scientific_status", "stage_outcomes"):
        if removed in payload:
            raise SystemExit(f"IRC still publishes an independent acceptance field: {removed}")
    requested = [d for d in ("forward", "backward") if payload.get(f"{d}_requested") is True]
    if not requested:
        raise SystemExit("IRC reported no requested direction")
    for direction in requested:
        if int(payload.get(f"n_frames_{direction}", 0)) < 1:
            raise SystemExit(f"{direction} produced no retained candidate frames")
        if f"{direction}_integration_stop_reason" not in payload:
            raise SystemExit(f"{direction} stop diagnostics are missing")


def check_irc_handoff(root: Path, hessian_file: Path) -> None:
    payload = json.loads((root / "result.json").read_text(encoding="utf-8"))
    require_finite(payload)
    projection = payload.get("rigid_projection") or {}
    if payload.get("status") != "completed":
        raise SystemExit(f"IRC did not complete: {payload.get('status')!r}")
    if projection.get("hessian_source") != "file":
        raise SystemExit(f"IRC did not consume --read-hess: {projection.get('hessian_source')!r}")
    if payload.get("never_stop") is not True:
        raise SystemExit("IRC never-stop flag did not reach runtime result")
    if int(payload.get("never_stop_energy_bypasses", 0)) < 1:
        raise SystemExit("IRC never-stop did not bypass an actual energy stop")
    if (
        payload.get("forward_integration_converged") is not False
        or payload.get("backward_integration_converged") is not False
    ):
        raise SystemExit(
            "IRC never-stop reported a stationarity stop it cannot reach"
        )
    check_irc_direction_status_contract(payload)
    if int(payload.get("n_frames_forward", 0)) < 2 or int(payload.get("n_frames_backward", 0)) < 2:
        raise SystemExit("IRC Hessian handoff did not produce both nontrivial branches")
    freq_result_path = hessian_file.parent / "result.json"
    freq_payload = json.loads(freq_result_path.read_text(encoding="utf-8"))
    reported_hessian = (freq_payload.get("files") or {}).get("hessian_npz")
    if Path(str(reported_hessian)).resolve() != hessian_file.resolve():
        raise SystemExit(
            f"freq result does not identify its Hessian dump: {reported_hessian!r}"
        )
    with np.load(hessian_file, allow_pickle=False) as data:
        active = np.asarray(data["wph_active_dofs"], dtype=int).reshape(-1)
        hessian = np.asarray(data["hessian"], dtype=float)
        full_n_dof = int(data["wph_full_n_dof"])
        coords = np.asarray(data["cart_coords_bohr"], dtype=float).reshape(-1)
        if not (active.size == hessian.shape[0] < full_n_dof == coords.size):
            raise SystemExit("freq dump was not a self-consistent partial Hessian")
        if list(projection.get("hessian_shape") or []) != list(hessian.shape):
            raise SystemExit("IRC result does not report the dumped partial-Hessian shape")


def check_tsopt_reference(root: Path) -> None:
    payload = json.loads((root / "result.json").read_text(encoding="utf-8"))
    require_finite(payload)
    safeguards = payload.get("safeguards") or {}
    if not payload.get("reference_mode_file"):
        raise SystemExit("TS optimization result omitted the reference-mode file")
    if payload.get("opt_mode_requested") != "hess" or payload.get("optimizer") != "rsprfo":
        raise SystemExit("TS optimization requested/effective optimizer provenance is wrong")
    index = safeguards.get("initial_reference_root_index")
    overlap = safeguards.get("initial_reference_root_overlap")
    eigenvalue = safeguards.get("initial_reference_root_eigenvalue")
    if not isinstance(index, int) or index < 0:
        raise SystemExit(f"reference-mode root selection was not recorded: {index!r}")
    if overlap is None or not 0.0 < float(overlap) <= 1.0 + 1.0e-12:
        raise SystemExit(f"invalid initial reference-root overlap: {overlap!r}")
    if eigenvalue is None or not math.isfinite(float(eigenvalue)):
        raise SystemExit(f"invalid initial reference-root eigenvalue: {eigenvalue!r}")


def check_tsopt_optimizer(root: Path, expected_mode: str, expected_optimizer: str) -> None:
    """Assert TS requested/effective optimizer provenance without requiring convergence.

    Bound to a deliberate max-cycle lane: the run reports its optimizer identity
    whether or not it converged, so this check must not gate on convergence.
    """
    payload = json.loads((root / "result.json").read_text(encoding="utf-8"))
    require_finite(payload)
    if payload.get("opt_mode_requested") != expected_mode:
        raise SystemExit(f"unexpected TS preset: {payload.get('opt_mode_requested')!r}")
    if payload.get("optimizer") != expected_optimizer:
        raise SystemExit(f"unexpected TS optimizer: {payload.get('optimizer')!r}")


def check_scan_optimizer(root: Path, expected_mode: str, expected_optimizer: str) -> None:
    payload = json.loads((root / "result.json").read_text(encoding="utf-8"))
    require_finite(payload)
    if payload.get("scan_opt_mode") != expected_mode:
        raise SystemExit(f"unexpected scan preset: {payload.get('scan_opt_mode')!r}")
    if payload.get("scan_optimizer") != expected_optimizer:
        raise SystemExit(f"unexpected scan optimizer: {payload.get('scan_optimizer')!r}")


def check_provenance(
    root: Path,
    *,
    expected_backend: str,
    expected_model: str,
    expected_precision: str,
) -> None:
    summary = json.loads((root / "summary.json").read_text(encoding="utf-8"))
    expected = {
        "mlip_backend": expected_backend,
        "mlip_model": expected_model,
        "mlip_precision": expected_precision,
    }
    for key, value in expected.items():
        if summary.get(key) != value:
            raise SystemExit(
                f"incorrect {key} provenance: {summary.get(key)!r}, expected {value!r}"
            )
    text = (root / "summary.log").read_text(encoding="utf-8")
    backend_label = {
        "uma": "UMA", "orb": "ORB", "mace": "MACE", "aimnet2": "AIMNet2",
    }.get(expected_backend, expected_backend)
    if f"MLIP backend        : {backend_label}" not in text:
        raise SystemExit("summary.log omits the resolved MLIP backend")
    model_label = summary.get("mlip_model_label") or mlip_model_label(
        expected_backend, expected_model, summary.get("mlip_task")
    )
    if f"MLIP model          : {model_label}" not in text:
        raise SystemExit("summary.log omits the resolved MLIP model")
    for stale in ("UMA model", "DFT//UMA", "energy_diagram_UMA"):
        if stale in text:
            raise SystemExit(f"summary.log contains stale backend label: {stale}")


def check_path_search_max_depth(root: Path) -> None:
    """Check that max_depth=0 retains one input interval without subdivision."""
    data = json.loads((root / "summary.json").read_text(encoding="utf-8"))
    n_segments = data.get("n_segments")
    segments = data.get("segments")
    if not isinstance(n_segments, int) or n_segments not in (0, 1):
        raise SystemExit(f"Invalid n_segments={n_segments!r}")
    if not isinstance(segments, list) or len(segments) != n_segments:
        raise SystemExit("Segment records disagree with n_segments")
    for segment in segments:
        if not segment.get("tag") or segment.get("barrier_kcal") is None:
            raise SystemExit("A reported segment lacks its tag or energetics")
    if n_segments == 0:
        reasons = " ".join(str(r) for r in (data.get("scientific_status_reasons") or []))
        if "endpoint_hei" not in reasons:
            raise SystemExit("No segment and no endpoint-HEI diagnostic")
    if data.get("search_max_depth") != 0:
        raise SystemExit("The effective zero recursion cap was not recorded")
    if not any(path.is_dir() for path in root.glob("seg_*_mep")):
        raise SystemExit("The raw MEP interval is missing")
    text = (root / "summary.log").read_text(encoding="utf-8")
    if "Recursion depth cap : 0 (subdivision disabled)" not in text:
        raise SystemExit("summary.log omits the disabled-subdivision setting")
    if any(str(segment.get("tag", "")).endswith("_maxdepth") for segment in segments):
        raise SystemExit("Disabled subdivision must retain the plain input-interval tag")
    optimizers = data.get("path_optimizers", [])
    if "lbfgs" in optimizers or "rfo" in optimizers:
        raise SystemExit("No single-structure optimizer is requested in this no-preopt case")
    if ("lbfgs" in optimizers) != ("Limited-memory BFGS (L-BFGS)" in text):
        raise SystemExit("L-BFGS citation does not match the executed optimizer record")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "kind",
        choices=("all", "tsopt-optimizer", "scan-optimizer", "dmf-freeze", "opt-config", "sp-hessian", "irc-handoff", "tsopt-reference", "path-search-max-depth", "provenance"),
    )
    parser.add_argument("root", type=Path)
    parser.add_argument("--require-thermo", action="store_true")
    parser.add_argument("--require-dft", action="store_true")
    parser.add_argument("--expected-model")
    parser.add_argument("--expected-max-cycles", type=int)
    parser.add_argument("--expected-precision")
    parser.add_argument("--expected-backend")
    parser.add_argument("--expected-mode")
    parser.add_argument("--frozen-atoms")
    parser.add_argument("--expected-optimizer")
    parser.add_argument("--expected-link-atom-method")
    parser.add_argument("--expected-thresh")
    parser.add_argument("--hessian-file", type=Path)
    args = parser.parse_args()
    if args.kind == "all":
        check_all(args.root, args.require_thermo, args.require_dft)
    elif args.kind == "dmf-freeze":
        if args.frozen_atoms is None:
            parser.error("dmf-freeze requires --frozen-atoms")
        check_dmf_frozen_atoms(args.root, args.frozen_atoms)
    elif args.kind == "tsopt-optimizer":
        if None in (args.expected_mode, args.expected_optimizer):
            parser.error("tsopt-optimizer requires --expected-mode and --expected-optimizer")
        check_tsopt_optimizer(args.root, args.expected_mode, args.expected_optimizer)
    elif args.kind == "scan-optimizer":
        if None in (args.expected_mode, args.expected_optimizer):
            parser.error("scan-optimizer requires --expected-mode and --expected-optimizer")
        check_scan_optimizer(args.root, args.expected_mode, args.expected_optimizer)
    elif args.kind == "opt-config":
        if None in (
            args.expected_model,
            args.expected_max_cycles,
            args.expected_precision,
            args.expected_link_atom_method,
            args.expected_thresh,
        ):
            parser.error(
                "opt-config requires --expected-model, --expected-max-cycles, "
                "--expected-precision, --expected-link-atom-method, and --expected-thresh"
            )
        check_opt_config(
            args.root,
            expected_model=args.expected_model,
            expected_max_cycles=args.expected_max_cycles,
            expected_precision=args.expected_precision,
            expected_link_atom_method=args.expected_link_atom_method,
            expected_thresh=args.expected_thresh,
        )
    elif args.kind == "path-search-max-depth":
        check_path_search_max_depth(args.root)
    elif args.kind == "sp-hessian":
        check_sp_hessian(args.root)
    elif args.kind == "irc-handoff":
        if args.hessian_file is None:
            parser.error("irc-handoff requires --hessian-file")
        check_irc_handoff(args.root, args.hessian_file)
    elif args.kind == "tsopt-reference":
        check_tsopt_reference(args.root)
    else:
        if None in (args.expected_backend, args.expected_model, args.expected_precision):
            parser.error(
                "provenance requires --expected-backend, --expected-model, "
                "and --expected-precision"
            )
        check_provenance(
            args.root,
            expected_backend=args.expected_backend,
            expected_model=args.expected_model,
            expected_precision=args.expected_precision,
        )


if __name__ == "__main__":
    main()
