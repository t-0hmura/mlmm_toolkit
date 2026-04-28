# mlmm/defaults.py

"""
Central configuration defaults for mlmm workflows.

All default dictionaries are defined here to avoid redundant definitions across modules.
Modules should import defaults from here instead of defining local copies.

Shared optimizer/IRC/path defaults — keep aligned with pdb2reaction/defaults.py.
Project-specific overrides are commented with rationale.
"""

from typing import Any, Dict

# -----------------------------------------------
# B-factor values for 3-layer ML/MM system
# -----------------------------------------------

# B-factor values encode atom layer membership in PDB files:
#   ML atoms: 0.0
#   Movable MM atoms: 10.0
#   Frozen MM atoms: 20.0
#
# NOTE:
#   BFACTOR_HESS_MM uses the same value as movable MM for unified 3-layer encoding.
#   Hessian-target MM atoms are selected dynamically (e.g., by hess_cutoff),
#   not by a dedicated B-factor layer.
BFACTOR_ML = 0.0
BFACTOR_MOVABLE_MM = 10.0
BFACTOR_HESS_MM = BFACTOR_MOVABLE_MM
BFACTOR_FROZEN = 20.0

# Tolerance for B-factor matching (allow small variations)
BFACTOR_TOLERANCE = 1.0

# -----------------------------------------------
# Output directory defaults
# -----------------------------------------------

OUT_DIR_OPT = "./result_opt/"
OUT_DIR_SCAN = "./result_scan/"
OUT_DIR_SCAN2D = "./result_scan2d/"
OUT_DIR_SCAN3D = "./result_scan3d/"
OUT_DIR_FREQ = "./result_freq/"
OUT_DIR_IRC = "./result_irc/"
OUT_DIR_TSOPT = "./result_tsopt/"
OUT_DIR_PATH_OPT = "./result_path_opt/"
OUT_DIR_PATH_SEARCH = "./result_path_search/"

# -----------------------------------------------
# Geometry defaults
# -----------------------------------------------

GEOM_KW_DEFAULT: Dict[str, Any] = {
    "coord_type": "cart",
    "freeze_atoms": [],
}

# -----------------------------------------------
# Calculator defaults (ML/MM with MLIP backend + hessian_ff MM)
# -----------------------------------------------

MLMM_CALC_KW: Dict[str, Any] = {
    "input_pdb": None,
    "real_parm7": None,
    "model_pdb": None,
    "model_charge": 0,
    "model_mult": 1,
    "link_mlmm": None,
    "link_atom_method": "scaled",  # "scaled" (g-factor, Gaussian ONIOM standard) | "fixed" (1.09/1.01 Å legacy)
    # ML backend selection: "uma" | "orb" | "mace" | "aimnet2"
    "backend": "uma",
    "uma_model": "uma-s-1p1",
    "uma_task_name": "omol",
    "orb_model": "orb_v3_conservative_omol",
    "orb_precision": "float32",
    "mace_model": "MACE-OMOL-0",
    "mace_dtype": "float64",
    "aimnet2_model": "aimnet2",
    # ML Hessian mode: "FiniteDifference" or "Analytical"
    "hessian_calc_mode": "FiniteDifference",
    "out_hess_torch": True,
    "H_double": True,
    "ml_device": "auto",
    "ml_cuda_idx": 0,
    "mm_backend": "hessian_ff",  # "hessian_ff" (analytical) | "openmm" (FD Hessian)
    "use_cmap": False,           # If False, disable CMAP terms in model parm7 (Gaussian ONIOM-compatible)
    "mm_device": "cpu",
    "mm_cuda_idx": 0,
    "mm_threads": 16,
    "mm_fd": True,
    "mm_fd_dir": None,
    "mm_fd_delta": 1e-3,         # Displacement step for OpenMM FD Hessian (Å)
    "symmetrize_hessian": True,  # Symmetrize final Hessian as 0.5*(H+H^T)
    "print_timing": True,        # Print ML/MM Hessian timing breakdown
    "print_vram": True,          # Print CUDA VRAM usage (peak) during Hessian
    "return_partial_hessian": True,
    "freeze_atoms": [],
    # 3-layer B-factor configuration:
    #   ML region (B=0)
    #   Movable MM (B=10)
    #   Frozen MM (B=20)
    # Hessian-target MM is selected by hess_cutoff (distance from ML),
    # independent of B-factor layer.
    "hess_cutoff": None,      # Å, MM atoms within this distance of ML get MM Hessian (None = all movable)
    "movable_cutoff": None,   # Å, MM atoms within this distance of ML are movable (None = use freeze_atoms)
    # CLI defaults to --detect-layer for layer-aware commands.
    # Keep the base default aligned so config merge semantics stay consistent.
    "use_bfactor_layers": True,  # If True, read layer assignments from input PDB B-factors
    # Explicit YAML-based layer specification (1-based indices, takes precedence over cutoffs/B-factors)
    "hess_mm_atoms": None,    # Explicit Hessian-target MM atom indices
    "movable_mm_atoms": None, # Explicit movable MM atom indices
    "frozen_mm_atoms": None,  # Explicit frozen MM atom indices
    # xTB point-charge embedding correction
    "embedcharge": False,           # Enable xTB-based point-charge embedding
    "embedcharge_step": 1.0e-3,     # Numerical Hessian step for embedding correction (Å)
    "embedcharge_cutoff": 12.0,     # Distance cutoff (Å) for MM point charges in xTB embedding
    "xtb_cmd": "xtb",              # xTB executable command
    "xtb_acc": 0.2,                # xTB accuracy parameter
    "xtb_workdir": "tmp",          # xTB working directory
    "xtb_keep_files": False,       # Keep xTB temporary files
    "xtb_ncores": 4,               # Number of cores for xTB
}

# -----------------------------------------------
# Optimizer base (common to LBFGS & RFO)
# -----------------------------------------------

OPT_BASE_KW: Dict[str, Any] = {
    "thresh": "gau",
    "max_cycles": 10000,
    "print_every": 100,
    "min_step_norm": 1e-8,
    "assert_min_step": True,
    "rms_force": None,
    "rms_force_only": False,
    "max_force_only": False,
    "force_only": False,
    "converge_to_geom_rms_thresh": 0.05,
    "overachieve_factor": 0.0,
    "check_eigval_structure": False,
    "energy_plateau": True,
    "energy_plateau_thresh": 1e-4,
    "energy_plateau_window": 50,
    "line_search": True,
    "dump": False,
    "dump_restart": False,
    "prefix": "",
    "out_dir": OUT_DIR_OPT,
}

# -----------------------------------------------
# LBFGS-specific
# -----------------------------------------------

LBFGS_KW: Dict[str, Any] = {
    **OPT_BASE_KW,
    "keep_last": 7,
    "beta": 1.0,
    "gamma_mult": False,
    "max_step": 0.30,
    "control_step": True,
    "double_damp": True,
    "mu_reg": None,
    "max_mu_reg_adaptions": 10,
}

# -----------------------------------------------
# RFO-specific
# -----------------------------------------------

RFO_KW: Dict[str, Any] = {
    **OPT_BASE_KW,
    "trust_radius": 0.10,
    "trust_update": True,
    "trust_min": 1e-4,
    "trust_max": 0.10,
    "max_energy_incr": None,
    "hessian_update": "bfgs",
    "hessian_init": "calc",
    "hessian_recalc": 500,
    "hessian_recalc_adapt": None,
    "small_eigval_thresh": 1e-8,
    "alpha0": 1.0,
    "max_micro_cycles": 50,
    "rfo_overlaps": False,
    "gediis": False,
    "gdiis": True,
    "gdiis_thresh": 2.5e-3,
    "gediis_thresh": 1.0e-2,
    "gdiis_test_direction": True,
    "adapt_step_func": True,
}

# -----------------------------------------------
# Bias (harmonic well) defaults
# -----------------------------------------------

BIAS_KW: Dict[str, Any] = {
    "k": 300,
}

# -----------------------------------------------
# Bond-change detection
# -----------------------------------------------

BOND_KW: Dict[str, Any] = {
    "device": "auto",
    "bond_factor": 1.20,
    "margin_fraction": 0.05,
    "delta_fraction": 0.05,
}

# -----------------------------------------------
# Optimizer mode aliases
# -----------------------------------------------

OPT_MODE_ALIASES = (
    (("grad", "light", "lbfgs"), "lbfgs"),
    (("hess", "heavy", "rfo"), "rfo"),
)

# -----------------------------------------------
# DMF (Direct Max Flux) defaults for path optimization
# -----------------------------------------------

DMF_KW: Dict[str, Any] = {
    "max_cycles": 300,
    "correlated": True,
    "sequential": True,
    "fbenm_only_endpoints": False,
    "fbenm_options": {
        "delta_scale": 0.2,
        "bond_scale": 1.25,
        "fix_planes": True,
    },
    "cfbenm_options": {
        "bond_scale": 1.25,
        "corr0_scale": 1.10,
        "corr1_scale": 1.50,
        "corr2_scale": 1.60,
        "eps": 0.05,
        "pivotal": True,
        "single": True,
        "remove_fourmembered": True,
    },
    "dmf_options": {
        "remove_rotation_and_translation": False,
        "mass_weighted": False,
        "parallel": False,
        "eps_vel": 0.01,
        "eps_rot": 0.01,
        "beta": 10.0,
        "update_teval": False,
    },
    "k_fix": 300.0,
}

# -----------------------------------------------
# GrowingString (path representation) defaults
# -----------------------------------------------

GS_KW: Dict[str, Any] = {
    "fix_first": True,
    "fix_last": True,
    "max_nodes": 20,
    "perp_thresh": 5e-3,
    "reparam_check": "rms",
    "reparam_every": 1,
    "reparam_every_full": 1,
    "param": "equi",
    "max_micro_cycles": 10,
    "reset_dlc": True,
    "climb": True,
    "climb_rms": 5e-4,
    "climb_lanczos": True,
    "climb_lanczos_rms": 5e-4,
    "climb_fixed": False,
    "scheduler": None,
}

# -----------------------------------------------
# StringOptimizer (optimization control) defaults
# -----------------------------------------------

STOPT_KW: Dict[str, Any] = {
    "type": "string",
    "thresh": "gau_loose",
    "stop_in_when_full": 300,
    "align": False,
    "scale_step": "global",
    "max_cycles": 300,
    "dump": False,
    "dump_restart": False,
    "reparam_thresh": 0.0,
    "coord_diff_thresh": 0.0,
    "out_dir": OUT_DIR_PATH_OPT,
    "print_every": 10,
}

# -----------------------------------------------
# Path search control defaults
# -----------------------------------------------

SEARCH_KW: Dict[str, Any] = {
    "max_depth": 10,
    "stitch_rmsd_thresh": 1.0e-4,
    "bridge_rmsd_thresh": 1.0e-4,
    "max_nodes_segment": 10,
    "max_nodes_bridge": 5,
    "kink_max_nodes": 3,
    "max_seq_kink": 2,
    "refine_mode": None,
}

# -----------------------------------------------
# IRC defaults
# -----------------------------------------------

IRC_KW: Dict[str, Any] = {
    "step_length": 0.10,
    "max_cycles": 125,
    "downhill": False,
    "forward": True,
    "backward": True,
    "root": 0,
    "hessian_init": "calc",
    "displ": "energy",
    "displ_energy": 1.0e-3,
    "displ_length": 0.10,
    "rms_grad_thresh": 1.0e-3,
    "hard_rms_grad_thresh": None,
    "energy_thresh": 1.0e-6,
    "imag_below": 0.0,
    "force_inflection": True,
    "check_bonds": False,
    "out_dir": OUT_DIR_IRC,
    "prefix": "",
    "hessian_update": "bofill",
    "hessian_recalc": None,
    "max_pred_steps": 500,
    "loose_cycles": 3,
    "corr_func": "mbs",
}

# -----------------------------------------------
# Microiteration defaults (opt heavy / tsopt heavy)
# -----------------------------------------------

MICROITER_KW: Dict[str, Any] = {
    "micro_thresh": None,           # Convergence threshold for MM relaxation (None → same as macro thresh)
    "micro_max_cycles": 10000,      # Max LBFGS cycles per micro iteration
}

# -----------------------------------------------
# Frequency analysis defaults
# -----------------------------------------------

FREQ_KW: Dict[str, Any] = {
    "amplitude_ang": 0.8,
    "n_frames": 20,
    "max_write": 10,
    "sort": "value",
    "out_dir": OUT_DIR_FREQ,
}

# -----------------------------------------------
# Thermochemistry defaults
# -----------------------------------------------

THERMO_KW: Dict[str, Any] = {
    "temperature": 298.15,
    "pressure_atm": 1.0,
    "dump": False,
}

# -----------------------------------------------
# TS optimization mode aliases
# -----------------------------------------------

TSOPT_MODE_ALIASES = (
    (("grad", "light", "dimer"), "dimer"),
    (("hess", "heavy", "rsirfo"), "rsirfo"),
)

# -----------------------------------------------
# Dimer defaults for TS optimization
# -----------------------------------------------

DIMER_KW: Dict[str, Any] = {
    "length": 0.0189,
    "rotation_max_cycles": 15,
    "rotation_method": "fourier",
    "rotation_thresh": 1e-4,
    "rotation_tol": 1,
    "rotation_max_element": 0.001,
    "rotation_interpolate": True,
    "rotation_disable": False,
    "rotation_disable_pos_curv": True,
    "rotation_remove_trans": True,
    "trans_force_f_perp": True,
    "bonds": None,
    "N_hessian": None,
    "bias_rotation": False,
    "bias_translation": False,
    "bias_gaussian_dot": 0.1,
    "seed": None,
    "write_orientations": True,
    "forward_hessian": True,
}

# -----------------------------------------------
# Hessian-dimer defaults for TS optimization
# -----------------------------------------------

HESSIAN_DIMER_KW: Dict[str, Any] = {
    "thresh_loose": "gau_loose",
    "thresh": "baker",
    "update_interval_hessian": 500,
    "neg_freq_thresh_cm": 5.0,
    "flatten_amp_ang": 0.10,
    "flatten_max_iter": 50,
    "flatten_sep_cutoff": 0.0,
    "flatten_k": 10,
    "flatten_loop_bofill": False,
    "mem": 100000,
    "device": "auto",
    "root": 0,
    # --- ONIOM-specific keys (not present in pdb2reaction) ---
    "partial_hessian_flatten": True,
    "ml_only_hessian_dimer": False,
}

# -----------------------------------------------
# RS-I-RFO defaults for TS optimization (heavy mode)
# -----------------------------------------------

# Inherit shared keys from RFO_KW, but exclude RFOptimizer-only params
# that RSIRFOptimizer (via TSHessianOptimizer → Optimizer) does not accept.
_RFO_ONLY_KEYS = {
    "gediis", "gdiis", "gdiis_thresh", "gediis_thresh",
    "gdiis_test_direction", "adapt_step_func", "rfo_overlaps",
}

RSIRFO_KW: Dict[str, Any] = {
    **{k: v for k, v in RFO_KW.items() if k not in _RFO_ONLY_KEYS},
    "thresh": "baker",
    "trust_radius": 0.10,
    "trust_max": 0.10,
    "max_energy_incr": None,
    "hessian_update": "bofill",
    "hessian_init": "calc",
    "hessian_recalc": 500,
    "small_eigval_thresh": 1e-8,
    "assert_neg_eigval": False,
    "track_mode_by_overlap": False,
    "out_dir": OUT_DIR_TSOPT,
}

# -----------------------------------------------
# DFT single-point defaults
# -----------------------------------------------

DFT_KW: Dict[str, Any] = {
    "func_basis": "wb97m-v/def2-tzvpd",
    "max_cycle": 100,
    "conv_tol": 1e-9,
    "grid_level": 3,
    "verbose": 4,
    "out_dir": "./result_dft/",
    "lowmem": True,  # Use gpu4pyscf rks_lowmem for closed-shell GPU runs (auto-fallback otherwise)
}

# Note: normalize_choice and deep_update are now in utils.py to avoid duplication
