# mlmm_toolkit/defaults.py

"""
Central configuration defaults for mlmm workflows.

All default dictionaries are defined here to avoid redundant definitions across modules.
Modules should import defaults from here instead of defining local copies.
"""

from typing import Any, Dict

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
# Geometry defaults
# -----------------------------------------------

GEOM_KW_DEFAULT: Dict[str, Any] = {
    "coord_type": "cart",
    "freeze_atoms": [],
}

# -----------------------------------------------
# Calculator defaults (ML/MM with UMA + hessian_ff MM)
# -----------------------------------------------

MLMM_CALC_KW: Dict[str, Any] = {
    "input_pdb": None,
    "real_parm7": None,
    "model_pdb": None,
    "model_charge": 0,
    "model_mult": 1,
    "link_mlmm": None,
    "uma_model": "uma-s-1p1",
    "uma_task_name": "omol",
    "ml_hessian_mode": "Analytical",
    "hessian_calc_mode": None,  # Alias for ml_hessian_mode
    "out_hess_torch": True,
    "H_double": False,
    "ml_device": "auto",
    "ml_cuda_idx": 0,
    "mm_backend": "hessian_ff",  # "hessian_ff" (analytical) | "openmm" (FD Hessian)
    "mm_device": "cpu",
    "mm_cuda_idx": 0,
    "mm_threads": 16,
    "mm_fd": True,
    "mm_fd_dir": None,
    "mm_fd_delta": 1e-3,         # Legacy parameter (retained)
    "symmetrize_hessian": True,  # Symmetrize final Hessian as 0.5*(H+H^T)
    "print_timing": True,        # Print ML/MM Hessian timing breakdown
    "print_vram": True,          # Print CUDA VRAM usage (peak) during Hessian
    "return_partial_hessian": False,
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
    # Explicit YAML-based layer specification (0-based indices, takes precedence over cutoffs/B-factors)
    "hess_mm_atoms": None,    # Explicit Hessian-target MM atom indices
    "movable_mm_atoms": None, # Explicit movable MM atom indices
    "frozen_mm_atoms": None,  # Explicit frozen MM atom indices
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
    "trust_min": 0.00,
    "trust_max": 0.10,
    "max_energy_incr": None,
    "hessian_update": "bfgs",
    "hessian_init": "calc",
    "hessian_recalc": 200,
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
    "k": 100,  # eV/Å²
}

# -----------------------------------------------
# Bond-change detection
# -----------------------------------------------

BOND_KW: Dict[str, Any] = {
    "device": "cuda",
    "bond_factor": 1.20,
    "margin_fraction": 0.05,
    "delta_fraction": 0.05,
}

# -----------------------------------------------
# Optimizer mode aliases
# -----------------------------------------------

OPT_MODE_ALIASES = (
    (("light", "lbfgs"), "lbfgs"),
    (("heavy", "rfo"), "rfo"),
)

TSOPT_MODE_ALIASES = (
    (("light", "lbfgs"), "light"),
    (("heavy", "rfo"), "heavy"),
)

# -----------------------------------------------
# GrowingString (path representation) defaults
# -----------------------------------------------

GS_KW: Dict[str, Any] = {
    "fix_first": True,
    "fix_last": True,
    "max_nodes": 10,
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
# Frequency analysis defaults
# -----------------------------------------------

# -----------------------------------------------
# Microiteration defaults (opt heavy / tsopt heavy)
# -----------------------------------------------

MICROITER_KW: Dict[str, Any] = {
    "micro_thresh": "gau_loose",    # Convergence threshold for MM relaxation
    "micro_max_cycles": 500,        # Max LBFGS cycles per micro iteration
}

# -----------------------------------------------
# Frequency analysis defaults
# -----------------------------------------------

FREQ_KW: Dict[str, Any] = {
    "amplitude_ang": 0.8,
    "n_frames": 20,
    "max_write": 20,
    "sort": "value",
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
    "flatten_amp_ang": 0.20,
    "flatten_max_iter": 50,
    "flatten_sep_cutoff": 0.0,
    "flatten_k": 10,
    "flatten_loop_bofill": False,
    "mem": 100000,
    "device": "auto",
    "root": 0,
    "partial_hessian_flatten": True,  # New: use partial Hessian for imaginary mode detection
    "ml_only_hessian_dimer": False,  # Use ML-region-only Hessian (no MM Hessian) for dimer orientation
}

# -----------------------------------------------
# RS-I-RFO defaults for TS optimization (heavy mode)
# -----------------------------------------------

RSIRFO_KW: Dict[str, Any] = {
    **OPT_BASE_KW,
    "thresh": "baker",
    "trust_radius": 0.10,
    "trust_update": True,
    "trust_min": 0.00,
    "trust_max": 0.30,
    "max_energy_incr": None,
    "hessian_update": "bofill",
    "hessian_init": "calc",
    "hessian_recalc": 200,
    "small_eigval_thresh": 1e-8,
    "out_dir": OUT_DIR_TSOPT,
}

# -----------------------------------------------
# DFT single-point defaults
# -----------------------------------------------

DFT_KW: Dict[str, Any] = {
    "func_basis": "wb97m-v/6-31g**",
    "max_cycle": 100,
    "conv_tol": 1e-9,
    "grid_level": 3,
}

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

# Note: normalize_choice and deep_update are now in utils.py to avoid duplication
