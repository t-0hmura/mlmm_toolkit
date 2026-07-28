# YAML Reference

## Overview


| Section | Description | Used by |
|---------|-------------|---------|
| [`geom`](#geom) | Geometry and coordinate settings | all, opt, scan, scan2d, scan3d, tsopt, freq, irc, path-opt, path-search |
| [`calc`](#calc) | ML/MM calculator settings (alias: `mlmm:`) | all, opt, scan, scan2d, scan3d, tsopt, freq, irc, path-opt, path-search |
| [`opt`](#opt) | Shared optimizer settings | opt, scan, scan2d, scan3d, tsopt, path-opt, path-search |
| [`lbfgs`](#lbfgs) | L-BFGS optimizer settings | opt, scan, scan2d, scan3d, path-opt, path-search |
| [`rfo`](#rfo) | RFO optimizer settings | opt |
| [`gs`](#gs) | Growing String Method settings | path-opt, path-search |
| [`dmf`](#dmf) | Direct Max Flux settings | path-opt, path-search |
| [`irc`](#irc-section) | IRC integration settings | irc |
| [`freq`](#freq-section) | Vibrational analysis settings | freq |
| [`thermo`](#thermo) | Thermochemistry settings | freq |
| [`dft`](#dft-section) | DFT calculation settings | dft |
| [`bias`](#bias) | Harmonic bias settings | scan, scan2d, scan3d |
| [`bond`](#bond) | Bond-change detection settings | scan, path-search |
| [`search`](#search) | Recursive path search settings | path-search |
| [`hessian_dimer`](#hessian_dimer) | Hessian Dimer TS optimization | tsopt |
| [`rsirfo`](#rsirfo) | RS-I-RFO TS optimization | tsopt |
| [`stopt`](#stopt) | String optimizer settings | path-opt, path-search |
| [`microiter`](#microiter) | Micro-iteration (MM relaxation) settings | opt, tsopt |

---

## Shared Sections

### `geom`

Geometry loading and coordinate handling.

```yaml
geom:
 coord_type: cart # Coordinate type: "cart" (Cartesian) or "dlc" (delocalized internals)
 freeze_atoms: [] # 1-based frozen-atom indices
 tr_projection: constrained # constrained (default) | legacy-active
```

**Notes:**
- Frozen atoms have zeroed forces; their Hessian columns are also zeroed
- `tr_projection: constrained` removes only full-system rigid motions that
  leave all frozen anchors fixed. Its generic effective rank is 6/3/1/0 for
  zero/one/two/at least three non-collinear anchors; realistic ML/MM boundaries
  therefore usually have rank 0. An all-frozen selection is an explicit error.
- `legacy-active` is deprecated, comparison-only, and must not be used for
  pass/HOSP transition-state certification. It uses the current common
  projection kernel and numerical rank handling; bitwise identity is not
  guaranteed for rank-degenerate geometries.
- `tr_projection` controls frozen-boundary PHVA used by `freq`, `irc`, `tsopt`,
  and `opt --flatten`. It is unrelated to the `tsopt --ref-mode` MEP tangent.
- For `irc`, `geom.coord_type` is forced to `cart` after YAML/CLI merging

---

(calc)=
### `calc` (section)

`mlmm:` is an accepted alias for this section — the loaders read `calc` first and fall back to
`mlmm`, so a file may use either name (not both).

```yaml
calc:
 # --- Input files ---
 input_pdb: null # Input PDB file path (usually set by CLI)
 real_parm7: null # Amber parm7 topology for the full (real) system
 model_pdb: null # PDB defining the ML (model) region atoms
 model_charge: 0 # Charge of the ML (model) region
 model_mult: 1 # Spin multiplicity of the ML (model) region
 link_mlmm: null # null: derive boundary pairs from parm7 bonds; list: explicit override
 link_atom_method: scaled    # Link atom placement: "scaled" (g-factor) or "fixed" (1.09/1.01 Å)

 # --- MLIP backend selection ---
 backend: uma # MLIP backend: "uma", "orb", "mace", or "aimnet2"

 # --- UMA backend settings ---
 uma_model: uma-s-1p2 # uma-s-1p2 | uma-s-1p1 | uma-m-1p1
 uma_task_name: omol # Task tag recorded in UMA batches (UMA backend only)
 uma_precision: fp32 # fp32 | fp64 (UMA backend numerical precision)

 # --- ORB backend settings ---
 orb_model: orb_v3_conservative_omol  # ORB model name (ORB backend only)
 orb_precision: float64  # ORB floating-point precision, default (ORB backend only; "float32-high" = TF32 matmul, also reachable via --precision fp32; "float32" accepted as a legacy alias)

 # --- MACE backend settings ---
 mace_model: MACE-OMOL-0 # MACE model name (MACE backend only)
 mace_dtype: float64      # MACE floating-point precision (MACE backend only)

 # --- AIMNet2 backend settings ---
 aimnet2_model: aimnet2   # AIMNet2 model name (AIMNet2 backend only)

 # --- ML device & Hessian ---
 ml_device: auto # Device for ML inference: "cuda", "cpu", or "auto"
 ml_cuda_idx: 0 # CUDA device index for ML inference
 hessian_calc_mode: FiniteDifference # ML Hessian mode: "Analytical" or "FiniteDifference"

 # --- Retired electronic-embedding compatibility keys ---
 embedcharge: false # Must remain false in v0.3.3; true is rejected
 embedcharge_cutoff: 12.0 # Inactive compatibility value
 embedcharge_step: 0.001 # Inactive compatibility value
 xtb_cmd: xtb # Inactive compatibility value
 xtb_acc: 0.2 # Inactive compatibility value
 xtb_workdir: tmp # Inactive compatibility value
 xtb_keep_files: false # Inactive compatibility value
 xtb_ncores: 4 # Inactive compatibility value

 # --- MM backend settings ---
 mm_backend: hessian_ff # MM backend; Hessian method is selected by mm_fd below
 use_cmap: true         # Preserve parm7 CMAP in both REAL and MODEL MM layers
 mm_device: cpu # Device for MM calculation: "cuda" or "cpu" (hessian_ff is CPU-only)
 mm_cuda_idx: 0 # CUDA device index for MM calculation (OpenMM only)
 mm_threads: 16 # Number of threads for MM calculation
 mm_fd: true # Use finite-difference for MM Hessian
 mm_hessian_mode: null # Explicit finite_difference/analytical mode; null maps mm_fd
 mm_fd_dir: null # Directory for MM finite-difference scratch files
 mm_fd_delta: 0.001 # Finite-difference step size for MM Hessian

 # --- Hessian output settings ---
 out_hess_torch: true # Return Hessian as torch.Tensor
 H_double: true # Assemble/return Hessian in float64
 symmetrize_hessian: true # Symmetrize final Hessian as 0.5*(H+H^T)
 return_partial_hessian: true # Active-block partial Hessian (CLI wrappers default to true)

 # --- Layer configuration ---
 freeze_atoms: [] # 1-based indices of atoms to freeze (Frozen layer)
 hess_cutoff: null # Å; null = all movable MM included in Hessian (default); >0.0 = only MM within this distance of ML
 movable_cutoff: null # Å; MM atoms within this distance of ML are movable (null = use freeze_atoms)
 use_bfactor_layers: true # If true, read layer assignments from input PDB B-factors
 hess_mm_atoms: null # Explicit Hessian-target MM atom indices (1-based; overrides cutoffs)
 movable_mm_atoms: null # Explicit movable MM atom indices (1-based; overrides cutoffs)
 frozen_mm_atoms: null # Explicit frozen MM atom indices (1-based; overrides cutoffs)

 # --- Diagnostics ---
 print_timing: true # Print ML/MM Hessian timing breakdown
 print_vram: true # Print CUDA VRAM usage during Hessian
```

**Notes:**
- The section name `calc:` is the canonical form; `mlmm:` is accepted as a legacy alias (recognized by `opt`, `tsopt`, `freq`, `irc`, `dft`, `path-opt`, `path-search`, `scan`, `scan2d`, `scan3d`). When both are present, `calc:` takes precedence.
- `backend` selects the MLIP backend: `uma` (default), `orb`, `mace`, or `aimnet2`. Alternative backends require optional dependencies (`pip install "mlmm-toolkit[orb]"`, etc.)
- Backend-specific model keys are only relevant when the corresponding backend is selected:
  - `uma_model`, `uma_task_name` — UMA backend only
  - `orb_model`, `orb_precision` — ORB backend only
  - `mace_model`, `mace_dtype` — MACE backend only
  - `aimnet2_model` — AIMNet2 backend only
- `embedcharge` must remain `false` in v0.3.3. `embedcharge: true`, `--embedcharge`, or an explicit CLI `--embedcharge-cutoff` is rejected before calculation because the retired path double-counted ML--MM electrostatics and used an inconsistent uncapped boundary model.
- YAML `embedcharge_cutoff` and the other `embedcharge_*`/`xtb_*` keys remain readable only as inert configuration-compatibility values; they do not activate a supported calculation path.
- `hessian_calc_mode: Analytical` requests the backend's analytical Hessian and is recommended when sufficient VRAM is available. UMA, ORB, MACE, and AIMNet2 implement this path; an incompatible installed backend version raises an error instead of silently changing methods. `workers > 1` cannot be combined with an analytical Hessian and also raises an error.
- `mm_fd: true` uses a finite-difference MM Hessian; `false` uses the
  analytical `hessian_ff` Hessian. `mm_hessian_mode` is the explicit
  `finite_difference`/`analytical` spelling; when it is `null`, `mm_fd`
  supplies the backward-compatible selection.
- `use_cmap: true` (default) preserves CMAP in both REAL and MODEL MM layers when the parm7 contains it. Set `false` only for an explicit modified-force-field calculation; the opt-out removes CMAP from both layers.
- `real_parm7` and `model_pdb` are required for ML/MM calculations
- `model_charge` and `model_mult` override `-q` and `-m` for the ML region specifically
- `opt`, `tsopt`, `irc`, and `freq` use partial Hessian by default when `calc.return_partial_hessian` is not explicitly set in YAML.
- To force full Hessian output in those commands, set `calc.return_partial_hessian: false` explicitly.
- `irc` forces `geom.coord_type = cart` regardless of YAML.

### `opt`

Shared optimizer controls used by both L-BFGS and RFO.

```yaml
opt:
 type: string # StringOptimizer-only (path-opt/path-search): optimizer type label
 thresh: gau # Convergence preset: gau_loose, gau, gau_tight, gau_vtight, baker, never
 stop_in_when_full: 300 # StringOptimizer-only: early stop threshold when string is full
 align: false # StringOptimizer-only: alignment toggle
 scale_step: global # StringOptimizer-only: step scaling mode
 max_cycles: 10000 # Maximum optimizer iterations
 print_every: 100 # Logging stride
 min_step_norm: 1.0e-08 # Minimum step norm for acceptance
 assert_min_step: true # Stop if steps fall below threshold
 rms_force: null # Explicit RMS force target
 rms_force_only: false # Rely only on RMS force convergence
 max_force_only: false # Rely only on max force convergence
 force_only: false # Skip displacement checks
 converge_to_geom_rms_thresh: 0.05 # RMS threshold when converging to reference geometry
 overachieve_factor: 0.0 # Factor to tighten thresholds
 check_eigval_structure: false # Validate Hessian eigenstructure
 energy_plateau: true # Fallback: stop the optimizer as stalled (status: "stalled", never converged) when the energy plateaus
 energy_plateau_thresh: 1.0e-4 # Energy range tolerance in au (~0.06 kcal/mol)
 energy_plateau_window: 50 # Number of trailing steps used for plateau detection
 line_search: true # Enable line search
 dump: false # Dump trajectory/restart data
 dump_restart: false # Dump restart checkpoints
 reparam_thresh: 0.0 # StringOptimizer-only: reparameterization threshold
 coord_diff_thresh: 0.0 # StringOptimizer-only: coordinate difference threshold
 prefix: "" # Filename prefix
 out_dir: ./result_opt/ # Output directory
```

**Convergence Presets:**

| Preset | Max Force | RMS Force | Max Step | RMS Step |
|--------|-----------|-----------|----------|----------|
| `gau_loose` | 2.5e-3 | 1.7e-3 | 1.0e-2 | 6.7e-3 |
| `gau` | 4.5e-4 | 3.0e-4 | 1.8e-3 | 1.2e-3 |
| `gau_tight` | 1.5e-5 | 1.0e-5 | 6.0e-5 | 4.0e-5 |
| `gau_vtight` | 2.0e-6 | 1.0e-6 | 6.0e-6 | 4.0e-6 |
| `baker` | 3.0e-4 | 2.0e-4 | 3.0e-4 | 2.0e-4 |

`baker` is the exception to the four-column rule: it converges when
`max(|force|) <= 3e-4` **and** (`|delta E| < 1e-6` **or**
`max(|step|) <= 3e-4`). Its RMS force and RMS step columns are diagnostics,
not additional terminal gates.

**Energy plateau fallback:**

When `energy_plateau: true` (default), a plateau stops the optimizer with
`status: "stalled"` (a distinct non-converged outcome, never `converged`) if the
energy range `max(E) - min(E)` over the last `energy_plateau_window` steps (default
50) falls below `energy_plateau_thresh` (default `1.0e-4` au, ~0.06 kcal/mol).

This is a safety net for ML/MM optimizations where the MLIP force noise floor
can exceed the standard gradient-based convergence thresholds, causing the
optimizer to wander indefinitely. Once the energy itself has plateaued within
the MLIP's numerical precision, further cycles cannot reduce the residual force
below the noise floor, so the plateau trigger stops the run cleanly.

Set `energy_plateau: false` to disable the fallback (the optimizer will then
rely solely on the `thresh` preset). The plateau check is automatically skipped
for chain-of-states (COS) optimizers (e.g., GS, DMF string optimizers).

---

### `lbfgs`

L-BFGS optimizer settings (extends `opt`).

```yaml
lbfgs:
 # Inherits all opt settings, plus:
 keep_last: 7 # History size for L-BFGS buffers
 beta: 1.0 # Initial damping beta
 gamma_mult: false # Multiplicative gamma update toggle
 max_step: 0.3 # Maximum step length
 control_step: true # Control step length adaptively
 double_damp: true # Double damping safeguard
 mu_reg: null # Regularization strength
 max_mu_reg_adaptions: 10 # Cap on mu adaptations
```

---

### `rfo`

Rational Function Optimizer settings (extends `opt`).

```yaml
rfo:
 # Inherits all opt settings, plus:
 trust_radius: 0.10 # Trust-region radius
 trust_update: true # Enable trust-region updates
 trust_min: 0.0001 # Minimum trust radius
 trust_max: 0.10 # Maximum trust radius (tuned for ML/MM stability)
 max_energy_incr: null # Allowed energy increase per step
 hessian_update: bfgs # Hessian update scheme: bfgs, bofill, etc.
 hessian_init: calc # Hessian initialization: calc, unit, etc.
 hessian_recalc: 500 # Rebuild Hessian every N steps
 hessian_recalc_adapt: null # Adaptive Hessian rebuild factor
 small_eigval_thresh: 1.0e-08 # Eigenvalue threshold for stability
 alpha0: 1.0 # Initial micro step
 max_micro_cycles: 50 # Micro-iteration limit
 rfo_overlaps: false # Enable RFO overlaps
 gediis: false # Enable GEDIIS
 gdiis: true # Enable GDIIS
 gdiis_thresh: 0.0025 # GDIIS acceptance threshold
 gediis_thresh: 0.01 # GEDIIS acceptance threshold
 gdiis_test_direction: true # Test descent direction before DIIS
 adapt_step_func: true # Adaptive step scaling
```

---

## Path Optimization Sections

### `gs`

Growing String Method settings.

```yaml
gs:
 fix_first: true # Keep first endpoint fixed
 fix_last: true # Keep last endpoint fixed
 max_nodes: 20 # Maximum string nodes (internal images)
 perp_thresh: 0.005 # Perpendicular displacement threshold
 reparam_check: rms # Reparameterization check metric
 reparam_every: 1 # Reparameterization stride
 reparam_every_full: 1 # Full reparameterization stride
 param: equi # Parametrization scheme
 max_micro_cycles: 10 # Micro-iteration limit
 reset_dlc: true # Rebuild delocalized coordinates each step
 climb: true # Enable climbing image
 climb_rms: 0.0005 # Climbing RMS threshold
 climb_lanczos: true # Lanczos refinement for climbing
 climb_lanczos_rms: 0.0005 # Lanczos RMS threshold
 climb_fixed: false # Keep climbing image fixed
 scheduler: null # Optional scheduler backend
```

---

### `dmf`

Direct Max Flux settings for MEP optimization.

```yaml
dmf:
 max_cycles: 300 # Maximum DMF/IPOPT iterations (overridden by --max-cycles)
 correlated: true # Correlated DMF propagation
 sequential: true # Sequential DMF execution
 fbenm_only_endpoints: false # Run FB-ENM beyond endpoints
 fbenm_options:
   delta_scale: 0.2 # FB-ENM displacement scaling
   bond_scale: 1.25 # Bond cutoff scaling
   fix_planes: true # Enforce planar constraints
 cfbenm_options:
   bond_scale: 1.25 # CFB-ENM bond cutoff scaling
   corr0_scale: 1.1 # Correlation scale for corr0
   corr1_scale: 1.5 # Correlation scale for corr1
   corr2_scale: 1.6 # Correlation scale for corr2
   eps: 0.05 # Correlation epsilon
   pivotal: true # Pivotal residue handling
   single: true # Single-atom pivots
   remove_fourmembered: true # Prune four-membered rings
 dmf_options:
   remove_rotation_and_translation: false # Keep rigid-body motions
   mass_weighted: false # Toggle mass weighting
   parallel: false # Enable parallel DMF
   eps_vel: 0.01 # Velocity tolerance
   eps_rot: 0.01 # Rotational tolerance
   beta: 10.0 # Beta parameter for DMF
   update_teval: false # Update transition evaluation
 k_fix: 300.0 # Harmonic constant for restraints
```

---

### `search`

Recursive path search settings (path-search only).

```yaml
search:
 max_depth: 10 # Recursion depth limit
 stitch_rmsd_thresh: 0.0001 # RMSD threshold for stitching segments
 bridge_rmsd_thresh: 0.0001 # RMSD threshold for bridging nodes
 max_nodes_segment: 20 # Max nodes per segment
 max_nodes_bridge: 5 # Max nodes per bridge
 kink_max_nodes: 3 # Max nodes for kink optimizations
 max_seq_kink: 2 # Max sequential kinks
 refine_mode: null # Refinement strategy: peak, minima, or null (auto)
```

---

### `stopt`

String optimizer settings for path-opt and path-search. Controls the GS/DMF string
optimization (not individual node optimization).

```yaml
stopt:
 type: string           # Optimizer type label (used by StringOptimizer)
 thresh: gau_loose      # Convergence preset for string optimization
 stop_in_when_full: 300 # Early stop threshold when string is full
 align: false           # Alignment toggle
 scale_step: global     # Step scaling mode
 max_cycles: 300        # Maximum string optimizer iterations
 dump: false            # Dump trajectory/restart data
 dump_restart: false    # Dump restart checkpoints
 reparam_thresh: 0.0    # Reparameterization threshold
 coord_diff_thresh: 0.0 # Coordinate difference threshold
 out_dir: ./result_path_opt/  # Output directory
 print_every: 10        # Logging stride
 lbfgs:
   # Same keys as lbfgs section (for single-structure optimizer)
   thresh: gau
   max_cycles: 10000
   #... (see lbfgs section)
```

**Notes:**
- `stopt.lbfgs` configures the single-structure L-BFGS optimizer used for
  HEI+/-1 endpoint optimization and kink node optimization within path-search.
  Only L-BFGS is consumed at this nested level; a `stopt.rfo:` block is not honored.
- The outer `stopt` keys control the string optimizer (GS or DMF wrapper)

---

## TS Optimization Sections

### `hessian_dimer`

Hessian Dimer TS optimization settings (`tsopt --opt-mode grad`).

```yaml
hessian_dimer:
 thresh_loose: gau_loose # Loose convergence preset
 thresh: baker # Main convergence preset
 update_interval_hessian: 500 # Hessian rebuild cadence
 neg_freq_thresh_cm: 5.0 # Negative frequency threshold (cm-1)
 flatten_amp_ang: 0.1 # Flattening amplitude (Å)
 flatten_max_iter: 50 # Flattening iteration cap (default 50; --no-flatten sets to 0)
 flatten_sep_cutoff: 0.0 # Minimum distance between representative atoms
 flatten_k: 10 # Representative atoms sampled per mode
 flatten_loop_bofill: false # Bofill update for flatten displacements
 partial_hessian_flatten: true # Use partial Hessian for imaginary mode detection
 ml_only_hessian_dimer: false # Use only ML-region atoms for Dimer rotation
 mem: 100000 # Memory limit for solver
 device: auto # Device selection for eigensolver
 root: 0 # Targeted TS root index
 dimer:
   length: 0.0189 # Dimer separation (Bohr)
   rotation_max_cycles: 15 # Max rotation iterations
   rotation_method: fourier # Rotation optimizer method
   rotation_thresh: 0.0001 # Rotation convergence threshold
   rotation_tol: 1 # Rotation tolerance factor
   rotation_max_element: 0.001 # Max rotation matrix element
   rotation_interpolate: true # Interpolate rotation steps
   rotation_disable: false # Disable rotations entirely
   rotation_disable_pos_curv: true # Disable when positive curvature detected
   rotation_remove_trans: true # Remove the selected rigid-null components
   trans_force_f_perp: true # Project forces perpendicular to translation
   bonds: null # Bond list for constraints
   N_hessian: null # Hessian size override
   bias_rotation: false # Bias rotational search
   bias_translation: false # Bias translational search
   bias_gaussian_dot: 0.1 # Gaussian bias dot product
   seed: null # RNG seed for rotations
   write_orientations: true # Write rotation orientations
   forward_hessian: true # Propagate Hessian forward
 lbfgs:
   # Same keys as lbfgs section
   thresh: baker
   max_cycles: 10000
```

**Notes:**
- `flatten_max_iter` controls the maximum number of imaginary-mode flattening iterations. The default value is 50.
- The CLI flags `--flatten` / `--no-flatten` (in `tsopt` and `all`) interact with this setting: `--flatten` enables the flattening loop with the default `flatten_max_iter` (50); `--no-flatten` forces `flatten_max_iter` to 0, effectively disabling the loop. An explicit YAML value for `flatten_max_iter` takes precedence when provided alongside `--flatten`.

---

### `rsirfo`

RS-I-RFO TS optimization settings (`tsopt --opt-mode hess`).

```yaml
rsirfo:
 thresh: baker # RS-I-RFO convergence preset
 max_cycles: 10000 # Iteration cap
 print_every: 100 # Logging stride
 min_step_norm: 1.0e-08 # Minimum accepted step norm
 assert_min_step: true # Assert when steps stagnate
 roots: [0] # Target root indices (pysisyphus default; not set by mlmm)
 hessian_ref: null # Reference Hessian
 rx_modes: null # Reaction-mode definitions
 prim_coord: null # Primary coordinates to monitor
 rx_coords: null # Reaction coordinates to monitor
 hessian_update: bofill # Hessian update scheme
 hessian_init: calc # Hessian initialization
 hessian_recalc_reset: true # Reset recalc counter after exact Hessian
 max_micro_cycles: 50 # Micro-iterations per macro cycle
 augment_bonds: false # Augment reaction path based on bond analysis
 min_line_search: false # Line search along imaginary mode (pysisyphus default)
 max_line_search: false # Line search in minimized subspace (pysisyphus default)
 assert_neg_eigval: false # Require negative eigenvalue at convergence
 track_mode_by_overlap: false # mlmm-specific: track the target mode by overlap
 trust_radius: 0.10 # Trust region radius
 trust_update: true # Trust region update
 trust_min: 0.0001 # Minimum trust radius
 trust_max: 0.10 # Maximum trust radius (tuned for ML/MM stability)
 hessian_recalc: 500 # Hessian rebuild cadence
 small_eigval_thresh: 1.0e-08 # Eigenvalue threshold for stability
 out_dir: ./result_tsopt/ # Output directory
```

---

## IRC Section

(irc-section)=
### `irc` (section)

IRC integration settings.

```yaml
irc:
 step_length: 0.1 # Integration step length
 max_cycles: 125 # Maximum steps along IRC
 downhill: false # Follow downhill direction only
 forward: true # Propagate in forward direction
 backward: true # Propagate in backward direction
 root: 0 # Normal-mode root index
 hessian_init: calc # Hessian initialization source
 hessian_update: bofill # Hessian update scheme
 hessian_recalc: null # Hessian rebuild cadence
 displ: energy # Displacement construction method
 displ_energy: 0.001 # Energy-based displacement scaling
 displ_length: 0.1 # Length-based displacement fallback
 rms_grad_thresh: 0.001 # RMS gradient convergence threshold
 hard_rms_grad_thresh: null # Hard RMS gradient stop
 energy_thresh: 0.000001 # Energy change threshold
 imag_below: 0.0 # Imaginary frequency cutoff
 force_inflection: true # Enforce inflection detection
 check_bonds: false # Check bonds during propagation
 out_dir: ./result_irc/ # Output directory
 prefix: "" # Filename prefix
 dump_fn: irc_data.h5 # IRC data filename
 dump_every: null # Disabled by default; set a positive dump stride to opt in
 max_pred_steps: 500 # Predictor-corrector max steps
 loose_cycles: 3 # Loose cycles before tightening
 corr_func: mbs # Correlation function choice
```

---

## Vibrational Analysis Sections

(freq-section)=
### `freq` (section)

Vibrational frequency analysis settings.

```yaml
freq:
 active_dof_mode: partial # Active-atom selection: "all" | "ml-only" | "partial" | "unfrozen"
 amplitude_ang: 0.8 # Displacement amplitude for modes (Å)
 n_frames: 20 # Number of frames per mode animation
 max_write: 10 # Maximum number of modes to write
 sort: value # Sort order: "value" or "abs"
 out_dir: ./result_freq/ # Output directory
```

**Notes:**
- `active_dof_mode` selects which atoms participate in the vibrational analysis. `all` uses every atom; `ml-only` restricts to ML-region atoms; `partial` (default) uses ML + Movable-MM atoms; `unfrozen` uses every non-frozen atom. The CLI flag `--active-dof-mode` overrides the YAML value when explicitly passed.

---

### `sp` (section)

Single-point settings. Read only by `mlmm sp`.

```yaml
sp:
 hess: false # Also compute the full ONIOM Hessian
 hessian_calc_mode: FiniteDifference # "FiniteDifference" | "Analytical"
 out_dir: ./result_sp/ # Output directory
```

**Notes:**
- The matching CLI flags (`--hess`, `--hessian-calc-mode`, `-o/--out-dir`) override these when passed explicitly.

---

### `thermo`

Thermochemistry settings.

```yaml
thermo:
 temperature: 298.15 # Thermochemistry temperature (K)
 pressure_atm: 1.0 # Thermochemistry pressure (atm)
 symmetry_number: 1 # External rotational symmetry number (integer >= 1)
 dump: false # Write thermoanalysis.yaml
```

---

### `microiter`

Micro-iteration settings for ML/MM optimization. When `--microiter` is enabled,
the MM region is relaxed (with frozen ML atoms) between each macro-step of the
ML-region optimizer. This can substantially reduce the number of expensive
ML Hessian evaluations needed.

```yaml
microiter:
 micro_thresh: null       # Convergence preset for MM relaxation (L-BFGS); null → same as macro thresh
 micro_max_cycles: 10000  # Maximum L-BFGS iterations per micro-iteration
```

**Notes:**
- Enabled via `--microiter` / `--no-microiter` CLI flag (default: on)
- Available in `opt` (with `--opt-mode hess`) and `tsopt` (with `--opt-mode hess`)
- Uses L-BFGS to minimize MM-region forces while ML atoms are frozen
- `micro_thresh` accepts the same presets as `opt.thresh` (gau_loose, gau, gau_tight, etc.); when `null` or omitted, defaults to the same threshold as the macro step

---

## DFT Section

(dft-section)=
### `dft` (section)

DFT calculation settings.

```yaml
dft:
 func_basis: wb97m-v/def2-tzvpd # Combined "FUNC/BASIS" string
 conv_tol: 1.0e-09 # SCF convergence tolerance (Hartree)
 max_cycle: 100 # Maximum SCF iterations
 grid_level: 3 # PySCF grid level
 engine: gpu # Compute engine: "gpu" (gpu4pyscf) or "cpu" (pyscf); CLI --engine takes precedence
 ecp: null # ECP basis name; null auto-derives from def2-* basis sets
 lowmem: true # Use gpu4pyscf rks_lowmem.RKS for closed-shell GPU runs
 verbose: 0 # PySCF verbosity level; CLI -v 2/3 raises runtime PySCF verbosity to >=4
 out_dir: ./result_dft/ # Output directory
```

**Notes:**
- `engine`: `gpu` runs through gpu4pyscf (closed-shell uses `rks_lowmem.RKS` when `lowmem: true`); `cpu` falls back to standard PySCF RKS/UKS. The CLI flag `--engine` overrides the YAML value when explicitly passed.
- `ecp`: when the basis name starts with `def2-` and `ecp` is `null`, the same basis name is used as the ECP automatically. Set explicitly to override.

---

## Scan Sections

### `bias`

Harmonic bias settings for scans.

```yaml
bias:
 k: 300.0 # Harmonic bias strength (eV/Å²)
```

---

### `bond`

MLIP-based bond-change detection.

```yaml
bond:
 device: auto # MLIP device for bond analysis
 bond_factor: 1.2 # Covalent-radius scaling for cutoff
 margin_fraction: 0.05 # Fractional tolerance for comparisons
 delta_fraction: 0.05 # Minimum relative change to flag bond formation/breaking
```

---

## Example: Complete Configuration File

Below is a full example combining multiple sections:

```yaml
# mlmm configuration example

geom:
 coord_type: cart
 freeze_atoms: []
 tr_projection: constrained

calc:
 model_charge: 0
 model_mult: 1
 backend: uma                  # MLIP backend: "uma", "orb", "mace", or "aimnet2"
 embedcharge: false            # Compatibility tombstone; true is rejected
 uma_model: uma-s-1p2          # uma-s-1p2 | uma-s-1p1 | uma-m-1p1
 ml_device: auto
 hessian_calc_mode: Analytical   # Recommended when VRAM permits
 mm_device: cpu
 mm_fd: true
 use_bfactor_layers: true # Read layers from PDB B-factors

gs:
 max_nodes: 20
 climb: true
 climb_lanczos: true

opt:
 thresh: gau
 max_cycles: 300
 dump: false
 out_dir: ./result_all/

stopt:
 thresh: gau_loose
 max_cycles: 300
 lbfgs:
   thresh: gau
   max_cycles: 10000

bond:
 bond_factor: 1.2
 delta_fraction: 0.05

search:
 max_depth: 10
 max_nodes_segment: 20

freq:
 max_write: 10
 amplitude_ang: 0.8

thermo:
 temperature: 298.15
 pressure_atm: 1.0
 symmetry_number: 1

dft:
 func_basis: wb97m-v/def2-tzvpd
 grid_level: 3
```

---

## See Also

- [all](all.md) — End-to-end workflow
- [opt](opt.md) — Single-structure optimization
- [tsopt](tsopt.md) — Transition state optimization
- [path-search](path-search.md) — Recursive MEP search
- [freq](freq.md) — Vibrational analysis
- [dft](dft.md) — DFT calculations
- [Concepts](concepts.md) — ML/MM 3-layer system and ONIOM energy decomposition
