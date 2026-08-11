# `mlmm tsopt`

## Purpose

Transition-state optimization. Two algorithms: the default full-Hessian
RS-I-RFO (`--opt-mode hess`/`rsirfo`), and the lighter Hessian-Guided
Dimer (`--opt-mode grad`/`dimer`). Use after
`path-search` or `scan` to refine
a HEI to a true first-order saddle, or as a standalone validator on an
externally-generated TS guess.

## Synopsis

```bash
mlmm tsopt -i ts_guess.{pdb,xyz} --parm real.parm7 \
    [-q 0 -m 1] [-l 'RES:Q,...'] \
    [--opt-mode grad|hess|light|heavy|dimer|rsirfo] \
    [--max-cycles 10000] \
    [-b uma|orb|mace|aimnet2] [-o ./result_tsopt/]
```


## ML/MM-aware flags (mlmm-toolkit specific)

In addition to the common flags below,
**`mlmm-toolkit` requires an Amber topology** and supports layer-aware
selection. Most subcommands accept:

| flag | purpose |
|---|---|
| `--parm FILE` | Amber `parm7` topology of the whole enzyme — **required** |
| `--model-pdb FILE` | Explicit ML-region PDB; takes precedence over B-factor ML membership |
| `--detect-layer` | Automatically read B-factor layers; explicit ML membership retains valid movable/frozen MM layers. Enabled by default. |
| `--model-indices` | Explicit ML atom indices used when `--model-pdb` is omitted; takes precedence over B-factor ML membership |
| `--ref-pdb FILE` | Full-enzyme PDB used as topology reference for XYZ inputs |
| `--link-atom-method [scaled\|fixed]` | g-factor (default) or fixed 1.09/1.01 Å |
| `--embedcharge / --no-embedcharge` | Unavailable in v0.3.3; use `--no-embedcharge` |
| `-q, --charge` | **ML-region** charge (not whole-system) |
| `-l, --ligand-charge` | Per-residue charge mapping for ML region |

Inspect via `mlmm <subcommand> --help` and `mlmm <subcommand> --help-advanced`.

## Key flags

| flag | type | default | description |
|---|---|---|---|
| `-i, --input` | path | required | TS candidate; `.pdb` / `.xyz` (XYZ requires `--ref-pdb`) |
| `-q` / `-l` / `-m` | — | — | Charge / spin (common conventions) |
| `--opt-mode` | str | `hess` | `grad`/`dimer` (Hessian-Guided Dimer) or `hess`/`rsirfo` (RS-I-RFO); also `trim` (TRIM/Helgaker) and `rsprfo` (RS-P-RFO/Banerjee); the mlmm-only `light` / `heavy` shortcuts are also accepted (light = Dimer, heavy = full-Hessian RS-I-RFO) |
| `--max-cycles` | int | 10000 | Optimization step cap |
| `--hessian-calc-mode` | str | `FiniteDifference` | `Analytical` or `FiniteDifference`; check `RSIRFO_KW` / `DIMER_KW` |
| `--ref-mode` | path | none | Advanced Cartesian 3N MEP tangent for initial-root selection and overlap tracking. `all` supplies it by default; with `all --no-tsopt-from-mep-tan`, TSOPT selects from the initial-structure Hessian modes. Ordinary standalone runs omit it. |
| `--precision` | str | backend-specific | UMA/AIMNet2 fp32; ORB/MACE fp64; AIMNet2 rejects fp64 |
| `--workers` | int | 1 | UMA predictor workers; `>1` requires `fairchem-core[extras]` and is incompatible with `Analytical` |
| `--allow-charge-mult-mismatch` | flag | off | Warn and skip ML-region charge/multiplicity electron-parity validation for an intentional mismatch |
| `-b, --backend` | str | `uma` | MLIP backend |
| `-o, --out-dir` | path | `./result_tsopt/` | Output directory |
| `--config` / `--show-config` / `--dry-run` / `--help-advanced` | — | — | Standard |

`tsopt` always forces `reject_uphill=False`, regardless of optimizer mode or
YAML. Uphill trial steps can be part of saddle-point mode following. The
`--reject-uphill/--no-reject-uphill` toggle belongs only to minimum
optimization (`opt`) and post-IRC endpoint refinement (`all`).

## Examples

### Default RS-I-RFO

```bash
mlmm tsopt -i hei.xyz --parm real.parm7 --ref-pdb enzyme_layered.pdb \
    -q 0 -m 1 -b uma -o result_tsopt
```

### Dimer mode

```bash
mlmm tsopt -i hei.xyz --parm real.parm7 --ref-pdb enzyme_layered.pdb -q 0 -m 1 \
    --opt-mode dimer -b uma -o result_tsopt_dimer
```

### Tighter convergence on an ill-conditioned saddle

```bash
mlmm tsopt -i hei.xyz --parm real.parm7 --ref-pdb enzyme_layered.pdb \
    -l 'SAM:1,GPP:-3' \
    --opt-mode rsirfo --max-cycles 200 -b mace \
    -o result_tsopt_rsirfo
```

## Output

```
result_tsopt/
├── result.json                     # when --out-json
├── final_geometry.{xyz,pdb}        # final geometry; check result.json status
├── optimization_trj.xyz            # macro-cycle trajectory
├── optimization_all_trj.xyz        # full per-step trajectory (when --dump)
└── vib/                            # imaginary-mode vibrations
    └── imag_*.{pdb,xyz}            # mode displacement visualization
```

`result.json` keys:

```python
import json
d = json.load(open("result_tsopt/result.json"))
print(d["status"])                      # "converged" / "stalled" / "not_converged" / "unverified"
print(d["energy_hartree"])
print(d["n_imaginary_modes"])           # should be 1 for a real TS
print(d["imaginary_frequencies_cm"])    # list of cm⁻¹
print(d["files"]["final_geometry_xyz"]) # final_geometry.xyz
print(d["rigid_projection"]["treatment"], d["rigid_projection"]["effective_rank"])
```

## `--opt-mode` choice

| Mode | Algorithm | When |
|---|---|---|
| `hess` / `rsirfo` (default) | RS-I-RFO with full Hessian | Direct curvature treatment; memory and runtime depend on active DOFs and backend |
| `grad` / `dimer` | Hessian-Guided Dimer | Uses initial and periodic orientation Hessians, which is more robust than a random initial direction for large systems; convergence remains seed- and system-dependent |

If Dimer stalls, inspect the followed mode and step diagnostics, then compare
RS-I-RFO on the same seed rather than using a universal cycle threshold.

## Validation: imaginary modes

A real TS has exactly one imaginary frequency that corresponds to the
reaction coordinate.

```python
import json
d = json.load(open("result_tsopt/result.json"))
if d["status"] != "converged":
    print("NOT CONVERGED:", d["status"])
elif d["n_imaginary_modes"] == 1:
    print("OK: single imaginary mode at", d["imaginary_frequencies_cm"][0], "cm-1")
elif d["n_imaginary_modes"] == 0:
    print("BAD: collapsed to a minimum during refinement")
elif d["n_imaginary_modes"] is not None and d["n_imaginary_modes"] > 1:
    print("AMBIGUOUS: multiple imaginary modes; inspect vib/imag_*.pdb")
```

For multi-imaginary cases, visualize the modes (`pymol vib/imag_*.pdb`)
to decide whether the extra modes are spurious (translation/rotation
of frozen residues) or real chemical second-order saddle points.

`n_imaginary_modes == 0` is a failed TS optimization, even when the force
optimizer stopped normally. `--flatten` can remove surplus negative modes but
cannot create a missing reaction direction. Improve the MEP/starting guess;
`all --refine-path` is opt-in because recursive refinement can split a poor
path into several costly segments.

`--ref-mode` is an advanced `all`-workflow handoff, not a routine standalone
requirement. Supply it manually only when the non-zero 3N vector uses exactly
the same atom ordering as the TS input.

Do not confuse `--ref-mode` with the fixed constrained rigid-mode treatment.
`--ref-mode` supplies an MEP tangent; the constrained treatment removes only
full-system rigid motions that leave
frozen anchors fixed (generic rank 6/3/1/0 for 0/1/2/3+ non-collinear
anchors; realistic boundaries normally rank 0). All-frozen input is an
explicit error. A stale non-constrained YAML value fails explicitly.
`result.json.rigid_projection` records treatment, rank, Hessian source, and
shape.

## Caveats

- A converged `tsopt` is **not** a complete validation; always follow
  with `irc.md` to confirm the TS connects the expected R and P.
- `--max-cycles` is a safety cap, not evidence of correctness. On repeated
  nonconvergence, inspect the TS seed, followed mode, optimizer diagnostics,
  and backend/model behavior.
- Backend/model choice changes the curvature surface. Validate every candidate
  by exactly one imaginary mode, its displacement, and the intended IRC
  connectivity.

## See also

- `path-search.md` — produces TS candidates for `tsopt`.
- `irc.md`, `freq.md` — downstream validation.
- `mlmm-install-backends/uma.md` / `mace.md` — TS-accurate
  backends.
- Defaults: `import mlmm.core.defaults as d; print(d.RSIRFO_KW, d.DIMER_KW, d.HESSIAN_DIMER_KW)`
