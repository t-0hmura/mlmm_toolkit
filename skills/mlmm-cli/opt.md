# `mlmm opt`

## Purpose

Single-structure geometry optimization with L-BFGS or RFO.
Use this to relax a starting geometry to its nearest local minimum
before feeding it to `path-search` / `path-opt`, or as a post-IRC
endpoint refinement.

## Synopsis

```bash
mlmm opt -i input.pdb --parm real.parm7 [-q 0 -m 1] \
    [--opt-mode grad|hess|light|heavy|lbfgs|rfo] \
    [-b uma|orb|mace|aimnet2] [-o ./result_opt/]
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
| `-i, --input` | path | required | `.pdb` / `.xyz` (XYZ requires `--ref-pdb`) |
| `-q` / `-l` / `-m` | — | — | Charge / spin (common conventions) |
| `--opt-mode` | str | `grad` | `grad` (L-BFGS) or `hess` (RFO); aliases `lbfgs` / `rfo` and the mlmm-only `light` / `heavy` shortcuts (light = L-BFGS, heavy = RFO with full Hessian) are accepted |
| `--reject-uphill / --no-reject-uphill` | toggle | off | Opt in to rejecting an energy-raising Hessian/RFO trial above `1e-4` Hartree, restoring the lower-energy geometry and shrinking the trust radius. At the emergency floor, run one final convergence check on the retained geometry. Ignored in L-BFGS mode. |
| `--mm-only` / `--no-mm-only` | flag | `False` | Skip the MLIP component and minimize on the MM force field only. Layers honored as usual; only `--opt-mode grad` supported (microiter auto-off). Useful as a cheap MM pre-relaxation before ML/MM ONIOM opt. |
| `--max-cycles` | int | (live default) | Stop after N cycles; check `OPT_BASE_KW` |
| `-b, --backend` | str | `uma` | MLIP backend |
| `-o, --out-dir` | path | `./result_opt/` | Output directory |
| `--config` / `--show-config` / `--dry-run` / `--help-advanced` | — | — | Standard |

With `--thresh baker`, convergence requires ALL of `max(|force|) <= 3e-4`,
`rms(force) <= 2e-4`, `max(|step|) <= 3e-4`, `rms(step) <= 2e-4` and
`|delta E| < 1e-6`. This is a deliberately tightened variant of the published
criterion (Bakken and Helgaker, J. Chem. Phys. 117, 9160 (2002)), which requires
only `max(|force|)` and (`|delta E|` or `max(|step|)`); the looser form accepts
geometries whose remaining RMS force still displaces the structure.

## Examples

### Default L-BFGS

```bash
mlmm opt -i my.pdb --parm real.parm7 -l 'SAM:1' -b uma -o result_opt
```

### RFO for stiffer convergence

```bash
mlmm opt -i my.xyz --parm real.parm7 --ref-pdb topology.pdb -q -1 -m 1 --opt-mode rfo -b mace -o result_opt_rfo
```

### Pre-relax endpoints before path-opt

```bash
mlmm opt -i 1.R.pdb --parm real.parm7 -l '...' -o /tmp/relax_R
mlmm opt -i 3.P.pdb --parm real.parm7 -l '...' -o /tmp/relax_P
mlmm path-opt -i /tmp/relax_R/final_geometry.xyz /tmp/relax_P/final_geometry.xyz \
    --parm real.parm7 --ref-pdb 1.R.pdb --ref-pdb 3.P.pdb \
    -l 'SAM:1,GPP:-3' -o result_path_opt
```

## Output

```
result_opt/
├── result.json                 # when --out-json
├── final_geometry.{xyz,pdb}    # final geometry
├── optimization_trj.xyz        # macro-cycle trajectory (when --dump)
└── optimization_all_trj.xyz    # full per-step trajectory (when --dump)
```

`result.json` keys: `status` (converged / stalled / not_converged), `n_opt_cycles`,
`energy_hartree`, `final_max_force`, `files.final_geometry_xyz`.
When `--flatten` runs, `rigid_projection` also records the selected treatment,
effective rank, Hessian source, and Hessian shape.

## `--opt-mode` choice

| Mode | Algorithm | When |
|---|---|---|
| `grad` / `lbfgs` | L-BFGS | Default, fast, robust for most well-conditioned minima |
| `hess` / `rfo` | RFO with Hessian updates | Stiffer convergence; useful when L-BFGS oscillates |

## Caveats

- Not a TS optimizer — for TS use `tsopt.md`.
- The fixed constrained treatment has no effect unless `--flatten` performs
  PHVA. It has generic rank 6/3/1/0 for 0/1/2/3+
  non-collinear frozen anchors; realistic boundaries normally rank 0 and
  all-frozen input is an error. A stale non-constrained YAML value fails
  explicitly.
- L-BFGS occasionally walks past a saddle on shallow surfaces; if the
  resulting geometry has imaginary frequencies (run `freq` to check),
  re-run with `--opt-mode rfo`.
- `--config` YAML is the way to override less-common settings (step
  limits, trust radius, etc.); inspect `OPT_BASE_KW` and `LBFGS_KW`
  in `mlmm.core.defaults`.

## See also

- `tsopt.md` — TS analog.
- `freq.md` — verify the optimized minimum (zero imaginary modes).
- Defaults: `import mlmm.core.defaults as d; print(d.OPT_BASE_KW, d.LBFGS_KW, d.RFO_KW)`
