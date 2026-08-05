# `mlmm irc`

## Purpose

Intrinsic Reaction Coordinate (IRC) integration from a TS geometry.
Default integrator: **EulerPC** (mass-weighted Cartesians). Forward
and backward branches are run; `forward_last` / `backward_last` are the
last raw IRC frames (the IRC endpoints), not optimized minima. Output:
a stitched IRC trajectory plus the forward/backward endpoint geometries.
Run `mlmm opt` separately to relax the endpoints to true minima.

## Synopsis

```bash
mlmm irc -i ts.{pdb,cif,mmcif,xyz} --parm real.parm7 \
    [-q 0 -m 1] [-l 'RES:Q,...'] \
    [--max-cycles 125] [--step-size 0.1] \
    [-b uma|orb|mace|aimnet2] [-o ./result_irc/]
```


## ML/MM-aware flags (mlmm-toolkit specific)

In addition to the common flags below,
**`mlmm-toolkit` requires an Amber topology** and supports layer-aware
selection. Most subcommands accept:

| flag | purpose |
|---|---|
| `--parm FILE` | Amber `parm7` topology of the whole enzyme — required unless provided in YAML as `calc.real_parm7` |
| `--model-pdb FILE` | Explicit ML-region PDB; takes precedence over B-factor ML membership |
| `--detect-layer` | Automatically read B-factor layers; explicit ML membership retains valid movable/frozen MM layers. Enabled by default. |
| `--model-indices` | Explicit ML atom indices used when `--model-pdb` is omitted; takes precedence over B-factor ML membership |
| `--ref-pdb FILE` | Full-enzyme PDB/mmCIF used as topology reference for XYZ inputs |
| `--link-atom-method [scaled\|fixed]` | g-factor (default) or fixed 1.09/1.01 Å |
| `--embedcharge / --no-embedcharge` | Unavailable in v0.3.3; use `--no-embedcharge` |
| `-q, --charge` | Net charge; overrides `calc.model_charge` from YAML |
| `-l, --ligand-charge` | Per-residue charge mapping for ML region |

Inspect via `mlmm <subcommand> --help` and `mlmm <subcommand> --help-advanced`.

## Key flags

| flag | type | default | description |
|---|---|---|---|
| `-i, --input` | path | required | Optimized TS geometry |
| `-q` / `-l` / `-m` | — | — | Charge / spin (common conventions) |
| `--max-cycles` | int | 125 | Max IRC steps per branch (forward + backward) |
| `--step-size` | float | 0.10 (Bohr) | Step in Bohr; maps to `IRC_KW['step_length']` |
| `--never-stop / --no-never-stop` | bool | off | Ignore gradient and energy endpoint criteria and trace to max cycles; propagation failures still stop |
| `--read-hess` | path | — | Identified NPZ from `freq --dump-hess`; geometry, atom order, active-DOF basis, and schema-2 charge/multiplicity must match |
| `--allow-unverified-hess-state` | bool | off | Permit a schema-1 Hessian whose charge/multiplicity cannot be verified. Requires `--read-hess` and independent state checking; schema-2 mismatches remain fatal. |
| `--workers` | int | 1 | UMA predictor workers; `>1` requires `fairchem-core[extras]` and is incompatible with `Analytical` |
| `-b, --backend` | str | `uma` | MLIP backend |
| `-o, --out-dir` | path | `./result_irc/` | Output directory |
| `--config` / `--show-config` / `--dry-run` / `--help-advanced` | — | — | Standard |

## Examples

### Default IRC from a tsopt'd geometry

```bash
mlmm irc -i result_tsopt/final_geometry.xyz --parm real.parm7 \
    --ref-pdb enzyme_layered.pdb -q 0 -m 1 -b uma -o result_irc
```

### Tighter step / longer integration for shallow surfaces

```bash
mlmm irc -i ts.xyz --parm real.parm7 --ref-pdb enzyme_layered.pdb -q -1 -m 1 \
    --max-cycles 250 --step-size 0.05 \
    -b uma -o result_irc_long
```

If a branch stops immediately, reduce `--step-size` first. Use
`--never-stop` when tracing to the maximum-cycle guard is intended. Numerical
or integration failure can still stop the branch.

## Output

```
result_irc/
├── result.json                     # written when --out-json
├── forward_irc_trj.xyz             # raw IRC forward trajectory
├── forward_irc.pdb                 # PDB companion when topology + conversion are available
├── forward_irc.cif                 # bridge-input companion with restored IDs
├── backward_irc_trj.xyz            # raw IRC backward trajectory
├── backward_irc.pdb                # PDB companion (same gating)
├── backward_irc.cif                # bridge-input companion with restored IDs
├── finished_irc_trj.xyz            # full stitched path (first endpoint -> TS -> last endpoint)
├── finished_irc.pdb                # PDB companion (same gating)
├── finished_irc.cif                # bridge-input companion with restored IDs
├── forward_last.{xyz,pdb,cif}      # single-frame forward IRC endpoint/companions
└── backward_last.{xyz,pdb,cif}     # single-frame backward IRC endpoint/companions
```

With a non-empty YAML `irc.prefix`, EulerPC inserts one underscore before each
filename (`prefix: trial` → `trial_finished_irc_trj.xyz`); read the normalized
names from `result.json.files`.

`result.json` keys:

```python
import json
d = json.load(open("result_irc/result.json"))
print(d["n_frames_forward"], d["n_frames_backward"])
print(d["energy_first_hartree"], d["energy_ts_hartree"], d["energy_last_hartree"])
print(d.get("bond_changes"))       # directed first -> last; may be omitted
print(d["status"])                  # "completed" (success path only; errors emit a separate error JSON)
print(d["never_stop"], d["never_stop_energy_bypasses"])
print(d["rigid_projection"]["electronic_state_verified"])  # False only for an opted-in schema-1 handoff
print(d["rigid_projection"]["treatment"], d["rigid_projection"]["effective_rank"])
```

Schema-2 Hessian handoffs fail closed on model charge or multiplicity
mismatch. Schema-1 files can be used only with
`--allow-unverified-hess-state`; this bypasses missing identity metadata, not a
known mismatch.

Standalone IRC does not know which endpoint is the chemical reactant or
product. Read `energy_first_hartree` / `energy_last_hartree`; the older
`energy_reactant_hartree` / `energy_product_hartree` keys are compatibility
aliases for first/last only. Assign R/P after inspecting or matching the
endpoint structures. `never_stop` records whether the opt-in mode was enabled;
`never_stop_energy_bypasses` is the observed bypass count.

The default `constrained` treatment removes only full-system rigid motions
that leave frozen anchors fixed. Generic ranks are 6/3/1/0 for
zero/one/two/at least three non-collinear anchors, and realistic ML/MM
boundaries normally have rank 0. All-frozen input is an explicit error.
A stale non-constrained YAML value fails explicitly. `result.json` records the
treatment, effective rank, initial-Hessian source, and Hessian shape.

## Forward / backward endpoints

Two forms of endpoint geometry are written:

| File | What |
|---|---|
| `forward_last.{xyz,pdb,cif}` / `backward_last.{xyz,pdb,cif}` | Single-frame raw IRC endpoints — **canonical** for downstream stages; companions depend on topology/bridge metadata |
| Last frame of `forward_irc_trj.xyz` / `backward_irc_trj.xyz` | Identical to `forward_last` / `backward_last` (same final IRC frame) |

The validator and bond-change detector use `forward_last` / `backward_last`
when present. Their direction is not a chemical R→P assignment. See
`mlmm-workflows-output/SKILL.md`.

## Bond-change check

`bond_changes` records the directed difference between the first and last
standalone IRC endpoints
according to a 1.20× covalent-radius cutoff (`bond_factor` default). This is the same algorithm
used by `bond-summary` and `path-search` segmentation.

```python
import json
bc = json.load(open("result_irc/result.json"))["bond_changes"]
for b in bc["formed"]: print("FORMED ", b)
for b in bc["broken"]: print("BROKEN ", b)
```

## Caveats

- IRC starts from a **single imaginary mode** TS. If `tsopt` produced
  multiple imaginary modes, IRC may follow the wrong one — re-tsopt
  first.
- `--max-cycles 125` is enough for most clusters. If forward / backward
  hits the cap, the surface is probably very shallow; try a smaller
  `--step-size`.
- The bond-change detector is geometry-based (covalent-radius cutoff),
  not physics-based. Metal–ligand bonds may flicker on the borderline.

## See also

- `tsopt.md` — produces the IRC starting geometry.
- `freq.md`, `dft.md` — downstream.
- `bond-summary.md` — same bond-change algorithm, standalone.
- `mlmm-workflows-output/SKILL.md` — R/TS/P path conventions.
- Defaults: `import mlmm.core.defaults as d; print(d.IRC_KW)`
