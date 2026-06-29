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
mlmm irc -i ts.{pdb,xyz,gjf} --parm real.parm7 \
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
| `--model-pdb FILE` | PDB defining the ML-region atoms (optional with `--detect-layer`) |
| `--detect-layer / --no-detect-layer` | Derive layer assignment from PDB B-factor (0.0=ML, 10.0=movable-MM, 20.0=frozen). Default on. |
| `--model-indices` | Comma-separated atom indices for ML region (e.g. `'1-50,75,100-110'`); used only when `--model-pdb` is omitted (`--model-pdb` takes precedence) |
| `--ref-pdb FILE` | Full-enzyme PDB used as topology reference for XYZ inputs |
| `--link-atom-method [scaled\|fixed]` | g-factor (default) or fixed 1.09/1.01 Å |
| `--embedcharge / --no-embedcharge` | xTB point-charge embedding for MM→ML environment (default off) |
| `-q, --charge` | Net charge; overrides `calc.charge` from YAML |
| `-l, --ligand-charge` | Per-residue charge mapping for ML region |

Inspect via `mlmm <subcommand> --help` and `mlmm <subcommand> --help-advanced`.

## Key flags

| flag | type | default | description |
|---|---|---|---|
| `-i, --input` | path | required | Optimized TS geometry |
| `-q` / `-l` / `-m` | — | — | Charge / spin (common conventions) |
| `--max-cycles` | int | 125 | Max IRC steps per branch (forward + backward) |
| `--step-size` | float | 0.10 (Bohr) | Step in Bohr; maps to `IRC_KW['step_length']` |
| `-b, --backend` | str | `uma` | MLIP backend |
| `-o, --out-dir` | path | `./result_irc/` | Output directory |
| `--config` / `--show-config` / `--dry-run` / `--help-advanced` | — | — | Standard |

## Examples

### Default IRC from a tsopt'd geometry

```bash
mlmm irc -i result_tsopt/final_geometry.xyz -q 0 -m 1 -b uma -o result_irc
```

### Tighter step / longer integration for shallow surfaces

```bash
mlmm irc -i ts.xyz -q -1 -m 1 \
    --max-cycles 250 --step-size 0.05 \
    -b uma -o result_irc_long
```

## Output

```
result_irc/
├── result.json                     # written when --out-json
├── forward_irc_trj.xyz             # raw IRC forward trajectory
├── forward_irc.pdb                 # PDB conversion when input is a PDB or --ref-pdb supplied
├── backward_irc_trj.xyz            # raw IRC backward trajectory
├── backward_irc.pdb                # PDB conversion when input is a PDB or --ref-pdb supplied
├── finished_irc_trj.xyz            # full stitched IRC trajectory (reactant -> TS -> product)
├── finished_irc.pdb                # PDB conversion when input is a PDB or --ref-pdb supplied
├── forward_last.{xyz,pdb}          # single-frame forward IRC endpoint
└── backward_last.{xyz,pdb}         # single-frame backward IRC endpoint
```

`result.json` keys:

```python
import json
d = json.load(open("result_irc/result.json"))
print(d["n_frames_forward"], d["n_frames_backward"])
print(d["energy_reactant_hartree"], d["energy_ts_hartree"], d["energy_product_hartree"])
print(d["bond_changes"])           # {"formed": [...], "broken": [...]}
print(d["status"])                  # "completed" (success path only; errors emit a separate error JSON)
```

## Forward / backward endpoints

Two forms of endpoint geometry are written:

| File | What |
|---|---|
| `forward_last.{xyz,pdb}` / `backward_last.{xyz,pdb}` | Single-frame raw IRC endpoints — **canonical** for downstream stages |
| Last frame of `forward_irc_trj.xyz` / `backward_irc_trj.xyz` | Identical to `forward_last` / `backward_last` (same final IRC frame) |

The validator and bond-change detector use `forward_last` / `backward_last`
when present. See `mlmm-workflows-output/SKILL.md`.

## Bond-change check

`bond_changes` records which bonds are different between R and P
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
