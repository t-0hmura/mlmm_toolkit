# `mlmm irc`

## Purpose

Intrinsic Reaction Coordinate (IRC) integration from a TS geometry.
Default integrator: **EulerPC** (mass-weighted Cartesians). Forward
and backward branches are run, then each endpoint is LBFGS-optimized
to the nearest minimum. Output: a stitched IRC trajectory plus the
optimized R and P geometries.

## Synopsis

```bash
mlmm irc -i ts.{pdb,xyz,gjf} \
    [-q 0 -m 1] [-l 'RES:Q,...'] \
    [--max-cycles 125] [--step-size 0.1] \
    [-b uma|orb|mace|aimnet2] [-o ./result_irc/]
```

## Key flags

| flag | type | default | description |
|---|---|---|---|
| `-i, --input` | path | required | Optimized TS geometry |
| `-q` / `-l` / `-m` | — | — | Charge / spin (common conventions) |
| `--max-cycles` | int | 125 | Max IRC steps per branch (forward + backward) |
| `--step-size` | float | (live default) | Step in Bohr; check `IRC_KW.step_size` |
| `--integrator` | str | `EulerPC` | EulerPC (default) or other pysisyphus integrator |
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
├── result.json
├── forward_irc_trj.xyz             # raw IRC forward trajectory
├── backward_irc_trj.xyz            # raw IRC backward trajectory
├── finished_irc_trj.xyz            # stitched: backward (reversed) + TS + forward
├── reactant.{xyz,pdb}              # backward endpoint after LBFGS
├── product.{xyz,pdb}               # forward endpoint after LBFGS
├── reactant_irc.xyz                # raw IRC backward end (pre-LBFGS)
├── product_irc.xyz                 # raw IRC forward end (pre-LBFGS)
└── irc.log
```

`result.json` keys:

```python
import json
d = json.load(open("result_irc/result.json"))
print(d["n_frames_forward"], d["n_frames_backward"])
print(d["energy_reactant_hartree"], d["energy_ts_hartree"], d["energy_product_hartree"])
print(d["bond_changes"])           # {"formed": [...], "broken": [...]}
print(d["status"])                  # "completed" / "diverged" / ...
```

## R/TS/P canonical geometries

Two flavours of endpoint geometry are written:

| File | What |
|---|---|
| `reactant.xyz`, `product.xyz` | Post-LBFGS-optimized minima — **canonical** for downstream stages |
| `reactant_irc.xyz`, `product_irc.xyz` | Raw IRC endpoint — useful for debugging IRC vs. LBFGS divergence |

The validator and bond-change detector use `reactant.xyz` / `product.xyz`
when present, falling back to `_irc.xyz` if not. See
`mlmm-workflows-output/SKILL.md`.

## Bond-change check

`bond_changes` records which bonds are different between R and P
according to a 1.3× covalent-radius cutoff. This is the same algorithm
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
- Defaults: `import mlmm.defaults as d; print(d.IRC_KW)`
