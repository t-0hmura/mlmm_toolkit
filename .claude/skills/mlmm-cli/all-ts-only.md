# `mlmm all` — TS-only mode

## When to use

You already have a **TS candidate** (typically from another QM code, an
older `mlmm-toolkit` run, or a manual guess) and want to run only the
validation + thermochemistry stages — `tsopt → irc → freq → (dft)` —
without the upstream extract / path-search.

This is the "I trust this geometry, just check it for me" mode.

## Synopsis

```bash
mlmm all -i ts_candidate.xyz \
    -q -1 -m 1 -b uma \
    --tsopt --thermo \
    [--dft --func-basis 'wb97m-v/def2-svp'] \
    -o result_ts_only
```

Or with a PDB that carries residue / charge info:

```bash
mlmm all -i ts_candidate.pdb \
    -l 'SAM:1,GPP:-3' \
    --tsopt --thermo \
    -o result_ts_only
```

## How it differs from the other two modes

`mlmm all` falls into TS-only mode when:

- exactly **one** `-i` input is given,
- **no** `--scan-lists` is provided.

The orchestrator skips path-search automatically and starts the
pipeline at `tsopt`. There is **no explicit "force TS-only" flag** — the
mode is selected purely from the input shape. If you also pass
`--no-tsopt`, `all` becomes a thin wrapper around `freq + irc + dft`
(rarely useful).

For finer control, run the underlying subcommands directly:

```bash
mlmm tsopt -i ts.xyz -q ... -m 1 -o result_tsopt -b uma
mlmm irc   -i result_tsopt/final_geometry.xyz -o result_irc -b uma
mlmm freq  -i result_tsopt/final_geometry.xyz -o result_freq -b uma
```

## Pipeline collapses to

```
ts_candidate.{xyz,pdb,gjf}
       │
       ▼
   [tsopt]            (Dimer or RS-I-RFO; default RS-I-RFO)
       │
       ▼
   [irc]              (forward + backward; LBFGS endpoint refinement)
       │
       ▼
   [freq]             (Hessian + thermo)
       │
       ▼
   [dft]              (optional)
```

`extract` and `path-search` are skipped entirely. The output tree
collapses to one segment:

```
result_ts_only/
├── summary.json
├── summary.log
├── post_seg_01/
│   ├── tsopt/         final_geometry.{xyz,pdb}, result.json
│   ├── irc/           forward_irc_trj.xyz, backward_irc_trj.xyz, finished_irc_trj.xyz
│   ├── freq/          frequencies_cm-1.txt, thermoanalysis.yaml
│   └── (dft/)
└── seg_01/
    ├── reactant.pdb   (from the IRC backward endpoint)
    ├── ts.pdb
    └── product.pdb    (from the IRC forward endpoint)
```

## Output keys

```python
import json
d = json.load(open("result_ts_only/summary.json"))
seg = d["segments"][0]
print(seg["barrier_kcal"], seg["delta_kcal"])
print(seg["bond_changes"])             # what bonds broke / formed along the IRC
print(seg["tsopt"]["n_imaginary"])     # should be 1
print(seg["irc"]["energy_reactant_hartree"], seg["irc"]["energy_ts_hartree"], seg["irc"]["energy_product_hartree"])
```

If `tsopt.n_imaginary != 1`, the geometry is **not a true first-order
saddle**; see "Distinctive failure modes" below.

## Distinctive failure modes

| Symptom | Likely cause | Fix |
|---|---|---|
| `tsopt.status == "not_converged"` | Initial Hessian misleading or step size too large | `mlmm tsopt -i ts.xyz --opt-mode rsirfo --max-cycles 200` standalone, then re-run downstream stages |
| `tsopt.n_imaginary == 0` | Geometry collapsed to a minimum during refinement | TS guess was not a real saddle; re-do `path-search` instead |
| `tsopt.n_imaginary == 2+` | Two near-degenerate negative modes | Normal for some metalloenzyme TSs; check whether the second imaginary mode is a translation / rotation residue (often resolved by tightening `freeze_atoms`) |
| `irc.bond_changes == {}` (no bonds change) | TS connects two essentially identical wells (numerical ringing) | Verify the imaginary mode visualization in `freq/`; this is sometimes a non-physical TS |

## When *not* to use TS-only mode

- You do not yet have a TS candidate. Run `path-search` (or the
  full `all` in endpoint-MEP / scan-list mode) instead.
- You have a candidate but suspect the connectivity is wrong (i.e.
  you're not sure whether your "TS" sits between the right reactant
  and product). Use `path-search` to discover the connectivity.

## Caveats

- The `--no-tsopt` flag would skip TS optimization entirely, which is
  rarely what you want in this mode.
- For an XYZ TS candidate, you must supply `-q` and `-m` explicitly
  (XYZ has no header). Use `--ref-pdb cluster.pdb` if you want
  `-l 'RES:Q'` to work.
- The IRC step here is the **canonical validation** that the TS
  connects the expected R and P. Always read `seg_01/{reactant,product}.pdb`
  to confirm the IRC ended up where you thought.

## See also

- `all.md` — base orientation.
- `tsopt.md`, `irc.md`, `freq.md`, `dft.md` — the underlying
  subcommands (which you can also run standalone if you want
  fine-grained control).
- `mlmm-workflows-output/SKILL.md` — IRC interpretation
  and bond-change conventions.
## ML/MM-aware flags (mlmm-toolkit specific)

Beyond the cluster-style flags below (inherited from `pdb2reaction`),
**`mlmm-toolkit` requires an Amber topology** and supports layer-aware
selection. Most subcommands accept:

| flag | purpose |
|---|---|
| `--parm FILE` | Amber `parm7` topology of the whole enzyme — **required** |
| `--model-pdb FILE` | PDB defining the ML-region atoms (optional with `--detect-layer`) |
| `--detect-layer / --no-detect-layer` | Pick layer assignment from PDB B-factor (0.0=ML, 10.0=movable-MM, 20.0=frozen). Default on. |
| `--model-indices` | Comma-separated atom indices for ML region (e.g. `'1-50,75,100-110'`); overrides `--model-pdb` |
| `--ref-pdb FILE` | Full-enzyme PDB used as topology reference for XYZ inputs |
| `--link-atom-method [scaled\|fixed]` | g-factor (default) or fixed 1.09/1.01 Å |
| `--embedcharge / --no-embedcharge` | xTB point-charge embedding for MM→ML environment (default off) |
| `-q, --charge` | **ML-region** charge (not whole-system) |
| `-l, --ligand-charge` | Per-residue charge mapping for ML region |

Inspect via `mlmm <subcommand> --help` and `mlmm <subcommand> --help-advanced`.
