---
name: mlmm-ts-strategy
description: Decision know-how for ML/MM enzyme reaction-barrier campaigns — precision (fp32 vs fp64 by GPU class), TS-candidate routes (path-search MEP vs distance-restrained scan), fixing wrong imaginary-frequency counts at TS-opt (fp64 / --coord-type dlc), reading a barrier when the scan started from Product (reverse direction), staged vs concerted scans, and the same-atom-set rule for controlled mutant-vs-WT comparisons (B-factor layer transplant + --detect-layer). TRIGGER on "barrier", "imaginary frequency", "wrong saddle", "fp64 / precision", "MEP vs restraint", "scan from product", "staged vs concerted scan", "mutant comparison", or "controlled experiment". SKIP for install / pure structure-format editing / MCP-transport questions.
---

# mlmm ts-strategy

Cross-cutting decisions for getting a *correct* reaction barrier out of an ML/MM
ONIOM campaign. Every flag below is verified against `mlmm/cli/common_options.py`,
`mlmm/workflows/{scan,opt,path_search,all}.py`, and `mlmm/core/defaults.py`.

## 1. Precision: pick by GPU class

| Hardware | Flag | Why |
|---|---|---|
| HPC datacenter GPU (H100 / H200 / A100) | `--precision fp64` | Deterministic-grade, low numerical noise; native fp64 throughput is affordable. Stabilises TS-opt / Hessian. |
| Consumer GPU (RTX 50xx / 40xx) | `--precision fp32` (default) | fp64 is much slower on consumer cards; fp32 is the speed/screening baseline. |

- `--precision` = `click.Choice(['fp32','fp64'])`, case-insensitive (`common_options.py` `add_precision_option`). Option default `None`; effective default `fp32` from `defaults.py` `MLMM_CALC_KW['uma_precision']='fp32'`.
- Backend routing: `uma`→`uma_precision`; `orb`→`orb_precision` (`float32-high`|`float64`); `mace`→`mace_dtype`; `aimnet2`→fp32 is a no-op and **fp64 is rejected** (inputs cast to float32 upstream).
- Accepted on `sp`, `opt`, `tsopt`, `freq`, `irc`, `scan`/`scan2d`/`scan3d`, `path-opt`, `path-search`, `all`.
- fp64 *reduces* GPU reduction-order drift but does **not** make a run bit-identical — only `--deterministic` does (see `reproducibility.md`).

## 2. Two routes to a TS candidate

| Route | Subcommand | Mechanism | Use when |
|---|---|---|---|
| (a) MEP / path-search | `path-search` (or `path-opt` for one segment) | Recursive GSM/DMF segmentation; brackets the TS between endpoints, auto-bridges gaps, one TS per segment | You have R (and optionally P / intermediates) and want the path discovered |
| (b) Distance-restrained build-up | `scan` (`scan2d`/`scan3d`) | Harmonic restraint `E = ½·k·(r_ij − target)²` (default `k=10`, `restraints.py` `HarmonicBiasCalculator`) drives the reacting distance(s) toward the barrier with L-BFGS relaxation | No usable second endpoint / TS guess — drive the reacting bond directly |

- There is **no `opt --restraint` flag**. The restrained-optimization route IS the `scan` subcommand (`scan.py` docstring "staged bond-length scan with harmonic restraints"); plain `opt` is un-restrained.
- `path-search` (`app.py`: "Search reaction pathways recursively.") auto-segments a multistep path; `path-opt` optimizes a single given segment.
- Feed a TS candidate from either route into `tsopt → irc → freq` (or `all --tsopt`).

## 3. Wrong imaginary-frequency count at TS-opt

A clean first-order saddle = **exactly one** dominant imaginary mode along the reaction coordinate.

| Symptom | Action |
|---|---|
| Spurious 2nd small imaginary mode, OR no dominant reaction mode | (i) raise precision: `--precision fp64`; AND/OR (ii) switch coordinates: `--coord-type dlc` |
| Still no clean saddle | combine both; verify the imaginary-mode visualization in `freq/` actually moves the reacting atoms |

- `--coord-type` = `click.Choice(['cart','redund','dlc','tric'])`, case-insensitive; `'dlc'` is present verbatim (`common_options.py` `add_coord_type_option`). Effective default `'cart'` (`defaults.py` `GEOM_KW_DEFAULT['coord_type']='cart'`).
- `dlc` = delocalized internal coordinates: slower but more robust convergence on torsion-rich systems.
- `dlc` requires a **Hessian-based optimizer**: in `opt.py`, `--coord-type dlc` with L-BFGS (`--opt-mode grad`) is forced back to `cart` with a warning. Use it on `tsopt` (RFO/RS-I-RFO) or `opt --opt-mode hess`.
- `path-opt`/`path-search` restrict choices to `('cart','dlc')` (pysisyphus ChainOfStates supports only those).
- mlmm caveat (help text): `DLC + link atom` and `DLC + 3-layer frozen MM` are numerically unverified — `cart` is the published-numbers default.

## 4. Reading the barrier when the scan started from Product

If the scan/path **starts from P**, the raw reported barrier is the **reverse** direction.

| Quantity | Formula |
|---|---|
| Forward barrier | `E(TS) − E(reactant)` |
| Reverse barrier (raw P-start number) | `E(TS) − E(product)` |

- This is a *read-time interpretation*, **not a CLI flag**. Always confirm which endpoint is R vs P (read `segments/seg_NN/{reactant,product}.pdb` from the IRC, not the scan direction).
- Concrete case (CM): the scan begins at P, so the campaign's "barrier" is reverse; the forward barrier = `E(TS) − E(R)`.

## 5. Staged vs concerted scan

`-s/--scan-lists` is `multiple=True` (`scan.py`). Help: "Multiple inline literals define sequential stages."

| Form | Invocation | Meaning | Needs mechanism up front? |
|---|---|---|---|
| Concerted | **single** `--scan-lists` literal with several `(i,j,target)` tuples | all coords driven together in one stage | No |
| Staged | **repeat** `--scan-lists` (one literal per stage) | each stage = its own restrained relaxation, written to `stage_NN/` | Yes — define the mechanism per stage |

```bash
# Concerted (one stage, two coords driven together):
mlmm scan -i r.pdb --parm e.parm7 -l 'LIG:Q' \
    --scan-lists '[(1,5,1.40),(7,9,1.60)]' -o result_concerted

# Staged (two sequential stages):
mlmm scan -i r.pdb --parm e.parm7 -l 'LIG:Q' \
    --scan-lists '[(1,5,1.40)]' \
    --scan-lists '[(7,9,0.95)]' -o result_staged
```

- `path-search` does multistep auto-segmentation, so a **concerted** scan needs no mechanism breakdown.
- A **staged** scan needs the mechanism defined up front, **but when the mechanism is known, staged gives cleaner per-step control and is generally preferred.**
- A 4-tuple expands into 2 stages (bidirectional scan).

## 6. Controlled mutant-vs-WT (or mechanism-vs-mechanism) comparison

**Rule: every compared model MUST use the SAME atom set — identical atom count and residues — or it is not a controlled experiment.** A different atom set changes the energy reference and invalidates the ΔE‡ comparison; a geometrically re-derived ML/movable/frozen partition on the mutant also produces spurious soft modes (`tsopt.n_imaginary ≥ 2`, both tiny → IRC aborts).

mlmm recipe (transplant the WT ML/MM layer encoding onto the mutant; run with `--detect-layer`):

| step | action |
|---|---|
| 1 | Build the mutant structure; keep the same residue set as WT (only the mutated residue's identity differs). |
| 2 | Transplant WT's B-factor layer encoding onto the mutant by `(resid, atom-name)`: **ML=0.0, MovableMM=10.0, FrozenMM=20.0** (`common_options.py` `add_ml_layer_detection_options`). |
| 3 | Run with `--detect-layer` (default on) so the SAME layer assignment is reused. |

```bash
mlmm all -i mutant_layered.pdb -l 'LIG:Q' \
    --tsopt True --thermo True [--dft True --dft-func-basis '...'] \
    -o result_mutant
```

| flag | do | why |
|---|---|---|
| `--detect-layer` | KEEP (default `True`) | reads transplanted B-factors → ML/movable/frozen byte-identical to WT |
| `-c/--center`, `-r/--radius` | **OMIT** | geometric extraction would re-derive a *different* pocket on the mutant (`all.py`: omitting `-c` skips extraction, uses the full structure as-is) |
| `-l/--ligand-charge` | KEEP | a non-standard ligand's charge isn't in the standard-AA table; `-l 'RES:Q'` auto-derives the total. Prefer over hardcoded `-q`. |
| `--movable-cutoff` | do NOT pass | it **disables** `--detect-layer` (`scan.py`) |

- Same principle applies to comparing two mechanisms on the same enzyme: identical atom set across both, only the reaction coordinate differs.

## See also

- `mlmm-cli/tsopt.md`, `scan.md`, `path-search.md`, `define-layer.md` — per-subcommand flags.
- `mlmm-cli/all-ts-only.md` — the full mutate→complete→transplant→run mutation recipe.
- `mlmm-workflows-output/SKILL.md` — IRC R/TS/P canonical paths and bond-change conventions.
- `mlmm-hpc/SKILL.md` — choosing the GPU class that determines §1.
- docs `reproducibility.md` — fp64 vs `--deterministic`.
