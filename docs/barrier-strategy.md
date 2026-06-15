# Reaction-barrier strategy

Practical decisions for getting a *correct* reaction barrier out of an ML/MM
ONIOM campaign: choosing precision for your GPU, building a transition-state
candidate, fixing a transition state that has the wrong number of imaginary
frequencies, reading a barrier that was scanned from the product side, choosing
between staged and concerted scans, and setting up a *controlled* comparison
between a mutant and the wild type.

---

## 1. Precision: choose by GPU class

`--precision` selects the MLIP backend floating-point precision (`fp32` or
`fp64`, case-insensitive). The effective default is `fp32`.

| Hardware | Recommended | Reasoning |
| --- | --- | --- |
| HPC datacenter GPU (H100 / H200 / A100) | `--precision fp64` | Deterministic-grade, low numerical noise; native fp64 throughput is affordable. Stabilises TS optimization and the Hessian. |
| Consumer GPU (RTX 50xx / 40xx) | `--precision fp32` (default) | fp64 is markedly slower on consumer cards. fp32 is the speed/screening baseline. |

```bash
# Datacenter H200 — full-precision base inference
mlmm tsopt -i ts.pdb --parm enzyme.parm7 -l 'LIG:Q' -b uma --precision fp64 -o result_ts

# Consumer RTX — fast screening with the default
mlmm scan -i r.pdb --parm enzyme.parm7 -l 'LIG:Q' -b uma --scan-lists '[(1,5,1.4)]' -o result_scan
```

`--precision` is accepted on every compute subcommand: `sp`, `opt`, `tsopt`,
`freq`, `irc`, `scan`/`scan2d`/`scan3d`, `path-opt`, `path-search`, and `all`.
The value is routed per backend (UMA precision, ORB precision, MACE
`default_dtype`).

```{note}
**AIMNet2 exception.** For `-b aimnet2`, `fp32` is a no-op and `fp64` is
*rejected* — its model inputs are cast to float32 upstream, so an "fp64" run
could not actually run in fp64. Use `uma`, `orb`, or `mace` when you need fp64.
```

```{note}
`--precision fp64` *reduces* GPU reduction-order drift but does **not** make a
run bit-identical. Only `--deterministic` gives bit-exactness — see
[Reproducibility](reproducibility.md).
```

## 2. Two routes to a transition-state candidate

There is no single "find the TS" button. Pick the route that matches the
information you already have.

**(a) MEP / path search** — when you have a reactant (and optionally a product
or intermediates) and want the path *discovered*. `path-search` runs a recursive
GSM/DMF minimum-energy-path search, brackets the TS between endpoints, bridges
gaps between segments, and emits one TS per segment. (`path-opt` optimizes a
single supplied segment.)

```bash
mlmm path-search -i r.pdb p.pdb --parm enzyme.parm7 -l 'LIG:Q' -o result_mep
```

**(b) Distance-restrained build-up** — when you have neither a usable second
endpoint nor a TS guess, drive the reacting bond directly. The `scan` subcommand
*is* the distance-restrained optimization route: it adds a harmonic restraint
`E = ½·k·(r_ij − target)²` to each reacting pair and relaxes everything else with
L-BFGS, walking the reacting distance toward the barrier.

```bash
mlmm scan -i r.pdb --parm enzyme.parm7 -l 'LIG:Q' \
    --scan-lists '[(1,5,1.40)]' -o result_scan
```

```{note}
There is no `opt --restraint` flag. Plain `opt` is an *un-restrained* optimizer;
the restrained build-up of a TS candidate is done with `scan` (drive the
distance) or with `path-search` (route a).
```

Feed the resulting candidate into `tsopt → irc → freq` (or `mlmm all --tsopt`)
to optimize and validate it.

## 3. Wrong number of imaginary frequencies at the TS

A clean first-order saddle has **exactly one** dominant imaginary mode along the
reaction coordinate. Two common failures are a spurious second small imaginary
mode, or no dominant reaction mode at all.

| Symptom | Fix |
| --- | --- |
| Spurious 2nd small imaginary mode, or no dominant reaction mode | Raise precision with `--precision fp64`, **and/or** switch coordinates with `--coord-type dlc`. |
| Still no clean saddle | Combine both, then verify in `freq/` that the imaginary mode actually moves the reacting atoms. |

```bash
mlmm tsopt -i ts_guess.pdb --parm enzyme.parm7 -l 'LIG:Q' -b uma \
    --precision fp64 --coord-type dlc -o result_ts
```

`--coord-type` selects the optimization coordinate system
(`cart` | `redund` | `dlc` | `tric`; default `cart`). `dlc` (delocalized
internal coordinates) is slower but converges more robustly on torsion-rich
systems and toward a clean first-order saddle.

```{warning}
`--coord-type dlc` needs a **Hessian-based** optimizer. On `opt` with the
default L-BFGS (`--opt-mode grad`) it is silently forced back to `cart`; use it
on `tsopt` (RFO / RS-I-RFO) or `opt --opt-mode hess`. `path-opt`/`path-search`
accept only `cart` and `dlc`. mlmm caveat: `DLC + link atom` and
`DLC + 3-layer frozen MM` are numerically unverified, so `cart` remains the
default behind the published numbers.
```

## 4. Reading a barrier scanned from the product side

The barrier you read depends on which endpoint the scan started from. If the
scan (or path) **starts at the product**, the raw reported barrier is the
**reverse** direction.

| Quantity | Formula |
| --- | --- |
| Forward barrier | `E(TS) − E(reactant)` |
| Reverse barrier (the raw product-start number) | `E(TS) − E(product)` |

This is a *read-time* interpretation, not a CLI flag. Always confirm which
endpoint is the reactant versus the product by reading
`segments/seg_NN/{reactant,product}.pdb` from the IRC, rather than trusting the
scan direction. For a product-start campaign, the forward barrier you want is
`E(TS) − E(reactant)`, not the number printed against the product start.

## 5. Staged versus concerted scans

`-s/--scan-lists` may be repeated. The number of times you pass it decides
whether the coordinates are driven together (concerted) or in sequence (staged).

| Form | How to invoke | Meaning | Mechanism needed up front? |
| --- | --- | --- | --- |
| Concerted | a **single** `--scan-lists` literal containing several `(i, j, target)` tuples | all coordinates driven together within one stage | No |
| Staged | **repeat** `--scan-lists`, one literal per stage | each stage is its own restrained relaxation, written to `stage_NN/` | Yes — define the mechanism per stage |

```bash
# Concerted: one stage, two distances driven together
mlmm scan -i r.pdb --parm enzyme.parm7 -l 'LIG:Q' \
    --scan-lists '[(1,5,1.40),(7,9,1.60)]' -o result_concerted

# Staged: two sequential stages
mlmm scan -i r.pdb --parm enzyme.parm7 -l 'LIG:Q' \
    --scan-lists '[(1,5,1.40)]' \
    --scan-lists '[(7,9,0.95)]' -o result_staged
```

A concerted scan needs no mechanism breakdown — `path-search` performs the
multistep auto-segmentation for you. A staged scan needs the mechanism defined
up front, **but when the mechanism is known, staged scans give cleaner per-step
control and are generally preferred.** (A four-tuple expands into two stages for
a bidirectional scan.)

## 6. Controlled mutant-vs-WT (or mechanism-vs-mechanism) comparison

```{important}
For a mutant-versus-wild-type (or mechanism-versus-mechanism) barrier
comparison, **all compared models must use the same atom set — identical atom
count and residues.** Otherwise the energy reference differs and the comparison
is not a controlled experiment. A geometrically re-derived ML/movable/frozen
partition on the mutant also produces spurious soft modes
(`tsopt.n_imaginary ≥ 2`, both tiny → the IRC aborts).
```

In mlmm, preserve the wild-type ML/MM layering by **transplanting the WT
B-factor layer encoding onto the mutated structure** and running with
`--detect-layer`:

| Step | Action |
| --- | --- |
| 1 | Build the mutant; keep the same residue set as WT (only the mutated residue's identity differs). |
| 2 | Copy WT's per-atom B-factor layer codes onto the mutant by `(resid, atom-name)`: **ML = 0.0, MovableMM = 10.0, FrozenMM = 20.0**. |
| 3 | Run with `--detect-layer` (default on) so the *same* layer assignment is reused. |

```bash
mlmm all -i mutant_layered.pdb -l 'LIG:Q' \
    --tsopt True --thermo True \
    -o result_mutant
```

| Flag | Action | Why |
| --- | --- | --- |
| `--detect-layer` | keep (default `True`) | reads the transplanted B-factors so the ML / movable / frozen layers are byte-identical to WT |
| `-c/--center`, `-r/--radius` | **omit** | geometric extraction would re-derive a *different* pocket on the mutant; omitting `-c` skips extraction and uses the structure as-is |
| `-l/--ligand-charge` | keep | a non-standard ligand's charge is not in the standard amino-acid table; `-l 'RES:Q'` auto-derives the total charge (preferred over a hardcoded `-q`) |
| `--movable-cutoff` | do not pass | it disables `--detect-layer` |

The same-atom-set principle applies equally when comparing two mechanisms on the
same enzyme: keep an identical atom set across both models and vary only the
reaction coordinate.

## See also

- [`tsopt`](tsopt.md) · [`scan`](scan.md) · [`path-search`](path-search.md) · [`define-layer`](define-layer.md) — per-subcommand references.
- [Reproducibility](reproducibility.md) — `--precision fp64` versus `--deterministic`.
- [Device & HPC Setup](device-hpc.md) — choosing the GPU class that drives §1.
