---
name: mlmm-workflows-output
description: Output parsing and multi-step workflow selection for mlmm-toolkit — `summary.json` schema, R/TS/P/IM canonical paths, `bond_changes` interpretation, energy-diagram conventions, and the cluster + 1-step / multistep / scan-list / endpoint-MEP / TS-only / DFT//ML/MM / stage-by-stage (subcommand-only, gate each stage) recipes used to extract numerical results. TRIGGER on output parsing (`summary.json`, `result.json`, `seg_NN/`), extracting barriers / ΔE / Gibbs for a paper, choosing between multi-input / scan-list / endpoint-MEP / TS-only modes, or running the pipeline subcommand-by-subcommand with a success check at each stage (instead of one `all` run). SKIP for single-subcommand syntax (CLI skill) or install / HPC questions.
---

# mlmm-toolkit Workflows and Output Parsing

## Purpose

This skill ties the per-subcommand mds in `mlmm-cli/` together
into **end-to-end recipes** and explains how to read the resulting
output trees, JSON, and figures. Use this when you have a goal
("compute the barrier of step 1 of this enzyme") and want the path
through the toolkit.

> **mlmm-specific notes**:
>
> - Every ML/MM-evaluating subcommand requires an Amber `--parm` plus
>   a layer-encoded PDB (`--detect-layer` from B-factor 0.0/10.0/20.0,
>   or explicit `--model-pdb`/`--model-indices`).
> - **Microiteration** alternates ML-region geometry steps with MM
>   relaxation (controlled by `MICROITER_KW` defaults; toggled via
>   `--microiter / --no-microiter` and `--opt-mode hess` driving the
>   outer RFO step). The MM block uses analytical-Hessian `hessian_ff`,
>   so it scales with CPU cores, not GPU.
> - The `mlmm dft` step (when `--dft` is set in `all`) computes a
>   single-point DFT energy on the **ML region only**, not the whole
>   enzyme. The MM contribution is taken from the parm7 force field.

## Six canonical workflows

### 1. Cluster + 1-step reaction (multi-input MEP)

You have R and P PDBs (from a published QM study). One step.

```bash
mlmm all -i 1.R.pdb 3.P.pdb \
    -c 'SAM,GPP,MG' -l 'SAM:1,GPP:-3' \
    --tsopt --thermo \
    -o result_mep
```

Result: `result_mep/seg_NN/{reactant,ts,product}.pdb`,
`summary.json["segments"][0]["barrier_kcal"]`.

### 2. Multi-step recursive (multi-input MEP, recursive segmentation)

You have R and P. The mechanism is multi-step, but you don't have
intermediates handy. `--refine-path` runs the recursive `path-search`,
which segments by detecting bond changes.

```bash
mlmm all -i 1.R.pdb 3.P.pdb \
    -c '...' -l '...' \
    --refine-path \
    --tsopt --thermo \
    -o result_mep
```

With `--refine-path` the output `summary.json["n_segments"]` may be > 1 —
that's the recursion finding intermediates the inputs didn't contain.
(Without it, single-pass `path-opt` yields one segment per adjacent input
pair.)

### 3. Single-input scan-list (when only R is available)

You have just the reactant. Articulate the reaction as a sequence of
distance-restraint scans.

```bash
mlmm all -i 1.R.pdb \
    -c '...' -l '...' \
    --scan-lists '[("CS1 SAM 320","GPP 321 C7",1.60)]' \
                 '[("GPP`321/H11","GLU`186/OE2",0.90)]' \
    --tsopt --thermo \
    -o result_scan
```

Each `--scan-lists` argument is one stage. See
`mlmm-cli/all-scan-list.md` for syntax details.

### 4. Endpoint-MEP with explicit intermediates

You have R, IM₁, IM₂, P from the literature. Pass them all in order:

```bash
mlmm all -i 1.R.pdb 2.IM1.pdb 3.IM2.pdb 4.P.pdb \
    -c '...' -l '...' \
    --tsopt --thermo \
    -o result_mep_4pt
```

Add `--refine-path` to run recursive sub-segmentation *between* adjacent
endpoints (then you don't have to provide every elementary step); default
single-pass `path-opt` treats the provided endpoints as the elementary
steps.

### 5. TS-only validation (existing TS candidate)

You have a TS guess from another code or a prior run. Skip
extract / path-search:

```bash
mlmm tsopt -i ts.xyz -q -1 -m 1 -b uma -o result_tsopt
mlmm freq  -i result_tsopt/final_geometry.xyz -q -1 -m 1 -b uma -o result_freq
mlmm irc   -i result_tsopt/final_geometry.xyz -q -1 -m 1 -b uma -o result_irc
```

Or use `mlmm all` with a single `-i` (collapses to TS-only
mode automatically; see `mlmm-cli/all-ts-only.md`).

### 6. DFT//ML/MM refinement

After any of the above, refine R / TS / P energies at DFT level:

```bash
mlmm dft -i seg_01/reactant.pdb \
    -l 'SAM:1,GPP:-3' \
    --func-basis 'wb97m-v/def2-tzvpd' \
    --engine gpu \
    -o dft_R
mlmm dft -i seg_01/ts.pdb       -l '...' --func-basis '...' -o dft_TS
mlmm dft -i seg_01/product.pdb  -l '...' --func-basis '...' -o dft_P
```

Composite the energies with `energy-diagram` (see below).

## Stage-by-stage execution (subcommand-only, gate each stage)

Run the pipeline as separate subcommands instead of one `mlmm all` when you want to
**judge each stage's success before spending GPU time on the next** — e.g. confirm
path-search found the right segments / bond changes before optimizing a TS, or validate
the TS (one imaginary mode + correct IRC connectivity) before thermo / DFT. `mlmm all`
runs this chain (the MEP stage is single-pass `path-opt` by default; recursive
`path-search` with `--refine-path`):

```
extract → [mm-parm] → path-opt → (per reactive seg) tsopt → freq → irc → [freq R/TS/P] → [dft] → energy-diagram
```

**mlmm carry-through**: every ML/MM-evaluating stage needs the *same* `--parm` + layer
definition (`--detect-layer` on a B-factor-encoded PDB, or `--model-pdb` / `--model-indices`)
and the *same* `-l` / `-q` / `-m`. Pass them on every command. After each stage, read its
`result.json` / `summary.json` `status` and gate before continuing.

**Stage 0 — prep** (only from a full enzyme PDB; most staged campaigns start from
already-prepared R/P cluster PDBs + parm7): generate topology, cut the pocket, encode
layers with `mm-parm` → `extract` → `define-layer` (flags in
`mlmm-cli/{mm-parm,extract,define-layer}.md`). **GATE**: `real.parm7` written and the
pocket PDB carries the intended ML atoms + layer B-factors (0/10/20).

**Stage 1 — MEP (`path-search`)**

```bash
mlmm path-search -i 1.R.pdb 3.P.pdb --parm real.parm7 --detect-layer \
    -l 'SAM:1,GPP:-3' -b uma -o ps/
```

**GATE** (`ps/summary.json`): `status == "success"`; inspect `n_segments` and EACH
segment's `bond_changes` — the intended bonds must be *formed AND broken* for the right
atoms (`mlmm-cli/bond-summary.md`, `mlmm-ts-strategy/SKILL.md`). Wrong segmentation /
spurious changes → fix chemistry or inputs **before** any TS work (don't optimize a TS
for the wrong step).

**Stage 2 — per reactive segment: TS → validate → connectivity** (seed = `ps/hei_seg_NN.xyz`)

```bash
mlmm tsopt -i ps/hei_seg_NN.xyz               --parm real.parm7 --detect-layer -l 'SAM:1,GPP:-3' -b uma -o seg_NN/tsopt
mlmm freq  -i seg_NN/tsopt/final_geometry.xyz --parm real.parm7 --detect-layer -l 'SAM:1,GPP:-3' -b uma -o seg_NN/freq
mlmm irc   -i seg_NN/tsopt/final_geometry.xyz --parm real.parm7 --detect-layer -l 'SAM:1,GPP:-3' -b uma -o seg_NN/irc
```

**GATE** in order: tsopt `result.json` `status` is `converged` (or `completed`; not
`not_converged`) → freq `result.json` `n_imaginary == 1` (exactly one imaginary frequency)
whose mode moves the reacting atoms (0 or >1 → fix via fp64 / `--coord-type dlc` /
`--flatten`, see `mlmm-ts-strategy/SKILL.md` §3, before trusting the barrier) → irc
`result.json` `status == "completed"` and forward/backward endpoints connect the **intended**
R and P (bond changes match this step). A TS that fails any gate is not this elementary step.

**Stage 3 — thermochemistry** (optional, = `all --thermo`): run `mlmm freq` on R / TS / P
for the Gibbs/QRRHO profile (`post_segments[i].gibbs_uma`).

**Stage 4 — DFT//ML/MM** (optional, = `all --dft`):

```bash
mlmm dft -i seg_NN/reactant.pdb --parm real.parm7 --detect-layer -l 'SAM:1,GPP:-3' --func-basis 'wb97m-v/def2-tzvpd' -o seg_NN/dft/R   # repeat for ts, product
```

**GATE**: each `dft/<state>/result.json` shows `"converged": true`.

**Stage 5 — energy diagram**

```bash
mlmm energy-diagram -i 0.0 -i 21.5 -i -0.7 --label-x R --label-x TS --label-x P -o diagram.png
```

> Resuming after a walltime hit uses these same commands — see `mlmm-cli/all.md`
> "Resume / restart". On any non-success status read `summary.log`, then
> `segments/seg_NN/<stage>/result.json`, before retrying. Large IRC/freq (n ≳ 4000
> atoms) on a 16 GB GPU can OOM — run pysisyphus `bofill_update` on CPU.

## summary.json schema (`mlmm all` output)

Top-level keys:

| Key | Description |
|---|---|
| `command` | The subcommand (`"all"`, `"tsopt"`, …) |
| `mlmm_toolkit_version` | Toolkit version that produced this output |
| `status` | `"success"` (all stages OK), `"partial"` (segments produced but diagrams missing), or `"failed"` |
| `charge` / `spin` | Resolved cluster charge / multiplicity |
| `environment` | `{device, gpu_name, gpu_vram_gb, cuda_version, cpu, n_cpus, ram_gb}` |
| `config` | Full effective config after CLI + YAML + defaults merge |
| `freeze_atoms` | Indices held fixed during optimization (link-H parents) |
| `n_segments` | Number of elementary steps detected |
| `n_segments_reactive` | Number of non-bridge (reactive, `kind != "bridge"`) segments |
| `rate_limiting_step` | Dict `{segment, barrier_kcal, method}` describing the highest-barrier segment (or `null` when no segments) |
| `overall_reaction_energy_kcal` | R → P total energy difference |
| `segments` | List, one per path-search segment (see below) |
| `post_segments` | List, one per post-processed segment (TS / IRC / freq / DFT details) |
| `key_output_files` | Map of role → path (mep_pdb, energy_diagrams, …) |
| `pipeline_mode` | Internal mode tag |
| `mlip_backend` | Which backend produced the energies |
| `energy_diagrams` | Paths to PNG / HTML diagrams |

Per-segment keys (`summary.json["segments"][i]`, from path-search output):

| Key | Description |
|---|---|
| `index` | Segment index (1, 2, …) |
| `tag` | Segment tag string |
| `kind` | `"seg"` or `"bridge"` (non-bridge reactive segments are `"seg"`) |
| `barrier_kcal` | MEP barrier — peak energy along the segment, kcal/mol (relative to segment reactant) |
| `delta_kcal` | MEP reaction energy — segment endpoint difference, kcal/mol |
| `bond_changes` | Multi-line **string** in the form `"Bond formed (k):\n  Cs-C : 3.17 Å -> 1.68 Å\n..."` (empty for `kind=="bridge"`). Default cutoff is 1.20× covalent radii (with internal margin 0.05) — see `mlmm-cli/bond-summary.md`. |

Per-segment keys in the post-processing list (`summary.json["post_segments"][i]`, populated when `--tsopt` / `--thermo` / `--dft` are active):

| Key | Description |
|---|---|
| `tag` | Matches the corresponding `segments[i].tag` |
| `post_dir` | `result/segments/seg_NN/` directory |
| `structures` | Map: `reactant`, `ts`, `product` → file paths |
| `irc_plot` / `irc_traj` | IRC-related artifact paths |
| `ts_imag` | `{n_imag}` |
| `uma` | `{energies_au, energies_kcal, barrier_kcal, delta_kcal, ...}` (or whichever MLIP backend was used) |
| `gibbs_uma` | `{energies, barrier_kcal, delta_kcal, ...}` (when `--thermo` is on) |
| `dft` | `{labels, energies_au, energies_kcal, diagram, structures, barrier_kcal, delta_kcal}` (when `--dft` is on) |
| `gibbs_dft_uma` | DFT//ML/MM Gibbs profile (when both `--dft` and `--thermo` are on) |
| `mep_barrier_kcal` / `mep_delta_kcal` | Plain-MEP energies (no Gibbs / DFT correction) |

## R/TS/P canonical paths

Two locations get written for each elementary step:

```
result_all/
└── segments/seg_NN/                        # CANONICAL — post-LBFGS optimized
    ├── reactant.{xyz,pdb}                  # IRC backward endpoint, then LBFGS-optimized
    ├── ts.{xyz,pdb}                        # tsopt'd transition state
    ├── product.{xyz,pdb}                   # IRC forward endpoint, then LBFGS-optimized
    └── structures/
        ├── reactant.{xyz,pdb}              # same as above (canonical) — nested copy
        ├── reactant_irc.{xyz,pdb}          # raw IRC backward end (pre-LBFGS)
        ├── ts.{xyz,pdb}                    # same as above
        ├── product.{xyz,pdb}              # same as above (canonical)
        └── product_irc.{xyz,pdb}           # raw IRC forward end (pre-LBFGS)
```

**Rule of thumb**: read from `segments/seg_NN/` for downstream stages. Use
`structures/reactant_irc.xyz` / `structures/product_irc.xyz` only when debugging
IRC vs. LBFGS divergence.

`bond_changes` are computed from `reactant.xyz` / `product.xyz`
(post-LBFGS), not from the raw IRC endpoints.

## Programmatic key extraction

```python
import json

d = json.load(open("result_all/summary.json"))

# Per-segment MEP barriers
for seg in d["segments"]:
    print(f"seg_{seg['index']:02d} ({seg['kind']}): "
          f"ΔE‡ = {seg['barrier_kcal']:.1f} kcal/mol, "
          f"ΔE = {seg['delta_kcal']:.1f} kcal/mol")

# Rate-limiting barrier (rate_limiting_step is a dict, not an int)
rls = d.get("rate_limiting_step")
if rls is not None:
    print(f"rate-limiting: seg_{rls['segment']:02d}, "
          f"barrier = {rls['barrier_kcal']:.1f} kcal/mol "
          f"({rls['method']})")

# Post-processed (tsopt / freq / dft) per-segment data
for ps in d.get("post_segments", []):
    n_imag = ps.get("ts_imag", {}).get("n_imag")
    if n_imag is not None and n_imag != 1:
        print(f"WARNING: {ps['tag']} has {n_imag} imaginary modes at TS")

# DFT//ML/MM energies (when --dft was used)
for ps in d.get("post_segments", []):
    if "dft" in ps:
        print(ps["tag"], ps["dft"]["energies_kcal"])
```

## Bond-change interpretation

`segments[i]["bond_changes"]` is a multi-line **string** assembled by
`mlmm.domain.bond_changes.summarize_changes`, e.g.:

```
Bond formed (2):
  CS1(SAM 320) -- C7(GPP 321)  : 3.17 Å -> 1.68 Å
  OE2(GLU 186) -- H11(GPP 321) : 2.41 Å -> 1.02 Å
Bond broken (2):
  CS1(SAM 320) -- S(SAM 320)    : 1.81 Å -> 3.05 Å
  C7(GPP 321) -- H11(GPP 321)  : 1.10 Å -> 2.94 Å
```

Reading rules:

- "Bond formed" / "Bond broken" headers are the structural change totals across the whole segment.
- For a single elementary step you usually expect 1–4 entries combined.
- Empty string for `kind=="bridge"` (bridge segments have no chemical bond change by construction).
- Default detection cutoff is 1.20× covalent radii (margin 0.05); see `mlmm-cli/bond-summary.md`.
- If a single segment shows > 8 bond changes, the recursive segmentation may have failed — inspect the geometries before trusting the barrier.

The flat `{"formed": [...], "broken": [...]}` dict shape is used in the **`irc` subcommand's** result, not in `path_search` / `summary.json`. Don't confuse the two.

## Failed-run output

When `summary.json["status"] != "success"`, look at:

1. `summary.log` — human-readable, prints the failure point first.
2. `segments/seg_NN/<stage>/result.json` — per-stage status (which step
   crashed).
3. `segments/seg_NN/<stage>/<stage>.log` — stack trace if any.

Even on failed runs, partial outputs are kept:

- `path_opt/seg_NN/` (`path_search/seg_NN/` under `--refine-path`) exists
  for any segment that completed the MEP stage (even if downstream stages
  failed).
- `segments/seg_NN/` is **only populated** for fully-successful
  segments.

## Energy diagrams

`mlmm all` writes `path_opt/energy_diagram_*.png` (`path_search/...` under
`--refine-path`):

- `energy_diagram_MEP.png` — bare MEP energies from the path-search
  string (MLIP, no thermochemistry).
- `energy_diagram_UMA_all.png` (etc.) — per-segment energies for
  whichever backend was used.
- `energy_diagram_G_UMA_all.png` — Gibbs free-energy diagram with
  QRRHO thermochemistry (when `--thermo`).

To compose a custom diagram from energies of multiple runs, use
`mlmm-cli/energy-diagram.md`:

```bash
mlmm energy-diagram \
    -i 0.0 -i 21.5 -i -0.7 -i 2.2 -i -18.2 \
    --label-x R --label-x TS1 --label-x IM --label-x TS2 --label-x P \
    -o my_diagram.png
```

## See also
- `mlmm-cli/all.md` and the three `all-*.md` mode files.
- `mlmm-cli/{tsopt,freq,irc,dft}.md` — per-stage
  `result.json` schemas.
- `mlmm-cli/bond-summary.md` — same bond-change algorithm,
  standalone.
- `mlmm-structure-io/SKILL.md` — input file formats that feed
  these workflows.
