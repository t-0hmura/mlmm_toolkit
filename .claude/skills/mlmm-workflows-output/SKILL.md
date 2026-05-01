---
name: mlmm-workflows-output
description: Output parsing and multi-step workflow selection for mlmm-toolkit — `summary.json` schema, R/TS/P/IM canonical paths, `bond_changes` interpretation, energy-diagram conventions, and the cluster + 1-step / multistep / scan-list / endpoint-MEP / TS-only / DFT//ML/MM recipes used to extract numerical results. TRIGGER on output parsing (`summary.json`, `result.json`, `seg_NN/`), extracting barriers / ΔE / Gibbs for a paper, or choosing between multi-input / scan-list / endpoint-MEP / TS-only modes. SKIP for single-subcommand syntax (CLI skill) or install / HPC questions.
---

# mlmm-toolkit Workflows and Output Parsing

## Overview

This skill ties the per-subcommand mds in `mlmm-cli/` together
into **end-to-end recipes** and explains how to read the resulting
output trees, JSON, and figures. Use this when you have a goal
("compute the barrier of step 1 of this enzyme") and want the path
through the toolkit.

> **mlmm-specific notes** (different from `pdb2reaction`):
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
intermediates handy. Path-search recursively segments by detecting
bond changes.

```bash
mlmm all -i 1.R.pdb 3.P.pdb \
    -c '...' -l '...' \
    --tsopt --thermo \
    -o result_mep
```

The output `summary.json["n_segments"]` may be > 1 — that's the
recursion finding intermediates the inputs didn't contain.

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

Recursive segmentation still runs *between* adjacent endpoints, so
you don't have to provide every elementary step.

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
| `n_segments_reactive` | Number of segments with non-empty bond changes |
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
| `kind` | `"reactive"` or `"bridge"` |
| `barrier_kcal` | MEP barrier — peak energy along the segment, kcal/mol (relative to segment reactant) |
| `delta_kcal` | MEP reaction energy — segment endpoint difference, kcal/mol |
| `bond_changes` | Multi-line **string** in the form `"Bond formed (k):\n  Cs-C : 3.17 Å -> 1.68 Å\n..."` (empty for `kind=="bridge"`). Default cutoff is 1.20× covalent radii (with internal margin 0.05) — see `mlmm-cli/bond-summary.md`. |

Per-segment keys in the post-processing list (`summary.json["post_segments"][i]`, populated when `--tsopt` / `--thermo` / `--dft` are active):

| Key | Description |
|---|---|
| `tag` | Matches the corresponding `segments[i].tag` |
| `post_dir` | `result/path_search/post_seg_NN/` directory |
| `structures` | Map: `reactant`, `ts`, `product` → file paths |
| `irc_plot` / `irc_traj` | IRC-related artifact paths |
| `ts_imag` | `{n_imag, frequencies_cm}` |
| `uma` | `{energies_au, energies_kcal, barrier_kcal, delta_kcal, ...}` (or whichever MLIP backend was used) |
| `gibbs_uma` | `{energies, barrier_kcal, delta_kcal, ...}` (when `--thermo` is on) |
| `dft` | `{labels, energies_au, energies_kcal, diagram, structures, barrier_kcal, delta_kcal}` (when `--dft` is on) |
| `gibbs_dft_uma` | DFT//ML/MM Gibbs profile (when both `--dft` and `--thermo` are on) |
| `mep_barrier_kcal` / `mep_delta_kcal` | Plain-MEP energies (no Gibbs / DFT correction) |

## R/TS/P canonical paths

Two locations get written for each elementary step:

```
result_all/
├── seg_NN/                                 # CANONICAL — top-level, post-LBFGS optimized
│   ├── reactant.{xyz,pdb}                  # IRC backward endpoint, then LBFGS-optimized
│   ├── ts.{xyz,pdb}                        # tsopt'd transition state
│   └── product.{xyz,pdb}                   # IRC forward endpoint, then LBFGS-optimized
└── result/path_search/post_seg_NN/structures/
    ├── reactant.{xyz,pdb}                  # same as above (canonical) — nested copy
    ├── reactant_irc.{xyz,pdb}              # raw IRC backward end (pre-LBFGS)
    ├── ts.{xyz,pdb}                        # same as above
    ├── product.{xyz,pdb}                   # same as above (canonical)
    └── product_irc.{xyz,pdb}               # raw IRC forward end (pre-LBFGS)
```

**Rule of thumb**: read from `seg_NN/` for downstream stages. Use
`reactant_irc.xyz` / `product_irc.xyz` only when debugging
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
`mlmm.bond_changes.format_changes`, e.g.:

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
2. `post_seg_NN/<stage>/result.json` — per-stage status (which step
   crashed).
3. `post_seg_NN/<stage>/<stage>.log` — stack trace if any.

Even on failed runs, partial outputs are kept:

- `path_search/seg_NN/` exists for any segment that completed
  path-search (even if downstream stages failed).
- `seg_NN/` (top-level) is **only populated** for fully-successful
  segments.

## Energy diagrams

`mlmm all` writes `path_search/energy_diagram_*.png`:

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
    --states 'R:0.0' 'TS1:21.5' 'IM:-0.7' 'TS2:2.2' 'P:-18.2' \
    -o my_diagram.png
```

## Cross-references

- `mlmm-cli/all.md` and the three `all-*.md` mode files.
- `mlmm-cli/{tsopt,freq,irc,dft}.md` — per-stage
  `result.json` schemas.
- `mlmm-cli/bond-summary.md` — same bond-change algorithm,
  standalone.
- `mlmm-structure-io/SKILL.md` — input file formats that feed
  these workflows.
