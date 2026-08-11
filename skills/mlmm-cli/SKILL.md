---
name: mlmm-cli
description: Per-subcommand reference for mlmm-toolkit's 22 CLI subcommands (extract / mm-parm / define-layer / opt / tsopt / freq / irc / dft / scan / path-search / all / …). SKILL.md is the orientation + cross-cutting flag conventions + canonical recipes; each subcommand has its own md (extract.md / mm-parm.md / tsopt.md / …) for flags, validation, and caveats. TRIGGER on questions about a specific subcommand, flag, or shell invocation. SKIP for install / HPC / output-parsing / structure-format-editing / overview questions.
---

# mlmm CLI

## Subcommand index

Each row points to the full per-subcommand md in this skill directory.

| md | subcommand | role (2 lines) |
|---|---|---|
| `all.md` | `all` | End-to-end pipeline: extract → MEP → TS → IRC → freq → (DFT) in one invocation.<br>Delegates to a base orientation; specific modes are in `all-{endpoint-mep,scan-list,ts-only}.md`. |
| `all-endpoint-mep.md` | `all` (mode 1) | Drives the pipeline from N reaction-ordered structures (R, optionally IM₁ … IMₙ, P).<br>Path search runs GSM/DMF between adjacent endpoints; recursion handles multi-step mechanisms. |
| `all-scan-list.md` | `all` (mode 2) | Drives the pipeline from a single reactant + a list of staged distance scans.<br>The scan list seeds the MEP; recursion handles intermediate states like in mode 1. |
| `all-ts-only.md` | `all` (mode 3) | Skips path search and starts from a TS candidate; runs `tsopt → irc`, with freq/DFT enabled by their flags.<br>Use when you already have a transition-state guess (from a different code or a prior run). |
| `extract.md` | `extract` | Selects and writes an active-site/model pocket around the substrate residues.<br>`define-layer` separately assigns B-factor layers and frozen atoms. |
| `mm-parm.md` | `mm-parm` | Generate Amber `parm7` + `rst7` from a PDB via tleap (and antechamber for non-standard ligands).<br>Required for any subcommand that needs MM gradients. |
| `define-layer.md` | `define-layer` | Assign / refine ML / movable-MM / frozen layers via the PDB B-factor field.<br>Standalone or post-`extract` adjustment without rebuilding parm7. |
| `oniom-export.md` | `oniom-export` | Export the layered system as a Gaussian g16 ONIOM input (or ORCA).<br>Useful for input-deck exchange or hand-comparing setup with a third-party DFT/MM run. |
| `oniom-import.md` | `oniom-import` | Reverse direction: read a g16 / ORCA ONIOM input and reconstruct an `mlmm-toolkit` PDB.<br>Use when adopting an existing Gaussian ONIOM workflow. |
| `path-search.md` | `path-search` | Recursive MEP search (GSM or DMF) across N endpoints with bond-change segmentation.<br>Splits multi-step paths into one-TS-per-segment automatically. |
| `path-opt.md` | `path-opt` | MEP optimization for a **single** segment between two endpoints.<br>Building block of `path-search`; also useful for refining one segment without re-running the whole search. |
| `opt.md` | `opt` | Single-structure geometry optimization with L-BFGS or RFO.<br>`--opt-mode grad` (L-BFGS, default) is fast; `--opt-mode hess` (RFO) is robust on tricky surfaces. |
| `tsopt.md` | `tsopt` | TS optimization: default RS-I-RFO (`--opt-mode hess/rsirfo`); Hessian-Guided Dimer is the lighter alternative (`--opt-mode grad/dimer`). |
| `freq.md` | `freq` | Vibrational analysis: Hessian, frequencies, normal-mode visualization, QRRHO thermochemistry.<br>Default temperature/pressure 298.15 K / 1 atm; partial-Hessian variant when `freeze_atoms` is non-empty. |
| `sp.md` | `sp` | ONIOM single-point energy + forces (and optional Hessian).<br>Cheapest stage; useful for spot-checking a geometry without running an optimization. |
| `irc.md` | `irc` | IRC integration with EulerPC in mass-weighted Cartesians.<br>Writes raw forward/backward endpoints; optimize them separately with `opt` or through `all`. |
| `dft.md` | `dft` | Single-point DFT through PySCF (CPU) or GPU4PySCF (CUDA, x86_64).<br>`--engine gpu` is the default; if the GPU backend is unavailable it raises an error — select CPU explicitly with `--engine cpu`. |
| `scan.md` | `scan` | 1D distance scan with harmonic restraints to seed a path search.<br>Useful when neither endpoint nor TS guess is available — drives the bond manually. |
| `scan2d.md` | `scan2d` | 2D analog of `scan` with two restrained distances.<br>Generates a grid; mlmm interpolates the MEP through the grid minima. |
| `scan3d.md` | `scan3d` | 3D analog with three restrained distances.<br>Rare but supported; output volume grows quickly, plan resources. |
| `trj2fig.md` | `trj2fig` | Plot an energy profile from an XYZ trajectory.<br>Reads ASE-style energies in the comment line and writes a figure or CSV (PNG/JPEG/HTML/SVG/PDF/CSV). |
| `energy-diagram.md` | `energy-diagram` | Build an ad-hoc energy diagram from a list of state names + energies.<br>For composing diagrams that combine multiple `mlmm-toolkit` runs. |
| `add-elem-info.md` | `add-elem-info` | Repair / add the element column (PDB cols 77-78).<br>Run before `extract` if your PDB came out of PyMOL or Maestro and elements are missing. |
| `fix-altloc.md` | `fix-altloc` | Resolve PDB alternate locations (`altloc` field).<br>Pick a single conformation per residue; needed before `extract` on raw RCSB downloads. |
| `bond-summary.md` | `bond-summary` | Detect bond changes between two structures (e.g. R vs P).<br>Uses the same bond-change algorithm `path-search` invokes for segmentation. |

## Pipeline at a glance

```
PDB(s) ──► extract ──► path-search ──┐
                       path-opt ─────┴──► tsopt ──► irc ──► freq ──► (dft)
                       (alternative)
```

`path-opt` and `path-search` are parallel MEP strategies (single-pass
vs recursive). `dft` is an optional terminal node and can be attached
to `freq` or directly after `irc`/`tsopt` for an ML-region DFT single-point
energy evaluation. `mlmm all` chains the whole pipeline; each box is
also available as its own subcommand.

## Common flag conventions

These flags appear on most subcommands (canonical list:
`mlmm <subcommand> --help`):

| Flag | Meaning |
|---|---|
| `-i, --input` | Input file(s); calculation workflows use `.pdb` / `.xyz` (Gaussian/ORCA input belongs to `oniom-import`) |
| `-q, --charge` | Net ML-region/model-system charge (integer) |
| `-l, --ligand-charge` | Unknown-ligand total or `'RES1:Q1,RES2:Q2'` mapping used to derive the ML-region charge |
| `-m, --multiplicity` | Spin multiplicity (2S+1), default 1 |
| `-b, --backend` | MLIP backend: `uma` / `orb` / `mace` / `aimnet2` |
| `--precision` | Unset defaults by backend: UMA/AIMNet2 fp32; ORB/MACE fp64 |
| `--workers` | UMA predictor workers; `>1` is incompatible with an analytical Hessian |
| `-o, --out-dir` | Output directory, subcommand-specific default |
| `--config` | YAML configuration file applied before CLI flags |
| `--show-config` | Print resolved merged config, then continue execution (`sp` prints and exits before evaluation) |
| `--dry-run` | Validate options and print the run plan without executing |
| `--help-advanced` | Reveal hidden / advanced flags |
| `--ref-pdb` | Reference PDB used to derive residue context for XYZ inputs |
Charge precedence: explicit `-q` > `-l 'RES:Q'` derivation > `--config` YAML > `defaults.py`.

## Canonical recipes

### Multi-input MEP for a 1-step reaction

```bash
mlmm all -i 1.R.pdb 3.P.pdb \
    -c 'SAM,GPP,MG' -l 'SAM:1,GPP:-3' \
    --tsopt --thermo \
    --out-dir result_mep
```

### Single-input scan-list (when only the reactant is available)

```bash
mlmm all -i 1.R.pdb \
    -c 'SAM,GPP,MG' -l 'SAM:1,GPP:-3' \
    --scan-lists '[("SAM 320 CS1","GPP 321 C7",1.60)]' \
                 '[("GPP`321/H11","GLU`186/OE2",0.90)]' \
    --tsopt --thermo \
    --out-dir result_scan
```

### Validate a TS candidate without re-running path search

```bash
mlmm tsopt -i ts_guess.xyz --parm real.parm7 --ref-pdb enzyme.pdb -q -1 -m 1 -b uma -o result_tsopt
mlmm freq  -i result_tsopt/final_geometry.xyz --parm real.parm7 --ref-pdb enzyme.pdb -q -1 -m 1 -b uma -o result_freq
mlmm irc   -i result_tsopt/final_geometry.xyz --parm real.parm7 --ref-pdb enzyme.pdb -q -1 -m 1 -b uma -o result_irc
```

### DFT//MLIP/MM single point on the highest-local-barrier TS candidate

```bash
mlmm dft -i seg_01/ts.pdb --parm real.parm7 \
    -l 'SAM:1,GPP:-3' \
    --func-basis 'wb97m-v/def2-tzvpd' \
    --engine gpu
```

### Bond-change report between R and P

```bash
mlmm bond-summary -i reactant.pdb product.pdb
```

## Cross-cutting caveats

| Pitfall | Fix |
|---|---|
| `--scan-lists` syntax error | The list is a Python literal-eval expression. Quote with single-quotes outside, double-quotes inside, and watch space- vs backtick-separated atom specs. |
| Wrong charge silently | Always run `--show-config` once before a long job; it prints the resolved charge. |
| Forgetting `-b` falls back to the default (`uma`) | Spell `-b uma` / `-b orb` / `-b mace` / `-b aimnet2` explicitly for production runs. |
| `--config` YAML ignored | YAML is read **after** built-in defaults but **before** explicit CLI flags. Anything also given on CLI overrides YAML. |
| `--help-advanced` flags differ between versions | They are subject to change; if a flag isn't in `--help`, check `--help-advanced` and version-pin if the workflow is shared. |
| OOM on the Hessian step | Reduce the active Hessian region and compare `Analytical` with `FiniteDifference` on the target backend/system; also keep `return_partial_hessian=True` where applicable. |
| `workers > 1` with `Analytical` | This is an intentional hard error. Use one UMA worker for an analytical Hessian, or request `FiniteDifference` before enabling the parallel predictor. |

## Defaults

Every default value is exported from `mlmm.core.defaults` (read with
`import` — the skill does not transcribe values that change between
releases):

```bash
python -c "import mlmm.core.defaults as d; print([n for n in dir(d) if n.endswith('_KW') or n.startswith('OUT_DIR')])"

# Examples:
python -c "import mlmm.core.defaults as d; print(d.LBFGS_KW)"
python -c "import mlmm.core.defaults as d; print(d.RSIRFO_KW)"
python -c "import mlmm.core.defaults as d; print(d.IRC_KW)"
python -c "import mlmm.core.defaults as d; print(d.MLMM_CALC_KW)"
```

Each per-subcommand md points at the relevant `_KW` dict in the
"See also" section.

## See also

- `mlmm-overview/SKILL.md` — what `mlmm-toolkit` is and when to
  use it.
- `mlmm-structure-io/` — input file formats and charge / spin.
- `mlmm-install-backends/` — AmberTools / backend installation.
- `mlmm-workflows-output/SKILL.md` — what comes out of each
  invocation, summary.json schema, R/TS/P canonical paths.
- `mlmm-hpc/SKILL.md` — running these recipes on PBS / SLURM.
