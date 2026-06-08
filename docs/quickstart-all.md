# Quickstart: `mlmm all`

## Goal

Run the end-to-end ML/MM ONIOM workflow once from a Reactant and Product PDB pair: MM parametrisation → ML-region selection → MEP search (path_search) → optional TS optimisation, frequencies, IRC, and DFT single-point.

## Prerequisites

- **Two complex PDBs** with matching atom order and identical residues (R/P; only coordinates differ — in PyMOL tick *Original atom order* on export).
- **Explicit hydrogens present** (`mlmm all` does not auto-protonate). Match `-l RES:CHARGE` to the H count actually in the file (e.g. SAM with 23 H = `SAM:1` cation, 22 H = `SAM:0` neutral; mismatch causes antechamber sqm odd-electron failure).
- **`mm-parm` is run automatically** by `mlmm all`. If you already have a `.parm7` / `.rst7`, pass `--parm real.parm7` (and optionally `--model-pdb ml_region.pdb`) to skip MM parametrisation.
- **GPU recommended.** Default backend is `uma` (one of `uma`, `orb`, `mace`, `aimnet2`).

## Minimal command

End-to-end MEP only (matches `examples/toy_system/run.sh` test18):

```bash
mlmm all -i r_complex.pdb p_complex.pdb -c PRE -r 6.0 \
  --ligand-charge 'PRE:0' -q -1 -m 1 --out-dir ./result_all
```

With TS optimisation, thermochemistry, and DFT single-point:

```bash
mlmm all -i r_complex.pdb p_complex.pdb -c PRE -r 6.0 \
  --ligand-charge 'PRE:0' -q -1 -m 1 \
  --tsopt --thermo --dft --out-dir ./result_all
```

(`-c` = ligand residue list for ML-region centring; `-r` = ML-region radius in Å; `-q`/`-m` = total ML charge / multiplicity.)

## Expected output tree

```
result_all/
├── summary.log                  # human-readable summary
├── summary.json                 # machine-readable: segments, ΔE‡, ΔE, bond_changes
├── mep.pdb / mep_trj.xyz        # full MEP (copied to the root)
├── energy_diagram_MEP.png       # R → TS → P energy diagram
├── mep_plot.png                 # smooth MEP energy profile
├── segments/
│   └── seg_NN/                  # per-reactive-segment deliverables (only when --tsopt/--thermo/--dft)
│       ├── reactant.pdb · ts.pdb · product.pdb   # canonical R/TS/P
│       ├── ts/final_geometry.pdb        # optimised TS
│       ├── structures/{R,TS,P}.pdb      # IRC endpoints
│       ├── freq/{R,TS,P}/               # thermochemistry per state
│       └── irc/                         # IRC trajectory & plot
└── _work/                       # pipeline scratch (safe to delete)
    └── path_search/             # raw MEP-engine output
        ├── seg_000_refine_mep/  # refined GSM/NEB path for segment 0
        └── hei_seg_01.pdb       # highest-energy image (TS approximation) per segment
```

(With `--no-refine-path` the engine directory is `_work/path_opt/` instead of `_work/path_search/`.)

## Inspecting the result

- **Per-segment barriers** — read `summary.json → segments[i].barrier_kcal` (ΔE‡, MEP electronic) and `delta_kcal` (ΔE, R→P). The `[1] Global MEP overview` and `[2] Segment-level MEP summary` blocks in `summary.log` show the same numbers and a segment-overview table.
- **Bond changes** — `segments[i].bond_changes` lists each formed / broken covalent bond as `Atom1-Atom2 : d_before Å --> d_after Å`. Use this to confirm the MEP is the chemistry you intended.
- **R / TS / P canonical paths** — when `--tsopt` is set the canonical R/TS/P are `segments/seg_NN/{reactant,ts,product}.pdb`, the optimised TS is `segments/seg_NN/ts/final_geometry.pdb`, and the IRC endpoints (chemically optimised R/P) are under `segments/seg_NN/structures/`. Without `--tsopt`, the closest TS proxy is `_work/path_search/hei_seg_NN.pdb` (highest-energy MEP image).
- **Post-processed energies** — with `--thermo` / `--dft`, `summary.log` `[3]` shows `UMA ΔG‡`, `DFT//UMA ΔE‡`, `DFT//UMA ΔG‡` per segment alongside the MEP values; the segment-overview table aggregates the same.
- **Visual** — `energy_diagram_MEP.png` is the R → TS → P stairstep figure (kcal/mol); `mep_plot.png` is the smooth MEP energy profile along the reaction coordinate.

## Tips

- `--dry-run` (shown in `--help-advanced`) validates parsing and execution plan without running heavy stages.
- `mlmm all --help` shows core options; `mlmm all --help-advanced` shows the full list.
- Switch MLIP backend with `-b orb` (or `mace`, `aimnet2`); default is `uma`.
- Add `--embedcharge` to enable xTB point-charge embedding for MM-to-ML environmental corrections.

## Troubleshoot

- **`antechamber failed (... odd-electron)` during mm-parm** — `-l RES:CHARGE` does not match the protonation state in the PDB. Count residue H atoms first and match (e.g. SAM 22 H = `:0`, 23 H = `:1`); do not blindly re-protonate.
- **`R and P atom order mismatch`** — re-export with *Original atom order* in PyMOL, or pass both files through `mlmm fix-altloc` then verify residue/atom names line up.
- **Missing `.parm7`** — `mlmm all` runs `mm-parm` automatically. If running a per-stage subcommand instead (`tsopt`, `freq`, …), prepare topology first with `mlmm mm-parm -i complex.pdb --ligand-charge '...'`.
- **GPU OOM during IRC / Hessian** — shrink the ML region (`-r 5.0`), or pass `--hessian-calc-mode FiniteDifference` (FD Hessian; slower but lower peak VRAM).
- **Link-atom non-convergence** — usually means the QM/MM boundary cuts a polar bond; reselect the ML region so the cut goes through a non-polar C–C bond (see [define-layer](define-layer.md)).
- **Memory not released between stages** — `mlmm all` already invokes `gc.collect()` and `torch.cuda.empty_cache()` between stages; if VRAM still climbs, lower `--max-cycles` or split the run into staged subcommands.

## Next step

- Single-structure scan route: [Quickstart: `mlmm scan` with `-s` (YAML spec)](quickstart-scan-spec.md)
- TS validation route: [Quickstart: `mlmm tsopt` -> `mlmm freq`](quickstart-tsopt-freq.md)
- Per-stage path control: [path-opt](path-opt.md), [path-search](path-search.md)
- Full option reference: [all](all.md)
