# ML/MM Workflow Tutorial

This tutorial walks through the full reaction–path workflow implemented in **ML/MM tools**. Commands refer to the YAML files and scripts in `examples/chorismate_mutase/` but the same procedure applies to other systems.

## 1. System preparation

1. Build an enzyme–substrate complex and generate AMBER parameters. Each complex was solvated, minimized, equilibrated and simulated for 40 ns of MD. Snapshots were extracted and trimmed to retain the protein, substrate and waters within 6 Å, then parameters and topologies were regenerated for these REAL systems.
2. The resulting snapshot provides `complex.pdb`, `complex.parm7` and `complex.rst7`. Define an ML region around the substrate and save it as `ml_region.pdb`.

## 2. ML region setup

Use the CLI tool `def_ml_region` (see `docs/code_and_cli.md`) or a molecular viewer to select residues near the substrate. Ensure the atom order matches the full `complex.pdb`.

## 3. Running the workflow

`examples/chorismate_mutase/run.sh` runs the following commands:
```bash
mlmm 1_opt.yaml
bond_scan 2_scan.yaml
mlmm 3_opt.yaml
mlmm 4_opt.yaml
mlmm 5_gs.yaml
ts_search 6_tsopt.yaml
mlmm 7_irc.yaml
mlmm 8_opt.yaml
mlmm 9_opt.yaml
energy_summary 10_energy_summary.yaml
```

Each step corresponds to the workflow described in the paper:
(i) geometry optimization of the initial structure,
(ii) bond scan to obtain candidate Reactant/Product,
(iii) optimization of those candidates,
(iv) Growing‑String MEP search,
(v) PH‑Dimer TS refinement,
(vi) IRC calculation from the TS,
(vii) optimization of the IRC endpoints,
(viii) summarizing ΔE/ΔG and ΔE‡/ΔG‡.

### 3.1 Optimization (`1_opt.yaml`)
Optimizes the MD snapshot while freezing distant atoms. Output is written under `./dump/opt1/` as `final_geometry.xyz` (copied to `coord/1_opt_final_geometry.xyz` for the next step).

### 3.2 Bond scan (`2_scan.yaml`)
Scans the forming bond and saves a trajectory in `./dump/cart_scan/`. Select representative frames from this scan as reactant and product (`2_scan_reac.pdb`, `2_scan_prod.pdb`).

### 3.3 Reactant/Product optimizations (`3_opt.yaml`, `4_opt.yaml`)
Optimize the structures selected above. Results are stored in `./dump/opt2/` and `./dump/opt3/`.

### 3.4 Growing String (`5_gs.yaml`)
Runs a GSM search between reactant and product. Afterward run `trj2fig -i gs.trj --output-peak 5_gs_peak.xyz` to extract the highest-energy geometry. Convert it to PDB with `xyz_geom2pdb 5_gs_peak.xyz --ref complex.pdb -o 5_gs_peak.pdb`. The TS search uses this geometry.

### 3.5 Transition‑state search (`6_tsopt.yaml`)
Starts from `5_gs_peak.xyz` and refines the TS with the PH‑Dimer method. The optimized TS is saved to `./dump/dimer/6_tsopt_final_geometry.xyz`.

### 3.6 IRC (`7_irc.yaml`)
Propagates the TS along the reaction coordinate in both directions. The last frames of the backward and forward paths (`7_irc_backward_last.xyz`, `7_irc_forward_last.xyz`) serve as refined reactant and product guesses.

### 3.7 Endpoint optimizations (`8_opt.yaml`, `9_opt.yaml`)
Relax the IRC endpoints to yield final reactant and product geometries stored in `./dump/opt4/` and `./dump/opt5/`.

### 3.8 Energy summary (`10_energy_summary.yaml`)
Computes electronic and vibrational contributions for the final reactant, TS and product structures. The command
```bash
energy_summary 10_energy_summary.yaml
```
produces ΔG and ΔG‡ tables and an optional energy profile.

All intermediate results are written to subdirectories under `./dump/`. Use these output files as the inputs for subsequent YAML jobs (the example `coord/` directory mirrors these files for convenience).

---