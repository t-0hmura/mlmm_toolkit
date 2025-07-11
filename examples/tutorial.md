# ML/MM Workflow Tutorial

This tutorial demonstrates **step by step** how to reproduce the
reaction‑path workflow described in our accompanying paper
(*ML/MM tools — towards Accelerated Mechanistic Investigation of
Enzymatic Reactions*). All commands refer to the files in
`examples/chorismate_mutase/`, but the same procedure applies
unchanged to any other system.
For a quick sanity check, you can also run `bash run.sh` in
`examples/toy_system/` to make sure the calculator works.

---

## 1  System preparation

1. **Build the solvated enzyme–substrate complex and generate AMBER
   parameters.**  
   * Minimize, equilibrate and run 40 ns of MD.  
   * From the trajectory, extract a snapshot and trim it to retain the
     protein, the substrate and every water molecule within 6 Å.  
   * Regenerate topology/parameter files for this trimmed *real* system.

   > If you use PYMOL, select complex and then `select byres (resn WAT within 6 of sele) or sele` to select water molecules.

2. **Collect the three required input files** (Example file name).

   * `complex.pdb`  
   * `complex.parm7`  
   * `complex.rst7`

3. **Define the ML region** around the substrate and save it as `ml_region.pdb` (Example file name.).

---

## 2  ML‑region setup

Run the CLI helper  
`def_ml_region` *(see `docs/code_and_cli.md`)* **or** use a molecular
viewer to select residues close to the substrate.  
Ensure that atom order, atom names, residue IDs and residue names are
identical to those in the full `complex.pdb`.

---

## 3  Directory layout

```text
examples/
├── chorismate_mutase/       # Example: chorismate‑mutase claisen rearrangement
│   ├── coord/               # Input + intermediate geometries (for example)
│   ├── parm/                # Topology/parameters, original complex PDB and user defined ML region (PDB)
│   │   └── freeze.txt       # Indices (>10 Å) to freeze; use as `freeze_atoms` parameter in YAML
│   ├── yaml/                # 1_opt.yaml … 10_energy_summary.yaml
│   └── run.sh               # Bash script to run the full workflow test
├── cm_mutation/             # Example: Arg90 → Ala mutant of chorismate‑mutase
└── methyltransferase/       # Example: Fur6 methyl‑transferase
```

The `cm_mutation` workflow reuses the original chorismate‑mutase setup but
mutates residue 90 from arginine to alanine. The side chain was removed in the
topology/parameter and the two PDB files used as GSM inputs, replacing the side chain except $C_{\beta}$ with hydrogens.
Freeze indices correspond to the original system, and the model charge is
decreased by one (`model_charge: -1`).  Running `run.sh` executes the steps from
`3_opt.yaml` to `10_energy_summary.yaml` to probe how this mutation affects the
reaction barrier.

---

## 4  Workflow overview

The canonical pipeline consists of **10 YAML job files**:

| Stage | YAML file(s)              | Purpose                                             | Typical wall‑time \* |
|------:|---------------------------|-----------------------------------------------------|----------------------|
| (i)   | `1_opt.yaml`              | Initial optimization of the snapshot                | 5 – 10 min |
| (ii)  | `2_scan.yaml`             | 1‑D bond‑length scan → seed Reactant/Product        | 5 – 10 min |
| (iii) | `3_opt.yaml`, `4_opt.yaml`| Relax *reactant* and *product* obtained by scan     | 1 – 2 min |
| (iv)  | `5_gs.yaml`               | Minimum Energy Path search by Growing String Method | 5 – 10 min |
| (v)   | `6_tsopt.yaml`            | **PH‑Dimer** transition‑state refinement            | **30 – 60 min** |
| (vi)  | `7_irc.yaml`              | IRC propagation                                     | 5 – 10 min |
| (vii) | `8_opt.yaml`, `9_opt.yaml`| Endpoint relaxation → Final Reactant/Product        | 2 – 3 min |
| (viii)| `10_energy_summary.yaml`  | $\Delta E$ and $\Delta G$ tables + optional energy plots         | 10 – 30 min |

\* Ryzen 7950X3D + RTX 5080; see Section 5 for full benchmarks.

The bash script `examples/chorismate_mutase/run.sh` executes the
stages in order:

```bash
mlmm          1_opt.yaml          # (i)
bond_scan     2_scan.yaml         # (ii)
mlmm          3_opt.yaml          # (iii‑R)
mlmm          4_opt.yaml          # (iii‑P)
mlmm          5_gs.yaml           # (iv)
ts_search     6_tsopt.yaml        # (v)
mlmm          7_irc.yaml          # (vi)
mlmm          8_opt.yaml          # (vii‑R)
mlmm          9_opt.yaml          # (vii‑P)
energy_summary 10_energy_summary.yaml   # (viii)
```

---

### 4.1  optimization – `1_opt.yaml`

optimizes the MD snapshot while freezing atoms farther than 10 Å from
the ML region.  
The final structure is written to
`./dump/opt1/final_geometry.xyz` and copied to
`coord/1_opt_final_geometry.xyz` for the next step.

### 4.2  Bond scan – `2_scan.yaml`

Performs a 1‑D scan along the forming bond.  
The trajectory is saved in `./dump/cart_scan/final_geometries.trj`.  
Afterwards:

```bash
xyz_geom2pdb -i final_geometries.trj -o scan_path.pdb -r ../../parm/complex.pdb
```

Pick two representative frames as the *reactant* and *product* from `scan_path.pdb` :
`2_scan_reac.pdb`, `2_scan_prod.pdb`.

### 4.3  Reactant / product optimizations – `3_opt.yaml`, `4_opt.yaml`

Refine the two structures chosen above.  
Results are stored as `./dump/opt2/final_geometry.xyz` and `./dump/opt3/final_geometry.xyz`.

### 4.4  Growing‑String search – `5_gs.yaml`

Runs a GSM search between the optimized reactant and product.  
Afterwards:

```bash
trj2fig -i current_geometries.trj -o gs.png --output-peak 5_gs_peak.xyz
```

`5_gs_peak.xyz` serves as the initial guess for the TS search.

### 4.5  Transition‑state search – `6_tsopt.yaml`

Starts from `5_gs_peak.xyz` and refines the TS with the
**PH‑Dimer** algorithm.  
The optimized TS is saved as `./dump/dimer/final_geometry.xyz`.

### 4.6  IRC – `7_irc.yaml`

Propagates the TS downhill in both directions.  
The final frames (`./dump/irc/backward_last.xyz` → `7_irc_backward_last.xyz`, 
`./dump/irc/forward_last.xyz` → `7_irc_forward_last.xyz`) provide improved reactant/product guesses.

### 4.7  Endpoint optimizations – `8_opt.yaml`, `9_opt.yaml`

Relax the IRC endpoints to obtain the *final* reactant and product
geometries (stored as `./dump/opt4/final_geometry.xyz` and `./dump/opt5/final_geometry.xyz`).

### 4.8  Energy summary – `10_energy_summary.yaml`

Compute electronic + vibrational contributions for the final
reactant, TS and product:

```bash
energy_summary 10_energy_summary.yaml
```

This command prints $\Delta G$, $\Delta G^{\ddagger}$, $\Delta E$ and $\Delta E^{\ddagger}$ and write figures of energy diagrams.

> All intermediate files live under `./dump/`.
> In this example, all input geometries are stored in `./coord/`.
---

## 5  Benchmarks  
*(UMA / ff14SB / GAFF2 / TIP3P)*

| Stage | Chorismate‑mutase | Fur6 |
|------:|------------------:|-----:|
| Geometry optimization      | 0 h 06 m 45 s | 0 h 08 m 23 s |
| Bond scan                  | 0 h 06 m 54 s | 0 h 05 m 46 s |
| R/P optimization           | 0 h 01 m 32 s | 0 h 02 m 18 s |
| GSM                        | 0 h 05 m 14 s | 0 h 06 m 01 s |
| **PH‑Dimer TS**            | **0 h 33 m 50 s** | **0 h 49 m 00 s** |
| IRC                        | 0 h 05 m 01 s | 0 h 07 m 31 s |
| Endpoint optimization      | 0 h 01 m 56 s | 0 h 02 m 27 s |
| Thermochemistry            | 0 h 13 m 35 s | 0 h 20 m 50 s |
| **Total**                  | **1 h 14 m 47 s** | **1 h 42 m 16 s** |

---

## 6  Energetic results

| Reaction | $\Delta E^{\ddagger}$ (kcal mol<sup>-1</sup>) | $\Delta E$ | $\Delta G^{\ddagger}$ | $\Delta G$ | $\Delta E^{\ddagger}$<sub>QM/MM</sub> | $\Delta G^{\ddagger}$<sub>Exp.</sub> |
|----------|-----------------:|----:|----:|----:|---------------------:|-------------------:|
| CM — Claisen rearrangement | 18.1 | −23.4 | 14.7 | −23.3 | 16.1 | 15.4 |
| CM — Arg90Ala Mutant       | 30.8 | −18.0 | 28.9 | −17.6 | — | — |
| Fur6 — methyl transfer     | 10.2 | −62.5 |  9.6 | −58.3 |  8.6 | — |

From these data it is evident that both the CM Claisen rearrangement and the Fur6 methyl-transfer reaction possess activation barriers low enough for the reactions to proceed at room temperature, with computed values lying within a few kcal mol<sup>-1</sup> of the available experimental and QM/MM results. Moreover, in the CM Arg90Ala mutant—where Arg90 is thought to contribute catalytically—the activation barrier rises sharply, consistent with the experimentally observed loss of activity (Kast P. et al., 2000, J. Biol. Chem. 275 (47), 36832–36838).

---

\* Timings were obtained on an AMD Ryzen 7950X3D (32 threads, 4.2 GHz)
workstation equipped with an NVIDIA RTX 5080 (16 GB VRAM).
