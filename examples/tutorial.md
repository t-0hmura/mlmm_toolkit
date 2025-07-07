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
   * Minimise, equilibrate and run 40 ns of MD.  
   * From the trajectory, extract a snapshot and trim it to retain the
     protein, the substrate and every water molecule within 6 Å.  
   * Regenerate topology/parameter files for this trimmed *real* system.

   > If you use PYMOL, select complex and then `select byres (resn WAT within 6 of sele) or sele` to select water molecules.

2. **Collect the three required input files.**

   * `complex.pdb`  
   * `complex.parm7`  
   * `complex.rst7`

3. **Define the ML region** around the substrate and save it as
   `ml_region.pdb` (see Section&nbsp;2).

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
├── chorismate_mutase/       # Example: chorismate‑mutase Claisen shift
│   ├── coord/               # Input + intermediate geometries
│   ├── parm/                # Topology/parameters, original complex PDB,
│   │   └── freeze.txt       # Indices (>10 Å) to freeze; paste into YAML
│   ├── yaml/                # 1_opt.yaml … 10_energy_summary.yaml
│   └── run.sh               # Bash script to run the full workflow
└── methyltransferase/       # Example: Fur6 methyl‑transferase
```

---

## 4  Workflow overview

The canonical pipeline consists of **ten YAML job files**:

| Stage | YAML file(s)              | Purpose                                        | Typical wall‑time \* |
|------:|---------------------------|------------------------------------------------|----------------------|
| (i)   | `1_opt.yaml`              | Freeze‑shell minimisation of the snapshot      | 5 – 10 min |
| (ii)  | `2_scan.yaml`             | 1‑D bond‑length scan → seed R/P                | 5 – 10 min |
| (iii) | `3_opt.yaml`, `4_opt.yaml`| Relax *reactant* and *product*                 | 1 – 2 min |
| (iv)  | `5_gs.yaml`               | Growing‑String MEP search                      | 5 – 10 min |
| (v)   | `6_tsopt.yaml`            | **PH‑Dimer** transition‑state refinement       | **30 – 60 min** |
| (vi)  | `7_irc.yaml`              | IRC propagation                                | 5 – 10 min |
| (vii) | `8_opt.yaml`, `9_opt.yaml`| Endpoint relaxation                            | 2 – 3 min |
| (viii)| `10_energy_summary.yaml`  | $\Delta E$ / $\Delta G$ tables + optional energy plots         | 10 – 30 min |

\* Ryzen 7950X3D + RTX 5080; see Section 5 for full benchmarks.

The helper script `examples/chorismate_mutase/run.sh` executes the
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

### 4.1  Optimisation – `1_opt.yaml`

Optimises the MD snapshot while freezing atoms farther than 10 Å from
the ML region.  
The final structure is written to
`./dump/opt1/final_geometry.xyz` and copied to
`coord/1_opt_final_geometry.xyz` for the next step.

### 4.2  Bond scan – `2_scan.yaml`

Performs a 1‑D scan along the forming bond.  
The trajectory is saved in `./dump/cart_scan/`.  
Pick two representative frames as the *reactant* and *product*:
`2_scan_reac.pdb`, `2_scan_prod.pdb`.

### 4.3  Reactant / product optimisations – `3_opt.yaml`, `4_opt.yaml`

Refine the two structures chosen above.  
Results are stored in `./dump/opt2/` and `./dump/opt3/`.

### 4.4  Growing‑String search – `5_gs.yaml`

Runs a GSM search between the optimised reactant and product.  
Afterwards:

```bash
trj2fig -i gs.trj --output-peak 5_gs_peak.xyz
```

`5_gs_peak.xyz` serves as the initial guess for the TS search.

### 4.5  Transition‑state search – `6_tsopt.yaml`

Starts from `5_gs_peak.xyz` and refines the TS with the
**PH‑Dimer** algorithm.  
The optimised TS is saved as
`./dump/dimer/6_tsopt_final_geometry.xyz`.

### 4.6  IRC – `7_irc.yaml`

Propagates the TS downhill in both directions.  
The final frames (`7_irc_backward_last.xyz`,
`7_irc_forward_last.xyz`) provide improved reactant/product guesses.

### 4.7  Endpoint optimisations – `8_opt.yaml`, `9_opt.yaml`

Relax the IRC endpoints to obtain the *final* reactant and product
geometries (stored in `./dump/opt4/` and `./dump/opt5/`).

### 4.8  Energy summary – `10_energy_summary.yaml`

Compute electronic + vibrational contributions for the final
reactant, TS and product:

```bash
energy_summary 10_energy_summary.yaml
```

This command prints $\Delta G$, $\Delta G^{\ddagger}$, $\Delta E$ and $\Delta E^{\ddagger}$ and can optionally write an
interactive energy diagram.

All intermediate files live under `./dump/`; the `coord/` directory
mirrors the key geometries for quick inspection.

---

## 5  Benchmarks  
*(UMA / ff14SB / GAFF2 / TIP3P)*

| Stage | Chorismate‑mutase | Fur6 |
|------:|------------------:|-----:|
| Geometry optimisation      | 0 h 06 m 45 s | 0 h 08 m 23 s |
| Bond scan                  | 0 h 06 m 54 s | 0 h 05 m 46 s |
| R/P optimisation           | 0 h 01 m 32 s | 0 h 02 m 18 s |
| GSM                        | 0 h 05 m 14 s | 0 h 06 m 01 s |
| **PH‑Dimer TS**            | **0 h 33 m 50 s** | **0 h 49 m 00 s** |
| IRC                        | 0 h 05 m 01 s | 0 h 07 m 31 s |
| Endpoint optimisation      | 0 h 01 m 56 s | 0 h 02 m 27 s |
| Thermochemistry            | 0 h 13 m 35 s | 0 h 20 m 50 s |
| **Total**                  | **1 h 14 m 47 s** | **1 h 42 m 16 s** |

---

## 6  Energetic results

| Reaction | $\Delta E^{\ddagger}$ (kcal mol⁻¹) | $\Delta E$ | $\Delta G^{\ddagger}$ | $\Delta G$ | $\Delta E^{\ddagger}$<sub>QM/MM</sub> | $\Delta G^{\ddagger}$<sub>Exp.</sub> |
|----------|-----------------:|----:|----:|----:|---------------------:|-------------------:|
| CM — Claisen rearrangement | 18.1 | −23.4 | 14.7 | −23.3 | 16.1 | 15.4 |
| Fur6 — methyl transfer     | 10.2 | −62.5 |  9.6 | −58.3 |  8.6 | — |

---

\* Timings were obtained on an AMD Ryzen 7950X3D (32 threads, 4.2 GHz)
workstation equipped with an NVIDIA RTX 5080 (16 GB VRAM).
