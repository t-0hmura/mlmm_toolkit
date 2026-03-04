# ML/MM Workflow Tutorial

This tutorial demonstrates **step by step** how to reproduce the
reaction-path workflow described in our accompanying paper
(*ML/MM toolkit — Toward Accelerated Mechanistic Investigation of
Enzymatic Reactions*).
For a quick sanity check, run `bash run.sh` in `examples/toy_system/`
to make sure the calculator works.

---

## 1  System preparation

1. **Build the solvated enzyme-substrate complex and generate AMBER
   parameters.**
   * Minimize, equilibrate and run 40 ns of MD.
   * From the trajectory, extract a snapshot and trim it to retain the
     protein, the substrate and every water molecule within 6 A.
   * Regenerate topology/parameter files for this trimmed *real* system.

   > With PyMOL: `select byres (resn WAT within 6 of sele) or sele`

2. **Collect the three required input files:**
   * `complex.pdb`
   * `complex.parm7`
   * `complex.rst7`

3. **Define the ML region** around the substrate and save it as `ml_region.pdb`.

---

## 2  ML-region setup

Use `mlmm define-layer` or a molecular viewer to select residues close
to the substrate.  Ensure that atom order, atom names, residue IDs and
residue names are identical to those in the full `complex.pdb`.

```bash
mlmm define-layer \
  -i parm/complex.pdb \
  --model-pdb parm/ml_region.pdb \
  --radius-freeze 10.0 \
  -o complex_layered.pdb
```

---

## 3  Directory layout

```text
examples/
├── chorismate_mutase/       # Chorismate-mutase Claisen rearrangement
│   ├── coord/               # Input + intermediate geometries
│   ├── parm/                # Topology, full PDB, ML region PDB
│   └── run.sh               # CLI workflow script
├── cm_mutation/             # Arg90 → Ala mutant (steps iii-viii)
├── methyltransferase/       # Fur6 methyl-transferase
├── toy_system/              # Python API examples + CLI smoke test
│   ├── opt_ase.py           # ASE interface
│   ├── opt_pysisyphus.py    # PySisyphus interface
│   ├── test_core.py         # Direct MLMMCore API
│   └── run.sh
└── tutorial.md              # This file
```

---

## 4  Workflow overview

The canonical pipeline consists of **8 stages**, each mapping to an
`mlmm` CLI subcommand:

| Stage | CLI command | Purpose | Typical wall-time \* |
|------:|-------------|---------|----------------------|
| (0)   | `mlmm define-layer` | Assign ML/MM layers | < 1 min |
| (i)   | `mlmm opt`     | Initial optimization of the snapshot | 5 - 10 min |
| (ii)  | `mlmm scan`    | 1-D bond-length scan → seed Reactant/Product | 5 - 10 min |
| (iii) | `mlmm opt`     | Relax *reactant* and *product* | 1 - 2 min |
| (iv)  | `mlmm path-opt`| GSM path search between R and P | 5 - 10 min |
| (v)   | `mlmm tsopt`   | Transition state refinement (Dimer) | **30 - 60 min** |
| (vi)  | `mlmm irc`     | IRC propagation | 5 - 10 min |
| (vii) | `mlmm opt`     | Endpoint relaxation → Final R/P | 2 - 3 min |
| (viii) | `mlmm freq`   | Vibrational analysis & thermochemistry | 10 - 30 min |

\* AMD Ryzen 7950X3D + NVIDIA RTX 5080.

The bash script `examples/chorismate_mutase/run.sh` executes the
stages in order:

```bash
mlmm define-layer  ...                         # (0)
mlmm opt           -i complex_layered.pdb ...  # (i)
mlmm scan          -i coord/1_opt_... ...      # (ii)
mlmm opt           -i coord/2_scan_reac.pdb ...  # (iii-R)
mlmm opt           -i coord/2_scan_prod.pdb ...  # (iii-P)
mlmm path-opt      -i coord/3_opt_... coord/4_opt_... ...  # (iv)
mlmm tsopt         -i coord/5_gs_peak.xyz ...  # (v)
mlmm irc           -i coord/6_tsopt_... ...    # (vi)
mlmm opt           -i coord/7_irc_forward_... ... # (vii-R)
mlmm opt           -i coord/7_irc_backward_... ... # (vii-P)
mlmm freq          -i coord/6_tsopt_... ...    # (viii)
```

---

### 4.1  Initial optimization — `mlmm opt`

Optimizes the MD snapshot using the LBFGS method while freezing atoms
in the outer layer (> 10 A from the ML region).

### 4.2  Bond scan — `mlmm scan`

Performs a 1-D scan along the forming bond.  Pick two representative
frames as the *reactant* and *product* seeds.

### 4.3  Reactant / product optimizations — `mlmm opt`

Refine the two structures obtained from the scan.

### 4.4  Growing String search — `mlmm path-opt`

Runs a GSM search between the optimized reactant and product.
The highest-energy image (HEI) serves as the initial TS guess.

### 4.5  Transition state search — `mlmm tsopt`

Refines the TS with the **Dimer** algorithm (`--opt-mode grad`).

### 4.6  IRC — `mlmm irc`

Propagates the TS downhill in both directions.  The endpoints provide
improved reactant/product geometries.

### 4.7  Endpoint optimizations — `mlmm opt`

Relax the IRC endpoints to obtain the *final* reactant and product.

### 4.8  Vibrational analysis — `mlmm freq`

Compute vibrational frequencies and thermochemical quantities
(ZPE, Delta G) for the TS, reactant and product.

---

## 5  Alternative: end-to-end with `mlmm all`

For routine use, `mlmm all` combines extraction, path search, TS
optimization, IRC, frequency analysis and DFT into a single command:

```bash
mlmm all \
  -i r_complex.pdb p_complex.pdb \
  -c PRE \
  -r 6.0 \
  --ligand-charge "PRE:0" \
  -q 0 -m 1 \
  --tsopt True \
  --thermo True \
  --dft True \
  --out-dir results
```

---

## 6  Python API

The calculator is also accessible from Python (see `examples/toy_system/`):

```python
# ASE interface
from mlmm_toolkit.mlmm_calc import MLMMASECalculator
atoms.calc = MLMMASECalculator(
    input_pdb="complex.pdb", real_parm7="complex.parm7",
    model_pdb="ml_region.pdb", model_charge=-1, model_mult=1,
)

# PySisyphus interface
from mlmm_toolkit.mlmm_calc import mlmm
geom.set_calculator(mlmm(
    input_pdb="complex.pdb", real_parm7="complex.parm7",
    model_pdb="ml_region.pdb", model_charge=-1, model_mult=1,
))

# Direct MLMMCore API
from mlmm_toolkit.mlmm_calc import MLMMCore
core = MLMMCore(
    input_pdb="complex.pdb", real_parm7="complex.parm7",
    model_pdb="ml_region.pdb", model_charge=-1, model_mult=1,
)
results = core.compute(coord_ang, return_forces=True, return_hessian=True)
```

---

## 7  Benchmarks
*(UMA / ff14SB / GAFF2 / TIP3P)*

| Stage | Chorismate-mutase | Fur6 |
|------:|------------------:|-----:|
| Geometry optimization      | 0 h 06 m 45 s | 0 h 08 m 23 s |
| Bond scan                  | 0 h 06 m 54 s | 0 h 05 m 46 s |
| R/P optimization           | 0 h 01 m 32 s | 0 h 02 m 18 s |
| GSM                        | 0 h 05 m 14 s | 0 h 06 m 01 s |
| **Dimer TS**               | **0 h 33 m 50 s** | **0 h 49 m 00 s** |
| IRC                        | 0 h 05 m 01 s | 0 h 07 m 31 s |
| Endpoint optimization      | 0 h 01 m 56 s | 0 h 02 m 27 s |
| Thermochemistry            | 0 h 13 m 35 s | 0 h 20 m 50 s |
| **Total**                  | **1 h 14 m 47 s** | **1 h 42 m 16 s** |

---

## 8  Energetic results

| Reaction | Delta E_act (kcal/mol) | Delta E | Delta G_act | Delta G | Delta E_act (QM/MM) | Delta G_act (Exp.) |
|----------|---:|---:|---:|---:|---:|---:|
| CM — Claisen rearrangement | 18.1 | -23.4 | 14.7 | -23.3 | 16.1 | 15.4 |
| CM — Arg90Ala Mutant       | 30.8 | -18.0 | 28.9 | -17.6 | — | — |
| Fur6 — methyl transfer     |  9.0 | -64.9 |  7.9 | -60.0 | 8.6 | — |

---

\* Timings were obtained on an AMD Ryzen 7950X3D (32 threads, 4.2 GHz)
workstation equipped with an NVIDIA RTX 5080 (16 GB VRAM).
