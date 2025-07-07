# ML/MM Workflow **Master Tutorial**  _(updated 2025-07-07)_

This document distills the **strongest points** of the three prior drafts  
(“Revised 2025” fileciteturn0file0, “Detailed v2025.4” fileciteturn0file2, and the “Quick‑Start” sheet fileciteturn0file4) into a single, self‑contained guide that serves:

* **Beginners** – get up‑and‑running in minutes;  
* **Power‑users** – learn every optimisation switch and recovery trick;  
* **Maintainers** – centralise parameters & paths so future updates touch one place.

---

## 0 . TL;DR Quick‑Start  *(five commands)*

```bash
# 1.  Clone examples   (or use your own complex)
git clone https://github.com/t-0hmura/mlmm_tools_examples.git
cd mlmm_tools_examples/chorismate_mutase

# 2.  Make/activate conda env (optional)
conda create -n mlmm python=3.11 -y && conda activate mlmm

# 3.  Install ML/MM tools (+ UMA backend)
pip install fairchem-core==2.3.0
pip install git+https://github.com/t-0hmura/mlmm_tools.git
huggingface-cli login               # UMA model access

# 4.  Pull pretrained weights (first run only)
mlmm pull-weights --backend uma

# 5.  Launch the 8‑stage workflow
bash run.sh                         # ≈ 75 min on 7950X3D + RTX 5080
```

Resulting tables and figures appear in `tables/` & `fig/`.

---

## 1 . Directory Layout (opinionated)

```
tutorial/
├── cm/                      # Chorismate Mutase example
│   ├── coord/               # Input & intermediate geometries
│   ├── yaml/                # 1_opt.yaml … 10_energy_summary.yaml
│   ├── freeze.txt           # Atoms to freeze (>10 Å)
│   ├── run.sh               # Convenience launcher
│   └── README.md            # System‑specific notes
├── fur6/                    # Fur6 methyl‑transferase example
└── common/                  # Helper scripts (trj2fig, xyz_geom2pdb …)
```

Change **exactly one** variable – `EXAMPLE=cm|fur6` – in `run.sh` to swap systems.

---

## 2 . System Preparation — 3 Key Steps

| # | Action | One‑liner |
|---|--------|-----------|
| **2‑1** | **Trim MD snapshot** (keep protein + substrate + waters ≤ 6 Å) | `cpptraj … strip !(:SUB<6.0)` |
| **2‑2** | **Define ML region** (residues ≤ 6 Å of substrate) | `def_ml_region -r REAL.pdb -c ligand.pdb -o ML.pdb --radius_ml 6` |
| **2‑3** | **Freeze outer shell** (atoms > 10 Å) | `get_freeze_indices REAL.pdb ML.pdb --cutoff 10 > freeze.txt` |

Tips for tricky systems are in § 9.

---

## 3 . Workflow Overview

Eight YAMLs correspond to the canonical mechanistic pipeline:

| Stage | YAML | Purpose | Typical time* |
|-------|------|---------|---------------|
| (i) | `1_opt.yaml` | Freeze‑shell minimisation of snapshot | 6–8 min |
| (ii)| `2_scan.yaml`| 1‑D bond scan → seed R/P | 6–7 min |
| (iii)| `3_opt.yaml` + `4_opt.yaml` | Relax Reactant / Product | 1–2 min |
| (iv)| `5_gs.yaml`  | Growing‑String MEP search | 5–6 min |
| (v) | `6_tsopt.yaml`| **PH‑Dimer** TS refinement | **34–49 min** |
| (vi)| `7_irc.yaml` | IRC propagation | 5–7 min |
| (vii)| `8_opt.yaml` + `9_opt.yaml` | Endpoint relax | 2–3 min |
| (viii)| `10_energy_summary.yaml` | ΔE / ΔG tables & plots | 14–21 min |

\* Timings on 7950X3D + RTX 5080 (see § 6 for full table).

Run everything:

```bash
bash run.sh            # or mlmm yaml/1_opt.yaml   (single stage)
```

Parallelism: `mm_threads` = physical cores (16 on 7950X3D). GPU is used **only** by the ML potential; MM stays on CPU.

---

## 4 . YAML Template Snippet

```yaml
# 5_gs.yaml  —  GSM search
mlmm:
  backend: uma                     # "uma" | "aimnet2"
  ml_device: auto                  # auto picks first CUDA GPU
  ml_cuda_idx: 0
  mm_threads: 16
  freeze_file: freeze.txt
gsm:
  n_nodes: 9
  max_iter: 100
  conv_E: 1e-4
  conv_F: 5e-3
```

Change `backend: aimnet2` to switch potentials (§ 8).

---

## 5 . Inspecting and Restarting

* **All raw trajectories** → `dump/`  
* **Geometry hand‑off files** → `coord/` (mirrors key XYZ/PDB so a stage can be rerun in isolation)  
* **Plots** (`*_profile.png`) & **ΔG tables** (`ΔG_summary.csv`) → `fig/` & `tables/`

To resume a failed PH‑Dimer (most common hic‑cup):

```bash
mlmm 6_tsopt.yaml --restart dump/dimer/last_state.chk
```

---

## 6 . Benchmarks (UMA / ff14SB / GAFF2 / TIP3P)

| Stage | CM | Fur6 |
|-------|-------|-------|
| Geometry opt. | 0 h 06 m 45 s | 0 h 08 m 23 s |
| Bond scan | 0 h 06 m 54 s | 0 h 05 m 46 s |
| R/P opt. | 0 h 01 m 32 s | 0 h 02 m 18 s |
| GSM | 0 h 05 m 14 s | 0 h 06 m 01 s |
| **PH‑Dimer TS** | **0 h 33 m 50 s** | **0 h 49 m 00 s** |
| IRC | 0 h 05 m 01 s | 0 h 07 m 31 s |
| Endpoint opt. | 0 h 01 m 56 s | 0 h 02 m 27 s |
| Thermochemistry | 0 h 13 m 35 s | 0 h 20 m 50 s |
| **Total** | **1 h 14 m 47 s** | **1 h 42 m 16 s** |

---

## 7 . Energetic Results

| Reaction | ΔE‡ | ΔE | ΔG‡ | ΔG | ΔE‡<sub>QM/MM</sub> | ΔG‡<sub>Exp.</sub> |
|----------|----:|---:|----:|---:|---------------------:|-------------------:|
| CM Claisen rearr. | 18.1 | −23.4 | 14.7 | −23.3 | 16.1 | 15.4 |
| Fur6 methyl transfer | 10.2 | −62.5 | 9.6 | −58.3 | 8.6 | — |

Agreement within **≤ 2 kcal mol⁻¹** validates the ML/MM approach.

---

## 8 . Switching to AIMNet2

1. `sed -i 's/backend: uma/backend: aimnet2/' yaml/*.yaml`  
2. Install backend:  
   ```bash
   pip install git+https://github.com/isayevlab/aimnetcentral.git
   ```  
3. Re‑run `bash run.sh` — barriers typically shift by ≤ 1 kcal mol⁻¹.

---

## 9 . Troubleshooting (top 4 issues)

| Symptom | Likely cause | Fix |
|---------|--------------|-----|
| **CUDA OOM** during PH‑Dimer | VRAM peak (numerical Hessian) | `ml_device=cpu` for TS step _only_ or reduce `ml_batch_size` |
| TS search oscillates | GSM peak too far from saddle | Decrease `n_nodes` or supply manual TS guess |
| `ValueError: atoms out of sync` | PDB atom order mismatch | Regenerate `ML.pdb` **after** final trimming |
| “Fake” imaginary modes (< 50 cm⁻¹) | Incomplete Hessian update | Increase `n_hess_update`; enable `mass_scale_flat` |

For deeper profiling see `docs/performance.md`.

---

## 10 . Next Steps & Extensibility

* **Large complexes (> 20 k atoms)** – shrink freeze radius to 8 Å and run PH‑Dimer with partial Hessians (see `examples/large_system/`).  
* **Periodic boundary support** – prototype branch `feature/pbc_md`.  
* **GPU MM layer** – OpenMM GPU build is recognised if `mm_device:auto` and CUDA available.

---

## 11 . References

Chook _et al._ 1994; Kast _et al._ 1996; Ranaghan & Mulholland 2004; Zhao _et al._ 2024.

---

> **Maintainer note:** edit hardware specs, freeze radius, and table values in the **header of `run.sh`**; all linked YAMLs inherit those via environment variables, centralising future updates.

*Happy modelling!*
