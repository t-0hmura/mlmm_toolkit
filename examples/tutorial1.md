# ML/MM Workflow Tutorial (Definitive Edition, July 2025)

This **definitive guide** merges the strongest points of the three previous drafts and the reviewers’ suggestions.  
It walks you through installing the software, preparing an enzyme–substrate system, running the automated reaction‑path workflow, and validating the results.

> **Benchmark systems used throughout this guide**  
> 1. **Chorismate Mutase (CM)** – Claisen rearrangement (chorismate → prephenate)  
> 2. **Fur6 methyl‑transferase** – SAM‑dependent CH₃ transfer

---

## 0  Prerequisites & Quick Install

```bash
# CUDA 12.6   (see README for CUDA 12.8 / RTX 50 series)
pip install fairchem-core==2.3.0
pip install git+https://github.com/isayevlab/aimnetcentral.git       # AIMNet2 backend
pip install git+https://github.com/t-0hmura/mlmm_tools.git           # ML/MM tools
huggingface-cli login                                                # pull UMA weights
```

* **GPU recommended.** The tutorial was benchmarked on a workstation with  
  _AMD Ryzen 7950X3D (32 threads), 128 GB RAM, NVIDIA RTX 5080 16 GB_.

---

## 1  Directory Layout

```
tutorial/
├── cm/                    # Chorismate Mutase example
│   ├── coord/             # PDB / XYZ files created on‑the‑fly
│   ├── yaml/              # Workflow YAML inputs (steps 1–10)
│   ├── run.sh             # One‑shot driver script
│   └── 00_prep.sh         # Optional structure‑prep helper
├── fur6/                  # Fur6 methyltransferase example
└── common/                # Re‑usable helper utilities
```

Clone `examples/` from the ML/MM tools repo into `tutorial/`, then work inside either `cm/` or `fur6/`.

---

## 2  System Preparation

1. **Obtain or build the complex**  
   * **CM:** PDB ID `1COM` (trimer)  
   * **Fur6:** coordinates from Zhao _et al._ 2024  
2. **Parameterise with AMBER** (ff14SB + GAFF2 + TIP3P)  
   Use `tleap` or `ParmEd` to create `complex.parm7 / complex.rst7`.
3. **40 ns MD** (minimisation → equilibration → production) with Amber24.  
4. **Trim the REAL system** – keep protein, substrate & waters **≤ 6 Å** from the substrate:  

   ```bash
   cpptraj << 'EOF'
   parm full_system.parm7
   trajin full_system.nc 2000
   strip !(:1-XYZ<6.0)
   trajout REAL.pdb pdb
   trajout REAL.parm7 parm7
   trajout REAL.rst7 restart
   EOF
   ```
5. **Define the ML region** (side chains within 6 Å of the substrate):

   ```bash
   def_ml_region -r REAL.pdb -c substrate_only.pdb \
                 --radius_ml 6.0 --exclude_backbone true \
                 --include_H2O false -o ML_region.pdb
   ```

---

## 3  Freezing Outer‑Shell Atoms

Speed up local optimisations by freezing atoms that are **> 10 Å** from the ML region:

```bash
get_freeze_indices REAL.pdb ML_region.pdb --cutoff 10.0 > freeze.txt
```

Reference this `freeze.txt` in every YAML file.

---

## 4  Automated Reaction‑Path Workflow

Execute **all** ten stages in one go:

```bash
bash run.sh        # inside cm/ or fur6/
```

| Step | YAML / CLI                   | Purpose |
|------|------------------------------|---------|
| (i)  | `mlmm 1_opt.yaml`            | Optimise MD snapshot (frozen shell) |
| (ii) | `bond_scan 2_scan.yaml`      | 1‑D scan → pick Reactant/Product |
| (iii)| `mlmm 3_opt.yaml` & `4_opt.yaml` | Relax candidates |
| (iv) | `mlmm 5_gs.yaml`             | Growing‑String MEP search |
| (v)  | `ts_search 6_tsopt.yaml`     | PH‑Dimer TS refinement |
| (vi) | `mlmm 7_irc.yaml`            | IRC from TS |
| (vii)| `mlmm 8_opt.yaml` & `9_opt.yaml` | Final Reactant/Product |
| (viii)| `energy_summary 10_energy_summary.yaml` | ΔE / ΔG tables & plots |

---

## 5  Validation Results

### 5.1 Wall‑Clock Timing

| Stage | CM | Fur6 |
|-------|-----------|-----------|
| Initial optimisation | 0 h 06 m 45 s | 0 h 08 m 23 s |
| Reaction scan        | 0 h 06 m 54 s | 0 h 05 m 46 s |
| Candidate opt.       | 0 h 01 m 32 s | 0 h 02 m 18 s |
| GSM (MEP) search     | 0 h 05 m 14 s | 0 h 06 m 01 s |
| **PH‑Dimer TS**      | **0 h 33 m 50 s** | **0 h 49 m 00 s** |
| IRC                  | 0 h 05 m 01 s | 0 h 07 m 31 s |
| Endpoint opt.        | 0 h 01 m 56 s | 0 h 02 m 27 s |
| Thermochemistry      | 0 h 13 m 35 s | 0 h 20 m 50 s |
| **Total**            | **1 h 14 m 47 s** | **1 h 42 m 16 s** |

### 5.2 Energetic Results (UMA/ff14SB/GAFF2/TIP3P)

| Reaction | ΔE‡ | ΔE | ΔG‡ | ΔG | QM/MM ΔE‡ | Exp. ΔG‡ |
|----------|----:|---:|----:|---:|----------:|---------:|
| CM Claisen rearr. | 18.1 | −23.4 | 14.7 | −23.3 | 16.1 | 15.4 |
| Fur6 methyl transfer | 10.2 | −62.5 |  9.6 | −58.3 |  8.6 | — |

The ML/MM protocol reproduces reference QM/MM and experimental barriers within **≈ 2 kcal mol⁻¹**.

---

## 6  Inspecting the Reaction Path

* **`dump/gs/`** – GSM trajectory  
* **`dump/dimer/`** – Optimised TS (`*_tsopt_final_geometry.xyz`)  
* **`dump/irc/`** – Forward / backward IRC paths  
* **`tables/ΔG_summary.csv`** – Energies  
* **`fig/energy_profile.png`** – Potential‑ & free‑energy plots  

Load `.xyz` files in **VMD** or **PyMOL** to visualise structural changes.

---

## 7  Switching Back‑Ends: UMA ↔ AIMNet2

Set **one line** in each YAML:

```yaml
backend: aimnet2          # default is "uma"
aimnet2_model: aimnet2-s1x   # optional; defaults to s1x
```

**CLI override** (single‑stage quick test):

```bash
mlmm 1_opt.yaml --backend aimnet2 --aimnet2_model aimnet2-s1x
```

In our tests, AIMNet2 changes ΔG‡ by ≤ 1 kcal mol⁻¹ while running 20‑30 % faster on RTX 5080 GPUs.

---

## 8  Troubleshooting at a Glance

| Symptom | Likely cause | Remedy |
|---------|--------------|--------|
| `CUDA out of memory` during TS search | VRAM exhausted | Set `ml_device=cpu` for (v) & (vi) *only*, or lower `ml_batch_size`. |
| Dimer fails to converge | TS guess too far from saddle | Re‑run GSM with finer `n_images` or increase `dimer_max_cycles`. |
| `atoms out of sync` error | Atom order mismatch between REAL & ML region | Regenerate `ML_region.pdb` **after** final trimming. |

---

## 9  Next Steps & Support

* **Periodic MD sampling** & **PBC** support are in the dev branch – follow the GitHub repo for updates.  
* For systems > 20 000 atoms tighten the freeze radius to 8 Å and use a partial‑Hessian Dimer (`phdimer_partial: true`).  
* Join `#mlmm-tools` on the **Materials Science Discord** for community help.

---

## 10  Citing this Workflow

```
@software{mlmm_tools_2025,
  title  = {{ML/MM tools}: Definitive Workflow Tutorial},
  author = {Ohmura, T. *et al.*},
  year   = {2025},
  url    = {https://github.com/t-0hmura/mlmm_tools}
}
```

&copy; 2025  The **ML/MM tools** team.  Released under the MIT License.

---

*Last updated 2025‑07‑07*
