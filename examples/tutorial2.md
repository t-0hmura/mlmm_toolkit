
# ML/MM Workflow Tutorial — **Concise Reproducible Guide**  
_Last updated: 2025‑07‑07_

This document integrates the **directory clarity** of *tutorial (3)*, the **robust troubleshooting & advanced tips** of *tutorial (2)*, and the **step‑wise explanatory depth** of *tutorial (1)* into **one self‑contained recipe**.  
It walks you from a prepared enzyme–substrate snapshot to activation‑ and reaction‑free energies (**ΔG‡**, **ΔG**) using the `mlmm_tools` package (AIMNet2 or UMA + OpenMM).

---

## 0 Prerequisites (one‑line recap)  

```bash
pip install git+https://github.com/t-0hmura/mlmm_tools.git  # plus AIMNet2 / UMA back‑ends
```

> Detailed CUDA‑specific wheels and dependency notes are in the project README.

---

## 1 Test Systems & Hardware

| ID | Enzyme | Reaction | Reference | Atoms (REAL / ML) |
|----|--------|----------|-----------|-------------------|
| **CM** | Chorismate Mutase (*B. subtilis*) | Claisen rearrangement | PDB `1COM` | 12 986 / 208 |
| **Fur6** | SAM‑dependent methyltransferase | CH₃ transfer (SAM → substrate) | Zhao *et al.* 2024 | 13 993 / 260 |

*Benchmark workstation:* Ryzen 7950X3D (32 threads), 128 GB RAM, RTX 5080 (16 GB VRAM).

---

## 2 Directory Layout (clone `examples/` here)

```
tutorial/
├── cm/           # Chorismate Mutase example
│   ├── coord/    # Pre‑made coordinates (MD snapshot, ML region, …)
│   ├── yaml/     # 10 YAML job files
│   └── run.sh    # Kick‑off script
├── fur6/         # Fur6 methyltransferase example
└── common/       # Helper scripts (def_ml_region, trj2fig, …)
```

---

## 3 System Preparation

1. **Trim the MD snapshot**  
   Keep protein, substrate and waters within **6 Å** of the substrate.

   ```bash
   cpptraj -p full.parm7 -y full.nc <<'EOF'
   autoimage
   strip !(:1-XYZ<6.0)
   trajout REAL.pdb pdb
   trajout REAL.parm7 parm7
   trajout REAL.rst7 restart
   EOF
   ```

2. **Define ML region**  

   ```bash
   def_ml_region -r REAL.pdb -c substrate_only.pdb \
                 --radius_ml 6.0 --exclude_backbone true \
                 -o ml_region.pdb
   ```

3. **Freeze outer shell (optional but speeds up Hessians)**  

   ```bash
   get_freeze_indices REAL.pdb ml_region.pdb --cutoff 10.0 > freeze.txt
   ```

---

## 4 Running the Workflow

### 4.1 One‑command run

```bash
bash run.sh          # executes 1_opt.yaml → … → 10_energy_summary.yaml
```

### 4.2 What happens under the hood?

| Stage | YAML / CLI | Description |
|-------|------------|-------------|
| (i)  | `1_opt.yaml`               | Optimise initial snapshot |
| (ii) | `2_scan.yaml` + `bond_scan`| 1‑D scan → pick Reactant/Product |
| (iii)| `3_opt.yaml`, `4_opt.yaml`| Optimise Reactant & Product |
| (iv) | `5_gs.yaml`               | GSM minimum‑energy path |
| (v)  | `6_tsopt.yaml` + `ts_search`| PH‑Dimer TS refinement |
| (vi) | `7_irc.yaml`              | IRC from TS |
| (vii)| `8_opt.yaml`, `9_opt.yaml`| Final Reactant / Product opt. |
| (viii)| `10_energy_summary.yaml`  | ΔE, ΔG tables & plots |

Intermediate geometries live in `dump/` (mirrored to `coord/` for chaining).

---

## 5 Wall‑Clock Timing (UMA backend)

| Stage | CM | Fur6 |
|-------|-----------|-----------|
| Initial opt. | 0 h 06 m 45 s | 0 h 08 m 23 s |
| Scan | 0 h 06 m 54 s | 0 h 05 m 46 s |
| Candidate opt. | 0 h 01 m 32 s | 0 h 02 m 18 s |
| GSM | 0 h 05 m 14 s | 0 h 06 m 01 s |
| **PH‑Dimer TS** | **0 h 33 m 50 s** | **0 h 49 m 00 s** |
| IRC | 0 h 05 m 01 s | 0 h 07 m 31 s |
| Endpoint opt. | 0 h 01 m 56 s | 0 h 02 m 27 s |
| Thermochem. | 0 h 13 m 35 s | 0 h 20 m 50 s |
| **Total** | **1 h 14 m 47 s** | **1 h 42 m 16 s** |

**Quick sanity‑check:** set `vib_run: false` in `6_tsopt.yaml` to slash runtime to ≈ 20 min.

---

## 6 Energetic Results

| Reaction | ΔE‡ | ΔE | ΔG‡ | ΔG | ΔE‡<sub>QM/MM</sub> | ΔG‡<sub>exp.</sub> |
|----------|----:|---:|----:|---:|---------------------:|-------------------:|
| CM | 18.1 | –23.4 | 14.7 | –23.3 | 16.1 | 15.4 |
| Fur6 | 10.2 | –62.5 | 9.6 | –58.3 | 8.6 | — |

The UMA‑based ML/MM workflow reproduces reference QM/MM and experimental barriers to **within ≈ 2 kcal mol⁻¹**.

---

## 7 Inspecting the Output

* **`dump/gs/`** – GSM trajectory (`gs.trj`)  
* **`dump/dimer/`** – Optimised TS (`*_tsopt_final_geometry.xyz`)  
* **`dump/irc/`** – IRC forward / backward paths  
* **`tables/ΔG_summary.csv`** – Collated energies  
* **`fig/energy_profile.png`** – Potential‑ & free‑energy plots  

Open XYZ files in PyMOL / VMD for structural insight.

---

## 8 Switching to AIMNet2

Change `backend: aimnet2` in every YAML, adjust `aimnet_model`, and rerun `run.sh`.  
ΔG‡ typically shifts by ≤ 1 kcal mol⁻¹.

---

## 9 Troubleshooting at a Glance

| Symptom | Cause | Fix |
|---------|-------|-----|
| **CUDA OOM** during TS or IRC | VRAM exceeded | Lower `ml_batch_size` or set `ml_device: cpu` for those steps |
| **TS search stalls** | Poor TS guess | Increase `n_hess_update`, reduce `dimer_step`, or supply a GSM peak geometry closer to the saddle |
| **“atoms out of sync”** | Atom order mismatch | Regenerate `ml_region.pdb` *after* final REAL trimming |
| **`torch.cuda.OutOfMemoryError` during vib analysis** | Hessian too large for GPU | Set `vib_run: false` or `ml_device: cpu` |

---

## 10 Next Steps & Support

* Periodic‑boundary MD and on‑the‑fly sampling are under active development—watch the GitHub repo.  
* Bug reports → GitHub Issues • Questions → `#mlmm-tools` Discord channel.

---

## 11 Citing the Workflow

```
@software{mlmm_tools_2025,
  title  = {{ML/MM tools}: v2025.4},
  author = {Ohmura, T. and colleagues},
  year   = 2025,
  url    = {https://github.com/t-0hmura/mlmm_tools}
}
```

Happy modelling! 🚀
