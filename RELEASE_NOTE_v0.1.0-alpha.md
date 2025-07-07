# ML/MM Tools — v0.1.0‑alpha Release Notes  
*(First public alpha)*  

**Release date:** 7 July 2025  

---

## Overview  
**ML/MM Tools** is an open‑source QM/MM framework that replaces the quantum region with machine‑learning potentials (UMA or AIMNet 2) while treating the environment with OpenMM. This alpha delivers a working proof‑of‑concept for enzymatic mechanism studies.

---

## Key Features  

* **Dual ML back‑ends – UMA & AIMNet 2**
* **Automatic link‑atom boundaries**
* **ASE calculator & PySiSyphus driver**  
* **Partial‑Hessian Dimer (PH‑Dimer) optimiser**  
* **CLI helpers**: `def_ml_region`, `ts_search`, `energy_summary`, and others.

---

## Requirements  

| Component | Minimum |
|-----------|---------|
| Python    | 3.11    |
| CUDA      | 12.6 (12.8 for RTX 50‑series) |
| OS        | Linux / Windows 11 + WSL 2 |

---

## Quick Installation  

For CUDA 12.6:
```bash
pip install fairchem-core==2.3.0
pip install git+https://github.com/isayevlab/aimnetcentral.git
pip install git+https://github.com/t-0hmura/mlmm_tools.git
huggingface-cli login
```

For CUDA 12.8 (required for RTX 50 series):
```bash
pip install fairchem-core==2.3.0
pip install git+https://github.com/isayevlab/aimnetcentral.git
pip install --force-reinstall torch==2.7.0 --index-url https://download.pytorch.org/whl/cu128
pip install git+https://github.com/t-0hmura/mlmm_tools.git
huggingface-cli login
```
---

## Limitations (Alpha)  

1. Periodic boundary conditions and full molecular‑dynamics drivers are not yet implemented.  
2. Single‑GPU execution only; multi‑GPU testing pending.  
3. The public API will evolve and may change without notice.  

---

## Licence & Repository  

* **Licence:** MIT  
* **Source:** <https://github.com/t-0hmura/mlmm_tools>  
* **Citation:** Methodology paper to appear on *ChemRxiv*.  

*End of notes.*  
