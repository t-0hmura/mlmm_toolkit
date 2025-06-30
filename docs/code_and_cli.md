# Code Overview and CLI Tools

This document provides an overview of the ML/MM calculator code features and detailed usage for the bundled command-line utilities.

## Code Features

The ML/MM hybrid calculator combines a machine-learning potential with classical molecular mechanics. It can compute:

- **Energies**
- **Forces**
- **Hessians** – analytical in the ML region and numerical or padded in the MM region

Covalent bonds cut at the QM/MM boundary are capped with hydrogen link atoms so that entire side chains can be included in the ML region.

Supported ML potentials:

| Potential | Repository |
|-----------|-----------|
| **AIMNet2** | <https://github.com/isayevlab/aimnetcentral> |
| **UMA** (uma-s-1)   | <https://github.com/facebookresearch/fairchem> |

The MM layer uses **OpenMM** and works with any force field capable of generating parameters.

## Command-Line Utilities

The package installs several small utilities which become available on your `$PATH`.

| Tool | Purpose | Typical use‑case |
|------|---------|------------------|
| `def_ml_region` | Build an ML region with residues around one or more substrate in a protein–substrate complex. | Preparing the subsystem for ML/MM calculator |
| `xyz_geom2pdb`  | Convert an XYZ geometry or trajectory to a multi‑model PDB while borrowing atom/residue metadata from a reference PDB. | Exporting optimization results for visualization |
| `add_elem_info` | Append element symbols (PDB columns 77–78) | Fixing element fields before running external tools |
| `get_freeze_indices` | List atom indices to *freeze* based on their distance from the ML region. | Constraining outer‑shell atoms during local relaxations |
| `bond_scan` | Scan a bond length with ML/MM optimization at each step. | Generating aligned structures along a reaction coordinate |
| `ts_search` | Dimer‑based TS search with partial Hessian updates. | Locating transition states in large systems |
| `energy_summary` | Compute ΔE/ΔG tables and plots from reactant, TS and product structures. | Summarizing reaction energetics |
| `trj2fig` | Plot ΔE from an XYZ trajectory and export the highest peak frame. | Visualizing optimization or scan profiles |

### `def_ml_region` – Automated ML‑region extractor

```bash
def_ml_region --real complex.pdb --center ligand.pdb --out ml_region.pdb --radius_ml 2.6 --radius_het 3.6 --include_H2O false --exclude_backbone true
```

| Option | Default | Meaning |
|--------|---------|---------|
| `-r, --real` | *(required)* | Full protein–substrate complex (PDB) from which the ML region will be carved out. |
| `-c, --center` | *(required)* | PDB containing the ligand/substrate whose heavy atoms act as distance centers. |
| `-o, --out` | `ml_region.pdb` | Output file containing only the atoms selected for the ML region. |
| `--radius_ml` | `2.6` Å | Residues with any atom ≤ radius from the substrate are included. |
| `--radius_het` | `3.6` Å | Cut‑off distance between hetero atoms in the center residues and those in other residues. |
| `--include_H2O` | `false` | Keep water molecules inside the ML region. |
| `--exclude_backbone` | `true` | When `true`, only side chains are kept; `false` preserves backbone connectivity. |

Example:

```bash
def_ml_region --real complex.pdb --center ligand.pdb --out ml_region.pdb --radius_ml 6.0
```

This generates `ml_region.pdb` containing all protein residues and HETATMs within 6.0 Å of any atom in `ligand.pdb`.

### `xyz_geom2pdb` – XYZ → PDB converter

```bash
xyz_geom2pdb --input traj.xyz --ref ref.pdb --output traj.pdb
```

If `--output` is omitted the file is written as `<input‑stem>.pdb`. Multi‑model XYZ files are automatically wrapped in MODEL/ENDMDL blocks for playback in graphics programs.

| Option | Default | Meaning |
|--------|---------|---------|
| `-i, --input` | *(required)* | Input XYZ geometry or trajectory |
| `-r, --ref` | *(required)* | Reference PDB file providing topology |
| `-o, --output` | `<input>.pdb` | Output PDB file |

### `add_elem_info` – Append element columns 77–78

```bash
add_elem_info protein.pdb       # over‑write in place
add_elem_info input.pdb output_with_elem.pdb
```

| Argument | Default | Meaning |
|----------|---------|---------|
| `input` | *(required)* | Input PDB file |
| `output` | same as `input` | Output file name (in-place if omitted) |

### `get_freeze_indices` – Distance‑based atom freeze list

```bash
get_freeze_indices --real complex.pdb --model ml_region.pdb --range 8.0 --out freeze.txt
```

| Option | Default | Meaning |
|--------|---------|---------|
| `--real` | *(required)* | Full complex PDB path |
| `--model` | *(required)* | ML‑region PDB |
| `--range` | *(required)* | Relaxed range in Å from ML-region |
| `--out` | `stdout` | Output text file of indices |

Indices are written in a line. Omit `--out` to print them to stdout and pipe the list directly into `ase.constraints.FixedAtoms` or other tools.

### `bond_scan` – Single‑bond scan

```bash
bond_scan scan.yml
```

Runs an inward and outward scan along a chosen bond. Each scan point is optimized with LBFGS and saved to `final_geometries.trj`.

| Key | Default | Meaning |
|-----|---------|---------|
| `geom.fn` | `./coord/gs.xyz` | Starting geometry (XYZ or PDB) |
| `scan.bond` | `[0, 1]` | Pair of 0‑based atom indices to scan |
| `scan.step` | `0.05` Å | Displacement per step |
| `scan.scan_range` | `[1.4, 4.0]` Å | Minimum and maximum bond length |
| `scan.init_thresh` | `gau` | Threshold for the first optimization |
| `scan.thresh` | `gau` | Threshold for subsequent steps |
| `scan.max_cycles` | `10000` | LBFGS cycle limit |
| `scan.out_dir` | `./dump/cart_scan/` | Directory for the trajectory |
| `scan.out_fn` | `final_geometries.trj` | Trajectory file name |
| `calc.*` | – | Options forwarded to the ML/MM calculator |

See `examples/cartesian_scan.yml` for a complete template.

### `ts_search` – Partial-Hessian ML/MM Dimer search   

```bash
ts_search config.yml
```

The script drives a **mixed ML/MM transition-state search** that combines  

1. *Partial* Hessians (ML region + nearby MM; distant MM rows/cols = 0)  
2. Periodic Dimer optimization (loose → final thresholds)  
3. A **mass-scaled multi-mode flatten loop** that flattens extra imaginary modes  
    until only the transition-state mode remains.

Before launching the Dimer loops the geometry of **MM region** is relaxed with LBFGS. The
cycle limit for this pre-optimization is controlled by `dimer.max_cycles_preopt`
(defaults to `dimer.max_cycles`).

---

## YAML configuration keys

| Section · Key | Default | Meaning |
|---------------|---------|---------|
| **geom.fn** | `./coord/gs_hei.xyz` | Input geometry (`.xyz` or `.pdb`) |
| **geom.freeze_atoms** | `[]` | List of indices kept fixed *throughout* |
| **geom.real_pdb** | `./parm/complex.pdb` | Template for PDB overlays |
| **misc.out_dir** | `./dump/dimer/` | optimization output directory |
| **misc.vib_dir** | `./dump/vib/` | Directory for vibrational analysis |
| **misc.dump** | `false` | Write intermediate trajectories |
| **dimer.thresh_loose** | `gau_loose` | Loose LBFGS tolerance |
| **dimer.thresh** | `baker` | Final LBFGS tolerance |
| **dimer.update_interval_hessian** | `50` | Steps between Hessian rebuilds |
| **dimer.lobpcg** | `true` | Use LOBPCG for lowest-mode search |
| **dimer.max_cycles** | `100000` | Global LBFGS cycle cap |
| **dimer.max_cycles_preopt** | same as `dimer.max_cycles` | LBFGS cap for the pre‑optimization step. 0 to no preopt. |
| **dimer.partial_mm_cutoff** | `0.0` | Å radius for *partial* Hessian |
| **dimer.neg_freq_thresh** | `10.0` | cm⁻¹ cutoff for “imaginary” |
| **flatten.amp_ang** | `0.20` | Carbon displacement amplitude (Å) |
| **flatten.max_iter** | `20` | Max flatten iterations |
| **flatten.sep_cutoff** | `2.0` | Å between representative atoms of modes |
| **flatten.k** | `10` | Number of representative atoms to displace |
| **calc.* ** | – | ML/MM calculator kwargs (must include `model_pdb`) |
| **dimer_kwargs.* ** | – | Extra kwargs forwarded to Pysisyphus `Dimer` |

---

## Minimal example `config.yml`

```yaml
geom:
    fn: "./coord/start.xyz"
    freeze_atoms: [0, 1, 2]
    real_pdb: "./parm/complex.pdb"
misc:
    out_dir: "./work/dimer/"
    vib_dir: "./work/vib/"
    dump: true
dimer:
    thresh_loose: "gau_loose"
    thresh: "baker"
    update_interval_hessian: 50
    partial_mm_cutoff: 0.0
flatten:
    amp_ang: 0.20        # optional – default shown
    max_iter: 20
    sep_cutoff: 2.0
calc:
    model_pdb: "./parm/model.pdb"
    # other MLMM-specific settings …
```

---

*For a richer template see* `examples/ts_search.yml`.

### `energy_summary` – Gibbs‑energy profiler

```bash
energy_summary reaction.yml
```

Computes electronic and vibrational contributions for reactant, transition state and product structures, producing a ΔE/ΔG table and optional energy diagrams.

| Key | Default | Meaning |
|-----|---------|---------|
| `geom.reac` | *(required)* | Reactant structure |
| `geom.ts` | *(required)* | Transition‑state structure |
| `geom.prod` | *(required)* | Product structure |
| `calc.*` | – | ML/MM calculator options |
| `--temp` | `298.15` K | Temperature for vibrational corrections |
| `--no_annotation` | *(flag)* | Do not label energy diagrams |

A ready‑to‑use YAML file is provided as `examples/energy_summary.yml`.

### `trj2fig` – Plot energy from an XYZ trajectory

```bash
trj2fig -i traj.xyz -o energy.png --output-peak ts.pdb --reverse-x
```

Extract energies from the two-line comments of an XYZ trajectory, plot the ΔE profile and optionally write the peak structure or a CSV table.

| Option | Default | Meaning |
|--------|---------|---------|
| `-i, --input` | *(required)* | XYZ trajectory file |
| `-o, --out` | `energy.html` | Output (.html/.png/.svg/.pdf/.csv) |
| `--unit` | `kcal` | ΔE unit (kcal/hartree) |
| `-r, --reference` | `0` | Reference frame index |
| `--output-peak` | – | Write the highest peak frame |
| `--reverse-x` | *(flag)* | Reverse the x-axis |

---
