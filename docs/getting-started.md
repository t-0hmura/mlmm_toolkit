# Getting Started

::::{container} p2r-intro
<img src="./overview.jpg" alt="mlmm-toolkit workflow overview" class="p2r-intro-figure">

:::{container} p2r-intro-text
`mlmm-toolkit` is a Python command-line toolkit that uses ML/MM (machine learning / molecular mechanics) to **search automatically for candidate enzyme reaction pathways, starting from PDB / mmCIF structures**.

ML/MM works like QM/MM, with a machine-learning interatomic potential (MLIP) in place of QM. The MLIPs are neural networks trained on DFT data; they approximate a DFT-level potential energy surface at a tiny fraction of the cost. The reacting part of the enzyme (the ML region) is computed with the MLIP, and the protein around it with an Amber force field (MM). The two are combined by the ONIOM subtraction:

```text
E_total = E_REAL_low + E_MODEL_high - E_MODEL_low
```

REAL is the full system and MODEL the ML region; high is the MLIP, and low is the MM backend. The full system is computed with MM, the ML region with both the MLIP and MM, and the MM energy of the ML region is subtracted so that it is not counted twice. Where the ML region cuts a covalent bond, a link hydrogen caps it.

The MM atoms form two layers: Movable-MM atoms relax during optimizations, and Frozen-MM atoms farther out stay fixed. The layers are stored in the B-factor column of the PDB. See [Building the ML region and layers](model-setup.md) for the layers and [ML/MM Calculator](mlmm-calc.md) for the energy, forces, and Hessian.

In many cases, a **single command** like this one gives a first draft of the reaction pathway:

```bash
mlmm -i 1.R.pdb 3.P.pdb -c 'SAM,GPP,MG' -l 'SAM:1,GPP:-3'
```
:::
::::

---

Add `--tsopt --thermo --dft` and the same run continues automatically through **minimum energy path (MEP) search → transition state (TS) optimization → intrinsic reaction coordinate (IRC) → vibrational analysis and thermochemical correction → DFT single points**.

```bash
mlmm -i 1.R.pdb 3.P.pdb -c 'SAM,GPP,MG' -l 'SAM:1,GPP:-3' --tsopt --thermo --dft
```

---

> **Examples:** the [`examples/beza/`](https://github.com/t-0hmura/mlmm_toolkit/tree/main/examples/beza) directory holds the structures used above (`1.R.pdb`, `3.P.pdb`) and a workflow script (MEP search and scan pipelines) built around the GPP C6-methyltransferase BezA ([Tsutsumi et al., *Angew. Chem. Int. Ed.* 2022, 61, e202111217](https://doi.org/10.1002/anie.202111217)). After [installation](installation.md), get it with `git clone https://github.com/t-0hmura/mlmm_toolkit && cd mlmm_toolkit/examples/beza` and run the commands above there.

---

## What it is for

* **Trial and error on reaction mechanisms** in the full enzyme, where QM/MM with DFT takes too long to check
* **Starting structures** for QM/MM (reactant, TS, and product of the full system; [`oniom-export`](oniom-export.md) turns them into Gaussian ONIOM or ORCA QM/MM input)
* **Many reaction-path calculations** across substrate variants and enzyme mutants

## What it automates

Provide one of three inputs: (1) several PDB structures in reaction order (R → … → P), (2) one structure plus a scan, or (3) one structure plus TS optimization. `mlmm-toolkit` then handles the following automatically.

1. **ML region**: cuts out the active site (binding pocket) around the specified substrates as the ML region
2. **MM topology and layers**: builds the Amber topology of the full system with `mm-parm` (AmberTools) and assigns the ML, Movable-MM, and Frozen-MM layers with `define-layer`
3. **Minimum energy path (MEP) search**: searches the pathway with the Growing String Method (GSM) or Direct Max Flux (DMF)
4. **High-accuracy checks**: TS optimization, IRC, vibrational analysis, and DFT single points

The ML region uses **UMA** (Meta) by default; [`-b/--backend`](backends.md) also selects **ORB**, **MACE**, and **AIMNet2**.

Once MLIP/MM has found a reasonable pathway, `mlmm-toolkit` can take its TS straight into a DFT/MM TS optimization. It runs the TS optimization → IRC → endpoint optimization → frequency workflow with GPU-accelerated DFT through GPU4PySCF. See [Refine an MLIP TS with DFT](dft-backend.md) for details.

> To run a [model you built yourself](model-setup.md#use-a-model-you-built-yourself) as is, omit `-c`.

---

## Workflow and pipeline

### The pipeline

The `all` subcommand (the default) runs the whole workflow in one go, stage by stage in this order:

```text
Input structure(s) (PDB / mmCIF)
  │
  ▼
[extract] extraction: cut out the ML region around the substrates (only with -c)
  │
  ▼
[mm-parm] MM topology: build the Amber parm7/rst7 of the full system (skipped with --parm7)
  │
  ▼
[define-layer] layers: write the ML / Movable-MM / Frozen-MM layers into the B-factor column
  │
  ▼
[scan] scan: staged scan of distances, angles, or dihedrals (only with -s)
  │
  ▼
[path-opt / path-search] path search: find the MEP (minimum energy path); skipped in TS-only mode
  │
  ▼
[tsopt] TS optimization: refine the transition state (only with --tsopt)
  │
  ▼
[irc] IRC: follow the intrinsic reaction coordinate and optimize its endpoints (only with --tsopt)
  │
  ▼
[freq] vibrational analysis: compute the thermochemical correction (only with --tsopt --thermo)
  │
  ▼
[dft] DFT single points: compute DFT/MM energies (only with --tsopt --dft)
```

Each stage also runs on its own as a subcommand ([`extract`](extract.md), [`tsopt`](tsopt.md), [`irc`](irc.md), and so on; see the [subcommand list](index.md#subcommands)).

---

## Where to start

For environment setup, see the [Installation guide](installation.md).

* **Try it in a web browser**: the [Colab GUI notebook](https://colab.research.google.com/github/t-0hmura/mlmm_toolkit/blob/main/examples/mlmm_colab.ipynb) (pick the ML region in 3D)
* **Start from several PDB structures**: [Quickstart: `mlmm all`](quickstart-all.md)
* **Explore from one PDB structure with a scan**: [Quickstart: `mlmm all --scan-lists`](quickstart-scan.md)
* **Optimize and check a TS candidate**: [Quickstart: TS-only mode](quickstart-tsopt.md)

---

## How the command works

Installation provides the `mlmm` command. Without a subcommand, `all` runs.

```bash
# These two do the same thing
mlmm [OPTIONS]...
mlmm all [OPTIONS]...
```

### Choosing an input mode

| Mode | Input | What happens |
| --- | --- | --- |
| **Multi-structure MEP search** | Two or more PDBs (`-i R.pdb P.pdb`) | Builds the ML region and layers from the structures and searches the MEP |
| **Single structure + scan** | One PDB + `--scan-lists` (`-s`) | Drives the chosen distances, angles, or dihedrals step by step to build the pathway |
| **TS-only mode** | One PDB + `--tsopt` | Skips the MEP search and goes straight to optimizing the TS candidate and running IRC |

> **Note:** a single-structure input needs either `--scan-lists/-s` or `--tsopt`.

### Choosing between all and individual commands

* **Use `all`** to run model setup → MEP search → TS optimization and IRC → frequencies and DFT in one command, or while you are still exploring and want one command to manage the outputs.
* **Use the individual commands** to run each stage in turn and check its result before the next; for a complex reaction, this often works better than one `all` run. They also fit a custom sequence and a run that reuses the parm7 and layered PDB of an earlier run.

The individual ML/MM commands need the full-system topology (`--parm7`) and the ML region (`--model-pdb`, `--model-indices`, or the B-factor layers of the input); `all` builds both. `-q` is the charge of the ML region, not of the whole system. See {ref}`ML/MM options <mlmm-options>`.

---

## Main CLI options

| Option | Example | Description |
| --- | --- | --- |
| `-i, --input` | `1.R.pdb 3.P.pdb` | Input structure files (PDB / mmCIF); accepts several |
| `-c, --center` | `'SAM,GPP'` / `'A:SAM:123'` | Extraction center (substrate residue names, residue IDs, or a PDB file); the ML region is cut out around it. Without it, no extraction runs, the whole structure is used, and the ML region comes from the B-factor layers or `--model-pdb` (with neither, the run stops with an error) |
| `-l, --ligand-charge` | `'SAM:1,GPP:-3'` | Formal charge of each ligand, as a mapping (standard residues and ions are counted automatically) |
| `-q, --charge` | `-2` | Total charge of the ML region, not of the whole system (set it to override the automatic value) |
| `-m, --multiplicity` | `1` | Spin multiplicity (default `1`, a singlet) |
| `--parm7` | `real.parm7` | Full-system Amber topology to reuse, such as one from an earlier run or from the MD that produced the input snapshot; skips `mm-parm` |
| `--model-pdb` | `ml_region.pdb` | ML region as a PDB; takes precedence over `-c` and the B-factor layers |
| `--tsopt/--no-tsopt` | (flag) | Turns on TS optimization and IRC |
| `--thermo/--no-thermo` | (flag) | Runs vibrational analysis and thermochemical correction with the QRRHO (quasi-rigid-rotor harmonic oscillator) model (with `--tsopt`) |
| `--dft/--no-dft` | (flag) | Runs DFT single points on the resulting structures (with `--tsopt`) |
| `-b, --backend` | `uma` / `orb` / `mace` | Backend for the ML region (default `uma`; `dft` is also available) |

For the syntax rules, see [Common options and selectors](cli-conventions.md); for every option, see the [`all` CLI reference](reference/commands/all.md).

---

## Before you run: the input structures

### 1. Add hydrogens (required)

Input structures must contain **every hydrogen atom**; `all` does not add them. When a structure lacks hydrogens (a crystal structure, for example), add them beforehand with a tool such as these:

| Recommended tool | Example command | Notes |
| --- | --- | --- |
| **reduce** (Richardson Lab) | `reduce input.pdb > output.pdb` | Fast; widely used to add hydrogens to crystal structures |
| **pdb2pqr** | `pdb2pqr --ff=AMBER input.pdb output.pqr` | Adds hydrogens and assigns partial charges |
| **Open Babel** | `obabel input.pdb -O output.pdb -h` | General-purpose cheminformatics toolkit |
| **mm-parm --add-h** | `mlmm mm-parm -i input.pdb --add-h` | Adds hydrogens with PDBFixer at `--ph` (default 7.0) |

`all` fills blank element columns (77–78) by itself; before a standalone command such as `extract`, fill them with [`add-elem-info`](add-elem-info.md). If the PDB has alternate locations (altLoc), keep one per residue with [`fix-altloc`](fix-altloc.md).

### 2. Keep the same atom order (multiple structures)

When the input has several structures, such as a reactant (R) and a product (P), **every structure must list the same atoms in the same order** (only the coordinates differ). Run the hydrogen tool on every structure with the same settings, and in PyMOL tick *Original atom order* when saving. An atom that moves to another residue keeps its residue and atom name from R: in the bundled example, the hydrogen that GPP passes to Glu186 is still `H11` of `GPP 321` in `3.P.pdb`.

### 3. Match `-l` to the hydrogens

Give each ligand the charge that matches the hydrogens in the file. In the bundled example, SAM has 23 hydrogens, so it is `SAM:1`; with 22 hydrogens it would be `SAM:0`. When the charge and the hydrogens do not match, `mm-parm` stops with an electron-count error before it runs `antechamber`.

mmCIF (`.cif`, `.mmcif`) and PDB files beyond the fixed-column limits of the PDB format work with `all` and the calculation commands; the standalone `mm-parm` reads PDB only. See {ref}`mmCIF input <mmcif-input>` for details.

---

## Output files

When the run finishes, the output directory (`-o`, default `./result_all/`) contains the following files. [Output Directory Layout](output-layout.md) lists the main files, and {ref}`JSON Output Reference <summary-json-path-search-all>` every key of `summary.json`.

| File / folder | Contents |
| --- | --- |
| `summary.log` | Text summary (directory layout and progress of each stage) |
| `summary.json` | Machine-readable results (barriers, energies of each state, bond changes) |
| `energy_diagram_*.png` | Energy profile plots (electronic energy / Gibbs-corrected) |
| `mep_trj.pdb` / `mep_trj.cif` | Animated trajectory of the minimum energy path (MEP) |
| `ml_region.pdb`, `mm_parm/`, `layered/` | The ML region, the full-system Amber topology, and the layered full-system PDBs; reuse them with `--model-pdb` and `--parm7` |
| `segments/seg_NN/` | Detailed results for each reaction segment (optimized R/TS/P structures, IRC trajectories, and more; with `--tsopt`) |

At the end of the terminal output, the `Scientific status:` line under `====== Pipeline summary ======` (`scientific_status` in `summary.json`) is `success` when every requested stage converged, otherwise `partial` or `failed` with the reasons in `scientific_status_reasons`. Whether the TS has n_imag = 1 and the endpoints are the intended R and P is for you to check; each quickstart lists the files to open.

---

## AI agent skills

`mlmm-toolkit` ships instructions for AI agents (Claude Code, Codex, Cursor, and others) in the `skills/` directory.

They cover the CLI subcommands, structure input and output, backend installation, TS search strategy, and HPC runs. Load `skills/` into an agent, and it can run and analyze calculations from plain-language instructions. For where to place the files and the full list of skills, see [`skills/README.md`](https://github.com/t-0hmura/mlmm_toolkit/blob/main/skills/README.md). To call the commands as tools from an MCP client, see [MCP server](mcp_server.md).

---

## Troubleshooting and support

If an error occurs during a run, see these pages:

* {ref}`Troubleshooting <troubleshooting-quick-table>`: fixes by error symptom, and solutions for installation and environment problems
* [MLIP Backends](backends.md): choosing a backend and running parallel workers; [Device Configuration & HPC Setup](device-hpc.md) for GPU memory, device settings, and job scripts on clusters

To see every option of a command, use the help options:

```bash
mlmm <subcommand> --help
mlmm all --help-advanced
```

Report unresolved problems and bugs on [GitHub Issues](https://github.com/t-0hmura/mlmm_toolkit/issues).
