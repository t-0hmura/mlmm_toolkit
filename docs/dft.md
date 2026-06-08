# `dft`

Run a single-point DFT calculation on the ML region using GPU4PySCF (or CPU PySCF), then recombine the high-level energy with MM evaluations to obtain the ML(dft)/MM total energy. The default functional/basis is `wb97m-v/def2-tzvpd`. Results include energy and population analysis (Mulliken, meta-Lowdin, IAO charges).

```
E_total = E_REAL_low + E_ML(DFT) - E_MODEL_low
```

## When to use

- Use when refining stationary-point energies (R / TS / P / IM) at DFT level after an MLIP path search, or sanity-checking an MLIP barrier against a benchmark functional / basis.

## Quick examples

```bash
# Minimal single-point DFT on the ML region
mlmm dft -i enzyme.pdb --parm real.parm7 --model-pdb ml_region.pdb \
 -q 0 -m 1 --out-dir ./result_dft
```

```bash
# Change functional/basis for a higher-level single point
mlmm dft -i enzyme.pdb --parm real.parm7 --model-pdb ml_region.pdb \
 -q 0 -m 1 --func-basis "wb97m-v/def2-tzvpd" --out-dir ./result_dft_tz
```

```bash
# Freeze selected atoms in the ML/MM setup before DFT
mlmm dft -i enzyme.pdb --parm real.parm7 --model-pdb ml_region.pdb \
 -q -1 -m 2 --freeze-atoms "1,3,5" --out-dir ./result_dft_freeze
# tighten SCF instead: -q 0 -m 1 --conv-tol 1e-10 --max-cycle 200 --out-dir ./result_dft_tight
```

## Inputs

Command form:

```bash
mlmm dft -i enzyme.pdb --parm real.parm7 --model-pdb ml_region.pdb -q 0 -m 1 [options]
```

`mlmm dft --help` shows core options; `mlmm dft --help-advanced` shows the full option list.

| Input | Required | Notes |
| --- | --- | --- |
| `-i, --input` | yes | Full enzyme structure (PDB or XYZ). If XYZ, use `--ref-pdb` for topology. |
| `--parm` | yes | Amber parm7 topology for the full system. |
| `--model-pdb` | optional | PDB defining the ML region (atom IDs must match the enzyme PDB). Optional when `--detect-layer` is enabled. |
| `--model-indices` | optional | Comma-separated atom indices for the ML region (ranges allowed, e.g. `1-5`). Used when `--model-pdb` is omitted. |
| `--ref-pdb` | for XYZ inputs | Reference PDB topology when input is XYZ. |
| `-q, --charge` | yes (unless `-l`) | Charge of the ML region. Required unless `-l/--ligand-charge` is given (PDB input or XYZ with `--ref-pdb`). |
| `-m, --multiplicity` | optional | Spin multiplicity (2S+1) for the ML region (default `1`). |

## Workflow

1. **Input handling** -- The full enzyme PDB (`-i`), Amber topology (`--parm`), and ML-region definition (`--model-pdb` or `--model-indices` or B-factor detection via `--detect-layer`) are loaded. Link hydrogens are appended automatically (C/N parents within 1.7 Å) unless explicit `link_mlmm` pairs are provided via YAML.
2. **SCF build** -- `--func-basis` is parsed into functional and basis. The GPU4PySCF backend is used when available; closed-shell GPU runs additionally use the low-memory `gpu4pyscf.dft.rks_lowmem.RKS` SCF when `--lowmem` is on (default). Use `--engine cpu` to force CPU mode. (For the SCF JK / `density_fit()` behaviour see the `--lowmem` row in the CLI options table.) When `--embedcharge` is enabled, MM point charges from the Amber topology are embedded into the QM Hamiltonian via `pyscf.qmmm.mm_charge()`, so the DFT wavefunction is self-consistently polarized by the MM environment.
3. **ML(dft)/MM recombination** -- After the DFT converges, MM evaluations of the full system (REAL-low) and the ML subset (MODEL-low) are computed. The combined energy is reported in Hartree and kcal/mol.
4. **Population analysis & outputs** -- Mulliken, meta-Lowdin, and IAO charges and spin densities (UKS only) are written alongside the combined energy block in `result.yaml`.

## Outputs

```
out_dir/ (default: ./result_dft/)
├── ml_region_with_linkH.xyz    # ML-region coordinates (with link-H) used for DFT
├── result.yaml                 # DFT + ML(dft)/MM energy summary, charges, spin densities
├── result.json                 # only when --out-json is passed
└── (stdout)                    # Pretty-printed configuration blocks and energies
```

- `result.yaml` expands to:
  - `energy`: Hartree/kcal/mol values, convergence flag, wall time, backend info (`engine`: `gpu4pyscf(rks_lowmem)` / `gpu4pyscf` / `pyscf(cpu)`; `used_gpu`; `used_lowmem`).
  - `mlmm_energy`: REAL-low / MODEL-low MM evaluations and the recombined `E_total = E_REAL_low + E_ML(DFT) - E_MODEL_low` in Hartree and kcal/mol.
  - `charges`: Mulliken, meta-Lowdin, and IAO atomic charges (`null` when a method fails).
  - `spin_densities`: Mulliken, meta-Lowdin, and IAO spin densities (UKS-only for spins).
- It also summarizes charge, multiplicity, spin (2S), functional, basis, convergence knobs, and resolved output directory.

## CLI options

The full flag list is in the generated [command reference](reference/commands/index.md); the table below covers the options that need explanation.

| Option | Description | Default |
| --- | --- | --- |
| `-i, --input PATH` | Full enzyme structure (PDB or XYZ). If XYZ, use `--ref-pdb` for topology. | Required |
| `--ref-pdb FILE` | Reference PDB topology when input is XYZ. | _None_ |
| `--parm PATH` | Amber parm7 topology for the full system. | Required |
| `--model-pdb PATH` | PDB defining the ML region (atom IDs must match the enzyme PDB). Optional when `--detect-layer` is enabled. | _None_ |
| `--model-indices TEXT` | Comma-separated atom indices for the ML region (ranges allowed, e.g. `1-5`). Used when `--model-pdb` is omitted. | _None_ |
| `--model-indices-one-based / --model-indices-zero-based` | Interpret `--model-indices` as 1-based or 0-based. | `True` (1-based) |
| `--detect-layer / --no-detect-layer` | Detect ML/MM layers from input PDB B-factors (B=0/10/20). | `True` |
| `-q, --charge INT` | Charge of the ML region. Required unless `-l/--ligand-charge` is given (PDB input or XYZ with `--ref-pdb`). | Required unless `-l/--ligand-charge` is provided |
| `-l, --ligand-charge TEXT` | Total charge or per-resname mapping (e.g. `SAM:1,GPP:-3`) used to derive the ML-region charge when `-q` is omitted (requires PDB input or `--ref-pdb`). | _None_ |
| `-m, --multiplicity INT` | Spin multiplicity (2S+1) for the ML region. | `1` |
| `--freeze-atoms TEXT` | Comma-separated 1-based indices to freeze (e.g. `"1,3,5"`). Merged with YAML `geom.freeze_atoms`. | _None_ |
| `--func-basis TEXT` | Functional/basis pair as `"FUNC/BASIS"`. | `wb97m-v/def2-tzvpd` |
| `--max-cycle INT` | Maximum SCF iterations. | `100` |
| `--conv-tol FLOAT` | SCF convergence tolerance (Hartree). | `1e-9` |
| `--grid-level INT` | DFT integration grid level (0=coarse, 3=default, 5=fine, 9=very fine). | `3` |
| `--engine {gpu,cpu}` | Force GPU4PySCF (`gpu`) or CPU PySCF (`cpu`); `gpu` raises an error if GPU4PySCF is unavailable. | `gpu` |
| `--lowmem/--no-lowmem` | Use `gpu4pyscf.dft.rks_lowmem.RKS` for closed-shell GPU runs (memory-efficient direct JK; `mlmm dft` does not call `density_fit()` on either path). Open-shell, CPU, or pre-`rks_lowmem` GPU4PySCF auto-fall back to standard RKS/UKS. | `True` |
| `-o, --out-dir DIR` | Output directory. | `./result_dft/` |
| `--config FILE` | Base YAML configuration file applied before explicit CLI options. | _None_ |
| `--show-config/--no-show-config` | Print resolved configuration and continue execution. | `False` |
| `-b, --backend CHOICE` | Metadata-only label recorded in output; does not select a calculator in `dft` (the ML region is computed with DFT/PySCF, and the MM low level is chosen via `--mm-backend`): `uma` (default), `orb`, `mace`, `aimnet2`. | `uma` |
| `--embedcharge/--no-embedcharge` | Electrostatic embedding: MM point charges polarize the DFT wavefunction (see Workflow step 2). | `False` |
| `--embedcharge-cutoff FLOAT` | Cutoff radius (Å) for embed-charge MM atoms. | `12.0` |
| `--link-atom-method {scaled,fixed}` | Link-atom placement: `scaled` (g-factor, Gaussian ONIOM standard) or `fixed` (legacy 1.09 Å for C, 1.01 Å for N). | `scaled` |
| `--mm-backend {hessian_ff,openmm}` | MM backend for the low-level ONIOM evaluation: `hessian_ff` (analytical, default) or `openmm` (finite-difference Hessian). | `hessian_ff` |
| `--cmap/--no-cmap` | Enable CMAP (backbone cross-map dihedral correction) in model parm7. Default: disabled (consistent with Gaussian ONIOM). | `--no-cmap` |
| `--out-json/--no-out-json` | Write a machine-readable `result.json` to `out_dir`. | `False` |
| `--dry-run/--no-dry-run` | Validate options and print execution plan without running DFT. Shown in `--help-advanced`. | `False` |
| `--convert-files/--no-convert-files` | Toggle XYZ/TRJ to PDB companions when a PDB template is available. | `True` |

## YAML configuration

Accepts a mapping root; the `dft` section (and optional `geom`, `calc`/`mlmm`) is applied when present. Merge order is:
- defaults
- `--config`
- explicit CLI options

```yaml
geom:
 coord_type: cart                  # optional geom_loader settings
calc:
 model_charge: 0                   # ML region charge
 model_mult: 1                     # spin multiplicity 2S+1
mlmm:
 real_parm7: real.parm7            # Amber parm7 topology
 model_pdb: ml_region.pdb          # ML-region definition
dft:
 func_basis: wb97m-v/def2-tzvpd      # exchange-correlation functional / basis set
 conv_tol: 1.0e-09                # SCF convergence tolerance (Hartree)
 max_cycle: 100                    # maximum SCF iterations
 grid_level: 3                     # PySCF grid level
 verbose: 0                        # PySCF verbosity (0-9); CLI -v 2/3 raises runtime PySCF verbosity to >=4
 out_dir: ./result_dft/            # output directory root
```

Full schema (every key and default): [YAML Reference](yaml-reference.md).

## Notes

- A matching def2 effective core potential is auto-attached whenever the basis name begins with `def2` (no element-presence check).
- **Blackwell-architecture GPUs** (RTX 50xx): GPU4PySCF may fail with out-of-memory errors even for small systems (~100 atoms). Use `--engine cpu` or an external DFT program (ORCA, Gaussian) for production calculations on these GPUs.
- **Out-of-memory with def2-TZVPD**: The default basis set `def2-tzvpd` is large and may cause OOM for systems with >150 atoms on 16–24 GB GPUs. Use `--func-basis 'wb97m-v/def2-svp'` as a practical alternative; barrier height errors between def2-SVP and def2-TZVPD are typically 1–3 kcal/mol.
- Compiled GPU4PySCF wheels may not support non-x86 systems; build from source in that case (see https://github.com/pyscf/gpu4pyscf).

## See Also

- [Common Error Recipes](recipes-common-errors.md) — Symptom-first failure routing
- [Troubleshooting](troubleshooting.md) — Detailed troubleshooting guide
- [freq](freq.md) — Vibrational frequency analysis (often precedes DFT refinement)
- [opt](opt.md) — Single-structure geometry optimization
- [all](all.md) — End-to-end workflow with `--dft`
- [YAML Reference](yaml-reference.md) — Full `dft` configuration options
- [Glossary](glossary.md) — Definitions of DFT, SP (Single Point)
