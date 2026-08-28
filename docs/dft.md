# `dft`

Run an energy-only single-point DFT calculation on the ML region using GPU4PySCF (or CPU PySCF), then recombine the high-level energy with MM evaluations to obtain the ML(dft)/MM total energy. DFT gradients and forces are not requested. Use it to evaluate stationary-point energies (R / TS / P / IM) at the DFT level after an MLIP path search, or to sanity-check an MLIP barrier against a benchmark functional / basis. The default functional/basis is `wb97m-v/def2-tzvpd`. Results include energy and population analysis (Mulliken, meta-Lowdin, IAO charges).

```
E_total = E_REAL_low + E_ML(DFT) - E_MODEL_low
```

## Examples

Minimal single-point DFT on the ML region:

```bash
# Minimal single-point DFT on the ML region
mlmm dft -i enzyme.pdb --parm real.parm7 --model-pdb ml_region.pdb \
 -q 0 -m 1 --out-dir ./result_dft
```

Change functional/basis for a higher-level single point:

```bash
# Change functional/basis for a higher-level single point
mlmm dft -i enzyme.pdb --parm real.parm7 --model-pdb ml_region.pdb \
 -q 0 -m 1 --func-basis "wb97m-v/def2-tzvpd" --out-dir ./result_dft_tz
```

Tighten the SCF convergence if needed:

```bash
mlmm dft -i enzyme.pdb --parm real.parm7 --model-pdb ml_region.pdb \
 -q 0 -m 1 --conv-tol 1e-10 --max-cycle 200 --out-dir ./result_dft_tight
```

## Workflow

1. **Input handling** -- `MLMMCore` loads the full enzyme PDB (`-i`), Amber topology (`--parm`), and ML-region definition (`--model-pdb` or `--model-indices` or B-factor detection via `--detect-layer`). Unless YAML supplies explicit `link_mlmm` pairs, it appends link hydrogens at parm7 bonds that cross the ML/MM selection; distance is not used to perceive those bonds.
2. **SCF build** -- `--func-basis` is parsed into functional and basis. The GPU4PySCF backend is used when available; closed-shell GPU runs additionally use the low-memory `gpu4pyscf.dft.rks_lowmem.RKS` SCF when `--lowmem` is on (default). Use `--engine cpu` to force CPU mode. The experimental `--embedcharge` option embeds MM point charges directly in the PySCF Hamiltonian; the DFT workflow does not apply the optional xTB correction used by MLIP workflows. (For the SCF JK / `density_fit()` behavior see the `--lowmem` row in the CLI options table.)
3. **ML(dft)/MM recombination** -- DFT replaces only `MLMMCore`'s high-level MODEL energy. `MLMMCore` evaluates REAL-low and MODEL-low with the selected MM backend and applies the subtractive expression. This workflow has no separate topology builder, MM calculator path, or DFT force evaluation.
4. **Population analysis & outputs** -- Mulliken, meta-Lowdin, and IAO charges and spin densities (UKS only) are written alongside the combined energy block in `result.yaml`.

## Outputs

```
out_dir/ (default: ./result_dft/)
├── ml_region_without_linkH.xyz # Exact ML selection before generated link-H
├── ml_region_with_linkH.xyz    # PySCF input snapshot after generated link-H
├── ml_region_without_linkH.pdb # PDB input with --convert-files; topology-bearing companion
├── ml_region_with_linkH.pdb    # PDB input with --convert-files; generated link-H as HL/LKH
├── result.yaml                 # DFT + ML(dft)/MM energy summary, charges, spin densities
├── result.json                 # only when --out-json is passed
└── (stdout)                    # Pretty-printed configuration blocks and energies
```

- `result.yaml` expands to:
  - `energy`: Hartree/kcal/mol values, convergence flag, wall time, backend info (`engine`: `gpu4pyscf(rks_lowmem)` / `gpu4pyscf` / `pyscf(cpu)`; `used_gpu`; `used_lowmem`).
  - `mlmm_energy`: REAL-low / MODEL-low MM evaluations and the recombined `E_total = E_REAL_low + E_ML(DFT) - E_MODEL_low` in Hartree and kcal/mol.
  - `charges [index, element, mulliken, lowdin, iao]`: Mulliken, meta-Lowdin, and IAO atomic charges (`null` when a method fails).
  - `spin_densities [index, element, mulliken, lowdin, iao]`: Mulliken, meta-Lowdin, and IAO spin densities (UKS-only for spins).
- It also summarizes charge, multiplicity, functional, basis, convergence knobs, and resolved output directory.

## CLI options

`mlmm dft --help` shows core options; `mlmm dft --help-advanced` shows the full option list. The full flag list is in the generated [command reference](reference/commands/index.md); the table below covers the options that need explanation.

| Option | Description | Default |
| --- | --- | --- |
| `-i, --input PATH` | Full enzyme structure (PDB/mmCIF, or XYZ with `--ref-pdb` topology). | Required |
| `--ref-pdb FILE` | Reference PDB topology when input is XYZ. | _None_ |
| `--parm PATH` | Amber parm7 topology for the full system. | Required |
| `--model-pdb PATH` | PDB defining the ML region (atom IDs must match the enzyme PDB). Optional when `--detect-layer` is enabled. | _None_ |
| `--model-indices TEXT` | Comma-separated atom indices for the ML region (ranges allowed, e.g. `1-5`). Used when `--model-pdb` is omitted. | _None_ |
| `--model-indices-one-based / --model-indices-zero-based` | Interpret `--model-indices` as 1-based or 0-based. | `True` (1-based) |
| `--detect-layer / --no-detect-layer` | Automatically detect ML/MM layers from input PDB B-factors (B=0/10/20). | Enabled |
| `-q, --charge INT` | Charge of the ML region. Required unless `-l/--ligand-charge` is given (PDB input or XYZ with `--ref-pdb`). | Required unless `-l/--ligand-charge` is provided |
| `-l, --ligand-charge TEXT` | Total charge or per-resname mapping (e.g. `SAM:1,GPP:-3`) used to derive the ML-region charge when `-q` is omitted (requires PDB input or `--ref-pdb`). | _None_ |
| `-m, --multiplicity INT` | Spin multiplicity (2S+1) for the ML region. | `1` |
| `--func-basis TEXT` | Functional/basis pair as `"FUNC/BASIS"`. | `wb97m-v/def2-tzvpd` |
| `--max-cycle INT` | SCF-iteration cap. | `100` |
| `--conv-tol FLOAT` | SCF convergence tolerance (Hartree). | `1e-9` |
| `--grid-level INT` | DFT integration grid level (0=coarse, 3=default, 5=fine, 9=very fine). | `3` |
| `--engine {gpu,cpu}` | Force GPU4PySCF (`gpu`) or CPU PySCF (`cpu`); `gpu` raises an error if GPU4PySCF is unavailable. | `gpu` |
| `--lowmem/--no-lowmem` | Use `gpu4pyscf.dft.rks_lowmem.RKS` for closed-shell GPU runs (memory-efficient direct JK; `mlmm dft` does not call `density_fit()` on either path). Open-shell, CPU, or pre-`rks_lowmem` GPU4PySCF auto-fall back to standard RKS/UKS. | `True` |
| `--embedcharge/--no-embedcharge` | Experimental direct PySCF electrostatic embedding of MM point charges. No xTB correction is used in `dft`. | `False` |
| `--embedcharge-cutoff FLOAT` | Include MM point charges within this distance of the ML region. | `12.0` Å |
| `-o, --out-dir DIR` | Output directory. | `./result_dft/` |
| `--config FILE` | Base YAML configuration file applied before explicit CLI options. | _None_ |
| `--show-config/--no-show-config` | Print resolved configuration and continue execution. | `False` |
| `--link-atom-method {scaled,fixed}` | Link-atom placement: `scaled` (g-factor, Gaussian ONIOM standard) or `fixed` (legacy 1.09 Å for C, 1.01 Å for N). | `scaled` |
| `--mm-backend {hessian_ff,openmm}` | MM backend for the low-level ONIOM evaluation. Hessians use finite differences by default; set `calc.mm_fd: false` for the `hessian_ff` analytical path. | `hessian_ff` |
| `--cmap/--no-cmap` | Preserve CMAP in both REAL and MODEL MM layers. | `--cmap` |
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
 real_parm7: real.parm7            # Amber parm7 topology
 model_pdb: ml_region.pdb          # ML-region definition
 embedcharge: false                # PySCF electrostatic embedding; no xTB in dft
 embedcharge_cutoff: 12.0          # MM point-charge cutoff from ML region (Å)
dft:
 func_basis: wb97m-v/def2-tzvpd      # exchange-correlation functional / basis set
 conv_tol: 1.0e-09                # SCF convergence tolerance (Hartree)
 max_cycle: 100                    # SCF iteration cap
 grid_level: 3                     # PySCF grid level
 verbose: 0                        # PySCF verbosity (0-9); CLI -v 2/3 raises runtime PySCF verbosity to >=4
 out_dir: ./result_dft/            # output directory root
```

Full schema (every key and default): [YAML Reference](yaml-reference.md).

## Notes

- A matching def2 effective core potential is auto-attached whenever the basis name begins with `def2` (no element-presence check).
- **Blackwell-architecture GPUs** (RTX 50xx): verify that the installed
  GPU4PySCF/CuPy stack supports the device. If the GPU path fails, use
  `--engine cpu` or an external DFT program.
- **Out-of-memory with def2-TZVPD**: memory depends on atom types, basis,
  functional, grid, and software stack. Pilot the target system and, if
  necessary, choose a smaller basis only after validating its effect on the
  quantities of interest.
- Compiled GPU4PySCF wheels may not support non-x86 systems; build from source in that case (see https://github.com/pyscf/gpu4pyscf).

## See Also

- [Common Error Recipes](recipes-common-errors.md) — Symptom-first failure routing
- [Troubleshooting](troubleshooting.md) — Detailed troubleshooting guide
- [freq](freq.md) — Vibrational frequency analysis (often precedes DFT single-point evaluation)
- [opt](opt.md) — Single-structure geometry optimization
- [all](all.md) — End-to-end workflow with `--dft`
- [YAML Reference](yaml-reference.md) — Full `dft` configuration options
- [Glossary](glossary.md) — Definitions of DFT, SP (Single Point)
