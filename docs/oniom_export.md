# `oniom-gaussian` / `oniom-orca`

## Overview

> **Summary:** Export an ML/MM system to Gaussian or ORCA ONIOM input format. The force field parameters are extracted directly from an Amber parm7 topology file, supporting any Amber-compatible force field (ff14SB, ff19SB, GAFF2, etc.).

### Quick reference
- **Two subcommands:** `mlmm oniom-gaussian` (Gaussian ONIOM `.com`/`.gjf`) and `mlmm oniom-orca` (ORCA QM/MM `.inp`).
- **Input:** Amber parm7 topology, optional coordinate file (PDB/XYZ), optional ML-region PDB.
- **Element validation:** Best-effort element sequence check between the parm7 and coordinate file (atom order is assumed invariant).
- **Boundary handling:** QM/MM covalent boundaries are detected from topology bonds and exported as link atoms (Gaussian) / auto-capping boundaries (ORCA).
- **Use when:** You want to run QM/MM calculations in Gaussian or ORCA using the same system prepared with mlmm_toolkit.

## Usage

### oniom-gaussian
```bash
mlmm oniom-gaussian --parm7 real.parm7 [-i coords.pdb] [--model-pdb ml_region.pdb] \
 -o output.com [-q CHARGE] [-m MULT] [--method "B3LYP/6-31G(d,p)"] \
 [--near FLOAT] [--nproc INT] [--mem TEXT] \
 [--element-check|--no-element-check]
```

### oniom-orca
```bash
mlmm oniom-orca --parm7 real.parm7 [-i coords.pdb] [--model-pdb ml_region.pdb] \
 -o output.inp [-q CHARGE] [-m MULT] [--method "B3LYP D3BJ def2-SVP"] \
 [--total-charge INT] [--total-mult INT] \
 [--nproc INT] [--near FLOAT] [--orcaff PATH] \
 [--convert-orcaff|--no-convert-orcaff] [--element-check|--no-element-check]
```

### Examples
```bash
# Generate Gaussian ONIOM input
mlmm oniom-gaussian --parm7 real.parm7 -i pocket.pdb --model-pdb ml_region.pdb \
 -o system.com -q 0 -m 1

# Generate ORCA QM/MM input
mlmm oniom-orca --parm7 real.parm7 -i pocket.pdb --model-pdb ml_region.pdb \
 -o system.inp -q 0 -m 1
```

## Description

Both subcommands read an Amber parm7 topology file and generate QM/MM input files for external programs. The workflow is:

1. **Load topology**: Read the parm7 file via ParmEd to extract atom types, charges, bonds, angles, dihedrals, and van der Waals parameters.
2. **Load coordinates** (optional): When `-i/--input` is provided, read a PDB or XYZ coordinate file. If `--element-check` is enabled (default), validate that the element sequence matches the parm7 topology.
3. **Identify QM region**: When `--model-pdb` is supplied, match atoms from the ML-region PDB to the full system to define the QM (high) layer.
4. **Generate input**: Write the appropriate input file:
 - **Gaussian**: ONIOM input with `%nproc`, `%mem`, method specification, atom coordinates with layer assignments, connectivity, and automatic link-atom annotations at QM/MM covalent boundaries. Movable atoms within `--near` angstroms of the QM region are identified.
 - **ORCA**: QM/MM input with `%pal nprocs`, method specification, QM atom selection, `Charge_Total`/`Mult_Total`, and a reference to the ORCAFF force-field file. ORCA handles link-atom placement from QM/MM boundary connectivity.

## CLI options

### oniom-gaussian
| Option | Description | Default |
| --- | --- | --- |
| `--parm7 PATH` | Amber parm7 topology file. | Required |
| `-i, --input PATH` | Coordinate file (PDB/XYZ); atom order must match parm7. | _None_ |
| `--element-check / --no-element-check` | Validate element sequence between input and parm7. | `True` |
| `--model-pdb PATH` | PDB file defining QM region atoms. | _None_ |
| `-o, --output PATH` | Output Gaussian input file (`.com` or `.gjf`). | Required |
| `--method TEXT` | QM method and basis set. | `B3LYP/6-31G(d,p)` |
| `-q, --charge INT` | Charge of QM region. | `0` |
| `-m, --mult INT` | Multiplicity of QM region. | `1` |
| `--near FLOAT` | Distance cutoff for movable atoms (angstrom). | `6.0` |
| `--nproc INT` | Number of processors. | `8` |
| `--mem TEXT` | Memory allocation. | `16GB` |

### oniom-orca
| Option | Description | Default |
| --- | --- | --- |
| `--parm7 PATH` | Amber parm7 topology file. | Required |
| `-i, --input PATH` | Coordinate file (PDB/XYZ); atom order must match parm7. | _None_ |
| `--element-check / --no-element-check` | Validate element sequence between input and parm7. | `True` |
| `--model-pdb PATH` | PDB file defining QM region atoms. | _None_ |
| `-o, --output PATH` | Output ORCA input file (`.inp`). | Required |
| `--method TEXT` | QM method and basis set. | `B3LYP D3BJ def2-SVP` |
| `-q, --charge INT` | Charge of QM region. | `0` |
| `-m, --mult INT` | Multiplicity of QM region. | `1` |
| `--total-charge INT` | Total charge of full QM+MM system (`Charge_Total`). | topology-derived |
| `--total-mult INT` | Total multiplicity of full QM+MM system (`Mult_Total`). | same as `--mult` |
| `--nproc INT` | Number of processors. | `8` |
| `--near FLOAT` | Distance cutoff (Å) used to define ActiveAtoms when no layer B-factors are present. | `6.0` |
| `--orcaff PATH` | Path to ORCAFF.prms; if omitted, exporter uses/creates `<parm7_stem>.ORCAFF.prms`. | auto |
| `--convert-orcaff/--no-convert-orcaff` | Auto-run `orca_mm -convff -AMBER` when ORCAFF.prms is missing. | `True` |

## Outputs
```
<output>.com / <output>.gjf # Gaussian ONIOM input (oniom-gaussian)
<output>.inp # ORCA QM/MM input (oniom-orca)
```
- Console summary with QM atom count and (for ORCA) a reminder to copy the parm7 file alongside the input.

## Notes
- Requires `ParmEd` (`parmed`) for topology parsing; an `ImportError` is raised if it is not installed.
- Element validation between the coordinate file and parm7 is best-effort; atom ordering is assumed to be invariant.
- For Gaussian, the `--near` cutoff determines which MM atoms are allowed to move during the ONIOM optimization.
- Gaussian `NonBon` is written with Amber-standard 1-4 electrostatic scaling (`-1.2`).

## Quick failure recovery
1. `ImportError: parmed` when starting export.

```bash
python -c "import parmed; print(parmed.__version__)"
pip install parmed
```

2. Element/order validation fails (`--element-check`).

```bash
mlmm oniom-gaussian --parm7 real.parm7 -i pocket.pdb --model-pdb ml_region.pdb \
 -o system.com --element-check
```

3. QM/MM boundary or layer assignment looks wrong.

```bash
mlmm define-layer -i pocket.pdb --real-parm7 real.parm7 -o layered.pdb
mlmm oniom-orca --parm7 real.parm7 -i layered.pdb --model-pdb ml_region.pdb \
 -o system.inp
```

---

## See Also

- [Common Error Recipes](recipes_common_errors.md) -- Symptom-first failure routing
- [Troubleshooting](troubleshooting.md) -- Detailed troubleshooting guide

- [mm_parm](mm_parm.md) -- Build AMBER topology (parm7/rst7)
- [define_layer](define_layer.md) -- Define ML/MM layers
- [opt](opt.md) -- ML/MM geometry optimization
