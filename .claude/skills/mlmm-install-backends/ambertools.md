# AmberTools (ambertools.md)

`mlmm-toolkit` uses **AmberTools** (specifically `tleap` and
`antechamber`) for the MM region parameterization step
(`mlmm mm-parm`). AmberTools is required if you want to drive the
parameter file (`parm7`) generation through the toolkit; if you have
hand-built `parm7` / `rst7` files, you can skip the install and feed
them directly to `mlmm extract` / `mlmm define-layer`.

## Install via conda (recommended)

```bash
conda activate <your_mlmm_env>
conda install -c conda-forge ambertools
```

This pulls `tleap`, `antechamber`, `parmchk2`, `cpptraj`, and the
standard force-field parameter sets (Amber14SB, GAFF2, etc.).

Verify:

```bash
which tleap antechamber parmchk2
tleap -h | head -3
```

## Install from the official tarball

If `conda-forge` is not available, AmberTools can also be installed
from the official tarball at https://ambermd.org/AmberTools.php
(free academic license required). After install:

```bash
source <amber_install>/amber.sh   # sets AMBERHOME, PATH, LD_LIBRARY_PATH
which tleap
```

Add `source <amber_install>/amber.sh` to every PBS / SLURM job that
calls `mlmm mm-parm`, or build it into your env-init shell hook.

## CLI usage (`mlmm mm-parm`)

```bash
mlmm mm-parm -i complex.pdb \
    --ligand 'SAM:1,GPP:-3' \
    --force-field amber14sb \
    --water tip3p \
    -o complex.parm7
```

`mm-parm` writes both `<basename>.parm7` and `<basename>.rst7`. Pass
both to downstream subcommands via the toolkit's standard flags.

## Common gotchas

| Symptom | Cause / fix |
|---|---|
| `tleap: command not found` | AmberTools not installed; or `<amber_install>/amber.sh` not sourced. |
| `Unknown residue name 'GPP'` | The ligand has no parameters; run `antechamber` to derive GAFF parameters first, then point `mm-parm` at the resulting `.frcmod`. |
| `parm7` written but `mm-parm` reports charge mismatch | Check `mlmm-structure-io/charge-multiplicity.md` — the `-l` mapping must agree with the protonation states in the PDB. |
| Aarch64 wheels not available | AmberTools provides aarch64 conda packages on conda-forge — `conda install -c conda-forge ambertools` works on ARM machines. |

## What AmberTools is **not** used for

- Running MD inside `mlmm-toolkit`. The toolkit only uses tleap /
  antechamber for parameter generation; the actual ML/MM gradient
  evaluation goes through bundled `hessian_ff`.
- Producing trajectories. The toolkit reads positions from the
  `parm7` + `rst7` pair (or from a fresh PDB) but does not call
  `sander` / `pmemd`.

## See also

- `core.md` — `mlmm-toolkit` install (the toolkit itself does not depend on AmberTools, only `mm-parm` does).
- `mlmm-cli/mm-parm.md` — full `mm-parm` flag reference.
- `mlmm-structure-io/parm7.md` — `parm7`/`rst7` file structure and editing.
- `mlmm-structure-io/charge-multiplicity.md` — when ligand charges
  conflict between PDB and the AmberTools parameter set.