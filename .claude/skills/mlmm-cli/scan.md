# `mlmm scan`

## Purpose

1D distance scan with harmonic restraints. Drives one bond between two
atoms from its current value to a target, generating a sequence of
relaxed geometries. Used standalone for surface exploration or as a
building block when you need to seed `path-search` manually.

For most workflows, prefer `mlmm all --scan-lists '...'`
(see `all-scan-list.md`); the standalone `scan` is for one-off
explorations.

## Synopsis

```bash
mlmm scan -i input.pdb \
    -c 'NAME RESNAME RESID' 'NAME RESNAME RESID' --target 1.6 \
    [-q / -l / -m] [-b uma] [-o ./result_scan/]
```

## Key flags

| flag | type | default | description |
|---|---|---|---|
| `-i, --input` | path | required | Reactant geometry |
| `-c, --bond` | str×2 | required | Two atom specs forming the bond to scan |
| `--target` | float | required | Target distance (Å) |
| `--n-steps` | int | (live default) | Number of intermediate snapshots |
| `--force-constant` | float | (live default) | Restraint stiffness (Hartree/Bohr²); check `BIAS_KW` |
| `-q` / `-l` / `-m` | — | — | Charge / spin |
| `-b, --backend` | str | `uma` | MLIP backend |
| `-o, --out-dir` | path | `./result_scan/` | Output directory |

## Examples

```bash
mlmm scan -i 1.R.pdb \
    -c 'CS1 SAM 320' 'C7 GPP 321' --target 1.6 \
    -l 'SAM:1,GPP:-3' -b uma -o result_scan
```

## Output

```
result_scan/
├── result.json
├── scan_NN.{xyz,pdb}      # per-step relaxed geometries
├── mep.xyz                # stitched scan trajectory
└── scan.log
```

`result.json` lists per-step energy + final distance; useful for
plotting with `trj2fig.md`.

## See also

- `scan2d.md`, `scan3d.md` — 2D / 3D analogs.
- `all-scan-list.md` — staged multi-bond scans inside the full pipeline.
- Defaults: `import mlmm.defaults as d; print(d.BIAS_KW, d.BOND_KW)`
