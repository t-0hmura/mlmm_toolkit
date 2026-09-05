# Quickstart: `mlmm scan`

## Goal

Generate restrained endpoint structures and trajectories from a single
structure using a YAML scan specification.

## Prerequisites

- Full-system structure (`-i`): `pocket.pdb`, with atom identity and order matching `real.parm7`
- MM topology (`--parm`): `real.parm7`
- ML subset (`--model-pdb`): `ml_region.pdb`, without link hydrogens. Explicit model indices or valid B-factor layers are also accepted.

One YAML `stages` entry defines one stage. Multiple distance tuples within an
entry are advanced concertedly; multiple entries form a sequential multistage scan.
Use `scan2d` when two distances must instead form independent grid axes.

## 1. Prepare `scan.yaml`

```yaml
one_based: true
stages:
 - [[12, 45, 2.20]]
 - [[10, 55, 1.35], [23, 34, 1.80]]
```

## 2. Run scan

```bash
mlmm scan -i pocket.pdb --parm real.parm7 --model-pdb ml_region.pdb \
 -q 0 -s scan.yaml -o ./result_scan
```

```{note}
To validate the spec without running (GPU-free), add `--print-parsed`. This prints the parsed targets and exits before any calculation, so it does **not** produce the scan outputs listed below.
```

## Output validation

- `result_scan/stage_01/result.pdb`
- `result_scan/stage_02/result.pdb`
- `result_scan/scan_trj.xyz` (always written); `result_scan/scan.pdb` when conversion and a reference topology are available

## Inline literal input (without YAML file)

Instead of a YAML spec file, you can pass scan targets directly on the command line:

```bash
mlmm scan -i layered.pdb --parm system.parm7 -q 0 \
  --scan-lists '[(1,5,1.4)]' --no-preopt --no-endopt
```

Or using PDB atom selectors:

```bash
mlmm scan -i layered.pdb --parm system.parm7 -q 0 \
  --scan-lists '[("TYR,285,CA","MMT,309,C10",2.20)]' --no-preopt --no-endopt
```

Both 1-based atom indices and PDB atom name strings are accepted. See [scan.md](scan.md) for full details.

For detailed options, run `mlmm scan --help-advanced`.

## Next step

- Feed scan results to path refinement with [all](all.md) or [path-search](path-search.md).
