# `scan`

Drive a reaction coordinate on a layered enzyme PDB by scanning bond distances with harmonic restraints using the ML/MM calculator. `mlmm scan` performs a staged, bond-length-driven scan using the ML/MM calculator (`mlmm.backends.mlmm_calc.mlmm`) with harmonic restraints: at each step the temporary targets are updated, restraint wells are applied, and the structure is relaxed with LBFGS. The ML/MM calculator couples an MLIP backend (selected via `-b/--backend`; default: UMA) and hessian_ff. Use `-s/--scan-lists` to define targets as a YAML/JSON spec file (recommended) or as inline Python literals.

## When to use

- Use when generating a coarse reaction trajectory from a single starting structure by driving one or more interatomic distances toward target values; provides intermediate/product candidates for downstream MEP refinement.

## Quick examples

```bash
mlmm scan -i pocket.pdb --parm real.parm7 --model-pdb ml_region.pdb \
 -q 0 -s scan.yaml -o ./result_scan
```
(Add `--print-parsed` to validate the parsed scan spec and exit without running the GPU calculation.)

```bash
# Inline Python literal
mlmm scan -i pocket.pdb --parm real.parm7 --model-pdb ml_region.pdb \
 -q 0 -s "[(12,45,2.20)]"
```

```bash
# Dump trajectories for stage-by-stage inspection
mlmm scan -i pocket.pdb --parm real.parm7 --model-pdb ml_region.pdb \
 -q 0 -s scan.yaml --dump -o ./result_scan_dump
```

## Inputs

Command form:

```bash
mlmm scan -i INPUT.pdb --parm real.parm7 --model-pdb ml_region.pdb \
 -q CHARGE [-m MULT] \
 [-s scan.yaml | -s "[(I,J,TARGET_ANG)]"] [options]
```

| Input | Required | Notes |
| --- | --- | --- |
| `-i, --input` | yes | Input PDB (or XYZ with `--ref-pdb` for topology). |
| `--parm` | yes | Amber prmtop for the full REAL system. |
| `--model-pdb` | optional | PDB defining the ML region; optional when `--detect-layer` is enabled or `--model-indices` is provided. |
| `-q, --charge` | yes (unless `-l`) | Net ML-region charge. |
| `-s, --scan-lists` | yes | Scan targets: a YAML/JSON spec file path (auto-detected) or inline Python literal(s). |

### Input syntax

**YAML/JSON spec format (recommended)**

`-s/--scan-lists` auto-detects YAML/JSON files. Pass a file path to use the spec format:

```yaml
one_based: true # optional; defaults to CLI --one-based/--zero-based
stages:
 - [[12, 45, 2.20]]
 - [[10, 55, 1.35], [23, 34, 1.80]]
```

- `stages` is required.
- Each stage is a list of `(i, j, target_A)` triples.
- Indices may be integers or PDB selectors (for PDB input), same as inline literals.

**Inline literal format**

When `-s/--scan-lists` receives a value that is not a file path, it is treated as a **Python literal** string evaluated by the CLI. Shell quoting matters.

Each literal is a Python list of triples `(atom1, atom2, target_A)`:

```
-s '[(atom1, atom2, target_A),...]'
```

- Wrap the entire literal in **single quotes** so the shell does not interpret parentheses or spaces.
- Each triple drives the distance between `atom1`--`atom2` toward `target_A`.
- One literal = one **stage**. For multiple stages, pass multiple literals after a **single** `-s/--scan-lists` flag (do not repeat the flag).

Atoms can be given as **integer indices** or **PDB selector strings**:

| Method | Example | Notes |
| --- | --- | --- |
| Integer index | `(1, 5, 2.0)` | 1-based by default (`--one-based`) |
| PDB selector | `("TYR,285,CA", "MMT,309,C10", 2.0)` | Residue name, residue number, atom name |

PDB selector tokens can be separated by any of: comma `,`, space, slash `/`, backtick `` ` ``, or backslash `\`. Token order is flexible.

```bash
# All of these specify the same atom:
"TYR,285,CA"
"TYR 285 CA"
"TYR/285/CA"
"285,TYR,CA" # order is flexible
```

Quoting rules:

```bash
# Correct: single-quote the list, double-quote selector strings inside
-s '[("TYR,285,CA","MMT,309,C10",1.35)]'

# Correct: integer indices need no inner quotes
-s '[(1, 5, 2.0)]'

# Avoid: double-quoting the outer literal requires escaping inner quotes
-s "[(\"TYR,285,CA\",\"MMT,309,C10\",1.35)]"
```

Pass multiple literals after a single `-s/--scan-lists` flag. Each literal becomes one stage:

```bash
# Stage 1: drive one bond to 1.35 Å
# Stage 2: drive two bonds simultaneously
-s \
 '[("TYR,285,CA","MMT,309,C10",1.35)]' \
 '[("TYR,285,CA","MMT,309,C10",2.20),("TYR,285,CB","MMT,309,C11",1.80)]'
```

Stages run sequentially; each starts from the previous stage's relaxed result. **Do not repeat the `-s/--scan-lists` flag** -- supply all stage literals after a single flag.

**Bidirectional scan (4-tuple)**

Instead of a 3-tuple `(i, j, target)`, you can pass a **4-tuple** `(i, j, start, end)` to scan in both directions from the current geometry. The CLI automatically expands each 4-tuple into two stages:

1. **Pass 1:** Drive `i`--`j` from the current distance toward `start`.
2. **Pass 2:** Restore the initial geometry and drive `i`--`j` toward `end`.

The concatenated trajectory is assembled as `start → initial → end`, giving a continuous path through the starting structure.

```bash
# Bidirectional scan: drive bond 12--45 from current geometry
# toward 1.35 Å (pass 1) and toward 2.50 Å (pass 2)
mlmm scan -i pocket.pdb --parm real.parm7 --model-pdb ml_region.pdb \
 -q 0 -s '[(12, 45, 1.35, 2.50)]'
```

This is equivalent to two manual stages with a geometry reset between them, but avoids the need to script it yourself. Mixed 3-tuples and 4-tuples are accepted in the same literal.

## Workflow

1. Load the structure through `geom_loader`, resolving charge/spin from the CLI
    or defaults. Provide `--parm`, `--model-pdb`, `-q/--charge`, and optionally
    `-m/--multiplicity` for the ML/MM calculator.
2. Optionally run an unbiased preoptimization (`--preopt`) before any
    biasing so the starting point is relaxed.
3. Parse stage targets from `-s/--scan-lists` (YAML/JSON spec file or inline literal), then normalize the
    `(i, j)` indices (1-based by default). When the input is a PDB, each entry
    may be either an integer index or an atom selector string like `'TYR,285,CA'`;
    selector fields can be separated by spaces, commas, slashes, backticks, or
    backslashes and may be in any order.
4. Compute the per-bond displacement and split into steps:
 - For scan tuples `[(i, j, target_A)]`, compute `delta = target - current_distance_A`.
 - With `--max-step-size = h`, the stage takes `N = ceil(max(|delta|) / h)` biased relaxations.
 - Each pair's incremental change is `delta_k = delta_k / N` (Å). At step `s`, the temporary
  target is `r_k(s) = r_k(0) + s * delta_k`.
5. March through all steps, applying the harmonic wells
    `E_bias = sum 1/2 * k * (|r_i - r_j| - target_k)^2` and minimizing with LBFGS.
    `k` comes from `--bias-k` (eV/Å²) and is converted once to Hartree/Bohr^2.
    Coordinates are stored in Bohr for PySisyphus and converted internally for reporting.
6. After the last step of each stage, optionally run an unbiased relaxation
    (`--endopt`) before reporting covalent bond changes and writing the
    `result.*` files.
7. Repeat for every stage; optional trajectories are dumped only when `--dump`
    is `True`.

## Outputs

Each stage writes its final geometry and biased-step trajectory under `stage_XX/`, with a combined trajectory at the root. The files you check first are the per-stage `result.pdb` (or `result.xyz`) and the always-generated `scan_trj.xyz` / `scan.pdb`.

```
out_dir/ (default: ./result_scan/)
├─ scan_trj.xyz              # Combined trajectory across all stages (always written)
├─ scan.pdb                  # Combined PDB companion (PDB inputs only; always written)
├─ preopt/                   # Present when --preopt is True
│  ├─ result.xyz
│  └─ result.pdb             # Only for PDB inputs
└─ stage_XX/                 # One folder per stage (k = 01..K)
   ├─ result.xyz             # Final (possibly endopt) geometry
   ├─ result.pdb             # If input was PDB
   ├─ scan_trj.xyz           # Per-stage biased step frames (always written)
   └─ scan.pdb               # PDB version of scan_trj.xyz (PDB inputs only; always written)
```

## CLI options

The full flag list is in the generated [command reference](reference/commands/index.md); the table below covers the options that need explanation.

| Option | Description | Default |
| --- | --- | --- |
| `-i, --input PATH` | Input PDB (or XYZ with `--ref-pdb` for topology). | Required |
| `--parm PATH` | Amber prmtop for the full REAL system. | Required |
| `--model-pdb PATH` | PDB defining the ML region (atom IDs). Optional when `--detect-layer` is enabled or `--model-indices` is provided. | _None_ |
| `--model-indices TEXT` | Comma-separated ML-region atom indices (ranges allowed). | _None_ |
| `--model-indices-one-based / --model-indices-zero-based` | Interpret `--model-indices` as 1-based or 0-based. | `True` (1-based) |
| `--detect-layer / --no-detect-layer` | Detect ML/MM layers from input PDB B-factors. | `True` |
| `-q, --charge INT` | Net ML-region charge. | _None_ (required unless `-l` is given) |
| `-l, --ligand-charge TEXT` | Per-resname charge mapping (e.g., `GPP:-3,SAM:1`). Derives net charge when `-q` is omitted. | _None_ |
| `-m, --multiplicity INT` | Spin multiplicity (2S+1). | `1` |
| `--freeze-atoms TEXT` | Comma-separated 1-based atom indices to freeze (merged with YAML `geom.freeze_atoms`). | _None_ |
| `--hess-cutoff FLOAT` | Distance cutoff (Å) from ML region for MM atoms to include in Hessian calculation. Can be combined with `--detect-layer`. | _None_ |
| `--movable-cutoff FLOAT` | Movable-MM distance cutoff (Å); providing this disables `--detect-layer`. | _None_ |
| `-s, --scan-lists TEXT` | Scan targets: a YAML/JSON spec file path (auto-detected) or inline Python literal(s) with `(i, j, target_A)` triples or `(i, j, start, end)` 4-tuples for bidirectional scans. Each literal is one stage; supply multiple literals after a single flag. `i`/`j` can be integer indices or PDB atom selectors like `"TYR,285,CA"`. | Required |
| `--one-based/--zero-based` | Interpret atom indices as 1-based (default) or 0-based. | `True` (1-based) |
| `--print-parsed/--no-print-parsed` | Print parsed stage tuples after `-s/--scan-lists` resolution. | `False` |
| `--max-step-size FLOAT` | Maximum change in any scanned bond per step (Å). Controls the number of biased relaxation steps. | `0.20` |
| `--bias-k FLOAT` | Harmonic bias strength `k` in eV/Å². | `300` |
| `--opt-mode {grad,hess,lbfgs,rfo,light,heavy}` | Compatibility option for `mlmm all` forwarding. Current scan relaxations use LBFGS regardless of mode. | _None_ |
| `--max-cycles INT` | Maximum LBFGS cycles per biased step and per pre/end optimization stage. | `10000` |
| `--relax-max-cycles INT` | Compatibility alias of `--max-cycles` (overrides it when provided). | _None_ |
| `--preopt/--no-preopt` | Run an unbiased optimization before scanning. | `False` |
| `--endopt/--no-endopt` | Run an unbiased optimization after each stage. | `False` |
| `--dump/--no-dump` | Dump per-step optimizer trajectory files. Note: `scan_trj.xyz`/`scan.pdb` are always written regardless of this flag. | `False` |
| `-o, --out-dir TEXT` | Output directory root. | `./result_scan/` |
| `--thresh TEXT` | Convergence preset (`gau_loose\|gau\|gau_tight\|gau_vtight\|baker\|never`). | _None_ (inherits `gau`) |
| `--config FILE` | Base YAML configuration file (applied first). | _None_ |
| `--ref-pdb FILE` | Reference PDB topology when `--input` is XYZ. | _None_ |
| `-b, --backend CHOICE` | MLIP backend for the ML region: `uma`, `orb`, `mace`, `aimnet2`. | `uma` |
| `--embedcharge/--no-embedcharge` | Enable xTB point-charge embedding correction for MM-to-ML environmental effects. | `False` |
| `--embedcharge-cutoff FLOAT` | Cutoff radius (Å) for embed-charge MM atoms. | `12.0` |
| `--cmap/--no-cmap` | Enable CMAP (backbone cross-map dihedral correction) in model parm7. Default: disabled (consistent with Gaussian ONIOM). | `--no-cmap` |
| `--mm-backend [hessian_ff\|openmm]` | MM backend (analytical Hessian vs OpenMM finite-difference). | `hessian_ff` |
| `--link-atom-method [scaled\|fixed]` | Link-atom placement: scaled ($g$-factor) or fixed 1.09/1.01 Å. | `scaled` |
| `--out-json/--no-out-json` | Write machine-readable `result.json` to `out_dir`. | `False` |
| `--dry-run/--no-dry-run` | Validate options and print the execution plan without running the scan. Shown in `--help-advanced`. | `False` |
| `--convert-files/--no-convert-files` | Toggle XYZ/TRJ to PDB companions when a PDB template is available. | `True` |

## YAML configuration

The scan reads the shared `geom` (`coord_type`, `freeze_atoms`), `calc` / `mlmm` (ML/MM calculator setup), and `opt` / `lbfgs` (optimizer) sections, plus `bias` (`k`, harmonic strength in eV/Å²) and a `bond` section for MLIP-based bond-change detection.

Full schema (every key and default): [YAML Reference](yaml-reference.md).

## See Also

- [Common Error Recipes](recipes-common-errors.md) — Symptom-first failure routing
- [Troubleshooting](troubleshooting.md) — Detailed troubleshooting guide
- [scan2d](scan2d.md) — 2D distance grid scan
- [scan3d](scan3d.md) — 3D distance grid scan
- [opt](opt.md) — Single-structure geometry optimization
- [all](all.md) — End-to-end workflow with `--scan-lists` for single-structure inputs
- [path-search](path-search.md) — MEP search using scan endpoints as intermediates
