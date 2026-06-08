# CLI Conventions

Conventions shared across all `mlmm` commands.

## Boolean options

Every boolean CLI flag accepts **all four forms**:

| Form | Example |
|---|---|
| Positive flag | `--tsopt` |
| Negative flag | `--no-tsopt` |
| Positive value | `--tsopt True`, `--tsopt yes`, `--tsopt 1`, `--tsopt on` |
| Negative value | `--tsopt False`, `--tsopt no`, `--tsopt 0`, `--tsopt off` |

```bash
--tsopt --thermo --no-dft                        # toggle form
--tsopt True --thermo yes --dft 0                # value form
--tsopt true --thermo on --dft off               # mix is fine
```

All four forms route through a single root-CLI `bool_compat` synthesizer; the toggle form (`--tsopt` / `--no-tsopt`) is canonical, and the value form (`--tsopt True`) is accepted as a legacy alias for backward compatibility. `tests/test_bool_compat_cli.py` walks every registered bool option against every form on every release, so a missing entry is caught by CI.

Common toggles: `--tsopt` / `--thermo` / `--dft` (post-processing stages) · `--dump` (write trajectory files) · `--preopt` / `--endopt` (pre/post optimisation) · `--climb` (climbing-image MEP).

### Contributing a new bool flag

When adding a boolean flag inside a subcommand, always route it through one of the `add_*_option()` factories in `mlmm/cli/common_options.py` and register the long name in the matching `_COMMAND_BOOL_*_OPTIONS` table in `mlmm/cli/app.py`. Avoid writing `@click.option("--foo/--no-foo", ...)` or `type=click.BOOL` directly in the subcommand body — that bypasses the registry, falls out of test coverage, and silently drops the value form.

## Progressive help

```bash
mlmm <subcmd> --help                # core options
mlmm <subcmd> --help-advanced       # full option set
```

`all`, every calc subcommand (`opt`, `tsopt`, `freq`, `irc`, `dft`, `scan` / `scan2d` / `scan3d`, `path-opt`, `path-search`), and the main utilities (`mm-parm`, `define-layer`, `add-elem-info`, `trj2fig`, `energy-diagram`, `oniom-export`, `extract`, `fix-altloc`) follow this pattern.

(verbosity-levels)=

## Verbosity levels

`-v/--verbose LEVEL` is an integer from 0 to 3 (**default 2**) that sets how much each command prints to the console. It is a per-command option, so write it with the subcommand, e.g. `mlmm opt -v 1 ...`. The same four levels apply to every command; command pages describe only their command-specific payload (e.g. the `opt` cycle table or the `freq` thermochemistry summary).

| Level | What you see |
|---|---|
| `-v 0` | Silent. Confirm success from the exit code and the output artifacts. |
| `-v 1` | Milestones only: version, input summary, key settings, output location, dry-run / final status. No banner, `[command]`, `[mode]`, or config dump. |
| `-v 2` | Default. Adds the banner, `[command]`, `[mode]`, stage progress, the main optimizer cycle table, terminal status, the one-line Hessian summary, thermo / DFT summaries, and elapsed time. |
| `-v 3` | Debug: resolved config / dry-run plan, backend DEBUG, raw optimizer and internal-coordinate chatter, `[HessianTiming]`, and `[HessianVRAM]`. |

A semantic failure is a failure at any level: a `Traceback` that appears only at `-v 3` still means the run failed.

## ML/MM required options

Per-stage subcommands (everything except `all`, `extract`, `mm-parm`, `define-layer`) need:

```bash
--parm real.parm7              # Amber parm7 topology of the full (real) system
--model-pdb model.pdb          # ML region (model) PDB
```

`mlmm all` generates both automatically; standalone subcommands require them explicitly.

```bash
mlmm path-search -i R.pdb P.pdb --parm real.parm7 --model-pdb model.pdb -q 0 -m 1
```

## Inspect resolved config

```bash
mlmm opt -i input.pdb --parm real.parm7 -q -1 --show-config --dry-run
```

## B-factor layer encoding

| Layer | B-factor | Meaning |
|---|---|---|
| ML | 0.0 | MLIP energy / force / Hessian |
| Movable-MM | 10.0 | MM atoms free to move |
| Frozen | 20.0 | Coordinates fixed (non-bonded interactions still included) |

Tolerance ±1.0 when reading B-factors. Inspect visually by colouring on B-factor. Four ways to assign layers:

```bash
# 1. define-layer subcommand (recommended)
mlmm define-layer -i system.pdb --model-pdb ml_region.pdb -o labeled.pdb
```

```yaml
# 2. Distance cutoffs (YAML)
calc:
  hess_cutoff: 3.6        # Hessian-target MM
  movable_cutoff: 8.0     # Movable-MM (beyond → Frozen)

# 3. Read existing B-factors
calc:
  use_bfactor_layers: true

# 4. Explicit index lists
calc:
  hess_mm_atoms:    [100, 101, 102, ...]
  movable_mm_atoms: [200, 201, 202, ...]
  frozen_mm_atoms:  [300, 301, 302, ...]
```

## Residue selectors

| Form | Example | Notes |
|---|---|---|
| By residue name | `-c 'SAM,GPP'` / `-c 'LIG'` | If multiple residues share a name, **all** matches are included (warning logged). |
| By residue ID | `-c '123,456'` / `-c 'A:123,B:456'` / `-c '123A'` / `-c 'A:123A'` | Optional chain prefix; trailing letter = insertion code. |
| By PDB file | `-c substrate.pdb` | Use coordinates from a separate PDB to locate substrates. |

## Charge specification

For PDB inputs, `--ligand-charge` lets you specify charges only for non-standard residues (substrates, cofactors, metal ions). The net system charge is **automatically derived** by summing standard amino-acid charges, ions, and your ligand charges.

```bash
-l 'SAM:1,GPP:-3'              # per-residue mapping (recommended)
-l 'SAM=1,GPP=-3'              # `=` separator also accepted
-l -3                           # single integer = total ligand charge
-q 0                            # explicit total system charge override
```

Residue-name matching is case-insensitive; unmapped non-standard residues default to charge 0.

**Resolution order** (highest priority first):

1. `-q/--charge` explicit CLI override.
2. ML-region determination summary (sum of amino acids, ions, `--ligand-charge`; only when `-c` is set and extraction runs).
3. `--ligand-charge` fallback when extraction is skipped (PDB input or `--ref-pdb` required).
4. `.gjf` template metadata (Gaussian-style charge / spin header; applies to `opt`, `tsopt`, `freq`, `irc`, `scan`, `dft`, `oniom-export` when the input is `.gjf`).
5. Default: abort if unresolved.

```{tip}
Always provide `--ligand-charge` for non-standard residues to ensure correct charge propagation.
```

## Spin multiplicity

```bash
-m 1    # singlet (default)
-m 2    # doublet
-m 3    # triplet
```

`mm-parm` uses the separate `--ligand-mult` for per-residue multiplicity metadata.

## Atom selectors

```bash
--scan-lists '[(1, 5, 2.0)]'                                          # 1-based integer indices
--scan-lists '[("TYR,285,CA", "MMT,309,C10", 2.20)]'                  # PDB-style selector strings
```

Selector field delimiters: space · comma · slash · backtick · backslash — e.g. `'TYR 285 CA'`, `'TYR,285,CA'`, `'TYR/285/CA'`, `` 'TYR`285`CA' ``, `'TYR\285\CA'`. The three tokens (residue name / residue number / atom name) may appear in any order — the parser falls back to a heuristic for non-standard orderings.

## Input file requirements

- **PDB** — must contain hydrogens (add via `reduce` / `pdb2pqr` / Open Babel / `mlmm mm-parm --add-h`) and element symbols in cols 77–78 (`mlmm add-elem-info` if missing). Multiple PDBs must share identical atoms in the same order.
- **XYZ** — accepted when ML-region determination is skipped (omit `-c/--center`).
- **Amber `--parm` (parm7)** — force-field parameters for the full system; atom ordering must match the input PDB exactly.

## Backend selection

All calc subcommands (`opt`, `tsopt`, `freq`, `irc`, `dft`, `scan` / `scan2d` / `scan3d`, `path-opt`, `path-search`, `all`) accept:

| Option | Description | Default |
|---|---|---|
| `-b, --backend` | MLIP backend: `uma`, `orb`, `mace`, `aimnet2` | `uma` |
| `--embedcharge` / `--no-embedcharge` | xTB point-charge embedding correction | off |

Install alternatives: `pip install "mlmm-toolkit[orb]"` / `"[aimnet]"` / `pip install --no-deps mace-torch` (MACE in a dedicated env).

## `--opt-mode` (subcommand-dependent)

```{warning}
Choices and defaults differ across subcommands.
```

| Subcommand | Gradient-only alias | Hessian-based alias | Default | Engine |
|---|---|---|---|---|
| `opt` | `grad` (`light`, `lbfgs`) | `hess` (`heavy`, `rfo`) | `grad` | L-BFGS vs RFO (with optional `--microiter`). |
| `tsopt` | `grad` (`light`, `dimer`) | `hess` (`heavy`, `rsirfo`) | `hess` | Dimer vs RS-I-RFO. |
| `path-search` | `grad` (= L-BFGS) | — (RFO/`hess` not yet wired) | `grad` | Inner optimiser for GSM / DMF nodes; `hess`/`rfo` rejected with a Click error. `path-opt` has no `--opt-mode` (use `--mep-mode gsm`/`dmf`). |
| `scan` / `scan2d` / `scan3d` | — | — | — | No `--opt-mode`; relaxation uses the YAML inner optimiser. |

`light` / `heavy` are accepted as legacy aliases for `grad` / `hess` for backward compatibility; prefer `grad` / `hess` in new scripts.

### YAML precedence

```
defaults < config < CLI options
```

Three tiers; `--config` is the only YAML layer (no separate "override YAML"). Full schema: [YAML Reference](yaml-reference.md).

## Output directory

`--out-dir ./my_results/` overrides. Defaults: `all → ./result_all/`, per-stage subcommand → `./result_<subcmd>/` (e.g. `result_opt/`, `result_tsopt/`, `result_path_search/`). `extract` / `mm-parm` / `define-layer` default to the current directory or the explicit `-o` / `--out-prefix`.

## See Also

[Getting Started](getting-started.md) · [Concepts & Workflow](concepts.md) · [Common Error Recipes](recipes-common-errors.md) · [Troubleshooting](troubleshooting.md) · [YAML Reference](yaml-reference.md).
