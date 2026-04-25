# `mlmm pysis`

## Purpose

Run a pysisyphus YAML workflow directly through `mlmm-toolkit`'s
bundled `pysisyphus` fork (v0.1.x compatibility shim). This is the
escape hatch for using pysisyphus features that don't have a
first-class `mlmm` subcommand — e.g. NEB, custom optimizers, growing
string with non-default kwargs.

For most workflows you don't need this — `mlmm path-search` /
`tsopt` / `irc` / `freq` cover the common cases through the same
bundled pysisyphus.

## Synopsis

```bash
mlmm pysis -i workflow.yaml [--ref-pdb reference.pdb] [-o ./pysis_out/]
```

## Key flags

| flag | type | default | description |
|---|---|---|---|
| `-i, --input` | path | required | pysisyphus YAML workflow file |
| `--ref-pdb` | path | none | Reference PDB for residue context |
| `-o, --out-dir` | path | `./pysis_out/` | Output directory |
| `--show-config` / `--dry-run` | flag | off | Print resolved config |

## YAML workflow

The YAML follows pysisyphus' v0.1 schema. A minimal example:

```yaml
geom:
  type: cart
  fn: ts_candidate.xyz
calc:
  type: mlmm                    # ML/MM coupled calculator
  ml_backend: uma
  charge: 0
  spin: 1
opt:
  type: rsirfo
  max_cycles: 200
  trust_radius: 0.3
freq:
  num_modes: 10
```

Run:

```bash
mlmm pysis -i workflow.yaml -o ./pysis_run/
```

`mlmm pysis` substitutes `calc.type: mlmm` with the toolkit's ONIOM
calculator (so `ml_backend`, `parm7`, `ref_pdb` are available);
otherwise the YAML works as in upstream pysisyphus.

## When to use this vs first-class subcommands

Use `pysis` when:

- The workflow needs NEB (not exposed by `mlmm path-search`).
- You want to chain multiple optimizers in non-standard order.
- You're prototyping a feature that isn't yet a subcommand.

Don't use `pysis` for routine TS / IRC / freq — the dedicated
subcommands handle MLIP / MM coupling and output schemas more
robustly.

## Output

```
pysis_out/
├── result.json
├── final_geom.xyz                # final geometry (if applicable)
├── opt.log                       # pysisyphus's own log
└── (engine-specific subdirs)
```

`result.json` schema is pysisyphus-native, not `mlmm`'s standard.
Don't expect `barrier_kcal` / `bond_changes` keys — parse the
pysisyphus log for those.

## Caveats

- Bundled `pysisyphus` is a **fork**: do **not** mix with an
  upstream-installed pysisyphus in the same env (will collide).
- Many YAML keys behave differently from mainline pysisyphus
  documentation (the fork has GPU-tensor extensions). Check
  `import pysisyphus; help(pysisyphus.<module>)` for current API.
- Output is stable but does not follow `mlmm`'s standard
  `summary.json` — downstream analysis scripts that expect
  `mlmm`-shaped output will need adapters.

## See also

- `path-search.md`, `tsopt.md`, `irc.md`, `freq.md` — first-class
  subcommands; prefer them when applicable.
- bundled pysisyphus: explore via
  `python -c "import pysisyphus; help(pysisyphus)"`.