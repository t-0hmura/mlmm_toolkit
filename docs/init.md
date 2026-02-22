# `init`

## Overview

`mlmm init` writes a starter YAML configuration file for `mlmm all`.

Use this when you want a reproducible config-first workflow:

```bash
mlmm init --out mlmm_all.config.yaml
mlmm all --config mlmm_all.config.yaml --dry-run
```

## Usage

```bash
mlmm init [--out PATH] [--force]
```

## Options

| Option | Description | Default |
| --- | --- | --- |
| `-o, --out PATH` | Destination YAML path. | `mlmm_all.config.yaml` |
| `--force/--no-force` | Overwrite existing file. | `False` |

## Notes

- The generated file is a starting point, not an exhaustive schema.
- Runtime precedence is: `defaults < --config < CLI < --override-yaml`.

---

## See Also

- [Common Error Recipes](recipes_common_errors.md) -- Symptom-first failure routing
- [Troubleshooting](troubleshooting.md) -- Detailed troubleshooting guide
- [all](all.md) -- Run with `--config` / `--override-yaml`
- [YAML Reference](yaml_reference.md) -- Full configurable keys
