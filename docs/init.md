# `init`

## Overview

> **Summary:** Write a starter YAML configuration file for `mlmm all`.

### At a glance
- **Use when:** You want a reproducible config-first workflow for `mlmm all`.
- **Method:** Generates a YAML template with default values for all configurable keys.
- **Outputs:** `mlmm_all.config.yaml` (or custom path via `--out`).
- **Defaults:** Output path `mlmm_all.config.yaml`; does not overwrite unless `--force`.
- **Next step:** Edit the generated YAML, then run `mlmm all --config mlmm_all.config.yaml`.

## Minimal example

```bash
mlmm init --out mlmm_all.config.yaml
mlmm all --config mlmm_all.config.yaml --dry-run
```

## CLI options

| Option | Description | Default |
| --- | --- | --- |
| `-o, --out PATH` | Destination YAML path. | `mlmm_all.config.yaml` |
| `--force/--no-force` | Overwrite existing file. | `False` |

## Notes
- For symptom-first diagnosis, start with [Common Error Recipes](recipes_common_errors.md), then use [Troubleshooting](troubleshooting.md) for detailed fixes.

- The generated file is a starting point, not an exhaustive schema.

---

## See Also

- [Common Error Recipes](recipes_common_errors.md) -- Symptom-first failure routing
- [Troubleshooting](troubleshooting.md) -- Detailed troubleshooting guide
- [YAML Reference](yaml_reference.md) -- Full configurable keys
