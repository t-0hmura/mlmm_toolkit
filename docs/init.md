# `init` (removed)

> **Note:** The `init` subcommand has been removed. Use `--show-config` on any subcommand to see the current configuration, or write a YAML file manually following the [YAML Reference](yaml-reference.md).

## Previous behavior

The `init` command used to generate a starter YAML template:

```text
mlmm init --out mlmm_all.config.yaml
mlmm all --config mlmm_all.config.yaml --dry-run
```

## CLI options

| Option | Description | Default |
| --- | --- | --- |
| `-o, --out PATH` | Destination YAML path. | `mlmm_all.config.yaml` |
| `--force/--no-force` | Overwrite existing file. | `False` |

## See Also

- [Common Error Recipes](recipes-common-errors.md) -- Symptom-first failure routing
- [Troubleshooting](troubleshooting.md) -- Detailed troubleshooting guide
- [YAML Reference](yaml-reference.md) -- Full configurable keys
