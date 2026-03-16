# `mlmm pysis`

```text

Usage: mlmm pysis [YAML_FILE] [ARGS]...

  Run a pysisyphus YAML workflow file.

  This subcommand provides compatibility with v0.1.x YAML-based workflows.
  The ``mlmm`` calculator type is automatically registered so that
  ``calc: type: mlmm`` works in pysisyphus YAML files.

  Usage:
      mlmm pysis opt.yaml
      mlmm pysis tsopt.yaml -- --clean
```

See [pysis](../../pysis.md) for YAML format details and migration guide.
