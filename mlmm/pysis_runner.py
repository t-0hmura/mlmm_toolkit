"""Thin wrapper to run pysisyphus YAML workflows with mlmm calculator registered.

This module provides the ``mlmm pysis`` subcommand, which is compatible with
v0.1.x YAML-based workflows (``mlmm opt.yaml`` in the old interface).

Usage::

    mlmm pysis opt.yaml
    mlmm pysis tsopt.yaml
"""

import sys

import click


@click.command(
    context_settings={"ignore_unknown_options": True, "allow_extra_args": True},
)
@click.argument("yaml_file", required=False)
@click.pass_context
def cli(ctx, yaml_file):
    """Run a pysisyphus YAML workflow file.

    This subcommand provides compatibility with v0.1.x YAML-based workflows.
    The ``mlmm`` calculator type is automatically registered so that
    ``calc: type: mlmm`` works in pysisyphus YAML files.

    \b
    Usage:
        mlmm pysis opt.yaml
        mlmm pysis tsopt.yaml -- --clean
    """
    from pysisyphus import run as _run

    from mlmm.mlmm_calc import mlmm as MLMMCalc

    _run.CALC_DICT["mlmm"] = MLMMCalc

    if yaml_file:
        sys.argv = ["pysisyphus", yaml_file] + ctx.args
    else:
        sys.argv = ["pysisyphus"] + ctx.args

    _run.run()
