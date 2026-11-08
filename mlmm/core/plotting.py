"""Per-invocation plotting policy shared by workflow helpers."""

import click

from mlmm.core.defaults import PLOT_DEFAULT
_KEY = "mlmm.plot_enabled"


def plot_enabled() -> bool:
    context = click.get_current_context(silent=True)
    return PLOT_DEFAULT if context is None else bool(context.meta.get(_KEY, PLOT_DEFAULT))


def _capture_plot(context, parameter, value):
    if not context.resilient_parsing:
        source = context.get_parameter_source(parameter.name)
        if source != click.core.ParameterSource.DEFAULT or _KEY not in context.meta:
            context.meta[_KEY] = bool(value)
    return value


def add_plot_option():
    return click.option(
        "--plot/--no-plot", default=PLOT_DEFAULT, show_default=True,
        expose_value=False, callback=_capture_plot,
        help="Write workflow plots. Use --no-plot to skip figures and Chrome rendering.",
    )
