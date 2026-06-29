"""Shared Click option decorators across mlmm subcommands.

Factories that collapse the identical
`--detect-layer / --no-detect-layer` and
`--model-indices-one-based / --model-indices-zero-based` pairs shared by the
ten subcommands that define an ML region
(`dft / freq / irc / opt / path_opt / path_search / scan / scan2d / scan3d /
tsopt`) into a single call site. Flag names, defaults, and help text are
byte-for-byte the same as the per-subcommand inline definitions they replaced;
only the order in `--help` shifts to the factory's call site.

`--model-pdb` and `--model-indices` are sibling options that vary their help
text per subcommand (ten variants), so they remain inline. If a future change
unifies that text deliberately, those can also be factored out.
"""

from __future__ import annotations

from typing import Callable, Sequence

import click


def add_print_every_option() -> Callable[[Callable], Callable]:
    """Attach `--print-every N` (debug verbosity throttle).

    Routes to pysisyphus `Optimizer.__init__(print_every=N)`. N=1 (pysis
    default) prints every cycle; larger N prints every N-th macro cycle.
    Useful when running long opts and only the periodic summary is needed.
    Default ``None`` so omission falls through to defaults / YAML.
    """
    def decorator(func: Callable) -> Callable:
        return click.option(
            "--print-every",
            "print_every",
            type=click.IntRange(min=1),
            default=None,
            show_default=False,
            help="Print optimizer status every N cycles (debug knob).",
        )(func)
    return decorator


def add_irc_pos_def_option() -> Callable[[Callable], Callable]:
    """Attach `--irc-pos-def/--no-irc-pos-def`.

    When enabled, IRC convergence additionally requires `eigvalsh(mw_hessian)[0] > 0`
    on top of rms(grad) <= threshold. Blocks the "shoulder" false-convergence
    where the IRC walker calls success on a downhill descent before reaching
    the local minimum. Default ``None`` falls through to the rms-only criterion.
    """
    def decorator(func: Callable) -> Callable:
        return click.option(
            "--irc-pos-def/--no-irc-pos-def",
            "irc_pos_def",
            default=None,
            show_default=False,
            help="Require pos-def Hessian at IRC convergence (blocks shoulder false-convergence).",
        )(func)
    return decorator


def add_coord_type_option(
    choices: Sequence[str] = ("cart", "redund", "dlc", "tric"),
) -> Callable[[Callable], Callable]:
    """Attach `--coord-type` to a Click command.

    Selects the optimization coordinate system passed through to pysisyphus'
    Geometry constructor. ``cart`` (default) preserves the historical
    Cartesian behavior used for the published paper data; ``dlc``
    (delocalised internal coordinates) often converges faster on torsion-
    rich systems but is brittle for bond-making, ML/MM link-atom paths and
    multi-fragment systems. ``redund`` and ``tric`` are accepted for
    single-structure optimizers (opt / tsopt / scan / freq) but NOT for
    Chain-of-States engines — ``path-opt`` and ``path-search`` pass
    ``choices=("cart", "dlc")`` here because pysisyphus' ChainOfStates only
    honours those two coordinate systems. Subcommands hard-coupled to
    Cartesian (irc, dft) skip the decorator entirely.

    Default is ``None`` so omission falls through to
    ``GEOM_KW_DEFAULT['coord_type']`` (cart) via the standard YAML override
    chain.

    NOTE: dest is `cli_coord_type` (not `coord_type`) because every
    downstream cli body already has a local `coord_type = geom_cfg.get(...)`
    right before the geom_loader call. Binding the Click param to
    `coord_type` would make Python treat the whole symbol as local to the
    closure and the assemble-block reference would UnboundLocalError.
    """
    options_str = "|".join(choices)
    def decorator(func: Callable) -> Callable:
        return click.option(
            "--coord-type",
            "cli_coord_type",
            type=click.Choice(list(choices), case_sensitive=False),
            default=None,
            show_default=False,
            help=(
                f"Optimization coordinate system ({options_str}). cart is the "
                f"robust default used in published numbers; dlc speeds up "
                f"torsion-rich opts. mlmm-specific caveats: DLC + link atom "
                f"and DLC + 3-layer frozen MM are numerically unverified."
            ),
        )(func)
    return decorator


def add_precision_option() -> Callable[[Callable], Callable]:
    """Attach `--precision fp32|fp64` to a Click command.

    Backend-agnostic precision flag. The CLI body routes the value into
    the backend-specific configuration key via
    ``mlmm.backends.apply_precision_to_calc_cfg``:

    - ``uma``  -> ``uma_precision`` ('fp32' | 'fp64')
    - ``orb``  -> ``orb_precision`` ('float32-high' | 'float64')
    - ``mace`` -> ``mace_dtype``    ('float32' | 'float64')
    - ``aimnet2`` -> fp32 no-op; fp64 rejected (model inputs are cast to
      float32 upstream, so fp64 cannot be honoured)

    fp64 base precision can have non-trivial TSopt/Hessian impact for
    OMol-trained UMA; for ORB/MACE the higher precision similarly costs
    throughput and can stabilise gradients/Hessians.

    Wire targets: every subcommand that constructs a backend calculator —
    currently ``sp``, ``opt``, ``tsopt``, ``freq``, ``irc``,
    ``scan`` / ``scan2d`` / ``scan3d``, ``path-opt``, ``path-search``, and
    ``all``.
    """
    def decorator(func: Callable) -> Callable:
        return click.option(
            "--precision",
            "precision",
            type=click.Choice(["fp32", "fp64"], case_sensitive=False),
            default=None,
            show_default=False,
            help=(
                "MLIP backend precision: fp32 (default) or fp64. Routed to "
                "backend-specific kwargs (UMA precision / ORB precision / "
                "MACE default_dtype). aimnet2: fp32 no-op; fp64 rejected."
            ),
        )(func)
    return decorator


# Deprecated alias for one cycle so external scripts that import
# `add_uma_precision_option` keep working. Internal call sites must use
# `add_precision_option`.
add_uma_precision_option = add_precision_option


def add_backend_model_option() -> Callable[[Callable], Callable]:
    """Attach ``--backend-model NAME`` to a Click command.

    Backend-agnostic model-variant override. The CLI body routes the value into
    the active backend's model kwarg via
    ``mlmm.backends.apply_backend_model_to_calc_cfg``:

    - ``uma``     -> ``uma_model``     (default ``uma-s-1p1``)
    - ``orb``     -> ``orb_model``     (default ``orb_v3_conservative_omol``)
    - ``mace``    -> ``mace_model``    (default ``MACE-OMOL-0``)
    - ``aimnet2`` -> ``aimnet2_model`` (default ``aimnet2``)

    Unset keeps the backend's built-in default model. Same wire targets as
    ``add_precision_option`` (every subcommand that constructs a backend
    calculator).
    """
    def decorator(func: Callable) -> Callable:
        return click.option(
            "--backend-model",
            "backend_model",
            type=str,
            default=None,
            show_default=False,
            help=(
                "Model variant for the selected --backend (e.g. "
                "uma-s-1p1 / uma-m-1p1 for uma, orb_v3_conservative_omol for orb, "
                "MACE-OMOL-0 / MACE-OFF23_small for mace). "
                "Default: the backend's built-in model."
            ),
        )(func)
    return decorator


def _deterministic_callback(ctx, param, value):
    """Eager callback: activate strict-deterministic mode when --deterministic
    is set. Process-global, so it covers every backend used in the run and
    every in-process child stage of ``all``. ``expose_value=False`` keeps it
    out of the command function signature (no body changes needed)."""
    if ctx.resilient_parsing:
        return value
    if value:
        from mlmm.backends._determinism import setup_deterministic
        setup_deterministic()
    return value


def add_deterministic_option() -> Callable[[Callable], Callable]:
    """Attach ``--deterministic/--no-deterministic`` to a Click command.

    Bit-reproducible mode: turns on ``torch.use_deterministic_algorithms`` plus
    an ``index_reduce_`` shim so repeated GPU runs are bit-identical. It is a
    process-global side effect applied via an eager, value-less callback, so it
    propagates to all backends and to the in-process child stages of ``all``
    without per-stage forwarding. Slower than the default, and raises (rather
    than silently degrading) if the torch build cannot honour strict mode.
    Default off; default runs carry ~1e-7 A scatter/atomic non-determinism that
    is chemically negligible. The env var ``MLMM_STRICT_DETERMINISTIC=1`` is the
    equivalent entry point for CI / the direct Python API.
    """
    def decorator(func: Callable) -> Callable:
        return click.option(
            "--deterministic/--no-deterministic",
            default=False,
            show_default=False,
            is_eager=True,
            expose_value=False,
            callback=_deterministic_callback,
            help=(
                "Strict bit-reproducible GPU runs (deterministic algorithms + "
                "index_reduce_ shim). Slower; raises if unsupported. Default off."
            ),
        )(func)
    return decorator


def _allow_charge_mult_mismatch_callback(ctx, param, value):
    """Eager callback: when --allow-charge-mult-mismatch is set, disable the ML-region
    electron-parity check process-globally (covers every backend + every in-process child
    stage of ``all`` without per-stage forwarding, like --deterministic). ``expose_value=False``
    keeps it out of the command signature."""
    if ctx.resilient_parsing:
        return value
    if value:
        from mlmm.core.utils import set_allow_charge_mult_mismatch
        set_allow_charge_mult_mismatch(True)
    return value


def add_allow_charge_mult_mismatch_option() -> Callable[[Callable], Callable]:
    """Attach ``--allow-charge-mult-mismatch`` to a Click command.

    Skips the ML-region charge/multiplicity electron-parity check (``validate_charge_spin``)
    and logs that it was skipped. For users who know their (charge, multiplicity) is intentional
    despite a parity warning -- e.g. a genuinely open-shell ML region, or a covalently-modified
    residue whose ML/MM cut leaves an unpaired electron. Process-global via an eager, value-less
    callback, so it propagates to every backend and child stage without per-stage forwarding.
    """
    def decorator(func: Callable) -> Callable:
        return click.option(
            "--allow-charge-mult-mismatch",
            is_flag=True,
            default=False,
            is_eager=True,
            expose_value=False,
            callback=_allow_charge_mult_mismatch_callback,
            help=(
                "Skip the ML-region charge/multiplicity electron-parity check (logs that it was "
                "skipped). For an intentional open-shell or covalently-cut ML region."
            ),
        )(func)
    return decorator


def add_ml_charge_spin_options() -> Callable[[Callable], Callable]:
    """Attach the standard ML region charge/spin triple to a Click command.

    Options: -q/--charge, -l/--ligand-charge, -m/--multiplicity (spin).
    All 3 wired subcommands (freq, opt, scan) share identical signature.
    scan2d/scan3d are NOT wired because their --charge help text differs
    ('ML-region total charge' with hyphen) while everything else matches;
    extracting them would require an extra `charge_help` parameter.
    """
    options = [
        click.option(
            "-q", "--charge",
            type=int,
            required=False,
            help="ML region charge. Required unless --ligand-charge is provided.",
        ),
        click.option(
            "-l", "--ligand-charge",
            type=str,
            default=None,
            show_default=False,
            help=(
                "Total charge or per-resname mapping (e.g., GPP:-3,SAM:1) used to derive "
                "charge when -q is omitted (requires PDB input or --ref-pdb)."
            ),
        ),
        click.option(
            "-m", "--multiplicity",
            "spin",
            type=int,
            default=None,
            show_default=False,
            help="Spin multiplicity (2S+1) for the ML region. Defaults to 1 when omitted.",
        ),
    ]

    def decorator(func: Callable) -> Callable:
        for opt in reversed(options):
            func = opt(func)
        return func

    return decorator


def add_ml_layer_detection_options() -> Callable[[Callable], Callable]:
    """Attach `--detect-layer` and `--model-indices-one-based` to a Click command.

    Both options have identical signature (default/help text) across all 10
    mlmm subcmds that use them, so the factory takes no parameters.
    """
    options = [
        click.option(
            "--detect-layer/--no-detect-layer",
            "detect_layer",
            default=True,
            show_default=True,
            help="Detect ML/MM layers from input PDB B-factors (ML=0, MovableMM=10, FrozenMM=20). "
                 "If disabled, you must provide --model-pdb or --model-indices.",
        ),
        click.option(
            "--model-indices-one-based/--model-indices-zero-based",
            "model_indices_one_based",
            default=True,
            show_default=True,
            help="Interpret --model-indices as 1-based (default) or 0-based.",
        ),
    ]

    def decorator(func: Callable) -> Callable:
        for opt in reversed(options):
            func = opt(func)
        return func

    return decorator
