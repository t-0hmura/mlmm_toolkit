"""Shared Click option decorators for mlmm subcommands."""

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
            show_default="100",
            help="Print optimizer status every N cycles.",
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
            show_default="no-irc-pos-def (rms-only criterion)",
            help="Require pos-def Hessian at IRC convergence (blocks shoulder false-convergence).",
        )(func)
    return decorator


def add_coord_type_option(
    choices: Sequence[str] = ("cart", "redund", "dlc", "tric"),
) -> Callable[[Callable], Callable]:
    """Attach `--coord-type` to a Click command.

    Selects the optimization coordinate system passed through to pysisyphus'
    Geometry constructor. ``cart`` is the default. ``redund`` and ``tric`` are accepted for
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
            show_default="cart",
            help=(
                f"Optimization coordinate system ({options_str}). cart is the "
                f"default; command-specific choices are listed here."
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

    Unset resolves per backend: UMA fp32, ORB and MACE fp64.

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
            show_default="per backend: uma fp32; orb, mace fp64",
            help=(
                "MLIP backend precision: fp32 or fp64. Unset defaults per "
                "backend (uma: fp32; orb, mace: fp64). Routed to "
                "backend-specific kwargs (UMA precision / ORB precision / "
                "MACE default_dtype). aimnet2: fp32 no-op; fp64 rejected."
            ),
        )(func)
    return decorator

def add_workers_options() -> Callable[[Callable], Callable]:
    """Attach ``--workers`` / ``--workers-per-node`` to a Click command.

    MLIP predictor parallelism. ``--workers > 1`` routes the UMA backend through
    ``ParallelMLIPPredictUnit`` (needs ``fairchem-core[extras]``); the parallel
    predictor exposes no autograd model, so analytical Hessians are unavailable and
    combining ``--workers >1`` with ``--hessian-calc-mode Analytical`` is an error.

    The CLI body routes the values via ``mlmm.backends.apply_workers_to_calc_cfg``.
    Wire targets: every subcommand that constructs a backend calculator — ``sp``,
    ``opt``, ``tsopt``, ``freq``, ``irc``, ``scan`` / ``scan2d`` / ``scan3d``,
    ``path-opt``, ``path-search``, and ``all``.
    """
    def decorator(func: Callable) -> Callable:
        func = click.option(
            "--workers-per-node",
            "workers_per_node",
            type=int,
            default=None,
            show_default="1",
            help="Workers per node when the parallel MLIP predictor is used (--workers > 1).",
        )(func)
        func = click.option(
            "--workers",
            "workers",
            type=int,
            default=None,
            show_default="1",
            help=(
                "MLIP predictor workers (UMA). >1 uses a parallel predictor "
                "(fairchem-core[extras]); combining it with an analytical Hessian "
                "is an error. Default 1."
            ),
        )(func)
        return func
    return decorator


def add_backend_model_option() -> Callable[[Callable], Callable]:
    """Attach ``--backend-model NAME`` to a Click command.

    Backend-agnostic model-variant override. The CLI body routes the value into
    the active backend's model kwarg via
    ``mlmm.backends.apply_backend_model_to_calc_cfg``:

    - ``uma``     -> ``uma_model``     (default ``uma-s-1p2``)
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
            show_default="the selected backend's own model",
            help=(
                "Model variant for the selected --backend (e.g. "
                "uma-s-1p2 / uma-m-1p1 for uma, orb_v3_conservative_omol for orb, "
                "MACE-OMOL-0 / off:small for mace). "
                "Default: the backend's built-in model."
            ),
        )(func)
    return decorator


def add_calc_file_option() -> Callable[[Callable], Callable]:
    """Attach ``--calc-file PATH`` (+ ``--calc-file-func-name NAME``) to a Click command.

    When set, mlmm drives the ML region with an arbitrary ASE Calculator loaded
    from the user Python file (overriding ``--backend``). The CLI body routes the
    value via ``mlmm.backends.apply_calc_file_to_calc_cfg``, which switches the
    backend to ``custom``. The file must expose a factory
    ``get_calculator(charge, spin, device, **kwargs) -> ase Calculator`` (rename
    via ``--calc-file-func-name``). Lets users couple GFN-xTB (tblite), DFTB+, ORCA, or
    any ASE-compatible engine for the ML region without modifying mlmm. Same wire
    targets as ``add_precision_option``.
    """
    def decorator(func: Callable) -> Callable:
        func = click.option(
            "--calc-file-func-name",
            "calc_factory",
            type=str,
            default=None,
            show_default="get_calculator",
            help=(
                "Name of the callable in --calc-file that returns an ASE "
                "Calculator (or a module-level Calculator instance). "
                "CLI overrides config YAML; otherwise defaults to get_calculator."
            ),
        )(func)
        func = click.option(
            "--calc-file",
            "calc_file",
            type=click.Path(exists=True, dir_okay=False),
            default=None,
            help=(
                "Python file exposing get_calculator(...) -> an ASE Calculator "
                "used as the ML-region backend (overrides --backend). Couples "
                "GFN-xTB / DFTB+ / any ASE engine. See --calc-file-func-name."
            ),
        )(func)
        return func
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

    Requests deterministic algorithms for controlled operations by enabling
    ``torch.use_deterministic_algorithms`` and an ``index_reduce_`` shim. Exact
    reproducibility must be validated on the complete target stack. The setting
    is process-global and propagates to in-process ``all`` stages. Default off;
    ``MLMM_STRICT_DETERMINISTIC=1`` is the direct-API equivalent.
    """
    def decorator(func: Callable) -> Callable:
        return click.option(
            "--deterministic/--no-deterministic",
            default=False,
            show_default=True,
            is_eager=True,
            expose_value=False,
            callback=_deterministic_callback,
            help=(
                "Request deterministic algorithms for controlled operations; "
                "verify exact reproducibility on the complete target stack."
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
    and logs that it was skipped. Every valid integer-electron spin state obeys the parity
    relation, so an open-shell ML region needs a matching multiplicity rather than this flag.
    It exists for genuinely intentional nonstandard inputs, such as a covalently-modified
    residue whose ML/MM cut leaves an unpaired electron. Process-global via an eager,
    value-less callback, so it propagates to every backend and child stage without
    per-stage forwarding.
    """
    def decorator(func: Callable) -> Callable:
        return click.option(
            "--allow-charge-mult-mismatch",
            is_flag=True,
            default=False, show_default=True,
            is_eager=True,
            expose_value=False,
            callback=_allow_charge_mult_mismatch_callback,
            help=(
                "Skip the ML-region charge/multiplicity electron-parity check (logs that it was "
                "skipped). An open-shell ML region needs a matching multiplicity; use this only "
                "for an intentional nonstandard input such as a covalently-cut region."
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
                "Total charge for unknown ligand residues or a per-resname mapping "
                "(e.g., GPP:-3,SAM:1), used to derive the ML-region charge when -q "
                "is omitted (requires PDB input or --ref-pdb)."
            ),
        ),
        click.option(
            "-m", "--multiplicity",
            "spin",
            type=click.IntRange(min=1),
            default=None,
            show_default="1",
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

    Both options have identical signatures (default/help text) across all 11
    mlmm subcmds that use them, so the factory takes no parameters.
    """
    options = [
        click.option(
            "--detect-layer/--no-detect-layer",
            "detect_layer",
            default=True,
            show_default=True,
            help="Automatically detect ML/MM layers from input PDB B-factors "
                 "(ML=0, MovableMM=10, FrozenMM=20) when explicit ML membership "
                 "is absent. With explicit membership, retain valid movable/frozen "
                 "MM B-factor layers.",
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
