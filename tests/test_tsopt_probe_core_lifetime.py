"""The TSOPT probe core must not survive into the IRC phase.

``_run_tsopt_on_hei`` builds one ML/MM core purely to put an energy on the
returned TS geometry.  ``_irc_and_match`` then builds the segment's own leased
core and adopts that geometry with ``if g_ts.calculator is None:
lease.attach(g_ts)``.  If the probe core is still attached, that guard is False,
the lease never owns the geometry, and two heavy cores stay resident across the
TSOPT->IRC handoff and can exhaust memory on a small GPU.

The release idiom is load-bearing and easy to "tidy" into a bug:
``Geometry.set_calculator`` defaults to ``clear=True`` and would drop the energy
that was just computed, so the detach must be a direct attribute assignment.
``CalculatorLease.release`` uses the same idiom for the same reason.

The end-to-end behaviour (at most one live core across the handoff) is exercised
by the GPU smoke lane; these tests pin the two contracts a unit test can hold.
"""

from __future__ import annotations

import gc
import weakref

import numpy as np
import pytest

from pysisyphus.Geometry import Geometry


class _ProbeCalc:
    """Minimal stand-in for a heavy ML/MM core."""

    def __init__(self) -> None:
        self.closed = False

    def get_energy(self, atoms, coords):
        return {"energy": -42.5}

    def close(self) -> None:
        self.closed = True


def _h2() -> Geometry:
    return Geometry(["H", "H"], np.array([0.0, 0.0, 0.0, 0.0, 0.0, 1.4]))


def test_direct_detach_preserves_the_energy_the_probe_computed() -> None:
    """The detach idiom used by _run_tsopt_on_hei must keep the cached energy."""
    geom = _h2()
    geom.set_calculator(_ProbeCalc())
    energy = float(geom.energy)

    geom.calculator = None

    assert geom.calculator is None
    assert float(geom.energy) == pytest.approx(energy)


def test_set_calculator_none_would_lose_the_energy() -> None:
    """Pin why the detach is a direct assignment and not set_calculator(None).

    This is the trap: set_calculator(None) clears the results, so a later
    ``float(g_ts.energy)`` has neither a cache nor a calculator to ask.
    """
    geom = _h2()
    geom.set_calculator(_ProbeCalc())
    _ = float(geom.energy)

    geom.set_calculator(None)

    with pytest.raises(AttributeError):
        float(geom.energy)


def test_detached_probe_core_is_closed_and_collectable() -> None:
    """A detached probe core must not outlive the handoff."""
    geom = _h2()
    calc = _ProbeCalc()
    geom.set_calculator(calc)
    _ = float(geom.energy)

    geom.calculator = None
    calc.close()
    ref = weakref.ref(calc)
    del calc
    gc.collect()

    assert ref() is None, "probe core survived the TSOPT->IRC handoff"


def test_lease_can_adopt_a_detached_geometry() -> None:
    """The IRC lease's adoption guard only fires on a detached geometry."""
    from mlmm.workflows._run_session import CalculatorLease

    geom = _h2()
    geom.set_calculator(_ProbeCalc())
    _ = float(geom.energy)
    geom.calculator = None

    lease = CalculatorLease(_ProbeCalc())
    if geom.calculator is None:
        lease.attach(geom)

    assert geom.calculator is lease.calculator, "lease failed to adopt the TS geometry"

    lease.release()
    assert geom.calculator is None
