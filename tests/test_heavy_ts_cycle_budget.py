from pathlib import Path

import pytest

from mlmm.workflows import tsopt
from mlmm.workflows.tsopt import _OptimizationCycleLedger


def test_heavy_ts_cycle_ledger_counts_every_trial():
    ledger = _OptimizationCycleLedger(3)

    ledger.debit(1)
    assert ledger.remaining == 2
    ledger.debit(1)  # rejected trials still consume command-level work
    ledger.debit(1)

    assert ledger.spent == 3
    assert ledger.remaining == 0


def test_heavy_ts_cycle_ledger_rejects_overspend():
    ledger = _OptimizationCycleLedger(2, spent=1)

    with pytest.raises(RuntimeError, match="exceeded"):
        ledger.debit(2)

    assert ledger.spent == 1


def test_heavy_ts_restarts_receive_only_remaining_cycle_budget():
    source = Path(tsopt.__file__).read_text(encoding="utf-8")

    assert 'restart_args["max_cycles"] = remaining_cycles' in source
    assert 'restart_opt_cfg["max_cycles"] = remaining_cycles' in source
    assert "_rsirfo_cycles_spent" not in source
    assert "_tsopt_n_opt_cycles = (" in source
    assert "_heavy_cycle_ledger.spent" in source
