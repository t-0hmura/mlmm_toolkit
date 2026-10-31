import pytest

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
