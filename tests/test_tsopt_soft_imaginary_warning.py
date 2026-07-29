"""A soft leading imaginary mode must warn without changing the saddle status.

Saddle certification counts imaginary modes without a magnitude cutoff. The
soft-mode warning is diagnostic only and does not alter terminal status.
"""

import pytest

from mlmm.core.defaults import TS_IMAG_SOFT_WARN_CM
from mlmm.workflows.tsopt import _warn_if_leading_imaginary_mode_is_soft


def test_soft_leading_mode_warns(capsys) -> None:
    _warn_if_leading_imaginary_mode_is_soft([-15.72])
    out = capsys.readouterr().out
    assert "WARNING" in out and "-15.72 cm^-1" in out
    assert str(int(TS_IMAG_SOFT_WARN_CM)) in out


def test_large_magnitude_mode_is_silent(capsys) -> None:
    _warn_if_leading_imaginary_mode_is_soft([-447.78])
    assert capsys.readouterr().out == ""


def test_leading_mode_is_the_most_negative_one(capsys) -> None:
    # The warning evaluates the most negative mode, not a soft companion.
    _warn_if_leading_imaginary_mode_is_soft([-447.78, -5.72])
    assert capsys.readouterr().out == ""


@pytest.mark.parametrize("ims", [None, []])
def test_no_imaginary_modes_is_silent(capsys, ims) -> None:
    _warn_if_leading_imaginary_mode_is_soft(ims)
    assert capsys.readouterr().out == ""


def test_threshold_boundary_is_not_warned(capsys) -> None:
    _warn_if_leading_imaginary_mode_is_soft([-TS_IMAG_SOFT_WARN_CM])
    assert capsys.readouterr().out == ""


def test_ts_imag_record_carries_the_frequency() -> None:
    """summary.json/summary.log retain the leading imaginary frequency.

    The record includes the leading frequency so downstream summaries can
    display its magnitude.
    """
    from mlmm.workflows.all import _ts_imag_record

    soft = _ts_imag_record(1, [-15.72])
    assert soft["n_imag"] == 1
    assert soft["nu_imag_max_cm"] == pytest.approx(-15.72)
    assert soft["min_abs_imag_cm"] == pytest.approx(15.72)

    large = _ts_imag_record(1, [-447.78])
    assert large["nu_imag_max_cm"] == pytest.approx(-447.78)

    # The certifying mode is the most negative one, not the softest companion.
    two = _ts_imag_record(2, [-447.78, -5.72])
    assert two["nu_imag_max_cm"] == pytest.approx(-447.78)
    assert two["min_abs_imag_cm"] == pytest.approx(5.72)

    # No frequency data available: keep the old shape rather than invent a value.
    assert _ts_imag_record(1, None) == {"n_imag": 1}


def test_summary_renders_the_frequency_and_its_warning() -> None:
    """The summary warning fires when the recorded mode is soft."""
    from mlmm.io.summary import _format_ts_imag_info
    from mlmm.workflows.all import _ts_imag_record

    soft = "\n".join(_format_ts_imag_info(_ts_imag_record(1, [-15.72])))
    assert "-15.7 cm^-1" in soft
    assert "WARNING" in soft

    large = "\n".join(_format_ts_imag_info(_ts_imag_record(1, [-447.78])))
    assert "-447.8 cm^-1" in large
    assert "WARNING" not in large


def test_thermo_branch_does_not_clobber_the_tsopt_frequency() -> None:
    """`--thermo` preserves the frequency published by the TS stage.

    `all` writes `ts_imag` twice: once from the tsopt result (which carries
    `imaginary_frequencies_cm`) and again from the thermo payload (which carries
    only `num_imag_freq`—thermoanalysis.yaml has no frequency list). The latter
    therefore retains the earlier list.
    """
    from mlmm.workflows.all import _ts_imag_record

    from_tsopt = _ts_imag_record(1, [-444.08301776728405])
    assert from_tsopt["nu_imag_max_cm"] == pytest.approx(-444.083017767)

    # The thermo branch re-derives n_imag but has no list of its own; it must
    # fall back to what is already published rather than dropping it.
    prior = from_tsopt.get("imag_freqs_cm")
    from_thermo = _ts_imag_record(1, None or prior)
    assert from_thermo["nu_imag_max_cm"] == pytest.approx(-444.083017767)
    assert from_thermo["imag_freqs_cm"] == from_tsopt["imag_freqs_cm"]


def test_thermo_fallback_is_wired_in_the_product() -> None:
    """Pin the wiring, not just the helper: the fallback must exist in `all`."""
    import inspect

    from mlmm.workflows import all as all_mod

    src = inspect.getsource(all_mod)
    assert src.count('_prior_freqs = (segment_log.get("ts_imag") or {}).get(') == 2
    assert src.count("or _prior_freqs,") == 2
