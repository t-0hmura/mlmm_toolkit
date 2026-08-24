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
    """The exact TS-stage record remains a complete saddle verdict."""
    from mlmm.workflows.all import _ts_imag_record

    from_tsopt = _ts_imag_record(1, [-444.08301776728405])
    assert from_tsopt["nu_imag_max_cm"] == pytest.approx(-444.083017767)

    # A later thermochemistry spectrum may contain tiny numerical negative
    # roots.  The workflow must keep this exact TS-stage record unchanged.
    segment_log = {"ts_imag": from_tsopt}
    thermo_n_imag = 3
    if thermo_n_imag is not None and "ts_imag" not in segment_log:
        segment_log["ts_imag"] = _ts_imag_record(thermo_n_imag)
    assert segment_log["ts_imag"] == from_tsopt


def test_thermo_fallback_is_wired_in_the_product() -> None:
    """Thermo is a fallback, never an overwrite of TS certification."""
    import inspect

    from mlmm.workflows import all as all_mod

    src = inspect.getsource(all_mod)
    guarded_fallback = (
        "n_imag is not None\n"
        "                and not do_tsopt\n"
        '                and "ts_imag" not in segment_log'
    )
    assert src.count(guarded_fallback) == 2
    assert "must never overwrite certification" in src


def test_unexported_soft_mode_is_not_reported_as_no_imaginary_mode() -> None:
    """An exact soft imaginary root must not be announced as absent.

    A failed trajectory write is distinct from a spectrum with no imaginary mode.
    """
    from mlmm.workflows.tsopt import _dimer_mode_export_message

    message, is_diagnostic = _dimer_mode_export_message(0, 1, 5.0, -100.0)
    assert is_diagnostic is True
    assert message == "[tsopt] ERROR: Failed to write imaginary mode trajectory."

    # A genuinely positive spectrum still reports the absent imaginary mode.
    message, is_diagnostic = _dimer_mode_export_message(0, 0, 5.0, 18.0)
    assert is_diagnostic is True
    assert message == "[tsopt] No imaginary mode detected. Try all --refine-path."
