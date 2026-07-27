"""A soft leading imaginary mode must warn without changing the saddle status.

Saddle certification counts imaginary modes and does not weigh them, so a
few-cm^-1 soft mode certifies exactly like a real reaction coordinate. Two
runs of the same input from bit-identical starting geometries were observed to
diverge and land on a real saddle (~-450 cm^-1) and on a soft-mode structure
(~-15 cm^-1) respectively, with the soft one still reported as n_imag=1. The
warning is diagnostic only: the terminal status must stay exactly as it was.
"""

import pytest

from mlmm.core.defaults import TS_IMAG_SOFT_WARN_CM
from mlmm.workflows.tsopt import _warn_if_leading_imaginary_mode_is_soft


def test_soft_leading_mode_warns(capsys) -> None:
    _warn_if_leading_imaginary_mode_is_soft([-15.72])
    out = capsys.readouterr().out
    assert "WARNING" in out and "-15.72 cm^-1" in out
    assert str(int(TS_IMAG_SOFT_WARN_CM)) in out


def test_real_reaction_coordinate_is_silent(capsys) -> None:
    _warn_if_leading_imaginary_mode_is_soft([-447.78])
    assert capsys.readouterr().out == ""


def test_leading_mode_is_the_most_negative_one(capsys) -> None:
    # n_imag=2 with a real coordinate plus a soft companion: the certifying
    # mode is the stiff one, so this must stay silent.
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
    """summary.json/summary.log must distinguish a soft mode from a real one.

    `mlmm/io/summary.py` has always been able to print `nu_imag (max)` and to
    warn when the magnitude is small, but the producer wrote a bare
    `{"n_imag": n}`, so the field rendered as `-` and the warning never fired.
    A -15 cm^-1 run and a -450 cm^-1 run were indistinguishable downstream.
    """
    from mlmm.workflows.all import _ts_imag_record

    soft = _ts_imag_record(1, [-15.72])
    assert soft["n_imag"] == 1
    assert soft["nu_imag_max_cm"] == pytest.approx(-15.72)
    assert soft["min_abs_imag_cm"] == pytest.approx(15.72)

    real = _ts_imag_record(1, [-447.78])
    assert real["nu_imag_max_cm"] == pytest.approx(-447.78)

    # The certifying mode is the most negative one, not the softest companion.
    two = _ts_imag_record(2, [-447.78, -5.72])
    assert two["nu_imag_max_cm"] == pytest.approx(-447.78)
    assert two["min_abs_imag_cm"] == pytest.approx(5.72)

    # No frequency data available: keep the old shape rather than invent a value.
    assert _ts_imag_record(1, None) == {"n_imag": 1}


def test_summary_renders_the_frequency_and_its_warning() -> None:
    """The dormant summary warning must actually fire once the data flows."""
    from mlmm.io.summary import _format_ts_imag_info
    from mlmm.workflows.all import _ts_imag_record

    soft = "\n".join(_format_ts_imag_info(_ts_imag_record(1, [-15.72])))
    assert "-15.7 cm^-1" in soft
    assert "WARNING" in soft

    real = "\n".join(_format_ts_imag_info(_ts_imag_record(1, [-447.78])))
    assert "-447.8 cm^-1" in real
    assert "WARNING" not in real


def test_thermo_branch_does_not_clobber_the_tsopt_frequency() -> None:
    """`--thermo` must not erase the frequency the tsopt branch published.

    `all` writes `ts_imag` twice: once from the tsopt result (which carries
    `imaginary_frequencies_cm`) and again from the thermo payload (which carries
    only `num_imag_freq` — thermoanalysis.yaml has no frequency list). The second
    write overwrote the first, so every `--thermo` run published a bare
    `{"n_imag": n}` and the smoke floor rejected it even though the run had found
    a -444 cm^-1 saddle.
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
