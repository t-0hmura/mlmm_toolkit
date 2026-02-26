"""Tests for defaults.py configuration constants."""

from mlmm_toolkit.defaults import (
    GEOM_KW_DEFAULT,
    MLMM_CALC_KW,
    OPT_BASE_KW,
    LBFGS_KW,
    RFO_KW,
    BIAS_KW,
    BOND_KW,
    GS_KW,
    STOPT_KW,
    SEARCH_KW,
    IRC_KW,
    FREQ_KW,
    THERMO_KW,
    DIMER_KW,
    HESSIAN_DIMER_KW,
    RSIRFO_KW,
    DFT_KW,
    BFACTOR_ML,
    BFACTOR_HESS_MM,
    BFACTOR_MOVABLE_MM,
    BFACTOR_FROZEN,
)


def test_geom_defaults():
    assert "coord_type" in GEOM_KW_DEFAULT
    assert GEOM_KW_DEFAULT["coord_type"] == "cart"


def test_mlmm_calc_defaults():
    assert isinstance(MLMM_CALC_KW, dict)
    assert "uma_model" in MLMM_CALC_KW
    assert "model_charge" in MLMM_CALC_KW


def test_opt_defaults():
    assert "thresh" in OPT_BASE_KW
    assert "max_cycles" in OPT_BASE_KW


def test_lbfgs_inherits_opt():
    for key in OPT_BASE_KW:
        assert key in LBFGS_KW, f"LBFGS_KW missing opt key: {key}"


def test_rfo_inherits_opt():
    for key in OPT_BASE_KW:
        assert key in RFO_KW, f"RFO_KW missing opt key: {key}"


def test_bfactors():
    assert BFACTOR_ML < BFACTOR_HESS_MM <= BFACTOR_MOVABLE_MM < BFACTOR_FROZEN


def test_all_defaults_are_dicts():
    for name, obj in [
        ("GEOM_KW_DEFAULT", GEOM_KW_DEFAULT),
        ("MLMM_CALC_KW", MLMM_CALC_KW),
        ("OPT_BASE_KW", OPT_BASE_KW),
        ("LBFGS_KW", LBFGS_KW),
        ("RFO_KW", RFO_KW),
        ("BIAS_KW", BIAS_KW),
        ("BOND_KW", BOND_KW),
        ("GS_KW", GS_KW),
        ("STOPT_KW", STOPT_KW),
        ("SEARCH_KW", SEARCH_KW),
        ("IRC_KW", IRC_KW),
        ("FREQ_KW", FREQ_KW),
        ("THERMO_KW", THERMO_KW),
        ("DIMER_KW", DIMER_KW),
        ("HESSIAN_DIMER_KW", HESSIAN_DIMER_KW),
        ("RSIRFO_KW", RSIRFO_KW),
        ("DFT_KW", DFT_KW),
    ]:
        assert isinstance(obj, dict), f"{name} should be a dict"
