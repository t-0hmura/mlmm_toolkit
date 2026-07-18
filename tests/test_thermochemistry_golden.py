"""M62: independent golden thermochemistry vectors.

The reference values below are HARD-CODED physical constants of the ideal-gas
RRHO/QRRHO model. They are NOT produced by the ``thermoanalysis`` helpers at
assertion time; they were derived from the closed-form statistical-mechanics
expressions and are frozen here as golden data. Each fixture is a fully
specified input (geometry, masses, wavenumbers, energy, multiplicity, T, p) with
a fixed expected output vector, so a formula or unit-conversion drift in
production breaks the test.

Derivation (Hartree / Hartree·K^-1 for energies & entropies, J·mol^-1·K^-1 for
heat capacities), with SI/atomic-unit constants:

    J2AU  = 1 / E_h(J)          Hartree energy in J
    KB    = Boltzmann (J/K)     KBAU = KB * J2AU
    PLANCK= Planck (J s)
    C     = speed of light (m/s)
    AMU2KG= unified atomic mass unit (kg)
    R     = molar gas constant (J/mol/K)

    nu[Hz]     = C * wavenumber[cm^-1] * 100
    ZPE        = J2AU * sum(PLANCK * nu / 2)
    U_trans    = 3/2 * KBAU * T
    U_rot      = {atom:0, linear:1, nonlinear:3/2} * KBAU * T
    U_vib      = KBAU * sum(theta*(1/2 + 1/(exp(theta/T)-1))),  theta = PLANCK*nu/KB
    U_tot      = E_scf + U_trans + U_rot + U_vib ;  H = U_tot + KBAU*T
    S_trans    = KBAU*(ln q_trans + 5/2),  q_trans = (2*pi*M*AMU2KG*KB*T/h^2)^(3/2)*KB*T/p
    S_rot      = KBAU*(ln q_rot + {linear:1, nonlinear:3/2}); atom: 0
    S_vib(RRHO)= sum KBAU*((theta/T)/(exp(theta/T)-1) - ln(1-exp(-theta/T)))
    S_vib(QRRHO)= Chai-Head-Gordon weighted mix of harmonic & free-rotor entropies
                  (rotor cutoff 100 cm^-1)
    c_trans=3/2 R, c_rot={atom:0,lin:R,nonlin:3/2 R},
    c_vib  = R*sum(exp(q)*(q/(exp(q)-1))^2), q = PLANCK*nu/(KB*T)
    G = H - T*S_tot ;  dG = G - E_scf

The full-3N wavenumber vectors carry rigid zeros; production drops the
(6 - is_linear) smallest-|w| entries before the analysis, so only the listed
positive modes contribute.
"""

from __future__ import annotations

import warnings

import numpy as np
import pytest

from thermoanalysis.QCData import QCData
from thermoanalysis.thermo import thermochemistry


T_K = 298.15
P_PA = 101325.0

# Atomic-unit energy/entropy terms are exact to ~1e-16 here; a small allowance
# covers the 3x3 inertia eigendecomposition that enters rotational terms on other
# platforms, while still catching any real formula/unit drift (>> 1e-9).
RTOL = 1.0e-9
ATOL = 1.0e-16

# ---- fixtures: (coords3d [Ang], masses [amu], wavenumbers [cm^-1], E_scf [Ha], mult) ----
FIXTURES = {
    "atom_He": dict(
        coords=[[0.0, 0.0, 0.0]],
        masses=[4.002602],
        wavenumbers=[],
        scf=-2.9,
        mult=1,
    ),
    "linear_CO": dict(
        coords=[[0.0, 0.0, 0.0], [1.128, 0.0, 0.0]],
        masses=[12.011, 15.999],
        wavenumbers=[0.0, 0.0, 0.0, 0.0, 0.0, 2143.0],  # full 3N; five rigid zeros
        scf=-113.3,
        mult=1,
    ),
    "nonlinear_H2O": dict(
        coords=[[0.0, 0.0, 0.0], [0.757, 0.586, 0.0], [-0.757, 0.586, 0.0]],
        masses=[15.999, 1.008, 1.008],
        # full 3N; six rigid zeros + three positive modes
        wavenumbers=[0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1595.0, 3657.0, 3756.0],
        scf=-76.4,
        mult=1,
    ),
}

# ---- golden output vectors (frozen; see module docstring for the derivation) ----
GOLDEN = {
    "atom_He": {
        "ZPE": 0.0,
        "U_trans": 0.0014162773014667915,
        "U_rot": 0.0,
        "U_vib": 0.0,
        "U_therm": 0.0014162773014667915,
        "U_tot": -2.898583722698533,
        "H": -2.8976395378308886,
        "S_trans": 4.8007401001081675e-05,
        "S_rot": 0.0,
        "S_vib": 0.0,
        "S_el": 0.0,
        "S_tot": 4.8007401001081675e-05,
        "c_trans": 12.47169392722986,
        "c_rot": 0.0,
        "c_vib": 0.0,
        "c_tot": 12.47169392722986,
        "G": -2.9119529444393613,
        "dG": -0.011952944439361346,
        "linear": True,
        "atom_num": 1,
        "wavenumbers": [],
        "U_rot_factor": 0.0,  # atom
        "c_rot_units_of_R": 0.0,
    },
    "linear_CO": {
        "ZPE": 0.00488211322349645,
        "U_trans": 0.0014162773014667915,
        "U_rot": 0.0009441848676445277,
        "U_vib": 0.004882428305571276,
        "U_therm": 0.007242890474682595,
        "U_tot": -113.29275710952531,
        "H": -113.29181292465766,
        "S_trans": 5.7249504334823436e-05,
        "S_rot": 1.797380862167058e-05,
        "S_vib": 1.1575359674431613e-09,
        "S_el": 0.0,
        "S_tot": 7.522447049246146e-05,
        "c_trans": 12.47169392722986,
        "c_rot": 8.31446261815324,
        "c_vib": 0.028694302145572923,
        "c_tot": 20.814850847528675,
        "G": -113.31424110053499,
        "dG": -0.014241100534988504,
        "linear": True,
        "atom_num": 2,
        "wavenumbers": [2143.0],
        "U_rot_factor": 1.0,  # linear: one rotational degree pair (U_rot = 1*kB*T)
        "c_rot_units_of_R": 1.0,
    },
    "nonlinear_H2O": {
        "ZPE": 0.020521733979120868,
        "U_trans": 0.0014162773014667915,
        "U_rot": 0.0014162773014667915,
        "U_vib": 0.020525036939745188,
        "U_therm": 0.023357591542678774,
        "U_tot": -76.37664240845733,
        "H": -76.37569822358968,
        "S_trans": 5.515296339394666e-05,
        "S_rot": 1.8843562408808126e-05,
        "S_vib": 1.2518110998048041e-08,
        "S_el": 0.0,
        "S_tot": 7.400904391375283e-05,
        "c_trans": 12.47169392722986,
        "c_rot": 12.47169392722986,
        "c_vib": 0.22402650858867354,
        "c_tot": 25.16741436304839,
        "G": -76.39776402003257,
        "dG": 0.002235979967437629,
        "linear": False,
        "atom_num": 3,
        "wavenumbers": [1595.0, 3657.0, 3756.0],
        "U_rot_factor": 1.5,  # nonlinear: three rotational degrees (U_rot = 3/2*kB*T)
        "c_rot_units_of_R": 1.5,
    },
}

def _build_qc(fix):
    return QCData(
        {
            "coords3d": np.asarray(fix["coords"], dtype=float),
            "wavenumbers": np.asarray(fix["wavenumbers"], dtype=float),
            "scf_energy": fix["scf"],
            "masses": np.asarray(fix["masses"], dtype=float),
            "mult": fix["mult"],
        },
        # DELIBERATE: symmetry number sigma=1 for every fixture (including H2O,
        # whose physical C2v group has sigma=2). The golden S_rot/S_tot/G below are
        # computed for sigma=1, so this must stay "c1" -- "fixing" it to the
        # physical point group would shift H2O's rotational entropy by R*ln2 and
        # break the golden vectors.
        point_group="c1",
        mult=fix["mult"],
    )


def _run(fix, **kw):
    qc = _build_qc(fix)
    with warnings.catch_warnings():
        # QCData.rot_temperatures divides by a zero inertia eigenvalue for the
        # atom/linear cases; the result is unused there (is_atom/is_linear guard).
        warnings.simplefilter("ignore", RuntimeWarning)
        return thermochemistry(qc, T_K, pressure=P_PA, **kw)


@pytest.mark.parametrize("name", list(FIXTURES))
def test_golden_thermochemistry_vectors(name):
    """Atom/linear/nonlinear RRHO golden vectors: ZPE, U, H, G, S, Cv."""
    fix = FIXTURES[name]
    gold = GOLDEN[name]
    tr = _run(fix)  # library defaults: QRRHO, rotor cutoff 100, no invert/floor

    # Scalar golden vector (energies & entropies in Ha; heat capacities in J/mol/K).
    for key in (
        "ZPE", "U_trans", "U_rot", "U_vib", "U_therm", "U_tot", "H",
        "S_trans", "S_rot", "S_vib", "S_el", "S_tot",
        "c_trans", "c_rot", "c_vib", "c_tot", "G", "dG",
    ):
        got = float(getattr(tr, key))
        np.testing.assert_allclose(
            got, gold[key], rtol=RTOL, atol=ATOL,
            err_msg=f"{name}.{key}: {got!r} != golden {gold[key]!r}",
        )

    # Rotational-degree structure and surviving vibrational wavenumbers.
    assert bool(tr.linear) is gold["linear"]
    assert tr.atom_num == gold["atom_num"]
    np.testing.assert_allclose(
        np.asarray(tr.wavenumbers, dtype=float),
        np.asarray(gold["wavenumbers"], dtype=float),
        rtol=RTOL, atol=ATOL,
    )


@pytest.mark.parametrize("name", list(FIXTURES))
def test_golden_rotational_and_capacity_degrees(name):
    """U_rot / c_rot expose the atom(0) / linear(1) / nonlinear(3/2) degree count."""
    import scipy.constants as spc

    kbau = spc.Boltzmann / spc.value("Hartree energy")
    R = spc.R
    fix = FIXTURES[name]
    gold = GOLDEN[name]
    tr = _run(fix)

    np.testing.assert_allclose(
        float(tr.U_rot), gold["U_rot_factor"] * kbau * T_K, rtol=RTOL, atol=ATOL
    )
    np.testing.assert_allclose(
        float(tr.c_rot), gold["c_rot_units_of_R"] * R, rtol=RTOL, atol=ATOL
    )


# ---- QRRHO vs RRHO on a shared molecule with one low-frequency mode ----
QRRHO_FIXTURE = dict(
    coords=[[0.0, 0.0, 0.0], [0.757, 0.586, 0.0], [-0.757, 0.586, 0.0]],
    masses=[15.999, 1.008, 1.008],
    wavenumbers=[50.0, 1595.0, 3657.0],  # one 50 cm^-1 mode below the 100 rotor cutoff
    scf=-76.4,
    mult=1,
)
QRRHO_GOLDEN = {
    "qrrho": {"S_tot": 7.977484016142647e-05, "S_vib": 5.778314358671693e-06, "G": -76.40709112887673},
    "rrho": {"S_tot": 8.168603473335565e-05, "S_vib": 7.689508930600874e-06, "G": -76.40766095153835},
}


@pytest.mark.parametrize("kind", ["qrrho", "rrho"])
def test_golden_qrrho_and_rrho_low_frequency(kind):
    """The QRRHO rotor damping shifts S/G of a low-frequency mode by a golden amount."""
    gold = QRRHO_GOLDEN[kind]
    tr = _run(QRRHO_FIXTURE, kind=kind)
    for key in ("S_tot", "S_vib", "G"):
        np.testing.assert_allclose(
            float(getattr(tr, key)), gold[key], rtol=RTOL, atol=ATOL,
            err_msg=f"{kind}.{key}",
        )


def test_qrrho_damps_low_frequency_entropy_below_rrho():
    """QRRHO caps the diverging low-frequency harmonic entropy -> smaller S, larger G."""
    tr_q = _run(QRRHO_FIXTURE, kind="qrrho")
    tr_r = _run(QRRHO_FIXTURE, kind="rrho")
    assert float(tr_q.S_vib) < float(tr_r.S_vib)
    assert float(tr_q.S_tot) < float(tr_r.S_tot)
    assert float(tr_q.G) > float(tr_r.G)


# ---- M16: zero-point-energy scale factor is applied exactly once ----
ZPE_SCALE_FIXTURE = FIXTURES["nonlinear_H2O"]


@pytest.mark.parametrize("f", [1.0, 0.5, 0.9])
def test_zpe_scale_factor_is_linear_single_scaling(f):
    """ZPE == f*raw and U/H/G shift by exactly (f-1)*raw (M16, non-squared).

    The f == 1.0 production result is the raw baseline; the relation asserted is
    the golden (it is not derived from a production partition/energy helper).
    """
    base = _run(ZPE_SCALE_FIXTURE, zpe_scale_factor=1.0)
    raw = float(base.ZPE)
    assert raw > 0.0

    tr = _run(ZPE_SCALE_FIXTURE, zpe_scale_factor=f)

    # Reported ZPE is exactly linear in f (a squared scaling would give f**2 * raw).
    np.testing.assert_allclose(float(tr.ZPE), f * raw, rtol=1.0e-12, atol=0.0)

    # U_tot / H / G / dG inherit the same additive (f - 1) * raw shift.
    for key in ("U_tot", "H", "G", "dG", "U_vib", "U_therm"):
        np.testing.assert_allclose(
            float(getattr(tr, key)),
            float(getattr(base, key)) + (f - 1.0) * raw,
            rtol=0.0, atol=1.0e-12,
            err_msg=f"{key} did not shift by (f-1)*raw at f={f}",
        )

    # Thermal excitation above the raw ZPE is invariant under the scale factor.
    np.testing.assert_allclose(
        float(tr.U_therm) - float(tr.ZPE),
        float(base.U_therm) - raw,
        rtol=0.0, atol=1.0e-12,
    )

    # Entropy and heat capacity are untouched by the ZPE scale factor.
    for key in ("S_tot", "S_vib", "c_tot", "c_vib"):
        np.testing.assert_allclose(
            float(getattr(tr, key)), float(getattr(base, key)), rtol=1.0e-12, atol=0.0
        )


# ---- imaginary-inversion / positive-frequency-floor policy (M62 pt.4) ----
IMAG_FIXTURE = dict(
    coords=[[0.0, 0.0, 0.0], [0.757, 0.586, 0.0], [-0.757, 0.586, 0.0]],
    masses=[15.999, 1.008, 1.008],
    wavenumbers=[-10.0, 10.0, 100.0],  # vibrational-only (size != 3N -> no drop)
    scf=-76.4,
    mult=1,
)


def test_imaginary_inversion_and_floor_policy_wavenumbers():
    """Workflow policy drops the -10; Geometry policy inverts it and floors to 25."""
    # Workflow/library policy: no inversion, no floor -> the -10 mode is discarded.
    tr_workflow = _run(IMAG_FIXTURE, invert_imags=0.0, cutoff=0.0)
    np.testing.assert_allclose(
        np.asarray(tr_workflow.wavenumbers, dtype=float),
        np.array([10.0, 100.0]),
        rtol=RTOL, atol=ATOL,
    )

    # Geometry policy: invert imaginaries down to -15, floor positives below 25.
    tr_geometry = _run(IMAG_FIXTURE, invert_imags=-15.0, cutoff=25.0)
    np.testing.assert_allclose(
        np.asarray(tr_geometry.wavenumbers, dtype=float),
        np.array([25.0, 25.0, 100.0]),
        rtol=RTOL, atol=ATOL,
    )
