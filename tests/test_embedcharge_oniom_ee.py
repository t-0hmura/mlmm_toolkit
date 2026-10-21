"""Electrostatic embedding follows ONIOM-EE (Vreven et al., JCTC 2006, 2, 815).

The embedded MM point charges enter the high level and, through the ML
force-field charges, the model-low term, so the ML-MM Coulomb that E_real_low
holds for them cancels. MM atoms one or two bonds from the ML region stay out
of the point charges. The Li+/Cl- ion-pair test is in
test_dft_embedding_charges.py.
"""

from __future__ import annotations

import shutil
from types import MethodType

import numpy as np
import pytest

pmd = pytest.importorskip("parmed")

BOHR = 0.529177210903


def _carbon(name: str):
    return pmd.Atom(name=name, type="CT", charge=0.0, mass=12.011, atomic_number=6)


def test_mm_atoms_one_or_two_bonds_from_the_ml_region_are_near() -> None:
    from mlmm.backends.mlmm_calc import EMBED_EXCLUDED_BOND_DEPTH, _atoms_within_bonds

    s = pmd.Structure()
    for i in range(6):
        s.add_atom(_carbon(f"C{i}"), "MOL", 1)
    s.add_atom(_carbon("X"), "ION", 2)  # not bonded to anything
    for i in range(5):
        s.bonds.append(pmd.Bond(s.atoms[i], s.atoms[i + 1]))

    assert EMBED_EXCLUDED_BOND_DEPTH == 2
    assert _atoms_within_bonds(s, [0, 1], EMBED_EXCLUDED_BOND_DEPTH) == {2, 3}


def test_embedded_point_charges_skip_near_atoms_and_stay_fixed() -> None:
    from ase import Atoms

    from mlmm.backends.mlmm_calc import MLMMCore

    core = MLMMCore.__new__(MLMMCore)
    core.selection_indices = [0]
    core._embedding_near_mm = {1}
    core.embedcharge_cutoff = 2.5
    first = Atoms("HHHH", positions=[[0, 0, 0], [1, 0, 0], [2, 0, 0], [5, 0, 0]])
    moved = Atoms("HHHH", positions=[[0, 0, 0], [1, 0, 0], [5, 0, 0], [2, 0, 0]])

    assert core.embedding_point_charge_indices(first) == (2,)
    assert core.embedding_point_charge_indices(moved) == (2,)


def _fd_forces(energy, x0: np.ndarray, step: float = 1.0e-5) -> np.ndarray:
    flat = x0.reshape(-1)
    out = np.zeros(flat.size)
    for k in range(flat.size):
        xp, xm = flat.copy(), flat.copy()
        xp[k] += step
        xm[k] -= step
        out[k] = -(energy(xp.reshape(x0.shape)) - energy(xm.reshape(x0.shape))) / (2 * step)
    return out


def test_point_charge_coulomb_derivatives_match_finite_differences() -> None:
    from mlmm.backends.mlmm_calc import point_charge_coulomb, point_charge_coulomb_hessian

    rng = np.random.default_rng(0)
    ml = rng.normal(size=(3, 3))
    pc = rng.normal(size=(4, 3)) + 4.0
    q_ml = rng.normal(size=3)
    q_pc = rng.normal(size=4)
    x0 = np.vstack([ml, pc])

    def energy(x):
        return point_charge_coulomb(x[:3], q_ml, x[3:], q_pc)[0]

    def forces(x):
        _, f_ml, f_pc = point_charge_coulomb(x[:3], q_ml, x[3:], q_pc, need_forces=True)
        return np.vstack([f_ml, f_pc]).reshape(-1)

    assert np.allclose(forces(x0), _fd_forces(energy, x0), atol=1.0e-6)

    hess = point_charge_coulomb_hessian(ml, q_ml, pc, q_pc, [0, 1, 2], [0, 1, 2, 3])
    fd_hess = np.column_stack([
        -(forces(xp) - forces(xm)) / 2.0e-5
        for xp, xm in (
            ((x0.reshape(-1) + e).reshape(x0.shape), (x0.reshape(-1) - e).reshape(x0.shape))
            for e in np.eye(x0.size) * 1.0e-5
        )
    ])
    assert np.allclose(hess, fd_hess, atol=1.0e-5)

    sub = point_charge_coulomb_hessian(ml, q_ml, pc, q_pc, [1], [2])
    rows = [3, 4, 5, 15, 16, 17]
    assert np.allclose(sub, hess[np.ix_(rows, rows)])


@pytest.mark.parametrize("method", ["scaled", "fixed"])
def test_link_hessian_mapping_matches_finite_differences(method) -> None:
    """A Hessian on [ML, link H, point charge] maps onto [ML, point charge, host]."""
    from ase import Atoms

    from mlmm.backends.mlmm_calc import MLMMCore

    rng = np.random.default_rng(1)
    A = rng.normal(size=(9, 9))
    A = A + A.T
    b = rng.normal(size=9)
    param = 0.7 if method == "scaled" else 1.09

    core = MLMMCore.__new__(MLMMCore)
    core.selection_indices = [0]
    core._idx_map_real_to_model = {0: 0}
    core.link_atom_method = method

    def link_position(q, m):
        if method == "scaled":
            return q + param * (m - q)
        return q + param * (m - q) / np.linalg.norm(m - q)

    def corr_coords(x):  # x rows: [ML, point charge, host]
        q, p, m = x
        return np.concatenate([q, link_position(q, m), p])

    def energy(x):
        y = corr_coords(x)
        return b @ y + 0.5 * y @ A @ y

    x0 = np.array([[0.0, 0.0, 0.0], [4.0, 1.0, 0.5], [1.5, 0.3, -0.2]])
    y0 = corr_coords(x0)
    dF = -(b + A @ y0).reshape(3, 3)
    atoms_real = Atoms("CHC", positions=[x0[0], x0[2], x0[1]])  # ML, host, point charge
    atoms_model_lh = Atoms("CH", positions=[x0[0], y0[3:6]])
    link_maps = core._link_jacobians([(1, 0, 1, param)], atoms_model_lh, atoms_real)
    H_real, basis = core._link_hessian_to_real(A, dF, link_maps, atoms_model_lh, atoms_real, [2])

    assert basis == [0, 2, 1]

    def forces(x):
        return _fd_forces(energy, x, 1.0e-6)

    step = 1.0e-4
    fd_hess = np.column_stack([
        -(forces((x0.reshape(-1) + e).reshape(3, 3)) - forces((x0.reshape(-1) - e).reshape(3, 3)))
        / (2 * step)
        for e in np.eye(9) * step
    ])
    assert np.allclose(H_real, fd_hess, atol=2.0e-3)


def test_xtb_correction_runs_on_ml_plus_link_hydrogens() -> None:
    from test_mlmm_result_constraints import _force_only_core

    def run(link_force):
        core = _force_only_core(with_link=True)
        core.embedcharge = True
        core.embedcharge_cutoff = None
        core._embedding_near_mm = set()
        core.model_charge = 0
        core.model_mult = 1
        seen = {}

        class Recorder:
            def compute_correction(self, **kwargs):
                seen.update(kwargs)
                n = len(kwargs["symbols"]) + len(kwargs["mm_charges"])
                dF = np.zeros((n, 3))
                dF[1] = link_force
                return 0.0, dF, None

        core._embed_correction = Recorder()
        core._get_mm_charges = MethodType(
            lambda self, indices: np.zeros(len(indices), dtype=float), core
        )
        result = core.compute(
            np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]]), return_forces=True
        )
        return seen, result

    seen, base = run([0.0, 0.0, 0.0])
    assert seen["symbols"] == ["H", "H"]  # the ML atom and its link hydrogen
    assert len(seen["coords_ml_ang"]) == 2
    _, pushed = run([2.0, 0.0, 0.0])
    # g = 0.5 in the fake core: half of the link-H force lands on the host.
    assert np.allclose(pushed["forces"][1] - base["forces"][1], [1.0, 0.0, 0.0])


def test_xtb_correction_minus_model_low_coulomb_removes_the_first_order_term(tmp_path) -> None:
    xtb = shutil.which("xtb")
    if xtb is None:
        pytest.skip("xtb is not installed")
    from pysisyphus.constants import AU2EV

    from mlmm.backends.mlmm_calc import point_charge_coulomb
    from mlmm.backends.xtb_embedcharge_correction import delta_embedcharge_minus_noembed

    for distance in (10.0, 20.0):
        d_e, _, _ = delta_embedcharge_minus_noembed(
            symbols=["Li"], coords_q_ang=np.zeros((1, 3)),
            mm_coords_ang=np.array([[distance, 0.0, 0.0]]), mm_charges=np.array([-1.0]),
            charge=1, multiplicity=1, need_forces=False, need_hessian=False,
            xtb_cmd=xtb, xtb_workdir=str(tmp_path), ncores=1,
        )
        e_coul, _, _ = point_charge_coulomb(
            np.zeros((1, 3)), [1.0], np.array([[distance, 0.0, 0.0]]), [-1.0]
        )
        assert e_coul == pytest.approx(-AU2EV * BOHR / distance, rel=1.0e-5)
        # xTB damps the Coulomb kernel at short range; the residual is far
        # below the -1/R term that a double count would leave.
        assert abs(d_e - e_coul) < 0.02 * abs(e_coul)
