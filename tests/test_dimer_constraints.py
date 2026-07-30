import numpy as np
import pytest
import torch

from pysisyphus.calculators.Dimer import Dimer
from pysisyphus.tr_projection import full_cartesian_tr_basis


COORDS = np.array(
    [
        [0.0, 0.0, 0.0],
        [1.2, 0.0, 0.0],
        [0.1, 1.1, 0.0],
        [0.2, 0.3, 1.3],
        [1.0, 0.8, 0.7],
    ],
    dtype=float,
)
MASSES = torch.tensor([12.0, 1.0, 16.0, 14.0, 32.0], dtype=torch.float64)


class RecordingCalculator:
    def __init__(self):
        self.force_coords = []

    def get_forces(self, atoms, coords):
        coords = np.asarray(coords, dtype=float).copy()
        self.force_coords.append(coords)
        forces = np.linspace(-0.5, 0.5, coords.size)
        return {"energy": 0.0, "forces": forces}

    def get_energy(self, atoms, coords):
        return {"energy": 0.0}


def _constraint_data(frozen):
    return _constraint_data_at(COORDS, frozen)


def _constraint_data_at(coords, frozen):
    active = [atom for atom in range(len(COORDS)) if atom not in frozen]
    return full_cartesian_tr_basis(
        torch.as_tensor(coords), MASSES, active, mode="constrained"
    )


def _assert_allowed(vector, frozen, basis, *, atol=2.0e-12):
    vector = np.asarray(vector).reshape(-1)
    frozen_dofs = [3 * atom + axis for atom in frozen for axis in range(3)]
    if frozen_dofs:
        np.testing.assert_array_equal(vector[frozen_dofs], 0.0)
    if basis.shape[1]:
        np.testing.assert_allclose(basis.T @ vector, 0.0, atol=atol)


def _assert_random_states_equal(left, right):
    assert left[0] == right[0]
    np.testing.assert_array_equal(left[1], right[1])
    assert left[2:] == right[2:]


@pytest.mark.parametrize("seed", [1729, None])
def test_random_orientation_does_not_mutate_global_numpy_state(tmp_path, seed):
    original_state = np.random.get_state()
    try:
        np.random.seed(314159)
        state = np.random.get_state()
        dimer = Dimer(
            calculator=None,
            write_orientations=False,
            seed=seed,
            out_dir=tmp_path,
        )
        dimer.set_N_raw(COORDS.reshape(-1))
        _assert_random_states_equal(np.random.get_state(), state)
    finally:
        np.random.set_state(original_state)


def test_random_orientation_is_reproducible_for_explicit_seed(tmp_path):
    dimers = [
        Dimer(
            calculator=None,
            write_orientations=False,
            seed=8675309,
            out_dir=tmp_path / str(index),
        )
        for index in range(2)
    ]
    for dimer in dimers:
        dimer.set_N_raw(COORDS.reshape(-1))

    np.testing.assert_array_equal(dimers[0].N_raw, dimers[1].N_raw)
    np.testing.assert_array_equal(dimers[0].N, dimers[1].N)


def test_explicit_seed_retains_legacy_randomstate_sequence(tmp_path):
    seed = 24680
    raw = np.random.RandomState(seed).rand(COORDS.size)
    expected = raw.reshape(-1, 3)
    expected = (expected - expected.mean(axis=0)).reshape(-1)
    expected /= np.linalg.norm(expected)

    dimer = Dimer(
        calculator=None,
        write_orientations=False,
        seed=seed,
        out_dir=tmp_path,
    )
    dimer.set_N_raw(COORDS.reshape(-1))

    np.testing.assert_array_equal(dimer.N_raw, expected)
    np.testing.assert_array_equal(dimer.N, expected)


def test_invalid_rotation_method_fails_before_output_setup(tmp_path):
    output_dir = tmp_path / "must-not-exist"
    with pytest.raises(ValueError) as exc_info:
        Dimer(
            calculator=None,
            rotation_method="invalid",
            write_orientations=False,
            seed=0,
            out_dir=output_dir,
        )

    message = str(exc_info.value)
    assert "invalid" in message
    assert "direct" in message
    assert "fourier" in message
    assert not output_dir.exists()


@pytest.mark.parametrize(
    ("frozen", "expected_rank"),
    [([], 6), ([0], 3), ([0, 1], 1), ([0, 1, 2], 0)],
)
def test_dimer_stays_on_constraint_manifold_for_anchor_ranks(
    tmp_path, frozen, expected_rank
):
    basis_t, info = _constraint_data(frozen)
    basis = basis_t.numpy()
    assert info.effective_rank == expected_rank

    recorder = RecordingCalculator()
    seed = np.arange(1, COORDS.size + 1, dtype=float)
    dimer = Dimer(
        calculator=recorder,
        N_raw=seed,
        frozen_atoms=frozen,
        rigid_basis=basis,
        write_orientations=False,
        rotation_disable=True,
        seed=0,
        out_dir=tmp_path,
    )
    _assert_allowed(dimer.N, frozen, basis)
    np.testing.assert_allclose(np.linalg.norm(dimer.N), 1.0, atol=1.0e-14)

    center = COORDS.reshape(-1).copy()
    dimer.atoms = ("C", "H", "O", "N", "S")
    dimer.coords0 = center
    endpoint = dimer.coords1
    _assert_allowed(endpoint - center, frozen, basis)
    if frozen:
        np.testing.assert_array_equal(
            endpoint.reshape(-1, 3)[frozen], COORDS[frozen]
        )

    _ = dimer.f1
    evaluated = recorder.force_coords[-1]
    if frozen:
        np.testing.assert_array_equal(
            evaluated.reshape(-1, 3)[frozen], COORDS[frozen]
        )

    dimer._f0 = np.linspace(0.4, -0.2, COORDS.size)
    dimer._f1 = np.linspace(-0.3, 0.6, COORDS.size)
    _assert_allowed(dimer.rot_force, frozen, basis)

    theta = np.linspace(-0.7, 0.9, COORDS.size)
    trial = dimer.rotate_coords1(0.31, theta)
    _assert_allowed(trial - center, frozen, basis)
    if frozen:
        np.testing.assert_array_equal(trial.reshape(-1, 3)[frozen], COORDS[frozen])

    dimer.restrict_step = lambda step: step
    dimer.direct_rotation(lambda force, previous: theta, None)
    _assert_allowed(dimer.N, frozen, basis)
    if frozen:
        np.testing.assert_array_equal(
            dimer.coords1.reshape(-1, 3)[frozen], COORDS[frozen]
        )


def test_positive_control_reproduces_legacy_frozen_orientation_bug(tmp_path):
    seed = np.zeros_like(COORDS)
    seed[1, 0] = 1.0
    seed[2, 1] = 1.0

    legacy = Dimer(
        calculator=None,
        N_raw=seed,
        write_orientations=False,
        seed=0,
        out_dir=tmp_path / "legacy",
    )
    assert np.linalg.norm(legacy.N[:3]) > 0.1

    basis, info = _constraint_data([0])
    constrained = Dimer(
        calculator=None,
        N_raw=seed,
        frozen_atoms=[0],
        rigid_basis=basis.numpy(),
        write_orientations=False,
        seed=0,
        out_dir=tmp_path / "constrained",
    )
    assert info.effective_rank == 3
    np.testing.assert_array_equal(constrained.N[:3], 0.0)
    np.testing.assert_allclose(basis.numpy().T @ constrained.N, 0.0, atol=2.0e-12)


def test_dimer_refreshes_rigid_basis_when_center_moves(tmp_path):
    frozen = [0]
    initial_basis, _ = _constraint_data_at(COORDS, frozen)
    moved = COORDS.copy()
    moved[1:] += np.array(
        [[0.4, 0.2, -0.3], [-0.2, 0.3, 0.5], [0.3, -0.4, 0.2], [-0.5, 0.1, -0.2]]
    )
    moved_basis, _ = _constraint_data_at(moved, frozen)
    initial_basis = initial_basis.numpy()
    moved_basis = moved_basis.numpy()

    raw = np.random.default_rng(4).normal(size=COORDS.size)
    stale = raw - initial_basis @ (initial_basis.T @ raw)
    stale /= np.linalg.norm(stale)
    assert np.linalg.norm(moved_basis.T @ stale) > 0.25

    seen_centers = []

    def basis_at(coords_flat):
        seen_centers.append(np.asarray(coords_flat).copy())
        basis, _ = _constraint_data_at(np.asarray(coords_flat).reshape(-1, 3), frozen)
        return basis.numpy()

    dimer = Dimer(
        calculator=None,
        N_raw=stale,
        frozen_atoms=frozen,
        rigid_basis=initial_basis,
        rigid_basis_getter=basis_at,
        write_orientations=False,
        seed=0,
        out_dir=tmp_path,
    )
    dimer.coords0 = COORDS.reshape(-1)
    dimer.N = stale
    _ = dimer.coords1
    assert len(seen_centers) == 1

    stale = dimer.N.copy()
    dimer.coords0 = moved.reshape(-1)
    assert len(seen_centers) == 2
    assert np.linalg.norm(moved_basis.T @ stale) > 0.25
    np.testing.assert_allclose(moved_basis.T @ dimer.N, 0.0, atol=2.0e-12)
    np.testing.assert_array_equal(dimer.N[:3], 0.0)


def test_unconstrained_dimer_retains_legacy_translation_removal(tmp_path):
    seed = np.arange(1, COORDS.size + 1, dtype=float)
    expected = seed.reshape(-1, 3)
    expected = (expected - expected.mean(axis=0)).reshape(-1)
    expected /= np.linalg.norm(expected)

    dimer = Dimer(
        calculator=None,
        N_raw=seed,
        write_orientations=False,
        seed=0,
        out_dir=tmp_path,
    )
    np.testing.assert_allclose(dimer.N, expected, rtol=0.0, atol=0.0)


def test_nearly_pure_rigid_orientation_is_rejected(tmp_path):
    basis, _ = _constraint_data([])
    rigid = basis[:, 0].numpy()
    allowed = np.random.default_rng(17).normal(size=rigid.size)
    allowed -= basis.numpy() @ (basis.numpy().T @ allowed)
    allowed /= np.linalg.norm(allowed)

    with pytest.raises(ValueError, match="orientation vanishes"):
        Dimer(
            calculator=None,
            N_raw=rigid + 1.0e-15 * allowed,
            rigid_basis=basis.numpy(),
            write_orientations=False,
            seed=0,
            out_dir=tmp_path,
        )


def test_rotation_disabled_can_evaluate_energy_before_force_cache(tmp_path):
    recorder = RecordingCalculator()
    dimer = Dimer(
        calculator=recorder,
        N_raw=np.arange(1, COORDS.size + 1, dtype=float),
        rotation_disable=True,
        rotation_disable_pos_curv=False,
        write_orientations=False,
        seed=0,
        out_dir=tmp_path,
    )

    result = dimer.get_forces(("C", "H", "O", "N", "S"), COORDS.reshape(-1))
    assert result["energy"] == 0.0
    assert np.all(np.isfinite(result["forces"]))


def test_hessian_dimer_passes_geometry_constraints_to_dimer(tmp_path, monkeypatch):
    import mlmm.workflows.tsopt as tsopt

    captured = {}

    class FakeGeometry:
        atomic_numbers = [6, 1, 8, 7, 16]
        cart_coords = COORDS.reshape(-1).copy()
        freeze_atoms = np.array([0], dtype=int)

        def set_calculator(self, calculator):
            self.calculator = calculator

    class FakeDimer:
        def __init__(self, **kwargs):
            captured.update(kwargs)

    class FakeLBFGS:
        cur_cycle = 0
        is_converged = True

        def __init__(self, geometry, **kwargs):
            self.geometry = geometry

        def run(self):
            return None

    runner = object.__new__(tsopt.HessianDimer)
    runner.geom = FakeGeometry()
    runner.calc_kwargs = {}
    runner.masses_au_t = MASSES
    runner.freeze_atoms = [0]
    runner.tr_projection = "constrained"
    runner.rigid_projection_info = {}
    runner.dimer_kwargs = {}
    runner.lbfgs_kwargs = {}
    runner.mode_path = tmp_path / "mode.dat"
    runner.mem = 100
    runner.out_dir = tmp_path
    runner.dump = False
    runner.optim_all_path = tmp_path / "all.xyz"

    monkeypatch.setattr(tsopt, "mlmm", lambda **kwargs: object())
    monkeypatch.setattr(tsopt, "Dimer", FakeDimer)
    monkeypatch.setattr(tsopt, "LBFGS", FakeLBFGS)

    steps, converged = runner._dimer_segment("baker", 2)
    basis = captured["rigid_basis"]
    assert steps == 1
    assert converged is True
    assert captured["frozen_atoms"] == [0]
    assert captured["out_dir"] == str(tmp_path)
    assert callable(captured["rigid_basis_getter"])
    assert basis.shape == (COORDS.size, 3)
    np.testing.assert_array_equal(basis[:3], 0.0)
    np.testing.assert_allclose(basis.T @ basis, np.eye(3), atol=2.0e-12)
    assert runner.rigid_projection_info["effective_rank"] == 3
