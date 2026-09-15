"""Opt-in Cartesian atom restriction; physical axes, mappings and actual steps.

All PESs/matrices here are analytic small fixtures. Molecular acceptance and
algorithmic speed are not inferred from these contracts.
"""
from types import SimpleNamespace

import numpy as np
import pytest
import torch

from pysisyphus.Geometry import Geometry
from pysisyphus.calculators.Calculator import Calculator
from pysisyphus.optimizers.RFOptimizer import RFOptimizer
from pysisyphus.tsoptimizers.RSIRFOptimizer import RSIRFOptimizer
from pysisyphus.tsoptimizers.RSPRFOptimizer import RSPRFOptimizer


def make_opt(tmp_path, *, norm="max_atom", radius=.03, partial=False, **extra):
    geometry = Geometry(
        ["H"] * (4 if partial else 2),
        np.array([[0., 0., 0.], [2., 0., 0.], [0., 2., 0.], [0., 0., 2.]])
        [:4 if partial else 2].ravel(),
        coord_type="cart", freeze_atoms=[1, 3] if partial else [],
    )
    options = dict(
        roots=[0], verify_saddle=False, flatten_enabled=False,
        hessian_init="unit", trust_radius=radius, trust_min=1e-8, trust_max=2.,
        max_micro_cycles=50, min_line_search=False, max_line_search=False,
        out_dir=tmp_path, dump=False, max_cycles=1,
    )
    if norm is not None:
        options["trust_norm"] = norm
    options.update(extra)
    optimizer = RSPRFOptimizer(geometry, **options)
    if partial:
        geometry.within_partial_hessian = dict(
            active_atoms=np.array([2, 0]), active_dofs=np.array([6, 7, 8, 0, 1, 2]),
            active_n_dof=6, full_n_dof=12,
        )
        optimizer._set_active_dofs(True)
    return optimizer


def proposal(optimizer, monkeypatch, *, rotation=None, line_step=None):
    # An index-one physical H with distributed gradient and a known root.
    lam = np.array([-1., 1., 1., 1., 1., 1.])
    vectors = np.eye(6) if rotation is None else rotation
    gradient = vectors @ np.ones(6)
    hessian = vectors @ np.diag(lam) @ vectors.T
    optimizer.H = optimizer.cur_H = hessian.copy()
    optimizer.forces = [-optimizer.full_from_active(gradient)]
    monkeypatch.setattr(
        optimizer, "housekeeping",
        lambda: (0., gradient.copy(), hessian.copy(), lam.copy(), vectors.copy(), False),
    )
    if line_step is not None:
        monkeypatch.setattr(
            optimizer, "step_and_grad_from_line_search",
            lambda energy, g, *args: (np.array(line_step), g),
        )
    step = optimizer.optimize()
    assert optimizer.geometry.calculator is None
    return step, gradient, hessian


def atom_norm(step):
    return float(np.linalg.norm(np.asarray(step).reshape(-1, 3), axis=1).max())


def test_native_prfo_distributed_cartesian_bound_and_default_parity(tmp_path, monkeypatch):
    atomic = make_opt(tmp_path / "atomic")
    default = make_opt(tmp_path / "default", norm=None)
    explicit = make_opt(tmp_path / "explicit", norm="l2")
    sa, _, _ = proposal(atomic, monkeypatch)
    sd, _, _ = proposal(default, monkeypatch)
    se, _, _ = proposal(explicit, monkeypatch)
    np.testing.assert_array_equal(sd, se)
    assert np.linalg.norm(sd) <= .03 * (1 + 1e-12)
    assert atom_norm(sa) == pytest.approx(.03, rel=1e-8)
    assert np.linalg.norm(sa) > 1.1 * .03
    assert atomic._last_atomic_trust["evaluations"] <= atomic.max_micro_cycles
    # The physical Hessian/root and sign-specific PRFO motion remain intact.
    np.testing.assert_array_equal(atomic.cur_H, np.diag([-1., 1., 1., 1., 1., 1.]))
    assert sa[0] > 0. and np.all(sa[1:] < 0.)


@pytest.mark.parametrize("partial", [False, True])
def test_native_prfo_atom_rotation_invariance(tmp_path, monkeypatch, partial):
    theta = .713
    rotation3 = np.array([[np.cos(theta), -np.sin(theta), 0.],
                          [np.sin(theta), np.cos(theta), 0.], [0., 0., 1.]])
    rotation = np.kron(np.eye(2), rotation3)
    first = make_opt(tmp_path / "a", partial=partial)
    second = make_opt(tmp_path / "b", partial=partial)
    original, _, _ = proposal(first, monkeypatch)
    rotated, _, _ = proposal(second, monkeypatch, rotation=rotation)
    np.testing.assert_allclose(
        rotated.reshape(-1, 3), original.reshape(-1, 3) @ rotation3.T,
        rtol=1e-8, atol=1e-11,
    )
    if partial:
        np.testing.assert_array_equal(rotated.reshape(4, 3)[[1, 3]], 0.)
        np.testing.assert_array_equal(first.active_dof_indices, [6, 7, 8, 0, 1, 2])
    assert atom_norm(rotated) <= .03 * (1 + 1e-12)


def test_combined_line_search_contribution_uses_same_cartesian_bound(tmp_path, monkeypatch):
    optimizer = make_opt(tmp_path)
    line = np.array([.004, 0., 0., 0., .002, 0.])
    actual, gradient, hessian = proposal(optimizer, monkeypatch, line_step=line)
    assert atom_norm(actual) == pytest.approx(.03, rel=1e-8)
    expected = (gradient @ actual + .5 * actual @ hessian @ actual) / (1 + actual @ actual)
    assert optimizer.predicted_energy_changes[-1] == pytest.approx(expected)


@pytest.mark.parametrize("backend", ["numpy", "torch"])
def test_partial_mapping_and_final_guard_preserve_distributed_step(tmp_path, backend):
    optimizer = make_opt(tmp_path, radius=.05, partial=True)
    compact = np.array([.03, .04, 0., 0., .03, .04])
    value = torch.tensor(compact) if backend == "torch" else compact.copy()
    bounded = optimizer._bound_to_trust_radius(value)
    observed = bounded.numpy() if backend == "torch" else bounded
    np.testing.assert_array_equal(observed, compact)
    assert np.linalg.norm(observed) > optimizer.trust_radius
    np.testing.assert_allclose(
        optimizer.full_from_active(observed), [0., .03, .04, 0., 0., 0., .03, .04, 0., 0., 0., 0.]
    )
    too_large = optimizer._bound_to_trust_radius(2 * value)
    np.testing.assert_allclose(too_large.numpy() if backend == "torch" else too_large, compact)


def test_actual_full_step_feedback_uses_atom_norm_under_partial_mapping(tmp_path):
    optimizer = make_opt(tmp_path, radius=.05, partial=True)
    # Stored steps are actual full displacements, not compact proposals.
    full = np.array([0., .03, .04, 0., 0., 0., .04, 0., 0., 0., 0., 0.])
    optimizer.steps = [full.copy()]
    optimizer.forces = [np.zeros(12), np.zeros(12)]
    optimizer.energies = [0., -1.]
    optimizer.predicted_energy_changes = [-1.]
    assert optimizer._trust_step_norm(full, full=True) == pytest.approx(.05)
    optimizer.update_trust_radius()
    assert optimizer.trust_radius == pytest.approx(.1)
    np.testing.assert_array_equal(optimizer.steps[0], full)


def test_cartesian_wrapper_frozen_coordinate_mapping_and_feedback(tmp_path, monkeypatch):
    optimizer = make_opt(tmp_path)
    optimizer.geometry = Geometry(
        ["H"] * 4, np.arange(12, dtype=float), coord_type="cartesian", freeze_atoms=[1, 3]
    )
    # The wrapper's native working basis excludes freezes without a Hessian map.
    assert not optimizer.using_active_dofs and optimizer.geometry.coords.size == 6
    step, _, _ = proposal(optimizer, monkeypatch)
    assert step.size == 6 and atom_norm(step) <= .03 * (1 + 1e-12)
    full = optimizer.geometry.internal.transform_int_step(step, pure=True)
    np.testing.assert_array_equal(full.reshape(4, 3)[[1, 3]], 0.)
    assert optimizer._trust_step_norm(step) == pytest.approx(atom_norm(full))
    optimizer.steps = [step]
    optimizer.forces.append(optimizer.forces[-1].copy())
    optimizer.energies = [0., -1.]
    optimizer.predicted_energy_changes = [-1.]
    optimizer.update_trust_radius()
    assert optimizer.trust_radius == pytest.approx(.06)


class AnalyticSaddle(Calculator):
    def __init__(self, reference, out_dir):
        super().__init__(out_dir=out_dir, check_mem=False)
        self.reference = np.array(reference)
        self.lam = np.array([-1., 1., 1., 1., 1., 1.])
        self.hessian_calls = 0

    def get_forces(self, atoms, coords, **kwargs):
        delta = np.asarray(coords) - self.reference
        return dict(energy=float(delta.sum() + .5 * delta @ (self.lam * delta)),
                    forces=-(np.ones(6) + self.lam * delta))

    get_energy = get_forces

    def get_hessian(self, atoms, coords, **kwargs):
        self.hessian_calls += 1
        return dict(self.get_forces(atoms, coords), hessian=np.diag(self.lam))


def test_native_run_applies_distributed_steps_without_extra_hessian_calls(tmp_path):
    optimizer = make_opt(tmp_path, hessian_init="calc", trust_update=False,
                         max_cycles=3)
    initial = optimizer.geometry.coords.copy()
    calculator = AnalyticSaddle(initial, tmp_path)
    optimizer.geometry.set_calculator(calculator)
    optimizer.run()
    assert calculator.hessian_calls == 1
    assert len(optimizer.steps) == 3 and not optimizer.is_converged
    points = [*optimizer.coords, optimizer.geometry.coords]
    for index, step in enumerate(optimizer.steps):
        np.testing.assert_array_equal(step, points[index + 1] - points[index])
        assert atom_norm(step) <= .03 * (1 + 1e-12)
        assert np.linalg.norm(step) > .03


def test_scalar_cap_retains_feasible_trial_and_does_not_commit_trial_root_history(tmp_path, monkeypatch):
    optimizer = make_opt(tmp_path, radius=.8, max_micro_cycles=2, rfo_overlaps=True)
    incoming = (np.array([7., 8.]), np.arange(6, dtype=float))
    optimizer.prev_eigvec_max, optimizer.prev_eigvec_min = incoming
    seen = []

    def family(values, gradient, alpha, kind, prev_eigvec):
        seen.append((kind, prev_eigvec))
        # A known smooth scalar family; exactly two evaluations bracket .8.
        step = np.zeros(len(values))
        if kind == "max":
            step[0] = 1 / np.sqrt(alpha)
        return step, 0., 1., np.full(len(values) + 1, alpha)

    monkeypatch.setattr(optimizer, "solve_rfo_secular", family)
    step = optimizer._max_atom_prfo_step(
        np.ones(6), np.eye(6), np.ones(6), np.zeros(6), [0], [1, 2, 3, 4, 5]
    )
    assert atom_norm(step) == pytest.approx(1 / np.sqrt(2))
    assert optimizer._last_atomic_trust["termination"] == "feasible_budget"
    assert optimizer._last_atomic_trust["evaluations"] == 2
    assert all(previous is incoming[0 if kind == "max" else 1] for kind, previous in seen)
    np.testing.assert_array_equal(optimizer.prev_eigvec_max, [2., 2.])


def test_scalar_cap_without_any_feasible_trial_raises(tmp_path):
    optimizer = make_opt(tmp_path, radius=1e-10, max_micro_cycles=1)
    with pytest.raises(ValueError, match="without a feasible point"):
        optimizer._max_atom_prfo_step(
            np.array([-1., 1., 1., 1., 1., 1.]), np.eye(6), np.ones(6), np.zeros(6),
            [0], [1, 2, 3, 4, 5],
        )


def test_full_size_permutation_is_distinct_from_full_actual_coordinates(tmp_path):
    optimizer = make_opt(tmp_path, radius=.1)
    optimizer._using_active_dofs = True
    optimizer._active_dof_indices = np.array([0, 3, 1, 4, 2, 5])
    step = np.array([.03, .04, 0., 0., 0., 0.])
    assert optimizer._trust_step_norm(step) == pytest.approx(.04)
    assert optimizer._trust_step_norm(step, full=True) == pytest.approx(.05)
    optimizer._active_dof_indices[1] = 0
    with pytest.raises(ValueError, match="ordered active DOFs"):
        optimizer._trust_step_norm(step)


def test_image_quadratic_substep_is_explicitly_conservative_and_bounded(tmp_path, monkeypatch):
    optimizer = make_opt(tmp_path, partial=True)
    # Reflecting root0 leaves root1 negative: this exercises the hard case.
    lam = np.array([-.1, -.05, 1., 2., 3., 4.])
    optimizer.H = optimizer.cur_H = np.diag(lam)
    optimizer.forces = [np.zeros(12)]
    monkeypatch.setattr(optimizer, "_hessian_system", lambda *args: (np.zeros(6), optimizer.H, lam, np.eye(6)))
    messages = []
    monkeypatch.setattr(optimizer, "log", messages.append)
    step, _ = optimizer._image_trust_step()
    assert np.linalg.norm(step) <= .03 * (1 + 1e-12)
    assert atom_norm(optimizer.full_from_active(step)) <= .03 * (1 + 1e-12)
    assert optimizer.atomic_trust_conservative_quadratic_steps == 1
    assert any("conservative L2" in message for message in messages)
    assert abs(step[1]) == pytest.approx(.03)


@pytest.mark.parametrize("coord_type", ["redund", "dlc", "mwcartesian"])
def test_unsupported_coordinates_fail_explicitly(tmp_path, coord_type):
    with pytest.raises(ValueError, match="Cartesian coordinates"):
        RSPRFOptimizer(SimpleNamespace(coord_type=coord_type), trust_norm="max_atom", out_dir=tmp_path)


@pytest.mark.parametrize("cls", [RFOptimizer, RSIRFOptimizer])
def test_other_step_owners_do_not_accept_an_ineffective_option(tmp_path, cls):
    geom = Geometry(["H"], [0., 0., 0.], coord_type="cart")
    with pytest.raises(ValueError, match="only by RS-PRFO"):
        cls(geom, trust_norm="max_atom", out_dir=tmp_path)


def test_yaml_constructor_plumbing_and_alias(tmp_path):
    from mlmm.workflows.tsopt import _build_rsirfo_kwargs
    kwargs = _build_rsirfo_kwargs(
        dict(trust_norm="max_atom", trust_radius=.188972612546,
             trust_max=.188972612546, trust_min=1e-4, hessian_init="unit"),
        max_cycles=1, out_dir=tmp_path, mode="rsprfo", dump=False,
    )
    geometry = Geometry(["H", "H"], [0., 0., 0., 2., 0., 0.], coord_type="cartesian")
    optimizer = RSPRFOptimizer(geometry, **kwargs)
    assert optimizer.trust_norm == "max_atom"
    assert optimizer.trust_radius == .188972612546


@pytest.mark.parametrize("stored", [None, "l2", "unknown"])
def test_max_atom_restart_rejects_legacy_or_mismatch_before_history_mutation(tmp_path, stored):
    optimizer = make_opt(tmp_path)
    original = optimizer.coords
    restart = {} if stored is None else dict(trust_norm=stored)
    with pytest.raises(ValueError, match="Restart trust_norm"):
        optimizer.set_restart_info(restart)
    assert optimizer.coords is original


def test_restart_norm_round_trip_and_reverse_mismatch(tmp_path):
    atomic = make_opt(tmp_path / "atomic")
    atomic.H = np.eye(6)
    saved = atomic._get_opt_restart_info()
    assert saved["trust_norm"] == "max_atom"
    restored = make_opt(tmp_path / "restored", radius=.02)
    restored._set_opt_restart_info(saved)
    assert restored.trust_radius == .03
    default = make_opt(tmp_path / "default", norm=None)
    with pytest.raises(ValueError, match="Restart trust_norm"):
        default.set_restart_info(saved)
    default.H = np.eye(6)
    legacy = default._get_opt_restart_info()
    legacy.pop("trust_norm")
    default._set_opt_restart_info(legacy)
    assert default.trust_norm == "l2"
