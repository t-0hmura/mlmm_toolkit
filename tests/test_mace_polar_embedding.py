"""Adapter lifetime and mixed ML/MM derivatives without optional model downloads."""
import sys
from types import SimpleNamespace

import numpy as np
import pytest
import torch

from mlmm.backends.mace_embedcharge import MACEPolarEmbedding


class Source(torch.nn.Module):
    def __init__(self, scale=1.):
        super().__init__()
        self.scale, self.external = scale, None
    @classmethod
    def from_energy(cls, energy):
        return cls()
    @classmethod
    def from_features(cls, features, external_scale):
        return cls(external_scale)
    def set_external_sources(self, features, positions, batch):
        self.external = (features, positions, batch)
    def clear_external_sources(self):
        self.external = None


class PolarMACE(torch.nn.Module):
    atomic_multipoles_max_l = 1
    def __init__(self):
        super().__init__()
        self.parameter = torch.nn.Parameter(torch.tensor(1., dtype=torch.float64))
        self.coulomb_energy = Source()
        self.electric_potential_descriptor = Source()
        self.fail = False
    def forward(self, graph, **kwargs):
        if self.fail:
            raise RuntimeError("model failure")
        x = graph["positions"]
        energy = .2*(x*x).sum()
        source = self.coulomb_energy.external
        if source is not None:
            features, positions, _ = source
            radius = torch.linalg.vector_norm(x[:,None]-positions[None,:], dim=2)
            energy = energy + (features[:,0]/radius).sum()
            assert self.electric_potential_descriptor.scale == .5
            assert self.electric_potential_descriptor.external is not None
        return {"energy": energy.reshape(1)}


@pytest.fixture
def adapter(monkeypatch):
    monkeypatch.setitem(sys.modules, "graph_longrange.external_source_energy",
                         SimpleNamespace(GTOElectrostaticExternalSourceEnergy=Source))
    monkeypatch.setitem(sys.modules, "graph_longrange.external_source_features",
                         SimpleNamespace(GTOElectrostaticExternalSourceFeatures=Source))
    model = PolarMACE()
    calculator = SimpleNamespace(models=[model], _atoms_to_batch=lambda atoms:
                                  SimpleNamespace(to_dict=lambda: {}))
    return MACEPolarEmbedding(calculator)


def test_mixed_hessian_reaction_forces_and_active_mm_mapping(adapter):
    symbols = ["O", "H"]
    x = np.array([[.1,.2,.3], [.9,.1,.2], [3.,1.,.4], [-2.,2.,.1]])
    charges = np.array([.3,-.2])
    def evaluate(coords, **kw):
        return adapter.compute_correction(symbols, coords[:2], coords[2:], charges,
                                            0, 1, need_forces=True, **kw)
    e, f, h = evaluate(x, need_hessian=True)
    fd = np.zeros_like(x)
    step = 1e-5
    for i in range(x.size):
        plus, minus = x.copy(), x.copy()
        plus.flat[i] += step
        minus.flat[i] -= step
        fd.flat[i] = -(evaluate(plus)[0]-evaluate(minus)[0])/(2*step)
    np.testing.assert_allclose(f, fd, atol=1e-9)
    np.testing.assert_allclose(f.sum(0), 0, atol=1e-12)
    assert np.any(np.abs(h[:6,6:]) > 1e-4)
    active = evaluate(x, need_hessian=True, hessian_mm_rows=[1])[2]
    dofs = [0,1,2,3,4,5,9,10,11]
    np.testing.assert_allclose(active, h[np.ix_(dofs,dofs)], atol=1e-12)
    assert e != 0


def test_zero_charges_and_translation_reproduce_vacuum(adapter):
    x = np.array([[.1,.2,.3], [3.,1.,.4]])
    args = (["He"], x[:1], x[1:], [0.], 0, 1)
    e, f, h = adapter.compute_correction(*args, need_hessian=True)
    assert e == 0
    np.testing.assert_array_equal(f, np.zeros((2,3)))
    np.testing.assert_array_equal(h, np.zeros((6,6)))
    e1, f1, _ = adapter.compute_correction(["He"],x[:1],x[1:],[.3],0,1,need_forces=True)
    e2, f2, _ = adapter.compute_correction(["He"],x[:1]+100,x[1:]+100,[.3],0,1,need_forces=True)
    assert e1 == pytest.approx(e2, abs=1e-12)
    np.testing.assert_allclose(f1,f2,atol=1e-12)


def test_modules_restored_and_sources_cleared_on_success_and_failure(adapter):
    originals = adapter.model.coulomb_energy, adapter.model.electric_potential_descriptor
    args = (["He"],np.zeros((1,3)),np.ones((1,3)),[.3],0,1)
    adapter.compute_correction(*args)
    assert (adapter.model.coulomb_energy,adapter.model.electric_potential_descriptor)==originals
    assert adapter._energy.external is adapter._features.external is None
    adapter.model.fail = True
    with pytest.raises(RuntimeError,match="model failure"):
        adapter.compute_correction(*args)
    assert (adapter.model.coulomb_energy,adapter.model.electric_potential_descriptor)==originals
    assert adapter._energy.external is adapter._features.external is None


def test_empty_mm_region_is_exact_zero(adapter):
    e, f, h = adapter.compute_correction(["He"],np.zeros((1,3)),[],[],0,1,need_hessian=True)
    assert e == 0
    np.testing.assert_array_equal(f,np.zeros((1,3)))
    np.testing.assert_array_equal(h,np.zeros((3,3)))
