"""Pins the default models and the routes for the newly supported models.

Defaults follow each backend's recommendation (UMA uma-s-1p2p1, ORB orbmol-v2,
MACE MACE-OMOL-0, AIMNet2 aimnet2). ORB resolves names through orb-models'
ORB_PRETRAINED_MODELS and accepts only conservative regressors; MACE-POLAR-1
loads through mace_polar. None of these tests needs a GPU, a network or a
backend package.
"""

from __future__ import annotations

import inspect
import sys
from types import ModuleType, SimpleNamespace

import numpy as np
import pytest
import torch
from ase import Atoms
from ase.calculators.calculator import Calculator, all_changes

from mlmm.backends import mlmm_calc
from mlmm.core import defaults


def test_default_models() -> None:
    assert defaults.DEFAULT_UMA_MODEL == "uma-s-1p2p1"
    assert defaults.DEFAULT_ORB_MODEL == "orbmol-v2"
    assert defaults.MLMM_CALC_KW["uma_model"] == "uma-s-1p2p1"
    assert defaults.MLMM_CALC_KW["orb_model"] == "orbmol-v2"
    assert defaults.MLMM_CALC_KW["mace_model"] == "MACE-OMOL-0"
    assert defaults.MLMM_CALC_KW["aimnet2_model"] == "aimnet2"
    for factory in (
        mlmm_calc._create_ml_backend,
        mlmm_calc.MLMMCore.__init__,
        mlmm_calc.mlmm.__init__,
    ):
        params = inspect.signature(factory).parameters
        assert params["uma_model"].default == "uma-s-1p2p1"
        assert params["orb_model"].default == "orbmol-v2"
    for backend_cls, key, model in (
        (mlmm_calc._UMABackend, "uma_model", "uma-s-1p2p1"),
        (mlmm_calc._OrbBackend, "orb_model", "orbmol-v2"),
    ):
        assert inspect.signature(backend_cls.__init__).parameters[key].default == model


class ConservativeForcefieldRegressor:
    pass


class DirectForcefieldRegressor:
    pass


def _install_fake_orb(monkeypatch, *, registry, attrs=None):
    """Register a fake ``orb_models`` whose loaders record their calls."""
    calls = []
    pretrained = ModuleType("orb_models.forcefield.pretrained")
    pretrained.ORB_PRETRAINED_MODELS = {
        key: _recording_loader(key, result, calls) for key, result in registry.items()
    }
    for name, result in (attrs or {}).items():
        setattr(pretrained, name, _recording_loader(name, result, calls))

    calculator = ModuleType("orb_models.forcefield.inference.calculator")

    class ORBCalculator:
        def __init__(self, model, adapter=None, *, device=None):
            self.model, self.adapter, self.device = model, adapter, device

    calculator.ORBCalculator = ORBCalculator
    inference = ModuleType("orb_models.forcefield.inference")
    inference.calculator = calculator
    forcefield = ModuleType("orb_models.forcefield")
    forcefield.pretrained = pretrained
    forcefield.inference = inference
    package = ModuleType("orb_models")
    package.forcefield = forcefield
    for name, module in (
        ("orb_models", package),
        ("orb_models.forcefield", forcefield),
        ("orb_models.forcefield.pretrained", pretrained),
        ("orb_models.forcefield.inference", inference),
        ("orb_models.forcefield.inference.calculator", calculator),
    ):
        monkeypatch.setitem(sys.modules, name, module)
    monkeypatch.setattr(mlmm_calc, "HAS_ORB", True)
    return calls


def _recording_loader(name, result, calls):
    def loader(*, device, precision):
        calls.append((name, device, precision))
        return result()

    return loader


def _orb(model: str):
    return mlmm_calc._OrbBackend(
        orb_model=model, orb_precision="fp64", ml_device=torch.device("cpu")
    )


@pytest.mark.parametrize(
    ("model", "key"),
    [
        ("orbmol-v2", "orbmol-v2"),
        ("OrbMol-v2", "orbmol-v2"),
        ("orb_v3_conservative_omol", "orb-v3-conservative-omol"),
        ("orb-v3-conservative-omol", "orb-v3-conservative-omol"),
    ],
)
def test_orb_names_resolve_through_the_pretrained_registry(monkeypatch, model, key) -> None:
    def with_adapter():
        return ConservativeForcefieldRegressor(), "adapter"

    calls = _install_fake_orb(
        monkeypatch,
        registry={"orbmol-v2": with_adapter, "orb-v3-conservative-omol": with_adapter},
    )
    backend = _orb(model)
    assert calls == [(key, "cpu", "float64")]
    assert isinstance(backend._model_obj, ConservativeForcefieldRegressor)
    assert backend._adapter == "adapter"
    assert backend._ase_calc.adapter == "adapter"


def test_orb_falls_back_to_the_loader_function(monkeypatch) -> None:
    calls = _install_fake_orb(
        monkeypatch,
        registry={},
        attrs={"orb_v3_conservative_omol": ConservativeForcefieldRegressor},
    )
    backend = _orb("orb-v3-conservative-omol")
    assert calls == [("orb_v3_conservative_omol", "cpu", "float64")]
    assert backend._adapter is None


def test_orbmol_v2_on_old_orb_models_names_the_requirement(monkeypatch) -> None:
    # orb-models 0.5.5 knows OrbMol-v1 only.
    _install_fake_orb(
        monkeypatch,
        registry={"orb-v3-conservative-omol": ConservativeForcefieldRegressor},
        attrs={"orb_v3_conservative_omol": ConservativeForcefieldRegressor},
    )
    with pytest.raises(ImportError, match=r"orb-models >= 0\.7\.0") as excinfo:
        _orb("orbmol-v2")
    assert "--backend-model orb_v3_conservative_omol" in str(excinfo.value)


def test_orb_rejects_a_direct_force_regressor(monkeypatch) -> None:
    _install_fake_orb(
        monkeypatch,
        registry={"orb-v3-direct-omol": lambda: (DirectForcefieldRegressor(), "adapter")},
    )
    with pytest.raises(ValueError, match="only conservative ORB models are supported") as excinfo:
        _orb("orb_v3_direct_omol")
    assert "orb_v3_direct_omol" in str(excinfo.value)


def test_orb_unknown_model_is_an_error(monkeypatch) -> None:
    _install_fake_orb(monkeypatch, registry={})
    with pytest.raises(ValueError, match="Unknown ORB model 'orb-nonexistent'"):
        _orb("orb-nonexistent")


class _InfoRecordingCalculator(Calculator):
    implemented_properties = ["energy", "forces"]

    def __init__(self):
        super().__init__()
        self.seen_info = []
        self.seen_positions = []

    def calculate(self, atoms=None, properties=("energy",), system_changes=all_changes):
        super().calculate(atoms, properties, system_changes)
        self.seen_info.append(dict(self.atoms.info))
        self.seen_positions.append(self.atoms.get_positions())
        self.results = {"energy": 0.0, "forces": np.zeros((len(self.atoms), 3))}


def _install_fake_mace(monkeypatch, *, mace_polar=None):
    calculators = ModuleType("mace.calculators")
    for name in ("mace_anicc", "mace_mp", "mace_off", "mace_omol"):
        setattr(calculators, name, lambda **_kwargs: pytest.fail("wrong MACE factory"))
    calculators.MACECalculator = object
    foundations = ModuleType("mace.calculators.foundations_models")
    if mace_polar is not None:
        foundations.mace_polar = mace_polar
    calculators.foundations_models = foundations
    package = ModuleType("mace")
    package.__path__ = []
    monkeypatch.setattr(mlmm_calc, "HAS_MACE", True)
    monkeypatch.setitem(sys.modules, "mace", package)
    monkeypatch.setitem(sys.modules, "mace.calculators", calculators)
    monkeypatch.setitem(sys.modules, "mace.calculators.foundations_models", foundations)


@pytest.mark.parametrize(
    ("model", "key"),
    [
        ("MACE-POLAR-1-S", "polar-1-s"),
        ("MACE-POLAR-1-M", "polar-1-m"),
        ("mace-polar-1-l", "polar-1-l"),
        ("polar-1-s", "polar-1-s"),
        ("POLAR-1-M", "polar-1-m"),
        ("polar-1-L", "polar-1-l"),
    ],
)
def test_mace_polar_names_route_to_mace_polar(monkeypatch, model, key) -> None:
    captured = {}
    ase_calc = _InfoRecordingCalculator()

    def mace_polar(**kwargs):
        captured.update(kwargs)
        return ase_calc

    _install_fake_mace(monkeypatch, mace_polar=mace_polar)
    backend = mlmm_calc._MACEBackend(
        mace_model=model,
        mace_dtype="float64",
        model_charge=-1,
        model_mult=2,
        ml_device=torch.device("cpu"),
    )
    assert captured == {"model": key, "device": "cpu", "default_dtype": "float64"}
    assert backend._ase_calc is ase_calc

    # Total charge and multiplicity reach the calculator as for MACE-OMOL-0.
    backend.eval(Atoms("H2", positions=[[0.0, 0.0, 0.0], [0.0, 0.0, 0.74]]))
    assert ase_calc.seen_info[-1]["charge"] == -1
    assert ase_calc.seen_info[-1]["spin"] == 2


def test_mace_sees_coordinates_centred_on_the_origin(monkeypatch) -> None:
    ase_calc = _InfoRecordingCalculator()
    _install_fake_mace(monkeypatch, mace_polar=lambda **_kwargs: ase_calc)
    backend = mlmm_calc._MACEBackend(
        mace_model="MACE-POLAR-1-M", mace_dtype="float64", ml_device=torch.device("cpu")
    )
    xyz = np.array([[40.0, 50.0, 60.0], [41.0, 50.0, 60.0], [40.0, 52.0, 61.0]])
    atoms = Atoms("OH2", positions=xyz)
    _, _, opaque = backend.eval(atoms)
    ase_calc.reset()  # same positions; drop the cached results
    backend.energy(atoms)
    assert len(ase_calc.seen_positions) == 2
    for seen in [*ase_calc.seen_positions, opaque.get_positions()]:
        assert np.allclose(seen.min(axis=0), -seen.max(axis=0))
        assert np.allclose(seen - seen[0], xyz - xyz[0])
    assert np.allclose(atoms.get_positions(), xyz)


def test_mace_polar_without_mace_polar_names_the_mace_torch_version(monkeypatch) -> None:
    _install_fake_mace(monkeypatch, mace_polar=None)
    with pytest.raises(ImportError, match=r"MACE-POLAR-1 needs mace-torch >= 0\.3\.16"):
        mlmm_calc._MACEBackend(
            mace_model="MACE-POLAR-1-M",
            mace_dtype="float64",
            ml_device=torch.device("cpu"),
        )


@pytest.mark.parametrize(
    ("backend", "model", "label"),
    [
        ("uma", "uma-s-1p2p1", "UMA-S-1.2.1 (OMol)"),
        ("uma", "uma-s-1p2p1-ft-d4rev", "UMA-S-1.2.1-FT-D4rev (OMol)"),
        ("orb", "orbmol-v2", "OrbMol-v2"),
        ("orb", "orb_v3_conservative_omol", "ORB-v3-conservative-OMol"),
    ],
)
def test_model_labels_match_pdb2reaction(backend, model, label) -> None:
    from mlmm.core.output import mlip_model_label

    assert mlip_model_label(backend, model, "omol") == label


@pytest.mark.parametrize(
    ("backend", "model", "expected"),
    [
        ("orb", "orbmol-v2", {"orb_v3", "omol25"}),
        ("mace", "MACE-POLAR-1-M", {"mace", "mace_design", "omol25"}),
    ],
)
def test_new_models_cite_their_model_and_dataset(backend, model, expected) -> None:
    from mlmm.io.summary import _method_citation_record_keys

    keys = set(_method_citation_record_keys({"mlip_backend": backend, "mlip_model": model}))
    assert expected <= keys


def test_mace_polar_without_graph_longrange_names_the_package(monkeypatch) -> None:
    def missing(**_kwargs):
        raise ModuleNotFoundError("No module named 'graph_longrange'", name="graph_longrange")

    _install_fake_mace(monkeypatch, mace_polar=missing)
    with pytest.raises(ImportError, match="MACE-POLAR-1 needs graph-longrange 0.4.0") as excinfo:
        mlmm_calc._MACEBackend(
            mace_model="MACE-POLAR-1-M",
            mace_dtype="float64",
            ml_device=torch.device("cpu"),
        )
    assert "graph-longrange==0.4.0" in str(excinfo.value)


def test_mace_polar_other_missing_modules_propagate(monkeypatch) -> None:
    def missing(**_kwargs):
        raise ModuleNotFoundError("No module named 'something_else'", name="something_else")

    _install_fake_mace(monkeypatch, mace_polar=missing)
    with pytest.raises(ModuleNotFoundError) as excinfo:
        mlmm_calc._MACEBackend(
            mace_model="MACE-POLAR-1-M",
            mace_dtype="float64",
            ml_device=torch.device("cpu"),
        )
    assert excinfo.value.name == "something_else"


def test_uma_reference_energies_skip_checkpoints_without_them(monkeypatch) -> None:
    registry = {
        "uma-s-1p2p1": SimpleNamespace(atom_refs={"f": 1}, form_elem_refs={"f": 2}),
        "uma-s-1p2p1-ft-d4rev": SimpleNamespace(atom_refs=None, form_elem_refs=None),
    }
    fetched = []
    fake = SimpleNamespace(
        _MODEL_CKPTS=SimpleNamespace(checkpoints=registry),
        get_reference_energies=lambda model, reference_type: fetched.append((model, reference_type)) or {"ok": 1},
    )
    monkeypatch.setattr(mlmm_calc, "pretrained_mlip", fake, raising=False)
    assert mlmm_calc._uma_reference_energies("uma-s-1p2p1-ft-d4rev", "atom_refs") is None
    assert mlmm_calc._uma_reference_energies("uma-s-1p2p1-ft-d4rev", "form_elem_refs") is None
    assert mlmm_calc._uma_reference_energies("uma-s-1p2p1", "form_elem_refs") == {"ok": 1}
    assert fetched == [("uma-s-1p2p1", "form_elem_refs")]
