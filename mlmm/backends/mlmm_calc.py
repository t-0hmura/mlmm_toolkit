"""
ONIOM-like ML/MM calculator coupling MLIP backends (UMA, ORB, MACE, AIMNet2)
with hessian_ff (MM).

Example:
    calc = mlmm(input_pdb="input.pdb", real_parm7="real.parm7", model_pdb="model.pdb", charge=0)

For backend configuration, see: docs/backends.md
"""
# DOMAIN_PURE

from __future__ import annotations

import abc
import logging
import os
from pathlib import Path
import warnings
import shutil
import tempfile
import time
from contextlib import contextmanager, nullcontext
from dataclasses import dataclass
from typing import Any, Dict, List, Optional, Sequence, Tuple
from concurrent.futures import ThreadPoolExecutor

logger = logging.getLogger(__name__)

import click
from mlmm.core.output import emit
import numpy as np
import torch
import torch.nn as nn

from ase import Atoms
from ase.io import read as _ase_read
from ase.calculators.calculator import Calculator, all_changes
from ase.constraints import FixAtoms

from pysisyphus._array import active_square
from mlmm.backends.methods import (
    normalize_hessian_calc_mode,
    normalize_link_atom_method,
    normalize_mm_hessian_mode,
)


def _gather_atom_hessian_square(H, atom_indices):
    """Gather an atom-indexed 4-D Hessian square with bounded workspace."""
    n_atoms = int(H.shape[0])
    atom_indices = torch.as_tensor(
        atom_indices, dtype=torch.long, device=H.device
    )
    dof_indices = (
        atom_indices[:, None] * 3
        + torch.arange(3, dtype=torch.long, device=H.device)[None, :]
    ).reshape(-1)
    square = active_square(
        H.reshape(3 * n_atoms, 3 * n_atoms), dof_indices
    )
    selected = int(atom_indices.numel())
    return square.reshape(selected, 3, selected, 3)


def read(filename, *args, **kwargs):
    """ASE read() that tolerates virtual-site atoms (element 'EP', e.g. OPC/TIP4P 4-point water
    'EPW', or lone pairs 'LP'). ASE's label_to_symbol cannot parse 'EP'/'EPW', so we rewrite the
    element column (PDB cols 77-78) of such atoms to 'X' (ASE dummy, Z=0) before reading. Virtual
    sites are MM point charges handled by OpenMM via the parm and never enter the MLIP region, so
    the dummy symbol is bookkeeping only. Non-PDB inputs pass straight through to ASE."""
    import io as _io
    try:
        is_pdb = isinstance(filename, str) and filename.lower().endswith((".pdb", ".pdb1", ".ent"))
    except Exception:
        is_pdb = False
    if not is_pdb:
        return _ase_read(filename, *args, **kwargs)
    out = []
    with open(filename) as fh:
        for ln in fh:
            if (ln.startswith("ATOM") or ln.startswith("HETATM")) and len(ln) >= 78:
                elem = ln[76:78].strip().upper()
                name = ln[12:16].strip().upper()
                if elem in ("EP", "LP") or name.startswith(("EPW", "LP")):
                    ln = ln[:76] + " X" + ln[78:]
            out.append(ln)
    kwargs.setdefault("format", "proteindatabank")
    return _ase_read(_io.StringIO("".join(out)), *args, **kwargs)

import parmed as pmd
from hessian_ff import ForceFieldTorch, load_coords, load_system
from hessian_ff.analytical_hessian import build_analytical_hessian

from mlmm.core.defaults import DEFAULT_UMA_MODEL  # noqa: E402
from mlmm.io.pdb_indexing import resolve_mlmm_atoms

# Optional OpenMM import
try:
    import openmm as mm
    from openmm import app, unit, Platform
    from openmm.unit import ScaledUnit, joule
    HAS_OPENMM = True
except ImportError:
    HAS_OPENMM = False

# Strict deterministic mode is opt-in via the `--deterministic` CLI flag or
# `MLMM_STRICT_DETERMINISTIC=1` and lives in `mlmm.backends._determinism`
# (the backend factory and CLI option callback drive it). `--precision fp64`
# reduces numerical drift but does not replace strict deterministic mode.


# Optional fairchem import (UMA backend)
try:
    from fairchem.core import pretrained_mlip
    from fairchem.core.datasets.atomic_data import AtomicData
    from fairchem.core.datasets import data_list_collater
    HAS_FAIRCHEM = True
except ImportError:
    HAS_FAIRCHEM = False

# Optional: parallel MLIP predictor, only needed when workers > 1.
try:
    from fairchem.core.units.mlip_unit.predict import ParallelMLIPPredictUnit
    from fairchem.core.units.mlip_unit.api.inference import guess_inference_settings
except Exception:
    ParallelMLIPPredictUnit = None
    guess_inference_settings = None

# fp64 base precision: switching OMol-trained UMA from default fp32 to
# fp64 can have non-trivial impact on TSopt + Hessian numerics. Available
# via InferenceSettings(base_precision_dtype="float64") in fairchem ≥ 2.0.
try:
    from fairchem.core.units.mlip_unit.api.inference import InferenceSettings as _UMAInferenceSettings
except Exception:
    _UMAInferenceSettings = None

# Importing orb_models registers the ORB backend with ASE/torch.
# Optional ORB backend
try:
    import orb_models  # noqa: F401
    HAS_ORB = True
except ImportError:
    HAS_ORB = False

# Importing mace registers the MACE backend with ASE/torch.
# Optional MACE backend
try:
    import mace  # noqa: F401
    HAS_MACE = True
except ImportError:
    HAS_MACE = False

# Importing aimnet registers the AIMNet2 backend with ASE/torch.
# Optional AIMNet2 backend
try:
    import aimnet  # noqa: F401
    HAS_AIMNET2 = True
except ImportError:
    HAS_AIMNET2 = False

# ---------- PySisyphus unit constants ----------
from pysisyphus.constants import BOHR2ANG, ANG2BOHR, AU2EV, AU2KCALPERMOL
from pysisyphus.elem_data import COVALENT_RADII as _COVALENT_RADII_BOHR


def _get_g_factor(qm1_elem: str, mm_elem: str, link_elem: str = "H") -> float:
    """Compute Morokuma/Dapprich g-factor for link atom placement.

    g = (CR(QM1) + CR(link)) / (CR(QM1) + CR(MM))

    The link atom is placed at: r_L = r_QM1 + g * (r_MM - r_QM1).
    Units cancel (covalent radii in Bohr).
    """
    cr_qm1 = _COVALENT_RADII_BOHR[qm1_elem.lower()]
    cr_mm = _COVALENT_RADII_BOHR[mm_elem.lower()]
    cr_link = _COVALENT_RADII_BOHR[link_elem.lower()]
    return (cr_qm1 + cr_link) / (cr_qm1 + cr_mm)
EV2AU = 1.0 / AU2EV  # eV → Hartree
KCALMOL2EV = AU2EV / AU2KCALPERMOL  # kcal/mol -> eV


def _prepare_model_for_autograd_hessian(model_obj: Any) -> Dict[str, Any]:
    """Temporarily make a torch model deterministic and double-backward safe."""
    state: Dict[str, Any] = {
        "was_training": bool(getattr(model_obj, "training", False)),
        "param_flags": [],
        "dropout_states": [],
    }
    if hasattr(model_obj, "parameters"):
        for param in model_obj.parameters():
            state["param_flags"].append((param, bool(param.requires_grad)))
            param.requires_grad_(False)
    if hasattr(model_obj, "train"):
        model_obj.train(True)
    dropout_types = tuple(
        cls
        for cls in (
            nn.Dropout,
            nn.Dropout1d,
            nn.Dropout2d,
            nn.Dropout3d,
            nn.AlphaDropout,
            nn.FeatureAlphaDropout,
        )
        if cls is not None
    )
    if hasattr(model_obj, "modules"):
        for module in model_obj.modules():
            if not isinstance(module, dropout_types):
                continue
            old_p = getattr(module, "p", None)
            state["dropout_states"].append(
                (module, bool(getattr(module, "training", False)), old_p)
            )
            if old_p is not None:
                module.p = 0.0
            module.train(False)
    return state


def _restore_model_after_autograd_hessian(
    model_obj: Any, state: Dict[str, Any]
) -> None:
    """Restore state saved by :func:`_prepare_model_for_autograd_hessian`."""
    for module, was_training, old_p in state.get("dropout_states", []):
        if old_p is not None:
            module.p = old_p
        module.train(was_training)
    if hasattr(model_obj, "train"):
        model_obj.train(state.get("was_training", False))
    for param, requires_grad in state.get("param_flags", []):
        param.requires_grad_(requires_grad)


def _autograd_hessian_with_mutation_guard(energy_fn, flat0: torch.Tensor) -> torch.Tensor:
    """Evaluate a Hessian while permitting ORB's saved-tensor mutations."""
    graph = getattr(torch.autograd, "graph", None)
    mutation_guard = getattr(graph, "allow_mutation_on_saved_tensors", nullcontext)
    with mutation_guard():
        return torch.autograd.functional.hessian(
            energy_fn,
            flat0,
            vectorize=False,
            create_graph=False,
        )




class _MLBackend(abc.ABC):
    """Internal abstraction for the ML part of the ONIOM ML/MM coupling.

    Each backend must provide energy/force evaluation and Hessian computation.
    All quantities are in eV and Angstrom.
    """

    @abc.abstractmethod
    def eval(
        self, atoms: Atoms, need_grad: bool = True
    ) -> Tuple[float, np.ndarray, Any]:
        """Evaluate energy and forces.

        Returns
        -------
        E : float
            Energy in eV.
        F : ndarray (N, 3)
            Forces in eV/Å.
        opaque : Any
            Backend-specific data needed for analytical Hessian (e.g., batch).
        """

    def energy(self, atoms: Atoms) -> float:
        """Evaluate only the high-level energy.

        Backends with a native energy-only route override this method.  The
        compatibility fallback keeps third-party backends working, while
        MLMMCore can give DFT and other energy-only evaluators a force-free
        execution path.
        """
        energy, _, _ = self.eval(atoms, need_grad=False)
        return float(energy)

    @abc.abstractmethod
    def hessian_analytical(self, opaque: Any, n_atoms: int, *, dtype: torch.dtype) -> torch.Tensor:
        """Compute analytical Hessian from the opaque batch returned by eval().

        Returns Hessian as a (N, 3, N, 3) torch Tensor in eV/Å².
        """

    def hessian_fd(
        self,
        atoms: Atoms,
        freeze_model: Sequence[int],
        *,
        eps_ang: float = 1.0e-3,
        dtype: torch.dtype = torch.float32,
        device: torch.device = torch.device("cpu"),
    ) -> torch.Tensor:
        """Compute Hessian via finite differences (central difference).

        Generic implementation that works for all backends.
        """
        n_atoms = len(atoms)
        dof = n_atoms * 3

        frozen_set = set(int(i) for i in freeze_model)
        active_atoms = [i for i in range(n_atoms) if i not in frozen_set]
        active_dof_idx = [3 * i + j for i in active_atoms for j in range(3)]

        # On the native in-process torch-model path, assemble columns directly
        # from device force tensors. ASE/MM calculators and the parallel
        # predictor use the NumPy path. Both paths cast displacement forces to
        # the requested Hessian dtype and device before central differencing.
        native_tensor = callable(getattr(self, "forces_tensor", None)) and not getattr(
            self, "parallel_predict", False
        )

        H = torch.zeros((dof, dof), device=device, dtype=dtype)
        coord0 = atoms.get_positions().copy()
        for k in active_dof_idx:
            a = k // 3
            c = k % 3

            atoms.positions = coord0.copy()
            atoms.positions[a, c] = coord0[a, c] + eps_ang
            if native_tensor:
                Fp_t = self.forces_tensor(atoms).reshape(-1).to(device, dtype=dtype)
            else:
                _, Fp, _ = self.eval(atoms, need_grad=False)
                Fp_t = torch.from_numpy(Fp.reshape(-1)).to(device, dtype=dtype)

            atoms.positions = coord0.copy()
            atoms.positions[a, c] = coord0[a, c] - eps_ang
            if native_tensor:
                Fm_t = self.forces_tensor(atoms).reshape(-1).to(device, dtype=dtype)
            else:
                _, Fm, _ = self.eval(atoms, need_grad=False)
                Fm_t = torch.from_numpy(Fm.reshape(-1)).to(device, dtype=dtype)

            col = -(Fp_t - Fm_t) / (2.0 * eps_ang)
            H[:, k] = col

        atoms.positions = coord0
        return H.view(n_atoms, 3, n_atoms, 3)

    @property
    @abc.abstractmethod
    def supports_analytical_hessian(self) -> bool:
        """Whether this backend supports analytical Hessian."""

    @property
    @abc.abstractmethod
    def device(self) -> torch.device:
        """The torch device this backend uses."""


@contextmanager
def _uma_analytical_head_scope(model):
    """Enable UMA's EFS derivative graph without training its backbone."""
    inner = getattr(model, "module", model)
    try:
        head = inner.output_heads["energyandforcehead"].head
        backbone = inner.backbone
    except (AttributeError, KeyError, TypeError) as exc:
        raise RuntimeError(
            "UMA analytical Hessian requires the energyandforcehead EFS head; "
            "use FiniteDifference with this unsupported model layout."
        ) from exc
    modules = list(model.modules())
    if (
        not isinstance(head, nn.Module)
        or not isinstance(backbone, nn.Module)
        or not any(head is module for module in modules)
        or any(head is module for module in backbone.modules())
        or head is model
        or head is inner
    ):
        raise RuntimeError("UMA analytical Hessian requires a separate EFS head.")

    module_flags = [(module, module.training) for module in modules]
    parameter_flags = [(parameter, parameter.requires_grad) for parameter in model.parameters()]
    dropout_types = (
        nn.Dropout, nn.Dropout1d, nn.Dropout2d, nn.Dropout3d,
        nn.AlphaDropout, nn.FeatureAlphaDropout,
    )
    dropout_flags = [
        (module, module.p) for module in head.modules()
        if isinstance(module, dropout_types)
    ]
    try:
        for parameter, _ in parameter_flags:
            parameter.requires_grad_(False)
        # Backbone training enables functional composition dropout in UMA.
        # Only the EFS head needs training=True to retain the force graph.
        model.eval()
        head.train(True)
        for module, _ in dropout_flags:
            module.p = 0.0
            module.training = False
        yield
    finally:
        for module, probability in dropout_flags:
            module.p = probability
        # Recursive train()/eval() would overwrite saved mixed child states.
        for module, training in module_flags:
            module.training = training
        for parameter, requires_grad in parameter_flags:
            parameter.requires_grad_(requires_grad)


class _UMABackend(_MLBackend):
    """UMA (FAIR-Chem) ML backend."""

    def __init__(
        self,
        *,
        uma_model: str = DEFAULT_UMA_MODEL,
        uma_task_name: str = "omol",
        model_charge: int = 0,
        model_mult: int = 1,
        ml_device: torch.device,
        precision: str = "fp32",
        workers: int = 1,
        workers_per_node: int = 1,
        analytical_hessian: bool = False,
    ):
        if not HAS_FAIRCHEM:
            raise ImportError(
                "fairchem-core is required for the UMA backend. "
                "Install with `pip install fairchem-core` "
                "and ensure Hugging Face authentication is configured."
            )
        self._device = ml_device
        device_str = "cuda" if ml_device.type == "cuda" else "cpu"
        self._AtomicData = AtomicData
        self._data_list_collater = data_list_collater

        # fp32 is the established baseline; fp64 enables full-precision
        # base inference and can change TSopt + Hessian numerics non-trivially.
        self.precision = str(precision or "fp32").lower()
        if self.precision not in ("fp32", "fp64"):
            raise ValueError(f"UMA precision must be 'fp32' or 'fp64', got {precision!r}")
        if self.precision == "fp64" and _UMAInferenceSettings is None:
            raise ImportError(
                "UMA precision='fp64' requires fairchem-core's InferenceSettings; "
                "upgrade fairchem-core (≥ 2.0) or pass precision='fp32'."
            )
        self.workers = max(int(workers or 1), 1)
        self.workers_per_node = max(int(workers_per_node or 1), 1)
        self.parallel_predict = self.workers > 1

        _uma_inference_settings = None
        if _UMAInferenceSettings is not None and (
            self.precision == "fp64" or analytical_hessian
        ):
            # FAIR-Chem 2.22's named default enables torch.compile. The
            # compiled backward does not support the double backward used by
            # ML/MM analytical Hessians, so that route requests the public
            # non-compiled settings object explicitly.
            _uma_inference_settings = _UMAInferenceSettings(
                compile=False,
                base_precision_dtype=(
                    "float64" if self.precision == "fp64" else "float32"
                )
            )
        if self.parallel_predict:
            # ParallelMLIPPredictUnit spreads inference over `workers` processes but
            # does NOT expose `.model`; analytical Hessians are therefore unavailable
            # in this mode. Automatic selection uses finite differences; an explicit
            # Analytical request is rejected by the preflight before model loading.
            if ParallelMLIPPredictUnit is None or guess_inference_settings is None:
                raise ImportError(
                    "workers>1 requested, but ParallelMLIPPredictUnit/guess_inference_settings "
                    "could not be imported from fairchem. Install `fairchem-core[extras]`."
                )
            ckpt_path = pretrained_mlip.pretrained_checkpoint_path_from_name(uma_model)
            inference_settings = _uma_inference_settings or guess_inference_settings("default")
            atom_refs = pretrained_mlip.get_reference_energies(uma_model, reference_type="atom_refs")
            form_elem_refs = pretrained_mlip.get_reference_energies(uma_model, reference_type="form_elem_refs")
            self.predictor = ParallelMLIPPredictUnit(
                inference_model_path=str(ckpt_path),
                device=device_str,
                inference_settings=inference_settings,
                atom_refs=atom_refs,
                form_elem_refs=form_elem_refs,
                num_workers=self.workers,
                num_workers_per_node=self.workers_per_node,
            )
        else:
            # Serial in-process predictor used when workers=1.
            _uma_kwargs = {"device": device_str}
            if _uma_inference_settings is not None:
                _uma_kwargs["inference_settings"] = _uma_inference_settings
            device_context = (
                torch.cuda.device(ml_device)
                if ml_device.type == "cuda"
                else nullcontext()
            )
            with device_context:
                self.predictor = pretrained_mlip.get_predict_unit(
                    uma_model, **_uma_kwargs
                )

        self.uma_task_name = uma_task_name
        self.model_charge = model_charge
        self.model_mult = model_mult
        # ParallelMLIPPredictUnit has no `.model`; guard every torch-model access.
        self._has_torch_model = hasattr(self.predictor, "model") and isinstance(
            getattr(self.predictor, "model", None), nn.Module
        )
        if self._has_torch_model:
            self.predictor.model.eval()
            for m in self.predictor.model.modules():
                if isinstance(m, nn.Dropout):
                    m.p = 0.0
            backbone = getattr(self.predictor.model, "module", self.predictor.model).backbone
            self._uma_max_neigh = getattr(backbone, "max_neighbors", None)
            self._uma_radius = getattr(backbone, "cutoff", None)
        else:
            # Graph-construction cutoffs live inside the worker processes; let
            # AtomicData.from_ase fall back to the checkpoint defaults (None).
            self._uma_max_neigh = None
            self._uma_radius = None

    @property
    def supports_analytical_hessian(self) -> bool:
        # ParallelMLIPPredictUnit (workers>1) exposes no `.model` for autograd, so
        # analytical Hessians require the in-process predictor (workers=1).
        return self._has_torch_model

    @property
    def device(self) -> torch.device:
        return self._device

    def _predict(self, atoms: Atoms, need_grad: bool) -> Tuple[Any, Any]:
        """Run the in-process/parallel predictor and return ``(res, batch)``.

        Shared by :meth:`eval` (reads energy and forces to NumPy) and
        :meth:`forces_tensor` (keeps the force tensor on its device).
        """
        # fairchem/OMol interprets atoms.info["spin"] as multiplicity (2S+1);
        # zero is its unspecified token. `model_mult` is already a multiplicity.
        atoms.info.update({"charge": self.model_charge, "spin": self.model_mult})
        # When the inference path uses fp64, hand AtomicData the matching
        # target dtype so it does not down-cast positions to fp32 only to be
        # re-upcasted (and emit the fairchem `Upcasting atomic coordinates`
        # WARNING on every call). fp32 path keeps fairchem's default.
        target_dtype = torch.float64 if self.precision == "fp64" else torch.float32
        # AtomicData.from_ase reads charge/spin from atoms.info only when the
        # corresponding keys are listed in r_data_keys.
        data = self._AtomicData.from_ase(
            atoms,
            max_neigh=self._uma_max_neigh,
            radius=self._uma_radius,
            r_edges=False,
            target_dtype=target_dtype,
            r_data_keys=["spin", "charge"],
        )
        data.dataset = self.uma_task_name
        if self.parallel_predict:
            # ParallelMLIPPredictUnit moves tensors to its worker devices itself and
            # exposes no autograd graph; forces come straight from the predictor
            # output (analytical Hessian is unavailable in this mode, so need_grad
            # is irrelevant here).
            batch = self._data_list_collater([data], otf_graph=True)
            res = self.predictor.predict(batch)
            return res, batch
        # FAIR-Chem owns the device transfer. Since 2.22 its first prediction
        # prepares the still-CPU model from this batch before moving both to
        # the execution device.
        # First-call mismatch analysis: https://github.com/t-0hmura/pdb2reaction/pull/298
        batch = self._data_list_collater([data], otf_graph=True)
        pos = batch.pos.detach().clone()
        pos.requires_grad_(need_grad)
        batch.pos = pos
        if need_grad:
            res = self.predictor.predict(batch)
        else:
            with torch.no_grad():
                res = self.predictor.predict(batch)
        return res, batch

    def eval(self, atoms: Atoms, need_grad: bool = True) -> Tuple[float, np.ndarray, Any]:
        res, batch = self._predict(atoms, need_grad)
        E = float(res["energy"].squeeze().detach().item())
        F = res["forces"].detach().cpu().numpy()
        return E, F, batch

    def energy(self, atoms: Atoms) -> float:
        res, _ = self._predict(atoms, need_grad=False)
        return float(res["energy"].squeeze().detach().item())

    def forces_tensor(self, atoms: Atoms) -> torch.Tensor:
        """Return device-native forces (eV/Å) for the FD Hessian assembler.

        Returns the force tensor on the model's execution device, skipping the
        ``.detach().cpu().numpy()`` round-trip and the scalar energy ``.item()``
        sync that :meth:`eval` performs for each displacement. Only valid on the
        native in-process torch-model path (``_has_torch_model`` with
        ``parallel_predict`` false); the parallel predictor has no in-process
        device tensor to hand back, so callers must gate on those flags and fall
        back to :meth:`eval`. The returned values are the same tensor ``eval``
        converts to NumPy.
        """
        res, _ = self._predict(atoms, need_grad=False)
        return res["forces"].detach()

    def hessian_analytical(self, opaque: Any, n_atoms: int, *, dtype: torch.dtype) -> torch.Tensor:
        batch = opaque
        try:
            with _uma_analytical_head_scope(self.predictor.model):
                pos = batch.pos

                def energy_fn(flat_pos: torch.Tensor):
                    batch.pos = flat_pos.view(-1, 3)
                    return self.predictor.predict(batch)["energy"].squeeze()

                # Analytical autograd Hessians allocate O(N²) tensors. Convert a
                # CUDA allocation failure into an actionable user-facing error.
                try:
                    H_flat = torch.autograd.functional.hessian(energy_fn, pos.view(-1), vectorize=False)
                except torch.cuda.OutOfMemoryError as _oom:
                    if self._device.type == "cuda":
                        torch.cuda.empty_cache()
                    raise RuntimeError(
                        f"Analytical Hessian (torch.autograd) ran out of GPU memory "
                        f"for {n_atoms} atoms ({_oom}). The analytical Hessian needs "
                        f"substantially more VRAM than finite differences; rerun with "
                        f"`--hessian-calc-mode FiniteDifference` (the default), or use "
                        f"a GPU with more memory / a smaller ML region."
                    ) from _oom
                H = H_flat.view(n_atoms, 3, n_atoms, 3).to(
                    device=self._device, dtype=dtype
                ).detach()
                # release the autograd-graph-bearing source tensor immediately so
                # the cast scratch peak (= 2× hessian) collapses to 1×.
                del H_flat
        finally:
            if self._device.type == "cuda":
                torch.cuda.empty_cache()
        return H


class _ASEMLBackend(_MLBackend):
    """Base class for ASE-calculator-based ML backends (ORB, MACE, AIMNet2).

    Subclasses must set ``self._ase_calc`` (an ASE Calculator) and
    ``self._device``.
    """

    _ase_calc: Calculator
    _device: torch.device
    _model_charge: int = 0
    _model_mult: int = 1

    @property
    def supports_analytical_hessian(self) -> bool:
        return False

    @property
    def device(self) -> torch.device:
        return self._device

    def eval(self, atoms: Atoms, need_grad: bool = True) -> Tuple[float, np.ndarray, Any]:
        atoms_copy = atoms.copy()
        atoms_copy.calc = self._ase_calc
        # Propagate charge/spin to ASE Atoms info for backends that use them.
        # AIMNet2 reads 'charge' + 'mult'; ORB/MACE (OMol) read 'charge' + 'spin'.
        # OMol-trained ORB/MACE models interpret "spin" as multiplicity (2S+1);
        # zero is an unspecified token. AIMNet2 uses its separate "mult" key.
        atoms_copy.info["charge"] = self._model_charge
        atoms_copy.info["mult"] = self._model_mult
        atoms_copy.info["spin"] = int(self._model_mult)
        E = float(atoms_copy.get_potential_energy())
        F = np.array(atoms_copy.get_forces(), dtype=np.float64)
        # Keep the prepared Atoms object as the opaque analytical-Hessian input.
        return E, F, atoms_copy

    def energy(self, atoms: Atoms) -> float:
        atoms_copy = atoms.copy()
        atoms_copy.calc = self._ase_calc
        atoms_copy.info["charge"] = self._model_charge
        atoms_copy.info["mult"] = self._model_mult
        atoms_copy.info["spin"] = int(self._model_mult)
        return float(atoms_copy.get_potential_energy())

    def hessian_analytical(self, opaque: Any, n_atoms: int, *, dtype: torch.dtype) -> torch.Tensor:
        raise NotImplementedError(
            f"Analytical Hessian is not supported by {self.__class__.__name__}. "
            "Use hessian_calc_mode='FiniteDifference'."
        )


def _normalize_orb_precision(precision: str) -> str:
    """Normalize unified/legacy ORB precision tokens case-insensitively."""
    token = str(precision).strip().lower()
    return _OrbBackend._PRECISION_ALIASES.get(token, token)


class _OrbBackend(_ASEMLBackend):
    """ORB (Orbital Materials) ML backend.

    `orb_precision` is mapped to ORB's pretrained `precision=` kwarg
    (e.g. ``"float32-high"`` or ``"float64"``). The historical alias
    ``"float32"`` is silently rewritten to ``"float32-high"`` for
    backward compatibility with older YAML defaults.
    """

    _PRECISION_ALIASES = {
        "float32": "float32-high",
        "fp32": "float32-high",
        "fp64": "float64",
    }

    def __init__(
        self,
        *,
        orb_model: str = "orb_v3_conservative_omol",
        orb_precision: str = "float64",
        model_charge: int = 0,
        model_mult: int = 1,
        ml_device: torch.device,
        **_kwargs,
    ):
        if not HAS_ORB:
            raise ImportError(
                "orb-models is required for the ORB backend. "
                "Install with `pip install orb-models`."
            )
        from orb_models.forcefield import pretrained

        device_str = str(ml_device)
        precision = _normalize_orb_precision(orb_precision)
        loaded = getattr(pretrained, orb_model)(
            device=device_str, precision=precision
        )
        if isinstance(loaded, tuple) and len(loaded) >= 2:
            self._model_obj, self._adapter = loaded[0], loaded[1]
        else:
            self._model_obj, self._adapter = loaded, None

        def _construct(calculator_cls):
            attempts = []
            if self._adapter is not None:
                attempts.extend(
                    [
                        ((self._model_obj, self._adapter), {"device": device_str}),
                        ((self._model_obj, self._adapter), {}),
                    ]
                )
            attempts.extend(
                [
                    ((self._model_obj,), {"device": device_str}),
                    ((self._model_obj,), {}),
                ]
            )
            for args, kwargs in attempts:
                try:
                    return calculator_cls(*args, **kwargs)
                except TypeError:
                    continue
            return None

        self._ase_calc = None
        try:
            from orb_models.forcefield.inference.calculator import ORBCalculator

            self._ase_calc = _construct(ORBCalculator)
        except ImportError:
            pass
        if self._ase_calc is None:
            from orb_models.forcefield.calculator import ORBCalculator

            self._ase_calc = _construct(ORBCalculator)
        if self._ase_calc is None:
            raise RuntimeError("Failed to build ORBCalculator.")
        self._device = ml_device
        self._model_charge = model_charge
        self._model_mult = model_mult
        self._orb_precision = precision

    @property
    def supports_analytical_hessian(self) -> bool:
        return hasattr(self._model_obj, "predict")

    def hessian_analytical(
        self, opaque: Any, n_atoms: int, *, dtype: torch.dtype
    ) -> torch.Tensor:
        """Compute ORB's analytical Hessian with double-backward safeguards."""
        atoms = opaque.copy()
        atoms.info.update(
            {"charge": self._model_charge, "spin": int(self._model_mult)}
        )
        device_str = str(self._device)
        if self._adapter is not None and hasattr(self._adapter, "from_ase_atoms"):
            base_graph = self._adapter.from_ase_atoms(
                atoms=atoms, device=device_str
            )
        else:
            try:
                from orb_models.forcefield import atomic_system
            except ImportError as exc:
                raise RuntimeError(
                    "ORB analytical Hessian requires orb_models.forcefield.atomic_system."
                ) from exc
            base_graph = atomic_system.ase_atoms_to_atom_graphs(
                atoms,
                getattr(self._model_obj, "system_config", None),
                device=device_str,
            )
        if (
            not hasattr(base_graph, "node_features")
            or "positions" not in base_graph.node_features
        ):
            raise RuntimeError(
                "Unexpected ORB graph: node_features['positions'] is missing."
            )

        flat0 = (
            base_graph.node_features["positions"]
            .detach()
            .clone()
            .reshape(-1)
            .to(device_str)
        )

        def energy_fn(flat_pos: torch.Tensor) -> torch.Tensor:
            base_graph.node_features["positions"] = flat_pos.view(n_atoms, 3)
            result = self._model_obj.predict(base_graph)
            if isinstance(result, dict):
                for key in ("energy", "free_energy", "total_energy", "E"):
                    if key in result:
                        return result[key].reshape(-1)[0]
                raise RuntimeError(
                    f"ORB predict() output has no energy key: {sorted(result)}"
                )
            return result.reshape(-1)[0]

        donated_before = None
        try:
            donated_before = torch._functorch.config.donated_buffer
            torch._functorch.config.donated_buffer = False
        except (AttributeError, RuntimeError):
            donated_before = None
        state = _prepare_model_for_autograd_hessian(self._model_obj)
        try:
            hessian = _autograd_hessian_with_mutation_guard(energy_fn, flat0)
        except (torch.cuda.OutOfMemoryError, RuntimeError) as exc:
            if "out of memory" in str(exc).lower():
                raise RuntimeError(
                    "ORB analytical Hessian ran out of GPU memory; use "
                    "hessian_calc_mode='FiniteDifference'."
                ) from exc
            raise RuntimeError(f"ORB analytical Hessian failed: {exc}") from exc
        finally:
            _restore_model_after_autograd_hessian(self._model_obj, state)
            if donated_before is not None:
                try:
                    torch._functorch.config.donated_buffer = donated_before
                except AttributeError:
                    pass
            if self._device.type == "cuda":
                torch.cuda.empty_cache()
        return hessian.detach().reshape(n_atoms, 3, n_atoms, 3).to(dtype)


class _MACEBackend(_ASEMLBackend):
    """MACE ML backend."""

    def __init__(
        self,
        *,
        mace_model: str = "MACE-OMOL-0",
        mace_dtype: str = "float64",
        model_charge: int = 0,
        model_mult: int = 1,
        ml_device: torch.device,
    ):
        if not HAS_MACE:
            raise ImportError(
                "mace-torch is required for the MACE backend. "
                "Install with `pip install mace-torch`."
            )
        from mace.calculators import (
            MACECalculator,
            mace_anicc,
            mace_mp,
            mace_off,
            mace_omol,
        )

        device_str = str(ml_device)
        model_lower = mace_model.lower()
        off23_aliases = {
            "mace-off23_small": "small",
            "mace-off23_medium": "medium",
            "mace-off23_large": "large",
        }

        # Resolve model name to the appropriate factory
        if model_lower in off23_aliases:
            self._ase_calc = mace_off(
                model=off23_aliases[model_lower],
                device=device_str,
                default_dtype=mace_dtype,
            )
        elif model_lower.startswith("mp:") or model_lower.startswith("mace-mp"):
            model_name = mace_model.split(":", 1)[-1] if ":" in mace_model else mace_model
            self._ase_calc = mace_mp(
                model=model_name, device=device_str, default_dtype=mace_dtype
            )
        elif model_lower.startswith("off:") or model_lower.startswith("mace-off"):
            model_name = mace_model.split(":", 1)[-1] if ":" in mace_model else mace_model
            self._ase_calc = mace_off(
                model=model_name, device=device_str, default_dtype=mace_dtype
            )
        elif model_lower.startswith("anicc") or model_lower.startswith("mace-anicc"):
            if mace_dtype == "float64":
                self._ase_calc = mace_anicc(device=device_str)
            else:
                raw_model = mace_anicc(
                    device=device_str,
                    return_raw_model=True,
                )
                self._ase_calc = MACECalculator(
                    models=raw_model,
                    device=device_str,
                    default_dtype=mace_dtype,
                )
        elif model_lower.startswith("omol") or model_lower.startswith("mace-omol"):
            # MACE-OMOL loads via the dedicated mace_omol factory. mace_off treats any non-preset,
            # non-URL string as a LOCAL file path, so the default "MACE-OMOL-0"
            # raises FileNotFoundError. mace_omol maps "extra_large"/None to the
            # published OMOL-0 checkpoint.
            self._ase_calc = mace_omol(
                model="extra_large", device=device_str, default_dtype=mace_dtype
            )
        else:
            # Treat as a local model file or direct mace_off model
            self._ase_calc = mace_off(
                model=mace_model, device=device_str, default_dtype=mace_dtype
            )

        self._device = ml_device
        self._model_charge = model_charge
        self._model_mult = model_mult

    @property
    def supports_analytical_hessian(self) -> bool:
        return hasattr(self._ase_calc, "get_hessian")

    def hessian_analytical(
        self, opaque: Any, n_atoms: int, *, dtype: torch.dtype
    ) -> torch.Tensor:
        """Return MACE's native analytical Hessian in eV/Å²."""
        atoms = opaque.copy()
        atoms.info.update(
            {"charge": self._model_charge, "spin": int(self._model_mult)}
        )
        calc = self._ase_calc
        if not hasattr(calc, "get_hessian"):
            raise RuntimeError(
                "Installed MACE calculator has no analytical Hessian API; "
                "upgrade mace-torch or use FiniteDifference."
            )
        try:
            internal = (
                hasattr(calc, "_atoms_to_batch")
                and hasattr(calc, "_clone_batch")
                and getattr(calc, "models", None) is not None
                and len(calc.models) == 1
            )
            if internal:
                batch = calc._atoms_to_batch(atoms)
                result = calc.models[0](
                    calc._clone_batch(batch).to_dict(),
                    compute_hessian=True,
                    compute_stress=False,
                    training=getattr(calc, "use_compile", False),
                )
                hessian = result["hessian"].detach()
                del result
            else:
                hessian = calc.get_hessian(atoms=atoms)
        except (torch.cuda.OutOfMemoryError, RuntimeError) as exc:
            if "out of memory" in str(exc).lower():
                raise RuntimeError(
                    "MACE analytical Hessian ran out of GPU memory; use "
                    "hessian_calc_mode='FiniteDifference'."
                ) from exc
            raise RuntimeError(f"MACE analytical Hessian failed: {exc}") from exc
        if not isinstance(hessian, torch.Tensor):
            hessian = torch.as_tensor(hessian, device=self._device)
        if hessian.ndim == 5 and hessian.shape[0] > 0:
            hessian = hessian[0]
        return hessian.detach().reshape(n_atoms, 3, n_atoms, 3).to(dtype)


class _AIMNet2Backend(_ASEMLBackend):
    """AIMNet2 ML backend."""

    def __init__(
        self,
        *,
        aimnet2_model: str = "aimnet2",
        model_charge: int = 0,
        model_mult: int = 1,
        ml_device: torch.device,
    ):
        if not HAS_AIMNET2:
            raise ImportError(
                "aimnet is required for the AIMNet2 backend. "
                "Install with `pip install aimnet`."
            )
        from aimnet.calculators import AIMNet2ASE, AIMNet2Calculator

        device_str = str(ml_device)
        self._aimnet_base_calc = AIMNet2Calculator(
            model=aimnet2_model,
            device=device_str,
        )
        self._ase_calc = AIMNet2ASE(
            base_calc=self._aimnet_base_calc,
            charge=model_charge,
            mult=model_mult,
        )
        self._device = ml_device
        self._model_charge = model_charge
        self._model_mult = model_mult

    @property
    def supports_analytical_hessian(self) -> bool:
        return callable(getattr(self._ase_calc, "get_hessian", None))

    def hessian_analytical(
        self, opaque: Any, n_atoms: int, *, dtype: torch.dtype
    ) -> torch.Tensor:
        """Return AIMNet2's native analytical Hessian in eV/Å²."""
        atoms = opaque
        try:
            hessian = self._ase_calc.get_hessian(atoms)
        except (torch.cuda.OutOfMemoryError, RuntimeError) as exc:
            if "out of memory" in str(exc).lower():
                raise RuntimeError(
                    "AIMNet2 analytical Hessian ran out of GPU memory; use "
                    "hessian_calc_mode='FiniteDifference'."
                ) from exc
            raise RuntimeError(f"AIMNet2 analytical Hessian failed: {exc}") from exc
        if hessian is None:
            raise RuntimeError(
                "AIMNet2 did not return an analytical Hessian; use "
                "hessian_calc_mode='FiniteDifference'."
            )
        if not isinstance(hessian, torch.Tensor):
            hessian = torch.as_tensor(hessian, device=self._device)
        return hessian.detach().reshape(n_atoms, 3, n_atoms, 3).to(dtype)


class _CustomBackend(_ASEMLBackend):
    """User-supplied ASE Calculator loaded from a Python file (``--calc-file``).

    Drives the ML region with any ASE-compatible engine (GFN-xTB, DFTB+, ORCA,
    …) via the standard ``_ASEMLBackend`` adapter. See ``mlmm/backends/custom.py``
    for the file/factory contract.
    """

    def __init__(
        self,
        *,
        calc_file: str,
        calc_factory: str = "get_calculator",
        model_charge: int = 0,
        model_mult: int = 1,
        ml_device: torch.device,
    ):
        from mlmm.backends.custom import load_ase_calculator

        device_str = "cuda" if ml_device.type == "cuda" else "cpu"
        self._ase_calc = load_ase_calculator(
            calc_file,
            calc_factory,
            charge=model_charge,
            spin=model_mult,
            device=device_str,
        )
        self._device = ml_device
        self._model_charge = model_charge
        self._model_mult = model_mult


_ANNOUNCED_MODEL_LOADS: set = set()


@contextmanager
def _announce_model_load(backend: str, model: str, task_name: str = ""):
    """Bracket the first load of each model so a download cannot look like a hang."""
    from mlmm.core.output import emit, mlip_model_label

    model = str(model or "").strip()
    task_name = str(task_name or "").strip()
    if backend == "uma" and not task_name:
        task_name = "omol"
    key = (backend, model, task_name)
    if key in _ANNOUNCED_MODEL_LOADS:
        yield
        return
    _ANNOUNCED_MODEL_LOADS.add(key)
    backend_label = {
        "uma": "UMA", "orb": "ORB", "mace": "MACE", "aimnet2": "AIMNet2",
    }.get(backend, backend)
    model_label = mlip_model_label(backend, model, task_name)
    label = f"{backend_label}{f' / {model_label}' if model else ''}"
    emit(f"\n[backend] Preparing MLIP model ({label})...", narrative=True)
    try:
        yield
    except BaseException:
        _ANNOUNCED_MODEL_LOADS.discard(key)
        raise
    emit("[backend] Done.", narrative=True)
    emit("", narrative=True)


def _create_ml_backend(
    backend: str,
    *,
    uma_model: str = DEFAULT_UMA_MODEL,
    uma_task_name: str = "omol",
    uma_precision: str = "fp32",
    workers: int = 1,
    workers_per_node: int = 1,
    orb_model: str = "orb_v3_conservative_omol",
    orb_precision: str = "float64",
    mace_model: str = "MACE-OMOL-0",
    mace_dtype: str = "float64",
    aimnet2_model: str = "aimnet2",
    calc_file: Optional[str] = None,
    calc_factory: str = "get_calculator",
    model_charge: int = 0,
    model_mult: int = 1,
    ml_device: torch.device,
    analytical_hessian: bool = False,
) -> _MLBackend:
    """Factory function to create the appropriate ML backend."""
    model_mult = int(model_mult)
    if model_mult < 1:
        raise ValueError(f"Spin multiplicity must be >= 1, got {model_mult}.")
    backend = backend.strip().lower()
    # Env-var / direct-API entry point for strict determinism (the CLI uses the
    # --deterministic flag callback). Idempotent; no-op unless requested.
    from mlmm.backends._determinism import (
        is_deterministic_active,
        is_deterministic_requested,
        setup_deterministic,
    )
    if is_deterministic_requested():
        setup_deterministic()
    if backend == "aimnet2" and is_deterministic_active():
        raise ValueError(
            "AIMNet2's custom CUDA force kernel is outside "
            "torch.use_deterministic_algorithms control."
        )
    if backend == "custom" and is_deterministic_active():
        raise ValueError(
            "--deterministic is not supported with --calc-file because an "
            "arbitrary ASE Calculator may use RNGs, external processes, or "
            "kernels outside mlmm-toolkit's control."
        )
    if backend == "uma":
        with _announce_model_load(backend, uma_model, uma_task_name):
            return _UMABackend(
                uma_model=uma_model,
                uma_task_name=uma_task_name,
                model_charge=model_charge,
                model_mult=model_mult,
                ml_device=ml_device,
                precision=uma_precision,
                workers=workers,
                workers_per_node=workers_per_node,
                analytical_hessian=analytical_hessian,
            )
    elif backend == "orb":
        with _announce_model_load(backend, orb_model):
            return _OrbBackend(
                orb_model=orb_model,
                orb_precision=orb_precision,
                model_charge=model_charge,
                model_mult=model_mult,
                ml_device=ml_device,
            )
    elif backend == "mace":
        with _announce_model_load(backend, mace_model):
            return _MACEBackend(
                mace_model=mace_model,
                mace_dtype=mace_dtype,
                model_charge=model_charge,
                model_mult=model_mult,
                ml_device=ml_device,
            )
    elif backend == "aimnet2":
        with _announce_model_load(backend, aimnet2_model):
            return _AIMNet2Backend(
                aimnet2_model=aimnet2_model,
                model_charge=model_charge,
                model_mult=model_mult,
                ml_device=ml_device,
            )
    elif backend == "custom":
        if not calc_file:
            raise ValueError(
                "ML backend 'custom' requires --calc-file pointing to a Python "
                "file that exposes get_calculator(...) -> an ASE Calculator."
            )
        return _CustomBackend(
            calc_file=calc_file,
            calc_factory=calc_factory,
            model_charge=model_charge,
            model_mult=model_mult,
            ml_device=ml_device,
        )
    else:
        raise ValueError(
            f"Unknown ML backend '{backend}'. "
            "Choose from: uma, orb, mace, aimnet2, custom (--calc-file)."
        )




class _EmbedChargeCorrection:
    """Experimental xTB point-charge correction for ML/MM calculations.

    The retained implementation evaluates:

        dE = E_xTB(ML + MM_charges) - E_xTB(ML_only)
        dF = F_xTB(ML + MM_charges) - F_xTB(ML_only)

    The energy, force, and Hessian corrections are added to the subtractive
    ML/MM expression when ``embedcharge`` is enabled.
    """

    def __init__(
        self,
        *,
        xtb_cmd: str = "xtb",
        xtb_acc: float = 0.2,
        xtb_workdir: str = "tmp",
        xtb_keep_files: bool = False,
        xtb_ncores: int = 4,
        hessian_step: float = 1.0e-3,
    ):
        self.xtb_cmd = xtb_cmd
        self.xtb_acc = xtb_acc
        self.xtb_workdir = xtb_workdir
        self.xtb_keep_files = xtb_keep_files
        self.xtb_ncores = xtb_ncores
        self.hessian_step = hessian_step

    def compute_correction(
        self,
        symbols: List[str],
        coords_ml_ang: np.ndarray,
        mm_coords_ang: np.ndarray,
        mm_charges: np.ndarray,
        charge: int,
        multiplicity: int,
        *,
        need_forces: bool = False,
        need_hessian: bool = False,
    ) -> Tuple[float, Optional[np.ndarray], Optional[np.ndarray]]:
        """Compute point-charge embedding correction.

        Parameters
        ----------
        symbols : list of str
            Element symbols for ML atoms.
        coords_ml_ang : ndarray (N_ML, 3)
            Coordinates of ML atoms in Angstrom.
        mm_coords_ang : ndarray (N_MM, 3)
            Coordinates of MM point charges in Angstrom.
        mm_charges : ndarray (N_MM,)
            Charges of MM point charges in atomic units.
        charge : int
            Total charge of the ML region.
        multiplicity : int
            Spin multiplicity of the ML region.
        need_forces : bool
            Whether to compute force corrections.
        need_hessian : bool
            Whether to compute Hessian corrections.

        Returns
        -------
        dE : float
            Energy correction in eV.
        dF : ndarray (N_ML + N_MM, 3) or None
            Force corrections in ``[ML atoms, MM charge sites]`` order, eV/Å.
        dH : ndarray (3*(N_ML + N_MM), 3*(N_ML + N_MM)) or None
            Hessian correction in the same ordered basis, eV/Å².
        """
        from mlmm.backends.xtb_embedcharge_correction import delta_embedcharge_minus_noembed

        n_ml = len(symbols)
        mm_coords = np.asarray(mm_coords_ang, dtype=np.float64).reshape(-1, 3)
        mm_q = np.asarray(mm_charges, dtype=np.float64).reshape(-1)
        n_mm = mm_q.shape[0]

        if n_mm == 0:
            dF = np.zeros((n_ml, 3), dtype=np.float64) if need_forces else None
            dH = np.zeros((3 * n_ml, 3 * n_ml), dtype=np.float64) if need_hessian else None
            return 0.0, dF, dH

        dE_ev, dF_full_ev, dH_full_ev = delta_embedcharge_minus_noembed(
            symbols=symbols,
            coords_q_ang=np.asarray(coords_ml_ang, dtype=np.float64).reshape(-1, 3),
            mm_coords_ang=mm_coords,
            mm_charges=mm_q,
            charge=charge,
            multiplicity=multiplicity,
            need_forces=need_forces or need_hessian,
            need_hessian=need_hessian,
            xtb_cmd=self.xtb_cmd,
            xtb_acc=self.xtb_acc,
            xtb_workdir=self.xtb_workdir,
            xtb_keep_files=self.xtb_keep_files,
            ncores=self.xtb_ncores,
            hessian_step=self.hessian_step,
        )

        dF = None
        if dF_full_ev is not None:
            dF = np.asarray(dF_full_ev, dtype=np.float64).reshape(-1, 3)
            if dF.shape != (n_ml + n_mm, 3):
                raise ValueError(
                    "Embedding correction force shape does not match its "
                    f"[ML, MM] basis: {dF.shape} vs {(n_ml + n_mm, 3)}."
                )

        dH = None
        if dH_full_ev is not None:
            dH = np.asarray(dH_full_ev, dtype=np.float64)
            expected = 3 * (n_ml + n_mm)
            if dH.shape != (expected, expected):
                raise ValueError(
                    "Embedding correction Hessian shape does not match its "
                    f"[ML, MM] basis: {dH.shape} vs {(expected, expected)}."
                )

        return float(dE_ev), dF, dH



def _fixed_indices_from_constraints(atoms: Atoms) -> set[int]:
    fixed: set[int] = set()
    for c in atoms.constraints or []:
        if isinstance(c, FixAtoms):
            fixed.update(int(i) for i in c.get_indices())
    return fixed


def apply_cmap_policy(top, use_cmap: bool) -> None:
    """Drop CMAP terms from an ONIOM layer topology unless ``use_cmap``.

    The single owner of the CMAP rule. Both layers must be stripped or kept
    together: subtractive ONIOM only cancels the model region's MM description
    when ``E_real_low`` and ``E_model_low`` come from the same MM Hamiltonian.
    Removing CMAP from the sliced model alone left every cross-map term lying
    entirely inside the ML region uncancelled, i.e. an empirical backbone
    potential stacked on top of the high-level (ML or DFT) description.

    CMAP is retained by default when it is present in the parm7.  Explicit
    ``use_cmap=False`` is a modified-force-field opt-out and removes the terms
    from both ONIOM layers.
    """
    if not use_cmap:
        top.cmaps[:] = []


def write_model_parm7(
    real_top,
    selection,
    real_parm7,
    real_rst7,
    model_parm7,
    model_rst7,
    use_cmap: bool,
) -> None:
    """Write the sliced ML-region topology and coordinates.

    The single owner of the model-parm7 contract: every workflow that needs an
    ONIOM ``E_model_low`` goes through here, so the ML/MM and DFT paths can never
    describe the same model region differently.

    A selection covering the whole system copies the real files verbatim.
    Otherwise the slice drops CMAP unless ``use_cmap``, and its LJ tables are
    normalized -- ParmEd leaves ``LENNARD_JONES_*COEF`` at the parent's length
    when the selection uses fewer atom types, which our own MM backend rejects.
    """
    if len(selection) == len(real_top.atoms):
        shutil.copyfile(str(real_parm7), str(model_parm7))
        shutil.copyfile(str(real_rst7), str(model_rst7))
        return

    model = real_top[selection]
    model.box = None
    apply_cmap_policy(model, use_cmap)
    model.save(str(model_parm7), overwrite=True)
    _normalize_prmtop_lj_tables(str(model_parm7))
    model.save(str(model_rst7), overwrite=True)


def _normalize_prmtop_lj_tables(parm7_path: str) -> None:
    """Normalize LJ table lengths in parm7 files generated from sliced structures.

    ParmEd slicing can leave ``LENNARD_JONES_*COEF`` longer than the ``POINTERS``
    ``NTYPES`` expectation. Trim only the trailing unused tail when detected.
    """
    from parmed.amber import AmberFormat, AmberParm

    try:
        AmberParm(parm7_path)
        return
    except Exception as exc:
        msg = str(exc)
        if (
            "FLAG LENNARD_JONES_ACOEF" not in msg
            and "FLAG LENNARD_JONES_BCOEF" not in msg
        ):
            raise

    af = AmberFormat(parm7_path)
    pointers = list(af.parm_data.get("POINTERS", []))
    if len(pointers) < 2:
        raise ValueError(f"Invalid POINTERS section in parm7: {parm7_path}")
    ntypes = int(pointers[1])
    expected = ntypes * (ntypes + 1) // 2

    changed = False
    for key in ("LENNARD_JONES_ACOEF", "LENNARD_JONES_BCOEF"):
        values = list(af.parm_data.get(key, []))
        if len(values) == expected:
            continue
        if len(values) < expected:
            raise ValueError(
                f"{key} has {len(values)} entries but expected at least {expected} "
                f"from NTYPES={ntypes} in {parm7_path}."
            )
        af.parm_data[key] = values[:expected]
        changed = True

    if changed:
        af.write_parm(parm7_path)

    # Validate normalized topology immediately.
    AmberParm(parm7_path)


#                    hessian_ff (MM) -> ASE calculator

def _expand_partial_hessian(
    h_sub: np.ndarray,
    active_atoms: np.ndarray,
    n_atoms: int,
    *,
    dtype: np.dtype,
) -> np.ndarray:
    h_full = np.zeros((3 * n_atoms, 3 * n_atoms), dtype=dtype)
    for i_local, i_atom in enumerate(active_atoms):
        i0 = 3 * int(i_atom)
        for j_local, j_atom in enumerate(active_atoms):
            j0 = 3 * int(j_atom)
            h_full[i0:i0 + 3, j0:j0 + 3] = h_sub[
                3 * i_local:3 * i_local + 3,
                3 * j_local:3 * j_local + 3,
            ]
    return h_full


class hessianffCalculator(Calculator):
    """Calculator for MM. hessian_ff-backed."""

    implemented_properties = ["energy", "forces"]

    def __init__(
        self,
        parm7: str,
        rst7: Optional[str] = None,
        *,
        device: str = "auto",
        cuda_idx: int = 0,
        threads: int = 16,
        **kwargs,
    ):
        super().__init__(**kwargs)

        requested = str(device).lower()
        if requested not in {"auto", "cpu"}:
            raise ValueError(
                "MM backend 'hessian_ff' is CPU-only. "
                f"Got device={device!r}. Use mm_device='cpu' or 'auto'."
            )

        self.device = "cpu"
        self.cuda_idx = int(cuda_idx)
        self.threads = int(threads)
        if self.threads > 0 and torch.get_num_threads() != self.threads:
            torch.set_num_threads(self.threads)

        self.system = load_system(parm7, device="cpu").to(dtype=torch.float64)
        self.ff = ForceFieldTorch(self.system)
        self.natom = int(self.system.natom)
        self._coords_dtype = torch.float64
        self._coords_device = torch.device("cpu")
        self._coord_buf = torch.empty((self.natom, 3), dtype=self._coords_dtype, device=self._coords_device)

        if rst7 is not None:
            xyz = load_coords(rst7, natom=self.natom, device=self._coords_device, dtype=self._coords_dtype)
            self._coord_buf.copy_(xyz)

    def _positions_to_tensor(self, positions_ang: np.ndarray) -> torch.Tensor:
        arr = np.asarray(positions_ang, dtype=np.float64)
        if arr.shape != (self.natom, 3):
            raise ValueError(
                f"Coordinate shape mismatch for '{type(self).__name__}': "
                f"got {arr.shape}, expected ({self.natom}, 3)."
            )
        self._coord_buf.copy_(torch.as_tensor(arr, dtype=self._coords_dtype, device=self._coords_device))
        return self._coord_buf

    def _energy_forces_from_positions(self, positions_ang: np.ndarray) -> Tuple[float, np.ndarray]:
        xyz = self._positions_to_tensor(positions_ang)
        out, force = self.ff.energy_force(xyz, force_calc_mode="Analytical")
        energy_ev = float(out["E_total"].detach().cpu()) * KCALMOL2EV
        forces_ev = force.detach().cpu().numpy().astype(np.float64, copy=False) * KCALMOL2EV
        return energy_ev, forces_ev

    def _energy_from_positions(self, positions_ang: np.ndarray) -> float:
        xyz = self._positions_to_tensor(positions_ang)
        out = self.ff(xyz)
        return float(out["E_total"].detach().cpu()) * KCALMOL2EV

    def calculate(self, atoms: Atoms = None, properties=None, system_changes=all_changes):
        super().calculate(atoms, properties, system_changes)
        if atoms is None:
            raise ValueError("ASE Atoms is required for MM evaluation.")
        need_forces = properties is None or "forces" in properties
        if need_forces:
            energy_ev, forces_ev = self._energy_forces_from_positions(
                atoms.get_positions()
            )
            self.results = {"energy": energy_ev, "forces": forces_ev}
        else:
            self.results = {
                "energy": self._energy_from_positions(atoms.get_positions())
            }

    def analytical_hessian(
        self,
        atoms: Atoms,
        *,
        info_path: Optional[str] = None,
        dtype: np.dtype = np.float64,
        return_partial_hessian: bool = False,
    ) -> Tuple[np.ndarray, Optional[np.ndarray]]:
        fixed = _fixed_indices_from_constraints(atoms)
        active_atoms = np.asarray([i for i in range(len(atoms)) if i not in fixed], dtype=int)

        if active_atoms.size == 0:
            if return_partial_hessian:
                return np.zeros((0, 0), dtype=dtype), active_atoms
            return np.zeros((3 * len(atoms), 3 * len(atoms)), dtype=dtype), None

        if info_path is not None:
            dir_ = os.path.dirname(info_path)
            if dir_:
                os.makedirs(dir_, exist_ok=True)
            with open(info_path, "w", encoding="utf-8") as log:
                log.write("Analytical Hessian (hessian_ff)\n")
                log.write("--------------------------------\n")
                log.write(f"n_active_atoms = {active_atoms.size}\n")
                log.flush()

        xyz = self._positions_to_tensor(atoms.get_positions())
        h_local, _ = build_analytical_hessian(
            system=self.system,
            coords=xyz,
            active_atoms=active_atoms.tolist(),
        )
        h_sub = h_local.detach().cpu().numpy().astype(np.float64, copy=False) * KCALMOL2EV
        h_sub = np.asarray(h_sub, dtype=dtype)
        # Release the GPU tensor after the required NumPy conversion.
        del h_local
        if torch.cuda.is_available():
            torch.cuda.empty_cache()

        if return_partial_hessian:
            return h_sub, active_atoms

        h_full = _expand_partial_hessian(h_sub, active_atoms, len(atoms), dtype=dtype)
        return h_full, None

    def finite_difference_hessian(
        self,
        atoms: Atoms,
        *,
        delta: float = 1e-3,
        info_path: Optional[str] = None,
        dtype: np.dtype = np.float64,
        return_partial_hessian: bool = False,
    ) -> Tuple[np.ndarray, Optional[np.ndarray]]:
        from mlmm.io.hessian_calc import hessian_calc

        H_full = hessian_calc(
            atoms,
            self,
            delta=float(delta),
            info_path=info_path,
            dtype=dtype,
        )
        if return_partial_hessian:
            fixed = _fixed_indices_from_constraints(atoms)
            active_atoms = np.asarray(
                [i for i in range(len(atoms)) if i not in fixed], dtype=int
            )
            idx3 = np.concatenate(
                [3 * active_atoms + component for component in range(3)]
            )
            idx3.sort()
            return H_full[np.ix_(idx3, idx3)], active_atoms
        return H_full, None



class OpenMMCalculator(Calculator):
    """
    ASE Calculator wrapper for OpenMM backend (finite-difference Hessian).

    This calculator uses OpenMM for MM force field evaluation and supports
    CUDA/CPU platforms. Unlike hessianffCalculator, it computes Hessians
    via numerical finite differences.

    Parameters
    ----------
    parm7 : str
        Path to Amber parm7 topology file.
    rst7 : str
        Path to Amber rst7 coordinate file.
    device : str, default "auto"
        Platform selection: "auto", "cuda", or "cpu".
    cuda_idx : int, default 0
        CUDA device index when device="cuda".
    threads : int, default 16
        Number of CPU threads when device="cpu".
    """

    implemented_properties = ["energy", "forces"]

    def __init__(
        self,
        parm7: str,
        rst7: str,
        *,
        device: str = "auto",
        cuda_idx: int = 0,
        threads: int = 16,
        **kwargs,
    ):
        super().__init__(**kwargs)

        if not HAS_OPENMM:
            raise ImportError(
                "OpenMM is required for OpenMMCalculator. "
                "Install with: conda install -c conda-forge openmm"
            )

        # Auto-detect device
        if device == "auto":
            registered = {
                Platform.getPlatform(index).getName().upper()
                for index in range(Platform.getNumPlatforms())
            }
            device = "cuda" if "CUDA" in registered else "cpu"

        # Expose the resolved device for MLMMCore.compute's ML/MM parallel gate.
        self.device = device

        # Platform selection
        if device == "cuda":
            platform = Platform.getPlatformByName("CUDA")
            properties = {
                "CudaDeviceIndex": str(cuda_idx),
                "CudaPrecision": "double",
                "DeterministicForces": "true",
                "CudaUseBlockingSync": "true",
            }
        else:
            platform = Platform.getPlatformByName("CPU")
            properties = {"Threads": str(threads)}

        # Load Amber topology and coordinates
        self.prmtop = app.AmberPrmtopFile(parm7)
        inpcrd = app.AmberInpcrdFile(rst7)

        # Create OpenMM system and context
        self.system = self.prmtop.createSystem(
            nonbondedMethod=app.NoCutoff,
            rigidWater=False
        )
        self.integrator = mm.VerletIntegrator(0 * unit.femtoseconds)
        self.context = mm.Context(self.system, self.integrator, platform, properties)
        self.context.setPositions(inpcrd.positions)
        # Virtual sites (4-point water OPC/TIP4P EPW, lone pairs): OpenMM's getForces() does NOT
        # redistribute their force to the parent atoms outside MD integration, so an external
        # optimizer would treat them as spurious DOF (stalls). We zero their force each eval and
        # rely on computeVirtualSites() to keep their positions correct from the parents.
        self._vsite_idx = [i for i in range(self.system.getNumParticles())
                           if self.system.isVirtualSite(i)]

    def calculate(self, atoms: Atoms = None, properties=None, system_changes=all_changes):
        """Compute energy and forces for the given atoms."""
        super().calculate(atoms, properties, system_changes)

        # Define eV unit for OpenMM
        ev_base_unit = ScaledUnit(1.602176634e-19, joule, "electron volt", "eV")
        eV = unit.Unit({ev_base_unit: 1.0})

        # Update positions and get state
        self.context.setPositions(atoms.get_positions() * unit.angstrom)
        # Place virtual sites (4-point water OPC/TIP4P "EPW", lone pairs) from their parent atoms,
        # so the MM point charge sits correctly even when the parent water moves. Without this the
        # optimizer's (inert) EP coordinate would be used verbatim. No-op if there are no v-sites.
        self.context.computeVirtualSites()
        need_forces = properties is None or "forces" in properties
        state = self.context.getState(
            getEnergy=True, getForces=need_forces
        )

        # Extract energy and, only when requested, forces in eV units.
        energy = state.getPotentialEnergy().value_in_unit(eV / unit.item)
        self.results = {"energy": energy}
        if need_forces:
            forces = state.getForces(asNumpy=True).value_in_unit(
                eV / unit.angstrom / unit.item
            )
            # Zero external virtual-site forces so the optimizer never moves
            # coordinates overwritten by computeVirtualSites(). OpenMM has
            # already redistributed their force to the parent atoms.
            if self._vsite_idx:
                forces[self._vsite_idx] = 0.0
            self.results["forces"] = forces

    def finite_difference_hessian(
        self,
        atoms: Atoms,
        *,
        delta: float = 0.01,
        info_path: Optional[str] = None,
        dtype: np.dtype = np.float64,
        return_partial_hessian: bool = False,
    ) -> Tuple[np.ndarray, Optional[np.ndarray]]:
        """
        Compute Hessian via finite differences using hessian_calc utility.

        Parameters
        ----------
        atoms : Atoms
            Structure to differentiate.
        delta : float, default 0.01
            Displacement size in Angstrom.
        info_path : str | None
            Progress log file path.
        dtype : numpy dtype, default float64
            Data type for the Hessian matrix.
        return_partial_hessian : bool, default False
            If True, return only the active sub-Hessian and active atom indices.

        Returns
        -------
        H_full : ndarray
            Full (3N, 3N) Hessian matrix in eV/Å².
        active_atoms : ndarray | None
            Active atom indices (only if return_partial_hessian=True).
        """
        from mlmm.io.hessian_calc import hessian_calc

        H_full = hessian_calc(atoms, self, delta=delta, info_path=info_path, dtype=dtype)

        if return_partial_hessian:
            fixed = _fixed_indices_from_constraints(atoms)
            excluded = fixed | set(self._vsite_idx)
            active_atoms = np.asarray(
                [i for i in range(len(atoms)) if i not in excluded],
                dtype=np.int64,
            )
            if active_atoms.size == 0:
                return np.zeros((0, 0), dtype=dtype), active_atoms
            # Extract active sub-Hessian to match hessian_ff convention
            idx3 = np.concatenate([3 * active_atoms + d for d in range(3)])
            idx3.sort()
            H_sub = H_full[np.ix_(idx3, idx3)]
            return H_sub, active_atoms

        return H_full, None


#                          ML/MM Core (Multi-Backend)

@dataclass(frozen=True)
class _MLHighOut:
    E: float
    F: Optional[np.ndarray]
    H: Optional[torch.Tensor]
    timing: Dict[str, float | str]


@dataclass(frozen=True)
class _MMLowOut:
    E_real: float
    F_real: Optional[np.ndarray]
    E_model: float
    F_model: Optional[np.ndarray]
    H_real: Optional[np.ndarray]
    H_model: Optional[np.ndarray]
    active_atoms_from_fd: Optional[np.ndarray]
    timing: Dict[str, float | str]


def validate_parmed_atom_order(
    input_structure,
    real_topology,
    *,
    input_label: str = "input structure",
    topology_label: str = "parm7 topology",
) -> None:
    """Validate count and available atom identity before positional assignment."""

    def amber_name(name: str) -> str:
        """Normalize the PDB/Amber placement of a leading hydrogen digit."""
        return name[1:] + name[0] if name[:1].isdigit() else name

    n_input = int(len(input_structure.atoms))
    n_top = int(len(real_topology.atoms))
    if n_input != n_top:
        raise ValueError(
            "Atom-count mismatch between input structure and real topology: "
            f"{input_label!r} has {n_input} atoms, {topology_label!r} expects "
            f"{n_top} atoms. Provide a full-system structure consistent with the parm7."
        )
    for index, (input_atom, top_atom) in enumerate(
        zip(input_structure.atoms, real_topology.atoms)
    ):
        input_z = int(getattr(input_atom, "atomic_number", 0) or 0)
        top_z = int(getattr(top_atom, "atomic_number", 0) or 0)
        # Unknown PDB elements cannot prove identity; known unequal elements do
        # prove that positional coordinate assignment would be unsafe.
        if input_z > 0 and top_z > 0 and input_z != top_z:
            from ase.data import chemical_symbols

            raise ValueError(
                "Atom-order mismatch between input structure and parm7 at "
                f"atom {index + 1} (0-based {index}): input has "
                f"{chemical_symbols[input_z]}, topology has "
                f"{chemical_symbols[top_z]}. Regenerate the topology from the "
                "same atom ordering or reorder the input structure."
            )
        input_name = str(getattr(input_atom, "name", "") or "").strip()
        top_name = str(getattr(top_atom, "name", "") or "").strip()
        input_residue = getattr(input_atom, "residue", None)
        top_residue = getattr(top_atom, "residue", None)
        input_resname = str(
            getattr(input_residue, "name", "") or ""
        ).strip()
        top_resname = str(getattr(top_residue, "name", "") or "").strip()
        input_residx = getattr(input_residue, "idx", None)
        top_residx = getattr(top_residue, "idx", None)
        residue_mismatch = bool(
            input_resname
            and top_resname
            and input_resname != top_resname
        )
        ordinal_mismatch = bool(
            input_residx is not None
            and top_residx is not None
            and int(input_residx) != int(top_residx)
        )
        names_match = bool(
            input_name == top_name
            or (
                input_z == 1
                and top_z == 1
                and amber_name(input_name) == amber_name(top_name)
            )
        )
        name_mismatch = bool(input_name and top_name and not names_match)
        if name_mismatch or residue_mismatch or ordinal_mismatch:
            raise ValueError(
                "Atom-order mismatch between input structure and parm7 at "
                f"atom {index + 1} (0-based {index}): input identity is "
                f"{input_resname or '?'}[{input_residx!s}]:{input_name or '?'}, "
                f"topology identity is "
                f"{top_resname or '?'}[{top_residx!s}]:{top_name or '?'}. "
                "Regenerate the topology from the same atom ordering or reorder "
                "the input structure."
            )


class MLMMCore:
    """ONIOM-like ML/MM engine supporting multiple MLIP backends.

    Supported ML backends: UMA (default), ORB, MACE, AIMNet2.
    Supported MM backends: hessian_ff (analytical), OpenMM (FD).
    The optional ``embedcharge`` path applies an experimental xTB point-charge
    correction for ML/MM environmental effects.
    """

    def __init__(
        self,
        *,
        input_pdb: Optional[str] = None,
        coordinate_path: Optional[str] = None,
        real_parm7: Optional[str] = None,
        model_pdb: Optional[str] = None,
        model_charge: Optional[int] = 0,
        model_mult: int = 1,
        link_mlmm: List[Tuple[str, str]] | None = None,
        link_atom_method: str = "scaled",
        # ML backend selection
        backend: str = "uma",
        uma_model: str = DEFAULT_UMA_MODEL,
        uma_task_name: str = "omol",
        uma_precision: str = "fp32",
        workers: int = 1,
        workers_per_node: int = 1,
        orb_model: str = "orb_v3_conservative_omol",
        orb_precision: str = "float64",
        mace_model: str = "MACE-OMOL-0",
        mace_dtype: str = "float64",
        aimnet2_model: str = "aimnet2",
        # Custom ML backend from a user Python file (--calc-file)
        calc_file: Optional[str] = None,
        calc_factory: str = "get_calculator",
        # MM settings
        mm_fd: bool = True,
        mm_hessian_mode: Optional[str] = None,
        mm_fd_dir: Optional[str] = None,
        mm_fd_delta: float = 1e-3,
        symmetrize_hessian: bool = True,
        print_timing: bool = True,
        print_vram: bool = True,
        H_double: bool = True,
        ml_device: str = "auto",
        ml_cuda_idx: int = 0,
        mm_backend: str = "hessian_ff",
        mm_device: str = "cpu",
        mm_cuda_idx: int = 0,
        mm_threads: int = 16,
        freeze_atoms: List[int] | None = None,
        hessian_calc_mode: str = "FiniteDifference",
        return_partial_hessian: bool = True,
        hess_cutoff: Optional[float] = None,
        movable_cutoff: Optional[float] = None,
        use_bfactor_layers: bool = True,       # matches MLMM_CALC_KW default
        hess_mm_atoms: Optional[List[int]] = None,
        movable_mm_atoms: Optional[List[int]] = None,
        frozen_mm_atoms: Optional[List[int]] = None,
        # Point-charge embedding correction
        embedcharge: bool = False,
        embedcharge_step: float = 1.0e-3,
        embedcharge_cutoff: float = 12.0,
        xtb_cmd: str = "xtb",
        xtb_acc: float = 0.2,
        xtb_workdir: str = "tmp",
        xtb_keep_files: bool = False,
        xtb_ncores: int = 4,
        use_cmap: bool = True,
        _high_level_backend: Optional[Any] = None,
        _skip_high_level_backend: bool = False,
        **kwargs,
    ):
        # --- v0.1.x backward compatibility aliases ---
        if "real_pdb" in kwargs:
            warnings.warn("'real_pdb' is deprecated; use 'input_pdb'.", DeprecationWarning, stacklevel=2)
            if input_pdb is None:
                input_pdb = kwargs.pop("real_pdb")
            else:
                kwargs.pop("real_pdb")
        for _old_name in ("real_rst7", "vib_run", "vib_dir"):
            if _old_name in kwargs:
                warnings.warn(f"'{_old_name}' is no longer used and will be ignored.", DeprecationWarning, stacklevel=2)
                kwargs.pop(_old_name)
        if kwargs:
            raise TypeError(f"MLMMCore.__init__() got unexpected keyword arguments: {', '.join(kwargs)}")
        if input_pdb is None:
            raise TypeError("MLMMCore.__init__() missing required keyword argument: 'input_pdb'")
        # Canonicalize the two numerical-method enums before any temporary
        # directory, file copy, topology preparation, or model allocation.
        hessian_calc_mode = normalize_hessian_calc_mode(hessian_calc_mode)
        mm_hessian_mode = normalize_mm_hessian_mode(
            mm_hessian_mode, mm_fd=mm_fd
        )
        link_atom_method = normalize_link_atom_method(link_atom_method)
        if (
            _high_level_backend is None
            and not _skip_high_level_backend
            and int(workers or 1) > 1
            and hessian_calc_mode == "Analytical"
        ):
            raise ValueError(
                "Analytical Hessian cannot be combined with workers>1: the "
                "parallel UMA predictor exposes no autograd model. Use workers=1 "
                "or select hessian_calc_mode='FiniteDifference'."
            )

        # ── Workspace setup ───────────────────────────────────────────────
        # The constructor copies input.pdb / real.parm7 / model.pdb into a
        # private TemporaryDirectory, prepares the ML region selection +
        # link-atom geometry, then loads the ML backend. Heavy ML weights
        # are only allocated AFTER the cheap structural checks
        # (atom-count match, charge/spin parity, layer detection) so a
        # bad input fails in milliseconds rather than seconds.
        self._tmpdir_obj = tempfile.TemporaryDirectory()
        self.tmpdir: str = self._tmpdir_obj.name
        for src, dst in [(input_pdb, "input.pdb"), (real_parm7, "real.parm7"), (model_pdb, "model.pdb")]:
            shutil.copy(src, os.path.join(self.tmpdir, dst))

        self.input_pdb = os.path.join(self.tmpdir, "input.pdb")
        self.real_parm7 = os.path.join(self.tmpdir, "real.parm7")
        self.real_rst7 = os.path.join(self.tmpdir, "real.rst7")
        self.model_pdb = os.path.join(self.tmpdir, "model.pdb")
        self.model_parm7 = os.path.join(self.tmpdir, "model.parm7")
        self.model_rst7 = os.path.join(self.tmpdir, "model.rst7")

        try:
            real_top = pmd.load_file(self.real_parm7)
        except pmd.exceptions.FormatNotFound as exc:
            raise ValueError(
                f"--parm '{self.real_parm7}' is not a valid Amber parm7 file "
                f"(parmed: {exc}). Regenerate with `mlmm mm-parm` or tleap."
            ) from exc
        start_struct = pmd.load_file(self.input_pdb)
        validate_parmed_atom_order(
            start_struct,
            real_top,
            input_label=f"input_pdb={input_pdb}",
            topology_label=f"real_parm7={real_parm7}",
        )
        coordinate_atoms: Optional[Atoms] = None
        if (
            coordinate_path is not None
            and Path(coordinate_path).resolve() != Path(input_pdb).resolve()
        ):
            coordinate_atoms = read(str(coordinate_path), index=0)
            if len(coordinate_atoms) != len(real_top.atoms):
                raise ValueError(
                    "Coordinate input and Amber topology have different atom "
                    f"counts: {len(coordinate_atoms)} != {len(real_top.atoms)}."
                )
            coordinate_numbers = np.asarray(
                coordinate_atoms.get_atomic_numbers(), dtype=int
            )
            topology_numbers = np.asarray(
                [int(atom.atomic_number or 0) for atom in real_top.atoms],
                dtype=int,
            )
            known = (coordinate_numbers > 0) & (topology_numbers > 0)
            mismatch = np.flatnonzero(
                known & (coordinate_numbers != topology_numbers)
            )
            if mismatch.size:
                index = int(mismatch[0])
                raise ValueError(
                    "Atom-order mismatch between coordinate input and parm7 at "
                    f"atom {index + 1} (0-based {index})."
                )
            real_top.coordinates = coordinate_atoms.get_positions()
        else:
            real_top.coordinates = start_struct.coordinates
        real_top.box = None
        apply_cmap_policy(real_top, use_cmap)
        real_top.save(self.real_parm7, overwrite=True)
        real_top.save(self.real_rst7, overwrite=True)

        self.link_mlmm = link_mlmm
        self.link_atom_method = link_atom_method
        self.use_cmap = use_cmap
        self.ml_ID, self.mlmm_links, self._link_elem_pairs = self._ml_prep(real_top)
        link_source = "manual link_mlmm" if self.link_mlmm is not None else "parm7 bonds"
        logger.info(
            "[MLMMCore] ML region = %d atoms from model_pdb; link H = %d "
            "boundary bond(s) from %s",
            len(self.ml_ID),
            len(self.mlmm_links),
            link_source,
        )
        if self.mlmm_links:
            logger.info(
                "[MLMMCore] Link-H boundary pairs (1-based full-system "
                "parm7 ML-MM): %s",
                ", ".join(f"{ml_idx}-{mm_idx}" for ml_idx, mm_idx in self.mlmm_links),
            )
        if self.link_atom_method == "scaled":
            self._link_g_factors = [
                _get_g_factor(qm_e, mm_e, "H") for qm_e, mm_e in self._link_elem_pairs
            ]
        else:
            self._link_g_factors = []
        self.selection_indices = self._mk_model_parm7()

        self.hess_cutoff = hess_cutoff
        self.movable_cutoff = movable_cutoff
        self.use_bfactor_layers = use_bfactor_layers
        self._original_input_pdb = input_pdb
        self._explicit_hess_mm_atoms = hess_mm_atoms
        self._explicit_movable_mm_atoms = movable_mm_atoms
        self._explicit_frozen_mm_atoms = frozen_mm_atoms
        self._compute_layer_indices(real_top.coordinates)

        self.freeze_atoms = [] if freeze_atoms is None else list(freeze_atoms)
        if self.frozen_layer_indices:
            self.freeze_atoms = sorted(set(self.freeze_atoms) | set(self.frozen_layer_indices))

        hess_set = set(self.hess_indices)
        all_atoms = set(range(len(real_top.atoms)))
        self.hess_freeze_atoms = sorted(all_atoms - hess_set)

        self.return_partial_hessian = bool(return_partial_hessian)

        self._n_real = len(real_top.atoms)
        self._idx_map_real_to_model = {idx: pos for pos, idx in enumerate(self.selection_indices)}
        self._update_active_dof_mappings()

        self.H_double = bool(H_double)
        self.H_dtype = torch.float64 if self.H_double else torch.float32
        self.H_np_dtype = np.float64 if self.H_double else np.float32

        self.mm_hessian_mode = mm_hessian_mode
        self.mm_fd = mm_hessian_mode == "finite_difference"
        self.mm_fd_dir = mm_fd_dir
        self.mm_fd_delta = mm_fd_delta
        self.symmetrize_hessian = symmetrize_hessian
        self.print_timing = bool(print_timing)
        self.print_vram = bool(print_vram)
        if self.mm_fd_dir and not os.path.exists(self.mm_fd_dir):
            os.makedirs(self.mm_fd_dir, exist_ok=True)

        if ml_device == "auto":
            ml_device = "cuda" if torch.cuda.is_available() else "cpu"
        if ml_device not in ("cuda", "cpu"):
            # Reject unknown devices to prevent unintended CPU execution.
            raise ValueError(
                "ml_device must be 'auto', 'cuda' or 'cpu' (choose the GPU with ml_cuda_idx), "
                f"got {ml_device!r}"
            )
        self.device_str = ml_device
        self.ml_device = torch.device(f"cuda:{ml_cuda_idx}" if ml_device == "cuda" else "cpu")

        self.model_charge = int(0 if model_charge is None else model_charge)
        self.model_mult = int(model_mult)
        if self.model_mult < 1:
            raise ValueError(
                f"Spin multiplicity must be >= 1, got {self.model_mult}."
            )
        logger.info(f"[MLMMCore] ML-region net charge = {self.model_charge}")
        if _skip_high_level_backend:
            self.backend_name = "none"
        elif _high_level_backend is not None:
            self.backend_name = str(
                getattr(_high_level_backend, "name", "external")
            ).strip().lower()
        elif backend is not None:
            self.backend_name = str(backend).strip().lower()
        else:
            self.backend_name = "uma"

        # Validate charge/spin parity before loading the ML model.
        # selection_indices is parmed 0-based (a.idx in _mk_model_parm7); index
        # real_top.atoms directly without subtracting 1.
        if not _skip_high_level_backend:
            from mlmm.core.utils import validate_charge_spin as _vcs
            from ase.data import chemical_symbols as _chem_sym
            _ml_elements = [
                _chem_sym[int(real_top.atoms[i].atomic_number)]
                for i in self.selection_indices
            ]
            # Each link atom adds one H (+1 electron); validate including them.
            _ml_elements.extend(["H"] * len(self.mlmm_links))
            _vcs(_ml_elements, self.model_charge, self.model_mult, source=model_pdb)
            self._charge_check_done = True
        else:
            self._charge_check_done = False

        # DFT and other energy-only high-level methods can replace only this
        # evaluator while retaining the ordinary ML/MM preparation, MM
        # calculators, and subtractive recombination.
        self._ml_backend = _high_level_backend
        if self._ml_backend is None and not _skip_high_level_backend:
            self._ml_backend = _create_ml_backend(
                self.backend_name,
                uma_model=uma_model,
                uma_task_name=uma_task_name,
                uma_precision=uma_precision,
                workers=workers,
                workers_per_node=workers_per_node,
                orb_model=orb_model,
                orb_precision=orb_precision,
                mace_model=mace_model,
                mace_dtype=mace_dtype,
                aimnet2_model=aimnet2_model,
                calc_file=calc_file,
                calc_factory=calc_factory,
                model_charge=self.model_charge,
                model_mult=self.model_mult,
                ml_device=self.ml_device,
                analytical_hessian=(hessian_calc_mode == "Analytical"),
            )

        # Point-charge embedding correction
        self.embedcharge = bool(embedcharge)
        self.embedcharge_cutoff = embedcharge_cutoff
        self._embed_correction: Optional[_EmbedChargeCorrection] = None
        if self.embedcharge:
            self._embed_correction = _EmbedChargeCorrection(
                xtb_cmd=xtb_cmd,
                xtb_acc=xtb_acc,
                xtb_workdir=xtb_workdir,
                xtb_keep_files=xtb_keep_files,
                xtb_ncores=xtb_ncores,
                hessian_step=embedcharge_step,
            )

        # MM backend selection: hessian_ff or openmm
        self.mm_backend = str(mm_backend).strip().lower()
        if self.mm_backend == "openmm":
            self.calc_real_low = OpenMMCalculator(
                parm7=self.real_parm7, rst7=self.real_rst7,
                device=mm_device, cuda_idx=mm_cuda_idx, threads=mm_threads
            )
            self.calc_model_low = OpenMMCalculator(
                parm7=self.model_parm7, rst7=self.model_rst7,
                device=mm_device, cuda_idx=mm_cuda_idx, threads=mm_threads
            )
        elif self.mm_backend == "hessian_ff":
            self.calc_real_low = hessianffCalculator(
                parm7=self.real_parm7, rst7=None,
                device=mm_device, cuda_idx=mm_cuda_idx, threads=mm_threads
            )
            self.calc_model_low = hessianffCalculator(
                parm7=self.model_parm7, rst7=None,
                device=mm_device, cuda_idx=mm_cuda_idx, threads=mm_threads
            )
        else:
            raise ValueError(
                f"Unknown mm_backend '{mm_backend}'. Choose 'hessian_ff' or 'openmm'."
            )

        self.hessian_calc_mode = hessian_calc_mode
        self._ml_hessian_mode = (
            "analytical" if hessian_calc_mode == "Analytical" else "fd"
        )

        self._atoms_real_tpl = (
            coordinate_atoms.copy()
            if coordinate_atoms is not None
            else read(self.input_pdb)
        )
        self._atoms_model_tpl = read(self.model_pdb)
        tmp = self._atoms_model_tpl.copy()
        for _ in self.mlmm_links:
            tmp += Atoms("H", positions=[[0.0, 0.0, 0.0]])
        self._atoms_model_LH_tpl = tmp

    def cleanup(self):
        """Clean up temporary directory."""
        if hasattr(self, '_tmpdir_obj') and self._tmpdir_obj is not None:
            try:
                self._tmpdir_obj.cleanup()
            except Exception:
                logger.debug("Failed to clean up tmpdir", exc_info=True)

    def __del__(self):
        self.cleanup()

    def _ml_prep(
        self,
        real_topology,
    ) -> Tuple[List[str], List[Tuple[int, int]], List[Tuple[str, str]]]:
        """Return (ml_ID, mlmm_links, link_elem_pairs)."""
        resolved = resolve_mlmm_atoms(
            self.input_pdb,
            self.model_pdb,
            self.link_mlmm,
            topology=real_topology,
        )
        return (
            [str(idx) for idx in resolved.model_indices],
            list(resolved.link_pairs),
            list(resolved.link_element_pairs),
        )

    def _mk_model_parm7(self) -> List[int]:
        real = pmd.load_file(self.real_parm7, self.real_rst7)
        real.box = None
        ml_atoms = [real.atoms[int(i) - 1] for i in self.ml_ID]
        selection = [a.idx for a in ml_atoms]
        write_model_parm7(
            real,
            selection,
            self.real_parm7,
            self.real_rst7,
            self.model_parm7,
            self.model_rst7,
            self.use_cmap,
        )
        return selection

    def _compute_layer_indices(self, coords: np.ndarray) -> None:
        self.ml_indices = sorted(self.selection_indices)

        n_atoms = int(coords.shape[0])
        all_indices = set(range(n_atoms))
        mm_indices = all_indices - set(self.ml_indices)

        def min_dist_to_ml(atom_idx: int) -> float:
            atom_coord = coords[atom_idx]
            dists = np.linalg.norm(coords[self.ml_indices] - atom_coord, axis=1)
            return float(np.min(dists))

        has_explicit = (
            self._explicit_hess_mm_atoms is not None
            or self._explicit_movable_mm_atoms is not None
            or self._explicit_frozen_mm_atoms is not None
        )
        if has_explicit:
            explicit_hess = set(self._explicit_hess_mm_atoms or [])
            explicit_movable = set(self._explicit_movable_mm_atoms or [])
            explicit_frozen = set(self._explicit_frozen_mm_atoms or [])

            for idx_set, name in [
                (explicit_hess, "hess_mm_atoms"),
                (explicit_movable, "movable_mm_atoms"),
                (explicit_frozen, "frozen_mm_atoms"),
            ]:
                for idx in idx_set:
                    if idx < 0 or idx >= n_atoms:
                        raise ValueError(f"Invalid atom index {idx} in {name}: must be 0 <= idx < {n_atoms}")
                    if idx in self.ml_indices:
                        raise ValueError(f"Atom index {idx} in {name} is also in ML region (model_pdb)")

            self.hess_mm_indices = sorted(explicit_hess & mm_indices)
            self.movable_mm_indices = sorted(explicit_movable & mm_indices)
            self.frozen_layer_indices = sorted(explicit_frozen & mm_indices)

            assigned_mm = explicit_hess | explicit_movable | explicit_frozen
            unassigned_mm = mm_indices - assigned_mm
            self.movable_mm_indices = sorted(set(self.movable_mm_indices) | unassigned_mm)

            self.hess_indices = sorted(self.ml_indices + self.hess_mm_indices)
            self.movable_indices = sorted(self.ml_indices + self.hess_mm_indices + self.movable_mm_indices)
            return

        if self.use_bfactor_layers:
            from mlmm.core.utils import read_bfactors_from_pdb, parse_layer_indices_from_bfactors, has_valid_layer_bfactors
            from pathlib import Path

            bfactors = read_bfactors_from_pdb(Path(self._original_input_pdb))
            if has_valid_layer_bfactors(bfactors):
                layer_info = parse_layer_indices_from_bfactors(bfactors)

                movable_from_layer = set(layer_info["movable_mm_indices"]) & mm_indices
                frozen_from_layer = set(layer_info["frozen_indices"]) & mm_indices
                hess_from_layer = set(layer_info["hess_mm_indices"]) & mm_indices

                if self.movable_cutoff is not None:
                    movable_pool = {
                        idx
                        for idx in mm_indices
                        if min_dist_to_ml(idx) <= float(self.movable_cutoff)
                    }
                    frozen_from_layer = mm_indices - movable_pool
                else:
                    # Unassigned MM atoms default to movable.
                    assigned_mm = (
                        movable_from_layer | frozen_from_layer | hess_from_layer
                    )
                    unassigned_mm = mm_indices - assigned_mm
                    movable_pool = set(movable_from_layer) | set(unassigned_mm)

                # Hessian-target MM selection:
                hess_mm: set[int]
                if self.hess_cutoff is not None:
                    hess_cut = float(self.hess_cutoff)
                    hess_mm = {idx for idx in movable_pool if min_dist_to_ml(idx) <= hess_cut}
                else:
                    hess_mm = set(hess_from_layer) & movable_pool

                movable_mm = movable_pool - hess_mm

                self.hess_mm_indices = sorted(hess_mm)
                self.movable_mm_indices = sorted(movable_mm)
                self.frozen_layer_indices = sorted(frozen_from_layer)

                self.hess_indices = sorted(self.ml_indices + self.hess_mm_indices)
                self.movable_indices = sorted(self.ml_indices + self.hess_mm_indices + self.movable_mm_indices)
                return

        if self.hess_cutoff is None and self.movable_cutoff is None:
            self.hess_mm_indices = sorted(mm_indices)
            self.movable_mm_indices = []
            self.frozen_layer_indices = []
            self.hess_indices = sorted(self.ml_indices + self.hess_mm_indices)
            self.movable_indices = sorted(self.ml_indices + self.hess_mm_indices)
            return

        mov_cut = self.movable_cutoff if self.movable_cutoff is not None else float("inf")

        movable_pool = {
            idx for idx in mm_indices if min_dist_to_ml(idx) <= mov_cut
        }
        frozen_mm = mm_indices - movable_pool
        if self.hess_cutoff is None:
            # The default Hessian target is every atom that remains movable;
            # it must not include atoms frozen by movable_cutoff.
            hess_mm = set(movable_pool)
        else:
            hess_cut = float(self.hess_cutoff)
            hess_mm = {
                idx
                for idx in movable_pool
                if min_dist_to_ml(idx) <= hess_cut
            }
        movable_mm = movable_pool - hess_mm

        self.hess_mm_indices = sorted(hess_mm)
        self.movable_mm_indices = sorted(movable_mm)
        self.frozen_layer_indices = sorted(frozen_mm)
        self.hess_indices = sorted(self.ml_indices + self.hess_mm_indices)
        self.movable_indices = sorted(self.ml_indices + self.hess_mm_indices + self.movable_mm_indices)

    def _update_active_dof_mappings(self) -> None:
        freeze_set = set(self.freeze_atoms)
        self.active_atoms_real = [i for i in range(self._n_real) if i not in freeze_set]
        self.n_active_real = len(self.active_atoms_real)
        self.full_to_active_real = {a: i for i, a in enumerate(self.active_atoms_real)}
        self.active_to_full_real = {i: a for i, a in enumerate(self.active_atoms_real)}

        # Hessian representation follows both the requested Hessian target and
        # the effective geometry constraints. An explicitly frozen atom must
        # never reappear merely because it lies inside the Hessian cutoff.
        hess_freeze_set = set(self.hess_freeze_atoms) | freeze_set
        self.effective_hess_freeze_atoms = sorted(hess_freeze_set)
        self.hess_active_atoms = [i for i in range(self._n_real) if i not in hess_freeze_set]
        self.n_hess_active = len(self.hess_active_atoms)
        self.full_to_hess_active = {a: i for i, a in enumerate(self.hess_active_atoms)}
        self.hess_active_to_full = {i: a for i, a in enumerate(self.hess_active_atoms)}

        self.ml_hess_active_indices = [
            self.full_to_hess_active[i] for i in self.selection_indices if i in self.full_to_hess_active
        ]

        self.freeze_model = [
            self._idx_map_real_to_model[i] for i in self.freeze_atoms if i in self._idx_map_real_to_model
        ]

    def _build_within_partial_hessian(self) -> Dict[str, np.ndarray | int | str]:
        """Build metadata for a partial (Hessian-target-only) Hessian."""
        n_real = int(self._n_real)
        active_atoms = np.asarray(self.hess_active_atoms, dtype=int)
        active_n_atoms = int(active_atoms.size)

        active_dofs = np.empty(active_n_atoms * 3, dtype=int)
        for i, a in enumerate(active_atoms):
            base = 3 * int(a)
            active_dofs[3 * i:3 * i + 3] = (base, base + 1, base + 2)

        full_to_active = -np.ones(n_real, dtype=int)
        if active_n_atoms:
            full_to_active[active_atoms] = np.arange(active_n_atoms, dtype=int)

        return {
            "kind": "hess-target-only",
            "active_atoms": active_atoms,
            "active_dofs": active_dofs,
            "active_to_full": active_atoms.copy(),
            "full_to_active": full_to_active,
            "full_n_atoms": n_real,
            "full_n_dof": int(3 * n_real),
            "active_n_atoms": active_n_atoms,
            "active_n_dof": int(3 * active_n_atoms),
        }

    def _finalize_result_constraints(self, results: Dict) -> Dict:
        """Apply the one final force/Hessian mask after all corrections."""

        force_constrained = sorted(
            {int(i) for i in self.freeze_atoms if 0 <= int(i) < self._n_real}
        )
        forces = results.get("forces")
        if forces is not None and force_constrained:
            forces_array = np.asarray(forces).reshape(self._n_real, 3)
            forces_array[np.asarray(force_constrained, dtype=int), :] = 0.0

        hessian = results.get("hessian")
        if hessian is None:
            return results

        effective_atoms = list(self.effective_hess_freeze_atoms)
        effective_dofs = [
            3 * atom + axis for atom in effective_atoms for axis in range(3)
        ]
        full_n_dof = 3 * self._n_real
        n_values = int(np.prod(tuple(hessian.shape)))
        if n_values == full_n_dof * full_n_dof:
            local_dofs = effective_dofs
        else:
            metadata = results.get("within_partial_hessian")
            if not isinstance(metadata, dict) or "active_dofs" not in metadata:
                raise RuntimeError(
                    "A compact ML/MM Hessian requires active_dofs metadata."
                )
            represented = np.asarray(metadata["active_dofs"], dtype=int).reshape(-1)
            if represented.size * represented.size != n_values:
                raise RuntimeError(
                    "Partial Hessian dimensions do not match active_dofs metadata."
                )
            constrained_set = set(effective_dofs)
            local_dofs = [
                local
                for local, global_dof in enumerate(represented.tolist())
                if global_dof in constrained_set
            ]
            if local_dofs:
                raise RuntimeError(
                    "Partial Hessian metadata includes effectively constrained atoms."
                )

        if local_dofs:
            if isinstance(hessian, torch.Tensor):
                square = hessian.reshape(full_n_dof, full_n_dof)
                idx = torch.as_tensor(
                    local_dofs, dtype=torch.long, device=square.device
                )
                square.index_fill_(0, idx, 0.0)
                square.index_fill_(1, idx, 0.0)
            else:
                square = np.asarray(hessian).reshape(full_n_dof, full_n_dof)
                idx = np.asarray(local_dofs, dtype=int)
                square[idx, :] = 0.0
                square[:, idx] = 0.0
        return results

    def _prep_3_layer_atoms(self, real_coord: np.ndarray):
        atoms_real = self._atoms_real_tpl.copy()
        atoms_real.set_positions(real_coord)

        atoms_model = self._atoms_model_tpl.copy()
        atoms_model_LH = self._atoms_model_LH_tpl.copy()

        for i, ridx in enumerate(self.ml_ID):
            pos = atoms_real[int(ridx) - 1].position
            atoms_model[i].position = pos
            atoms_model_LH[i].position = pos

        added_link_atoms = []
        base_model_len = len(self._atoms_model_tpl)
        for k, (ml_idx, mm_idx) in enumerate(self.mlmm_links):
            ml_i = ml_idx - 1
            mm_i = mm_idx - 1
            vec = atoms_real[mm_i].position - atoms_real[ml_i].position
            R = np.linalg.norm(vec)
            if R < 1e-6:
                # A degenerate ML/MM parent bond cannot define a link-H
                # position and indicates an invalid input geometry.
                raise ValueError(
                    f"Degenerate link distance |r(MM={mm_idx}) - r(ML={ml_idx})| "
                    f"= {R:.3e} A < 1e-6 A. Check the input geometry: the link "
                    f"parent and boundary MM atom share a position."
                )

            if self.link_atom_method == "scaled":
                g = self._link_g_factors[k]
                H_pos = atoms_real[ml_i].position + g * vec
                param = g  # g-factor stored as 4th element
            else:
                ml_elem = atoms_real[ml_i].symbol
                if ml_elem == "C":
                    dist = 1.09
                elif ml_elem == "N":
                    dist = 1.01
                else:
                    raise ValueError(
                        f"Unsupported link parent element: {ml_elem}. Only C and N are supported."
                    )
                u = vec / R
                H_pos = atoms_real[ml_i].position + u * dist
                param = dist  # fixed bond length stored as 4th element

            link_idx_in_model_LH = base_model_len + k
            atoms_model_LH[link_idx_in_model_LH].position = H_pos
            added_link_atoms.append((link_idx_in_model_LH, ml_i, mm_i, param))

        freeze_model: List[int] = []
        if self.freeze_atoms:
            atoms_real.set_constraint(FixAtoms(indices=self.freeze_atoms))
            real_to_model = self._idx_map_real_to_model
            freeze_model = [real_to_model[i] for i in self.freeze_atoms if i in real_to_model]
            if freeze_model:
                atoms_model.set_constraint(FixAtoms(indices=freeze_model))
                atoms_model_LH.set_constraint(FixAtoms(indices=freeze_model))

        if not getattr(self, "_charge_check_done", False):
            from mlmm.core.utils import validate_charge_spin
            validate_charge_spin([a.symbol for a in atoms_model_LH],
                                 self.model_charge, self.model_mult,
                                 source=getattr(self, "model_pdb", None))
            self._charge_check_done = True

        return atoms_real, atoms_model, atoms_model_LH, added_link_atoms, freeze_model

    def prepare_atoms(
        self, coord_ang: Optional[np.ndarray] = None
    ) -> Tuple[Atoms, Atoms, Atoms]:
        """Return the REAL, MODEL, and MODEL+link-H structures prepared by this core."""

        coordinates = (
            self._atoms_real_tpl.get_positions()
            if coord_ang is None
            else np.asarray(coord_ang, dtype=float).reshape(-1, 3)
        )
        atoms_real, atoms_model, atoms_model_lh, _, _ = (
            self._prep_3_layer_atoms(coordinates)
        )
        return atoms_real, atoms_model, atoms_model_lh

    @staticmethod
    def _jacobian_blocks_numpy(r_ml: np.ndarray, r_mm: np.ndarray, dist: float) -> Optional[np.ndarray]:
        """Returns J shape (6, 3): rows=[Q_xyz, M_xyz], cols=L_xyz."""
        vec = r_mm - r_ml
        R = np.linalg.norm(vec)
        if R < 1e-12:
            return None
        u = vec / R
        I = np.eye(3)
        du_dQ = (I - np.outer(u, u)) / R
        dR_dQ = I - dist * du_dQ
        dR_dM = dist * du_dQ
        return np.hstack([dR_dQ, dR_dM]).T

    @staticmethod
    def _jacobian_blocks_torch(
        r_ml: torch.Tensor,
        r_mm: torch.Tensor,
        dist: float,
        *,
        dtype: torch.dtype,
        device: torch.device,
    ) -> Optional[torch.Tensor]:
        """Returns K shape (3, 6): rows=L_xyz, cols=[Q_xyz, M_xyz]."""
        vec = r_mm - r_ml
        Rlen = torch.norm(vec)
        if float(Rlen) < 1e-12:
            return None
        u = vec / Rlen
        I3 = torch.eye(3, dtype=dtype, device=device)
        du_dQ = (I3 - torch.outer(u, u)) / Rlen
        dR_dQ = I3 - dist * du_dQ
        dR_dM = dist * du_dQ
        return torch.hstack([dR_dQ, dR_dM])

    @staticmethod
    def _jacobian_blocks_numpy_scaled(g: float) -> np.ndarray:
        """Jacobian for scaled (g-factor) link atoms. Shape (6, 3)."""
        I3 = np.eye(3)
        return np.vstack([(1.0 - g) * I3, g * I3])

    @staticmethod
    def _jacobian_blocks_torch_scaled(
        g: float, *, dtype: torch.dtype, device: torch.device,
    ) -> torch.Tensor:
        """Jacobian for scaled (g-factor) link atoms. Shape (3, 6)."""
        I3 = torch.eye(3, dtype=dtype, device=device)
        return torch.hstack([(1.0 - g) * I3, g * I3])

    def _get_mm_charges(self, atom_indices: Sequence[int]) -> np.ndarray:
        """Retrieve MM partial charges for the given atom indices.

        Works with both hessian_ff (AmberSystem) and OpenMM backends.
        """
        calc = self.calc_real_low
        # hessian_ff: AmberSystem with .charge tensor
        if isinstance(calc, hessianffCalculator) and hasattr(calc, "system"):
            return np.array(
                [calc.system.charge[i].item() for i in atom_indices],
                dtype=np.float64,
            )
        # OpenMM: extract charges from NonbondedForce
        if isinstance(calc, OpenMMCalculator) and HAS_OPENMM:
            sys_omm = calc.system
            for fi in range(sys_omm.getNumForces()):
                force = sys_omm.getForce(fi)
                if force.__class__.__name__ == "NonbondedForce":
                    charges = np.array(
                        [force.getParticleParameters(i)[0].value_in_unit(
                            unit.elementary_charge)
                         for i in atom_indices],
                        dtype=np.float64,
                    )
                    return charges
        # Fallback: zero charges
        warnings.warn(
            "Could not extract MM charges from the calculator; returning zeros. "
            "Embedcharge correction will have no effect.",
            RuntimeWarning,
        )
        return np.zeros(len(atom_indices), dtype=np.float64)

    def _eval_ml_high(
        self,
        atoms_model_LH: Atoms,
        freeze_model: Sequence[int],
        *,
        need_forces: bool,
        return_hessian: bool,
    ) -> _MLHighOut:
        local_timing: Dict[str, float | str] = {}
        if need_forces or return_hessian:
            E_model_high, F_model_high, opaque = self._ml_backend.eval(
                atoms_model_LH, need_grad=True
            )
        else:
            E_model_high = self._ml_backend.energy(atoms_model_LH)
            F_model_high = None
            opaque = None
        local_timing["ml_backend"] = self.backend_name

        H_high = None
        if return_hessian:
            n_mlLH = len(atoms_model_LH)
            if self._ml_hessian_mode == "analytical" and self._ml_backend.supports_analytical_hessian:
                t0 = time.perf_counter()
                H_high = self._ml_backend.hessian_analytical(opaque, n_mlLH, dtype=self.H_dtype)
                local_timing["ml_hessian_mode"] = "Analytical"
                local_timing["ml_hessian_s"] = time.perf_counter() - t0
            else:
                if self._ml_hessian_mode == "analytical" and not self._ml_backend.supports_analytical_hessian:
                    raise RuntimeError(
                        "The requested analytical Hessian is unavailable for "
                        f"the {self.backend_name} backend. Select "
                        "hessian_calc_mode='FiniteDifference' explicitly."
                    )
                t0 = time.perf_counter()
                H_high = self._ml_backend.hessian_fd(
                    atoms_model_LH, freeze_model,
                    eps_ang=1.0e-3, dtype=self.H_dtype, device=self.ml_device,
                )
                local_timing["ml_hessian_mode"] = "FiniteDifference"
                local_timing["ml_hessian_s"] = time.perf_counter() - t0

        return _MLHighOut(E=E_model_high, F=F_model_high, H=H_high, timing=local_timing)

    def _eval_mm_low(
        self,
        atoms_real: Atoms,
        atoms_model: Atoms,
        *,
        need_forces: bool,
        return_hessian: bool,
    ) -> _MMLowOut:
        local_timing: Dict[str, float | str] = {}

        atoms_real.calc = self.calc_real_low
        atoms_model.calc = self.calc_model_low

        if need_forces:
            # Request forces first: both supported MM calculators return energy
            # and forces together, so the following energy reads hit ASE's
            # cache instead of evaluating each layer twice.
            F_real_low = np.double(atoms_real.get_forces())
            F_model_low = np.double(atoms_model.get_forces())
            E_real_low = atoms_real.get_potential_energy()
            E_model_low = atoms_model.get_potential_energy()
        else:
            E_real_low = atoms_real.get_potential_energy()
            E_model_low = atoms_model.get_potential_energy()
            F_real_low = None
            F_model_low = None

        H_real_np = None
        H_model_np = None
        active_atoms_from_fd = None

        mm_hessian_mode = getattr(
            self,
            "mm_hessian_mode",
            "finite_difference" if self.mm_fd else "analytical",
        )
        if return_hessian and mm_hessian_mode != "none":
            info_real = os.path.join(self.mm_fd_dir, "real.log") if self.mm_fd_dir else None
            info_model = os.path.join(self.mm_fd_dir, "model.log") if self.mm_fd_dir else None

            atoms_real_for_hess = atoms_real.copy()
            # Clear any inherited constraints before applying hess-specific ones
            atoms_real_for_hess.set_constraint()
            atoms_real_for_hess.calc = self.calc_real_low
            if self.effective_hess_freeze_atoms:
                atoms_real_for_hess.set_constraint(
                    FixAtoms(indices=self.effective_hess_freeze_atoms)
                )

            method_name = (
                "finite_difference_hessian"
                if mm_hessian_mode == "finite_difference"
                else "analytical_hessian"
            )
            real_method = getattr(self.calc_real_low, method_name, None)
            model_method = getattr(self.calc_model_low, method_name, None)
            if real_method is None or model_method is None:
                raise RuntimeError(
                    f"MM backend '{type(self.calc_real_low).__name__}' does not "
                    f"support mm_hessian_mode='{mm_hessian_mode}'."
                )
            method_kwargs = (
                {"delta": self.mm_fd_delta}
                if mm_hessian_mode == "finite_difference"
                else {}
            )

            t0 = time.perf_counter()
            H_real_np, active_atoms_from_fd = real_method(
                atoms_real_for_hess,
                info_path=info_real,
                dtype=self.H_np_dtype,
                return_partial_hessian=True,
                **method_kwargs,
            )
            local_timing["mm_hessian_real_s"] = time.perf_counter() - t0

            t0 = time.perf_counter()
            H_model_np, _ = model_method(
                atoms_model,
                info_path=info_model,
                dtype=self.H_np_dtype,
                return_partial_hessian=False,
                **method_kwargs,
            )
            local_timing["mm_hessian_model_s"] = time.perf_counter() - t0
            local_timing["mm_hessian_total_s"] = (
                float(local_timing["mm_hessian_real_s"])
                + float(local_timing["mm_hessian_model_s"])
            )
            local_timing["mm_hessian_mode"] = mm_hessian_mode

        return _MMLowOut(
            E_real=E_real_low,
            F_real=F_real_low,
            E_model=E_model_low,
            F_model=F_model_low,
            H_real=H_real_np,
            H_model=H_model_np,
            active_atoms_from_fd=active_atoms_from_fd,
            timing=local_timing,
        )

    def compute(
        self,
        coord_ang: np.ndarray,
        *,
        return_forces: bool = False,
        return_hessian: bool = False,
    ) -> Dict:
        """Run an ML/MM ONIOM single point and (optionally) forces / Hessian.

        Parameters
        ----------
        coord_ang:
            Cartesian coordinates (Å) of the full real system in the order
            defined by ``real_parm7``.
        return_forces:
            When True, include analytical forces in the returned dict.
        return_hessian:
            When True, run the 3-layer 5-pass partial Hessian assembly
            (see CHEMISTRY-RULE:8 in this module). VRAM peak is reported
            via ``timing['hess_vram_*']`` when ``print_vram`` is enabled.

        Returns
        -------
        dict
            ``{"energy": float, ...}`` plus the optional ``"forces"`` /
            ``"hessian"`` keys depending on the flags. ``"timing"`` is a
            dict of per-stage elapsed seconds (and VRAM bytes when active).
        """
        # ``freeze_atoms`` is a mutable public list on MLMMCore. Refresh the
        # derived active maps at the evaluation boundary so direct API users
        # receive the same mask semantics as the higher-level adapters.
        self._update_active_dof_mappings()
        timing: Dict[str, float | str] = {}
        hess_total_start: Optional[float] = time.perf_counter() if return_hessian else None
        hess_vram_base_alloc: Optional[float] = None
        hess_vram_base_reserved: Optional[float] = None
        hess_vram_total: Optional[float] = None
        if return_hessian and self.print_vram and self.ml_device.type == "cuda":
            torch.cuda.synchronize(device=self.ml_device)
            hess_vram_base_alloc = float(torch.cuda.memory_allocated(device=self.ml_device))
            hess_vram_base_reserved = float(torch.cuda.memory_reserved(device=self.ml_device))
            hess_vram_total = float(torch.cuda.get_device_properties(self.ml_device).total_memory)
            torch.cuda.reset_peak_memory_stats(device=self.ml_device)

        # CHEMISTRY-RULE:8 3-layer 5-pass partial Hessian assembly entry.
        atoms_real, atoms_model, atoms_model_LH, added_link_atoms, freeze_model = self._prep_3_layer_atoms(coord_ang)
        atoms_real.set_pbc(False)
        atoms_model.set_pbc(False)
        atoms_model_LH.set_pbc(False)

        need_forces = return_forces or return_hessian
        use_parallel = (self.ml_device.type == "cuda") and (getattr(self.calc_real_low, "device", None) == "cpu")
        if use_parallel:
            with ThreadPoolExecutor(max_workers=2) as executor:
                fut_ml = executor.submit(
                    self._eval_ml_high,
                    atoms_model_LH,
                    freeze_model,
                    need_forces=need_forces,
                    return_hessian=return_hessian,
                )
                fut_mm = executor.submit(
                    self._eval_mm_low,
                    atoms_real,
                    atoms_model,
                    need_forces=need_forces,
                    return_hessian=return_hessian,
                )
                ml_out = fut_ml.result()
                mm_out = fut_mm.result()
        else:
            ml_out = self._eval_ml_high(
                atoms_model_LH,
                freeze_model,
                need_forces=need_forces,
                return_hessian=return_hessian,
            )
            mm_out = self._eval_mm_low(
                atoms_real,
                atoms_model,
                need_forces=need_forces,
                return_hessian=return_hessian,
            )

        timing.update(ml_out.timing)
        timing.update(mm_out.timing)
        ml_energy = ml_out.E
        ml_forces = ml_out.F
        H_high = ml_out.H
        mm_real_energy = mm_out.E_real
        mm_real_forces = mm_out.F_real
        mm_model_energy = mm_out.E_model
        mm_model_forces = mm_out.F_model
        mm_real_hessian = mm_out.H_real
        mm_model_hessian = mm_out.H_model
        mm_active_atoms_from_fd = mm_out.active_atoms_from_fd
        del ml_out, mm_out

        # CHEMISTRY-RULE:1 Subtractive ONIOM formula. Do NOT alter sign/sum order.
        total_E = mm_real_energy + ml_energy - mm_model_energy
        results: Dict = {
            "energy": total_E,
            "energy_components": {
                "real_low": mm_real_energy,
                "model_high": ml_energy,
                "model_low": mm_model_energy,
            },
        }

        if return_forces or return_hessian:
            if (
                mm_real_forces is None
                or mm_model_forces is None
                or ml_forces is None
            ):
                raise RuntimeError(
                    "A force/Hessian evaluation returned energy-only data."
                )
            F_combined = np.copy(mm_real_forces)
            for i, ridx in enumerate(self.selection_indices):
                F_combined[ridx] += ml_forces[i] - mm_model_forces[i]

            real_to_model = self._idx_map_real_to_model
            # CHEMISTRY-RULE:2 Link-atom Hessian B-matrix (= scaled Jacobian projection). Do NOT skip param.
            for link_idx, ml_idx, mm_idx, param in added_link_atoms:
                ml_model_idx = real_to_model[ml_idx]
                r_ml = atoms_model_LH[ml_model_idx].position
                r_mm = atoms_real[mm_idx].position
                grad_link = ml_forces[link_idx]
                if self.link_atom_method == "scaled":
                    J = self._jacobian_blocks_numpy_scaled(param)
                else:
                    J = self._jacobian_blocks_numpy(r_ml, r_mm, param)
                if J is None:
                    continue
                redistributed = J @ grad_link
                F_combined[ml_idx] += redistributed[:3]
                F_combined[mm_idx] += redistributed[3:]
            results["forces"] = F_combined

        # Point-charge embedding correction (optional)
        embed_dH_info = None
        if self.embedcharge and self._embed_correction is not None:
            t0_embed = time.perf_counter()
            # ML atom symbols and coordinates
            ml_symbols = [atoms_model_LH[i].symbol for i in range(len(self._atoms_model_tpl))]
            ml_coords = np.array([atoms_model_LH[i].position for i in range(len(self._atoms_model_tpl))])
            # MM atom coordinates and charges from the real topology
            ml_set = set(self.selection_indices)
            mm_atom_indices = [i for i in range(len(atoms_real)) if i not in ml_set]
            if mm_atom_indices and self.embedcharge_cutoff is not None:
                from scipy.spatial.distance import cdist
                _ml_ref_coords = atoms_real.get_positions()[sorted(ml_set)]
                mm_coords_all = atoms_real.get_positions()[mm_atom_indices]
                dists = cdist(mm_coords_all, _ml_ref_coords).min(axis=1)
                n_before = len(mm_atom_indices)
                mask = dists <= self.embedcharge_cutoff
                mm_atom_indices = [mm_atom_indices[j] for j in range(n_before) if mask[j]]
                if self.print_timing and not getattr(self, '_embedcharge_logged', False):
                    emit(f"[embedcharge] {len(mm_atom_indices)}/{n_before} MM atoms within {self.embedcharge_cutoff:.1f} Å cutoff.",
                         narrative=True)
                    self._embedcharge_logged = True
            if mm_atom_indices:
                mm_coords = atoms_real.get_positions()[mm_atom_indices]
                # Get MM partial charges from the topology
                mm_charges = self._get_mm_charges(mm_atom_indices)

                dE_embed, dF_embed, dH_embed = self._embed_correction.compute_correction(
                    symbols=ml_symbols,
                    coords_ml_ang=ml_coords,
                    mm_coords_ang=mm_coords,
                    mm_charges=mm_charges,
                    charge=self.model_charge,
                    multiplicity=self.model_mult,
                    need_forces=return_forces or return_hessian,
                    need_hessian=return_hessian,
                )

                # Add energy correction
                results["energy"] += dE_embed

                correction_real_indices = [
                    *[int(i) for i in self.selection_indices],
                    *[int(i) for i in mm_atom_indices],
                ]

                # Scatter the complete conservative correction, including
                # movable MM point-charge sites.
                if dF_embed is not None and (return_forces or return_hessian):
                    if len(dF_embed) != len(correction_real_indices):
                        raise ValueError(
                            "Embedding force basis length does not match the "
                            "real-system index map."
                        )
                    for row, ridx in enumerate(correction_real_indices):
                        results["forces"][ridx] += dF_embed[row]

                # Store Hessian correction for later assembly
                if dH_embed is not None:
                    embed_dH_info = (dH_embed, correction_real_indices)

                timing["embedcharge_s"] = time.perf_counter() - t0_embed
                del dE_embed, dF_embed, dH_embed

        if return_hessian:
            n_real = len(atoms_real)
            n_ml = len(self.selection_indices)
            n_hess_active = self.n_hess_active

            mm_hessian_mode = getattr(
                self,
                "mm_hessian_mode",
                "finite_difference" if self.mm_fd else "analytical",
            )
            if mm_hessian_mode != "none":
                if mm_real_hessian is None or mm_model_hessian is None:
                    raise RuntimeError("MM Hessians were not computed as expected.")

                if mm_active_atoms_from_fd is not None:
                    expected = set(self.hess_active_atoms)
                    got = set(mm_active_atoms_from_fd.tolist())
                    if expected != got:
                        raise RuntimeError(
                            f"Hessian active atoms mismatch: expected {len(expected)} atoms, got {len(got)}"
                        )

                H = torch.from_numpy(mm_real_hessian).to(
                    self.ml_device, self.H_dtype
                )
                H = H.view(n_hess_active, 3, n_hess_active, 3)
                del mm_real_hessian

                H_model = torch.from_numpy(mm_model_hessian).to(
                    self.ml_device, self.H_dtype
                )
                H_model = H_model.view(n_ml, 3, n_ml, 3)
                del mm_model_hessian
            else:
                H = torch.zeros((n_hess_active, 3, n_hess_active, 3), dtype=self.H_dtype, device=self.ml_device)
                H_model = torch.zeros((n_ml, 3, n_ml, 3), dtype=self.H_dtype, device=self.ml_device)

            ml_pairs = [
                (i, self.full_to_hess_active[gi_real])
                for i, gi_real in enumerate(self.selection_indices)
                if gi_real in self.full_to_hess_active
            ]
            if ml_pairs:
                ml_sel_idx = torch.as_tensor([p[0] for p in ml_pairs], dtype=torch.long, device=self.ml_device)
                ml_active_idx = torch.as_tensor([p[1] for p in ml_pairs], dtype=torch.long, device=self.ml_device)
            else:
                ml_sel_idx = torch.empty((0,), dtype=torch.long, device=self.ml_device)
                ml_active_idx = torch.empty((0,), dtype=torch.long, device=self.ml_device)

            if H_high is not None and ml_sel_idx.numel() > 0:
                t_asm = time.perf_counter()
                H_high_mm = _gather_atom_hessian_square(
                    H_high, ml_sel_idx
                )
                H_model_mm = _gather_atom_hessian_square(
                    H_model, ml_sel_idx
                )
                H_high_mm.sub_(H_model_mm)
                H[ml_active_idx[:, None], :, ml_active_idx[None, :], :] += H_high_mm.permute(0, 2, 1, 3)
                del H_high_mm, H_model_mm
                timing["hess_asm_mlml_s"] = time.perf_counter() - t_asm
            del H_model

            real_to_model = self._idx_map_real_to_model
            link_data: List[
                Tuple[
                    int,
                    int,
                    int,
                    Optional[int],
                    Optional[int],
                    float,
                    torch.Tensor,
                ]
            ] = []
            for link_idx, ml_idx, mm_idx, param in added_link_atoms:
                ml_model_idx = real_to_model[ml_idx]
                if self.link_atom_method == "scaled":
                    K = self._jacobian_blocks_torch_scaled(param, dtype=self.H_dtype, device=self.ml_device)
                else:
                    r_ml_t = torch.tensor(atoms_model_LH[ml_model_idx].position, dtype=self.H_dtype, device=self.ml_device)
                    r_mm_t = torch.tensor(atoms_real[mm_idx].position, dtype=self.H_dtype, device=self.ml_device)
                    K = self._jacobian_blocks_torch(r_ml_t, r_mm_t, param, dtype=self.H_dtype, device=self.ml_device)
                if K is None:
                    continue
                ml_active = self.full_to_hess_active.get(ml_idx)
                mm_active = self.full_to_hess_active.get(mm_idx)
                if ml_active is None and mm_active is None:
                    continue
                link_data.append((link_idx, ml_idx, mm_idx, ml_active, mm_active, param, K))

            def _add_endpoint_blocks(
                block6: torch.Tensor,
                left_ml: Optional[int],
                left_mm: Optional[int],
                right_ml: Optional[int],
                right_mm: Optional[int],
            ) -> None:
                left = ((left_ml, slice(0, 3)), (left_mm, slice(3, 6)))
                right = ((right_ml, slice(0, 3)), (right_mm, slice(3, 6)))
                for left_active, left_slice in left:
                    if left_active is None:
                        continue
                    for right_active, right_slice in right:
                        if right_active is None:
                            continue
                        H[left_active, :, right_active, :].add_(
                            block6[left_slice, right_slice]
                        )

            F_high_t = torch.as_tensor(
                ml_forces, dtype=self.H_dtype, device=self.ml_device
            )
            has_link_force = bool((F_high_t.abs() > 1e-12).any().item())
            if link_data and (H_high is not None or has_link_force):
                t_asm = time.perf_counter()
                I3 = torch.eye(3, dtype=self.H_dtype, device=self.ml_device)
                for link_idx, ml_idx, mm_idx, ml_active, mm_active, param, K in link_data:

                    if H_high is not None:
                        H_l = H_high[link_idx, :, link_idx, :]
                        H_self = K.T @ H_l @ K
                        _add_endpoint_blocks(
                            H_self,
                            ml_active,
                            mm_active,
                            ml_active,
                            mm_active,
                        )

                    # B-matrix constraint correction: only needed for fixed-distance
                    # link atoms. For scaled (g-factor) link atoms, d²L/dQ² = 0
                    # (position is linear in QM1 and MM), so no correction needed.
                    if self.link_atom_method != "scaled":
                        f_L = -F_high_t[link_idx]
                        dist = param

                        r_ml_t = torch.as_tensor(atoms_model_LH[real_to_model[ml_idx]].position,
                                                 dtype=self.H_dtype, device=self.ml_device)
                        r_mm_t = torch.as_tensor(atoms_real[mm_idx].position,
                                                 dtype=self.H_dtype, device=self.ml_device)
                        v = r_mm_t - r_ml_t
                        R_sq = torch.dot(v, v)
                        inv_R = torch.rsqrt(torch.clamp(R_sq, min=1.0e-24))
                        inv_R2 = inv_R * inv_R
                        u = v * inv_R

                        alpha = torch.dot(u, f_L)
                        uuT = torch.outer(u, u)
                        ufT = torch.outer(u, f_L)
                        fTu = torch.outer(f_L, u)
                        B = (alpha * (3.0 * uuT - I3) - (ufT + fTu)) * inv_R2

                        H_corr6 = torch.zeros((6, 6), dtype=self.H_dtype, device=self.ml_device)
                        H_corr6[0:3, 0:3] = B
                        H_corr6[3:6, 3:6] = B
                        H_corr6[0:3, 3:6] = -B
                        H_corr6[3:6, 0:3] = -B
                        H_corr6.mul_(dist)

                        _add_endpoint_blocks(
                            H_corr6,
                            ml_active,
                            mm_active,
                            ml_active,
                            mm_active,
                        )
                timing["hess_asm_link_self_s"] = time.perf_counter() - t_asm
            del F_high_t, ml_forces

            if H_high is not None and link_data and ml_sel_idx.numel() > 0:
                t_asm = time.perf_counter()
                for link_idx, _ml_idx, _mm_idx, ml_active, mm_active, _param, K in link_data:
                    H_coup = H_high[link_idx].index_select(1, ml_sel_idx).permute(1, 0, 2).contiguous()  # (K,3,3)
                    H_row = torch.einsum("ac,bcd->bad", K.T, H_coup)  # (K,6,3)
                    H_col = torch.einsum("bca,cd->bad", H_coup, K)    # (K,3,6)

                    # Mixed scalar/tensor indexing in PyTorch returns (3, K, 3) for
                    # H[scalar, :, tensor, :], so align H_row blocks explicitly.
                    # Write back by ASSIGNMENT (as the ML-ML block above does), never `.add_()`:
                    # mixing a scalar and a tensor index is advanced indexing, which returns a
                    # COPY, so an in-place add on it silently discards the coupling block.
                    if ml_active is not None:
                        H[ml_active, :, ml_active_idx, :] = (
                            H[ml_active, :, ml_active_idx, :]
                            + H_row[:, 0:3, :].permute(1, 0, 2)
                        )
                        H[ml_active_idx, :, ml_active, :] = (
                            H[ml_active_idx, :, ml_active, :] + H_col[:, :, 0:3]
                        )
                    if mm_active is not None:
                        H[mm_active, :, ml_active_idx, :] = (
                            H[mm_active, :, ml_active_idx, :]
                            + H_row[:, 3:6, :].permute(1, 0, 2)
                        )
                        H[ml_active_idx, :, mm_active, :] = (
                            H[ml_active_idx, :, mm_active, :] + H_col[:, :, 3:6]
                        )
                timing["hess_asm_link_ml_s"] = time.perf_counter() - t_asm

            if H_high is not None and link_data:
                t_asm = time.perf_counter()
                n_links = len(link_data)
                for a in range(n_links):
                    link_idx_a, _ml_a, _mm_a, ml_a_active, mm_a_active, _param_a, K_a = link_data[a]

                    for b in range(a + 1, n_links):
                        link_idx_b, _ml_b, _mm_b, ml_b_active, mm_b_active, _param_b, K_b = link_data[b]

                        H_ab = H_high[link_idx_a, :, link_idx_b, :]
                        HAB = K_a.T @ H_ab @ K_b
                        _add_endpoint_blocks(
                            HAB,
                            ml_a_active,
                            mm_a_active,
                            ml_b_active,
                            mm_b_active,
                        )
                        _add_endpoint_blocks(
                            HAB.T,
                            ml_b_active,
                            mm_b_active,
                            ml_a_active,
                            mm_a_active,
                        )
                timing["hess_asm_link_link_s"] = time.perf_counter() - t_asm
            del H_high
            if self.ml_device.type == "cuda":
                torch.cuda.empty_cache()

            # Add point-charge embedding Hessian correction
            if embed_dH_info is not None:
                t_asm = time.perf_counter()
                embed_dH, correction_real_indices = embed_dH_info
                n_correction_atoms = len(correction_real_indices)
                correction_positions: List[int] = []
                hessian_active_positions: List[int] = []
                for position, real_idx in enumerate(correction_real_indices):
                    active_position = self.full_to_hess_active.get(real_idx)
                    if active_position is None:
                        continue
                    correction_positions.append(position)
                    hessian_active_positions.append(active_position)
                if correction_positions:
                    active_idx = torch.as_tensor(
                        hessian_active_positions,
                        dtype=torch.long,
                        device=self.ml_device,
                    )
                    correction_dofs = np.asarray(
                        [
                            3 * position + axis
                            for position in correction_positions
                            for axis in range(3)
                        ],
                        dtype=np.intp,
                    )
                    embed_matrix = np.asarray(embed_dH).reshape(
                        3 * n_correction_atoms,
                        3 * n_correction_atoms,
                    )
                    if (
                        correction_dofs.size == embed_matrix.shape[0]
                        and np.array_equal(
                            correction_dofs,
                            np.arange(embed_matrix.shape[0]),
                        )
                    ):
                        dH_sub_cpu = embed_matrix
                    else:
                        dH_sub_cpu = active_square(
                            embed_matrix, correction_dofs
                        )
                    dH_sub = torch.as_tensor(
                        dH_sub_cpu,
                        dtype=self.H_dtype,
                        device=self.ml_device,
                    ).reshape(
                        len(correction_positions),
                        3,
                        len(correction_positions),
                        3,
                    )
                    H[
                        active_idx[:, None],
                        :,
                        active_idx[None, :],
                        :,
                    ] += dH_sub.permute(0, 2, 1, 3)
                    del (
                        active_idx,
                        correction_dofs,
                        dH_sub,
                        dH_sub_cpu,
                        embed_matrix,
                    )
                embed_dH_info = None
                del embed_dH
                if self.ml_device.type == "cuda":
                    torch.cuda.empty_cache()
                timing["hess_asm_embed_s"] = time.perf_counter() - t_asm

            if self.symmetrize_hessian:
                t_asm = time.perf_counter()
                # Bounded-peak symmetrization via shared helper (writes both
                # triangles; peak temp <= chunk^2 instead of full N×N clone).
                # H_flat is a view of H's storage — symmetrize_inplace mutates
                # in place, so the 4D H view automatically reflects the result
                # without rebinding.
                from mlmm.core.utils import symmetrize_inplace
                H_flat = H.view(3 * n_hess_active, 3 * n_hess_active)
                symmetrize_inplace(H_flat)
                timing["hess_asm_sym_s"] = time.perf_counter() - t_asm

            if self.return_partial_hessian:
                results["hessian"] = H.detach()
                results["within_partial_hessian"] = self._build_within_partial_hessian()
            else:
                t_asm = time.perf_counter()
                H_full = torch.zeros((n_real, 3, n_real, 3), dtype=self.H_dtype, device=self.ml_device)
                active_idx = torch.as_tensor(self.hess_active_atoms, dtype=torch.long, device=self.ml_device)
                if active_idx.numel() > 0:
                    # Scatter-assign accepts a strided source, avoiding the
                    # additional full-size copy created by `.contiguous()`.
                    H_full[active_idx[:, None], :, active_idx[None, :], :] = H.permute(0, 2, 1, 3)
                results["hessian"] = H_full.detach()
                timing["hess_asm_full_expand_s"] = time.perf_counter() - t_asm
                del H_full

            if hess_total_start is not None:
                timing["hessian_total_s"] = time.perf_counter() - hess_total_start
                results["timing"] = timing
                hessian_vram_summary = ""
                hessian_vram_detail = None
                if self.print_vram and self.ml_device.type == "cuda":
                    torch.cuda.synchronize(device=self.ml_device)
                    base_alloc = float(hess_vram_base_alloc or 0.0)
                    base_reserved = float(hess_vram_base_reserved or 0.0)
                    peak_alloc = max(
                        float(torch.cuda.max_memory_allocated(device=self.ml_device)) - base_alloc,
                        0.0,
                    ) / 1e9
                    peak_reserved_abs = float(torch.cuda.max_memory_reserved(device=self.ml_device))
                    peak_reserved = max(
                        peak_reserved_abs - base_reserved,
                        0.0,
                    ) / 1e9
                    total_vram = float(hess_vram_total or torch.cuda.get_device_properties(self.ml_device).total_memory) / 1e9
                    remaining_vram = max((total_vram * 1e9) - peak_reserved_abs, 0.0) / 1e9
                    hessian_vram_summary = f", peak VRAM {peak_alloc:.2f} GB"
                    hessian_vram_detail = (
                        f"[HessianVRAM] total={total_vram:.3f} GB | "
                        f"peak_allocated={peak_alloc:.3f} GB | "
                        f"peak_reserved={peak_reserved:.3f} GB | "
                        f"remaining={remaining_vram:.3f} GB"
                    )
                from mlmm.core.utils import verbose_level
                if self.print_timing:
                    ml_mode = timing.get("ml_hessian_mode")
                    ml_time = timing.get("ml_hessian_s")
                    mode_label = str(ml_mode) if ml_mode is not None else "ML/MM"
                    emit(
                        f"[hessian] Completed {mode_label} Hessian: "
                        f"{timing['hessian_total_s']:.2f} s{hessian_vram_summary}",
                        detail=True,
                    )
                    if verbose_level() >= 3:
                        if ml_mode is not None and ml_time is not None:
                            click.echo(f"[HessianTiming] ML Hessian ({ml_mode}): {ml_time:.2f} s")
                        if "mm_hessian_total_s" in timing:
                            click.echo(
                                "[HessianTiming] MM Hessian "
                                f"({timing.get('mm_hessian_mode', 'unknown')}): "
                                f"REAL {timing['mm_hessian_real_s']:.2f} s | "
                                f"MODEL {timing['mm_hessian_model_s']:.2f} s | "
                                f"total {timing['mm_hessian_total_s']:.2f} s"
                            )
                        asm_parts = []
                        for key, label in (
                            ("hess_asm_mlml_s", "ML-ML"),
                            ("hess_asm_link_self_s", "link-self"),
                            ("hess_asm_link_ml_s", "link-ML"),
                            ("hess_asm_link_link_s", "link-link"),
                            ("hess_asm_sym_s", "sym"),
                            ("hess_asm_full_expand_s", "full-expand"),
                        ):
                            if key in timing:
                                asm_parts.append(f"{label} {float(timing[key]):.2f} s")
                        if asm_parts:
                            click.echo(f"[HessianTiming] Assembly: {' | '.join(asm_parts)}")
                        click.echo(f"[HessianTiming] Hessian total: {timing['hessian_total_s']:.2f} s")
                if hessian_vram_detail is not None and verbose_level() >= 3:
                    click.echo(hessian_vram_detail)

            del H
            if self.ml_device.type == "cuda":
                torch.cuda.empty_cache()

        return self._finalize_result_constraints(results)


#                ASE Calculator wrapper for ML/MM (ONIOM)

class MLMMASECalculator(Calculator):
    """ASE Calculator wrapping MLMMCore for use with DMF and other ASE-based methods.

    The underlying MLMMCore takes full-system coordinates (Angstrom) and
    returns energy in eV and forces in eV/Angstrom, which matches ASE conventions.
    """

    implemented_properties = ["energy", "forces"]

    def __init__(self, core: "MLMMCore", **kwargs):
        super().__init__(**kwargs)
        self.core = core

    def calculate(self, atoms=None, properties=("energy",), system_changes=all_changes):
        super().calculate(atoms, properties, system_changes)
        coord_ang = atoms.get_positions().astype(float)
        want_forces = "forces" in properties
        res = self.core.compute(coord_ang, return_forces=want_forces, return_hessian=False)
        self.results = {
            "energy": float(res["energy"]),
        }
        if want_forces:
            self.results["forces"] = res["forces"].reshape(-1, 3)


#                     PySisyphus Calculator (ML/MM)

from pysisyphus.calculators.Calculator import Calculator as PySiCalc


class mlmm(PySiCalc):
    """Pysisyphus ``Calculator`` adapter for :class:`MLMMCore`.

    Wraps the lower-level MLMMCore engine so existing pysisyphus
    optimizers / IRC / NEB drivers can consume ML/MM ONIOM energies +
    forces + Hessian through the standard ``Calculator`` interface
    (``get_energy`` / ``get_forces`` / ``get_hessian``). All constructor
    arguments map 1:1 onto the underlying MLMMCore; see that class for
    the per-argument contract.
    """

    implemented_properties = ["energy", "forces", "hessian"]

    def __init__(
        self,
        input_pdb: Optional[str] = None,
        real_parm7: Optional[str] = None,
        model_pdb: Optional[str] = None,
        *,
        coordinate_path: Optional[str] = None,
        model_charge: int = 0,
        model_mult: int = 1,
        link_mlmm: List[Tuple[str, str]] | None = None,
        link_atom_method: str = "scaled",
        # ML backend selection
        backend: str = "uma",
        uma_model: str = DEFAULT_UMA_MODEL,
        uma_task_name: str = "omol",
        uma_precision: str = "fp32",
        workers: int = 1,
        workers_per_node: int = 1,
        orb_model: str = "orb_v3_conservative_omol",
        orb_precision: str = "float64",
        mace_model: str = "MACE-OMOL-0",
        mace_dtype: str = "float64",
        aimnet2_model: str = "aimnet2",
        # Custom ML backend from a user Python file (--calc-file)
        calc_file: Optional[str] = None,
        calc_factory: str = "get_calculator",
        # MM settings
        mm_fd: bool = True,
        mm_hessian_mode: Optional[str] = None,
        mm_fd_dir: Optional[str] = None,
        mm_fd_delta: float = 1e-3,
        symmetrize_hessian: bool = True,
        out_hess_torch: bool = True,
        H_double: bool = True,
        hessian_calc_mode: str = "FiniteDifference",
        ml_device: str = "auto",
        ml_cuda_idx: int = 0,
        mm_device: str = "cpu",
        mm_cuda_idx: int = 0,
        mm_threads: int = 16,
        mm_backend: str = "hessian_ff",
        freeze_atoms: List[int] | None = None,
        return_partial_hessian: bool = True,
        print_timing: bool = True,
        print_vram: bool = True,
        hess_cutoff: Optional[float] = None,
        movable_cutoff: Optional[float] = None,
        use_bfactor_layers: bool = True,       # matches MLMM_CALC_KW default
        hess_mm_atoms: Optional[List[int]] = None,
        movable_mm_atoms: Optional[List[int]] = None,
        frozen_mm_atoms: Optional[List[int]] = None,
        # Point-charge embedding correction
        embedcharge: bool = False,
        embedcharge_step: float = 1.0e-3,
        embedcharge_cutoff: float = 12.0,
        xtb_cmd: str = "xtb",
        xtb_acc: float = 0.2,
        xtb_workdir: str = "tmp",
        xtb_keep_files: bool = False,
        xtb_ncores: int = 4,
        use_cmap: bool = True,
        _high_level_backend: Optional[Any] = None,
        _skip_high_level_backend: bool = False,
        **kwargs,
    ):
        # --- v0.1.x backward compatibility aliases ---
        if "real_pdb" in kwargs:
            warnings.warn("'real_pdb' is deprecated; use 'input_pdb'.", DeprecationWarning, stacklevel=2)
            if input_pdb is None:
                input_pdb = kwargs.pop("real_pdb")
            else:
                kwargs.pop("real_pdb")
        for _old_name in ("real_rst7", "vib_run", "vib_dir"):
            if _old_name in kwargs:
                warnings.warn(f"'{_old_name}' is no longer used and will be ignored.", DeprecationWarning, stacklevel=2)
                kwargs.pop(_old_name)

        self._freeze_atoms = [] if freeze_atoms is None else list(freeze_atoms)
        super().__init__(charge=model_charge, mult=model_mult, **kwargs)

        self.core = MLMMCore(
            input_pdb=input_pdb,
            coordinate_path=coordinate_path,
            real_parm7=real_parm7,
            model_pdb=model_pdb,
            model_charge=model_charge,
            model_mult=model_mult,
            link_mlmm=link_mlmm,
            link_atom_method=link_atom_method,
            backend=backend,
            uma_model=uma_model,
            uma_task_name=uma_task_name,
            uma_precision=uma_precision,
            workers=workers,
            workers_per_node=workers_per_node,
            orb_model=orb_model,
            orb_precision=orb_precision,
            mace_model=mace_model,
            mace_dtype=mace_dtype,
            aimnet2_model=aimnet2_model,
            calc_file=calc_file,
            calc_factory=calc_factory,
            mm_fd=mm_fd,
            mm_hessian_mode=mm_hessian_mode,
            mm_fd_dir=mm_fd_dir,
            mm_fd_delta=mm_fd_delta,
            symmetrize_hessian=symmetrize_hessian,
            H_double=H_double,
            ml_device=ml_device,
            ml_cuda_idx=ml_cuda_idx,
            mm_device=mm_device,
            mm_cuda_idx=mm_cuda_idx,
            mm_threads=mm_threads,
            mm_backend=mm_backend,
            freeze_atoms=self._freeze_atoms,
            hessian_calc_mode=hessian_calc_mode,
            return_partial_hessian=return_partial_hessian,
            print_timing=print_timing,
            print_vram=print_vram,
            hess_cutoff=hess_cutoff,
            movable_cutoff=movable_cutoff,
            use_bfactor_layers=use_bfactor_layers,
            hess_mm_atoms=hess_mm_atoms,
            movable_mm_atoms=movable_mm_atoms,
            frozen_mm_atoms=frozen_mm_atoms,
            embedcharge=embedcharge,
            embedcharge_step=embedcharge_step,
            embedcharge_cutoff=embedcharge_cutoff,
            xtb_cmd=xtb_cmd,
            xtb_acc=xtb_acc,
            xtb_workdir=xtb_workdir,
            xtb_keep_files=xtb_keep_files,
            xtb_ncores=xtb_ncores,
            use_cmap=use_cmap,
            _high_level_backend=_high_level_backend,
            _skip_high_level_backend=_skip_high_level_backend,
        )

        self.out_hess_torch = bool(out_hess_torch)
        self.hess_torch_double = bool(H_double)
        self._hess_scale = EV2AU / ANG2BOHR / ANG2BOHR

    @property
    def freeze_atoms(self) -> List[int] | None:
        return self.core.freeze_atoms

    @freeze_atoms.setter
    def freeze_atoms(self, indices: List[int] | None):
        requested = set([] if indices is None else map(int, indices))
        requested.update(map(int, getattr(self.core, "frozen_layer_indices", [])))
        self._freeze_atoms = sorted(requested)
        self.core.freeze_atoms = self._freeze_atoms
        self.core._update_active_dof_mappings()

    def _run_core(self, coords, *, want_forces: bool, want_hessian: bool):
        coord_ang = np.asarray(coords).reshape(-1, 3) * BOHR2ANG
        res = self.core.compute(coord_ang, return_forces=want_forces or want_hessian, return_hessian=want_hessian)
        out = {"energy": res["energy"] * EV2AU}
        if want_forces or want_hessian:
            out["forces"] = (res["forces"] * (EV2AU / ANG2BOHR)).flatten()
        if want_hessian:
            H = res.pop("hessian")
            H = H.view(H.size(0) * 3, H.size(2) * 3)
            H.mul_(self._hess_scale)
            if self.out_hess_torch:
                target_dtype = torch.float64 if self.hess_torch_double else torch.float32
                out["hessian"] = H.to(target_dtype).detach().requires_grad_(False)
            else:
                out["hessian"] = H.detach().cpu().numpy()
                # Release GPU storage after the required NumPy conversion.
                del H
                if torch.cuda.is_available():
                    torch.cuda.empty_cache()
            if "within_partial_hessian" in res:
                out["within_partial_hessian"] = res["within_partial_hessian"]
        return out

    def get_energy(self, elem, coords):
        return self._run_core(coords, want_forces=False, want_hessian=False)

    def get_forces(self, elem, coords):
        return self._run_core(coords, want_forces=True, want_hessian=False)

    def get_hessian(self, elem, coords):
        return self._run_core(coords, want_forces=True, want_hessian=True)


#                     PySisyphus Calculator (MM-only)


class mlmm_mm_only(PySiCalc):
    """PySisyphus calculator that returns MM-only energy and forces (F_real_mm).

    Used for microiteration: relaxes the MM region without ML computation.
    Shares the MLMMCore from an existing ``mlmm`` calculator to avoid
    re-initializing topology and force field objects.
    """

    implemented_properties = ["energy", "forces"]

    def __init__(self, core: "MLMMCore", *, freeze_atoms: list[int] | None = None, **kwargs):
        super().__init__(charge=core.model_charge, mult=core.model_mult, **kwargs)
        self.core = core
        self._freeze_atoms = list(freeze_atoms) if freeze_atoms else []

    def _run_core(self, coords, *, want_forces: bool):
        coord_ang = np.asarray(coords).reshape(-1, 3) * BOHR2ANG
        atoms_real = self.core._atoms_real_tpl.copy()
        atoms_real.set_positions(coord_ang)
        atoms_real.set_pbc(False)
        atoms_real.calc = self.core.calc_real_low
        E_real = float(atoms_real.get_potential_energy())
        out = {"energy": E_real * EV2AU}
        if want_forces:
            F_real = np.double(atoms_real.get_forces())
            # Zero forces on frozen atoms
            for i in self._freeze_atoms:
                if 0 <= i < F_real.shape[0]:
                    F_real[i, :] = 0.0
            out["forces"] = (F_real * (EV2AU / ANG2BOHR)).flatten()
        return out

    def get_energy(self, elem, coords):
        return self._run_core(coords, want_forces=False)

    def get_forces(self, elem, coords):
        return self._run_core(coords, want_forces=True)

    def get_hessian(self, elem, coords):
        raise NotImplementedError("MM-only calculator does not support Hessian computation.")


#               v0.1.x compatibility: mlmm_ase() factory


def mlmm_ase(**kwargs):
    """v0.1.x compatibility wrapper.

    Accepts all MLMMCore parameters as keyword arguments and returns
    an MLMMASECalculator.  Equivalent to::

        MLMMASECalculator(MLMMCore(**kwargs))
    """
    warnings.warn(
        "mlmm_ase() is deprecated; use MLMMASECalculator(MLMMCore(...)) instead.",
        DeprecationWarning,
        stacklevel=2,
    )
    core = MLMMCore(**kwargs)
    return MLMMASECalculator(core)
