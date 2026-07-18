"""Falsifiers for M40: the pure normal-mode kernel is lowered below the product.

The mass-weighting / rigid-projection / diagonalization / mode-conversion kernel
now lives in ``pysisyphus.normal_modes`` (a sibling of
``pysisyphus.tr_projection``). The bundled engine must not import the product
namespaces, and ``mlmm.workflows.freq`` must keep byte-identical
compatibility re-exports so existing callers are unaffected.
"""

from __future__ import annotations

import ast
from pathlib import Path

import numpy as np
import pytest
import torch

import pysisyphus
import pysisyphus.normal_modes as nm
from mlmm.workflows import freq as product_freq
from mlmm.core import utils as product_utils


PRODUCT_NAMESPACES = ("pdb2reaction", "mlmm")

# The five pure functions lowered in M40 plus the bounded symmetrizer.
LOWERED_KERNEL_NAMES = (
    "_safe_masses_amu",
    "_mw_projected_hessian",
    "_mass_weighted_hessian",
    "_frequencies_cm_and_modes",
    "_mw_mode_to_cart",
)


def _iter_pysisyphus_py_files():
    root = Path(pysisyphus.__file__).resolve().parent
    for path in sorted(root.rglob("*.py")):
        if "__pycache__" in path.parts:
            continue
        yield path


def _imports_a_product_namespace(tree: ast.AST) -> list[str]:
    """Return the offending fully-qualified module names, if any."""
    offenders: list[str] = []
    for node in ast.walk(tree):
        if isinstance(node, ast.Import):
            for alias in node.names:
                top = alias.name.split(".", 1)[0]
                if top in PRODUCT_NAMESPACES:
                    offenders.append(alias.name)
        elif isinstance(node, ast.ImportFrom):
            # level > 0 is a package-relative import inside pysisyphus; skip.
            if node.level == 0 and node.module is not None:
                top = node.module.split(".", 1)[0]
                if top in PRODUCT_NAMESPACES:
                    offenders.append(node.module)
    return offenders


def test_no_product_imports_anywhere_in_pysisyphus():
    """AST falsifier: no file under pysisyphus/** imports a product namespace.

    This is the headline invariant of M40 (and the base case C12 must preserve):
    the bundled engine, including ``TSHessianOptimizer``, is free of any upward
    ``pdb2reaction``/``mlmm`` import.
    """
    violations: dict[str, list[str]] = {}
    scanned = 0
    for path in _iter_pysisyphus_py_files():
        scanned += 1
        tree = ast.parse(path.read_text(encoding="utf-8"), filename=str(path))
        offenders = _imports_a_product_namespace(tree)
        if offenders:
            violations[str(path)] = offenders
    assert scanned > 0, "AST scan found no pysisyphus source files"
    assert not violations, (
        "Bundled-engine files import a product namespace (upward dependency): "
        f"{violations}"
    )


def test_tshessianoptimizer_imports_the_lower_kernel_directly():
    """TSHessianOptimizer sources the kernel from pysisyphus.normal_modes."""
    src = Path(
        pysisyphus.__file__
    ).resolve().parent.joinpath("tsoptimizers", "TSHessianOptimizer.py").read_text(
        encoding="utf-8"
    )
    assert "from pysisyphus.normal_modes import" in src
    assert "from mlmm.workflows.freq import" not in src
    assert "from mlmm" not in src and "import mlmm" not in src


def test_product_freq_reexports_are_the_lower_objects():
    """The old workflow symbols ARE the lowered implementations (same object)."""
    for name in LOWERED_KERNEL_NAMES:
        lower = getattr(nm, name)
        reexport = getattr(product_freq, name)
        assert reexport is lower, f"{name}: re-export is not the lowered object"


def test_symmetrizer_reexport_is_the_lower_object():
    assert product_utils.symmetrize_inplace is nm.symmetrize_inplace


COORDS = np.array(
    [
        [0.0, 0.0, 0.0],
        [1.2, 0.0, 0.0],
        [0.1, 1.1, 0.0],
        [0.2, 0.3, 1.3],
        [1.0, 0.8, 0.7],
    ],
    dtype=np.float64,
)
ATOMIC_NUMBERS = [6, 1, 8, 7, 16]


def _random_spd(n: int, seed: int) -> torch.Tensor:
    gen = torch.Generator().manual_seed(seed)
    raw = torch.randn((n, n), generator=gen, dtype=torch.float64)
    return raw.T @ raw + 0.2 * torch.eye(n, dtype=torch.float64)


@pytest.mark.parametrize("tr_projection", ["constrained", "legacy-active"])
def test_full_and_partial_phva_identical_through_both_surfaces(tr_projection):
    """Controlled full/partial PHVA vectors are identical via freq and normal_modes.

    The product ``workflows.freq`` surface and the lower ``normal_modes`` surface
    must return bit-identical frequencies (they are the same object, so this also
    guards against a future divergent shim).
    """
    device = torch.device("cpu")
    full_hessian = _random_spd(15, seed=41)
    active_atoms = [2, 3, 4]
    active_dofs = [3 * a + axis for a in active_atoms for axis in range(3)]
    active_hessian = full_hessian[active_dofs][:, active_dofs].clone()
    freeze_idx = [0, 1]

    results = {}
    for surface_name, surface in (("freq", product_freq), ("normal_modes", nm)):
        full_info: dict = {}
        full_freqs, _ = surface._frequencies_cm_and_modes(
            full_hessian.clone(),
            ATOMIC_NUMBERS,
            COORDS,
            device,
            freeze_idx=freeze_idx,
            tr_projection=tr_projection,
            projection_info=full_info,
        )
        active_info: dict = {}
        active_freqs, _ = surface._frequencies_cm_and_modes(
            active_hessian.clone(),
            ATOMIC_NUMBERS,
            COORDS,
            device,
            freeze_idx=freeze_idx,
            tr_projection=tr_projection,
            projection_info=active_info,
        )
        # Full-Hessian and active-block PHVA agree within one surface.
        np.testing.assert_allclose(full_freqs, active_freqs, rtol=1e-10, atol=1e-8)
        results[surface_name] = (full_freqs, active_freqs, full_info, active_info)

    # The two import surfaces produce identical numbers and provenance.
    np.testing.assert_array_equal(results["freq"][0], results["normal_modes"][0])
    np.testing.assert_array_equal(results["freq"][1], results["normal_modes"][1])
    assert results["freq"][2]["treatment"] == results["normal_modes"][2]["treatment"]
    assert (
        results["freq"][2]["effective_rank"]
        == results["normal_modes"][2]["effective_rank"]
    )


def test_lowered_module_imports_no_product_namespace():
    """The lowered module's own AST is free of product imports (self-check)."""
    tree = ast.parse(
        Path(nm.__file__).read_text(encoding="utf-8"), filename=nm.__file__
    )
    assert _imports_a_product_namespace(tree) == []
