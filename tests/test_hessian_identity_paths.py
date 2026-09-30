"""Hessian reuse accepts aliases of one input while keeping distinct inputs apart."""

import importlib.util
from pathlib import Path

import numpy as np
import pytest

from mlmm.core.utils import prepare_input_structure
from mlmm.io import hessian_cache


def _reuses_hessian(first_config, second_config):
    def identity(config):
        return hessian_cache.build_identity(
            atoms=["C"], cart_coords=np.zeros(3), run_id="path-reuse",
            potential=hessian_cache._potential_identity(config),
        )
    hessian_cache.clear()
    try:
        hessian_cache.store("ts", np.eye(3), identity=identity(first_config))
        return hessian_cache.load_matching("ts", identity(second_config)) is not None
    finally:
        hessian_cache.clear()


@pytest.mark.parametrize("key", ["input_pdb", "real_parm7", "model_pdb", "calc_file"])
@pytest.mark.parametrize("alias_kind", ["relative", "symlink"])
def test_input_path_aliases_reuse_the_hessian(tmp_path, monkeypatch, key, alias_kind):
    source = tmp_path / "input.pdb"
    source.write_text("END\n")
    alias = tmp_path / "alias.pdb"
    alias.symlink_to(source)
    monkeypatch.chdir(tmp_path)
    other_path = "input.pdb" if alias_kind == "relative" else str(alias)
    assert _reuses_hessian({key: str(source)}, {key: other_path})


@pytest.mark.parametrize("source_format", ["cif", "overflow-pdb"])
def test_normalized_structures_reuse_the_source_hessian(tmp_path, source_format):
    source = tmp_path / ("large-id.cif" if source_format == "cif" else "large.pdb")
    if source_format == "cif":
        spec = importlib.util.spec_from_file_location(
            "structure_fixture", Path(__file__).with_name("test_structure_formats.py")
        )
        fixture = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(fixture)
        fixture._write_minimal_cif(source)
    else:
        source.write_text(
            "ATOM      1  C   MOL A   1       0.000   0.000   0.000  1.00  0.00           C\n"
            * 100000 + "END\n"
        )
    first = prepare_input_structure(source)
    second = prepare_input_structure(source)
    try:
        assert first.source_path != second.source_path
        assert _reuses_hessian(
            {"input_pdb": str(first.source_path)}, {"input_pdb": str(second.source_path)},
        )
    finally:
        first.cleanup()
        second.cleanup()


def test_distinct_topology_files_still_reject_reuse(tmp_path):
    first = tmp_path / "first.pdb"
    second = tmp_path / "second.pdb"
    first.write_text("END\n")
    second.write_text("END\n")
    assert not _reuses_hessian({"input_pdb": str(first)}, {"input_pdb": str(second)})
