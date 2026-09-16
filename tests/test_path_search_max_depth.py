"""Verify subdivision limits and numerical MEP outcomes.

Depth N terminates subdivision at cap N; cap zero runs one input interval.
Chemical diagnostic failures must not overwrite the solver's convergence fact.
"""

from __future__ import annotations


import click
import pytest

from mlmm.core.defaults import SEARCH_KW
from mlmm.workflows import path_search


def test_max_depth_option_is_declared_on_both_entry_points() -> None:
    from mlmm.workflows.all import cli as all_cli

    for command in (path_search.cli, all_cli):
        options = [p for p in command.params if "--max-depth" in getattr(p, "opts", ())]
        assert len(options) == 1, command.name
        option = options[0]
        # A declared default of None keeps `cli_param_overridden` able to tell an
        # explicit value from an omission, so YAML stays the middle layer.
        assert option.default is None
        assert option.show_default == str(SEARCH_KW["max_depth"])
        assert isinstance(option.type, click.IntRange)
        assert option.type.min == 0

@pytest.mark.parametrize("mep_mode", ["gsm", "dmf"])
@pytest.mark.parametrize(
    "max_depth,depth,primary_hei,refined_hei,n_segments,single_opt_executed",
    [
        (0, 0, 1, 1, 1, False),
        (0, 0, 0, 1, 0, False),
        (0, 0, 2, 1, 0, False),
        (0, 1, 1, 1, 1, False),
        (2, 1, 1, 1, 1, True),
        (2, 2, 1, 1, 1, False),
        (2, 3, 1, 1, 1, False),
        (1, 0, 0, 1, 0, False),
        (1, 0, 1, 0, 0, True),
        (1, 0, 1, 2, 0, True),
        (1, 0, 1, 1, 1, True),
    ],
)
def test_single_optimizer_provenance_and_refined_hei_boundary(
    tmp_path, monkeypatch, max_depth, depth, primary_hei, refined_hei,
    n_segments, single_opt_executed, mep_mode,
):
    from types import SimpleNamespace

    def mep(hei):
        energies = [0.0, 0.0, 0.0]
        energies[hei] = 1.0
        return SimpleNamespace(
            images=[SimpleNamespace(energy=energy) for energy in energies],
            energies=energies, hei_idx=hei, is_converged=True,
        )

    primary, refined = mep(primary_hei), mep(refined_hei)
    optimizer_calls = []

    def optimize(geometry, *_args, **kwargs):
        optimizer_calls.append(kwargs["tag"])
        return geometry, True

    # A reactive refined interval with no further changes on either side.
    changes = iter([True, True, False, False])
    monkeypatch.setattr(path_search, "_run_mep_between", lambda *_a, **_k: primary)
    if hasattr(path_search, "_run_dmf_between"):
        monkeypatch.setattr(path_search, "_run_dmf_between", lambda *_a, **_k: primary)
    monkeypatch.setattr(path_search, "_refine_between", lambda *_a, **_k: refined)
    monkeypatch.setattr(path_search, "_optimize_single", optimize)
    monkeypatch.setattr(
        path_search, "_has_bond_change",
        lambda *_a, **_k: (next(changes), "Bond formed"),
    )
    monkeypatch.setattr(path_search, "_stitch_paths", lambda parts, **_k: parts[0])

    result = path_search._build_multistep_path(
        primary.images[0], primary.images[-1], None,
        geom_cfg={}, gs_cfg={}, stopt_cfg={}, single_opt_cfg={}, bond_cfg={},
        search_cfg={"max_depth": max_depth, "stitch_rmsd_thresh": 1e-4,
                    "bridge_rmsd_thresh": 1e-4},
        refine_mode_kind="peak", mep_mode_kind=mep_mode,
        out_dir=tmp_path, ref_pdb_path=None, depth=depth, seg_counter=[0],
        branch_tag="pair_00",
    )

    assert len(optimizer_calls) == (2 if single_opt_executed else 0)
    assert result.single_opt_executed is single_opt_executed
    assert len(result.segments) == n_segments
    if max_depth == 0 and primary_hei in (0, 2):
        assert len(result.required_outcomes) == 1
        assert len(result.images) == 3
    if result.segments:
        assert result.segments[0].tag.endswith("_maxdepth") is (max_depth > 0 and depth >= max_depth)

    assert len(result.required_outcomes) == (0 if n_segments else 1)
    if not n_segments:
        raw = result.required_outcomes[0]
        assert raw.executed is True and raw.converged is True
        assert raw.usable is False and raw.reason == "endpoint_hei"
        leaves, expected = path_search._path_leaves_and_expected(
            result.segments, required_outcomes=result.required_outcomes)
        assert raw in leaves and raw.item_id in expected


@pytest.mark.parametrize("route", ["depth", "kink"])
@pytest.mark.parametrize("solver_converged", [True, False, None])
@pytest.mark.parametrize("bond_failure", [False, True])
def test_terminal_path_keeps_solver_fact_and_diagnostics(
    tmp_path, monkeypatch, capsys, route, solver_converged, bond_failure,
):
    from types import SimpleNamespace

    def image(energy):
        return SimpleNamespace(energy=energy, set_calculator=lambda _calc: None)

    def mep(converged):
        energies = [0.0, 1.0, 0.2]  # finite, interior HEI
        return SimpleNamespace(
            images=[image(e) for e in energies], energies=energies,
            hei_idx=1, is_converged=converged,
        )

    primary = mep(True)
    terminal = mep(solver_converged)
    responses = iter([primary, terminal] if route == "kink" else [terminal])
    mep_tags = []

    def run_mep(*_args, **kwargs):
        mep_tags.append(kwargs["tag"])
        return next(responses)

    monkeypatch.setattr(path_search, "_run_mep_between", run_mep)
    monkeypatch.setattr(path_search, "_optimize_single",
                        lambda geometry, *_args, **_kwargs: (geometry, True))
    monkeypatch.setattr(path_search, "_make_linear_interpolations",
                        lambda *_args: [image(1.0)])

    n_bond = 0
    def compare(*_args, **_kwargs):
        nonlocal n_bond
        n_bond += 1
        if route == "kink" and n_bond <= 4:
            # Detect a kink; no left/right recursion; then hit max_seq_kink=1.
            return False, "(no covalent changes detected)"
        if bond_failure:
            raise ValueError("injected bond diagnostic failure")
        return True, "Bond formed"

    is_p = path_search.__name__.startswith("pdb2reaction.")
    monkeypatch.setattr(path_search,
                        "has_bond_change" if is_p else "_has_bond_change", compare)
    kwargs = dict(
        geom_cfg={}, gs_cfg={}, stopt_cfg={}, single_opt_cfg={}, bond_cfg={},
        search_cfg={"max_depth": 0 if route == "depth" else 2,
                    "max_seq_kink": 1, "kink_max_nodes": 1,
                    "stitch_rmsd_thresh": 1e-4, "bridge_rmsd_thresh": 1e-4},
        refine_mode_kind="peak", mep_mode_kind="gsm", out_dir=tmp_path,
        ref_pdb_path=None, depth=1 if route == "depth" else 0,
        seg_counter=[0], branch_tag="pair_00", calc_cfg={}, dmf_cfg={},
    )
    if is_p:
        kwargs.update(single_opt_kind="lbfgs", prepared_inputs=[], prepared_input=None)
    result = path_search._build_multistep_path(
        primary.images[0], primary.images[-1], None, **kwargs
    )
    assert len(mep_tags) == (2 if route == "kink" else 1)
    assert n_bond == (5 if route == "kink" else 1)
    assert len(result.segments) == 1
    assert not result.required_outcomes
    segment = result.segments[0]
    assert segment.converged is solver_converged
    assert segment.summary == (
        "(bond-change evaluation failed)" if bond_failure else "Bond formed"
    )
    segment.seg_index = 1
    leaves, expected = path_search._path_leaves_and_expected(result.segments)
    assert expected == ["segment_1"]
    assert leaves[0].converged is solver_converged
    assert leaves[0].usable is (solver_converged is True)
    captured = capsys.readouterr()
    text = captured.out + captured.err
    if bond_failure:
        assert "injected bond diagnostic failure" in text
    if route == "kink":
        assert "Consecutive kink segments" in text
