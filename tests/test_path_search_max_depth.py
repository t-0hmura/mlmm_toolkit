"""`search.max_depth` counts LEVELS of recursive subdivision.

The cap therefore compares ``depth >= max_depth``: ``0`` performs no
subdivision at all and reproduces a single-segment MEP, which is the setting a
user reaches for when recursive splitting is not wanted. A deliberate ``0``
also keeps the ordinary ``seg_NNN`` tag, because ``_maxdepth`` means "the
recursion was cut off while covalent changes remained" and that segment is not
guaranteed to be a single elementary step.
"""

from __future__ import annotations

import inspect
from pathlib import Path

import click
import pytest

from mlmm.core.defaults import SEARCH_KW
from mlmm.workflows import path_search


SRC = Path(inspect.getfile(path_search)).read_text(encoding="utf-8")


def test_cap_counts_levels_so_zero_subdivides_nothing() -> None:
    assert 'max_depth = int(search_cfg.get("max_depth", SEARCH_KW["max_depth"]))' in SRC
    assert "if depth >= max_depth:" in SRC
    # `depth > max_depth` allows one split even at 0, leaving no way to switch
    # subdivision off and making the name mean "N+1 levels".
    assert 'if depth > int(search_cfg.get("max_depth"' not in SRC


def test_deliberate_zero_keeps_the_ordinary_segment_tag() -> None:
    assert "if max_depth > 0:" in SRC
    assert '_terminate_with_single_segment(\n                f"seg_{seg_counter[0]:03d}_maxdepth"' in SRC
    assert '_terminate_with_single_segment(f"seg_{seg_counter[0]:03d}")' in SRC
    # Every branch that abandons recursion publishes a segment record, so the
    # interval cannot vanish from `n_segments` / `segments[]`.
    assert "def _terminate_with_single_segment(" in SRC
    assert 'f"seg_{seg_counter[0]:03d}_kinklimit",' in SRC
    assert "segments=[]" in SRC  # the endpoint-HEI guard still refuses a report


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

@pytest.mark.parametrize(
    "max_depth,primary_hei,refined_hei,n_segments,single_opt_executed",
    [
        (0, 1, 1, 1, False),
        (1, 0, 1, 0, False),
        (1, 1, 0, 0, True),
        (1, 1, 2, 0, True),
        (1, 1, 1, 1, True),
    ],
)
def test_single_optimizer_provenance_and_refined_hei_boundary(
    tmp_path, monkeypatch, max_depth, primary_hei, refined_hei,
    n_segments, single_opt_executed,
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
        refine_mode_kind="peak", mep_mode_kind="gsm",
        out_dir=tmp_path, ref_pdb_path=None, depth=0, seg_counter=[0],
        branch_tag="pair_00",
    )

    assert len(optimizer_calls) == (2 if single_opt_executed else 0)
    assert result.single_opt_executed is single_opt_executed
    assert len(result.segments) == n_segments

    assert len(result.required_outcomes) == (0 if n_segments else 1)
    if not n_segments:
        raw = result.required_outcomes[0]
        assert raw.executed is True and raw.converged is True
        assert raw.usable is False and raw.reason == "endpoint_hei"
        leaves, expected = path_search._path_leaves_and_expected(
            result.segments, required_outcomes=result.required_outcomes)
        assert raw in leaves and raw.item_id in expected
