from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest

from mlmm.workflows import path_search
from mlmm.workflows.path_opt import _select_hei_index as select_path_opt_hei
from mlmm.workflows.path_search import _select_hei_index as select_path_search_hei
from pysisyphus.constants import AU2KCALPERMOL


@pytest.mark.parametrize("selector", [select_path_opt_hei, select_path_search_hei])
@pytest.mark.parametrize(
    ("energies", "expected"),
    [
        ([0.0, 1.0, 2.0], 2),
        ([0.0, 2.0, 1.0, 3.0], 1),
        ([3.0, 2.0, 1.0], 0),
        ([0.0, 3.0, 1.0], 1),
    ],
)
def test_hei_prefers_a_local_peak_over_higher_endpoints(selector, energies, expected):
    assert selector(energies) == expected


@pytest.mark.parametrize("drop_shared_start", [False, True])
def test_segment_peak_survives_stitching_a_shared_first_frame(drop_shared_start):
    from mlmm.workflows.path_search import SegmentReport, _segment_hei_index

    images = [SimpleNamespace() for _ in range(4)]
    energies = [0.0, 2.0, 1.0, 3.0]
    report = SegmentReport("rxn", 2.0, 3.0, "Bond formed", peak_image=images[1])
    start = int(drop_shared_start)
    assert _segment_hei_index(images[start:], energies[start:], report) == 1 - start
    assert report.barrier_kcal == 2.0


def test_monotone_bridge_does_not_overwrite_a_reactive_segment_barrier(tmp_path, monkeypatch):
    from mlmm.workflows import path_search

    left = SimpleNamespace(coords3d=np.zeros((1, 3)))
    right = SimpleNamespace(coords3d=np.ones((1, 3)))
    bridge_images = [left, SimpleNamespace(coords3d=np.full((1, 3), 0.5)), right]
    bridge = SimpleNamespace(images=bridge_images, energies=[0.0, 1.0, 2.0], hei_idx=2, is_converged=True)
    monkeypatch.setattr(path_search, "_run_mep_between", lambda *_a, **_k: bridge)
    report = path_search.SegmentReport("rxn", 5.0, -1.0, "Bond formed")
    reports = [report]
    path_search._stitch_paths(
        [([left], [0.0]), ([right], [2.0])],
        stitch_rmsd_thresh=1e-6, bridge_rmsd_thresh=1e-6,
        shared_calc=None, gs_cfg={}, stopt_cfg={}, out_dir=tmp_path,
        tag="pair", ref_pdb_path=None, segments_out=reports,
    )
    assert report.barrier_kcal == 5.0
    assert len(reports) == 2
    assert reports[1].kind == "bridge" and reports[1].barrier_kcal is None
    assert reports[1].peak_image is right


@pytest.mark.parametrize("selector", [select_path_opt_hei, select_path_search_hei])
@pytest.mark.parametrize("energies", [[], [0.0, np.nan], [0.0, np.inf]])
def test_hei_selection_rejects_missing_or_nonfinite_energies(selector, energies):
    with pytest.raises(ValueError):
        selector(energies)


def test_frame_ranges_preserve_bridge_and_disjoint_segment_provenance():
    images = [
        *[SimpleNamespace(mep_seg_index=1, mep_seg_kind="seg") for _ in range(3)],
        *[SimpleNamespace(mep_seg_index=2, mep_seg_kind="bridge") for _ in range(2)],
        *[SimpleNamespace(mep_seg_index=3, mep_seg_kind="seg") for _ in range(3)],
        SimpleNamespace(mep_seg_index=1, mep_seg_kind="seg"),
    ]

    ranges = path_search._frame_ranges_by_segment(images)

    assert ranges[1] == {
        "kind": "seg",
        "frame_ranges": [[0, 3], [8, 9]],
    }
    assert ranges[2] == {
        "kind": "bridge",
        "frame_ranges": [[3, 5]],
        "frame_start": 3,
        "frame_stop": 5,
    }
    assert ranges[3]["frame_ranges"] == [[5, 8]]


def test_kink_step_is_reported_as_kink_not_as_a_ts_segment(monkeypatch, tmp_path):
    """A refined step whose two ends keep one bonding pattern is a kink.

    The kink keeps its frames in the path, but its report, its frame tags and
    the segment entries that ``all`` reads say ``kink`` (with no MEP barrier),
    so ``all --tsopt`` does not count it as a TS segment.
    """
    from mlmm.workflows.all import _is_reactive_segment

    def image(state, energy):
        # ``state`` stands for the bonding pattern that _has_bond_change compares.
        return SimpleNamespace(
            state=state, energy=energy, coords3d=np.zeros((1, 3)),
            set_calculator=lambda _calc: None,
        )

    def mep(images, hei):
        return SimpleNamespace(
            images=images, energies=[im.energy for im in images],
            hei_idx=hei, is_converged=True,
        )

    reactant, product = image("R", 0.0), image("P", -0.002)
    # Both neighbours of the outer HEI already have the product bonds, so
    # End1 -> End2 is a kink and the bond change lies between R and End1.
    outer = mep([reactant, image("P", 0.004), image("P", 0.006),
                 image("P", 0.003), product], hei=2)
    inner = mep([reactant, image("R", 0.002), image("R", 0.012),
                 image("P", 0.005), outer.images[1]], hei=2)
    refined = mep([inner.images[1], image("R", 0.013), inner.images[3]], hei=1)

    monkeypatch.setattr(
        path_search, "_run_mep_between",
        lambda *_a, tag, **_k: {"seg_000": outer, "seg_001": inner}[tag],
    )
    monkeypatch.setattr(path_search, "_refine_between", lambda *_a, **_k: refined)
    monkeypatch.setattr(path_search, "_optimize_single",
                        lambda geometry, *_a, **_k: (geometry, True))
    monkeypatch.setattr(path_search, "_make_linear_interpolations",
                        lambda *_a: [image("P", 0.008)])
    monkeypatch.setattr(
        path_search, "_has_bond_change",
        lambda a, b, _cfg: (
            (True, "Bond formed") if a.state != b.state
            else (False, "(no covalent changes detected)")
        ),
    )

    result = path_search._build_multistep_path(
        reactant, product, None,
        geom_cfg={}, gs_cfg={}, stopt_cfg={}, single_opt_cfg={}, bond_cfg={},
        search_cfg={"max_depth": 3, "kink_max_nodes": 1,
                    "stitch_rmsd_thresh": 1e-4, "bridge_rmsd_thresh": 1e-4},
        refine_mode_kind="peak", mep_mode_kind="gsm", out_dir=tmp_path,
        ref_pdb_path=None, depth=0, seg_counter=[0], branch_tag="pair_00",
        single_opt_kind="lbfgs", calc_cfg={}, dmf_cfg={},
    )

    reaction, kink = result.segments
    assert (reaction.tag, reaction.kind) == ("seg_001_refine", "seg")
    assert reaction.barrier_kcal == pytest.approx((0.013 - 0.002) * AU2KCALPERMOL)
    assert (kink.tag, kink.kind) == ("seg_000_kink", "kink")
    assert kink.barrier_kcal is None
    assert kink.delta_kcal == pytest.approx((0.003 - 0.004) * AU2KCALPERMOL)

    # The cli numbers the segments along the path and the frames by their tags;
    # each summary.json entry takes ``kind`` last from its frame range.
    for index, report in enumerate(result.segments, 1):
        report.seg_index = index
    for im in result.images:
        im.mep_seg_index = {"seg_001_refine": 1, "seg_000_kink": 2}[im.mep_seg_tag]
    ranges = path_search._frame_ranges_by_segment(result.images)
    entries = [
        {"index": s.seg_index, "kind": s.kind, **ranges[s.seg_index]}
        for s in result.segments
    ]
    assert [entry["kind"] for entry in entries] == ["seg", "kink"]
    assert [entry["index"] for entry in entries if _is_reactive_segment(entry)] == [1]

    leaves, expected = path_search._path_leaves_and_expected(result.segments)
    assert expected == ["segment_1", "segment_2"]
    assert all(leaf.usable for leaf in leaves)


def test_recursive_refinement_preserves_resolved_no_climb(monkeypatch, tmp_path):
    captured = {}
    sentinel = object()

    def fake_run(*args, **kwargs):
        captured["gs_cfg"] = args[3]
        return sentinel

    monkeypatch.setattr(path_search, "_run_mep_between", fake_run)
    result = path_search._refine_between(
        object(),
        object(),
        object(),
        {"climb": False, "climb_lanczos": False},
        {},
        Path(tmp_path),
        "seg_000",
        None,
    )

    assert result is sentinel
    assert captured["gs_cfg"]["climb"] is False
    assert captured["gs_cfg"]["climb_lanczos"] is False


@pytest.mark.parametrize("level", [0, 1, 2, 3])
def test_mep_diagram_inputs_print_only_at_level_three(monkeypatch, capsys, level):
    from mlmm.core import utils

    monkeypatch.setattr(utils, "_VERBOSE_LEVEL", level)
    path_search._echo_mep_diagram_inputs(["R", "TS1", "P"], [0.0, 12.5, -3.25])

    out = capsys.readouterr().out
    if level >= 3:
        assert out == (
            '[diagram] build_energy_diagram.labels = ["R", "TS1", "P"]\n'
            "[diagram] build_energy_diagram.energies_kcal = "
            "[0.000000, 12.500000, -3.250000]\n"
        )
    else:
        assert out == ""
