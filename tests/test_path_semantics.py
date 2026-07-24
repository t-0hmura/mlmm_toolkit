from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest

from mlmm.workflows import path_search
from mlmm.workflows.path_opt import _select_hei_index as select_path_opt_hei
from mlmm.workflows.path_search import _select_hei_index as select_path_search_hei


@pytest.mark.parametrize("selector", [select_path_opt_hei, select_path_search_hei])
@pytest.mark.parametrize(
    ("energies", "expected"),
    [
        ([0.0, 1.0, 2.0], 2),
        ([0.0, 2.0, 1.0, 3.0], 3),
        ([3.0, 2.0, 1.0], 0),
        ([0.0, 3.0, 1.0], 1),
    ],
)
def test_hei_is_the_global_energy_maximum(selector, energies, expected):
    assert selector(energies) == expected


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
