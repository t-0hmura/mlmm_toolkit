"""``all --hessian-calc-mode`` reaches the IRC child, as it does tsopt and freq."""

from __future__ import annotations

from pathlib import Path
from types import SimpleNamespace

import pytest

from mlmm.workflows import all as all_workflow


class ReachedChild(BaseException):
    """Stop at the IRC child dispatch."""


def _ts_pdb(path: Path) -> Path:
    path.write_text(
        "HETATM    1  C1  MOL A   1       0.000   0.000   0.000  1.00  0.00           C\n"
        "HETATM    2  O1  MOL A   1       1.200   0.000   0.000  1.00  0.00           O\n"
        "END\n",
        encoding="utf-8",
    )
    return path


@pytest.mark.parametrize("mode", ["Analytical", None])
def test_irc_child_receives_hessian_calc_mode(tmp_path, monkeypatch, mode) -> None:
    captured = {}

    def child(name, _cli, args, **kwargs):
        captured["name"] = name
        captured["args"] = list(args)
        raise ReachedChild

    monkeypatch.setattr(all_workflow, "_run_cli_main", child)
    ts_pdb = _ts_pdb(tmp_path / "ts.pdb")

    with pytest.raises(ReachedChild):
        all_workflow._irc_and_match(
            seg_idx=1,
            seg_dir=tmp_path / "seg",
            ref_pdb_for_seg=ts_pdb,
            seg_pocket_pdb=ts_pdb,
            g_ts=None,
            q_int=0,
            spin=1,
            resolved_calc_template=SimpleNamespace(materialize=lambda: {}),
            real_parm7=tmp_path / "system.parm7",
            model_pdb=ts_pdb,
            hessian_calc_mode=mode,
        )

    args = captured["args"]
    assert captured["name"] == "irc"
    if mode is None:
        assert "--hessian-calc-mode" not in args
    else:
        assert args[args.index("--hessian-calc-mode") + 1] == mode
