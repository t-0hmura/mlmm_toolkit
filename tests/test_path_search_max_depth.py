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
