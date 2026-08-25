"""GSM progress tables retain their segment section headings."""

from __future__ import annotations

import inspect


def test_recursive_and_standalone_gsm_emit_sections_before_optimization() -> None:
    import mlmm.workflows.path_opt as path_opt
    import mlmm.workflows.path_search as path_search

    recursive = inspect.getsource(path_search._run_gsm_between)
    tagged = 'emit(f"\\n====== [{tag}] GSM ======\\n", narrative=True)'
    assert tagged in recursive
    assert recursive.index(tagged) < recursive.index("optimizer.run()")

    standalone = inspect.getsource(path_opt.cli.callback)
    generic = 'emit("\\n====== Growing String optimization ======\\n", narrative=True)'
    assert generic in standalone
    assert standalone.index(generic) < standalone.index(
        "optimizer.run()", standalone.index(generic)
    )
