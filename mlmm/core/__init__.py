"""L5 Foundation — default values, shared utilities, logging.

Modules:
- ``defaults`` — single source of truth for keyword arguments (``MLMM_CALC_KW``,
  ``OPT_BASE_KW``, ``LBFGS_KW``, ``RSIRFO_KW``, ``BIAS_KW``, ``IRC_KW``, ``FREQ_KW``,
  ``DFT_KW``, ``MM_BACKEND_KW``, etc.) and B-factor constants.
- ``utils`` — pure helpers (YAML parse, atom-index parsing, freeze-atom resolution,
  PDB metadata, format helpers, pretty-printers).
- ``logging`` — ``setup_logging(verbose)`` for the per-subcommand
  ``-v/--verbose LEVEL`` control (0--3); level 3 enables DEBUG records.

Dependency direction (measured): this layer **does not import the L1 (cli) or
L2 (workflows) layers** — that one-way direction, plus the absence of any import
cycle, is what ``.github/scripts/check_import_graph.py`` enforces. A few core
utilities still reach *down* into ``backends`` / ``domain`` / ``io`` as
compatibility back-edges (e.g. ``core.utils`` uses ``domain.add_elem_info`` and
``io.structure_formats``); these introduce no cycle and are retired in the major
rewrite. The leaf ideal is a design intent, not a current fact.
"""
