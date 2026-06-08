"""L5 Foundation — default values, shared utilities, logging.

Modules:
- ``defaults`` — single source of truth for keyword arguments (``MLMM_CALC_KW``,
  ``OPT_BASE_KW``, ``LBFGS_KW``, ``RSIRFO_KW``, ``BIAS_KW``, ``IRC_KW``, ``FREQ_KW``,
  ``DFT_KW``, ``MM_BACKEND_KW``, etc.) and B-factor constants.
- ``utils`` — pure helpers (YAML parse, atom-index parsing, freeze-atom resolution,
  PDB metadata, format helpers, pretty-printers).
- ``logging`` — ``setup_logging(verbose)`` driven by the ``-v`` / ``-vv`` root flag;
  configures stdlib logging level (WARNING / INFO / DEBUG).

Dependency direction: this layer **does not import any other in-package layer**
(L1 / L2 / L3 / L4). It is the leaf.
"""
