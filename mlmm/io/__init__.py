"""L4 Infra — output / I/O modules (no chemistry logic).

Modules:
- ``energy_diagram`` — render energy diagrams from numeric values (Plotly).
- ``hessian_cache`` — ``--read-hess`` / ``--dump-hess`` Hessian cache I/O.
- ``hessian_calc`` — Hessian computation helpers (CPU/GPU dispatch).
- ``pdb_fix`` — ``mlmm fix-altloc`` subcommand backend (PDB altloc resolution).
- ``summary`` — per-run ``summary.json`` / ``summary.log`` writer.
- ``trj2fig`` — trajectory-to-figure plotting (``mlmm trj2fig``).

Note: harmonic-restraint setup lives in ``mlmm/workflows/restraints.py``
(L2 layer), not here.
"""
