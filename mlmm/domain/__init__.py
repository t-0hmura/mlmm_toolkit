"""L3 Domain — chemistry helper logic (bond detection, element repair, structure summary).

Modules:
- ``add_elem_info`` — repair/add element columns in PDB files.
- ``bond_changes`` — detect bond changes between two structures (R↔P, R↔TS, TS↔P).
- ``bond_summary`` — summarise detected bond changes for human-readable reports.

These modules may use numeric backends (``torch`` / ``numpy`` /
``pysisyphus.constants``) but must **not** import MLIP SDKs
(``fairchem`` / ``orb_models`` / ``mace`` / ``aimnet``). The ``check_domain_pure``
gate enforces this isolation.
"""
