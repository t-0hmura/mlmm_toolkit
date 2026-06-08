"""L2 Application — per-subcommand orchestration (one module per ``mlmm`` subcommand).

Each module exposes a ``cli`` callable wired to the corresponding CLI subcommand
(``mlmm <name>``). The ``_LAZY_SUBCOMMANDS`` registry in ``mlmm.cli.app`` lazy-loads
these modules on demand to keep ``mlmm --help`` startup fast.

Modules:
- ``all`` — end-to-end pipeline (extract → MEP → TS → IRC → freq → DFT).
- ``opt`` — optimization (LBFGS grad / RFO hess, with optional microiteration).
- ``tsopt`` — TS optimization (RS-I-RFO / Dimer / Bofill flatten loop).
- ``freq`` — vibrational analysis + thermochemistry.
- ``irc`` — IRC integration (EulerPC).
- ``path_opt`` / ``path_search`` — reaction path search and optimization.
- ``scan`` / ``scan2d`` / ``scan3d`` — bond-length scans (1D/2D/3D).
- ``extract`` — pocket / model-region extractor.
- ``dft`` — single-point DFT on the ML region (PySCF / GPU4PySCF).
- ``define_layer`` — assign ML / MM layers via B-factor encoding.
- ``mm_parm`` — Amber parm7/rst7 topology generation.
- ``oniom_export`` / ``oniom_import`` — Gaussian g16 / ORCA ONIOM input I/O.

Internal helper module (not a subcommand):
- ``align_freeze`` — coordinate alignment + freeze-atom selection.
"""
