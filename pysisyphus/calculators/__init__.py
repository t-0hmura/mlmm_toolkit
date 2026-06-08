"""Public calculator surface of the bundled pysisyphus fork.

Only the abstract ``Calculator`` base and the orientation-projected
``Dimer`` TS-force calculator are exposed. The historical 30+ QM-chemistry
backends (Gaussian / ORCA / MOPAC / Psi4 / Turbomole / ...) and the
``pysisyphus run`` CLI driver are not part of this fork.
"""
import logging

__all__ = [
    "Calculator",
    "Dimer",
]


from pysisyphus.calculators.Calculator import Calculator
from pysisyphus.calculators.Dimer import Dimer


logger = logging.getLogger("dimer")
logger.setLevel(logging.DEBUG)
handler = logging.FileHandler("dimer.log", mode="w", delay=True)
fmt_str = "%(message)s"
formatter = logging.Formatter(fmt_str)
handler.setFormatter(formatter)
# logger.addHandler(handler)
