import logging

__all__ = [
    "LBFGS",
    "RFOptimizer",
    "StringOptimizer",
]

logger = logging.getLogger("optimizer")
logger.setLevel(logging.DEBUG)
handler = logging.FileHandler("optimizer.log", mode="w", delay=True)
fmt_str = "%(message)s"
formatter = logging.Formatter(fmt_str)
handler.setFormatter(formatter)
# logger.addHandler(handler)
