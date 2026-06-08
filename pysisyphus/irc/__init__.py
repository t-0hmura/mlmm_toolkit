import logging

__all__ = [
    "EulerPC",
]

from pysisyphus.irc.EulerPC import EulerPC

logger = logging.getLogger("irc")
logger.setLevel(logging.DEBUG)
handler = logging.FileHandler("irc.log", mode="w", delay=True)
fmt_str = "%(message)s"
formatter = logging.Formatter(fmt_str)
handler.setFormatter(formatter)
# logger.addHandler(handler)
