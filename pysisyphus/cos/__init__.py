import logging

__all__ = [
    "ChainOfStates",
    "GrowingChainOfStates",
    "GrowingString",
]

logger = logging.getLogger("cos")
logger.setLevel(logging.DEBUG)
handler = logging.FileHandler("cos.log", mode="w", delay=True)
fmt_str = "%(message)s"
formatter = logging.Formatter(fmt_str)
handler.setFormatter(formatter)
# logger.addHandler(handler)
