import logging
import sys

logger = logging.getLogger("pysisyphus")
logger.setLevel(logging.DEBUG)

file_handler = logging.FileHandler("pysisyphus.log", mode="w", delay=True)
# logger.addHandler(file_handler)

# silenced: leaked pysisyphus logger output outside the CLI -v gate
# stdout_handler = logging.StreamHandler(sys.stdout)
# stdout_handler.setLevel(logging.INFO)
# logger.addHandler(stdout_handler)
