# log_utils.py
"""
Author : Ryan Sonderman
Date   : 2025-06-20
Version: 1.0.0
Purpose: Logging script
"""

import logging

logging.basicConfig(
    format='%(asctime)s [%(levelname)s]: %(message)s',
    level=logging.INFO
)
logger = logging.getLogger()

def log(msg, level="INFO"):
    """
    Logging wrapper to log messages with different levels.
    
    Args:
        msg (str): Message to be saved in the log
        level (str): Level of the log (e.g. INFO, WARNING, ERROR)
    """
    if level == "INFO":
        logger.info(msg)
    elif level == "WARNING":
        logger.warning(msg)
    elif level == "ERROR":
        logger.error(msg)