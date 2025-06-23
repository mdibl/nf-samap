# log_utils.py
"""
Author : Ryan Sonderman
Date   : 2025-06-20
Version: 1.0.0
Purpose: Logging script
"""

import logging
from contextlib import redirect_stdout
import io

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
    elif level == "DEBUG":
        logger.debug(msg)


def capture_output(func, *args, **kwargs):
    """
    Captures and logs stdout from a function and passes output through.
    
    Args:
        func (func): Function to capture output from
    """
    # Create a StringIO buffer to capture the stdout
    captured_output = io.StringIO()

    # Redirect stdout to the StringIO buffer and call the function
    with redirect_stdout(captured_output):
        func(*args, **kwargs)

    # Log captured output line by line
    for line in captured_output.getvalue().splitlines():
        logging.info(line)