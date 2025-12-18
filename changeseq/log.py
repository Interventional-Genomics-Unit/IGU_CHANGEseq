"""
log.py
=====

Setup logging utils for nested module logging

Adapted from the accepted answer here: http://stackoverflow.com/questions/7621897/python-logging-module-globally
"""

import logging
from datetime import datetime

# def createCustomLogger(name):
#     formatter = logging.Formatter(fmt='[%(asctime)s][%(levelname)s][%(module)s] %(message)s', datefmt='%m/%d %I:%M:%S%p')
#
#     handler = logging.StreamHandler()
#     handler.setFormatter(formatter)
#
#     logger = logging.getLogger(name)
#     logger.setLevel(logging.DEBUG)
#     logger.addHandler(handler)
#     return logger
#
# import logging

def createCustomLogger(name):
    formatter = logging.Formatter(
        fmt='[%(asctime)s][%(levelname)s][%(module)s] %(message)s',
        datefmt='%m/%d %I:%M:%S%p'
    )

    # StreamHandler for console output
    stream_handler = logging.StreamHandler()
    stream_handler.setFormatter(formatter)


    # Create and configure the logger
    logger = logging.getLogger(name)
    logger.setLevel(logging.DEBUG)

    # Avoid adding multiple handlers if already added
    if not any(isinstance(h, logging.StreamHandler) for h in logger.handlers):
        logger.addHandler(stream_handler)

    return logger


def addFileHandlerToLogger(logger, base_log_filename='cs_pipeline.log', date_format='%Y-%m-%d_%H:%M'):
    """
    Adds a file handler to an existing logger.
    Useful when you want to start logging to a file later.
    """


    date_str = datetime.now().strftime(date_format)
    name_parts = base_log_filename.rsplit('.', 1)
    if len(name_parts) == 2:
        filename_with_date = f"{name_parts[0]}_{date_str}.{name_parts[1]}"
    else:
        filename_with_date = f"{base_log_filename}_{date_str}"

    formatter = logging.Formatter(
        fmt='[%(asctime)s][%(levelname)s][%(module)s] %(message)s',
        datefmt='%m/%d %I:%M:%S%p'
    )

    file_handler = logging.FileHandler(filename_with_date)
    file_handler.setFormatter(formatter)

    # Prevent adding duplicate file handlers
    if not any(isinstance(h, logging.FileHandler) and h.baseFilename == file_handler.baseFilename for h in logger.handlers):
        logger.addHandler(file_handler)