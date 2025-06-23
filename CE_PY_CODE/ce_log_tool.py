import logging, os, sys
from contextlib import contextmanager


@contextmanager
def printto(log: str, lst: str):
    """Redirect standard output and errors to files and yield a logger.

    Parameters
    ----------
    log : str
        Path to the debug log file.
    lst : str
        Path to the listing file capturing stdout.

    Yields
    ------
    logging.Logger
        Logger configured with handlers for the two files.
    """
    original_stdout = sys.stdout
    original_stderr = sys.stderr

    os.makedirs(os.path.dirname(log), exist_ok=True)
    os.makedirs(os.path.dirname(lst), exist_ok=True)

    logger = logging.getLogger("CE_Pipeline")
    logger.setLevel(logging.DEBUG)

    for h in logger.handlers[:]:
        logger.removeHandler(h)

    log_fh = logging.FileHandler(log, mode="w")
    log_fh.setLevel(logging.DEBUG)
    lst_fh = logging.FileHandler(lst, mode="w")
    lst_fh.setLevel(logging.INFO)

    formatter = logging.Formatter("%(asctime)s %(levelname)s: %(message)s")
    log_fh.setFormatter(formatter)
    lst_fh.setFormatter(formatter)

    logger.addHandler(log_fh)
    logger.addHandler(lst_fh)

    sys.stdout = open(lst, "a")
    sys.stderr = open(log, "a")

    try:
        yield logger
    finally:
        sys.stdout.close()
        sys.stderr.close()
        sys.stdout = original_stdout
        sys.stderr = original_stderr
        logger.removeHandler(log_fh)
        logger.removeHandler(lst_fh)
        log_fh.close()
        lst_fh.close()
