import os
import sys
import logging

def init_log(logger, filename=None, loglevel='INFO'):
    """
    Initializes the log file in the proper format.
    
    filename (str): Path to a file. Or None if logging is to
                     be disabled.
    loglevel (str): Determines the level of the log output.
    """
    
    formatter = logging.Formatter(
                    '[%(asctime)s] %(levelname)s: %(name)s: %(message)s'
                )
    logger.setLevel(getattr(logging, loglevel))
    
    # We will allways print warnings and hogher to stderr
    ch = logging.StreamHandler()
    ch.setLevel('WARNING')
    ch.setFormatter(formatter)
    logger.addHandler(ch)
    
    if filename:
        fi = logging.FileHandler(filename)
        fi.setLevel(getattr(logging, loglevel))
        fi.setFormatter(formatter)
        logger.addHandler(fi)
        logger.svsim_init = True
    

# Note: I could not find a clean way to find the stream
# of a log file, so this function is currently a bit of
# a hack. It will only return the stream of the root
# logger if init_log has been called. However, if logging
# is enabled for logging.INFO then sys.stdout will be
# returned. In all other cases os.devnull will be returned.
#
def get_log_stream():
    """
    Returns a stream to the root log file.
    
    Returns:
        A stream to the root log file.
    """
    log = logging.getLogger()
    try:
        if log.svsim_init and not log.disabled:
            return log.handlers[0].stream
        else:
            return open(os.devnull, "w")
    except AttributeError:
        if log.isEnabledFor(logging.INFO):
            return sys.stdout
        else:
            return open(os.devnull, "w")