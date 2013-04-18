import logging
import os
import sys

##
# Initializes the log file in the proper format.
#
# @param filename Path to a file. Or None if logging is to
#                 be disabled.
#
def init_log(filename):
    if filename:
        logging.basicConfig( filename = filename,
                             filemode = "w",
                             format = '%(asctime)s: %(message)s',
                             level = logging.INFO )
        logging.getLogger( ).svsim_init = True
    else:
        logging.getLogger( ).disabled = True

##
# Returns a stream to the root log file.
#
# Note: I could not find a clean way to find the stream
# of a log file, so this function is currently a bit of
# a hack. It will only return the stream of the root
# logger if init_log has been called. However, if logging 
# is enabled for logging.INFO then sys.stdout will be
# returned. In all other cases os.devnull will be returned.
#
# @return A stream to the root log file.
#
def get_log_stream():
    log = logging.getLogger( )
    try:
        if log.svsim_init and not log.disabled:
            return log.handlers[ 0 ].stream
        else:
            return open( os.devnull, "w" )
    except AttributeError:
        if log.isEnabledFor( logging.INFO ):
            return sys.stdout
        else:
            return open( os.devnull, "w" )

