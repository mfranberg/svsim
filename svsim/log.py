from logbook import FileHandler, StderrHandler

def init_log(filename):
    """
    Initializes the log file in the proper format.
    
    filename (str): Path to a file. Or None if logging is to
                     be disabled.
    """
    if filename:
        log_handler = FileHandler(filename)
    else:
        log_handler = StderrHandler()
    
    return log_handler

