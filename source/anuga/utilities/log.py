#!/usr/bin/env python

'''
A simple logging module that logs to the console and a logfile, and has a
configurable threshold loglevel for each of console and logfile output.

Use it this way:
    import anuga.utilities.log as log
    log.debug('A message at DEBUG level')
    log.info('Another message, INFO level')

This class uses the 'borg' pattern - there is never more than one instance
of log data.  See
<http://www.suttoncourtenay.org.uk/duncan/accu/pythonpatterns.html>
for the basic idea used here: modules *are* singletons!

Until the first call to log() the user is free to play with the module data
to configure the logging.
'''

import sys
import os.path
import traceback
import logging


################################################################################
# Module variables - only one copy of these, ever.
################################################################################

# flag variable to determine if logging set up or not
_setup = False

# logging level for the console
console_logging_level = logging.INFO

# logging level for the logfile
log_logging_level = logging.DEBUG

# The name of the file to log to.
log_filename = './anuga.log'


################################################################################
# Module code.
################################################################################

##
# @brief Log a message at a specified level.
# @param level The loglevel to log with (logging.DEBUG, etc).
# @param msg Message string to log.
# @note First call of this method initializes the logging system.
def log(level, msg):
    '''Log a message at a particular loglevel.

    The first call to this method (by anybody) initializes logging and
    then logs the message.  Subsequent calls just log the message.
    '''

    global _setup

    # have we been setup?
    if not _setup:
        # setup the file logging system
        fmt = '%(asctime)s %(levelname)-8s %(mname)25s:%(lnum)-4d|%(message)s'
        logging.basicConfig(level=log_logging_level, format=fmt,
                            filename=log_filename, filemode='w')

        # define a console handler which writes to sys.stdout
        console = logging.StreamHandler(sys.stdout)
        console.setLevel(console_logging_level)
        formatter = logging.Formatter('%(message)s')
        console.setFormatter(formatter)
        logging.getLogger('').addHandler(console)

        # tell the world how we are set up
        start_msg = ("Logfile is '%s' with logging level of %s, "
                     "console logging level is %s"
                     % (log_filename,
                        logging.getLevelName(log_logging_level),
                        logging.getLevelName(console_logging_level)))
        logging.log(logging.CRITICAL, start_msg,
                    extra={'mname': __name__, 'lnum': 0})

        # mark module as *setup*
        _setup = True

    # get caller information - look back for first module != <this module name>
    frames = traceback.extract_stack()
    frames.reverse()
    for (mname, lnum, _, _) in frames:
        mname = os.path.basename(mname).rsplit('.', 1)[0]
        if mname != __name__:
            break

    logging.log(level, msg, extra={'mname': mname, 'lnum': lnum})

################################################################################
# Shortcut routines to make for simpler user code.
################################################################################

##
# @brief Shortcut for log(DEBUG, msg).
# @param msg Message string to log at logging.DEBUG level.
def debug(msg):
    log(logging.DEBUG, msg)

##
# @brief Shortcut for log(INFO, msg).
# @param msg Message string to log at logging.INFO level.
def info(msg):
    log(logging.INFO, msg)

##
# @brief Shortcut for log(WARNING, msg).
# @param msg Message string to log at logging.WARNING level.
def warning(msg):
    log(logging.WARNING, msg)

##
# @brief Shortcut for log(ERROR, msg).
# @param msg Message string to log at logging.ERROR level.
def error(msg):
    log(logging.ERROR, msg)

##
# @brief Shortcut for log(CRITICAL, msg).
# @param msg Message string to log at logging.CRITICAL level.
def critical(msg):
    log(logging.CRITICAL, msg)

