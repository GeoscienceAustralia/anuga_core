#!/usr/bin/env python

'''
A simple logging module that logs to the console and a logfile, and has a
configurable threshold loglevel for each of console and logfile output.

Use it this way:
    import anuga.utilities.log as log
    log.console_logging_level = log.DEBUG
    log.debug('A message at DEBUG level')
    log.info('Another message, INFO level')

This class uses the 'borg' pattern - there is never more than one instance
of log data.  See
<http://www.suttoncourtenay.org.uk/duncan/accu/pythonpatterns.html>
for the basic idea used here: modules *are* singletons!

Until the first call to log() the user is free to play with the module data
to configure the logging.

Note that this module uses features of the logging package that were introduced
in python2.5.  If running on earlier versions, these features are disabled:
    . Calling module name + line number
'''

import os
import sys
import traceback
import logging


################################################################################
# Module variables - only one copy of these, ever.
#
# The console logging level is set to a high level, like CRITICAL.  The logfile
# logging is set lower, between DEBUG and CRITICAL.  The idea is to log least to
# the console, but ensure that everything that goes to the console *will* also
# appear in the log file.  There is code to ensure log <= console levels.
################################################################################

# flag variable to determine if logging set up or not
_setup = False

# logging level for the console
console_logging_level = logging.CRITICAL

# logging level for the logfile
log_logging_level = logging.INFO

# The default name of the file to log to.
log_filename = os.path.join('.', 'anuga.log')

# set module variables so users don't have to do 'import logging'.
CRITICAL = logging.CRITICAL
ERROR = logging.ERROR
WARNING = logging.WARNING
INFO = logging.INFO
DEBUG = logging.DEBUG
NOTSET = logging.NOTSET

# set True if python version 2.5 or later
(version_major, version_minor, _, _, _) = sys.version_info
new_python = ((version_major == 2 and version_minor >= 5) or version_major > 2)


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

    global _setup, log_logging_level

    # have we been setup?
    if not _setup:
        # sanity check the logging levels, require console >= file
        if log_logging_level > console_logging_level:
            log_logging_level = console_logging_level

        # setup the file logging system
        if new_python:
            fmt = '%(asctime)s %(levelname)-8s %(mname)25s:%(lnum)-4d|%(message)s'
        else:
            fmt = '%(asctime)s %(levelname)-8s|%(message)s'
        logging.basicConfig(level=log_logging_level, format=fmt,
                            filename=log_filename, filemode='w')

        # define a console handler which writes to sys.stdout
        console = logging.StreamHandler(sys.stdout)
        console.setLevel(console_logging_level)
        formatter = logging.Formatter('%(message)s')
        console.setFormatter(formatter)
        logging.getLogger('').addHandler(console)

        # catch exceptions
        sys.excepthook = log_exception_hook

        # tell the world how we are set up
        start_msg = ("Logfile is '%s' with logging level of %s, "
                     "console logging level is %s"
                     % (log_filename,
                        logging.getLevelName(log_logging_level),
                        logging.getLevelName(console_logging_level)))
        if new_python:
            logging.log(logging.CRITICAL, start_msg,
                        extra={'mname': __name__, 'lnum': 0})
        else:
            logging.log(logging.CRITICAL, start_msg)

        # mark module as *setup*
        _setup = True

    # get caller information - look back for first module != <this module name>
    frames = traceback.extract_stack()
    frames.reverse()
    try:
        (_, mod_name) = __name__.rsplit('.', 1)
    except ValueError:
        mod_name = __name__
    for (fpath, lnum, mname, _) in frames:
        (fname, _) = os.path.basename(fpath).rsplit('.', 1)
        if fname != mod_name:
            break

    if new_python:
        logging.log(level, msg, extra={'mname': fname, 'lnum': lnum})
    else:
        logging.log(level, msg)

##
# @brief Hook function to process uncaught exceptions.
# @param type
# @param value
# @param traceback
# @note Same interface as sys.excepthook()
def log_exception_hook(type, value, tb):
    msg = ''.join(traceback.format_exception(type, value, tb))
    critical(msg)

    
################################################################################
# Shortcut routines to make for simpler user code.
################################################################################

##
# @brief Shortcut for log(DEBUG, msg).
# @param msg Message string to log at logging.DEBUG level.
def debug(msg=''):
    log(logging.DEBUG, msg)

##
# @brief Shortcut for log(INFO, msg).
# @param msg Message string to log at logging.INFO level.
def info(msg=''):
    log(logging.INFO, msg)

##
# @brief Shortcut for log(WARNING, msg).
# @param msg Message string to log at logging.WARNING level.
def warning(msg=''):
    log(logging.WARNING, msg)

##
# @brief Shortcut for log(ERROR, msg).
# @param msg Message string to log at logging.ERROR level.
def error(msg=''):
    log(logging.ERROR, msg)

##
# @brief Shortcut for log(CRITICAL, msg).
# @param msg Message string to log at logging.CRITICAL level.
def critical(msg=''):
    log(logging.CRITICAL, msg)

##
# @brief Log memory usage at time of call.
# @param level Override the default INFO logging level.
# @note From http://code.activestate.com/recipes/286222/.
def resource_usage(level=logging.INFO):
    '''Log memory usage at given log level.'''

    if sys.platform != 'win32':
        _proc_status = '/proc/%d/status' % os.getpid()
        _scale = {'KB': 1024, 'MB': 1024*1024, 'GB': 1024*1024*1024,
                  'kB': 1024, 'mB': 1024*1024, 'gB': 1024*1024*1024}

        def _VmB(VmKey):
            '''Get number of virtual bytes used.'''

            # get pseudo file /proc/<pid>/status
            try:
                t = open(_proc_status)
                v = t.read()
                t.close()
            except IOError:
                return 0.0

            # get VmKey line, eg: 'VmRSS: 999 kB\n ...
            i = v.index(VmKey)
            v = v[i:].split(None, 3)
            if len(v) < 3:
                return 0.0

            # convert Vm value to bytes
            return float(v[1]) * _scale[v[2]]

        def memory(since=0.0):
            '''Get virtual memory usage in bytes.'''

            return _VmB('VmSize:') - since

        def resident(since=0.0):
            '''Get resident memory usage in bytes.'''

            return _VmB('VmRSS:') - since

        def stacksize(since=0.0):
            '''Get stack size in bytes.'''

            return _VmB('VmStk:') - since

        msg = ('Resource usage: memory=%.1fMB resident=%.1fMB stacksize=%.1fMB'
               % (memory()/_scale['MB'], resident()/_scale['MB'],
                  stacksize()/_scale['MB']))
        log(level, msg)
    else:
        msg = ('Sorry, no memory statistics for Windows (yet).')
        log(level, msg)


if __name__ == '__main__':
##    critical('Testing exception capturing')
    def test_it(num=100):
        if num > 0:
            test_it(num-1)
        else:
            resource_usage()

    import numpy as num
    
    a = num.zeros((1000,1000), num.float)

    test_it()
