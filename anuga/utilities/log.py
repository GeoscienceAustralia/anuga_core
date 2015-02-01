#!/usr/bin/env python

'''
A simple logging module that logs to the console and a logfile, and has a
configurable threshold loglevel for each of console and logfile output.

Use it this way:
    import anuga.utilities.log as log

    # configure my logging
    log.console_logging_level = log.INFO
    log.log_logging_level = log.DEBUG
    log.log_filename = './my.log'

    # log away!
    log.debug('A message at DEBUG level')
    log.info('Another message, INFO level')

This class uses the 'borg' pattern - there is never more than one instance
of log data.  See the URL for the basic idea used here: modules *are*
singletons!

<http://www.suttoncourtenay.org.uk/duncan/accu/pythonpatterns.html>

Until the first call to log() the user is free to play with the module data
to configure the logging.

Note that this module uses some features of the logging package that were
introduced in python2.5.  If running on earlier versions, the following
features are disabled:
    . Calling module name + line number
'''

import os
import sys
import traceback
import logging
import datetime

DefaultConsoleLogLevel = logging.CRITICAL
DefaultFileLogLevel = logging.INFO
TimingDelimiter ='#@# '

################################################################################
# Module variables - only one copy of these, ever.
#
# The console logging level is set to a high level, like CRITICAL.  The logfile
# logging is set lower, between DEBUG and CRITICAL.  The idea is to log least to
# the console, but ensure that everything that goes to the console *will* also
# appear in the log file.  There is code to ensure log <= console levels.
#
# If console logging level is set to CRITICAL+1 then nothing will print on the
# console.
################################################################################

# flag variable to determine if logging set up or not
_setup = False

# logging level for the console
console_logging_level = DefaultConsoleLogLevel

# logging level for the logfile
log_logging_level = DefaultFileLogLevel

# The default name of the file to log to.


log_filename = os.path.join('.', 'anuga.log')

# set module variables so users don't have to do 'import logging'.
CRITICAL = logging.CRITICAL
ERROR = logging.ERROR
WARNING = logging.WARNING
INFO = logging.INFO
DEBUG = logging.DEBUG
NOTSET = logging.NOTSET

# set _new_python to True if python version 2.5 or later
_new_python = (sys.version_info[0]*10 + sys.version_info[1] >= 25)      # 2.5.x.x


################################################################################
# Module code.
################################################################################

def log(msg, level=None):
    '''Log a message at a particular loglevel.

    msg:    The message string to log.
    level:  The logging level to log with (defaults to console level).

    The first call to this method (by anybody) initializes logging and
    then logs the message.  Subsequent calls just log the message.
    '''

    global _setup, log_logging_level
    fname = '' # default to no frame name if it cannot be found
    lnum = 0

    # have we been setup?
    if not _setup:
        # sanity check the logging levels, require console >= file
        if log_logging_level > console_logging_level:
            log_logging_level = console_logging_level

        # setup the file logging system
        if _new_python:
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
        if _new_python:
            logging.log(logging.INFO, start_msg,
                        extra={'mname': __name__, 'lnum': 0})
        else:
            logging.log(logging.INFO, start_msg)

        # mark module as *setup*
        _setup = True

    # if logging level not supplied, assume console level
    if level is None:
        level = console_logging_level

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

    # why are we here? ... Oh yes! Log the message!
    if _new_python:
        logging.log(level, msg, extra={'mname': fname, 'lnum': lnum})
    else:
        logging.log(level, msg)


def log_exception_hook(type, value, tb):
    '''Hook function to process uncaught exceptions.

    type:   Type of exception.
    value:  The exception data.
    tb:     Traceback object.

    This has the same interface as sys.excepthook().
    '''

    msg = '\n' + ''.join(traceback.format_exception(type, value, tb))
    critical(msg)

    
################################################################################
# Shortcut routines to make for simpler user code.
################################################################################

def debug(msg=''):
    '''Shortcut for log(DEBUG, msg).'''

    log(msg, logging.DEBUG)


def info(msg=''):
    '''Shortcut for log(INFO, msg).'''

    log(msg, logging.INFO)


def warning(msg=''):
    '''Shortcut for log(WARNING, msg).'''

    log(msg, logging.WARNING)


def error(msg=''):
    '''Shortcut for log(ERROR, msg).'''

    log(msg, logging.ERROR)


def critical(msg=''):
    '''Shortcut for log(CRITICAL, msg).'''

    log(msg, logging.CRITICAL)

def timingInfo(msg=''):
    '''Shortcut for log(timingDelimiter, msg).'''

    log(TimingDelimiter + msg, logging.INFO)


def resource_usage(level=logging.INFO):
    '''Log memory usage at given log level.'''

    _scale = {'KB': 1024, 'MB': 1024*1024, 'GB': 1024*1024*1024,
              'kB': 1024, 'mB': 1024*1024, 'gB': 1024*1024*1024}

    if sys.platform != 'win32':
        _proc_status = '/proc/%d/status' % os.getpid()
        
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
        log(msg, level)
    else:
        # Windows code from: http://code.activestate.com/recipes/511491/
        try:
            import ctypes
            import _winreg
        except:
            log(level, 'Windows resource usage not available')
            return

        kernel32 = ctypes.windll.kernel32
        c_ulong = ctypes.c_ulong
        c_ulonglong = ctypes.c_ulonglong
        class MEMORYSTATUSEX(ctypes.Structure):
            _fields_ = [('dwLength', c_ulong),
                        ('dwMemoryLoad', c_ulong),
                        ('ullTotalPhys', c_ulonglong),
                        ('ullAvailPhys', c_ulonglong),
                        ('ullTotalPageFile', c_ulonglong),
                        ('ullAvailPageFile', c_ulonglong),
                        ('ullTotalVirtual', c_ulonglong),
                        ('ullAvailVirtual', c_ulonglong),
                        ('ullAvailExtendedVirtual', c_ulonglong)
                       ]

        memoryStatusEx = MEMORYSTATUSEX()
        memoryStatusEx.dwLength = ctypes.sizeof(MEMORYSTATUSEX)
        kernel32.GlobalMemoryStatusEx(ctypes.byref(memoryStatusEx))

        msg = ('Resource usage: total memory=%.1fMB free memory=%.1fMB'
               % (memoryStatusEx.ullTotalPhys/_scale['MB'],
                  memoryStatusEx.ullAvailPhys/_scale['MB']))
        log(msg, level)

def CurrentDateTime():
    return datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")

def TimeStamp():
    return datetime.datetime.now().strftime('%Y%m%d_%H%M%S')


def resource_usage_timing(level=logging.INFO, prefix =""):
    '''Log memory usage at given log level.'''

    _scale = {'KB': 1024, 'MB': 1024*1024, 'GB': 1024*1024*1024,
              'kB': 1024, 'mB': 1024*1024, 'gB': 1024*1024*1024}

    if sys.platform != 'win32':
        _proc_status = '/proc/%d/status' % os.getpid()
        
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
        log(msg, level)
        timingInfo('sys_platform, ' + sys.platform)
        timingInfo(prefix + 'memory, ' + str(memory()/_scale['MB']))
        timingInfo(prefix + 'resident, ' + str(resident()/_scale['MB']))
        timingInfo(prefix + 'stacksize, ' + str(stacksize()/_scale['MB']))
    else:
        # Windows code from: http://code.activestate.com/recipes/511491/
        try:
            import ctypes
            import _winreg
        except:
            log(level, 'Windows resource usage not available')
            return

        kernel32 = ctypes.windll.kernel32
        c_ulong = ctypes.c_ulong
        c_ulonglong = ctypes.c_ulonglong
        class MEMORYSTATUSEX(ctypes.Structure):
            _fields_ = [('dwLength', c_ulong),
                        ('dwMemoryLoad', c_ulong),
                        ('ullTotalPhys', c_ulonglong),
                        ('ullAvailPhys', c_ulonglong),
                        ('ullTotalPageFile', c_ulonglong),
                        ('ullAvailPageFile', c_ulonglong),
                        ('ullTotalVirtual', c_ulonglong),
                        ('ullAvailVirtual', c_ulonglong),
                        ('ullAvailExtendedVirtual', c_ulonglong)
                       ]

        memoryStatusEx = MEMORYSTATUSEX()
        memoryStatusEx.dwLength = ctypes.sizeof(MEMORYSTATUSEX)
        kernel32.GlobalMemoryStatusEx(ctypes.byref(memoryStatusEx))

        msg = ('Resource usage: total memory=%.1fMB free memory=%.1fMB'
               % (memoryStatusEx.ullTotalPhys/_scale['MB'],
                  memoryStatusEx.ullAvailPhys/_scale['MB']))
        log(msg, level)
        timingInfo('sys_platform, ' + sys.platform)
        timingInfo(prefix + 'total_memory, ' + str(memoryStatusEx.ullTotalPhys/_scale['MB']))
        timingInfo(prefix + 'free_memory, ' + str(memoryStatusEx.ullAvailPhys/_scale['MB']))

    
################################################################################
if __name__ == '__main__':
    critical('#' * 80)
    warning('Test of logging...')
    log('CRITICAL+1', CRITICAL+1)
    log('CRITICAL', CRITICAL)
    log('CRITICAL-1', CRITICAL-1)
    log('CRITICAL-2', CRITICAL-2)
    log('default - CRITICAL?')

    def test_it(num=100):
        if num > 0:
            test_it(num-1)
        else:
            resource_usage()

    import numpy as num
    
    a = num.zeros((1000,1000), num.float)

    info('sys.version_info=%s, _new_python=%s'
         % (str(sys.version_info), str(_new_python)))
    test_it()
