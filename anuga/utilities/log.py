#!/usr/bin/env python

"""
A simple logging module that logs to the console and a logfile.

Basic usage (print + log both go to terminal and file):

    import anuga.utilities.log as log

    log.set_logfile('./my.log')   # activates tee to file

    log.debug('A message at DEBUG level')
    log.info('Another message, INFO level')
    print('This also goes to both terminal and file')

Or via the public API:

    import anuga
    anuga.set_logfile('./my.log')

Level defaults when a logfile is active:
    console: INFO   (info/warning/error/critical visible on terminal)
    file:    DEBUG  (everything recorded in the file)

This module uses the 'borg' pattern — modules are singletons.
"""

import os
import sys
import traceback
import logging
from contextlib import contextmanager


DefaultConsoleLogLevel = logging.INFO
DefaultFileLogLevel = logging.DEBUG
TimingDelimiter = '#@# '

################################################################################
# TeeStream — write to both terminal and a log file simultaneously
################################################################################

class TeeStream:
    """Tee sys.stdout to a file: every write goes to both terminal and file.

    Usage:
        sys.stdout = TeeStream('run.log')
        # From now on, print() and sys.stdout.write() go to both places.
        sys.stdout.close()   # when done (optional)
    """

    def __init__(self, logfile_path, mode='a'):
        self._terminal = sys.__stdout__
        self._log = open(logfile_path, mode, encoding='utf-8')

    def write(self, message):
        self._terminal.write(message)
        self._terminal.flush()
        self._log.write(message)
        self._log.flush()

    def flush(self):
        self._terminal.flush()
        self._log.flush()

    def close(self):
        self._log.close()

    # Proxy attribute reads to the underlying terminal so code that inspects
    # sys.stdout (e.g. checks for .encoding) still works.
    def __getattr__(self, name):
        return getattr(self._terminal, name)


class _FileOnlyStream:
    """Write to the log file only — used by the file_only() context manager."""

    def __init__(self, log_fh):
        self._log = log_fh

    def write(self, message):
        self._log.write(message)
        self._log.flush()

    def flush(self):
        self._log.flush()

    def __getattr__(self, name):
        return getattr(self._log, name)


@contextmanager
def file_only():
    """Context manager: send all print() output to the log file only.

    Terminal output is suppressed for the duration of the block.
    Requires set_logfile() to have been called first; if no logfile
    is active, output is suppressed entirely.

    Typical use — capture verbose internal output without cluttering
    the terminal::

        with log.file_only():
            anuga.create_pmesh_from_regions(..., verbose=True, ...)
    """
    original = sys.stdout
    if isinstance(sys.stdout, TeeStream):
        sys.stdout = _FileOnlyStream(sys.stdout._log)
    else:
        # No logfile active — discard output
        import io
        sys.stdout = io.StringIO()
    try:
        yield
    finally:
        sys.stdout = original


################################################################################
# Module variables — only one copy, ever.
################################################################################

# flag: has logging been set up yet?
_setup = False

# logging level for the console handler
console_logging_level = DefaultConsoleLogLevel

# logging level for the file handler
log_logging_level = DefaultFileLogLevel

# Path to the log file.  None = file logging disabled (no file created).
log_filename = None

# set module variables so users don't have to do 'import logging'.
CRITICAL = logging.CRITICAL
ERROR    = logging.ERROR
WARNING  = logging.WARNING
INFO     = logging.INFO
DEBUG    = logging.DEBUG
NOTSET   = logging.NOTSET


################################################################################
# set_logfile — the main entry point for enabling file+tee logging
################################################################################

VERBOSE = logging.DEBUG  # level used by log.verbose() — file only by default


def set_logfile(path,
                console_level=DefaultConsoleLogLevel,
                file_level=DefaultFileLogLevel,
                verbose_to_screen=False):
    """Enable logging to *path*, tee-ing all print() output as well.

    After this call:
    - sys.stdout is replaced with a TeeStream so every print() goes to
      both the terminal and *path*.
    - log.info() writes to both terminal and file.
    - log.verbose() / log.debug() write to the file only (unless
      verbose_to_screen=True).
    - The previous log file (if any) is closed.

    Parameters
    ----------
    path : str
        File path for the log file.
    console_level : int
        Logging level for console output (default INFO).
        log.verbose() and log.debug() are below this threshold and go
        to the file only.
    file_level : int
        Logging level for file output (default DEBUG — everything).
    verbose_to_screen : bool
        If True, lower the console threshold to DEBUG so that
        log.verbose() output also appears on the terminal.  Useful
        when debugging without needing a clean screen.
    """
    if verbose_to_screen:
        console_level = logging.DEBUG
    global log_filename, console_logging_level, log_logging_level, _setup

    # Close any existing TeeStream
    if isinstance(sys.stdout, TeeStream):
        sys.stdout.close()

    log_filename = path
    console_logging_level = console_level
    log_logging_level = file_level
    _setup = False  # force re-initialisation on next log() call

    # Tee stdout so print() goes to both terminal and file
    sys.stdout = TeeStream(path)

    # Trigger logging setup now
    log('Logfile opened: ' + path, INFO)


################################################################################
# Module code.
################################################################################

def log(msg, level=None):
    '''Log a message at a particular loglevel.

    msg:    The message string to log.
    level:  The logging level to log with (defaults to console level).

    The first call to this method initialises the logging.FileHandler if a
    log_filename has been configured.
    '''

    global _setup, log_logging_level

    fname = ''
    lnum = 0

    if not _setup:
        # File logging: only if a filename has been configured
        if log_filename is not None:
            fmt = '%(asctime)s %(levelname)-8s %(mname)25s:%(lnum)-4d|%(message)s'
            file_handler = logging.FileHandler(log_filename, mode='a')
            file_handler.setLevel(log_logging_level)
            file_handler.setFormatter(logging.Formatter(fmt))

            root = logging.getLogger('')
            root.setLevel(min(log_logging_level, console_logging_level))

            # Remove any pre-existing handlers to avoid duplicates on re-init
            for h in root.handlers[:]:
                root.removeHandler(h)

            root.addHandler(file_handler)

            console = logging.StreamHandler(sys.__stdout__)
            console.setLevel(console_logging_level)
            console.setFormatter(logging.Formatter('%(message)s'))
            root.addHandler(console)
        else:
            # No file configured: just console at console_logging_level
            root = logging.getLogger('')
            root.setLevel(console_logging_level)
            for h in root.handlers[:]:
                root.removeHandler(h)
            console = logging.StreamHandler(sys.__stdout__)
            console.setLevel(console_logging_level)
            console.setFormatter(logging.Formatter('%(message)s'))
            root.addHandler(console)

        sys.excepthook = log_exception_hook
        _setup = True

    if level is None:
        level = console_logging_level

    # get caller information
    frames = traceback.extract_stack()
    frames.reverse()

    try:
        (_, mod_name) = __name__.rsplit('.', 1)
    except ValueError:
        mod_name = __name__

    for (fpath, lnum, mname, _) in frames:
        try:
            (fname, _) = os.path.basename(fpath).rsplit('.', 1)
        except ValueError:
            fname = __name__
        if fname != mod_name:
            break

    logging.log(level, msg, extra={'mname': fname, 'lnum': lnum})


def log_exception_hook(type, value, tb):
    '''Hook function to process uncaught exceptions.'''
    msg = '\n' + ''.join(traceback.format_exception(type, value, tb))
    critical(msg)


################################################################################
# Shortcut routines
################################################################################

def verbose(msg=''):
    """Log a verbose/internal message — goes to file only (not screen).

    Use this instead of print() inside ANUGA code that has a verbose flag.
    Output appears on screen only when set_logfile(..., verbose_to_screen=True).
    """
    log(msg, logging.DEBUG)

def debug(msg=''):
    """Log a DEBUG-level message (file only by default)."""
    log(msg, logging.DEBUG)

def info(msg=''):
    """Log an INFO-level message (terminal and file)."""
    log(msg, logging.INFO)

def warning(msg=''):
    """Log a WARNING-level message (terminal and file)."""
    log(msg, logging.WARNING)

def error(msg=''):
    """Log an ERROR-level message (terminal and file)."""
    log(msg, logging.ERROR)

def critical(msg=''):
    """Log a CRITICAL-level message (terminal and file)."""
    log(msg, logging.CRITICAL)

def timingInfo(msg=''):
    log(TimingDelimiter + msg, logging.INFO)


def resource_usage(level=logging.INFO):
    '''Log memory usage at given log level.'''

    _scale = {'KB': 1024, 'MB': 1024*1024, 'GB': 1024*1024*1024,
              'kB': 1024, 'mB': 1024*1024, 'gB': 1024*1024*1024}

    if sys.platform != 'win32':
        _proc_status = '/proc/%d/status' % os.getpid()

        def _VmB(VmKey):
            try:
                t = open(_proc_status)
                v = t.read()
                t.close()
            except IOError:
                return 0.0
            i = v.index(VmKey)
            v = v[i:].split(None, 3)
            if len(v) < 3:
                return 0.0
            return float(v[1]) * _scale[v[2]]

        def memory(since=0.0):
            return _VmB('VmSize:') - since

        def resident(since=0.0):
            return _VmB('VmRSS:') - since

        def stacksize(since=0.0):
            return _VmB('VmStk:') - since

        msg = ('Resource usage: memory=%.1fMB resident=%.1fMB stacksize=%.1fMB'
               % (memory() / _scale['MB'],
                  resident() / _scale['MB'],
                  stacksize() / _scale['MB']))
        log(msg, level)
    else:
        try:
            import ctypes
            import winreg
        except ImportError:
            log('Windows resource usage not available', level)
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
                        ('ullAvailExtendedVirtual', c_ulonglong)]

        memoryStatusEx = MEMORYSTATUSEX()
        memoryStatusEx.dwLength = ctypes.sizeof(MEMORYSTATUSEX)
        kernel32.GlobalMemoryStatusEx(ctypes.byref(memoryStatusEx))

        msg = ('Resource usage: total memory=%.1fMB free memory=%.1fMB'
               % (memoryStatusEx.ullTotalPhys / _scale['MB'],
                  memoryStatusEx.ullAvailPhys / _scale['MB']))
        log(msg, level)


def current_datetime():
    from datetime import datetime
    return datetime.now().strftime("%Y%m%d_%H%M%S%z")

def CurrentDateTime():
    from datetime import datetime
    return datetime.now().strftime("%Y-%m-%d_%H:%M:%S")

def TimeStamp():
    from datetime import datetime
    return datetime.now().strftime('%Y%m%d_%H%M%S')


def resource_usage_timing(level=logging.INFO, prefix=''):
    '''Log memory usage with timing info.'''

    _scale = {'KB': 1024, 'MB': 1024*1024, 'GB': 1024*1024*1024,
              'kB': 1024, 'mB': 1024*1024, 'gB': 1024*1024*1024}

    if sys.platform != 'win32':
        _proc_status = '/proc/%d/status' % os.getpid()

        def _VmB(VmKey):
            try:
                t = open(_proc_status)
                v = t.read()
                t.close()
            except IOError:
                return 0.0
            i = v.index(VmKey)
            v = v[i:].split(None, 3)
            if len(v) < 3:
                return 0.0
            return float(v[1]) * _scale[v[2]]

        memory   = lambda since=0.0: _VmB('VmSize:') - since
        resident = lambda since=0.0: _VmB('VmRSS:')  - since
        stacksize= lambda since=0.0: _VmB('VmStk:')  - since

        msg = ('Resource usage: memory=%.1fMB resident=%.1fMB stacksize=%.1fMB'
               % (memory() / _scale['MB'],
                  resident() / _scale['MB'],
                  stacksize() / _scale['MB']))
        log(msg, level)
        timingInfo('sys_platform, ' + sys.platform)
        timingInfo(prefix + 'memory, '    + str(memory()    / _scale['MB']))
        timingInfo(prefix + 'resident, '  + str(resident()  / _scale['MB']))
        timingInfo(prefix + 'stacksize, ' + str(stacksize() / _scale['MB']))
    else:
        try:
            import ctypes
            import winreg
        except ImportError:
            log('Windows resource usage not available', level)
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
                        ('ullAvailExtendedVirtual', c_ulonglong)]

        memoryStatusEx = MEMORYSTATUSEX()
        memoryStatusEx.dwLength = ctypes.sizeof(MEMORYSTATUSEX)
        kernel32.GlobalMemoryStatusEx(ctypes.byref(memoryStatusEx))

        msg = ('Resource usage: total memory=%.1fMB free memory=%.1fMB'
               % (memoryStatusEx.ullTotalPhys / _scale['MB'],
                  memoryStatusEx.ullAvailPhys / _scale['MB']))
        log(msg, level)
        timingInfo('sys_platform, ' + sys.platform)
        timingInfo(prefix + 'total_memory, ' + str(memoryStatusEx.ullTotalPhys / _scale['MB']))
        timingInfo(prefix + 'free_memory, '  + str(memoryStatusEx.ullAvailPhys / _scale['MB']))


################################################################################
if __name__ == '__main__':
    set_logfile('/tmp/anuga_test.log')
    critical('#' * 80)
    warning('Test of logging...')
    info('An info message')
    debug('A debug message (file only if console level is INFO)')
    print('This print() goes to both terminal and /tmp/anuga_test.log')
