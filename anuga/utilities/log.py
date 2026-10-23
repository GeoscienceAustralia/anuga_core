#!/usr/bin/env python

"""
A simple logging module that logs to the console and a logfile.

Basic usage (print + log both go to terminal and file):

    import anuga.utilities.log as log

    log.set_logfile('./my.log')   # truncates my.log, activates tee to file

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

import atexit
import ctypes
import io
import os
import sys
import threading
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


################################################################################
# _StdoutTee — OS-level tee of file descriptor 1 (captures C output too)
################################################################################

# Control markers pushed through the pipe to switch the terminal side of the
# tee off and on again.  Sending them in-band (rather than just flipping a flag)
# means the switch happens at exactly the right point in the byte stream, with
# no race against output still in flight.  \x00 never occurs in real output.
#
# They nest: OFF increments a depth counter and ON decrements it, so wrapping an
# individual write (see _TeeFileStream) still behaves inside a file_only block.
_MARK_OFF = '\x00\x01ANUGA_ECHO_OFF\x01\x00'
_MARK_ON = '\x00\x01ANUGA_ECHO_ON\x01\x00'


def _force_c_stdout_line_buffered():
    """Ask the C runtime to line-buffer its stdout.

    Redirecting fd 1 into a pipe makes libc switch stdout from line-buffered
    to block-buffered (it is no longer a tty), which would hold C printf()
    output back until 4 kB had accumulated or the process exited — so the GPU
    banner would land in the log far away from the Python lines around it.

    setvbuf() is only defined before any output has been written to the stream,
    which is why set_logfile() wants calling early in a run.  Failure is not
    fatal (non-glibc platforms simply keep block buffering).
    """
    try:
        libc = ctypes.CDLL(None)
        c_stdout = ctypes.c_void_p.in_dll(libc, 'stdout')
        libc.setvbuf(c_stdout, None, 1, 0)      # 1 == _IOLBF
    except Exception:
        pass


class _LockedStream:
    """Serialise writes to the terminal between the pump thread and logging."""

    def __init__(self, stream, lock):
        self._stream = stream
        self._lock = lock

    def write(self, message):
        with self._lock:
            self._stream.write(message)
            self._stream.flush()

    def flush(self):
        with self._lock:
            self._stream.flush()

    def __getattr__(self, name):
        return getattr(self._stream, name)


class _TeeFileStream:
    """Stream for the logging file handler that writes *through* the tee.

    Having the file handler open its own handle would give the log file two
    independent writers — the handler and the tee's pump thread — so a log
    record could overtake print()/C output emitted earlier.  Pushing the record
    through the same pipe, wrapped in echo-off markers so it is not echoed to
    the terminal (the console handler already does that), keeps the file in
    exact write order with a single writer.
    """

    def __init__(self, tee):
        self._tee = tee

    def write(self, message):
        # One write() call keeps the record and its markers contiguous in the
        # pipe (atomic up to PIPE_BUF; log records are far smaller).
        os.write(1, (_MARK_OFF + message + _MARK_ON).encode('utf-8',
                                                            errors='replace'))

    def flush(self):
        pass


class _StdoutTee:
    """Tee file descriptor 1 to both the terminal and the log file.

    TeeStream only sees writes that pass through Python's sys.stdout object.
    Anything the C/Cython extensions print — the GPU domain banner, Triangle's
    mesh generation reporting, the many printf()s in the kernels — is written
    straight to file descriptor 1 and bypasses it completely, which is why such
    output used to reach the terminal but never the logfile.

    This class redirects fd 1 into a pipe and pumps everything arriving there
    to the terminal and the log file alike, so C and Python output are captured
    together and in true write order.

    Caveat: the *log file* is exactly ordered (single writer — everything,
    including the log records, goes through the pipe).  The *terminal* is not
    quite: the logging console handler writes to it synchronously while
    print()/C output takes the pump thread's detour, so a log record can appear
    just ahead of a print issued microseconds earlier.  Real runs space these
    seconds apart, so it does not show in practice.
    """

    def __init__(self, logfile_path, mode='a'):
        self._log = open(logfile_path, mode, encoding='utf-8', errors='replace')
        self._lock = threading.Lock()
        self._echo_off_depth = 0
        self._switched = threading.Event()
        self.file_stream = _TeeFileStream(self)

        sys.stdout.flush()

        # Keep the real terminal reachable after fd 1 has been taken over.
        self._saved_fd = os.dup(1)
        self._terminal = open(self._saved_fd, 'w', encoding='utf-8',
                              errors='replace', buffering=1, closefd=True)
        self.terminal_stream = _LockedStream(self._terminal, self._lock)

        read_fd, write_fd = os.pipe()
        os.dup2(write_fd, 1)
        os.close(write_fd)
        self._pipe_r = read_fd

        self._thread = threading.Thread(target=self._pump, daemon=True,
                                        name='anuga-log-tee')
        self._thread.start()

        _force_c_stdout_line_buffered()

        # Python buffers its own stdout in blocks whenever fd 1 is not a tty —
        # which it no longer is, and never was for a redirected batch run.  Left
        # alone, print() output would reach the pipe late and interleave wrongly
        # with the C output.  Line buffering keeps the log in true write order.
        self._old_line_buffering = None
        if hasattr(sys.stdout, 'reconfigure'):
            try:
                self._old_line_buffering = sys.stdout.line_buffering
                sys.stdout.reconfigure(line_buffering=True)
            except (AttributeError, ValueError, OSError):
                self._old_line_buffering = None

    # -- pump ---------------------------------------------------------------

    def _emit(self, text):
        if not text:
            return
        with self._lock:
            if self._echo_off_depth == 0:
                self._terminal.write(text)
                self._terminal.flush()
            self._log.write(text)
            self._log.flush()

    def _pump(self):
        buf = ''
        while True:
            try:
                chunk = os.read(self._pipe_r, 65536)
            except OSError:
                break
            if not chunk:
                break
            buf += chunk.decode('utf-8', errors='replace')

            # Act on any complete control markers, in stream order.
            while True:
                off, on = buf.find(_MARK_OFF), buf.find(_MARK_ON)
                if off < 0 and on < 0:
                    break
                if on < 0 or (0 <= off < on):
                    mark, step, at = _MARK_OFF, 1, off
                else:
                    mark, step, at = _MARK_ON, -1, on
                self._emit(buf[:at])
                self._echo_off_depth = max(0, self._echo_off_depth + step)
                buf = buf[at + len(mark):]
                self._switched.set()

            # Hold back a trailing fragment that may be a partial marker.
            cut = buf.rfind('\x00')
            if cut >= 0 and (_MARK_OFF.startswith(buf[cut:])
                             or _MARK_ON.startswith(buf[cut:])):
                self._emit(buf[:cut])
                buf = buf[cut:]
            else:
                self._emit(buf)
                buf = ''
        self._emit(buf)

    # -- terminal echo control ---------------------------------------------

    def _send(self, mark):
        self._switched.clear()
        os.write(1, mark.encode('utf-8'))
        # Wait for the pump to reach the marker so the switch is exact.
        self._switched.wait(timeout=1.0)

    @contextmanager
    def echo_off(self):
        """Send output to the log file only for the duration of the block."""
        sys.stdout.flush()
        self._send(_MARK_OFF)
        try:
            yield
        finally:
            sys.stdout.flush()
            self._send(_MARK_ON)

    # -- teardown -----------------------------------------------------------

    def stop(self):
        try:
            sys.stdout.flush()
        except ValueError:
            pass
        if self._old_line_buffering is not None:
            try:
                sys.stdout.reconfigure(line_buffering=self._old_line_buffering)
            except (AttributeError, ValueError, OSError):
                pass
        # Restoring fd 1 drops the last reference to the pipe's write end, so
        # the pump sees EOF and drains what is left.
        os.dup2(self._saved_fd, 1)
        self._thread.join(timeout=5.0)
        try:
            os.close(self._pipe_r)
        except OSError:
            pass
        with self._lock:
            self._terminal.flush()
            self._log.close()


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
    if _fd_tee is not None:
        # Mute the terminal side of the fd-level tee.  This also hides
        # C-level output (Triangle's mesh reporting, the GPU kernels) that
        # never passed through sys.stdout in the first place.
        with _fd_tee.echo_off():
            yield
        return

    original = sys.stdout
    if isinstance(sys.stdout, TeeStream):
        sys.stdout = _FileOnlyStream(sys.stdout._log)
    else:
        # No logfile active — discard output
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

# Active OS-level stdout tee (captures C printf as well), or None.
_fd_tee = None

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

    - *path* is truncated: each run starts a fresh log rather than
      appending to the previous run's.
    - file descriptor 1 is tee'd, so print() *and* output from the C
      extensions go to both the terminal and *path*.
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
    global _fd_tee

    # Close any existing TeeStream
    if isinstance(sys.stdout, TeeStream):
        sys.stdout.close()
        sys.stdout = sys.__stdout__

    if _fd_tee is not None:
        _fd_tee.stop()
        _fd_tee = None

    log_filename = path
    console_logging_level = console_level
    log_logging_level = file_level
    _setup = False  # force re-initialisation on next log() call

    # Start a fresh log for this run.  Truncating here, once, rather than
    # opening the writers in 'w' mode means a later handler rebuild (after
    # close_logfile(), say) reopens the file without wiping what the run has
    # already written.
    try:
        open(path, 'w').close()
    except OSError:
        pass  # unwritable path — the open() below reports it properly

    # Tee file descriptor 1 so that print() *and* output from the C extensions
    # reach both the terminal and the file.  sys.stdout is deliberately left
    # alone: it already writes to fd 1, so the tee picks it up.
    _fd_tee = _StdoutTee(path)
    atexit.register(close_logfile)

    # Trigger logging setup now
    log('Logfile opened: ' + path, INFO)


def close_logfile():
    """Stop tee-ing output, restore fd 1 and close the log file.

    Safe to call more than once.  Registered with atexit by set_logfile(), so
    an ordinary run needs no explicit call.
    """
    global _fd_tee, _setup

    if _fd_tee is None:
        return

    tee, _fd_tee = _fd_tee, None

    # Drop the console handler bound to the tee's terminal stream before the
    # stream goes away; the next log() call rebuilds the handlers.
    root = logging.getLogger('')
    for h in root.handlers[:]:
        root.removeHandler(h)
    _setup = False

    tee.stop()


def _console_stream():
    """The stream the logging console handler should write to.

    With the fd-level tee running this must be the *saved* terminal, not fd 1:
    writing log records to fd 1 would send them through the tee and duplicate
    every one of them in the file — once formatted by the file handler, and
    again as raw text.
    """
    if _fd_tee is not None:
        return _fd_tee.terminal_stream
    return sys.__stdout__


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
            if _fd_tee is not None:
                # Write through the tee so the file has a single writer and
                # stays in true write order with print()/C output.
                file_handler = logging.StreamHandler(_fd_tee.file_stream)
            else:
                file_handler = logging.FileHandler(log_filename, mode='a')
            file_handler.setLevel(log_logging_level)
            file_handler.setFormatter(logging.Formatter(fmt))

            root = logging.getLogger('')
            root.setLevel(min(log_logging_level, console_logging_level))

            # Remove any pre-existing handlers to avoid duplicates on re-init
            for h in root.handlers[:]:
                root.removeHandler(h)

            root.addHandler(file_handler)

            console = logging.StreamHandler(_console_stream())
            console.setLevel(console_logging_level)
            console.setFormatter(logging.Formatter('%(message)s'))
            root.addHandler(console)
        else:
            # No file configured: just console at console_logging_level
            root = logging.getLogger('')
            root.setLevel(console_logging_level)
            for h in root.handlers[:]:
                root.removeHandler(h)
            console = logging.StreamHandler(_console_stream())
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
            except OSError:
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
            except OSError:
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
