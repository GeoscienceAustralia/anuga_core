"""
Pytest driver: launches run_parallel_kv_unit_tests.py serially and under
mpiexec, checks exit codes.
"""

import os
import sys
import subprocess
import unittest

import pytest

try:
    import mpi4py
    MPI4PY_AVAILABLE = True
except ImportError:
    MPI4PY_AVAILABLE = False

path = os.path.dirname(os.path.abspath(__file__))
run_filename = os.path.join(path, 'run_parallel_kv_unit_tests.py')


def _extra_mpiexec_options():
    """Return '--oversubscribe' if mpiexec supports it, else ''."""
    import platform
    if platform.system() == 'Windows':
        return ''
    result = subprocess.run(
        ['mpiexec', '-np', '2', '--oversubscribe', 'echo', 'ok'],
        capture_output=True)
    return '--oversubscribe' if result.returncode == 0 else ''


@pytest.mark.slow
@pytest.mark.skipif(not MPI4PY_AVAILABLE, reason='requires mpi4py')
class Test_parallel_kv_unit_tests(unittest.TestCase):

    def _run(self, cmd):
        result = subprocess.run(cmd, capture_output=True)
        if result.returncode != 0:
            raise AssertionError(
                f'Command failed: {" ".join(cmd)}\n'
                + result.stdout.decode()
                + result.stderr.decode())

    def test_serial_primitives(self):
        """All KV parallel primitives pass in serial (numprocs == 1)."""
        self._run([sys.executable, run_filename])

    def test_parallel_primitives_3proc(self):
        """All KV parallel primitives pass under mpiexec -np 3."""
        extra = _extra_mpiexec_options()
        cmd = (['mpiexec', '-np', '3']
               + (extra.split() if extra else [])
               + [sys.executable, run_filename])
        self._run(cmd)


if __name__ == '__main__':
    unittest.main()
