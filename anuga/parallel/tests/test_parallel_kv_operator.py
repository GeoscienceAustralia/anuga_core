"""
Test that the kinematic viscosity operator gives the same xvelocity field
when run serially and when run in parallel across 3 processes.
"""

import os
import sys
import subprocess
import unittest
import numpy as num

import pytest

try:
    import mpi4py
    MPI4PY_AVAILABLE = True
except ImportError:
    MPI4PY_AVAILABLE = False

path = os.path.dirname(os.path.abspath(__file__))
run_filename = os.path.join(path, 'run_parallel_kv_operator.py')

serial_file   = os.path.join(path, 'kv_serial_xvel.npy')
parallel_file = os.path.join(path, 'kv_parallel_xvel.npy')


def _remove_if_exists(f):
    if os.path.exists(f):
        os.remove(f)


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
class Test_parallel_kv_operator(unittest.TestCase):

    def setUp(self):
        _remove_if_exists(serial_file)
        _remove_if_exists(parallel_file)

        # --- serial run ---
        result = subprocess.run(
            [sys.executable, run_filename],
            capture_output=True)
        if result.returncode != 0:
            raise RuntimeError('Serial run failed:\n' +
                               result.stdout.decode() +
                               result.stderr.decode())

        extra = _extra_mpiexec_options()
        cmd = ['mpiexec', '-np', '3'] + (extra.split() if extra else []) + \
              [sys.executable, run_filename]

        # --- parallel run ---
        result = subprocess.run(cmd, capture_output=True)
        if result.returncode != 0:
            raise RuntimeError('Parallel run failed:\n' +
                               result.stdout.decode() +
                               result.stderr.decode())

    def tearDown(self):
        _remove_if_exists(serial_file)
        _remove_if_exists(parallel_file)

    def test_xvelocity_matches_serial(self):
        """Parallel xvelocity centroids match serial to within CG tolerance."""
        assert os.path.exists(serial_file),   'serial output file missing'
        assert os.path.exists(parallel_file), 'parallel output file missing'

        serial   = num.load(serial_file)
        parallel = num.load(parallel_file)

        self.assertEqual(serial.shape, parallel.shape,
                         'Serial and parallel result arrays have different shapes')

        # Allow tolerance commensurate with CG convergence (tol = 1e-5)
        self.assertTrue(
            num.allclose(serial, parallel, rtol=1e-3, atol=1e-5),
            f'Max absolute difference: {num.max(num.abs(serial - parallel)):.3e}')


if __name__ == '__main__':
    unittest.main()
