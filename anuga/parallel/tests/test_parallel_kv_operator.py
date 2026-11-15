"""
Test that the kinematic viscosity operator gives the same xvelocity field
when run serially and when run in parallel across 3 processes, in its
velocity and momentum forms and with the explicit scheme.
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

FORMS = ('velocity', 'momentum', 'explicit')


def serial_file(form):
    return os.path.join(path, 'kv_serial_xvel_%s.npy' % form)


def parallel_file(form):
    return os.path.join(path, 'kv_parallel_xvel_%s.npy' % form)


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
        extra = _extra_mpiexec_options()
        for form in FORMS:
            _remove_if_exists(serial_file(form))
            _remove_if_exists(parallel_file(form))

            # --- serial run ---
            result = subprocess.run(
                [sys.executable, run_filename, form],
                capture_output=True)
            if result.returncode != 0:
                raise RuntimeError('Serial run failed:\n' +
                                   result.stdout.decode() +
                                   result.stderr.decode())

            cmd = ['mpiexec', '-np', '3'] + (extra.split() if extra else []) + \
                  [sys.executable, run_filename, form]

            # --- parallel run ---
            result = subprocess.run(cmd, capture_output=True)
            if result.returncode != 0:
                raise RuntimeError('Parallel run failed:\n' +
                                   result.stdout.decode() +
                                   result.stderr.decode())

    def tearDown(self):
        for form in FORMS:
            _remove_if_exists(serial_file(form))
            _remove_if_exists(parallel_file(form))

    def test_xvelocity_matches_serial(self):
        """Parallel xvelocity centroids match serial to within CG tolerance."""
        for form in FORMS:
            assert os.path.exists(serial_file(form)),   'serial output file missing'
            assert os.path.exists(parallel_file(form)), 'parallel output file missing'

            serial   = num.load(serial_file(form))
            parallel = num.load(parallel_file(form))

            self.assertEqual(serial.shape, parallel.shape,
                             'Serial and parallel result arrays have different shapes')

            # Allow tolerance commensurate with CG convergence (tol = 1e-5);
            # the explicit scheme has no solve, so it should match to round-off
            rtol, atol = (0.0, 1e-10) if form == 'explicit' else (1e-3, 1e-5)
            self.assertTrue(
                num.allclose(serial, parallel, rtol=rtol, atol=atol),
                f'{form}: max absolute difference: '
                f'{num.max(num.abs(serial - parallel)):.3e}')


if __name__ == '__main__':
    unittest.main()
