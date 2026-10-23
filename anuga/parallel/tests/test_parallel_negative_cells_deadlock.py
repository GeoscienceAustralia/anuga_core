"""Regression test: update_conserved_quantities() must not deadlock in parallel
when negative-depth cells are clamped on some ranks but not others.

The "possible loss of conservation" warning once computed the total volume via
the collective get_water_volume() inside a per-rank `if num_negative_ids > 0`
branch, so ranks that clamped cells entered an MPI_Allreduce that the other
ranks never joined -> deadlock (seen as `mpirun -np 2 run_small_towradgi.py`
hanging at the first sub-timestep while `-np 1` worked).

This launches run_parallel_negative_cells_deadlock.py under mpiexec with a
bounded timeout. A regression re-introduces the hang, the subprocess times out,
and the test fails; the fixed code returns promptly and prints the sentinel.
"""

import os
import subprocess
import sys
import unittest

try:
    import mpi4py  # noqa: F401
except ImportError:
    pass

import pytest

path = os.path.dirname(__file__)
run_filename = os.path.join(path, 'run_parallel_negative_cells_deadlock.py')

# Generous relative to the ~seconds the run needs, tight enough that a genuine
# deadlock is caught quickly rather than hanging CI.
TIMEOUT_SECONDS = 120


@pytest.mark.skipif('mpi4py' not in sys.modules,
                    reason="requires the mpi4py module")
class Test_parallel_negative_cells_deadlock(unittest.TestCase):

    def test_no_deadlock_with_divergent_negative_cells(self):
        # --oversubscribe where the MPI supports it (matches the other parallel
        # tests); harmless to omit otherwise.
        extra = '--oversubscribe'
        probe = subprocess.run(('mpiexec -np 2 ' + extra + ' echo').split(),
                               capture_output=True)
        if probe.returncode != 0:
            extra = ''

        cmd = ('mpiexec -np 2 ' + extra + ' python ' + run_filename).split()
        try:
            result = subprocess.run(cmd, capture_output=True, text=True,
                                    timeout=TIMEOUT_SECONDS)
        except subprocess.TimeoutExpired:
            self.fail(
                'update_conserved_quantities() deadlocked in parallel with '
                'negative-depth cells on only one rank (timed out after '
                '%ds) — a collective was called under a per-rank condition.'
                % TIMEOUT_SECONDS)

        if result.returncode != 0:
            print(result.stdout)
            print(result.stderr)
            self.fail('parallel runner exited %d' % result.returncode)

        self.assertIn('NEGATIVE_CELLS_DEADLOCK_OK', result.stdout)


if __name__ == '__main__':
    runner = unittest.TextTestRunner()
    suite = unittest.TestLoader().loadTestsFromTestCase(
        Test_parallel_negative_cells_deadlock)
    runner.run(suite)
