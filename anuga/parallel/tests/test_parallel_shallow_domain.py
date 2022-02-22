"""Test a run of the sequential shallow water domain against
a run of the parallel shallow water domain.

WARNING: This assumes that the command to run jobs is mpiexec.
Tested with MPICH and LAM (Ole)
"""

# ------------------------
# Import necessary modules
# ------------------------
import platform
import unittest
import numpy as num
import os
import subprocess

verbose = False

path = os.path.dirname(__file__)  # Get folder where this script lives
run_filename = os.path.join(path, 'run_parallel_shallow_domain.py')

# These must be the same as given in the run_file.
sequential_sww_file = 'shallow_water_sequential.sww'
parallel_sww_file = 'shallow_water_parallel.sww'

class Test_parallel_shallow_domain(unittest.TestCase):
    def setUp(self):
        # Run the sequential and parallel simulations to produce sww files for comparison.

        # ----------------------
        # First run sequentially
        # ----------------------
        cmd = 'python ' + run_filename
        if verbose:
            print(cmd)

        result = subprocess.run(cmd.split(), capture_output=True)
        if result.returncode != 0:
            print(result.stdout)
            print(result.stderr)
            raise Exception(result.stderr)

        # --------------------
        # Then run in parallel
        # --------------------
        if platform.system() == 'Windows':
            extra_options = ' '
        else:
            # E.g. for Ubuntu Linux
            extra_options = '--oversubscribe'

        cmd = 'mpiexec -np 3 ' + extra_options + ' python ' + run_filename
        if verbose:
            print(cmd)

        result = subprocess.run(cmd.split(), capture_output=True)
        if result.returncode != 0:
            print(result.stdout)
            print(result.stderr)
            raise Exception(result.stderr)

    def tearDown(self):
        os.remove(sequential_sww_file)
        os.remove(parallel_sww_file)

    def test_that_sequential_and_parallel_outputs_are_identical(self):
        from anuga.file.sww import sww_files_are_equal
        assert sww_files_are_equal(sequential_sww_file, parallel_sww_file) 

if __name__ == "__main__":
    runner = unittest.TextTestRunner()
    suite = unittest.makeSuite(Test_parallel_shallow_domain, 'test')
    runner.run(suite)
