"""
Test that parallel and sequential results of riverwall simulation are identical
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
run_filename = os.path.join(path, 'run_parallel_sw_flow.py')

# These must be the same as given in the run_file.
sequential_sww_file = 'sw_flow_sequential.sww'
parallel_sww_file = 'sw_flow_parallel.sww'

class Test_parallel_sw_flow(unittest.TestCase):
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
            extra_options = ' '

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
    suite = unittest.makeSuite(Test_parallel_sw_flow, 'test')
    runner.run(suite)
