"""Test a run of the sequential shallow water domain against
a run of the parallel shallow water domain.

WARNING: This assumes that the command to run jobs is mpiexec.
Tested with MPICH and LAM (Ole)
"""

# ------------------------
# Import necessary modules
# ------------------------
import platform
import re
import unittest
import numpy as num
import os
import subprocess

verbose = False

path = os.path.dirname(__file__)  # Get folder where this script lives
run_filename = os.path.join(path, 'run_parallel_distribute_mesh.py')



class Test_parallel_distribute_mesh(unittest.TestCase):
    def setUp(self):

        # --------------------
        # Calculate extra_options
        # --------------------
        extra_options = '--oversubscribe'
        cmd = 'mpiexec -np 3 ' + extra_options + ' echo '

        # See if we get an error with --oversubscribe option
        result = subprocess.run(cmd.split(), capture_output=True)
        if result.returncode != 0:
            extra_options = ' '

        import platform
        if platform.system() == 'Windows':
            extra_options = ' '

        # --------------------
        # Run in parallel on 3 processes
        # --------------------

        cmd = 'mpiexec -np 3 ' + extra_options + ' python ' + run_filename
        if verbose:
            print(cmd)

        result = subprocess.run(cmd.split(), capture_output=True)
        if result.returncode != 0:
            print(result.stdout)
            print(result.stderr)
            raise Exception(result.stderr)

        if verbose:
            lines = result.stdout.decode('UTF-8')
            print(lines)

    def tearDown(self):
        pass

    def test_that_sequential_and_parallel_outputs_are_identical(self):

        pass
           
if __name__ == "__main__":
    runner = unittest.TextTestRunner()
    suite = unittest.makeSuite(Test_parallel_distribute_mesh, 'test')
    runner.run(suite)
