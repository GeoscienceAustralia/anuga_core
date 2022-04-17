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
run_filename = os.path.join(path, 'run_parallel_distribute_domain.py')

# These must be the same as given in the run_file.
sequential_file = 'distribute_domain_sequential.txt'
parallel_file = 'distribute_domain_parallel.txt'

class Test_parallel_distribute_domain(unittest.TestCase):
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
        os.remove(sequential_file)
        os.remove(parallel_file)

    def test_that_sequential_and_parallel_outputs_are_identical(self):
        fid_seq = open(sequential_file)
        fid_par = open(parallel_file)
        
        seq_values = fid_seq.readlines()        
        par_values = fid_par.readlines()      
        
        N = len(seq_values)
        assert len(par_values) == N
        
        for i in range(N):
            seq_nums = [float(x) for x in seq_values[i].split()]          
            par_nums = [float(x) for x in par_values[i].split()]     
            assert num.allclose(seq_nums, par_nums)
           
if __name__ == "__main__":
    runner = unittest.TextTestRunner()
    suite = unittest.makeSuite(Test_parallel_distribute_domain, 'test')
    runner.run(suite)
