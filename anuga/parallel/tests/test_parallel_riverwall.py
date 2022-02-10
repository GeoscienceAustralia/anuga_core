"""
Test parallel and sequential results of riverwall procedure
"""
#------------------------------------------------------------------------------
# Import necessary modules
#------------------------------------------------------------------------------
from future.utils import raise_
import unittest

verbose = False

# Test an nprocs-way run of the shallow water equations
# against the sequential code.

class Test_parallel_riverwall(unittest.TestCase):
    def test_parallel_riverwall(self):
        
        import os.path

        run_filename = os.path.abspath(__file__).replace('test_parallel', 'run_parallel')

        extra_options = '--oversubscribe'

        import platform
        if platform.system() == 'Windows':
            extra_options = ' '

        cmd = 'mpiexec -np 3 ' + extra_options + ' python ' + run_filename
        if verbose : print(cmd)
        
        import os
        result = os.system(cmd)
        assert_(result == 0)

# Because we are doing assertions outside of the TestCase class
# the PyUnit defined assert_ function can't be used.
def assert_(condition, msg="Assertion Failed"):
    if condition == False:
        #pypar.finalize()
        raise_(AssertionError, msg)

if __name__=="__main__":
    runner = unittest.TextTestRunner()
    suite = unittest.makeSuite(Test_parallel_riverwall, 'test')
    runner.run(suite)
