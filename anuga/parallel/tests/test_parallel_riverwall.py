"""
Test that parallel and sequential results of riverwall simulation are identical
"""

#------------------------------------------------------------------------------
# Import necessary modules
#------------------------------------------------------------------------------
import platform
import unittest
import numpy as num
import os
import subprocess

verbose = False

class Test_parallel_riverwall(unittest.TestCase):
    def setUp(self):
        # Run the sequential and parallel simulations to produce sww files for comparison.
        
        path = os.path.dirname(__file__)  # Get folder where this script lives
        run_filename = os.path.join(path, 'run_parallel_riverwall.py')

        #-----------------------        
        # First run sequentially
        #-----------------------
        cmd = 'python ' + run_filename
        if verbose: 
            print(cmd)        
            
        result = subprocess.run(cmd.split(), capture_output=True)
        if result.returncode != 0:
            print(result.stdout)
            print(result.stderr)
            raise Exception(result.stderr)        

        #---------------------
        # Then run in parallel
        #---------------------
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
        os.remove('s_riverwall.sww')
        os.remove('p_riverwall.sww')
        os.remove('runup.msh')

        
    def test_parallel_riverwall(self):

        import anuga.utilities.plot_utils as util # Note if this is imported at the top level
                                                  # it'll interfere with running the subprocesses.       
        
        # Assuming both sequential and parallel simulations have been run, compare the 
        # merged sww files from the parallel run to the sequential output.
        if verbose: 
            print('Comparing SWW files')

        sdomain_v = util.get_output('s_riverwall.sww')
        sdomain_c = util.get_centroids(sdomain_v)

        pdomain_v = util.get_output('p_riverwall.sww')
        pdomain_c = util.get_centroids(pdomain_v)
        

        # Test some values against the original ordering

        if verbose:
            order = 0
            print('Centroid values')
            print(num.linalg.norm(sdomain_c.x-pdomain_c.x, ord=order))
            print(num.linalg.norm(sdomain_c.y-pdomain_c.y, ord=order))
            print(num.linalg.norm(sdomain_c.stage[-1] - pdomain_c.stage[-1], ord=order))
            print(num.linalg.norm(sdomain_c.xmom[-1] - pdomain_c.xmom[-1], ord=order))
            print(num.linalg.norm(sdomain_c.ymom[-1] - pdomain_c.ymom[-1], ord=order))
            print(num.linalg.norm(sdomain_c.xvel[-1] - pdomain_c.xvel[-1], ord=order))
            print(num.linalg.norm(sdomain_c.yvel[-1] - pdomain_c.yvel[-1], ord=order))        

        assert num.allclose(sdomain_c.stage, pdomain_c.stage)
        assert num.allclose(sdomain_c.xmom, pdomain_c.xmom)
        assert num.allclose(sdomain_c.ymom, pdomain_c.ymom)
        assert num.allclose(sdomain_c.xvel, pdomain_c.xvel)
        assert num.allclose(sdomain_c.yvel, pdomain_c.yvel)

        assert num.allclose(sdomain_v.x, pdomain_v.x)
        assert num.allclose(sdomain_v.y, pdomain_v.y)


if __name__ == "__main__":
    runner = unittest.TextTestRunner()
    suite = unittest.makeSuite(Test_parallel_riverwall, 'test')
    runner.run(suite)
