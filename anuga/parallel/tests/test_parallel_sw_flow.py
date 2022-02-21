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

class Test_parallel_sw_flow(unittest.TestCase):
    def setUp(self):
        # Run the sequential and parallel simulations to produce sww files for comparison.
        
        path = os.path.dirname(__file__)  # Get folder where this script lives
        run_filename = os.path.join(path, 'run_parallel_sw_flow.py')

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
        pass
        #os.remove('sw_flow_sequential.sww')
        #os.remove('sw_flow_parallel.sww')


        
    def test_parallel_sw_flow(self):

        import anuga.utilities.plot_utils as util  # Note if this is imported at the top level
                                                   # it'll interfere with running the subprocesses.       
        
        # Assuming both sequential and parallel simulations have been run, compare the 
        # merged sww files from the parallel run to the sequential output.
        if verbose: 
            print('Comparing SWW files')

        sdomain_v = util.get_output('sw_flow_sequential.sww')
        sdomain_c = util.get_centroids(sdomain_v)

        pdomain_v = util.get_output('sw_flow_parallel.sww')
        pdomain_c = util.get_centroids(pdomain_v)
        

        # Test some values against the original ordering

        if verbose:
            order = 0
            print('Centroid values')
            print(num.linalg.norm(sdomain_c.x - pdomain_c.x, ord=order))
            print(num.linalg.norm(sdomain_c.y - pdomain_c.y, ord=order))
            print(num.linalg.norm(sdomain_c.stage[-1] - pdomain_c.stage[-1], ord=order))
            print(num.linalg.norm(sdomain_c.xmom[-1] - pdomain_c.xmom[-1], ord=order))
            print(num.linalg.norm(sdomain_c.ymom[-1] - pdomain_c.ymom[-1], ord=order))
            print(num.linalg.norm(sdomain_c.xvel[-1] - pdomain_c.xvel[-1], ord=order))
            print(num.linalg.norm(sdomain_c.yvel[-1] - pdomain_c.yvel[-1], ord=order))        

        print('T, S stage,     P stage')    
        for i in range(len(sdomain_c.stage)):
            print(i, sdomain_c.stage[-1, i], pdomain_c.stage[-1, i])    

        
        msg = 'Values not identical'
        assert num.allclose(sdomain_c.stage, pdomain_c.stage), msg
        assert num.allclose(sdomain_c.xmom, pdomain_c.xmom)
        assert num.allclose(sdomain_c.ymom, pdomain_c.ymom)
        assert num.allclose(sdomain_c.xvel, pdomain_c.xvel)
        assert num.allclose(sdomain_c.yvel, pdomain_c.yvel)

        assert num.allclose(sdomain_v.x, pdomain_v.x)
        assert num.allclose(sdomain_v.y, pdomain_v.y)


if __name__ == "__main__":
    runner = unittest.TextTestRunner()
    suite = unittest.makeSuite(Test_parallel_sw_flow, 'test')
    runner.run(suite)
