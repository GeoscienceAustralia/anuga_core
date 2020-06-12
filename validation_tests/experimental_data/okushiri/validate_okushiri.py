"""Automatic verification of ANUGA flows.
See functions exercised by this wrapper for more details
"""

import unittest
import os
import numpy
import anuga

indent = anuga.indent

args = anuga.get_args()
alg = args.alg
verbose = args.verbose

class Test_results(unittest.TestCase):
    def setUp(self):
        for file in os.listdir('.'):    
            if file.endswith('.stdout') or\
                    file.endswith('.sww') or\
                    file.endswith('.msh') or\
                    file.endswith('.png'):
                os.remove(file)
                
        
    def tearDown(self):
        pass

    def test_simulation(self):
    
        if verbose:
            print(indent+'Running simulation script')

        s = 'run_okushiri.py'
        res = anuga.run_anuga_script(s, args=args)

        # Test that script runs ok
        assert res == 0


        if verbose:
            print(indent+'Running test script')

        s = 'test_results.py'
        res = os.system('python %s' %s)
        # Test that script runs ok
        assert res == 0

        


#-------------------------------------------------------------
if __name__ == '__main__':
    suite = unittest.makeSuite(Test_results, 'test')
    runner = unittest.TextTestRunner(verbosity=2)
    runner.run(suite)
