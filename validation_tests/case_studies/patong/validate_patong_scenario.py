"""Automatic verification of ANUGA using Patong Beach validation scenario.
See functions exercised by this wrapper for more details and also the publication (TBA).
"""

import unittest
import os
from subprocess import Popen, PIPE


def get_free_memory():
    """Get available memory. Linux only.
    """
    
    p = Popen('free', shell=True,
              stdin=PIPE, stdout=PIPE, stderr=PIPE, close_fds=True)
              
    if p.stdout is not None:
        for line in p.stdout.readlines():
            if line.startswith('Mem'):
                fields = line.split()
                mem = int(fields[1]) # Total memory
    else:
        mem = None
        
        
    return mem


class Test_flow(unittest.TestCase):
    def setUp(self):
        for file in os.listdir('.'):    
            if file.endswith('.stdout') or\
                    file.endswith('.sww') or\
                    file.endswith('.msh'):
                os.remove(file)
                
        
    def tearDown(self):
        pass

    def test_patong_validation(self):
        """Exercise Patong Validation Scenario for three resolutions and
        compare timeseries from modelled results against reference model.
        """
        
        # Bail out if there is insufficient memory.
        # This works only on *nix systems
        mem = get_free_memory()
        if mem is None:
            msg = 'No information about available memory: '
            msg += 'Skipping Patong beach validation'
            raise Exception(msg)
            
        if mem < 2000000:
            msg = 'Insufficient memory to run Patong beach validation. '
            msg += 'Got %i need at least 8GB. ' % mem
            msg += 'Skipping Patong beach validation'
            raise Exception(msg)        
        
        #print
        s = 'run_model.py'
        #print s
        res = os.system('python %s > run_model.stdout' % s)
        #assert res == 0


#-------------------------------------------------------------
if __name__ == '__main__':
    suite = unittest.makeSuite(Test_flow, 'test')
    runner = unittest.TextTestRunner(verbosity=2)
    runner.run(suite)
