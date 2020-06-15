
"""Test a run of the sequential shallow water domain against
a run of the parallel shallow water domain.

WARNING: This assumes that the command to run jobs is mpiexec.
Tested with MPICH and LAM (Ole)
"""
from __future__ import print_function

#------------------------------------------------------------------------------
# Import necessary modules
#------------------------------------------------------------------------------

from builtins import object
from future.utils import raise_
import unittest
import os
import sys



import numpy as num

from anuga.utilities import parallel_abstraction as pypar

#------------------------------------------
# anuga imports
#------------------------------------------
import anuga 

from anuga.utilities.numerical_tools import ensure_numeric
from anuga.utilities.util_ext        import double_precision
from anuga.utilities.norms           import l1_norm, l2_norm, linf_norm
from anuga.utilities.system_tools    import get_pathname_from_package

from anuga import Domain
from anuga import Reflective_boundary
from anuga import Dirichlet_boundary
from anuga import Time_boundary
from anuga import Transmissive_boundary

from anuga import rectangular_cross
from anuga import create_domain_from_file


from anuga import distribute, myid, numprocs, finalize


#--------------------------------------------------------------------------
# Setup parameters
#--------------------------------------------------------------------------

mod_path = get_pathname_from_package('anuga.parallel')

mesh_filename = os.path.join(mod_path,'data','merimbula_10785_1.tsh')
#mesh_filename = os.path.join('..','data','test-100.tsh')
yieldstep = 1
finaltime = 1
quantity = 'stage'
nprocs = 3
verbose = False

#--------------------------------------------------------------------------
# Setup procedures
#--------------------------------------------------------------------------
class Set_Stage(object):
    """Set an initial condition with constant water height, for x<x0
    """

    def __init__(self, x0=0.25, x1=0.5, h=1.0):
        self.x0 = x0
        self.x1 = x1
        self.h  = h

    def __call__(self, x, y):
        return self.h*((x>self.x0)&(x<self.x1))

#--------------------------------------------------------------------------
# Setup test
#--------------------------------------------------------------------------
def run_simulation(parallel=False):


    domain = create_domain_from_file(mesh_filename)
    domain.set_quantity('stage', Set_Stage(756000.0, 756500.0, 2.0))

    #--------------------------------------------------------------------------
    # Create parallel domain if requested
    #--------------------------------------------------------------------------

    if parallel:
        if myid == 0 and verbose: print('DISTRIBUTING PARALLEL DOMAIN')
        domain = distribute(domain)

    #------------------------------------------------------------------------------
    # Setup boundary conditions
    # This must currently happen *after* domain has been distributed
    #------------------------------------------------------------------------------
    domain.store = False
    Br = Reflective_boundary(domain)      # Solid reflective wall

    domain.set_boundary({'outflow' :Br, 'inflow' :Br, 'inner' :Br, 'exterior' :Br, 'open' :Br})

    #------------------------------------------------------------------------------
    # Evolution
    #------------------------------------------------------------------------------
    if parallel:
        if myid == 0 and verbose: print('PARALLEL EVOLVE')
    else:
        if verbose: print('SEQUENTIAL EVOLVE')

    for t in domain.evolve(yieldstep = yieldstep, finaltime = finaltime):
        pass

    #domain.dump_triangulation()

# Test an nprocs-way run of the shallow water equations
# against the sequential code.

class Test_parallel_shallow_domain(unittest.TestCase):
    def test_parallel_shallow_domain(self):
        #print "Expect this test to fail if not run from the parallel directory."
        
        cmd = anuga.mpicmd(os.path.abspath(__file__))
        result = os.system(cmd)
        
        assert_(result == 0)


# Because we are doing assertions outside of the TestCase class
# the PyUnit defined assert_ function can't be used.
def assert_(condition, msg="Assertion Failed"):
    if condition == False:
        #pypar.finalize()
        raise_(AssertionError, msg)

if __name__=="__main__":
    if numprocs == 1: 
        runner = unittest.TextTestRunner()
        suite = unittest.makeSuite(Test_parallel_shallow_domain, 'test')
        runner.run(suite)
    else:

        from anuga.utilities.parallel_abstraction import global_except_hook
        import sys
        sys.excepthook = global_except_hook

        pypar.barrier()
        if myid ==0:
            if verbose: print('PARALLEL START')

        run_simulation(parallel=True)
        
        if myid == 0:     
            if verbose: print('Parallel test OK')



    finalize()
