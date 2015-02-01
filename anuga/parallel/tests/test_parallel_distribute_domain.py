
"""Test a run of the sequential shallow water domain against
a run of the parallel shallow water domain.

WARNING: This assumes that the command to run jobs is mpirun.
Tested with MPICH and LAM (Ole)
"""

#------------------------------------------------------------------------------
# Import necessary modules
#------------------------------------------------------------------------------

import unittest
import os
import sys

from anuga.utilities.system_tools import get_pathname_from_package

import numpy as num

#------------------------------------------
# Import pypar without the initial output
#------------------------------------------
class NullStream:
    def write(self,text):
        pass
sys.stdout = NullStream()
import pypar
sys.stdout = sys.__stdout__


#------------------------------------------
# anuga imports
#------------------------------------------

from anuga.utilities.numerical_tools import ensure_numeric
from anuga.utilities.util_ext        import double_precision
from anuga.utilities.norms           import l1_norm, l2_norm, linf_norm

from anuga import Domain
from anuga import Reflective_boundary
from anuga import Dirichlet_boundary
from anuga import Time_boundary
from anuga import Transmissive_boundary

from anuga import rectangular_cross
from anuga import create_domain_from_file


from anuga.parallel import distribute, myid, numprocs, finalize


#--------------------------------------------------------------------------
# Setup parameters
#--------------------------------------------------------------------------

mod_path = get_pathname_from_package('anuga.parallel')

mesh_filename = os.path.join(mod_path,'data','merimbula_10785_1.tsh')
#mesh_filename = os.path.join(mod_path,'data','test-100.tsh')
yieldstep = 1
finaltime = 20
quantity = 'stage'
nprocs = 4
verbose = False

#--------------------------------------------------------------------------
# Setup procedures
#--------------------------------------------------------------------------
class Set_Stage:
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
        if myid == 0 and verbose: print 'DISTRIBUTING PARALLEL DOMAIN'
        domain = distribute(domain)

    #------------------------------------------------------------------------------
    # Setup boundary conditions
    # This must currently happen *after* domain has been distributed
    #------------------------------------------------------------------------------
    domain.store = False
    Br = Reflective_boundary(domain)      # Solid reflective wall

    domain.set_boundary({'outflow' :Br, 'inflow' :Br, 'inner' :Br, 'exterior' :Br, 'open' :Br})

    #------------------------------------------------------------------------------
    # Setup diagnostic arrays
    #------------------------------------------------------------------------------
    l1list = []
    l2list = []
    linflist = []
    l1norm = num.zeros(3, num.float)
    l2norm = num.zeros(3, num.float)
    linfnorm = num.zeros(3, num.float)
    recv_norm = num.zeros(3, num.float)

    #------------------------------------------------------------------------------
    # Evolution
    #------------------------------------------------------------------------------
    if parallel:
        if myid == 0 and verbose: print 'PARALLEL EVOLVE'
    else:
        if verbose: print 'SEQUENTIAL EVOLVE'
        
    for t in domain.evolve(yieldstep = yieldstep, finaltime = finaltime):
        edges = domain.quantities[quantity].edge_values.take(num.flatnonzero(domain.tri_full_flag),axis=0)
        l1norm[0] = l1_norm(edges[:,0])
        l1norm[1] = l1_norm(edges[:,1])
        l1norm[2] = l1_norm(edges[:,2])
        l2norm[0] = l2_norm(edges[:,0])
        l2norm[1] = l2_norm(edges[:,1])
        l2norm[2] = l2_norm(edges[:,2])
        linfnorm[0] = linf_norm(edges[:,0])
        linfnorm[1] = linf_norm(edges[:,1])
        linfnorm[2] = linf_norm(edges[:,2])
        if parallel:
            l2norm[0] = pow(l2norm[0], 2)
            l2norm[1] = pow(l2norm[1], 2)
            l2norm[2] = pow(l2norm[2], 2)
            if myid == 0:
                #domain.write_time()

                #print edges[:,1]            
                for p in range(1, numprocs):
                    recv_norm = pypar.receive(p)
                    l1norm += recv_norm
                    recv_norm = pypar.receive(p)
                    l2norm += recv_norm
                    recv_norm = pypar.receive(p)
                    linfnorm[0] = max(linfnorm[0], recv_norm[0])
                    linfnorm[1] = max(linfnorm[1], recv_norm[1])
                    linfnorm[2] = max(linfnorm[2], recv_norm[2])

                l2norm[0] = pow(l2norm[0], 0.5)
                l2norm[1] = pow(l2norm[1], 0.5)
                l2norm[2] = pow(l2norm[2], 0.5)

                l1list.append(l1norm)                
                l2list.append(l2norm)
                linflist.append(linfnorm)                
            else:
                pypar.send(l1norm, 0)
                pypar.send(l2norm, 0)
                pypar.send(linfnorm, 0)
        else:
            #domain.write_time()
            l1list.append(l1norm)                
            l2list.append(l2norm)
            linflist.append(linfnorm)
            

    return (l1list, l2list, linflist)

# Test an nprocs-way run of the shallow water equations
# against the sequential code.

class Test_parallel_distribute_domain(unittest.TestCase):
    def test_parallel_distribute_domain(self):
        #print "Expect this test to fail if not run from the parallel directory."

        abs_script_name = os.path.abspath(__file__)
        cmd = "mpirun -np %d python %s" % (nprocs, abs_script_name)
        result = os.system(cmd)

        assert_(result == 0)


# Because we are doing assertions outside of the TestCase class
# the PyUnit defined assert_ function can't be used.
def assert_(condition, msg="Assertion Failed"):
    if condition == False:
        #pypar.finalize()
        raise AssertionError, msg

if __name__=="__main__":
    if numprocs == 1: 
        runner = unittest.TextTestRunner()
        suite = unittest.makeSuite(Test_parallel_distribute_domain, 'test')
        runner.run(suite)
    else:

        pypar.barrier()
        if myid == 0:
            if verbose: print 'SEQUENTIAL START'
            l1norm_seq, l2norm_seq, linfnorm_seq = run_simulation(parallel=False)

        pypar.barrier()
        if myid ==0:
            if verbose: print 'PARALLEL START'
        
        l1norm_par, l2norm_par, linfnorm_par = run_simulation(parallel=True)
        
        if myid == 0:
            assert_(len(l1norm_seq) == len(l1norm_par))
            assert_(len(l2norm_seq) == len(l2norm_par))
            assert_(len(linfnorm_seq) == len(linfnorm_par))
            assert_(len(l1norm_seq) == len(l2norm_seq))
            assert_(len(l2norm_seq) == len(linfnorm_seq))
            # Anything smaller than tol we consider to be 0.
            # This usualy comes out to be 10^-14 (DBL_DIG = 10^-15)
            # * 10 to cater for rounding error in computation.
            tol = pow(10, -1 * (double_precision() - 1))
            for x in range(len(l1norm_seq)):
                for y in range(3):
                    # Calculate relative difference in the norms
                    assert_(abs(l1norm_seq[x][y] - l1norm_par[x][y])/l1norm_seq[x][y] < tol)
                    assert_(abs(l2norm_seq[x][y] - l2norm_par[x][y])/l2norm_seq[x][y] < tol)
                    assert_(abs(linfnorm_seq[x][y] - linfnorm_par[x][y])/linfnorm_seq[x][y] < tol)
                    if x > 0:
                        # Verify that the quantity is being conserved across iterations.
                        assert_(abs(l1norm_seq[x][y] - l1norm_seq[x-1][y]) < tol)
                        assert_(abs(l2norm_seq[x][y] - l2norm_seq[x-1][y]) < tol)
                        assert_(abs(linfnorm_seq[x][y] - linfnorm_seq[x-1][y]) < tol)
                        assert_(abs(l1norm_par[x][y] - l1norm_par[x-1][y]) < tol)
                        assert_(abs(l2norm_par[x][y] - l2norm_par[x-1][y]) < tol)
                        assert_(abs(linfnorm_par[x][y] - linfnorm_par[x-1][y]) < tol)
                
            if verbose: print 'Parallel test OK'



    finalize()
