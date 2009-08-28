
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
import pypar

#from Numeric import allclose, array, zeros, Float, take, nonzero

import numpy as num

from anuga.pmesh.mesh_interface import create_mesh_from_regions

from anuga.interface import rectangular_cross
from anuga.abstract_2d_finite_volumes.pmesh2domain import pmesh_to_domain_instance

from anuga.utilities.numerical_tools import ensure_numeric
from anuga.utilities.util_ext        import double_precision
from anuga.utilities.norms           import l1_norm, l2_norm, linf_norm

from anuga.interface import Domain
from anuga.interface import Reflective_boundary
from anuga.interface import Dirichlet_boundary
from anuga.interface import Time_boundary
from anuga.interface import Transmissive_boundary


from anuga_parallel.parallel_api import distribute, myid, numprocs


#--------------------------------------------------------------------------
# Setup parameters
#--------------------------------------------------------------------------

mesh_filename = "merimbula_10785_1.tsh"
#mesh_filename = "test-100.tsh"
yieldstep = 1
finaltime = 20
quantity = 'stage'
nprocs = 4

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
def evolution_test(parallel=False):


    domain = pmesh_to_domain_instance(mesh_filename, Domain)
    domain.set_quantity('stage', Set_Stage(756000.0, 756500.0, 2.0))

    #--------------------------------------------------------------------------
    # Create parallel domain if requested
    #--------------------------------------------------------------------------

    if parallel:
        if myid == 0: print 'DISTRIBUTING PARALLEL DOMAIN'
        domain = distribute(domain, verbose=False)

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
        if myid == 0: print 'PARALLEL EVOLVE'
    else:
        print 'SEQUENTIAL EVOLVE'
        
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
                domain.write_time()

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
            domain.write_time()
            l1list.append(l1norm)                
            l2list.append(l2norm)
            linflist.append(linfnorm)
            

    return (l1list, l2list, linflist)

# Test an 8-way run of the shallow water equations
# against the sequential code.

class Test_Parallel_Sw(unittest.TestCase):
    def testParallelSw(self):
        print "Expect this test to fail if not run from the parallel directory."
        result = os.system("mpirun -np %d python test_parallel_sw.py" % nprocs)
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
        suite = unittest.makeSuite(Test_Parallel_Sw, 'test')
        runner.run(suite)
    else:

        pypar.barrier()
        if myid == 0:
            print 'SEQUENTIAL START'
            l1norm_seq, l2norm_seq, linfnorm_seq = evolution_test(parallel=False)

        pypar.barrier()
        if myid ==0:
            print 'PARALLEL START'
        
        l1norm_par, l2norm_par, linfnorm_par = evolution_test(parallel=True)
        
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
                
            print 'Parallel test OK'


