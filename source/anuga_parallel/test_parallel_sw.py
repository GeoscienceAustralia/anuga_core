#!/usr/bin/env python
# Test a run of the sequential shallow water domain against
# a run of the parallel shallow water domain.
# WARNING: This assumes that the command to run jobs is mpirun.
# Tested with MPICH.

#mesh_filename = "test-100.tsh"
mesh_filename= "merimbula_10785_1.tsh"
yieldstep = 1
finaltime = 90
quantity = 'stage'
nprocs = 8

import unittest
import os
import sys
import print_stats
import pypar

from Numeric import array, zeros, Float, take, nonzero
from anuga.shallow_water import Domain
from anuga.shallow_water import Reflective_boundary as sw_reflective_boundary
from anuga.shallow_water import Transmissive_boundary as sw_transmissive_boundary
from parallel_shallow_water import Parallel_Domain
from parallel_shallow_water import Reflective_boundary as par_reflective_boundary
from parallel_shallow_water import Transmissive_boundary as par_transmissive_boundary
from anuga.abstract_2d_finite_volumes.pmesh2domain\
     import pmesh_to_domain_instance

from anuga.utilities.norms import *
from anuga.utilities.util_ext import double_precision
from print_stats import print_test_stats, build_full_flag
from pmesh_divide import pmesh_divide_metis
from build_submesh import build_submesh
from build_local import build_local_mesh
from build_commun import send_submesh, rec_submesh, extract_hostmesh

class Set_Stage:
    """Set an initial condition with constant water height, for x<x0
    """

    def __init__(self, x0=0.25, x1=0.5, h=1.0):
        self.x0 = x0
        self.x1 = x1
        self.h  = h

    def __call__(self, x, y):
        return self.h*((x>self.x0)&(x<self.x1))


def parallel_test():
    myid = pypar.rank()
    numprocs = pypar.size()
    proc_name = pypar.Get_processor_name()

    rect = zeros(4,Float) # Buffer for results
    if myid == 0:
        # Partition
        domain_full = pmesh_to_domain_instance(mesh_filename, Domain)
        rect = array(domain_full.xy_extent, Float)
        
        domain_full.set_quantity('stage', Set_Stage(756000.0, 756500.0, 2.0))
        #domain_full.set_quantity('stage', Set_Stage(200.0,300.0,1.0))
        domain_full.check_integrity()
        domain_full.smooth = False
        domain_full.reduction = min
        domain_full.checkpoint = False
        domain_full.visualise = False

        nodes, triangles, boundary, triangles_per_proc, quantities = \
               pmesh_divide_metis(domain_full, numprocs)

        submesh = build_submesh(nodes, triangles, boundary,\
                                quantities, triangles_per_proc)

        for p in range(1, numprocs):
            send_submesh(submesh, triangles_per_proc, p)

        points, vertices, boundary, quantities, ghost_recv_dict, full_send_dict = \
                extract_hostmesh(submesh, triangles_per_proc)

    else:
        points, vertices, boundary, quantities, ghost_recv_dict, full_send_dict = \
                rec_submesh(0)

    pypar.broadcast(rect, 0)
    domain = Parallel_Domain(points, vertices, boundary,
                             full_send_dict = full_send_dict,
                             ghost_recv_dict = ghost_recv_dict)
    tri_full_flag = build_full_flag(domain, ghost_recv_dict)

    domain.default_order = 1
    R = par_reflective_boundary(domain)
    domain.set_boundary({'outflow' : R, 'inflow' : R, 'inner' : R, 'exterior' : R, 'open' : R, 'ghost' : None})
    domain.set_quantity('stage', quantities['stage'])
    domain.set_quantity('elevation', quantities['elevation'])
    domain.store = False

    l1list = []
    l2list = []
    linflist = []
    l1norm = zeros(3, Float)
    l2norm = zeros(3, Float)
    linfnorm = zeros(3, Float)
    recv_norm = zeros(3, Float)
    # Evolution
    for t in domain.evolve(yieldstep = yieldstep, finaltime = finaltime):
        edges = take(domain.quantities[quantity].edge_values, nonzero(tri_full_flag))
        l1norm[0] = l1_norm(edges[:, 0])
        l1norm[1] = l1_norm(edges[:, 1])
        l1norm[2] = l1_norm(edges[:, 2])
        l2norm[0] = pow(l2_norm(edges[:,0]), 2)
        l2norm[1] = pow(l2_norm(edges[:,1]), 2)
        l2norm[2] = pow(l2_norm(edges[:,2]), 2)
        linfnorm[0] = linf_norm(edges[:,0])
        linfnorm[1] = linf_norm(edges[:,1])
        linfnorm[2] = linf_norm(edges[:,2])
        if myid == 0:
            domain.write_time()
            for p in range(1, numprocs):
                pypar.receive(p, recv_norm)
                l1norm += recv_norm
                pypar.receive(p, recv_norm)
                l2norm += recv_norm
                pypar.receive(p, recv_norm)
                linfnorm[0] = max(linfnorm[0], recv_norm[0])
                linfnorm[1] = max(linfnorm[1], recv_norm[1])
                linfnorm[2] = max(linfnorm[2], recv_norm[2])
            l1list.append(l1norm)
            l2norm[0] = pow(l2norm[0], 0.5)
            l2norm[1] = pow(l2norm[1], 0.5)
            l2norm[2] = pow(l2norm[2], 0.5)
            l2list.append(l2norm)
            linflist.append(linfnorm)
        else:
            pypar.send(l1norm, 0)
            pypar.send(l2norm, 0)
            pypar.send(linfnorm, 0)
    return (l1list, l2list, linflist)

def sequential_test():
    domain_full = pmesh_to_domain_instance(mesh_filename, Domain)

    domain_full.set_quantity('stage', Set_Stage(756000.0, 756500.0, 2.0))
    domain_full.check_integrity()
    domain_full.default_order = 1
    domain_full.smooth = False
    domain_full.reduction = min
    domain_full.store = False
    domain_full.checkpoint = False
    domain_full.visualise = False
    R = sw_reflective_boundary(domain_full)
    domain_full.set_boundary({'outflow' : R, 'inflow' : R, 'inner' : R, 'exterior' : R, 'open' : R})
    l1list = []
    l2list = []
    linflist = []
    l1norm = zeros(3, Float)
    l2norm = zeros(3, Float)
    linfnorm = zeros(3, Float)
    for t in domain_full.evolve(yieldstep = yieldstep, finaltime = finaltime):
        domain_full.write_time()
        edge = domain_full.quantities[quantity].edge_values
        l1norm[0] = l1_norm(edge[:,0])
        l1norm[1] = l1_norm(edge[:,1])
        l1norm[2] = l1_norm(edge[:,2])
        l2norm[0] = l2_norm(edge[:,0])
        l2norm[1] = l2_norm(edge[:,1])
        l2norm[2] = l2_norm(edge[:,2])
        linfnorm[0] = linf_norm(edge[:,0])
        linfnorm[1] = linf_norm(edge[:,1])
        linfnorm[2] = linf_norm(edge[:,2])
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
        pypar.finalize()
        raise AssertionError, msg

if __name__=="__main__":
    if pypar.size() == 1: 
        runner = unittest.TextTestRunner()
        suite = unittest.makeSuite(Test_Parallel_Sw, 'test')
        runner.run(suite)
    else:
        if pypar.rank() == 0:
            l1norm_seq, l2norm_seq, linfnorm_seq = sequential_test()
        l1norm_par, l2norm_par, linfnorm_par = parallel_test()
        if pypar.rank() == 0:
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
        pypar.finalize()

