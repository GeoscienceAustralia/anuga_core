"""
Simple water flow example using ANUGA

Water driven up a linear slope and time varying boundary,
similar to a beach environment

This is a very simple test of the parallel algorithm using the simplified parallel API
"""
from __future__ import print_function
from __future__ import division


#------------------------------------------------------------------------------
# Import necessary modules
#------------------------------------------------------------------------------
from builtins import range
from past.utils import old_div
from future.utils import raise_
import unittest
import os
import sys
#import pypar
import numpy as num

import anuga

from anuga import Domain
from anuga import Reflective_boundary
from anuga import Dirichlet_boundary
from anuga import Time_boundary
from anuga import Transmissive_boundary
from anuga import rectangular_cross_domain

from anuga import distribute, myid, numprocs, send, receive, barrier, finalize

#--------------------------------------------------------------------------
# Setup parameters
#--------------------------------------------------------------------------
yieldstep = 0.25
finaltime = 1.0
nprocs = 4
N = 29
M = 29
verbose = False

#---------------------------------
# Setup Functions
#---------------------------------
def topography(x,y):
    return old_div(-x,2)

###########################################################################
# Setup Test
##########################################################################
def run_simulation(parallel=False, G = None, seq_interpolation_points=None, verbose=False):

    #--------------------------------------------------------------------------
    # Setup computational domain and quantities
    #--------------------------------------------------------------------------
    domain = rectangular_cross_domain(M, N)


    domain.set_quantity('elevation', topography) # Use function for elevation
    domain.set_quantity('friction', 0.0)         # Constant friction
    domain.set_quantity('stage', expression='elevation') # Dry initial stage

    domain.set_low_froude(0)

    domain.set_name('runup')                    # Set sww filename
    domain.set_datadir('.')                     # Set output dir

    #--------------------------------------------------------------------------
    # Create the parallel domain
    #--------------------------------------------------------------------------
    if parallel:
        if myid == 0 and verbose : print('DISTRIBUTING PARALLEL DOMAIN')
        domain = distribute(domain, verbose=False)

    #--------------------------------------------------------------------------
    # Setup domain parameters
    #--------------------------------------------------------------------------


    domain.set_quantities_to_be_stored(None)


    #------------------------------------------------------------------------------
    # Setup boundary conditions
    # This must currently happen *AFTER* domain has been distributed
    #------------------------------------------------------------------------------

    Br = Reflective_boundary(domain)      # Solid reflective wall
    Bd = Dirichlet_boundary([-0.2,0.,0.]) # Constant boundary values

    # Associate boundary tags with boundary objects
    domain.set_boundary({'left': Br, 'right': Bd, 'top': Br, 'bottom': Br})

    #------------------------------------------------------------------------------
    # Find which sub_domain in which the interpolation points are located
    #
    # Sometimes the interpolation points sit exactly
    # between two centroids, so in the parallel run we
    # reset the interpolation points to the centroids
    # found in the sequential run
    #------------------------------------------------------------------------------
    interpolation_points = [[0.4,0.5], [0.6,0.5], [0.8,0.5], [0.9,0.5]]


    gauge_values = []
    tri_ids = []
    for i, point in enumerate(interpolation_points):
        gauge_values.append([]) # Empty list for timeseries

        #if is_inside_polygon(point, domain.get_boundary_polygon()):
        #print "Point ", myid, i, point
        try:
            k = domain.get_triangle_containing_point(point)
            if domain.tri_full_flag[k] == 1:
                tri_ids.append(k)
            else:
                tri_ids.append(-1)
        except:
            tri_ids.append(-2)

        #print "  tri_ids ",myid, i, tri_ids[-1]

    if verbose: print('P%d has points = %s' %(myid, tri_ids))


    c_coord = domain.get_centroid_coordinates()
    interpolation_points = []
    for id in tri_ids:
        if id<1:
            if verbose: print('WARNING: Interpolation point not within the domain!')
        interpolation_points.append(c_coord[id,:])

    #------------------------------------------------------------------------------
    # Evolve system through time
    #------------------------------------------------------------------------------
    time = []

    if parallel:
        if myid == 0 and verbose: print('PARALLEL EVOLVE')
    else:
        if myid == 0 and verbose: print('SEQUENTIAL EVOLVE')


    for t in domain.evolve(yieldstep = yieldstep, finaltime = finaltime):
        if myid == 0 and verbose : domain.write_time()

        # Record time series at known points
        time.append(domain.get_time())

        stage = domain.get_quantity('stage')

        for i in range(4):
            if tri_ids[i] > -1:
                gauge_values[i].append(stage.centroid_values[tri_ids[i]])


    #----------------------------------------
    # Setup test arrays during sequential run
    #----------------------------------------
    if not parallel:
        G = []
        for i in range(4):
            G.append(gauge_values[i])

    success = True

    for i in range(4):
        if tri_ids[i] > -1:
            #print num.max(num.array(gauge_values[i])- num.array(G[i]))
            success = success and num.allclose(gauge_values[i], G[i])

    assert_(success)

    return G, interpolation_points

# Test an nprocs-way run of the shallow water equations
# against the sequential code.

class Test_parallel_sw_flow(unittest.TestCase):
    def test_parallel_sw_flow(self):
        if verbose : print("Expect this test to fail if not run from the parallel directory.")

        cmd = anuga.mpicmd(os.path.abspath(__file__))
        result = os.system(cmd)

        # Just use the normal Python assert
        msg = 'Result == %i, expected 0' % result
        assert result == 0, msg

# Because we are doing assertions outside of the TestCase class
# the PyUnit defined assert_ function can't be used.
# FIXME (Ole): Why not use the normal Python assert?
def assert_(condition, msg="Assertion Failed"):
    if condition == False:
        #pypar.finalize()
        raise_(AssertionError, msg)

if __name__=="__main__":
    if numprocs == 1:
        runner = unittest.TextTestRunner()
        suite = unittest.makeSuite(Test_parallel_sw_flow, 'test')
        runner.run(suite)
    else:

        #------------------------------------------
        # Run the sequential code on each processor
        # and save results at 4 gauge stations to
        # array G
        #------------------------------------------
        barrier()
        if myid == 0 and verbose: print('SEQUENTIAL START')

        G , interpolation_points = run_simulation(parallel=False,verbose=verbose)
        G = num.array(G,float)

        barrier()

        #------------------------------------------
        # Run the code code and compare sequential
        # results at 4 gauge stations
        #------------------------------------------
        if myid ==0 and verbose: print('PARALLEL START')

        from anuga.utilities.parallel_abstraction import global_except_hook
        import sys
        sys.excepthook = global_except_hook

        run_simulation(parallel=True, G=G, seq_interpolation_points = interpolation_points, verbose= verbose)

        finalize()
