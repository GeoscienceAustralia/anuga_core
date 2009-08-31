
"""Simple water flow example using ANUGA

Water driven up a linear slope and time varying boundary,
similar to a beach environment

This is a very simple test of the parallel algorithm using the simplified parallel API
"""


#------------------------------------------------------------------------------
# Import necessary modules
#------------------------------------------------------------------------------
import unittest
import os
import sys
#import pypar
import numpy as num

from anuga.interface import Domain
from anuga.interface import Reflective_boundary
from anuga.interface import Dirichlet_boundary
from anuga.interface import Time_boundary
from anuga.interface import Transmissive_boundary

from anuga.interface import rectangular_cross_domain


from anuga_parallel.interface import distribute, myid, numprocs, send, receive, barrier

#--------------------------------------------------------------------------
# Setup parameters
#--------------------------------------------------------------------------
yieldstep = 0.25
finaltime = 6.0
nprocs = 4
N = 25
M = 25
verbose =  True

#---------------------------------
# Setup Functions
#---------------------------------
def topography(x,y): 
    return -x/2    

###########################################################################
# Setup Test
##########################################################################
def evolution_test(parallel=False, G = None, seq_interpolation_points=None):

    #--------------------------------------------------------------------------
    # Setup computational domain and quantities
    #--------------------------------------------------------------------------
    domain = rectangular_cross_domain(N, M) 
    domain.set_quantity('elevation', topography) # Use function for elevation
    domain.set_quantity('friction', 0.0)         # Constant friction 
    domain.set_quantity('stage', expression='elevation') # Dry initial stage

    #--------------------------------------------------------------------------
    # Create the parallel domain
    #--------------------------------------------------------------------------
    if parallel:
        if myid == 0 and verbose : print 'DISTRIBUTING PARALLEL DOMAIN'
        domain = distribute(domain, verbose=False)

    #--------------------------------------------------------------------------
    # Setup domain parameters
    #--------------------------------------------------------------------------
    domain.set_name('runup')                    # Set sww filename
    domain.set_datadir('.')                     # Set output dir

    domain.set_default_order(1)        
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
    # between to centroids, so in the parallel run we
    # reset the interpolation points to the centroids
    # found in the sequential run
    #------------------------------------------------------------------------------
    interpolation_points = [[0.4,0.5], [0.6,0.5], [0.8,0.5], [0.9,0.5]]

    if parallel:
        interpolation_points = seq_interpolation_points


    gauge_values = []
    tri_ids = []
    for i, point in enumerate(interpolation_points):
        gauge_values.append([]) # Empty list for timeseries

        #if is_inside_polygon(point, domain.get_boundary_polygon()):
        try:
            k = domain.get_triangle_containing_point(point)
            if domain.tri_full_flag[k] == 1:
                tri_ids.append(k)
            else:
                tri_ids.append(-1)            
        except:
            tri_ids.append(-2)


    if verbose: print 'P%d has points = %s' %(myid, tri_ids)

    if not parallel:
        c_coord = domain.get_centroid_coordinates()
        interpolation_points = []
        for id in tri_ids:
            if id<1:
                print 'ERROR: All interpolation points be within the sequential domain!'
            interpolation_points.append(c_coord[id,:])
            
    #------------------------------------------------------------------------------
    # Evolve system through time
    #------------------------------------------------------------------------------
    time = []

    if parallel:
        if myid == 0 and verbose: print 'PARALLEL EVOLVE'
    else:
        if myid == 0 and verbose: print 'SEQUENTIAL EVOLVE'
    

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
            success = success and num.allclose(gauge_values[i], G[i])

    assert_(success)
    
    return G, interpolation_points

# Test an nprocs-way run of the shallow water equations
# against the sequential code.

class Test_parallel_sw_flow(unittest.TestCase):
    def test_parallel_sw_flow(self):
        print "Expect this test to fail if not run from the parallel directory."
        result = os.system("mpirun -np %d python test_parallel_sw_flow.py" % nprocs)
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
        suite = unittest.makeSuite(Test_parallel_sw_flow, 'test')
        runner.run(suite)
    else:

        #------------------------------------------
        # Run the sequential code on each processor
        # and save results at 4 gauge stations to
        # array G
        #------------------------------------------
        barrier()
        if myid == 0 and verbose: print 'SEQUENTIAL START'

        G , interpolation_points = evolution_test(parallel=False)
        G = num.array(G,num.float)

        barrier()
        
        #------------------------------------------
        # Run the code code and compare sequential
        # results at 4 gauge stations
        #------------------------------------------
        if myid ==0 and verbose: print 'PARALLEL START'

        evolution_test(parallel=True, G=G, seq_interpolation_points = interpolation_points)
        



