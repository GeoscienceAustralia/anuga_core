"""
Simple water flow example using ANUGA

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



from anuga import Domain
from anuga import Reflective_boundary
from anuga import Dirichlet_boundary
from anuga import Time_boundary
from anuga import Transmissive_boundary

from anuga import rectangular_cross_domain


from anuga_parallel import distribute, myid, numprocs, send, receive, barrier, finalize

from anuga_parallel.sequential_distribute import sequential_distribute_dump
from anuga_parallel.sequential_distribute import sequential_distribute_load


#--------------------------------------------------------------------------
# Setup parameters
#--------------------------------------------------------------------------
yieldstep = 0.25
finaltime = 3.0
nprocs = 4
N = 29
M = 29 
verbose = False

#---------------------------------
# Setup Functions
#---------------------------------
def topography(x,y): 
    return -x/2    

###########################################################################
# Setup Test
##########################################################################
def evolution_test(parallel=False, verbose=False):

    #--------------------------------------------------------------------------
    # Setup computational domain and quantities
    #--------------------------------------------------------------------------
    if myid == 0:
        domain = rectangular_cross_domain(M, N)
        domain.set_name('sequential_dist_runup')                    # Set sww filename
        domain.set_datadir('.')   
        domain.set_quantity('elevation', topography) # Use function for elevation
        domain.set_quantity('friction', 0.0)         # Constant friction 
        domain.set_quantity('stage', expression='elevation') # Dry initial stage
    else:
        domain = None
        
    #--------------------------------------------------------------------------
    # Create pickled partition
    #--------------------------------------------------------------------------
    if myid == 0:
        if verbose: print 'DUMPING PARTITION DATA'
        sequential_distribute_dump(domain, numprocs, verbose=verbose)    

    #--------------------------------------------------------------------------
    # Create the parallel domains
    #--------------------------------------------------------------------------
    if parallel:
        
        if myid == 0 and verbose : print 'DISTRIBUTING TO PARALLEL DOMAIN'
        pdomain = distribute(domain, verbose=verbose)
        
        if myid == 0 and verbose : print 'LOADING IN PARALLEL DOMAIN'
        sdomain = sequential_distribute_load(filename='sequential_dist_runup', verbose = verbose)
        
        
    if myid == 0 and verbose: print 'EVOLVING pdomain'    
    setup_and_evolve(pdomain, verbose=verbose)
 
    if myid == 0 and verbose: print 'EVOLVING sdomain'   
    setup_and_evolve(sdomain, verbose=verbose)
    
    
    assert num.allclose(pdomain.quantities['stage'].centroid_values, sdomain.quantities['stage'].centroid_values)
    assert num.allclose(pdomain.quantities['stage'].vertex_values, sdomain.quantities['stage'].vertex_values)
    
    assert num.allclose(pdomain.vertex_coordinates, sdomain.vertex_coordinates)
    assert num.allclose(pdomain.centroid_coordinates, sdomain.centroid_coordinates)
    
    





def setup_and_evolve(domain, verbose=False):

    #--------------------------------------------------------------------------
    # Setup domain parameters
    #--------------------------------------------------------------------------
    domain.set_flow_algorithm('DE0')
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
    # Evolve
    #------------------------------------------------------------------------------
    for t in domain.evolve(yieldstep = yieldstep, finaltime = finaltime):
        if myid == 0 and verbose : domain.write_time()



# Test an nprocs-way run of the shallow water equations
# against the sequential code.

class Test_parallel_sw_flow(unittest.TestCase):
    def test_parallel_sw_flow(self):
        if verbose : print "Expect this test to fail if not run from the parallel directory."
        result = os.system("mpirun -np %d python test_sequential_dist_sw_flow.py" % nprocs)
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
        # Run the code code and compare sequential
        # results at 4 gauge stations
        #------------------------------------------
        if myid ==0 and verbose: print 'PARALLEL START'

        evolution_test(parallel=True, verbose=verbose)
        
        finalize()


