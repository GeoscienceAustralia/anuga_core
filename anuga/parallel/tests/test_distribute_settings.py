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

from anuga import Geo_reference

from anuga import rectangular_cross_domain


from anuga import distribute, myid, numprocs, send, receive, barrier, finalize

from anuga.parallel.sequential_distribute import sequential_distribute_dump
from anuga.parallel.sequential_distribute import sequential_distribute_load

import anuga.utilities.plot_utils as util


#--------------------------------------------------------------------------
# Setup parameters
#--------------------------------------------------------------------------
yieldstep = 0.25
finaltime = 3.0
nprocs = 4
N = 29
M = 29 
verbose = False


new_parameters = {}
new_parameters['ghost_layer_width'] = 2

#---------------------------------
# Setup Functions
#---------------------------------
def topography(x,y): 
    return -x/2    

###########################################################################
# Setup Test
##########################################################################
def run_simulation(parallel=False, verbose=False):

    #--------------------------------------------------------------------------
    # Setup computational domain and quantities
    #--------------------------------------------------------------------------
    domain = rectangular_cross_domain(M, N)
    domain.set_name('odomain')                    # Set sww filename
    domain.set_datadir('.')   
    domain.set_quantity('elevation', topography) # Use function for elevation
    domain.set_quantity('friction', 0.0)         # Constant friction 
    domain.set_quantity('stage', expression='elevation') # Dry initial stage
    
    domain.set_quantities_to_be_stored({'elevation': 2,'stage': 2,'xmomentum': 2,'ymomentum': 2})
    domain.set_store_vertices_uniquely(False)
    domain.set_flow_algorithm('DE1')
    georef = Geo_reference(zone=56,xllcorner=100000.0,yllcorner=200000.0)
    domain.set_georeference(georef)
        
    #--------------------------------------------------------------------------
    # Create pickled partition
    #--------------------------------------------------------------------------
    if myid == 0:
        if verbose: print 'DUMPING PARTITION DATA'
        sequential_distribute_dump(domain, numprocs, verbose=verbose, parameters=new_parameters)    

    #--------------------------------------------------------------------------
    # Create the parallel domains
    #--------------------------------------------------------------------------
    if parallel:
        
        if myid == 0 and verbose : print 'DISTRIBUTING TO PARALLEL DOMAIN'
        pdomain = distribute(domain, verbose=verbose, parameters=new_parameters)
        pdomain.set_name('pdomain')
        
        if myid == 0 and verbose : print 'LOADING IN PARALLEL DOMAIN'
        sdomain = sequential_distribute_load(filename='odomain', verbose = verbose)
        sdomain.set_name('sdomain')
        

        assert domain.get_datadir() == pdomain.get_datadir()
        assert domain.get_store() == pdomain.get_store()
        assert domain.get_store_centroids() == pdomain.get_store_centroids()
        assert domain.smooth == pdomain.smooth
        assert domain.reduction == pdomain.reduction
        assert domain.minimum_storable_height == pdomain.minimum_storable_height
        assert domain.get_flow_algorithm() == pdomain.get_flow_algorithm()
        assert domain.get_minimum_allowed_height() == pdomain.get_minimum_allowed_height()
        assert domain.geo_reference == pdomain.geo_reference
        
        assert domain.get_datadir() == sdomain.get_datadir()
        assert domain.get_store() == sdomain.get_store()
        assert domain.get_store_centroids() == sdomain.get_store_centroids()
        assert domain.smooth == sdomain.smooth
        assert domain.reduction == sdomain.reduction
        assert domain.minimum_storable_height == sdomain.minimum_storable_height
        assert domain.get_flow_algorithm() == sdomain.get_flow_algorithm()
        assert domain.get_minimum_allowed_height() == sdomain.get_minimum_allowed_height()
        assert domain.geo_reference == sdomain.geo_reference
        


    

# Test an nprocs-way run of the shallow water equations
# against the sequential code.

class Test_parallel_sw_flow(unittest.TestCase):
    def test_parallel_sw_flow(self):
        if verbose : print "Expect this test to fail if not run from the parallel directory."

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
        suite = unittest.makeSuite(Test_parallel_sw_flow, 'test')
        runner.run(suite)
    else:


        
        #------------------------------------------
        # Run the codel and compare sequential
        # results at 4 gauge stations
        #------------------------------------------
        if myid ==0 and verbose: print 'PARALLEL START'

        run_simulation(parallel=True, verbose=verbose)
        
        finalize()


