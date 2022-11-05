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

import anuga

from anuga import Domain
from anuga import Reflective_boundary
from anuga import Dirichlet_boundary
from anuga import Time_boundary
from anuga import Transmissive_boundary

from anuga import rectangular_cross_domain


from anuga import distribute, myid, numprocs, send, receive, barrier, finalize

from anuga.parallel.sequential_distribute import sequential_distribute_dump
from anuga.parallel.sequential_distribute import sequential_distribute_load

import anuga.utilities.plot_utils as util

# Setup to skip test if mpi4py not available
import sys
try:
    import mpi4py
except ImportError:
    pass

import pytest


#--------------------------------------------------------------------------
# Setup parameters
#--------------------------------------------------------------------------
yieldstep = 0.25
finaltime = 1.0
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
    if myid == 0:
        domain = rectangular_cross_domain(M, N)
        domain.set_name('odomain')                    # Set sww filename
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
        if verbose: print('DUMPING PARTITION DATA')
        sequential_distribute_dump(domain, numprocs, verbose=verbose, parameters=new_parameters)    

    #--------------------------------------------------------------------------
    # Create the parallel domains
    #--------------------------------------------------------------------------
    if parallel:
        
        if myid == 0 and verbose : print('DISTRIBUTING TO PARALLEL DOMAIN')
        pdomain = distribute(domain, verbose=verbose, parameters=new_parameters)
        pdomain.set_name('pdomain')
        
        if myid == 0 and verbose : print('LOADING IN PARALLEL DOMAIN')
        sdomain = sequential_distribute_load(filename='odomain', verbose = verbose)
        sdomain.set_name('sdomain')
        
    if myid == 0 and verbose: print('EVOLVING pdomain')    
    setup_and_evolve(pdomain, verbose=verbose)
 
    if myid == 0 and verbose: print('EVOLVING sdomain')   
    setup_and_evolve(sdomain, verbose=verbose)
    
    if myid == 0:
        if verbose: print('EVOLVING odomain')   
        setup_and_evolve(domain, verbose=verbose)
    

    if myid == 0 and verbose:
        parameter_file=open('odomain.txt', 'w')
        from pprint import pprint
        pprint(domain.get_algorithm_parameters(),parameter_file,indent=4)
        parameter_file.close()

        parameter_file=open('sdomain.txt', 'w')
        from pprint import pprint
        pprint(sdomain.get_algorithm_parameters(),parameter_file,indent=4)
        parameter_file.close()

        parameter_file=open('pdomain.txt', 'w')
        from pprint import pprint
        pprint(pdomain.get_algorithm_parameters(),parameter_file,indent=4)
        parameter_file.close()        
    
    assert num.allclose(pdomain.quantities['stage'].centroid_values, sdomain.quantities['stage'].centroid_values)
    assert num.allclose(pdomain.quantities['stage'].vertex_values, sdomain.quantities['stage'].vertex_values)
    
    assert num.allclose(pdomain.vertex_coordinates, sdomain.vertex_coordinates)
    assert num.allclose(pdomain.centroid_coordinates, sdomain.centroid_coordinates)
    
    

    #---------------------------------
    # Now compare the merged sww files
    #---------------------------------
    if myid == 0:
        if verbose: print('COMPARING SWW FILES')
        
        odomain_v = util.get_output('odomain.sww')
        odomain_c = util.get_centroids(odomain_v)

        pdomain_v = util.get_output('pdomain.sww')
        pdomain_c = util.get_centroids(pdomain_v)
        
        sdomain_v = util.get_output('sdomain.sww')
        sdomain_c = util.get_centroids(sdomain_v)

        # Test some values against the original ordering
        
        if verbose:
            
            order = 2
            print('PDOMAIN CENTROID VALUES')
            print(num.linalg.norm(odomain_c.x-pdomain_c.x,ord=order))
            print(num.linalg.norm(odomain_c.y-pdomain_c.y,ord=order))
            print(num.linalg.norm(odomain_c.stage[-1]-pdomain_c.stage[-1],ord=order))
            print(num.linalg.norm(odomain_c.xmom[-1]-pdomain_c.xmom[-1],ord=order))
            print(num.linalg.norm(odomain_c.ymom[-1]-pdomain_c.ymom[-1],ord=order))
            print(num.linalg.norm(odomain_c.xvel[-1]-pdomain_c.xvel[-1],ord=order))
            print(num.linalg.norm(odomain_c.yvel[-1]-pdomain_c.yvel[-1],ord=order))        
            
             
            print('SDOMAIN CENTROID VALUES')        
            print(num.linalg.norm(odomain_c.x-sdomain_c.x,ord=order))
            print(num.linalg.norm(odomain_c.y-sdomain_c.y,ord=order))
            print(num.linalg.norm(odomain_c.stage[-1]-sdomain_c.stage[-1],ord=order))
            print(num.linalg.norm(odomain_c.xmom[-1]-sdomain_c.xmom[-1],ord=order))
            print(num.linalg.norm(odomain_c.ymom[-1]-sdomain_c.ymom[-1],ord=order))
            print(num.linalg.norm(odomain_c.xvel[-1]-sdomain_c.xvel[-1],ord=order))
            print(num.linalg.norm(odomain_c.yvel[-1]-sdomain_c.yvel[-1],ord=order))
            
            print('PDOMAIN VERTEX VALUES')        
            print(num.linalg.norm(odomain_v.stage[-1]-pdomain_v.stage[-1],ord=order))
            print(num.linalg.norm(odomain_v.xmom[-1]-pdomain_v.xmom[-1],ord=order))
            print(num.linalg.norm(odomain_v.ymom[-1]-pdomain_v.ymom[-1],ord=order))
            print(num.linalg.norm(odomain_v.xvel[-1]-pdomain_v.xvel[-1],ord=order))
            print(num.linalg.norm(odomain_v.yvel[-1]-pdomain_v.yvel[-1],ord=order))
            
            print('SDOMAIN VERTEX VALUES')     
            print(num.linalg.norm(odomain_v.stage[-1]-sdomain_v.stage[-1],ord=order))
            print(num.linalg.norm(odomain_v.xmom[-1]-sdomain_v.xmom[-1],ord=order))
            print(num.linalg.norm(odomain_v.ymom[-1]-sdomain_v.ymom[-1],ord=order))
            print(num.linalg.norm(odomain_v.xvel[-1]-sdomain_v.xvel[-1],ord=order))
            print(num.linalg.norm(odomain_v.yvel[-1]-sdomain_v.yvel[-1],ord=order))
            
            
            

        assert num.allclose(odomain_c.stage,pdomain_c.stage)
        assert num.allclose(odomain_c.xmom,pdomain_c.xmom)
        assert num.allclose(odomain_c.ymom,pdomain_c.ymom)
        assert num.allclose(odomain_c.xvel,pdomain_c.xvel)
        assert num.allclose(odomain_c.yvel,pdomain_c.yvel)
        
        assert num.allclose(odomain_v.x,pdomain_v.x)
        assert num.allclose(odomain_v.y,pdomain_v.y)
                
        assert num.linalg.norm(odomain_v.x-pdomain_v.x,ord=0) == 0
        assert num.linalg.norm(odomain_v.y-pdomain_v.y,ord=0) == 0
        assert num.linalg.norm(odomain_v.stage[-1]-pdomain_v.stage[-1],ord=0) < 100
        assert num.linalg.norm(odomain_v.xmom[-1]-pdomain_v.xmom[-1],ord=0) < 100 
        assert num.linalg.norm(odomain_v.ymom[-1]-pdomain_v.ymom[-1],ord=0) < 100
        assert num.linalg.norm(odomain_v.xvel[-1]-pdomain_v.xvel[-1],ord=0) < 100
        assert num.linalg.norm(odomain_v.yvel[-1]-pdomain_v.yvel[-1],ord=0) < 100     
        
        assert num.allclose(odomain_c.x,sdomain_c.x)
        assert num.allclose(odomain_c.y,sdomain_c.y)
        assert num.allclose(odomain_c.stage,sdomain_c.stage)
        assert num.allclose(odomain_c.xmom,sdomain_c.xmom)
        assert num.allclose(odomain_c.ymom,sdomain_c.ymom)
        assert num.allclose(odomain_c.xvel,sdomain_c.xvel)
        assert num.allclose(odomain_c.yvel,sdomain_c.yvel)
        
        assert num.allclose(odomain_v.x,sdomain_v.x)
        assert num.allclose(odomain_v.y,sdomain_v.y)
        
        order = 0
        assert num.linalg.norm(odomain_v.x-sdomain_v.x,ord=order) == 0
        assert num.linalg.norm(odomain_v.y-sdomain_v.y,ord=order) == 0
        assert num.linalg.norm(odomain_v.stage[-1]-sdomain_v.stage[-1],ord=order) < 100
        assert num.linalg.norm(odomain_v.xmom[-1]-sdomain_v.xmom[-1],ord=order) < 100
        assert num.linalg.norm(odomain_v.ymom[-1]-sdomain_v.ymom[-1],ord=order) < 100
        assert num.linalg.norm(odomain_v.xvel[-1]-sdomain_v.xvel[-1],ord=order) < 100
        assert num.linalg.norm(odomain_v.yvel[-1]-sdomain_v.yvel[-1],ord=order) < 100        
                
        # COMPARE CENTROID PDOMAIN SDOMAIN  
        assert num.allclose(pdomain_c.x,sdomain_c.x)
        assert num.allclose(pdomain_c.y,sdomain_c.y)
        assert num.allclose(pdomain_c.stage[-1],sdomain_c.stage[-1])
        assert num.allclose(pdomain_c.xmom[-1],sdomain_c.xmom[-1])
        assert num.allclose(pdomain_c.ymom[-1],sdomain_c.ymom[-1])
        assert num.allclose(pdomain_c.xvel[-1],sdomain_c.xvel[-1])
        assert num.allclose(pdomain_c.yvel[-1],sdomain_c.yvel[-1])
            
            
        # COMPARE VERTEX PDOMAIN SDOMAIN
        assert num.allclose(pdomain_v.x,sdomain_v.x)
        assert num.allclose(pdomain_v.y,sdomain_v.y)
        assert num.allclose(pdomain_v.stage[-1],sdomain_v.stage[-1])
        assert num.allclose(pdomain_v.xmom[-1],sdomain_v.xmom[-1])
        assert num.allclose(pdomain_v.ymom[-1],sdomain_v.ymom[-1])
        assert num.allclose(pdomain_v.xvel[-1],sdomain_v.xvel[-1])
        assert num.allclose(pdomain_v.yvel[-1],sdomain_v.yvel[-1])   
        
        
        import os
        os.remove('odomain.sww')
        os.remove('pdomain.sww')
        os.remove('sdomain.sww')
        os.remove('odomain_P3_0.pickle')
        os.remove('odomain_P3_1.pickle')
        os.remove('odomain_P3_2.pickle')
        #os.remove('odomain_P4_3.pickle')
        import glob
        [ os.remove(fl) for fl in glob.glob('*.npy') ]
        
        
def setup_and_evolve(domain, verbose=False):

    #--------------------------------------------------------------------------
    # Setup domain parameters
    #--------------------------------------------------------------------------
    domain.set_flow_algorithm('DE0')
    #domain.set_store_vertices_uniquely()

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
        #if myid == 0 and verbose : print domain.quantities['stage'].get_maximum_value()


    domain.sww_merge(delete_old=True)
    

# Test an nprocs-way run of the shallow water equations
# against the sequential code.

@pytest.mark.skipif('mpi4py' not in sys.modules,
                    reason="requires the mpi4py module")
class Test_parallel_sw_flow(unittest.TestCase):
    def test_parallel_sw_flow(self):
        if verbose : print("Expect this test to fail if not run from the parallel directory.")

        cmd = anuga.mpicmd(os.path.abspath(__file__))
        result = os.system(cmd)
        
        assert_(result == 0)

# Because we are doing assertions outside of the TestCase class
# the PyUnit defined assert_ function can't be used.
def assert_(condition, msg="Assertion Failed"):
    if condition == False:
        #pypar.finalize()
        raise (AssertionError, msg)

if __name__=="__main__":
    if numprocs == 1: 
        runner = unittest.TextTestRunner()
        suite = unittest.makeSuite(Test_parallel_sw_flow, 'test')
        runner.run(suite)
    else:

        from anuga.utilities.parallel_abstraction import global_except_hook
        import sys
        sys.excepthook = global_except_hook
        
        #------------------------------------------
        # Run the codel and compare sequential
        # results at 4 gauge stations
        #------------------------------------------
        if myid ==0 and verbose: print('PARALLEL START')

        run_simulation(parallel=True, verbose=verbose)
        
        finalize()


