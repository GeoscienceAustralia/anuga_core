"""
Test parallel and sequential results of riverwall procedure
"""


#------------------------------------------------------------------------------
# Import necessary modules
#------------------------------------------------------------------------------
import unittest
import os
import sys
import anuga
import numpy as num

from anuga import myid, finalize, distribute, barrier, numprocs

import anuga.utilities.plot_utils as util

from math import exp


alg = 'DE0'
verbose = False

scale_me=1.0
nprocs = 3

## Set up mesh
boundaryPolygon=[ [0., 0.], [0., 100.], [100.0, 100.0], [100.0, 0.0]]
midResPolygon=[ [30., 30.], [30., 70.], [70., 70.], [70., 30.]]
higherResPolygon=[ [40., 40.], [40., 60.], [60., 60.], [60., 40.]]
# Riverwall = list of lists, each with a set of x,y,z (and optional QFactor) values
riverWall={ 'centralWall':
                [ [50., 0.0, -0.0],
                  [50., 45., -0.0],
                  [50., 46., -0.2],
                  [50., 54., -0.2],
                  [50., 55., -0.0], 
                  [50., 100.0, -0.0]] 
          }

riverWall_Par={'centralWall':{'Qfactor':1.0}}
# Try to avoid any shallow-water type solution -- becomes unstable
#riverWall_Par={'centralWall':{'Qfactor':1.0, 's1': 0.999, 's2':0.9999, 'h1':100, 'h2':150}}

# The boundary polygon + riverwall breaks the mesh into multiple regions
# Must define the resolution in these areas with an xy point + maximum area
# Otherwise triangle.c gets confused
regionPtAreas=[ [99., 99., 10.0*10.0*0.5],
                [1., 1., 10.0*10.0*0.5],
                [45., 45., 1.0*1.0*0.5],
                [55., 55., 1.0*1.0*0.5],
                [65., 65., 3.0*3.0*0.5],
                [35., 35., 3.0*3.0*0.5] ]




  

###########################################################################
# Setup Test
##########################################################################
def run_simulation(parallel=False, verbose=False):

    #--------------------------------------------------------------------------
    # Setup computational domain and quantities
    #--------------------------------------------------------------------------
    if myid == 0:
        anuga.create_mesh_from_regions(boundaryPolygon, 
                                 boundary_tags={'left': [0],
                                                'top': [1],
                                                'right': [2],
                                                'bottom': [3]},
                                   maximum_triangle_area = 1.0e+20,
                                   minimum_triangle_angle = 28.0,
                                   filename = 'runup.msh',
                                   interior_regions = [ [higherResPolygon, 1.*1.*0.5],
                                                        [midResPolygon, 3.0*3.0*0.5]],
                                   breaklines=riverWall.values(),
                                   use_cache=False,
                                   verbose=verbose,
                                   regionPtArea=regionPtAreas)
        
        sdomain=anuga.create_domain_from_file('runup.msh')
        
        
        sdomain.set_flow_algorithm(alg)
        
        
        sdomain.set_name('s_riverwall')                         
        sdomain.set_datadir('.')                         
        sdomain.set_store_vertices_uniquely()
        
        #------------------
        # Define topography
        #------------------
        
        def topography(x,y):
            return -x/150.*scale_me 
        
        def stagefun(x,y):
            stg=-0.5*scale_me
            return stg 
        
        sdomain.set_quantity('elevation',topography)     # Use function for elevation
        sdomain.set_quantity('friction',0.03)             # Constant friction
        sdomain.set_quantity('stage', stagefun)              # Constant negative initial stage
    else:
        sdomain = None
        
    #--------------------------------------------------------------------------
    # Create the parallel domains
    #--------------------------------------------------------------------------
    if parallel:
        
        if myid == 0 and verbose : print 'DISTRIBUTING TO PARALLEL DOMAIN'
        pdomain = distribute(sdomain, verbose=verbose)
        pdomain.set_name('p_riverwall')
        pdomain.set_store_vertices_uniquely()
        
        
    if myid == 0 and verbose: 
        print 60*'='
        print 'EVOLVING pdomain'
        print 60*'='
            
    setup_and_evolve(pdomain, verbose=verbose)
 
    barrier()
   
    if myid == 0:
        if verbose: 
            print 60*'='
            print 'EVOLVING sdomain'
            print 60*'='  
        setup_and_evolve(sdomain, verbose=verbose)
      
    barrier()
    
    #---------------------------------
    # Now compare the merged sww files
    #---------------------------------
    if myid == 0:
        if verbose: print 'COMPARING SWW FILES'
        
        sdomain_v = util.get_output('s_riverwall.sww')
        sdomain_c = util.get_centroids(sdomain_v)

        pdomain_v = util.get_output('p_riverwall.sww')
        pdomain_c = util.get_centroids(pdomain_v)
        

        # Test some values against the original ordering
        
        if verbose:
            
            order = 0
            print 'PDOMAIN CENTROID VALUES'
            print num.linalg.norm(sdomain_c.x-pdomain_c.x,ord=order)
            print num.linalg.norm(sdomain_c.y-pdomain_c.y,ord=order)
            print num.linalg.norm(sdomain_c.stage[-1]-pdomain_c.stage[-1],ord=order)
            print num.linalg.norm(sdomain_c.xmom[-1]-pdomain_c.xmom[-1],ord=order)
            print num.linalg.norm(sdomain_c.ymom[-1]-pdomain_c.ymom[-1],ord=order)
            print num.linalg.norm(sdomain_c.xvel[-1]-pdomain_c.xvel[-1],ord=order)
            print num.linalg.norm(sdomain_c.yvel[-1]-pdomain_c.yvel[-1],ord=order)        
            
        assert num.allclose(sdomain_c.stage,pdomain_c.stage)
        assert num.allclose(sdomain_c.xmom,pdomain_c.xmom)
        assert num.allclose(sdomain_c.ymom,pdomain_c.ymom)
        assert num.allclose(sdomain_c.xvel,pdomain_c.xvel)
        assert num.allclose(sdomain_c.yvel,pdomain_c.yvel)
        
        assert num.allclose(sdomain_v.x,pdomain_v.x)
        assert num.allclose(sdomain_v.y,pdomain_v.y)
        
        
        import os
        os.remove('s_riverwall.sww')
        os.remove('p_riverwall.sww')
        os.remove('runup.msh')
        
        
        
def setup_and_evolve(domain, verbose=False):
    

    domain.riverwallData.create_riverwalls(riverWall, riverWall_Par, verbose=verbose)
    
    #--------------------------
    # Setup boundary conditions
    #--------------------------
    Br=anuga.Reflective_boundary(domain)            # Solid reflective wall
    
    def boundaryFun(t):
        output=-0.4*exp(-t/100.)-0.1
        output=min(output,-0.11)
        #output=min(output,-0.3)
        return output  
    

    Bin_tmss = anuga.Transmissive_momentum_set_stage_boundary(domain=domain, function = boundaryFun) 
    
    #----------------------------------------------
    # Associate boundary tags with boundary objects
    #----------------------------------------------
    domain.set_boundary({'left': Br, 'right': Bin_tmss, 'top': Br, 'bottom':Br})
    

    
    #------------------------------
    #Evolve the system through time
    #------------------------------
    if verbose: print 'Evolve'
    for t in domain.evolve(yieldstep=10.0,finaltime=150.0):
        if myid == 0 and verbose: print domain.timestepping_statistics()
    
    domain.sww_merge(delete_old=True)

    

# Test an nprocs-way run of the shallow water equations
# against the sequential code.

class Test_parallel_riverwall(unittest.TestCase):
    def test_parallel_riverwall(self):
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
        suite = unittest.makeSuite(Test_parallel_riverwall, 'test')
        runner.run(suite)
    else:

        #------------------------------------------
        # Run the code code and compare sequential
        # and parallel values
        #------------------------------------------
        if myid ==0 and verbose: print 'PARALLEL START'
        run_simulation(parallel=True, verbose=verbose)
        finalize()


