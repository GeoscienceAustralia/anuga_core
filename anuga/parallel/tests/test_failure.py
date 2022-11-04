
import os.path
import sys

from anuga.utilities.system_tools import get_pathname_from_package
from anuga.geometry.polygon_function import Polygon_function
        
from anuga.abstract_2d_finite_volumes.mesh_factory import rectangular_cross
from anuga.abstract_2d_finite_volumes.quantity import Quantity
from anuga.abstract_2d_finite_volumes.util import file_function

import anuga

from anuga.structures.boyd_box_operator import Boyd_box_operator
from anuga.structures.inlet_operator import Inlet_operator
                            
#from anuga.culvert_flows.culvert_routines import boyd_generalised_culvert_model
     
from math import pi, pow, sqrt

import numpy as num
#from parallel_inlet_operator import Parallel_Inlet_operator
from anuga.parallel import distribute, myid, numprocs, finalize
from anuga.geometry.polygon import inside_polygon, is_inside_polygon, line_intersect

#from parallel_operator_factory import Inlet_operator, Boyd_box_operator
from anuga.utilities import parallel_abstraction as pypar
import random


"""test_that_culvert_runs_rating

This test exercises the culvert and checks values outside rating curve
are dealt with       
"""
verbose = False
path = get_pathname_from_package('anuga.culvert_flows')    

length = 40.
width = 15.

dx = dy = 0.5          # Resolution: Length of subdivisions on both axes

#----------------------------------------------------------------------
# Setup initial conditions
#----------------------------------------------------------------------

def topography(x, y):
    """Set up a weir
    
    A culvert will connect either side
    """
    # General Slope of Topography
    z=-x/1000
    
    N = len(x)
    for i in range(N):

       # Sloping Embankment Across Channel
        if 5.0 < x[i] < 10.1:
            # Cut Out Segment for Culvert face                
            if  1.0+(x[i]-5.0)/5.0 <  y[i]  < 4.0 - (x[i]-5.0)/5.0: 
               z[i]=z[i]
            else:
               z[i] +=  0.5*(x[i] -5.0)    # Sloping Segment  U/S Face
        if 10.0 < x[i] < 12.1:
           z[i] +=  2.5                    # Flat Crest of Embankment
        if 12.0 < x[i] < 14.5:
            # Cut Out Segment for Culvert face                
            if  2.0-(x[i]-12.0)/2.5 <  y[i]  < 3.0 + (x[i]-12.0)/2.5:
               z[i]=z[i]
            else:
               z[i] +=  2.5-1.0*(x[i] -12.0) # Sloping D/S Face
                   
        
    return z

filename=os.path.join(path, 'example_rating_curve.csv')


mod_path = get_pathname_from_package('anuga.parallel')

line0 = [[10.0, 10.0], [30.0, 10.0]]
#line0 = [[29.0, 10.0], [30.0, 10.0]]
line1 = [[0.0, 10.0], [0.0, 15.0]]
Q0 = file_function(os.path.join(mod_path,'data','test_hydrograph.tms'), quantities=['hydrograph'])
Q1 = 5.0

samples = 50

def run_simulation(parallel = False, control_data = None, test_points = None, verbose = False):
    success = True

##-----------------------------------------------------------------------
## Setup domain
##-----------------------------------------------------------------------

    points, vertices, boundary = rectangular_cross(int(length/dx),
                                                   int(width/dy),
                                                   len1=length, 
                                                   len2=width)

    domain = anuga.Domain(points, vertices, boundary)   
    domain.set_store(False)
    #domain.set_name('output_parallel_failure')                 # Output name
    domain.set_default_order(2)

##-----------------------------------------------------------------------
## Distribute domain
##-----------------------------------------------------------------------

    if parallel: domain = distribute(domain)
    

##-----------------------------------------------------------------------
## Setup boundary conditions
##-----------------------------------------------------------------------
    
    domain.set_quantity('elevation', topography) 
    domain.set_quantity('friction', 0.01)         # Constant friction 
    domain.set_quantity('stage',
                        expression='elevation')   # Dry initial condition

    
    Bi = anuga.Dirichlet_boundary([5.0, 0.0, 0.0])
    Br = anuga.Reflective_boundary(domain)              # Solid reflective wall
    domain.set_boundary({'left': Bi, 'right': Br, 'top': Br, 'bottom': Br})


##-----------------------------------------------------------------------
## Determine triangle index coinciding with test points
##-----------------------------------------------------------------------

    assert(test_points is not None)
    assert(len(test_points) == samples)

    tri_ids = []

    for point in test_points:
        try:
            k = domain.get_triangle_containing_point(point)
            if domain.tri_full_flag[k] == 1:
                tri_ids.append(k)
            else:
                tri_ids.append(-1)            
        except:
            tri_ids.append(-2)

    if verbose: print('P%d has points = %s' %(myid, tri_ids))

    if not parallel: control_data = []

    ################ Define Fractional Operators ##########################

    inlet0 = None
    inlet1 = None
    boyd_box0 = None
    
    #inlet0 = Inlet_operator(domain, line0, Q0, debug = True)
    #inlet1 = Inlet_operator(domain, line1, Q1, debug = True)
    
    ## # Enquiry point [ 19.    2.5] is contained in two domains in 4 proc case
    ## if myid == 5 and parallel:
    ##     boyd_box0 = Boyd_box_operator(domain,
    ##                                   end_points=[[9.0, 2.5],[13.0, 2.5]],
    ##                                   losses=1.5,
    ##                                   width=1.0,
    ##                                   apron=0.5,
    ##                                   use_momentum_jet=True,
    ##                                   use_velocity_head=False,
    ##                                   manning=0.013,
    ##                                   verbose=False)
    ## elif not parallel:
    ##     boyd_box0 = Boyd_box_operator(domain,
    ##                                   end_points=[[9.0, 2.5],[13.0, 2.5]],
    ##                                   losses=1.5,
    ##                                   width=1.0,
    ##                                   apron=0.5,
    ##                                   use_momentum_jet=True,
    ##                                   use_velocity_head=False,
    ##                                   manning=0.013,
    ##                                   verbose=False)

#    if parallel:
#        factory = Parallel_operator_factory(domain, debug = True)
#
#        inlet0 = factory.inlet_operator_factory(line0, Q0)
#        inlet1 = factory.inlet_operator_factory(line1, Q1)
#        
#        boyd_box0 = factory.boyd_box_operator_factory(end_points=[[9.0, 2.5],[19.0, 2.5]],
#                                          losses=1.5,
#                                          width=1.5,
#                                          apron=5.0,
#                                          use_momentum_jet=True,
#                                          use_velocity_head=False,
#                                          manning=0.013,
#                                          verbose=False)
#
#    else:
#        inlet0 = Inlet_operator(domain, line0, Q0)
#        inlet1 = Inlet_operator(domain, line1, Q1)
#
#        # Enquiry point [ 19.    2.5] is contained in two domains in 4 proc case
#        boyd_box0 = Boyd_box_operator(domain,
#                          end_points=[[9.0, 2.5],[19.0, 2.5]],
#                          losses=1.5,
#                          width=1.5,
#                          apron=5.0,
#                          use_momentum_jet=True,
#                          use_velocity_head=False,
#                          manning=0.013,
#                          verbose=False)
    
    #######################################################################

    ##-----------------------------------------------------------------------
    ## Evolve system through time
    ##-----------------------------------------------------------------------

    for t in domain.evolve(yieldstep = 1.0, finaltime = 4):
        if verbose: domain.write_time()

        #print domain.volumetric_balance_statistics()
    
        stage = domain.get_quantity('stage')

        #for i in range(samples):
        #    if tri_ids[i] >= 0:                
        #        if verbose: print 'P%d tri %d, value = %s' %(myid, i, stage.centroid_values[tri_ids[i]])
                    
        sys.stdout.flush()
 
        pass

 
    domain.sww_merge(delete_old=True)
    
    success = True

##-----------------------------------------------------------------------
## Assign/Test Control data
##-----------------------------------------------------------------------

    if not parallel:
        stage = domain.get_quantity('stage')

        for i in range(samples):
            assert(tri_ids[i] >= 0)
            control_data.append(stage.centroid_values[tri_ids[i]])

        if inlet0 is not None:
            control_data.append(inlet0.inlet.get_average_stage())
            control_data.append(inlet0.inlet.get_average_xmom())
            control_data.append(inlet0.inlet.get_average_ymom())
            control_data.append(inlet0.inlet.get_total_water_volume())
            control_data.append(inlet0.inlet.get_average_depth())

        if verbose: print('P%d control_data = %s' %(myid, control_data))
    else:
        stage = domain.get_quantity('stage')
        
        for i in range(samples):
            if tri_ids[i] >= 0:
                local_success = num.allclose(control_data[i], stage.centroid_values[tri_ids[i]])
                success = success and local_success
                if verbose and not local_success: 
                    print('P%d tri %d, control = %s, actual = %s, Success = %s' %(myid, i, control_data[i], stage.centroid_values[tri_ids[i]], local_success))
                    sys.stdout.flush()
                    
                
                
        if inlet0 is not None:
            inlet_master_proc = inlet0.inlet.get_master_proc()
            average_stage = inlet0.inlet.get_global_average_stage()
            average_xmom = inlet0.inlet.get_global_average_xmom()
            average_ymom = inlet0.inlet.get_global_average_ymom()
            average_volume = inlet0.inlet.get_global_total_water_volume()
            average_depth = inlet0.inlet.get_global_average_depth()

            if myid == inlet_master_proc:
                if verbose: 
                    print('P%d average stage, control = %s, actual = %s' %(myid, control_data[samples], average_stage))

                    print('P%d average xmom, control = %s, actual = %s' %(myid, control_data[samples+1], average_xmom))

                    print('P%d average ymom, control = %s, actual = %s' %(myid, control_data[samples+2], average_ymom))

                    print('P%d average volume, control = %s, actual = %s' %(myid, control_data[samples+3], average_volume))

                    print('P%d average depth, control = %s, actual = %s' %(myid, control_data[samples+4], average_depth))


        #assert(success)

    return control_data


if __name__=="__main__":
    
    test_points = []

    if myid == 0:

        for i in range(samples):
            x = random.randrange(0,1000)/1000.0 * length
            y = random.randrange(0,1000)/1000.0 * width
            point = [x, y]
            test_points.append(point)
        
        for i in range(1,numprocs):
            pypar.send(test_points, i)
    else:
        test_points = pypar.receive(0)

    print("Test Points::")
    print(test_points)

    if myid == 0:
        control_data = run_simulation(parallel=False, test_points = test_points, verbose = True)
        
        for proc in range(1,numprocs):
            pypar.send(control_data, proc)
    else:
        control_data = pypar.receive(0)

    pypar.barrier()
    run_simulation(parallel=True, control_data = control_data, test_points = test_points, verbose = True)


finalize()
