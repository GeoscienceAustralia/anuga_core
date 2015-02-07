import os.path
import sys


import anuga

from anuga import myid, numprocs, finalize
from anuga import Inlet_operator, Boyd_box_operator




"""
This test exercises the parallel culvert
"""
verbose = True

    

length = 40.
width = 16.

dx = dy = 2           # Resolution: Length of subdivisions on both axes

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

##-----------------------------------------------------------------------
## Setup domain
##-----------------------------------------------------------------------
if myid == 0:
    points, vertices, boundary = anuga.rectangular_cross(int(length/dx),
                                                   int(width/dy),
                                                   len1=length, 
                                                   len2=width)

    domain = anuga.Domain(points, vertices, boundary)   
    domain.set_name()                 # Output name output_script_name.sww
    domain.set_flow_algorithm('1_5')
else:
    domain = None
    
##-----------------------------------------------------------------------
## Distribute domain
##-----------------------------------------------------------------------


domain = anuga.distribute(domain)
#domain.dump_triangulation("run_parallel_boyd_box_op_domain.png")
    

##-----------------------------------------------------------------------
## Setup boundary conditions
##-----------------------------------------------------------------------
    
domain.set_quantity('elevation', topography) 
domain.set_quantity('friction', 0.01)         # Constant friction 
domain.set_quantity('stage',
                        expression='elevation')   # Dry initial condition

    
Bi = anuga.Dirichlet_boundary([5.0, 0.0, 0.0])
Br = anuga.Reflective_boundary(domain)              # Solid reflective wall
domain.set_boundary({'left': Br, 'right': Br, 'top': Br, 'bottom': Br})



    ################ Define Fractional Operators ##########################
line0 = [[10.0, 10.0], [30.0, 10.0]]
#line0 = [[29.0, 10.0], [30.0, 10.0]]
poly1 = [[0.0, 10.0], [0.0, 15.0], [5.0, 15.0], [5.0, 10.0]]
Q0 = anuga.file_function('../data/test_hydrograph.tms', quantities=['hydrograph'])
Q1 = 5.0

samples = 50


inlet0 = None
inlet1 = None
boyd_box0 = None

inlet0 = Inlet_operator(domain, line0, Q0, logging=True, description='inlet0', verbose = False)
inlet1 = Inlet_operator(domain, poly1, Q1, logging=True, description='inlet1', verbose = False)

# Enquiry point [ 19.    2.5] is contained in two domains in 4 proc case

boyd_box0 = Boyd_box_operator(domain,
                              end_points=[[9.0, 2.5],[19.0, 2.5]],
                              losses=1.5,
                              width=5.0,
                              apron=5.0,
                              use_momentum_jet=True,
                              use_velocity_head=False,
                              manning=0.013,
                              logging=True,
                              description='boyd_box_0',
                              verbose=False)

if inlet0 is not None and verbose: inlet0.print_statistics()
if inlet1 is not None and verbose: inlet1.print_statistics()
if boyd_box0 is not None and verbose: boyd_box0.print_statistics()


##-----------------------------------------------------------------------
## Evolve system through time
##-----------------------------------------------------------------------

for t in domain.evolve(yieldstep = 2.0, finaltime = 20.0):
    if myid == 0 and verbose:
        domain.write_time()

    #print domain.volumetric_balance_statistics()

    stage = domain.get_quantity('stage')


    if boyd_box0 is not None and verbose : boyd_box0.print_timestepping_statistics()

    #for i in range(samples):
    #    if tri_ids[i] >= 0:                
    #        if verbose: print 'P%d tri %d, value = %s' %(myid, i, stage.centroid_values[tri_ids[i]])

    sys.stdout.flush()

    pass

##-----------------------------------------------------------------------
## Assign/Test Control data
##-----------------------------------------------------------------------

domain.sww_merge(delete_old=True)

finalize()

    

