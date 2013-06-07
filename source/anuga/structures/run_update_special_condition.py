import os.path
import sys

from anuga.utilities.system_tools import get_pathname_from_package
from anuga.geometry.polygon_function import Polygon_function
        
from anuga.abstract_2d_finite_volumes.mesh_factory import rectangular_cross
from anuga.abstract_2d_finite_volumes.quantity import Quantity

import anuga

from anuga.structures.boyd_pipe_operator import Boyd_pipe_operator
from anuga.structures.inlet_operator import Inlet_operator
                            
#from anuga.culvert_flows.culvert_routines import boyd_generalised_culvert_model
     
from math import pi, pow, sqrt

import numpy as num





"""test_that_culvert_runs_rating

This test exercises the culvert and checks values outside rating curve
are dealt with       
"""

path = get_pathname_from_package('anuga.culvert_flows')    

length = 40.
width = 15.

dx = dy = 0.5          # Resolution: Length of subdivisions on both axes

points, vertices, boundary = rectangular_cross(int(length/dx),
                                               int(width/dy),
                                               len1=length, 
                                               len2=width)
domain = anuga.Domain(points, vertices, boundary)   
domain.set_name('Test_culvert')                 # Output name
domain.set_default_order(2)
#domain.set_beta(1.5)


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


line = [[0.01, 5.0], [0.01, 10.0]]
Q = 50.0

import anuga

inlet = anuga.structures.inlet.Inlet(domain, line)

domain.special_inlet = inlet

def my_update_special_conditions(domain):

    #print 'inside my domain update_special condition'

    domain.special_inlet.set_depths(5.0)

import types
domain.update_special_conditions = types.MethodType(my_update_special_conditions, domain)

domain.set_quantity('elevation', topography) 
domain.set_quantity('friction', 0.01)         # Constant friction 
domain.set_quantity('stage',
                    expression='elevation')   # Dry initial condition

filename=os.path.join(path, 'example_rating_curve.csv')



#Boyd_pipe_operator(domain,
#                            end_point0=[9.0, 2.5],
#                            end_point1=[13.0, 2.5],
#                            losses=1.5,
#                            diameter=1.5,
#                            apron=5.0,
#                            use_momentum_jet=True,
#                            use_velocity_head=False,
#                            manning=0.013,
#                            verbose=False)


line = [[0.0, 5.0], [0.0, 10.0]]
Q = 50.0
#Inlet_operator(domain, line, Q)







##-----------------------------------------------------------------------
## Setup boundary conditions
##-----------------------------------------------------------------------

## Inflow based on Flow Depth and Approaching Momentum
Bi = anuga.Dirichlet_boundary([2.0, 0.0, 0.0])
Br = anuga.Reflective_boundary(domain)              # Solid reflective wall
#Bo = anuga.Dirichlet_boundary([-5, 0, 0])           # Outflow

## Upstream and downstream conditions that will exceed the rating curve
## I.e produce delta_h outside the range [0, 10] specified in the the 
## file example_rating_curve.csv
#Btus = anuga.Time_boundary(domain, \
            #lambda t: [100*num.sin(2*pi*(t-4)/10), 0.0, 0.0])
#Btds = anuga.Time_boundary(domain, \
            #lambda t: [-5*(num.cos(2*pi*(t-4)/20)), 0.0, 0.0])
#domain.set_boundary({'left': Btus, 'right': Btds, 'top': Br, 'bottom': Br})
domain.set_boundary({'left': Br, 'right': Br, 'top': Br, 'bottom': Br})


#--------------------------------------------------------------------------
# Turn on the visualisation
#--------------------------------------------------------------------------
visualise = True
if visualise:
    from anuga.visualiser import RealtimeVisualiser
    vis = RealtimeVisualiser(domain)
    vis.render_quantity_height("elevation", offset=0.001, dynamic=False)
    vis.render_quantity_height("stage", dynamic=True)
    vis.colour_height_quantity('stage', (0.2, 0.2, 0.8))
    vis.start()
    import time
    time.sleep(2.0)



##-----------------------------------------------------------------------
## Evolve system through time
##-----------------------------------------------------------------------

#min_delta_w = sys.maxint 
#max_delta_w = -min_delta_w
for t in domain.evolve(yieldstep = 1.0, finaltime = 30):
    domain.write_time()

    #if domain.get_time() > 150.5 and domain.get_time() < 151.5 :
        #Bi = anuga.Dirichlet_boundary([0.0, 0.0, 0.0])
        #domain.set_boundary({'left': Bi, 'right': Br, 'top': Br, 'bottom': Br})

    #delta_w = culvert.inlet.stage - culvert.outlet.stage
    
    #if delta_w > max_delta_w: max_delta_w = delta_w
    #if delta_w < min_delta_w: min_delta_w = delta_w

    #print domain.volumetric_balance_statistics()

    if visualise:
        vis.update()

    pass



if visualise: vis.evolveFinished()


## Check that extreme values in rating curve have been exceeded
## so that we know that condition has been exercised
#assert min_delta_w < 0
#assert max_delta_w > 10        


#os.remove('Test_culvert.sww')