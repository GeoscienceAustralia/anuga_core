import os.path
   
from anuga import rectangular_cross
import anuga

from anuga import Boyd_box_operator
from anuga import Inlet_operator
                            

import numpy as num

  

length = 40.
width = 15.

dx = dy = 0.5          # Resolution: Length of subdivisions on both axes

points, vertices, boundary = rectangular_cross(int(length/dx),
                                               int(width/dy),
                                               len1=length, 
                                               len2=width)
domain = anuga.Domain(points, vertices, boundary)   
domain.set_name()                 # Output name



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


domain.set_quantity('elevation', topography) 
domain.set_quantity('friction', 0.01)         # Constant friction 
domain.set_quantity('stage',
                    expression='elevation')   # Dry initial condition



gate = Boyd_box_operator(domain,
                            end_points=[[9.0, 2.5],[13.0, 2.5]],
                            losses=1.5,
                            width=1.5,
                            height = 10.0,
                            apron=5.0,
                            use_momentum_jet=True,
                            use_velocity_head=False,
                            manning=0.013,
                            verbose=False)

line = [[0.0, 5.0], [0.0, 10.0]]
Q = 1.0
Inlet_operator(domain, line, Q)




##-----------------------------------------------------------------------
## Setup boundary conditions
##-----------------------------------------------------------------------

## Inflow based on Flow Depth and Approaching Momentum
Br = anuga.Reflective_boundary(domain)              # Solid reflective wall

domain.set_boundary({'left': Br, 'right': Br, 'top': Br, 'bottom': Br})


##-----------------------------------------------------------------------
## Evolve system through time
##-----------------------------------------------------------------------

#min_delta_w = sys.maxint 
#max_delta_w = -min_delta_w
for t in domain.evolve(yieldstep = 1.0, finaltime = 50):
    domain.write_time()

    if num.allclose(t, 10.0):
        gate.set_culvert_height(0.0)

    Q, velocity, depth = gate.discharge_routine()

    print gate.inlets[0].get_enquiry_stage()
    print gate.inlets[1].get_enquiry_stage()
    print gate.get_culvert_height()

    print Q
    print velocity
    print depth
        
    gate.print_timestepping_statistics()
    




