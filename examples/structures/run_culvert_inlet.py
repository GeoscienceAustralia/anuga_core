import os.path


from anuga.utilities.system_tools import get_pathname_from_package   
from anuga import rectangular_cross
import anuga

from anuga import Boyd_pipe_operator
from anuga import Inlet_operator
                            
     
from math import pi, pow, sqrt


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


Boyd_pipe_operator(domain,
                            end_points=[[9.0, 2.5],[13.0, 2.5]],
                            losses=1.5,
                            diameter=1.5,
                            apron=5.0,
                            use_momentum_jet=True,
                            use_velocity_head=False,
                            manning=0.013,
                            verbose=False)


line = [[0.0, 5.0], [0.0, 10.0]]
Q = 5.0
Inlet_operator(domain, line, Q)




##-----------------------------------------------------------------------
## Setup boundary conditions
##-----------------------------------------------------------------------

## Inflow based on Flow Depth and Approaching Momentum
Bi = anuga.Dirichlet_boundary([2.0, 0.0, 0.0])
Br = anuga.Reflective_boundary(domain)              # Solid reflective wall

domain.set_boundary({'left': Br, 'right': Br, 'top': Br, 'bottom': Br})


##-----------------------------------------------------------------------
## Evolve system through time
##-----------------------------------------------------------------------

#min_delta_w = sys.maxint 
#max_delta_w = -min_delta_w
for t in domain.evolve(yieldstep = 1.0, finaltime = 200):
    domain.write_time()

    print domain.volumetric_balance_statistics()
    
    pass


