import os.path
import sys

from anuga.utilities.system_tools import get_pathname_from_package

        
from anuga import rectangular_cross
from anuga import file_function

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


path = get_pathname_from_package('anuga.structures') 
filename=os.path.join(path, 'tests', 'data', 'test_hydrograph.tms')


Q = file_function(filename, quantities=['hydrograph'])

line1 = [[1.0, 5.0], [2.0, 10.0]]
poly1 = [[1.0, 5.0], [2.0, 5.0], [2.0, 10.0], [1.0, 10.0]]
inlet1 = Inlet_operator(domain, poly1, Q)



##-----------------------------------------------------------------------
## Setup boundary conditions
##-----------------------------------------------------------------------

## Inflow based on Flow Depth and Approaching Momentum
Bi = anuga.Dirichlet_boundary([2.0, 0.0, 0.0])
Br = anuga.Reflective_boundary(domain)              # Solid reflective wall
#Bo = anuga.Dirichlet_boundary([-5, 0, 0])           # Outflow


domain.set_boundary({'left': Br, 'right': Br, 'top': Br, 'bottom': Br})


##-----------------------------------------------------------------------
## Evolve system through time
##-----------------------------------------------------------------------

for t in domain.evolve(yieldstep = 1.0, finaltime = 38):
    domain.write_time()
    print domain.volumetric_balance_statistics()
    inlet1.print_timestepping_statistics()
    pass


