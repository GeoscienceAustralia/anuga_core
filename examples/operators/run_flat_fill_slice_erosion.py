"""Simple water flow example using ANUGA

Water flowing down a channel with a topography that varies with time
"""

#------------------------------------------------------------------------------
# Import necessary modules
#------------------------------------------------------------------------------
from anuga import rectangular_cross
from anuga import Domain
from anuga import Reflective_boundary
from anuga import Dirichlet_boundary
from anuga import Time_boundary

import numpy as num

#===============================================================================
# Setup Functions
#===============================================================================

#------------------------------------------------------------------------------
# Setup initial conditions
#------------------------------------------------------------------------------
def topography(x,y):
    """Complex topography defined by a function of vectors x and y."""
    print (' Create topography....')
    z = -x/100

    N = len(x)
    for i in range(N):
        # Step
        if 2 < x[i] < 4:
            z[i] += 0.4 - 0.05*y[i]

        # Permanent pole
        if (x[i] - 8)**2 + (y[i] - 2)**2 < 0.4**2:
            z[i] += 1

        # Dam
        #if 12 < x[i] < 13:
        #    z[i] += 1.4
    # Sloping Embankment Across Channel
        if 10.6 < x[i] < 12.1:
            # Cut Out Segment for Culvert face                
            z[i] +=  1.0*(x[i] -10.6)    # Sloping Segment  U/S Face
        if 12.0 < x[i] < 13.0:
           z[i] +=  1.4                    # Flat Crest of Embankment
        if 12.9 < x[i] < 13.7:
            z[i] +=  1.4-2.0*(x[i] -12.9) # Sloping D/S Face
            
    return z



#===============================================================================
# Setup and Run Model
#===============================================================================


#------------------------------------------------------------------------------
# Setup computational domain
#------------------------------------------------------------------------------
print (' Set up Domain first...')
length = 24.
width = 5.
dx = dy = 0.2 #.1           # Resolution: Length of subdivisions on both axes

points, vertices, boundary = rectangular_cross(int(length/dx), int(width/dy),
                                               len1=length, len2=width)
domain = Domain(points, vertices, boundary)
domain.set_flow_algorithm('DE1')
domain.set_name() # Output name based on script name
print (domain.statistics())

domain.set_quantities_to_be_stored({'elevation': 2,
                                    'stage': 2,
                                    'xmomentum': 2,
                                    'ymomentum': 2})



domain.set_quantity('elevation', topography)           # elevation is a function
domain.set_quantity('friction', 0.01)                  # Constant friction
domain.set_quantity('stage', expression='elevation')   # Dry initial condition

#------------------------------------------------------------------------------
# Setup boundary conditions
#------------------------------------------------------------------------------
Bi = Dirichlet_boundary([1.5, 0, 0])          # Inflow
Br = Reflective_boundary(domain)              # Solid reflective wall
Bo = Dirichlet_boundary([-5, 0, 0])           # Outflow

domain.set_boundary({'left': Bi, 'right': Bo, 'top': Br, 'bottom': Br})

#------------------------------------------------------------------------------
# Setup erosion operator in the middle of dam
#------------------------------------------------------------------------------
print ('Set up Erosion Area to test...')

from anuga import Flat_fill_slice_erosion_operator as Default_erosion_operator
from anuga import Set_elevation_operator as Default_erosion_operator

polygon1 = [ [10.6, 1.0], [13.7, 1.0], [13.7, 4.0], [10.6, 4.0] ]

# setup interpolated function
time = num.array([ 5.0 , 8.0])
elev = num.array([ 1.4,  2.0])

from scipy.interpolate import interp1d
flevel = interp1d(time,elev)

# create operator
op1 = Default_erosion_operator(domain, elevation=flevel, polygon=polygon1)





#------------------------------------------------------------------------------
# Evolve system through time
#------------------------------------------------------------------------------
for t in domain.evolve(yieldstep=0.1, finaltime=15.0):
    domain.print_timestepping_statistics()
    #domain.print_operator_timestepping_statistics()










