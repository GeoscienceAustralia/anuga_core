"""Simple water flow example using ANUGA

Water driven up a linear slope and time varying boundary,
similar to a beach environment
"""

#------------------------------------------------------------------------------
# Import necessary modules
#------------------------------------------------------------------------------

import sys
import anuga


from math import cos
from numpy import zeros, float, where
import numpy
from time import localtime, strftime, gmtime


#------------------------------------------------------------------------------
# Setup domain
#------------------------------------------------------------------------------
dx = 1000.
dy = dx
L = 100000.
W = 10*dx
#W = dx

# structured mesh
domain = anuga.rectangular_cross_domain(int(L/dx), int(W/dy), L, W, (0.0, -W/2))

domain.set_name() # based on script name

print (domain.mesh.statistics(nbins=50))

#------------------------------------------------------------------------------
# Setup Algorithm
#------------------------------------------------------------------------------
domain.set_flow_algorithm('2_0')

#------------------------------------------------------------------------------
# Setup initial conditions
#------------------------------------------------------------------------------
h0 = 10.0
h1 = 0.0

def elevation(x,y):

    return where(x<=95000.0, 0.0, h0)

def height(x,y):

    return where(x<=50000.0, h0, h1)

domain.set_quantity('elevation',elevation)
domain.set_quantity('friction', 0.0)
domain.add_quantity('stage', height)

#-----------------------------------------------------------------------------
# Setup boundary conditions
#------------------------------------------------------------------------------
from math import sin, pi, exp
Br = anuga.Reflective_boundary(domain)      # Solid reflective wall
Bt = anuga.Transmissive_boundary(domain)    # Continue all values on boundary 
Bd = anuga.Dirichlet_boundary([1,0.,0.])    # Constant boundary values


# Associate boundary tags with boundary objects
domain.set_boundary({'left': Br, 'right': Br, 'top': Br, 'bottom': Br})

#-----------------------------------------------------------------------------
# Setup operators that will be applied each inner timestep
#------------------------------------------------------------------------------
from anuga.operators.collect_max_stage_operator import Collect_max_stage_operator
max_operator = Collect_max_stage_operator(domain)

#------------------------------------------------------------------------------
# Evolve system through time
#------------------------------------------------------------------------------
for t in domain.evolve(yieldstep = 100.0, finaltime = 60*60.):
    #print domain.timestepping_statistics(track_speeds=True)
    domain.print_timestepping_statistics()
    domain.print_operator_timestepping_statistics()


# save the max_stage centroid data to a text file
max_operator.save_centroid_data_to_csv()

# Let's have a look at the max_stage
#max_operator.plot_quantity()




