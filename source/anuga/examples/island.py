"""Example of shallow water wave equation.

Island surrounded by water.
This example investigates onshore 'creep'

"""


#------------------------------------------------------------------------------
# Import necessary modules
#------------------------------------------------------------------------------

# Standard modules
from math import exp

# Application specific imports
from anuga.abstract_2d_finite_volumes.mesh_factory import rectangular
from anuga.shallow_water import Domain, Reflective_boundary, Dirichlet_boundary
from anuga.pmesh.mesh_interface import create_mesh_from_regions
from anuga.utilities.polygon import Polygon_function, read_polygon

from anuga.abstract_2d_finite_volumes.quantity import Quantity
from Numeric import allclose

#------------------------------------------------------------------------------
# Setup computational domain
#------------------------------------------------------------------------------

#Create basic mesh
create_mesh_from_regions( [[0,0], [100,0], [100,100], [0,100]],
                          boundary_tags = {'bottom': [0],
                                           'right': [1],
                                           'top': [2],
                                           'left': [3]},
                          maximum_triangle_area = 10.0,
                          filename = 'island.msh' ,
                          interior_regions=[ ([[50,25], [70,25], [70,75], [50,75]], 10.0)]
                          #interior_holes=[[[50,25], [70,25], [70,75], [50,75]]],
                          )



#Create shallow water domain
domain = Domain(mesh_filename = 'island.msh')
domain.smooth = False
domain.set_name('island')
domain.default_order = 2


#I tried to introduce this parameter top control the h-limiter,
#but it doesn't remove the 'lapping water'
#NEW (ON): I think it has been fixed (or at least reduced significantly) now
#
# beta_h == 1.0 means that the largest gradients (on h) are allowed
# beta_h == 0.0 means that constant (1st order) gradients are introduced
# on h. This is equivalent to the constant depth used previously.
#domain.beta_h     = 0.5
#domain.beta_w_dry = 0.0
#domain.alpha_balance = 10.0
#domain.minimum_allowed_height = 1.0e-4 
#domain.maximum_allowed_speed = 100.0 
#domain.minimum_storable_height = 1.0e-4 

domain.beta_h     = 0.0
domain.limit2007 = 1

#------------------------------------------------------------------------------
# Setup initial conditions
#------------------------------------------------------------------------------

def island(x, y):
    z = 0*x
    for i in range(len(x)):
        z[i] = 20*exp( -((x[i]-50)**2 + (y[i]-50)**2)/100 )

        #z[i] += 0.5*exp( -((x[i]-10)**2 + (y[i]-10)**2)/50 )

    return z

def slump(x, y):
    z = 0*x
    for i in range(len(x)):
        z[i] -= 0.7*exp( -((x[i]-10)**2 + (y[i]-10)**2)/200 )

    return z

stage_value = 15.0
#domain.set_quantity('friction', 0.1)  #Honky dory
domain.set_quantity('friction', 0.01)     #Creep
domain.set_quantity('elevation', island)
domain.set_quantity('stage', stage_value)
domain.max_timestep = 0.01



#------------------------------------------------------------------------------
# Setup boundary conditions (all reflective)
#------------------------------------------------------------------------------

Br = Reflective_boundary(domain)
Bd = Dirichlet_boundary([stage_value, 0, 0])

domain.set_boundary({'left': Bd, 'right': Bd, 'top': Bd, 'bottom': Bd, 'exterior': Br})
domain.check_integrity()

#------------------------------------------------------------------------------
# Evolve system through time
#------------------------------------------------------------------------------

import time
for t in domain.evolve(yieldstep = 1, finaltime = 10):
    domain.write_time()
    #if allclose(t, 100):
    #    Q = domain.get_quantity('stage')
    #    Q_s = Quantity(domain)
    #    Q_s.set_values(slump)
    #    domain.set_quantity('stage', Q + Q_s)
    #print '    Volume: ', domain.get_quantity('stage').get_integral()

print 'Done'
