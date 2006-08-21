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
from pyvolution.mesh_factory import rectangular
from pyvolution.shallow_water import Domain, Reflective_boundary, Dirichlet_boundary
from pmesh.mesh_interface import create_mesh_from_regions
from utilities.polygon import Polygon_function, read_polygon

from pyvolution.quantity import Quantity
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
                          maximum_triangle_area = 5,
                          filename = 'island.msh' ,
                          interior_regions=[ ([[50,25], [70,25], [70,75], [50,75]], 3)]
                          )



#Create shallow water domain
domain = Domain(mesh_filename = 'island.msh')
domain.smooth = False
domain.set_name('island')
domain.set_default_order(2)


#I tried to introduce this parameter top control the h-limiter,
#but it doesn't remove the 'lapping water'
#NEW (ON): I think it has been fixed (or at least reduced significantly) now
#
# beta_h == 1.0 means that the largest gradients (on h) are allowed
# beta_h == 0.0 means that constant (1st order) gradients are introduced
# on h. This is equivalent to the constant depth used previously.
domain.beta_h = 0.9


#------------------------------------------------------------------------------
# Setup initial conditions
#------------------------------------------------------------------------------

def island(x, y):
    z = 0*x
    for i in range(len(x)):
        z[i] = 8*exp( -((x[i]-50)**2 + (y[i]-50)**2)/100 )

        #z[i] += 0.5*exp( -((x[i]-10)**2 + (y[i]-10)**2)/50 )

    return z

def slump(x, y):
    z = 0*x
    for i in range(len(x)):
        z[i] -= 0.7*exp( -((x[i]-10)**2 + (y[i]-10)**2)/200 )

    return z

#domain.set_quantity('friction', 0.1)  #Honky dory
domain.set_quantity('friction', 100)     #Creep
domain.set_quantity('elevation', island)
domain.set_quantity('stage', 1)
domain.max_timestep = 0.01


#------------------------------------------------------------------------------
# Setup boundary conditions (all reflective)
#------------------------------------------------------------------------------

Br = Reflective_boundary(domain)
Bd = Dirichlet_boundary([1, 0, 0])

domain.set_boundary({'left': Bd, 'right': Bd, 'top': Bd, 'bottom': Bd})
domain.check_integrity()


#------------------------------------------------------------------------------
# Evolve system through time
#------------------------------------------------------------------------------

import time
for t in domain.evolve(yieldstep = 1, finaltime = 100):
    domain.write_time()
    #if allclose(t, 100):
    #    Q = domain.get_quantity('stage')
    #    Q_s = Quantity(domain)
    #    Q_s.set_values(slump)
    #    domain.set_quantity('stage', Q + Q_s)
    #print '    Volume: ', domain.get_quantity('stage').get_integral()

print 'Done'
