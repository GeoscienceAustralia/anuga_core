"""Example of shallow water wave equation
using time boundary from a file
"""

######################
# Module imports 
#
from mesh_factory import rectangular
from shallow_water import Domain, Reflective_boundary, Dirichlet_boundary,\
     Constant_height, Time_boundary, File_boundary
from Numeric import array, ones

#Create basic mesh (100m x 100m)
points, vertices, boundary = rectangular(50, 50, 100, 100)

#Create shallow water domain
domain = Domain(points, vertices, boundary)
domain.smooth = False
domain.visualise = False
domain.default_order = 2
domain.store = True     #Store for visualisation purposes
domain.format = 'sww'   #Native netcdf visualisation format



#######################
#Bed-slope and friction
domain.set_quantity('elevation', 0.0) #Flat
#domain.set_quantity('elevation', lambda x,y: -(x+y)/3) #Slopy
domain.set_quantity('friction', 0.1)


######################
# Boundary conditions


#Write file
import os, time
from anuga.config import time_format
from math import sin, pi

filename = 'Eden_Australia_31082004.txt'

Br = Reflective_boundary(domain)
Bd = Dirichlet_boundary([1. ,0.,0.])
Bw = Time_boundary(domain=domain,
                   f=lambda t: [(1*sin(t*pi/20)), 0.0, 0.0])

Bf = File_boundary(filename, domain)

#domain.set_boundary({'left': Bw, 'right': Br, 'top': Br, 'bottom': Br})
domain.set_boundary({'left': Bf, 'right': Br, 'top': Br, 'bottom': Br})


######################
#Initial condition (constant)
domain.set_quantity('stage', 0.0)
domain.check_integrity()

  
######################
#Evolution
for t in domain.evolve(yieldstep = 60, finaltime = 385*15*60):
    domain.write_time()

