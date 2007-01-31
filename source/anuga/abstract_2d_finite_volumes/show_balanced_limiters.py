"""Example of shallow water wave equation.

Specific methods pertaining to the 2D shallow water equation
are imported from shallow_water
for use with the generic finite volume framework

Conserved quantities are h, uh and vh stored as elements 0, 1 and 2 in the
numerical vector named conserved_quantities.
"""

######################
# Module imports 
#
from anuga.shallow_water import Domain,\
     Reflective_boundary, Dirichlet_boundary,\
     Transmissive_boundary, Time_boundary
from anuga.shallow_water.shallow_water_domain import Weir_simple as Weir

from mesh_factory import rectangular
from Numeric import array
    

######################
# Domain
#

N = 12

print 'Creating domain'
#Create basic mesh
points, vertices, boundary = rectangular(N, N/2, len1=1.2,len2=0.6,
                                         origin=(-0.07, 0))

print 'Number of elements', len(vertices)
#Create shallow water domain
domain = Domain(points, vertices, boundary)
domain.smooth = False
domain.default_order = 2
domain.set_name('show_balanced_limiters')
domain.store = True
domain.format = 'sww'   #Native netcdf visualisation format

#Set bed-slope and friction
inflow_stage = 0.1
manning = 0.1
Z = Weir(inflow_stage)

print 'Field values'
domain.set_quantity('elevation', Z)
domain.set_quantity('friction', manning)


######################
# Boundary conditions
#
print 'Boundaries'
Br = Reflective_boundary(domain)
Bt = Transmissive_boundary(domain)

#Constant inflow
Bd = Dirichlet_boundary([inflow_stage, 0.0, 0.0])

#Time dependent inflow
from math import sin, pi
Bw = Time_boundary(domain=domain,
                   f=lambda x: [(1 + sin(x*pi/4))*\
                                (inflow_stage*(sin(2.5*x*pi)+0.7)),0,0])

#Set boundary conditions
domain.set_boundary({'left': Bd, 'right': Br, 'bottom': Br, 'top': Br})

                    

######################
#Initial condition
#
print 'Initial condition'
domain.set_quantity('stage', Z)

from Numeric import allclose

#Evolve
for t in domain.evolve(yieldstep = 0.1, finaltime = 30):
    domain.write_time(track_location=True)
    domain.write_boundary_statistics(['stage'],'left')

print 'Done'    
    

