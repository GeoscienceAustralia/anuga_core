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
from shallow_water import Domain, Reflective_boundary, Dirichlet_boundary,\
     Transmissive_boundary, Time_boundary, Constant_height, Weir

from mesh_factory import rectangular
from pmesh2domain import pmesh_to_domain

from Numeric import array

######################
# Domain

import sys
if len(sys.argv) > 1: 
    filename = sys.argv[1]
else:
    filename = 'weir_domain_refined.tsh'

print 'Creating domain from', filename
domain_list = pmesh_to_domain(filename)
vertex_coordinates = domain_list[0]
volumes = domain_list[1]
marker_dict = domain_list[2]
vertex_quantity_dict = domain_list[3]

domain = Domain(vertex_coordinates, volumes, marker_dict)
print "Number of triangles = ", len(domain)

domain.store = False #True
domain.format = 'sww'
domain.filename = 'weir'
domain.checkpoint = False #True
domain.visualise = True #False
domain.default_order = 2

#Set bed-slope and friction
inflow_stage = 0.15
manning = 0.07
W = Weir(inflow_stage)

print 'Field values'

domain.set_quantity('elevation', W)
domain.set_quantity('friction', manning)



######################
# Boundary conditions
#
print 'Boundaries'
Br = Reflective_boundary(domain)
Bt = Transmissive_boundary(domain)

#Constant inflow
Bd = Dirichlet_boundary(array([inflow_stage, 0.0, 0.0]))

#Time dependent inflow
from math import sin, pi
Bw = Time_boundary(domain=domain,
                   f=lambda x: array([(1 + sin(x*pi/4))*\
                                      (inflow_stage*(sin(2.5*x*pi)+0.7)),0,0]))


print 'Available boundary tags are', domain.get_boundary_tags()

#Set boundary conditions
domain.set_boundary({'left': Bw, '0': Br, '1':Bw, 'external':Br})


#print domain.quantities['elevation'].vertex_values
#print domain.quantities['stage'].vertex_values
         
######################
#Initial condition
print 'Initial condition'
domain.set_quantity('stage', Constant_height(W, 0.))
domain.check_integrity()

######################
#Evolution
for t in domain.evolve(yieldstep = 0.01, finaltime = 5):
    domain.write_time()
    
print 'Done'

    
