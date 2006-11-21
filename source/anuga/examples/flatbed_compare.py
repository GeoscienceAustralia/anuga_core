"""Example of shallow water wave equation.

Generate slope

"""

######################
# Module imports 
#
from mesh_factory import rectangular
from shallow_water import Domain, Reflective_boundary, Dirichlet_boundary,\
     Constant_height
from Numeric import array

#Create basic mesh
N = 8
points, vertices, boundary = rectangular(N, N)

#Create shallow water domain
domain = Domain(points, vertices, boundary)
domain.smooth = False
domain.default_order=2


######################
# Boundary conditions
Br = Reflective_boundary(domain)
Bd = Dirichlet_boundary([0.2,0.,0.])

domain.set_boundary({'left': Bd, 'right': Br, 'top': Br, 'bottom': Br})


domain.check_integrity()

######################
#Evolution
#for t in domain.evolve(yieldstep = 0.05, finaltime = 1.):
#    domain.write_time()


#import sys; sys.exit()

######################
#Evolution
for t in domain.evolve(yieldstep = 0.01, finaltime = 0.03):
    domain.write_time()

print domain.quantities['stage'].centroid_values[:4]
print domain.quantities['xmomentum'].centroid_values[:4]
print domain.quantities['ymomentum'].centroid_values[:4]
domain.distribute_to_vertices_and_edges()
print
print domain.quantities['stage'].vertex_values[:4,0]
print domain.quantities['xmomentum'].vertex_values[:4,0]
print domain.quantities['ymomentum'].vertex_values[:4,0]  
