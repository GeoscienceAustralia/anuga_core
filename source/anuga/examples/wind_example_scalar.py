"""Example of shallow water wave equation.

Flat bed with constant wind stress
"""

######################
# Module imports 
from mesh_factory import rectangular
from shallow_water import Domain, Reflective_boundary, Constant_height, Wind_stress

#Create basic mesh (100m x 100m)
N = 100     
len = 100
points, vertices, boundary = rectangular(N, N, len, len)

#Create shallow water domain
domain = Domain(points, vertices, boundary)
domain.default_order = 2
domain.store = True
domain.set_name('wind_laminar')

#Set initial conditions
domain.set_quantity('elevation', 0.0)
domain.set_quantity('stage', 1.0)
domain.set_quantity('friction', 0.03)

#Constant (quite extreme :-) windfield: 9000 m/s, bearing 120 degrees
domain.forcing_terms.append( Wind_stress(s=9000, phi=120) )

######################
# Boundary conditions
Br = Reflective_boundary(domain)
domain.set_boundary({'left': Br, 'right': Br, 'top': Br, 'bottom': Br})

######################
#Evolution
for t in domain.evolve(yieldstep = 0.5, finaltime = 1000):
    domain.write_time()

print 'Done'    

