"""Simple example of shallow water wave equation using Abstract_2d_finite_volumes

Water driven up a linear slope with a time varying boundary,
similar to beach environment

"""

######################
# Module imports 
from anuga.abstract_2d_finite_volumes.mesh_factory import rectangular
from anuga.abstract_2d_finite_volumes.shallow_water import Domain, Reflective_boundary,\
     Dirichlet_boundary, Time_boundary, Transmissive_boundary

#Create basic triangular mesh
points, vertices, boundary = rectangular(10, 10)

#Create shallow water domain
domain = Domain(points, vertices, boundary)
domain.set_name('bedslope')
domain.set_datadir('.')                      #Use current directory for output
domain.set_quantities_to_be_stored('stage')  #See shallow_water.py

print domain.statistics()

#######################
# Initial conditions
def f(x,y):
    return -x/2

domain.set_quantity('elevation', f)
domain.set_quantity('friction', 0.1)

h = 0.00               # Constant depth over elevation
domain.set_quantity('stage', -.4)


######################
# Boundary conditions
from math import sin, pi
Br = Reflective_boundary(domain)
Bt = Transmissive_boundary(domain)
Bd = Dirichlet_boundary([0.2,0.,0.])
Bw = Time_boundary(domain=domain,
                   f=lambda t: [(0.1*sin(t*2*pi)-0.3), 0.0, 0.0])


print 'Tags are ', domain.get_boundary_tags()

domain.set_boundary({'left': Br, 'right': Bw, 'top': Br, 'bottom': Br})


######################
#Evolution
domain.check_integrity()

for t in domain.evolve(yieldstep = 0.1, finaltime = 4.0):
    domain.write_time()

