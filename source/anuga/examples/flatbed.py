"""Example of shallow water wave equation.

Generate slope

"""

######################
# Module imports 
#
from os import sep, path
from mesh_factory import rectangular
from shallow_water import Domain, Reflective_boundary, Dirichlet_boundary,\
     Constant_height
from Numeric import array
from anuga.pyvolution.util import Polygon_function, read_polygon


#Create basic mesh
N = 50
points, vertices, boundary = rectangular(N, N, 100, 100)

#Create shallow water domain
domain = Domain(points, vertices, boundary)
domain.store = True
domain.set_name('polygons')
print "Output being written to " + domain.get_datadir() + sep + \
	      domain.get_name() + "_size%d." %len(domain) + domain.format
	

domain.default_order=2

#Set driving forces
manning = 0.07
manning = 0.0
inflow_stage = 10.0
domain.set_quantity('friction', manning)

def wiggle(x, y):
    from Numeric import sin
    from math import pi
    return 10 + sin(2*pi*x/10)

def slope(x, y):
    return 20*(x/100+y/100)    

#Define polynomials
p0 = [[20,27], [30,25], [40,40], [20,40]]         
p1 = [[80,19], [90,20], [85,50], [80,55], [75,58], [70,60], [60,24]]         
p2 = read_polygon('testpoly.txt')


#Set elevation  
domain.set_quantity('elevation', 
	Polygon_function([(p0,slope), (p1,wiggle), (p2,15)]))
#domain.set_quantity('stage', 
#	Polygon_function([(p0,slope), (p1,wiggle), (p2,15)]))	


 	  



######################
# Boundary conditions
Br = Reflective_boundary(domain)
Bd = Dirichlet_boundary([inflow_stage,0.,0.])

domain.set_boundary({'left': Bd, 'right': Br, 'top': Br, 'bottom': Br})
domain.check_integrity()


######################
#Evolution
for t in domain.evolve(yieldstep = 1, finaltime = 1):
    domain.write_time()

print 'Done'    

