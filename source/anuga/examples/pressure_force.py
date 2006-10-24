"""Example of shallow water wave equation.

This example applies a dummy forcing function mimicking a
static cyclone
"""

######################
# Module imports 
#
from shallow_water import Transmissive_boundary, Reflective_boundary,\
     Dirichlet_boundary, Time_boundary
from shallow_water import Domain, Constant_height
from mesh_factory import rectangular

from math import sin, pi


#######################
# Domain
#
N = 50
elements, triangles, boundary = rectangular(N, N)
domain = Domain(elements, triangles, boundary)

domain.store = True
domain.set_name('pressure')

#######################
#Bed-slope and friction
def x_slope(x, y):
    #return -x/3
    return 0.0*x

domain.set_quantity('elevation', x_slope)
domain.set_quantity('friction', 0.1)


#######################
#Forcing terms
def pressure(x,y,t):
    r=0.1
    x0 = (1+sin(t*pi/5))/2
    y0 = (1+sin(t*pi/5))/2
    if x >= x0-r and x <= x0+r and y >= y0-r and y <= y0+r:
        return 1000*(1+sin(t*pi))
    else:
        return 800

def cyclone(domain):
    from anuga.config import rho_w 
    from anuga.pyvolution.util import gradient

    xmom = domain.quantities['xmomentum'].explicit_update
    ymom = domain.quantities['ymomentum'].explicit_update

    Stage = domain.quantities['stage']
    Elevation = domain.quantities['elevation']    
    h = Stage.vertex_values - Elevation.vertex_values
    x = domain.get_vertex_coordinates()

    for k in range(len(domain)):
        avg_h = sum( h[k,:] )/3
        
        #Compute bed slope
        x0, y0, x1, y1, x2, y2 = x[k,:]    

        p0 = pressure(x0,y0,t)
        p1 = pressure(x1,y1,t)
        p2 = pressure(x2,y2,t)        
        
        px, py = gradient(x0, y0, x1, y1, x2, y2, p0, p1, p2)

        #Update momentum
        xmom[k] += -px*avg_h/rho_w
        ymom[k] += -py*avg_h/rho_w        


domain.forcing_terms.append(cyclone)

######################
# Boundary conditions
Br = Reflective_boundary(domain)

domain.set_boundary({'left': Br, 'right': Br, 'top': Br, 'bottom': Br})


######################
#Initial condition
h = 0.05
domain.set_quantity('stage', Constant_height(x_slope, h))


######################
#Evolution
for t in domain.evolve(yieldstep = 0.05, finaltime = 7):
    domain.write_time()


    
