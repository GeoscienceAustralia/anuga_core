"""Simple water flow example using ANUGA

Water flowing along a spiral wall and draining into a hole in the centre.
"""

#------------------------------------------------------------------------------
# Import necessary modules
#------------------------------------------------------------------------------
from math import acos, cos, sin, sqrt, pi
import anuga


#------------------------------------------------------------------------------
# Setup computational domain
#------------------------------------------------------------------------------
length = 10.
width = 10.
dx = dy = 0.04           # Resolution: Length of subdivisions on both axes

# Create a domain with named boundaries "left", "right", "top" and "bottom"
domain = anuga.rectangular_cross_domain(int(length/dx), int(width/dy),
                                        len1=length, len2=width)

domain.set_name('spiral_wall')               # Output name

#------------------------------------------------------------------------------
# Setup initial conditions
#------------------------------------------------------------------------------


# Define wall polygon - circle wall

def circle_wall_polygon():

    N = 100
    c = (5, 5)     # Centre    
    r_inner = 1.5
    r_outer = 2     
    
    vertices = []
    
    # Outer wall edge
    for i in range(N):
        theta = i * 2 * pi / N
        
        x = r_outer * cos(theta) + c[0]
        y = r_outer * sin(theta) + c[1]       
        vertices.append((x, y))
        
    # Inner wall edge
    for i in range(N):
        theta = i * 2 * pi / N
        
        x = r_inner * cos(theta) + c[0]
        y = r_inner * sin(theta) + c[1]       
        vertices.append((x, y))        
        
    return vertices    
    

def topography(x, y):

    # Define wall for given polygon
    
    P = circle_wall_polygon()
    z = x * 0.   # Flat surface
    c = (5, 5)     # Centre            
    N = len(x)
    for i in range(N):

        # Raise elevation for points in polygon
        if anuga.geometry.polygon.is_inside_polygon((x[i], y[i]), P, closed=True, verbose=False):
            z[i] = 2                
            
        # Mark the centre
        if (x[i] - c[0])**2 + (y[i] - c[1])**2 < 0.02:
            z[i] = 4                        
        
    return z            
            

domain.set_quantity('elevation', topography) # Use function for elevation
domain.set_quantity('friction', 0.01)        # Constant friction 
domain.set_quantity('stage',                 # Dry bed
                    expression='elevation')  

#------------------------------------------------------------------------------
# Setup boundary conditions
#------------------------------------------------------------------------------
Bi = anuga.Dirichlet_boundary([0.4, 0, 0])         # Inflow
Br = anuga.Reflective_boundary(domain)             # Solid reflective walls

domain.set_boundary({'left': Br, 'right': Bi, 'top': Br, 'bottom': Br})

#------------------------------------------------------------------------------
# Evolve system through time
#------------------------------------------------------------------------------
for t in domain.evolve(yieldstep=0.2, finaltime=0):
    domain.print_timestepping_statistics()
