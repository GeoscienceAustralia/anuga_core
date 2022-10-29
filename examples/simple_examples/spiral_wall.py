"""Simple water flow example using ANUGA

Water flowing along a spiral wall and draining into a hole in the centre.
"""

#------------------------------------------------------------------------------
# Import necessary modules
#------------------------------------------------------------------------------
from math import acos, cos, sin, sqrt, pi
import anuga
import matplotlib.pyplot as plt
import numpy as np

#------------------------------------------------------------------------------
# Setup computational domain
#------------------------------------------------------------------------------
plotting = False
length = 10.
width = 10.
dx = dy = 0.05           # Resolution: Length of subdivisions on both axes
center = (length/2 * 0.7, width/2)

# Create a domain with named boundaries "left", "right", "top" and "bottom"
domain = anuga.rectangular_cross_domain(int(length/dx), int(width/dy),
                                        len1=length, len2=width)

domain.set_name('spiral_wall')               # Output name

#------------------------------------------------------------------------------
# Setup initial conditions
#------------------------------------------------------------------------------


# Define wall polygon - spiral wall

def wall_polygon():

    N = 50
    c = center
    
    r_outer = 2
    r_inner = 1.8    
    width = 0.2
    
    outer_vertices = []
    inner_vertices = []    
    
    # Outer wall edge
    for i in range(1, N):
        theta = i * (2+0.3) * pi / N
        a = theta * 0.5  # Spiral expansion term
        
        x = r_outer * a * cos(theta) + c[0]
        y = r_outer * a * sin(theta) + c[1]       
        outer_vertices.append((x, y))
        
        vector = (x - c[0], y - c[1])
        distance = sqrt(vector[0]**2 + vector[1]**2)
        if distance > 0 and i > 6:
            x = (distance - width) * vector[0]/distance + c[0]
            y = (distance - width) * vector[1]/distance + c[1]
            
            inner_vertices.append((x, y))        
    
    if plotting:
        # Diagnostic plotting only
        xos = [x[0] for x in outer_vertices]
        yos = [x[1] for x in outer_vertices]    
        
        xis = [x[0] for x in inner_vertices]
        yis = [x[1] for x in inner_vertices]        
        plt.plot(xos, yos, 'bo', xis, yis, 'g*')        
        plt.show()

    return outer_vertices + inner_vertices[::-1]  # Reverse inner points to make polygon sensible
    

def topography(x, y):

    # Define wall for given polygon
    
    P = wall_polygon()
    z = y * 0.0    # Flat surface # Sloping surface in the y direction
    c = center     # Center            
    N = len(x)
    
    # Identify points inside polygon
    x = x.reshape(-1, 1)    
    y = y.reshape(-1, 1)
    points = np.concatenate((x, y), axis=1)    
    indices = anuga.geometry.polygon.inside_polygon(points, P, closed=True, verbose=False)

    # Raise elevation for points in polygon        
    for i in indices:
        z[i] += 1.0                
            
    return z            
            

domain.set_quantity('elevation', topography) # Use function for elevation
domain.set_quantity('friction', 0.01)        # Constant friction 
domain.set_quantity('stage',                 # Dry bed
                    expression='elevation')  

                    
#------------------------------------------------------------------------------
# Setup inflow and outflow operators
#------------------------------------------------------------------------------                    
drain_region = anuga.Region(domain, center=center, radius=0.2)
drain = anuga.Inlet_operator(domain, region=drain_region, Q=0.0)

source_region = anuga.Region(domain, center=(9.4, 6.0), radius=0.2)
source = anuga.Inlet_operator(domain, region=source_region, Q=0.1)
                    
#------------------------------------------------------------------------------
# Setup boundary conditions
#------------------------------------------------------------------------------
#Bi = anuga.Dirichlet_boundary([0.4, 0, 0])         # Inflow
Br = anuga.Reflective_boundary(domain)             # Solid reflective walls

domain.set_boundary({'left': Br, 'right': Br, 'top': Br, 'bottom': Br})

#------------------------------------------------------------------------------
# Evolve system through time
#------------------------------------------------------------------------------
for t in domain.evolve(yieldstep=0.2, finaltime=40):
    domain.print_timestepping_statistics()

    if domain.get_time() >= 14 and drain.get_Q() == 0.0:
        print('    Turning drain on')
        drain.set_Q(-2.5)       
