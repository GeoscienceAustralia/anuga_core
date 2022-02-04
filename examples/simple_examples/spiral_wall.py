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
length = 10.
width = 10.
dx = dy = 0.02           # Resolution: Length of subdivisions on both axes
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
    
    # Diagnostic plotting only
    xos = [x[0] for x in outer_vertices]
    yos = [x[1] for x in outer_vertices]    
    
    xis = [x[0] for x in inner_vertices]
    yis = [x[1] for x in inner_vertices]        
    plt.plot(xos, yos, 'bo', xis, yis, 'g*')        
    #plt.show()

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
# Setup forcing functions
#------------------------------------------------------------------------------                    
# FIXME: Let's use the built in Inflow class from ANUGA
class Inflow:
    """Class Inflow - general 'rain and drain' forcing term.
    
    Useful for implementing flows in and out of the domain.
    
    Inflow(center, radius, flow)
    
    center [m]: Coordinates at center of flow point
    radius [m]: Size of circular area
    flow [m/s]: Rate of change of quantity over the specified area.  
                This parameter can be either a constant or a function of time. 
                Positive values indicate inflow, 
                negative values indicate outflow.
    
    Examples
    
    Inflow((0.7, 0.4), 0.07, -0.2) # Constant drain at 0.2 m/s.
                                   # This corresponds to a flow of 
                                   # 0.07**2*pi*0.2 = 0.00314 m^3/s
                                   
    Inflow((0.5, 0.5), 0.001, lambda t: min(4*t, 5)) # Tap turning up to 
                                                     # a maximum inflow of
                                                     # 5 m/s over the 
                                                     # specified area 
    """
    

    def __init__(self, 
                 center=None, radius=None,
                 flow=0.0,
                 quantity_name = 'stage'):
                 
        if center is not None and radius is not None:
            assert len(center) == 2
        else:
            msg = 'Both center and radius must be specified'
            raise Exception(msg)
    
        self.center = center
        self.radius = radius
        self.flow = flow
        self.quantity = domain.quantities[quantity_name].explicit_update
        
    
    def __call__(self, domain):
    
        # Determine indices in flow area
        if not hasattr(self, 'indices'):
            center = self.center
            radius = self.radius
            
            N = len(domain)    
            self.indices = []
            coordinates = domain.get_centroid_coordinates()     
            for k in range(N):
                x, y = coordinates[k,:] # Centroid
                if ((x - center[0])**2 + (y - center[1])**2) < radius**2:
                    self.indices.append(k)    

        # Update inflow
        if callable(self.flow):
            flow = self.flow(domain.get_time())
        else:
            flow = self.flow
                   
        for k in self.indices:
            self.quantity[k] += flow                    
                    

drain = Inflow(center=center, radius=0.2, flow=0.0)  # Zero initially             
domain.forcing_terms.append(drain)

source = Inflow(center=(9.4, 6.0), radius=0.2, flow=1.0)             
domain.forcing_terms.append(source)
                    
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

    if domain.get_time() >= 14 and drain.flow == 0.0:
        print('Turning drain on')
        drain.flow = -2.5        
