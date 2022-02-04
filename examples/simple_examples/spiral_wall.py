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
center = (length/2, width/2)

# Create a domain with named boundaries "left", "right", "top" and "bottom"
domain = anuga.rectangular_cross_domain(int(length/dx), int(width/dy),
                                        len1=length, len2=width)

domain.set_name('spiral_wall')               # Output name

#------------------------------------------------------------------------------
# Setup initial conditions
#------------------------------------------------------------------------------


# Define wall polygon - spiral wall

def wall_polygon():

    N = 10
    c = center
    r_inner = 1.5
    r_outer = 2     
    
    outer_vertices = []
    inner_vertices = []    
    
    # Outer wall edge
    for i in range(N):
        theta = i * 2 * pi / N
        a = theta * 0.5
        
        x = r_outer * a * cos(theta) + c[0]
        y = r_outer * a * sin(theta) + c[1]       
        outer_vertices.append((x, y))
        
    # Inner wall edge
    for i in range(N):
        theta = i * 2 * pi / N
        a = theta * 0.5        
        
        x = r_inner * a * cos(theta) + c[0]
        y = r_inner * a * sin(theta) + c[1]       
        inner_vertices.append((x, y))        
        
    return outer_vertices + inner_vertices    
    

def topography(x, y):

    # Define wall for given polygon
    
    P = wall_polygon()
    z = x * 0.   # Flat surface
    c = (5, 5)     # Center            
    N = len(x)
    for i in range(N):

        # Raise elevation for points in polygon
        if anuga.geometry.polygon.is_inside_polygon((x[i], y[i]), P, closed=True, verbose=False):
            z[i] = 2                
            
        # Mark the center
        #if (x[i] - c[0])**2 + (y[i] - c[1])**2 < 0.02:
        #    z[i] = 4                        
        
    return z            
            

domain.set_quantity('elevation', topography) # Use function for elevation
domain.set_quantity('friction', 0.01)        # Constant friction 
domain.set_quantity('stage',                 # Dry bed
                    expression='elevation')  

                    
#------------------------------------------------------------------------------
# Setup forcing functions
#------------------------------------------------------------------------------                    
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
    
    # FIXME (OLE): Add a polygon as an alternative.
    # FIXME (OLE): Generalise to all quantities

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
                    

drain = Inflow(center=center, radius=0.0707, flow=-3.0)             
domain.forcing_terms.append(drain)
                    
#------------------------------------------------------------------------------
# Setup boundary conditions
#------------------------------------------------------------------------------
Bi = anuga.Dirichlet_boundary([0.4, 0, 0])         # Inflow
Br = anuga.Reflective_boundary(domain)             # Solid reflective walls

domain.set_boundary({'left': Br, 'right': Bi, 'top': Br, 'bottom': Br})

#------------------------------------------------------------------------------
# Evolve system through time
#------------------------------------------------------------------------------
for t in domain.evolve(yieldstep=0.2, finaltime=5):
    domain.print_timestepping_statistics()
