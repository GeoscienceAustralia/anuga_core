"""Simple water flow example using ANUGA

Water flowing down a channel with a topography that varies with time
"""

#------------------------------------------------------------------------------
# Import necessary modules
#------------------------------------------------------------------------------
from anuga import rectangular_cross
from anuga import Domain
from anuga import Reflective_boundary
from anuga import Dirichlet_boundary
from anuga import Time_boundary

#------------------------------------------------------------------------------
# Setup computational domain
#------------------------------------------------------------------------------
length = 24.
width = 5.
dx = dy = 0.1           # Resolution: Length of subdivisions on both axes

points, vertices, boundary = rectangular_cross(int(length/dx), int(width/dy),
                                               len1=length, len2=width)
domain = Domain(points, vertices, boundary)
domain.set_name('channel_variable_bed_dx=%.2f_dy=%.2f' % (dx, dy)) # Output name

print(domain.statistics())
domain.set_quantities_to_be_stored({'elevation': 2,
                                    'stage': 2})

#------------------------------------------------------------------------------
# Setup initial conditions
#------------------------------------------------------------------------------
def topography(x,y):
    """Complex topography defined by a function of vectors x and y."""

    z = -x/100
    
    N = len(x)
    for i in range(N):
        # Step
        if 2 < x[i] < 4:
            z[i] += 0.4 - 0.05*y[i]
    
        # Permanent pole
        if (x[i] - 8)**2 + (y[i] - 2)**2 < 0.4**2:
            z[i] += 1
    return z
    

def pole_increment(x,y):
    """This provides a small increment to a pole located mid stream
    For use with variable elevation data
    """
    
    z = 0.0*x

    N = len(x)
    for i in range(N):
        # Pole 1
        if (x[i] - 12)**2 + (y[i] - 3)**2 < 0.4**2:
            z[i] += 0.01
            
    for i in range(N):
        # Pole 2
        if (x[i] - 14)**2 + (y[i] - 2)**2 < 0.4**2:
            z[i] += 0.005
                        
    return z
    
    
domain.set_quantity('elevation', topography)           # elevation is a function
domain.set_quantity('friction', 0.01)                  # Constant friction
domain.set_quantity('stage', expression='elevation')   # Dry initial condition

#------------------------------------------------------------------------------
# Setup boundary conditions
#------------------------------------------------------------------------------
Bi = Dirichlet_boundary([0.4, 0, 0])          # Inflow
Br = Reflective_boundary(domain)              # Solid reflective wall
Bo = Dirichlet_boundary([-5, 0, 0])           # Outflow

domain.set_boundary({'left': Bi, 'right': Bo, 'top': Br, 'bottom': Br})

#------------------------------------------------------------------------------
# Evolve system through time
#------------------------------------------------------------------------------

growing = False
shrinking = False
done = False
for t in domain.evolve(yieldstep=0.1, finaltime=40.0):
    domain.print_timestepping_statistics()

    #w = domain.get_quantity('stage').\
    #    get_values(interpolation_points=[[18, 2.5]])
    #print 'Level at gauge point = %.2fm' % w
           
    #z = domain.get_quantity('elevation').\
    #    get_values(interpolation_points=[[12, 3]])           
    #print 'Elevation at pole location = %.2fm' % z           

    # Start variable elevation after 10 seconds    
    if t > 10 and not (shrinking or growing or done):
        growing = True
           
    # Stop growing when pole has reached a certain height
    if t > 16 and growing:
        growing = False
        shrinking = False
        
    # Start shrinking
    if t > 20:
        shrinking = True
        growing = False
        
    # Stop changing when pole has shrunk to original level
    if t > 25 and shrinking:
        done = True
        shrinking = growing = False
        domain.set_quantity('elevation', topography)

    # Grow or shrink               
    if growing:       
        domain.add_quantity('elevation', pole_increment)
        
    if shrinking:    
        domain.add_quantity('elevation', lambda x,y: -2*pole_increment(x,y))    
        
