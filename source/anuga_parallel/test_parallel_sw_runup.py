#!/usr/bin/env python


"""Simple water flow example using ANUGA

Water driven up a linear slope and time varying boundary,
similar to a beach environment

This is a very simple test of the parallel algorithm using the simplified parallel API
"""


#------------------------------------------------------------------------------
# Import necessary modules
#------------------------------------------------------------------------------

from anuga.pmesh.mesh_interface import create_mesh_from_regions
from anuga.abstract_2d_finite_volumes.mesh_factory import rectangular_cross
from anuga.shallow_water import Domain
from anuga.shallow_water import Reflective_boundary
from anuga.shallow_water import Dirichlet_boundary
from anuga.shallow_water import Time_boundary
from anuga.shallow_water import Transmissive_boundary

from parallel_api import distribute, myid


#--------------------------------------------------------------------------
# Setup computational domain
#--------------------------------------------------------------------------
points, vertices, boundary = rectangular_cross(10, 10) # Basic mesh
domain = Domain(points, vertices, boundary) # Create domain
domain.set_name('runup')                    # Set sww filename
domain.set_datadir('.')                     # Set output dir


#--------------------------------------------------------------------------
# Setup initial conditions
#--------------------------------------------------------------------------

def topography(x,y): 
    return -x/2                              # linear bed slope

domain.set_quantity('elevation', topography) # Use function for elevation
domain.set_quantity('friction', 0.0)         # Constant friction 
domain.set_quantity('stage', -.5)            # Constant initial stage


#------------------------------------------------------------------------------
# Setup boundary conditions
#------------------------------------------------------------------------------

Br = Reflective_boundary(domain)      # Solid reflective wall
Bd = Dirichlet_boundary([-0.2,0.,0.]) # Constant boundary values

# Associate boundary tags with boundary objects
domain.set_boundary({'left': Br, 'right': Bd, 'top': Br, 'bottom': Br})


#--------------------------------------------------------------------------
# Create the parallel domain
#--------------------------------------------------------------------------
domain = distribute(domain, verbose=True)


#------------------------------------------------------------------------------
# Evolve system through time
#------------------------------------------------------------------------------

interpolation_points = [[0.4,0.5], [0.6,0.5], [0.8,0.5], [0.9,0.5]]
gauge_values = []
for _ in interpolation_points:
    gauge_values.append([])

time = []

for t in domain.evolve(yieldstep = 0.1, finaltime = 5.0):
    domain.write_time()

    
    # Record time series at known points
    time.append(domain.get_time())
    
    stage = domain.get_quantity('stage')
    w = stage.get_values(interpolation_points=interpolation_points)
    
    for i, _ in enumerate(interpolation_points):
        gauge_values[i].append(w[i])


print
print time
print
for i, (x,y) in enumerate(interpolation_points):
    print i, gauge_values[i]
    print 

        
    try:
        from pylab import *
    except:
        pass
    else:
        ion()
        hold(False)
        plot(time, gauge_values[i], 'r.-')
        #time, predicted_gauge_values[i], 'k-')
        
        title('Gauge %d (%f,%f)' %(i,x,y))
        xlabel('time(s)')
        ylabel('stage (m)')    
        #legend(('Observed', 'Modelled'), shadow=True, loc='upper left')
        #savefig('Gauge_%d.png' %i, dpi = 300)
    
        raw_input('Next')
        


