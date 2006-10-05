"""Simple water flow example using ANUGA

Water driven up a linear slope and time varying boundary,
similar to a beach environment.

"""


#------------------------------------------------------------------------------
# Import necessary modules
#------------------------------------------------------------------------------

import sys

from anuga.pmesh.mesh_interface import create_mesh_from_regions
from abstract_2d_finite_volumes.mesh_factory import rectangular_cross
from anuga.config import g
from anuga.shallow_water import Domain
from anuga.shallow_water import Reflective_boundary
from anuga.shallow_water import Dirichlet_boundary
from anuga.shallow_water import Time_boundary
from anuga.shallow_water import Transmissive_Momentum_Set_Stage_boundary
from abstract_2d_finite_volumes.util import file_function
from pylab import plot, xlabel, ylabel, title, ion, close, savefig,\
     figure, axis, legend, grid, hold



#------------------------------------------------------------------------------
# Model constants

slope = -0.02       # 1:50 Slope, reaches h=20m 1000m from western bndry,
                    # and h=0 (coast) at 300m
highest_point = 6   # Highest elevation (m)
sea_level = 0       # Mean sea level
min_elevation = -20 # Lowest elevation (elevation of offshore flat part)
offshore_depth = sea_level-min_elevation # offshore water depth

amplitude = 0.5     # Solitary wave height H
normalized_amplitude = amplitude/offshore_depth 
simulation_name = 'runup_convergence'   
coastline_x = -highest_point/slope

# Basin dimensions (m)
west = 0          # left boundary
east = 1500       # right boundary 
south = 0         # lower boundary
north = 100       # upper boundary


#------------------------------------------------------------------------------
# Setup computational domain all units in meters
#------------------------------------------------------------------------------

# Structured mesh
dx = 30           # Resolution: Length of subdivisions on x axis (length)
dy = 30           # Resolution: Length of subdivisions on y axis (width)

length = east-west
width = north-south
points, vertices, boundary = rectangular_cross(length/dx, width/dy,
                                               len1=length, len2=width,
                                               origin = (west, south)) 

domain = Domain(points, vertices, boundary) # Create domain


# Unstructured mesh
polygon = [[east,north],[west,north],[west,south],[east,south]]
interior_polygon = [[400,north-10],[west+10,north-10],
                    [west+10,south+10],[400,south+10]]
meshname = simulation_name + '.msh'
create_mesh_from_regions(polygon,
                         boundary_tags={'top': [0], 'left': [1],
                                        'bottom': [2], 'right': [3]},
                         maximum_triangle_area=dx*dy/4,
                         filename=meshname,
                         interior_regions=[[interior_polygon,dx*dy/32]])

domain = Domain(meshname, use_cache=True, verbose = True)
domain.set_minimum_storable_height(0.01)

domain.set_name(simulation_name)


#------------------------------------------------------------------------------
# Setup initial conditions
#------------------------------------------------------------------------------

#def topography(x,y):
#    return slope*x+highest_point  # Return linear bed slope (vector)

def topography(x,y):
    """Two part topography - slope and flat part
    """

    from Numeric import zeros, Float

    z = zeros(len(x), Float) # Allocate space for return vector
    for i in range(len(x)):

        z[i] = slope*x[i]+highest_point # Linear bed slope bathymetry

        if z[i] < min_elevation:        # Limit depth
            z[i] = min_elevation

    return z       
        

        

domain.set_quantity('elevation', topography) # Use function for elevation
domain.set_quantity('friction', 0.0 )        # Constant friction 
domain.set_quantity('stage', sea_level)      # Constant initial stage


#------------------------------------------------------------------------------
# Setup boundary conditions
#------------------------------------------------------------------------------

from math import sin, pi, cosh, sqrt
Br = Reflective_boundary(domain)      # Solid reflective wall
Bd = Dirichlet_boundary([0.,0.,0.])   # Constant boundary values


def waveform(t): 
    return sea_level +\
           amplitude/cosh(((t-50)/offshore_depth)*(0.75*g*amplitude)**0.5)**2

# Time dependent boundary for stage, where momentum is set automatically
Bts = Transmissive_Momentum_Set_Stage_boundary(domain, waveform)

# Associate boundary tags with boundary objects
domain.set_boundary({'left': Br, 'right': Bts, 'top': Br, 'bottom': Br})


# Find initial runup location and height (coastline)
w0 = domain.get_maximum_inundation_elevation()
x0, y0 = domain.get_maximum_inundation_location()
print
print 'Coastline elevation = %.2f at (x,y)=(%.2f, %.2f)' %(w0, x0, y0)

# Sanity check
w_i = domain.get_quantity('stage').get_values(interpolation_points=[[x0,y0]])
print 'Interpolated elevation at (x,y)=(%.2f, %.2f) is %.2f' %(x0, y0, w_i) 


#------------------------------------------------------------------------------
# Evolve system through time
#------------------------------------------------------------------------------

w_max = w0
for t in domain.evolve(yieldstep = 1, finaltime = 300):
    domain.write_time()

    w = domain.get_maximum_inundation_elevation()
    x, y = domain.get_maximum_inundation_location()
    print '  Coastline elevation = %.2f at (x,y)=(%.2f, %.2f)' %(w, x, y)
    print    

    if w > w_max:
        w_max = w
        x_max = x
        y_max = y


print '**********************************************'
print 'Max coastline elevation = %.2f at (%.2f, %.2f)' %(w_max, x_max, y_max)
print 'Run up distance = %.2f' %sqrt( (x_max-x0)**2 + (y_max-y0)**2 )
print '**********************************************'


#-----------------------------------------------------------------------------
# Interrogate further
#---------------------------------------------------------------

# Generate time series of one "gauge" situated at right hand boundary
from anuga.abstract_2d_finite_volumes.util import sww2timeseries
production_dirs = {'.': 'test'}
swwfiles = {}
for label_id in production_dirs.keys():
    
    swwfile = simulation_name + '.sww'
    swwfiles[swwfile] = label_id
    
texname, elev_output = sww2timeseries(swwfiles,
                                      'boundary_gauge.xya',
                                      production_dirs,
                                      report = False,
                                      reportname = 'test',
                                      plot_quantity = ['stage', 'speed'],
                                      surface = False,
                                      time_min = None,
                                      time_max = None,
                                      title_on = True,
                                      verbose = True)




