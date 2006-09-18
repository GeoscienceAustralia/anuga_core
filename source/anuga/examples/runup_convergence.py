"""Simple water flow example using ANUGA

Water driven up a linear slope and time varying boundary,
similar to a beach environment.

The study area is discretised as a regular triangular grid 100m x 100m
"""


#------------------------------------------------------------------------------
# Import necessary modules
#------------------------------------------------------------------------------

import sys

from anuga.pmesh.mesh_interface import create_mesh_from_regions
from abstract_2d_finite_volumes.mesh_factory import rectangular_cross
from anuga.shallow_water import Domain
from anuga.shallow_water import Reflective_boundary
from anuga.shallow_water import Dirichlet_boundary
from anuga.shallow_water import Time_boundary
from anuga.shallow_water import Transmissive_Momentum_Set_Stage_boundary
from abstract_2d_finite_volumes.util import file_function
from pylab import plot, xlabel, ylabel, title, ion, close, savefig, figure, axis, legend, grid, hold



#------------------------------------------------------------------------------
# Model constants

slope = -0.02       # 1:50 Slope, reaches h=20m 1000m from western bndry, and h=0 (coast) at 300m
highest_point = 6   # Highest elevation (m)
sea_level = 0       # Mean sea level
min_elevation = -20 # Lowest elevation (elevation of offshore flat part)
offshore_depth=sea_level-min_elevation # offshore water depth
amplitude = 0.5       # Solitary wave height H
normalized_amplitude = amplitude/offshore_depth 
simulation_name = 'runup_convergence'   


# Basin dimensions (m)
west = 0          # left boundary
east = 1500        # right boundary 
south = 0         # lower boundary
north = 100       # upper bourdary


#------------------------------------------------------------------------------
# Setup computational domain all units in meters
#------------------------------------------------------------------------------

# Structured mesh
dx = 20           # Resolution: Lenght of subdivisions on x axis (length)
dy = 20           # Resolution: Lenght of subdivisions on y axis (width)

length = east-west
width = north-south
points, vertices, boundary = rectangular_cross(length/dx, width/dy, len1=length, len2=width,
                                               origin = (west, south)) 

domain = Domain(points, vertices, boundary) # Create domain



# Unstructured mesh
polygon = [[east,north],[west,north],[west,south],[east,south]]
interior_polygon = [[400,north-10],[west+10,north-10],[west+10,south+10],[400,south+10]]
meshname = simulation_name + '.msh'
create_mesh_from_regions(polygon,
                         boundary_tags={'top': [0], 'left': [1], 'bottom': [2], 'right': [3]},
                         maximum_triangle_area=dx*dy/4, # Triangle area commensurate with structured mesh
                         filename=meshname,
                         interior_regions=[[interior_polygon,dx*dy/32]])
domain = Domain(meshname, use_cache=True, verbose = True)



domain.set_name(simulation_name)


#------------------------------------------------------------------------------
# Setup initial conditions
#------------------------------------------------------------------------------

#def topography(x,y):
#    return slope*x+highest_point  # Return linear bed slope bathymetry as vector

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
    return sea_level + amplitude/cosh(((t-50)/offshore_depth)*(0.75*9.81*amplitude)**0.5)**2

Bw = Time_boundary(domain=domain,     # Time dependent boundary  
                   f=lambda t: [waveform(t), 0.0, 0.0])

# Time dependent boundary for stage, where momentum is set automatically
Bts = Transmissive_Momentum_Set_Stage_boundary(domain, waveform)


# Associate boundary tags with boundary objects
domain.set_boundary({'left': Br, 'right': Bts, 'top': Br, 'bottom': Br})


#------------------------------------------------------------------------------
# Evolve system through time
#------------------------------------------------------------------------------

stagestep = []
for t in domain.evolve(yieldstep = 1, finaltime = 300):
    domain.write_time()


#-----------------------------------------------------------------------------
# Interrogate solution
#-----------------------------------------------------------------------------


# Define line of gauges through center of domain
def gauge_line(west,east,north,south):
    from Numeric import arange
    x_vector = arange(west,600, 1) # Gauges every 1 meter from west to 600m from western bdry
    y = (north+south)/2.

    gauges = []
    for x in x_vector:
        gauges.append([x,y])
        
    return gauges, x_vector

gauges, x_vector = gauge_line(west,east,north,south)

# Obtain interpolated timeseries at gauges
f = file_function(domain.get_name()+'.sww',
                  quantities = ['stage', 'elevation', 'xmomentum', 'ymomentum'],
                  interpolation_points = gauges,
                  verbose = True,
                  use_cache = True)


# Find runup distance from western boundary through a linear search
max_stage = []
min_stage = []
runup_point = west
coastline = east        
for k, g in enumerate(gauges):
    z = f(0, point_id=k)[1] # Elevation

    min_w = sys.maxint
    max_w = -min_w
    for i, t in enumerate(f.get_time()):
        w = f(t, point_id = k)[0]
        if w > max_w: max_w = w
        if w < min_w: min_w = w        

    if max_w-z <= 0.01:  # Find first gauge where max depth > eps (runup)
        runup_point = g[0]

    if min_w-z <= 0.01:  # Find first gauge where min depth > eps (coastline)
        coastline = g[0]        
        
    max_stage.append(max_w)
    min_stage.append(min_w)    


# Print
print 'wave height [m]:                    ', amplitude
runup_height = topography([runup_point], [(north+south)/2.])[0]
print 'run up height [m]:                  ', runup_height 

runup_distance = runup_point-coastline
print 'run up distance from coastline [m]: ', runup_distance

print 'Coastline (meters form west):       ', coastline



# Take snapshots and plot
ion()
figure(1)
plot(x_vector, topography(x_vector,(north+south)/2.), 'r-')
xlabel('x')
ylabel('Elevation')
#legend(('Max stage', 'Min stage', 'Elevation'), shadow=True, loc='upper right')
title('Stage snapshots (t=0, 10, ...) for gauge line')
grid()
hold(True)

for i, t in enumerate(f.get_time()):
    if i % 10 == 0:
        # Take only some timesteps to avoid clutter
        stages = []    
        for k, g in enumerate(gauges):
            w = f(t, point_id = k)[0]        
            stages.append(w)

        plot(x_vector, stages, 'b-')
         
savefig('snapshots')



# Store
filename = 'maxrunup'+str(amplitude)+'.csv'
fid = open(filename,'w')    
s = 'Waveheight,Runup distance,Runup height\n'
fid.write(s)

s = '%.2f,%.2f,%.2f\n' %(amplitude, runup_distance, runup_height)
fid.write(s)

fid.close()

# Plot max runup etc
ion()
figure(1)
plot(x_vector, max_stage, 'g+',
     x_vector, min_stage, 'b+',     
     x_vector, topography(x_vector,(north+south)/2.), 'r-')
xlabel('x')
ylabel('stage')
legend(('Max stage', 'Min stage', 'Elevation'), shadow=True, loc='upper right')
title('Maximum stage for gauge line')
grid()
#axis([33000, 47000, -1000, 3000])
savefig('max_stage')

close('all')
    


