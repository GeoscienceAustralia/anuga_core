"""Example of shallow water wave equation.

Generate slope

"""

######################
# Module imports 
#
from mesh_factory import rectangular
from shallow_water import Domain, Reflective_boundary, Dirichlet_boundary,\
     Constant_height, Time_boundary, File_boundary
from Numeric import array

#Create basic mesh
points, vertices, boundary = rectangular(10, 10, 100, 100)

#Create shallow water domain
domain = Domain(points, vertices, boundary)
domain.smooth = False
domain.visualise = True
domain.default_order=2

#######################
#Bed-slope and friction
def x_slope(x, y):
    return -x/3

#domain.set_quantity('elevation', x_slope)
domain.set_quantity('elevation', lambda x,y: -x/3)
domain.set_quantity('friction', 0.1)


######################
# Boundary conditions


#Write file
import os, time
from anuga.config import time_format
from math import sin, pi

finaltime = 100
filename = 'bed_w_boundary'
fid = open(filename + '.txt', 'w')
start = time.mktime(time.strptime('2000', '%Y'))
dt = 5  #Five second intervals
t = 0.0
while t <= finaltime:
    t_string = time.strftime(time_format, time.gmtime(t+start))    
    fid.write('%s, %f %f %f\n' %(t_string, 10*sin(t*0.1*pi), 0.0, 0.0))
    
    t += dt
    
fid.close()


#Convert ASCII file to NetCDF (Which is what we really like!)
from anuga.pyvolution.data_manager import timefile2swww        
timefile2swww(filename, quantity_names = domain.conserved_quantities)


Br = Reflective_boundary(domain)
Bd = Dirichlet_boundary([0.2,0.,0.])
Bw = Time_boundary(domain=domain,
                   f=lambda t: [(10*sin(t*0.1*pi)), 0.0, 0.0])

Bf = File_boundary(filename + '.sww', domain)

#domain.set_boundary({'left': Bw, 'right': Br, 'top': Br, 'bottom': Br})
domain.set_boundary({'left': Bf, 'right': Br, 'top': Br, 'bottom': Br})


######################
#Initial condition
h = 0.5
h = 0.0
domain.set_quantity('stage', Constant_height(x_slope, h))

domain.check_integrity()

  
######################
#Evolution
for t in domain.evolve(yieldstep = 1, finaltime = 100.0):
    domain.write_time()

