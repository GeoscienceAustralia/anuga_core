import os
from math import sqrt, pow, pi 

import numpy as np

from anuga_1d.channel.channel_domain import *
from anuga_1d.config import g, epsilon
from anuga_1d.base.generic_mesh import uniform_mesh


print "Channel Flow 1 Padarn Test"

# Define functions for initial quantities
def initial_area(x):
    return 1.4691*width(x)

def width(x):
    x1=(x/1000)*(x/1000)
    x2=x1*(x/1000)
    x3=x2*(x/1000)
    return 10-64*(x1-2*x2+x3)

def bed(x):
    y = np.zeros(len(x),np.float)
    for i in range(len(x)):
        if x[i]<525 and x[i]>475:
            y[i]=0
        else:
            y[i]=0
    return y
 

def initial_discharge(x):
    return 20
 
import time



# Length of channel (m)
L = 1000.0   
# Define the number of cells
number_of_cells = [200]

# Define cells for finite volume and their size
N = int(number_of_cells[0])     
print "Evaluating domain with %d cells" %N
cell_len = L/N # Origin = 0.0
points = np.zeros(N+1,np.float)

# Define the centroid points
for j in range(N+1):
    points[j] = j*cell_len

# Create domain with centroid points as defined above        
domain = Domain(points)

# Set initial values of quantities - default to zero
domain.set_quantity('area', initial_area)
domain.set_quantity('width',width)
domain.set_quantity('elevation',bed)
domain.set_quantity('discharge',initial_discharge)

# Set boundry type, order, timestepping method and limiter
#domain.set_boundary({'exterior': Dirichlet_boundary([14,20,0,1.4,20/14,9])})
domain.set_boundary({'exterior': Reflective_boundary(domain)})
domain.order = 2
domain.set_timestepping_method('rk2')
domain.set_CFL(1.0)
domain.set_limiter("vanleer")
#domain.h0=0.0001

# Start timer
t0 = time.time()

# Set final time and yield time for simulation
finaltime = 10.0
yieldstep = finaltime

#===================================================================
# Time loop
#===================================================================
for t in domain.evolve(yieldstep = yieldstep, finaltime = finaltime):
    domain.write_time()

N = float(N)
HeightC = domain.quantities['height'].centroid_values
DischargeC = domain.quantities['discharge'].centroid_values
C = domain.centroids
print 'That took %.2f seconds' %(time.time()-t0)
X = domain.vertices
HeightQ = domain.quantities['height'].vertex_values
VelocityQ = domain.quantities['velocity'].vertex_values
x = X.flat
z = domain.quantities['elevation'].vertex_values
stage=HeightQ+z
stage = stage.flat
z = z.flat

from pylab import plot,title,xlabel,ylabel,legend,savefig,show,hold,subplot

hold(False)
   
plot1 = subplot(211)

plot(x,z,x,stage)
 
plot1.set_ylim([-1,11])
xlabel('Position')
ylabel('Stage')
legend(('Analytical Solution', 'Numerical Solution'),
           'upper right', shadow=True)
plot2 = subplot(212)
plot(x,VelocityQ.flat)
plot2.set_ylim([-10,10])
    
xlabel('Position')
ylabel('Velocity')
   
show()
    
