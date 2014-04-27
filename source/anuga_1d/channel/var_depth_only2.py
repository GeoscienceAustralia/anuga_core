import os
import random
from math import sqrt, pow, pi 
from channel_domain_Ab import *
from Numeric import allclose, array, zeros, ones, Float, take, sqrt
from config import g, epsilon


print "Variable Width Only Test"

# Define functions for initial quantities
def initial_area(x):
    y = zeros(len(x),Float)
    for i in range(len(x)):
        #if x[i]>=500:
        #    y[i]=0
        #else:
       # if x[i]<900:
        y[i]=(5-bed([x[i]]))*bed([x[i]])
       # else:
         #   y[i]=10-bed([x[i]])
    return y

def bed(x):
    y = zeros(len(x),Float)
    for i in range(len(x)):
        if x[i]>100:
            y[i]=2
        else:
            y[i]=1
    return y

def width(x):
    return 1


 
import time

# Set final time and yield time for simulation
finaltime =100.0
yieldstep = finaltime

# Length of channel (m)
L = 1000.0   
# Define the number of cells
number_of_cells = [10]

# Define cells for finite volume and their size
N = int(number_of_cells[0])     
print "Evaluating domain with %d cells" %N
cell_len = L/N # Origin = 0.0
points = zeros(N+1,Float)

# Define the centroid points
for j in range(N+1):
    points[j] = j*cell_len

# Create domain with centroid points as defined above        
domain = Domain(points)


# Set initial values of quantities - default to zero
domain.set_quantity('area', initial_area)
domain.set_quantity('elevation',bed)
domain.set_quantity('width',bed)

# Set boundry type, order, timestepping method and limiter
domain.set_boundary({'exterior':Reflective_boundary(domain)})
domain.order = 2
domain.set_timestepping_method('rk2')
domain.set_CFL(1.0)
domain.set_limiter("vanleer")
#domain.h0=0.0001

# Start timer
t0 = time.time()

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
z = domain.quantities['elevation'].vertex_values.flat
stage=HeightQ.flat+z
b = domain.quantities['width'].vertex_values.flat


from pylab import plot,title,xlabel,ylabel,legend,savefig,show,hold,subplot

hold(False)
   
plot1 = subplot(211)

plot(x,z,x,stage,x,b)
 
plot1.set_ylim([-1,11])
xlabel('Position')
ylabel('Stage')
## legend(('Analytical Solution', 'Numerical Solution'),
##            'upper right', shadow=True)
plot2 = subplot(212)
plot(x,VelocityQ.flat)
plot2.set_ylim([-10,10])
    
xlabel('Position')
ylabel('Velocity')
   
show()
