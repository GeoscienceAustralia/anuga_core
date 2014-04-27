import os
import random
from math import sqrt, pow, pi ,sin, cos
from channel_domain import *
from Numeric import allclose, array, zeros, ones, Float, take, sqrt
from config import g, epsilon


print "Variable Width Only Test"

# Define functions for initial quantities
def initial_area(x):   
    z_infty = 10.0       ## max equilibrium water depth at lowest point.
    L_x = 2500.0         ## width of channel
    A0 = 0.5*L_x                  ## determines amplitudes of oscillations
    omega = sqrt(2*g*z_infty)/L_x ## angular frequency of osccilation
    t=0.0
    y = zeros(len(x),Float)
    for i in range(len(x)):
        y[i] = max(z_infty+2*A0*z_infty/L_x*cos(omega*t)*(x[i]/L_x-0.5*A0/(L_x)*cos(omega*t))-z_infty*(x[i]**2/L_x**2),0.0)
    return y

def width(x):
    y = zeros(len(x),Float)
    for i in range(len(x)):
        y[i]=1
    return y

def bed(x):
    N = len(x)
    z_infty = 10.0
    z = zeros(N,Float)
    L_x = 2500.0
    A0 = 0.5*L_x
    omega = sqrt(2*g*z_infty)/L_x
    for i in range(N):
        z[i] = z_infty*(x[i]**2/L_x**2)
    return z

import time

# Set final time and yield time for simulation
finaltime = 10.0
yieldstep = finaltime

# Length of channel (m)
L = 2500.0   
# Define the number of cells
number_of_cells = [200]

# Define cells for finite volume and their size
N = int(number_of_cells[0])     
print "Evaluating domain with %d cells" %N
cell_len = L/N*4 # Origin = 0.0
points = zeros(N+1,Float)

# Define the centroid points
points = zeros(N+1,Float)
for i in range(N+1):
        points[i] = -2*L +i*cell_len

# Create domain with centroid points as defined above        
domain = Domain(points)

# Define random array for width
randomarray=zeros(len(points),Float)
for j in range(N+1):
    randomarray[j]=random.normalvariate(5,1)


# Set initial values of quantities - default to zero
domain.set_quantity('area', initial_area)
domain.set_quantity('width',width)
domain.set_quantity('elevation',bed)

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
b = domain.quantities['elevation'].vertex_values.flat
stage=HeightQ.flat+b

from pylab import plot,title,xlabel,ylabel,legend,savefig,show,hold,subplot

hold(False)
   
plot1 = subplot(211)

plot(x,b,x,stage)
 
plot1.set_ylim([-1,20])
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
