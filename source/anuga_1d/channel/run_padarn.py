import os
import random
from math import sqrt, pow, pi ,sin, cos
import numpy

from anuga_1d.channel.channel_domain import *

from anuga_1d.config import g, epsilon



print "Variable Width Only Test"

z_infty = 10.0       ## max equilibrium water depth at lowest point.
L_x = 2500.0         ## width of channel
A0 = 0.5*L_x         ## determines amplitudes of oscillations
omega = sqrt(2*g*z_infty)/L_x ## angular frequency of osccilation

# Define functions for initial quantities
def initial_stage(x):
    t=0.0
    y = numpy.zeros(len(x),numpy.float)
    for i in range(len(x)):
        y[i] = z_infty+2*A0*z_infty/L_x*cos(omega*t)*(x[i]/L_x-0.5*A0/(L_x)*cos(omega*t))
        #y[i] = 12.0
    return y

def bed(x):
    z = numpy.zeros(len(x),numpy.float)
    for i in range(len(x)):
        z[i] = z_infty*(x[i]**2/L_x**2)
    return z

def width(x):
    return 1.0

def initial_area(x):
    y = numpy.zeros(len(x),numpy.float)
    for i in range(len(x)):
        y[i]=(initial_stage([x[i]])-bed([x[i]]))*width(x[i])
    return y


import time



# Set final time and yield time for simulation
finaltime = 1000.0
yieldstep = 10.0

# Length of channel (m)
L = 2500.0
# Define the number of cells
number_of_cells = [100]

# Define cells for finite volume and their size
N = int(number_of_cells[0])     
print "Evaluating domain with %d cells" %N
cell_len = 4*L/N # Origin = 0.0
points = numpy.zeros(N+1,numpy.float)
for i in range(N+1):
        points[i] = -2*L +i*cell_len



# Create domain with centroid points as defined above        
domain = Domain(points)


# Set initial values of quantities - default to zero
domain.set_quantity('area', initial_area)
domain.set_quantity('elevation',bed)
domain.set_quantity('width',width)

# Set boundry type, order, timestepping method and limiter
domain.set_boundary({'exterior':Reflective_boundary(domain)})
domain.order = 2
domain.set_timestepping_method('euler')
domain.set_CFL(1.0)
domain.set_limiter("vanleer")
#domain.h0=0.0001

# Start timer
t0 = time.time()





for t in domain.evolve(yieldstep = yieldstep, finaltime = finaltime):
    domain.write_time()

N = float(N)
HeightC    = domain.quantities['height'].centroid_values
DischargeC = domain.quantities['discharge'].centroid_values
C          = domain.centroids
print 'That took %.2f seconds' %(time.time()-t0)
X          = domain.vertices
HeightQ    = domain.quantities['area'].vertex_values
VelocityQ  = domain.quantities['velocity'].vertex_values
Z          = domain.quantities['elevation'].vertex_values
Stage      = HeightQ + Z



from pylab import plot,title,xlabel,ylabel,legend,savefig,show,hold,subplot

hold(False)
   
plot1 = subplot(211)
plot(X.flat,Z.flat, X.flat,Stage.flat)
 
plot1.set_ylim([-1,35])
xlabel('Position')
ylabel('Stage')
legend(('Analytical Solution', 'Numerical Solution'),
           'upper right', shadow=True)
plot2 = subplot(212)
plot(X.flat,VelocityQ.flat)
plot2.set_ylim([-10,10])
    
xlabel('Position')
ylabel('Velocity')
   
show()
