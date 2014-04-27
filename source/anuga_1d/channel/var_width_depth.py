import os
import random
from math import sqrt, pow, pi
from channel_domain import *
from numpy import allclose, array, zeros, ones, take, sqrt
from anuga_1d.config import g, epsilon


random.seed(111)

print "Variable Width Only Test"

# Define functions for initial quantities
def initial_area(x):
    y = zeros(len(x),'f')
    for i in range(len(x)):
        y[i]=(10-randomarray[i])

   
    return y

def width(x):

    y = zeros(len(x),'f')
    for i in range(len(x)):
        y[i]=randomarray[i]
        randomarray[i]=random.normalvariate(3,1)
    return y

def bed(x):

    y = zeros(len(x),'f')
    for i in range(len(x)):
        y[i]=randomarray2[i]
        randomarray2[i]=random.normalvariate(3,1)
    return y



def stage(x):

    y = zeros(len(x),'f')
    for i in range(len(x)):
        if x[i]<200 and x[i]>150:
            y[i]=8.0
        else:
            y[i]=8.0
    return y
  


 
import time

# Set final time and yield time for simulation
finaltime = 10.0
yieldstep = finaltime

# Length of channel (m)
L = 1000.0   
# Define the number of cells
number_of_cells = [100]

# Define cells for finite volume and their size
N = int(number_of_cells[0])     
print "Evaluating domain with %d cells" %N
cell_len = L/N # Origin = 0.0
points = zeros(N+1,'f')

# Define the centroid points
for j in range(N+1):
    points[j] = j*cell_len

# Create domain with centroid points as defined above        
domain = Domain(points)

# Define random array for width
randomarray=zeros(len(points),'f')
for j in range(N+1):
    randomarray[j]=random.normalvariate(3,1)
randomarray2=zeros(len(points),'f')
for j in range(N+1):
    randomarray2[j]=random.normalvariate(3,1)


# Set initial values of quantities - default to zero
domain.set_quantity('stage',stage,'centroids')
domain.set_quantity('elevation', bed)
domain.set_quantity('width',width)
domain.setstageflag = True
# Set boundry type, order, timestepping method and limiter
domain.set_boundary({'exterior':Reflective_boundary(domain)})
domain.order = 2
domain.set_timestepping_method('rk2')
domain.set_CFL(1.0)
domain.set_limiter("vanleer")
#domain.h0=0.0001

# Start timer
t0 = time.time()


AreaC = domain.quantities['area'].centroid_values
BedC = domain.quantities['elevation'].centroid_values
WidthC = domain.quantities['width'].centroid_values
#
AreaC[:] = (10.0 - BedC)* WidthC


#print domain.quantities['elevation'].vertex_values

finaltime = 100.0
yieldstep = 1.0

domain.initialize_plotting(stage_lim = [-1.0, 20.0],
                           width_lim = [0.0, 6.0],
                           velocity_lim = [-15.0, 15.0])


for t in domain.evolve(yieldstep = yieldstep, finaltime = finaltime):
    domain.write_time()

    domain.update_plotting()

#print domain.quantities['elevation'].vertex_values

domain.hold_plotting(save='not_well_balanced_random_depth_and_width')

