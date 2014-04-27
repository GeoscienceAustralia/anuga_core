import os
import random
from math import sqrt, pow, pi 
from channel_domain import *
from numpy import allclose, array, zeros, ones, take, sqrt
from anuga_1d.config import g, epsilon


print "Variable Width Only Test"

# Define functions for initial quantities



#def initial_area(x):
#    y = zeros(len(x),'f')
#    for i in range(len(x)):
#        y[i]=(10-randomarray[i])
#    return y



def bed(x):

    y = zeros(len(x),'f')
    for i in range(len(x)):
        y[i]=randomarray[i]
        randomarray[i]=random.normalvariate(3,1)
    return y


  


def width(x):
    y = zeros(len(x),'f')
    return y+1


stage = 6.0

def initial_area(x):

    a_bed = bed(x)
    a_width = width(x)

    a_height = 6.0 - a_bed

    y = a_height*a_width

    return y


import time



# Length of channel (m)
L = 1000.0   
# Define the number of cells
number_of_cells = [50]

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


# Set initial values of quantities - default to zero
#domain.set_quantity('stage',6.0)
domain.set_quantity('elevation',bed, location='centroids')
domain.set_quantity('width',width, location='centroids')
domain.set_quantity('area', initial_area, location='centroids')


#domain.setstageflag = True
# Set boundry type, order, timestepping method and limiter
domain.set_boundary({'exterior':Reflective_boundary(domain)})
domain.order = 2
domain.set_timestepping_method('rk2')
domain.set_CFL(1.0)
domain.set_limiter("vanleer")
#domain.h0=0.0001


#domain.distribute_to_vertices_and_edges()


AreaC = domain.quantities['area'].centroid_values
BedC = domain.quantities['elevation'].centroid_values
WidthC = domain.quantities['width'].centroid_values
#
AreaC[:] = (8.0 - BedC)* WidthC

#domain.set_quantity('area', initial_area)

#domain.distribute_to_vertices_and_edges()

# Start timer
t0 = time.time()
i=0

print 'elevation vertex values'
print domain.quantities['elevation'].vertex_values
print 'stage vertex values'
print domain.quantities['stage'].vertex_values
print 'area vertex values'
print domain.quantities['area'].vertex_values
print 'width vertex values'
print domain.quantities['width'].vertex_values


domain.distribute_to_vertices_and_edges()



# Set final time and yield time for simulation
finaltime = 100.0
yieldstep = 1.0

domain.initialize_plotting(stage_lim = [-2.0, 12.0],
                           width_lim = [0.0, 2.0],
                           velocity_lim = [-10.0, 10.0])

for t in domain.evolve(yieldstep = yieldstep, finaltime = finaltime):
    domain.write_time()

    domain.update_plotting()


print domain.quantities['elevation'].vertex_values

domain.hold_plotting(save="not_well_balanced_random_depth")


