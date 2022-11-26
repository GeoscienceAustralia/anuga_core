
import matplotlib.pyplot as plt
import matplotlib.tri as tri
from pylab import *

# The abstract Python-MPI interface
from anuga.utilities.parallel_abstraction import size, rank, get_processor_name
from anuga.utilities.parallel_abstraction import finalize, send, receive
from anuga.utilities.parallel_abstraction import pypar_available, barrier

import random

import numpy as num

'''
Adds random noise to stage quantities (centroid) in the domain.
domain - domain where noise is to be added
samples - number of centroids to add noise to
dthr - minimum depth threshold, stages corresponding to depths less than this value are not modified
mag - (relative) magnitude of the noise, uniformly distributed in [-mag, mag]
'''
def addRandomNoiseToStage(domain, samples = 50, dthr = 0.1, mag = 1E-15):

    # Calculate depth
    depth = domain.get_quantity('stage').centroid_values[domain.tri_full_flag == 1] \
        - domain.get_quantity('elevation').centroid_values[domain.tri_full_flag == 1]

    stage = domain.get_quantity('stage')

    # Sample centroid values with depth greater than dthr
    n_tri = int(num.sum(domain.tri_full_flag))
    tri_index = num.array(list(range(0,n_tri)))
    submerged = tri_index[depth > dthr]
    submerged = random.sample(submerged, min(samples, len(submerged)))

    # Add random noise to sampled centroid values
    if len(submerged) > 0:
        for i in range(len(submerged)):
            new_val = (1.0 + random.uniform(-1.0, 1.0)*mag)*stage.centroid_values[submerged[i]]
            domain.get_quantity('stage').centroid_values.put(submerged[i], new_val)

'''
Plot errors in centroid quantities, triangles with errors are denoted in blue, while
triangles without errors are denoted in green.

domain - domain for which the error plots will be displayed
control_data - control data to compare quantity values against (must have serial indexing)
rthr - relative error threshold
athr - absolute error threshold
quantity - quantity e.g stage
filename - output filename
'''

def plotCentroidError(domain, control_data, rthr = 1E-7, athr = 1E-12, 
                      quantity = 'stage', filename = 'centroid_error.png'):

    n_triangles = num.sum(domain.tri_full_flag)
    
    if size() > 1:
        # If parallel, translate control data to parallel indexing
        local_control_data = num.zeros(n_triangles)
        inv_tri_map = domain.get_inv_tri_map()
        
        for i in range(n_triangles):
            local_control_data[i] = control_data[inv_tri_map[(rank(), i)]]
    else:
        local_control_data = control_data

    # Evaluate absolute and relative difference between control and actual values
    stage = domain.get_quantity(quantity)
    actual_data = stage.centroid_values[:n_triangles]
    adiff = num.fabs((actual_data - local_control_data))
    rdiff = adiff/num.fabs(local_control_data)

    # Compute masks for error (err_mask) and non-error (acc_mask) vertex indices based on thresholds
    vertices = domain.get_vertex_coordinates()
    err_mask = rdiff > rthr
    err_mask[adiff <= athr] = False    
    err_mask = num.repeat(err_mask, 3)

    acc_mask = ~err_mask
    inv_tri_map = domain.get_inv_tri_map()

    # Plot error and non-error triangle
    if rank() == 0:
        fx = {}
        fy = {}
        gx = {}
        gy = {}

        fx[0] = vertices[acc_mask,0]
        fy[0] = vertices[acc_mask,1]
        gx[0] = vertices[err_mask,0]
        gy[0] = vertices[err_mask,1]

        # Receive vertex indices of non-error triangles (fx, fy) and error triangles (gx, gy)
        for i in range(1,size()):
            fx[i] = receive(i)
            fy[i] = receive(i)
            gx[i] = receive(i)
            gy[i] = receive(i)

        # Plot non-error triangles in green
        for i in range(0,size()):
            n = int(len(fx[i])//3)
            triang = num.array(list(range(0,3*n)))
            triang.shape = (n, 3)

            if len(fx[i]) > 0:
                plt.triplot(fx[i], fy[i], triang, 'g-')

        # Plot error triangles in blue
        for i in range(0,size()):
            n = int(len(gx[i])//3)
                            
            triang = num.array(list(range(0,3*n)))
            triang.shape = (n, 3)

            if len(gx[i]) > 0: 
                plt.triplot(gx[i], gy[i], triang, 'b--')
                
        # Save plot
        plt.savefig(filename)
        
    else:
        # Send error and non-error vertex indices to Proc 0
        send(vertices[acc_mask,0], 0)
        send(vertices[acc_mask,1], 0)
        send(vertices[err_mask,0], 0)
        send(vertices[err_mask,1], 0)      
