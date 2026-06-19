"""  Test environmental forcing - rain, wind, etc.
"""

import unittest, os

import anuga

from anuga import Reflective_boundary
from anuga import rectangular_cross_domain

from anuga import Domain

import numpy as num
import warnings
import time
import math

from anuga.shallow_water.sw_domain_cuda import nvtxRangePush, nvtxRangePop


nx = 50
ny = 50

def create_domain(name='domain'):

    domain = anuga.rectangular_cross_domain(nx, ny, len1=1., len2=1.)

    domain.set_flow_algorithm('DE0')
    domain.set_low_froude(0)

    domain.set_name(name)  
    domain.set_datadir('.')

    #------------------
    # Define topography
    #------------------
    scale_me=1.0

    #def topography(x,y):
    #    return (-x/2.0 +0.05*num.sin((x+y)*50.0))*scale_me

    def topography(x,y):
        return 0.0

    #def stagefun(x,y):
    #    stage=-0.2*scale_me #+0.01*(x>0.9)
    #    return stage

    def stagefun(x,y):
        stage=1.0-0.5*x
        return stage

    domain.set_quantity('elevation',topography)     # Use function for elevation
    domain.set_quantity('friction',0.03)            # Constant friction
    domain.set_quantity('stage', stagefun)          # Constant negative initial stage

    #--------------------------
    # Setup boundary conditions
    #--------------------------
    Br=anuga.Reflective_boundary(domain)                 # Solid reflective wall
    Bd=anuga.Dirichlet_boundary([-0.1*scale_me,0.,0.])   # Constant boundary values -- not used in this example

    #----------------------------------------------
    # Associate boundary tags with boundary objects
    #----------------------------------------------
    domain.set_boundary({'left': Br, 'right': Bd, 'top': Br, 'bottom':Br})

    return domain

print('')
print(70*'=')
print('Test Runup')
print(70*'=')


nvtxRangePush('create domain1')
domain1 = create_domain('domain_original')
domain1.set_multiprocessor_mode(1)
nvtxRangePop()

nvtxRangePush('create domain1')
domain2 = create_domain('domain_cuda')
domain2.set_multiprocessor_mode(1) # will change to 2 once burn in
nvtxRangePop()

#------------------------------
#Evolve the system through time
#------------------------------
yieldstep = 0.0002
finaltime = 0.0002
nvtxRangePush('evolve domain1')
print('Evolve domain1')
print('domain1 number of triangles ',domain1.number_of_elements)
for t in domain1.evolve(yieldstep=yieldstep,finaltime=finaltime):
    domain1.print_timestepping_statistics()
nvtxRangePop()

nvtxRangePush('evolve domain2')
print('Evolve domain2')
print('domain2 number of triangles ',domain2.number_of_elements)
for t in domain2.evolve(yieldstep=yieldstep,finaltime=finaltime):
    domain2.print_timestepping_statistics()
nvtxRangePop()



#---------------------------------------
# run domain1 using standard routine
#---------------------------------------
timestep = 0.1

#nvtx marker
nvtxRangePush('distribute_to_vertices_and_edges')
# From centroid values calculate edge and vertex values
domain1.distribute_to_vertices_and_edges()
#nvtx marker
nvtxRangePop()
#nvtx marker
nvtxRangePush('update_boundary')
# Apply boundary conditions
domain1.update_boundary()
#nvtx marker
nvtxRangePop()
#nvtx marker
nvtxRangePush('compute_fluxes')
# Compute fluxes across each element edge
domain1.compute_fluxes()
#nvtx marker
nvtxRangePop()
#nvtx marker
nvtxRangePush('compute_forcing_terms')
# Compute forcing terms
domain1.compute_forcing_terms()
#nvtx marker
nvtxRangePop()



##########################################################
# GPU INTERFACE
##########################################################

#nvtx marker
nvtxRangePush('distribute_to_vertices_and_edges')
# From centroid values calculate edge and vertex values
domain2.distribute_to_vertices_and_edges()
#nvtx marker
nvtxRangePop()
#nvtx marker
nvtxRangePush('update_boundary')
# Apply boundary conditions
domain2.update_boundary()
#nvtx marker
nvtxRangePop()
#nvtx marker
nvtxRangePush('compute_fluxes')
# Compute fluxes across each element edge
domain2.compute_fluxes()
#nvtx marker
nvtxRangePop()

from anuga.shallow_water.sw_domain_cuda import GPU_interface
gpu_interface2 = GPU_interface(domain2)

# Some of these arrays are "static" and so only
# need to be allocated and set once per simulation
nvtxRangePush('allocate gpu arrays for domain2')
gpu_interface2.allocate_gpu_arrays()
nvtxRangePop()



nvtxRangePush('compile gpu kernels for domain2')
gpu_interface2.compile_gpu_kernels()
nvtxRangePop()

nvtxRangePush('compute forcing terms on gpu for domain2')
from anuga.shallow_water.shallow_water_domain import manning_friction_implicit
domain2.set_multiprocessor_mode(2)
manning_friction_implicit(domain2)
nvtxRangePop()


quantities1 = domain1.quantities
stage1 = quantities1["stage"]
xmom1 = quantities1["xmomentum"]
ymom1 = quantities1["ymomentum"]
max_speed_1 = domain1.max_speed


quantities2 = domain2.quantities
stage2 = quantities2["stage"]
xmom2 = quantities2["xmomentum"]
ymom2 = quantities2["ymomentum"]
max_speed_2 = domain2.max_speed

N = domain1.number_of_elements
import math
sqrtN = 1.0/N

print('xmom semi implicit update diff L2-norm    ', num.linalg.norm(xmom1.semi_implicit_update-xmom2.semi_implicit_update)*sqrtN)
print('ymom semi implicit update diff L2-norm    ', num.linalg.norm(ymom1.semi_implicit_update-ymom2.semi_implicit_update)*sqrtN)


