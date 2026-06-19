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

nvtxRangePush('distribute domain1')
domain1.distribute_to_vertices_and_edges()
nvtxRangePop()

nvtxRangePush('update boundary domain1')
domain1.update_boundary()
nvtxRangePop()

nvtxRangePush('compute fluxes domain1')
domain1.compute_fluxes()
timestep1 = domain1.flux_timestep
boundary_flux1 = domain1.boundary_flux_sum[0]
nvtxRangePop()


#-----------------------------------------
# Test the kernel versions of
# distribute_to_vertices_and_edges
# update_boundary
# compute_fluxes
# as found in evolve_one_euler_step in
# generic_domain.py (abstract_2d_finite_volume)
#----------------------------------------

nvtxRangePush('distribute on cpu for domain2')
domain2.distribute_to_vertices_and_edges()
nvtxRangePop()

nvtxRangePush('update boundary domain2')
domain2.update_boundary()
nvtxRangePop()

# nvtxRangePush('initialise gpu_interface domain2')
# domain2.set_multiprocessor_mode(2)
# nvtxRangePop()

# nvtxRangePush('compute fluxes domain2')
# domain2.compute_fluxes()
# timestep2 = domain2.flux_timestep
# boundary_flux2 = domain2.boundary_flux_sum[0]
# nvtxRangePop()


# Setup gpu interface (if multiprocessor_mode == 4 and cupy available)
import pdb
pdb.set_trace()
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

nvtxRangePush('compute fluxes on gpu for domain2')
timestep2 = domain2.evolve_max_timestep 
timestep2 = gpu_interface2.compute_fluxes_ext_central_kernel(timestep2)
nvtxRangePop()

boundary_flux2 = domain2.boundary_flux_sum[0]

# Compare update arrays and timestep


print('domain1 timestep ', timestep1)
print('domain2 timestep ', timestep2)

print('domain1 boundary_flux ', boundary_flux1)
print('domain2 boundary_flux ', boundary_flux2)




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
# scale linalg.norm by number of elements
import math
sqrtN = 1.0/N


print('timestep diff                ', abs(timestep1-timestep2))
print('boundary_flux diff           ', abs(boundary_flux1-boundary_flux2))
print('max_speed diff L2-norm       ', num.linalg.norm(max_speed_1-max_speed_2)*sqrtN)
print('stage update diff L2-norm    ', num.linalg.norm(stage1.explicit_update-stage2.explicit_update)*sqrtN)
print('xmom  update diff L2-norm    ', num.linalg.norm(xmom1.explicit_update-xmom2.explicit_update)*sqrtN)
print('ymom  update diff L2-norm    ', num.linalg.norm(ymom1.explicit_update-ymom2.explicit_update)*sqrtN)

print('stage update diff Linf-norm  ', num.linalg.norm(stage1.explicit_update-stage2.explicit_update,num.inf))
print('xmom  update diff Linf-norm  ', num.linalg.norm(xmom1.explicit_update-xmom2.explicit_update,num.inf))
print('ymom  update diff Linf-norm  ', num.linalg.norm(ymom1.explicit_update-ymom2.explicit_update,num.inf))
#print('edge timestep error         ', num.linalg.norm(domain1.edge_timestep-domain2.edge_timestep))
#print('pressure work error         ', num.linalg.norm(domain1.pressuregrad_work-domain2.pressuregrad_work))
#print('edge flux work error        ', num.linalg.norm(domain1.edge_flux_work-domain2.edge_flux_work))


#(num.abs(stage1.explicit_update-stage2.explicit_update)/num.abs(stage1.explicit_update)).max()

from pprint import pprint

if False:
    pprint(stage1.explicit_update.reshape(2*nx,2*ny))
    pprint(stage2.explicit_update.reshape(2*nx,2*ny))
    pprint((stage1.explicit_update-stage2.explicit_update).reshape(2*nx,2*ny))
    pprint(max_speed_1.reshape(2*nx,2*ny))
    pprint(max_speed_2.reshape(2*nx,2*ny))
    pprint((max_speed_1-max_speed_2).reshape(2*nx,2*ny))


    stage_ids = num.argsort(num.abs(stage1.explicit_update-stage2.explicit_update))

    print('stage max diff values')
    pprint(stage_ids[-10:])
    pprint(stage1.explicit_update[stage_ids[-10:]])
    pprint(stage2.explicit_update[stage_ids[-10:]])
    pprint(num.abs(stage1.explicit_update-stage2.explicit_update)[stage_ids[-10:]])
    print(num.abs(stage1.explicit_update-stage2.explicit_update).max())


    xmom_ids = num.argsort(num.abs(xmom1.explicit_update-xmom2.explicit_update))

    print('xmom max diff values')
    pprint(xmom_ids[-10:])
    pprint(xmom1.explicit_update[xmom_ids[-10:]])
    pprint(xmom2.explicit_update[xmom_ids[-10:]])
    pprint(num.abs(xmom1.explicit_update-xmom2.explicit_update)[xmom_ids[-10:]])
    print(num.abs(xmom1.explicit_update-xmom2.explicit_update).max())








#assert num.allclose(timestep1,timestep2)
#assert num.allclose(boundary_flux1,boundary_flux2)
#assert num.allclose(stage1.explicit_update,stage2.explicit_update)
#assert num.allclose(xmom1.explicit_update,xmom2.explicit_update)
#assert num.allclose(ymom1.explicit_update,ymom2.explicit_update)
#assert num.allclose(domain1.edge_timestep,domain2.edge_timestep)
#assert num.allclose(domain1.pressuregrad_work,domain2.pressuregrad_work)
#assert num.allclose(domain1.edge_flux_work,domain2.edge_flux_work)



