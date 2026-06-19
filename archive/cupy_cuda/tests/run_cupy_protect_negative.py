import anuga
import os
from anuga import Reflective_boundary
from anuga import rectangular_cross_domain

from anuga import Domain

import numpy as num
import warnings
import time
import math
from pprint import pprint

from anuga.shallow_water.sw_domain_cuda import nvtxRangePush, nvtxRangePop


nx = 10
ny = 10

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
    #print(domain.__dict__)
    #print(dir(domain))
    return domain


print('')
print(70*'=')
print('Test Runup for update_conserved_quantities')
print(70*'=')
import time
start = time.time()

nvtxRangePush('create domain1')
domain1 = create_domain('domain_original')
domain1.set_multiprocessor_mode(1)
nvtxRangePop()

nvtxRangePush('create domain1')
domain2 = create_domain('domain_cuda')
domain2.set_multiprocessor_mode(1) # will change to 2 once burn in
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

end_time =time.time()
print('total creation time', (end_time - start))

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

#---------------------------------------
# run domain1 using standard routine
#---------------------------------------
timestep = 0.1



#------------------------------
#Evolve the system through time
#------------------------------
nvtxRangePush('evolve domain2')
print('Evolve domain2')
print('domain2 number of triangles ',domain2.number_of_elements)
for t in domain2.evolve(yieldstep=yieldstep,finaltime=finaltime):
    domain2.print_timestepping_statistics()
nvtxRangePop()

#---------------------------------------
# run domain2 using standard routine
#---------------------------------------
timestep = 0.1

# nvtxRangePush('distribute domain1')
# domain2.distribute_to_vertices_and_edges()
# nvtxRangePop()

# nvtxRangePush('update boundary domain1')
# domain2.update_boundary()
# nvtxRangePop()


nvtxRangePush('initialise gpu_interface : domain2')
domain2.set_multiprocessor_mode(2)
nvtxRangePop()

# import pdb; pdb.set_trace()
timestep2 = domain2.timestep
from anuga.shallow_water.sw_domain_cuda import GPU_interface
gpu_domain2 = GPU_interface(domain2)

nvtxRangePush('allocate gpu arrays for domain2')
gpu_domain2.allocate_gpu_arrays()
nvtxRangePop()

nvtxRangePush('compile gpu kernels for domain2')
gpu_domain2.compile_gpu_kernels()
nvtxRangePop()


stage1_centroid_values_before = num.copy(stage1.centroid_values)
stage2_centroid_values_before = num.copy(stage2.centroid_values)
xmom1_centroid_values_before = num.copy(xmom1.centroid_values)
xmom2_centroid_values_before = num.copy(xmom2.centroid_values)
ymom1_centroid_values_before = num.copy(ymom1.centroid_values)
ymom2_centroid_values_before = num.copy(ymom2.centroid_values)


nvtxRangePush('protect_kernal : domain2')

protect_new = gpu_domain2.protect_against_infinitesimal_and_negative_heights_kernal
mass_error = protect_new()

nvtxRangePop()



N = domain1.number_of_elements
# scale linalg.norm by number of elements
import math
sqrtN = 1.0/N


# print('stage edge      diff L2 norm ', num.linalg.norm(stage1.edge_values-stage2.edge_values)/N)
# print('xmom  edge      diff L2 norm ', num.linalg.norm(xmom1.edge_values-xmom2.edge_values)/N)
# print('ymom  edge      diff L2 norm ', num.linalg.norm(ymom1.edge_values-ymom2.edge_values)/N)

print('stage centroid diff L2 norm ', num.linalg.norm(stage1.centroid_values-stage2.centroid_values)/N)
print('stage vertex diff L2 norm ', num.linalg.norm(stage1.vertex_values-stage2.vertex_values)/N)


# print('stage semi implicit update diff L2 norm ', num.linalg.norm(stage1.semi_implicit_update-stage2.semi_implicit_update)/N)
# print('xmom  semi implicit update diff L2 norm ', num.linalg.norm(xmom1.semi_implicit_update-xmom2.semi_implicit_update)/N)
# print('ymom  semi implicit update diff L2 norm ', num.linalg.norm(ymom1.semi_implicit_update-ymom2.semi_implicit_update)/N)

# print('stage explicit update diff L2 norm ', num.linalg.norm(stage1.explicit_update-stage2.explicit_update)/N)
# print('stage explicit update diff L2 norm ', num.linalg.norm(stage1.explicit_update-stage2.explicit_update)/N)
# print('stage explicit update diff L2 norm ', num.linalg.norm(stage1.explicit_update-stage2.explicit_update)/N)


stage1_centroid_values_after = num.copy(stage1.centroid_values)
stage2_centroid_values_after = num.copy(stage2.centroid_values)
xmom1_centroid_values_after = num.copy(xmom1.centroid_values)
xmom2_centroid_values_after = num.copy(xmom2.centroid_values)
ymom1_centroid_values_after = num.copy(ymom1.centroid_values)
ymom2_centroid_values_after = num.copy(ymom2.centroid_values)


# print("change stage1.centroid_values")
# pprint(stage1_centroid_values_after - stage1_centroid_values_before)

# print("change stage2.centroid_values")
# pprint(stage2_centroid_values_after - stage2_centroid_values_before)

# print("change xmom1.centroid_values")
# pprint(xmom1_centroid_values_after - xmom1_centroid_values_before)

# print("change xmom2.centroid_values")
# pprint(xmom2_centroid_values_after - xmom2_centroid_values_before)

# print()
# print()
# print()

# print('Stage 1 and stage 2 cenrtoid values before')
# pprint(stage1_centroid_values_before)
# pprint(stage2_centroid_values_before)
