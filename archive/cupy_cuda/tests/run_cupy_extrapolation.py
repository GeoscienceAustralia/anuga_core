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
    #print(domain.__dict__)
    #print(dir(domain))
    return domain


print('')
print(70*'=')
print('Test Runup for extrapolate')
print(70*'=')

nvtxRangePush('create domain1')
domain1 = create_domain('domain_original')
domain1.set_multiprocessor_mode(1)

quantities1 = domain1.quantities
stage1 = quantities1["stage"]
xmom1 = quantities1["xmomentum"]
ymom1 = quantities1["ymomentum"]
nvtxRangePop()

nvtxRangePush('create domain2')
domain2 = create_domain('domain_cuda')
domain2.set_multiprocessor_mode(1) # will change to 4 once burn in

quantities2 = domain2.quantities
stage2 = quantities2["stage"]
xmom2 = quantities2["xmomentum"]
ymom2 = quantities2["ymomentum"]
nvtxRangePop()

import time
start = time.time()

#------------------------------
#Evolve the system through time
#------------------------------
yieldstep = 0.02
finaltime = 0.02
nvtxRangePush('evolve domain1')
print('Evolve domain1')
print('domain1 number of triangles ',domain1.number_of_elements)
for t in domain1.evolve(yieldstep=yieldstep,finaltime=finaltime):
    domain1.print_timestepping_statistics()
nvtxRangePop()


#---------------------------------------
# run domain1 using standard routine
#---------------------------------------
#timestep = 0.1

nvtxRangePush('distribute domain1')
domain1.distribute_to_vertices_and_edges()
nvtxRangePop()

#nvtxRangePush('update boundary domain1')
#domain1.update_boundary()
#nvtxRangePop()

# nvtxRangePush('extrapolate domain1')
# vol1=domain1.get_water_volume()
# boundaryFluxInt1=domain1.get_boundary_flux_integral()
# timestep1 = domain1.flux_timestep
# boundary_flux1 = domain1.boundary_flux_sum[0]
# nvtxRangePop()

end = time.time()
print('DOMAIN 1 time ' + str(end - start))



#-----------------------------------------
# Test the kernel version of compute fluxes
#----------------------------------------

start = time.time()


nvtxRangePush('evolve domain2')
print('Evolve domain2')
print('domain2 number of triangles ',domain2.number_of_elements)
for t in domain2.evolve(yieldstep=yieldstep,finaltime=finaltime):
    domain2.print_timestepping_statistics()
nvtxRangePop()

stage1_vertex_values_before = num.copy(stage1.vertex_values)
stage2_vertex_values_before = num.copy(stage2.vertex_values)
stage1_edge_values_before = num.copy(stage1.edge_values)
stage2_edge_values_before = num.copy(stage2.edge_values)

nvtxRangePush('distribute domain2')
# Now run the distribute procedure on the GPU
domain2.set_multiprocessor_mode(2)
domain2.distribute_to_vertices_and_edges()
nvtxRangePop()


# from anuga.shallow_water.sw_domain_cuda import GPU_interface
# gpu_domain2 = GPU_interface(domain2)

# nvtxRangePush('allocate gpu arrays domain2')
# gpu_domain2.allocate_gpu_arrays()
# nvtxRangePop()

# nvtxRangePush('compile gpu kernels domain2')
# gpu_domain2.compile_gpu_kernels()
# nvtxRangePop()

# nvtxRangePush('distribute_to_vertices_and_edges on gpu domain2')
# gpu_domain2.extrapolate_second_order_edge_sw_kernel()
# nvtxRangePop()



end = time.time()
print('DOMAIN 2 time ' + str(end - start))



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

print('stage edge      diff L2 norm ', num.linalg.norm(stage1.edge_values-stage2.edge_values)/N)
print('xmom  edge      diff L2 norm ', num.linalg.norm(xmom1.edge_values-xmom2.edge_values)/N)
print('ymom  edge      diff L2 norm ', num.linalg.norm(ymom1.edge_values-ymom2.edge_values)/N)

print('stage centroid diff L2 norm ', num.linalg.norm(stage1.centroid_values-stage2.centroid_values)/N)
print('xmom  centroid diff L2 norm ', num.linalg.norm(xmom1.centroid_values-xmom2.centroid_values)/N)
print('ymom  centroid diff L2 norm ', num.linalg.norm(ymom1.centroid_values-ymom2.centroid_values)/N)

print('stage vertex diff L2 norm ', num.linalg.norm(stage1.vertex_values-stage2.vertex_values)/N)
print('xmom  vertex diff L2 norm ', num.linalg.norm(xmom1.vertex_values-xmom2.vertex_values)/N)
print('ymom  vertex diff L2 norm ', num.linalg.norm(ymom1.vertex_values-ymom2.vertex_values)/N)



# FIXME SR: Why are these equal? I didn't think the vertex values had been copied back to the cpu


stage1_vertex_values_after = num.copy(stage1.vertex_values)
stage2_vertex_values_after = num.copy(stage2.vertex_values)
stage1_edge_values_after = num.copy(stage1.edge_values)
stage2_edge_values_after = num.copy(stage2.edge_values)



print("change stage1.vertex_values")
pprint(stage1_vertex_values_after - stage1_vertex_values_before)

print("change stage2.vertex_values")
pprint(stage2_vertex_values_after - stage2_vertex_values_before)

print("change stage1.edge_values")
pprint(stage1_edge_values_after - stage1_edge_values_before)

print("change stage2.edge_values")
pprint(stage2_edge_values_after - stage2_edge_values_before)


# print('timestep error                ', abs(timestep1-timestep2))
# print('boundary_flux error           ', abs(boundary_flux1-boundary_flux2))
# print('max_speed L2error             ', num.linalg.norm(max_speed_1-max_speed_2)*sqrtN)
# print('stage explicit update L2error ', num.linalg.norm(stage1.explicit_update-stage2.explicit_update)*sqrtN)
# print('xmom  explicit update L2error ', num.linalg.norm(xmom1.explicit_update-xmom2.explicit_update)*sqrtN)
# print('ymom  explicit update L2error ', num.linalg.norm(ymom1.explicit_update-ymom2.explicit_update)*sqrtN)

# print('stage update inferror         ', num.linalg.norm(stage1.explicit_update-stage2.explicit_update,num.inf))
# print('xmom  update inferror         ', num.linalg.norm(xmom1.explicit_update-xmom2.explicit_update,num.inf))
# print('ymom  update inferror         ', num.linalg.norm(ymom1.explicit_update-ymom2.explicit_update,num.inf))
#print('edge timestep error         ', num.linalg.norm(domain1.edge_timestep-domain2.edge_timestep))
#print('pressure work error         ', num.linalg.norm(domain1.pressuregrad_work-domain2.pressuregrad_work))
#print('edge flux work error        ', num.linalg.norm(domain1.edge_flux_work-domain2.edge_flux_work))


#(num.abs(stage1.explicit_update-stage2.explicit_update)/num.abs(stage1.explicit_update)).max()

# from pprint import pprint

# if False:
#     pprint(stage1.explicit_update.reshape(2*nx,2*ny))
#     pprint(stage2.explicit_update.reshape(2*nx,2*ny))
#     pprint((stage1.explicit_update-stage2.explicit_update).reshape(2*nx,2*ny))
#     pprint(max_speed_1.reshape(2*nx,2*ny))
#     pprint(max_speed_2.reshape(2*nx,2*ny))
#     pprint((max_speed_1-max_speed_2).reshape(2*nx,2*ny))


# stage_ids = num.argsort(num.abs(stage1.explicit_update-stage2.explicit_update))

# print('stage max diff values')
# pprint(stage_ids[-10:])
# pprint(stage1.explicit_update[stage_ids[-10:]])
# pprint(stage2.explicit_update[stage_ids[-10:]])
# pprint(num.abs(stage1.explicit_update-stage2.explicit_update)[stage_ids[-10:]])
# print(num.abs(stage1.explicit_update-stage2.explicit_update).max())


# xmom_ids = num.argsort(num.abs(xmom1.explicit_update-xmom2.explicit_update))

# print('xmom max diff values')
# pprint(xmom_ids[-10:])
# pprint(xmom1.explicit_update[xmom_ids[-10:]])
# pprint(xmom2.explicit_update[xmom_ids[-10:]])
# pprint(num.abs(xmom1.explicit_update-xmom2.explicit_update)[xmom_ids[-10:]])
# print(num.abs(xmom1.explicit_update-xmom2.explicit_update).max())


# print('vol1 ', vol1)
# print('vol2 ', vol2)

# import numpy as np

# # Check if 'vol' and 'vol2' are close with a relative tolerance of 5%
# if np.allclose(vol1, vol2, rtol=0.05):
#     print("vol1 and vol2 are close.")
# else:
#     print("vol1 and vol2 are not close.")

# # Check if 'vol' and 'boundaryFluxInt' are close
# if np.allclose(vol1, boundaryFluxInt1):
#     print("vol1 and boundaryFluxInt1 are close.")
# else:
#     print("vol1 and boundaryFluxInt1 are not close.")

# # Check if 'vol2' and 'boundaryFluxInt2' are close
# if np.allclose(vol2, boundaryFluxInt2):
#     print("vol2 and boundaryFluxInt2 are close.")
# else:
#     print("vol2 and boundaryFluxInt2 are not close.")

# # Check the absolute difference between 'domain.quantities['stage'].centroid_values' and 'domain2.quantities['stage'].centroid_values'
# abs_diff = np.abs(domain1.quantities['stage'].centroid_values - domain2.quantities['stage'].centroid_values)
# if np.all(abs_diff < 0.02):
#     print("Absolute difference is less than 0.02.")
# else:
#     print("Absolute difference is greater than or equal to 0.02.")
