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

#-----------------------------------------------------
# Code for profiling cuda version
#-----------------------------------------------------
def nvtxRangePush(*arg):
    pass
def nvtxRangePop(*arg):
    pass

try:
    from cupy.cuda.nvtx import RangePush as nvtxRangePush
    from cupy.cuda.nvtx import RangePop  as nvtxRangePop
except:
    pass

try:
    from nvtx import range_push as nvtxRangePush
    from nvtx import range_pop  as nvtxRangePop
except:
    pass

nx = 2
ny = 2

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

#----------------------------------------
# Now just run the cuda code on domain2
#----------------------------------------
domain2.set_multiprocessor_mode(2)
timestep = 0.1

nvtxRangePush('distribute domain1')
domain1.distribute_to_vertices_and_edges()
nvtxRangePop()

nvtxRangePush('compute fluxes domain1')
domain1.compute_fluxes()
timestep1 = domain1.flux_timestep
boundary_flux1 = domain1.boundary_flux_sum[0]
nvtxRangePop()



def compute_fluxes_ext_central_kernel(domain, timestep):
    

    #--------------------------------
    # create alias to domain variables
    #--------------------------------
    number_of_elements =  domain.number_of_elements
    boundary_length =  domain.boundary_length 
    number_of_riverwall_edges =  domain.number_of_riverwall_edges
    epsilon =  domain.epsilon
    H0 =  domain.H0
    g =  domain.g
    optimise_dry_cells =  domain.optimise_dry_cells
    evolve_max_timestep =  domain.evolve_max_timestep
    timestep_fluxcalls =  domain.timestep_fluxcalls
    low_froude =  domain.low_froude

    # minimum_allowed_height =  domain.minimum_allowed_height
    # maximum_allowed_speed =  domain.maximum_allowed_speed
    # extrapolate_velocity_second_order =  domain.extrapolate_velocity_second_order
    # beta_w =  domain.beta_w
    # beta_w_dry =  domain.beta_w_dry
    # beta_uh =  domain.beta_uh
    # beta_uh_dry =  domain.beta_uh_dry
    # beta_vh =  domain.beta_vh
    # beta_vh_dry =  domain.beta_vh_dry
    # max_flux_update_frequency =  domain.max_flux_update_frequency
    limiting_threshold = 10.0 * domain.H0

    # Quantity structures
    quantities = domain.quantities
    stage  = quantities["stage"]
    xmom   = quantities["xmomentum"]
    ymom   = quantities["ymomentum"]
    bed    = quantities["elevation"]
    height = quantities["height"]
    riverwallData = domain.riverwallData

    riverwall_ncol_hydraulic_properties = riverwallData.ncol_hydraulic_properties

    #-------------------------------------
    # FIXME SR: Need to calc substep_count
    # properly as used to be static in C code
    #-------------------------------------
    timestep_fluxcalls = 1
    call = 1

    timestep_fluxcalls = domain.timestep_fluxcalls
    base_call = call

    substep_count = (call - base_call) % domain.timestep_fluxcalls

    #----------------------------------------
    # Arrays to collect local timesteps and boundary fluxes
    # and do reduction after kernel call
    # FIXME SR: Maybe we only need a cupy array and use cupy amin 
    # and sum to collect global timestep and boundary flux
    #---------------------------------------
    local_boundary_flux_sum = num.zeros(number_of_elements, dtype=float) 
    timestep_array = num.zeros(number_of_elements, dtype=float) 

    #------------------------------------
    # create cupy arrays
    #------------------------------------
    import cupy as cp
    
    nvtxRangePush('to gpu')

    gpu_timestep_array          = cp.array(timestep_array)           #InOut
    gpu_local_boundary_flux_sum = cp.array(local_boundary_flux_sum ) #InOut

    gpu_max_speed              = cp.array(domain.max_speed)         #InOut
    gpu_stage_explicit_update  = cp.array(stage.explicit_update)    #InOut
    gpu_xmom_explicit_update   = cp.array(xmom.explicit_update)     #InOut
    gpu_ymom_explicit_update   = cp.array(ymom.explicit_update)     #InOut

    gpu_stage_centroid_values  = cp.array(stage.centroid_values) 
    gpu_stage_edge_values      = cp.array(stage.edge_values)
    gpu_xmom_edge_values       = cp.array(xmom.edge_values)
    gpu_ymom_edge_values       = cp.array(ymom.edge_values)
    gpu_bed_edge_values        = cp.array(bed.edge_values)
    gpu_height_edge_values     = cp.array(height.edge_values)
    gpu_height_centroid_values = cp.array(height.centroid_values)
    gpu_bed_centroid_values    = cp.array(bed.centroid_values)
    gpu_stage_boundary_values  = cp.array(stage.boundary_values)  
    gpu_xmom_boundary_values   = cp.array(xmom.boundary_values) 
    gpu_ymom_boundary_values   = cp.array(ymom.boundary_values) 
    gpu_areas                  = cp.array(domain.areas)
    gpu_normals                = cp.array(domain.normals)
    gpu_edgelengths            = cp.array(domain.edgelengths)
    gpu_radii                  = cp.array(domain.radii)
    gpu_tri_full_flag          = cp.array(domain.tri_full_flag)
    gpu_neighbours             = cp.array(domain.neighbours)
    gpu_neighbour_edges        = cp.array(domain.neighbour_edges)
    gpu_edge_flux_type         = cp.array(domain.edge_flux_type)  
    gpu_edge_river_wall_counter = cp.array(domain.edge_river_wall_counter)

    gpu_riverwall_elevation    = cp.array(riverwallData.riverwall_elevation)
    gpu_riverwall_rowIndex     = cp.array(riverwallData.hydraulic_properties_rowIndex)
    gpu_riverwall_hydraulic_properties = cp.array(riverwallData.hydraulic_properties)

    nvtxRangePop()

    #----------------------------------------
    # Read in precompiled kernel function
    #----------------------------------------
    # create a Module object in python
    #mod = cp.cuda.function.Module()

    # load the cubin created by comiling ../cuda_anuga.cu 
    # with 
    # nvcc -arch=sm_70 -cubin -o cuda_anuga.cubin cuda_anuga.cu
    #mod.load_file("../cuda_anuga.cubin")


    # fetch the kernel to make it a Python function object
    #_cuda_compute_fluxes_loop_1 = mod.get_function("_cuda_compute_fluxes_loop_1")



    with open('../cuda_anuga.cu') as f:
        code = f.read()

    mod  = cp.RawModule(code=code, options=("--std=c++17",), name_expressions=("_cuda_compute_fluxes_loop_1",))

    kernel = mod.get_function("_cuda_compute_fluxes_loop_1")

    # call the function with a tuple of grid size, a tuple of block size, 
    # and a tuple of all arguments required by the kernel
    # if the kernel requires shared memory, append `shared_mem=n_bytes` to the function call


    import math

    THREADS_PER_BLOCK = 128
    NO_OF_BLOCKS = int(math.ceil(number_of_elements/THREADS_PER_BLOCK))


    kernel( (NO_OF_BLOCKS, 0, 0), 
                                 (THREADS_PER_BLOCK, 0, 0), 
                                  ( 
                                    gpu_timestep_array, 
                                    gpu_local_boundary_flux_sum, 

                                    gpu_max_speed, 
                                    gpu_stage_explicit_update,
                                    gpu_xmom_explicit_update,
                                    gpu_ymom_explicit_update,

                                    gpu_stage_centroid_values,
                                    gpu_stage_edge_values,
                                    gpu_xmom_edge_values, 
                                    gpu_ymom_edge_values,
                                    gpu_bed_edge_values,
                                    gpu_height_edge_values,
                                    gpu_height_centroid_values,
                                    gpu_bed_centroid_values,
                                    gpu_stage_boundary_values, 
                                    gpu_xmom_boundary_values, 
                                    gpu_ymom_boundary_values, 
                                    gpu_areas,
                                    gpu_normals,
                                    gpu_edgelengths,
                                    gpu_radii,
                                    gpu_tri_full_flag,
                                    gpu_neighbours,
                                    gpu_neighbour_edges,
                                    gpu_edge_flux_type, 
                                    gpu_edge_river_wall_counter,

                                    gpu_riverwall_elevation,
                                    gpu_riverwall_rowIndex,
                                    gpu_riverwall_hydraulic_properties,

                                    number_of_elements,
                                    substep_count,
                                    riverwall_ncol_hydraulic_properties,
                                    epsilon,
                                    g,
                                    low_froude,
                                    limiting_threshold 
                                 ) 
                                 )


    #-------------------------------------
    # Recover values from gpu
    #-------------------------------------

    #print('=================')
    #print('boundary_flux_sum', boundary_flux_sum)
    #print('gpu_boundary_flux_sum', gpu_boundary_flux_sum)


    nvtxRangePush('calculate flux: from gpu')

    timestep_array[:]        = cp.asnumpy(gpu_timestep_array)          #InOut
    local_boundary_flux_sum[:] = cp.asnumpy(gpu_local_boundary_flux_sum) #InOut
    domain.max_speed[:]      = cp.asnumpy(gpu_max_speed)               #InOut
    stage.explicit_update[:] = cp.asnumpy(gpu_stage_explicit_update)   #InOut
    xmom.explicit_update[:]  = cp.asnumpy(gpu_xmom_explicit_update)    #InOut
    ymom.explicit_update[:]  = cp.asnumpy(gpu_ymom_explicit_update)    #InOut

    nvtxRangePop()


    nvtxRangePush('calculate flux: reduction operations')
    
    if substep_count == 0:
        timestep = timestep_array.min()

    domain.boundary_flux_sum[substep_count] = local_boundary_flux_sum.sum()

    nvtxRangePop()


    #print('boundary_flux_sum', boundary_flux_sum)
    #print('gpu_boundary_flux_sum', gpu_boundary_flux_sum)
    #print('=================')

    return timestep

#-----------------------------------------
# Test the kernel version of compute fluxes
#----------------------------------------
nvtxRangePush('distribute domain2')
domain2.distribute_to_vertices_and_edges()
nvtxRangePop()

nvtxRangePush('compute fluxes domain2')
#domain2.compute_fluxes()
timestep = domain2.evolve_max_timestep 
domain2.flux_timestep = compute_fluxes_ext_central_kernel(domain2, timestep)
nvtxRangePop()


## 
timestep2 = domain2.flux_timestep
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
sqrtN = 1.0/math.sqrt(N)


print('timestep error              ', abs(timestep1-timestep2))
print('boundary_flux error         ', abs(boundary_flux1-boundary_flux2))
print('max_speed error             ', num.linalg.norm(max_speed_1-max_speed_2))
print('stage explicit update error ', num.linalg.norm(stage1.explicit_update-stage2.explicit_update)*sqrtN)
print('xmom  explicit update error ', num.linalg.norm(xmom1.explicit_update-xmom2.explicit_update)*sqrtN)
print('ymom  explicit update error ', num.linalg.norm(ymom1.explicit_update-ymom2.explicit_update)*sqrtN)
#print('edge timestep error         ', num.linalg.norm(domain1.edge_timestep-domain2.edge_timestep))
#print('pressure work error         ', num.linalg.norm(domain1.pressuregrad_work-domain2.pressuregrad_work))
#print('edge flux work error        ', num.linalg.norm(domain1.edge_flux_work-domain2.edge_flux_work))


from pprint import pprint

if False:
    pprint(stage1.explicit_update.reshape(2*nx,2*ny))
    pprint(stage2.explicit_update.reshape(2*nx,2*ny))
    pprint((stage1.explicit_update-stage2.explicit_update).reshape(2*nx,2*ny))
    pprint(max_speed_1.reshape(2*nx,2*ny))
    pprint(max_speed_2.reshape(2*nx,2*ny))
    pprint((max_speed_1-max_speed_2).reshape(2*nx,2*ny))
#assert num.allclose(timestep1,timestep2)
#assert num.allclose(boundary_flux1,boundary_flux2)
#assert num.allclose(stage1.explicit_update,stage2.explicit_update)
#assert num.allclose(xmom1.explicit_update,xmom2.explicit_update)
#assert num.allclose(ymom1.explicit_update,ymom2.explicit_update)
#assert num.allclose(domain1.edge_timestep,domain2.edge_timestep)
#assert num.allclose(domain1.pressuregrad_work,domain2.pressuregrad_work)
#assert num.allclose(domain1.edge_flux_work,domain2.edge_flux_work)



