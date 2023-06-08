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


""" Run a version of the validation test runup_sinusoid
to ensure limiting solution has small velocity
"""

def create_domain(name='domain'):
    domain = anuga.rectangular_cross_domain(2,2, len1=1., len2=1.)

    domain.set_flow_algorithm('DE0')
    domain.set_low_froude(0)

    domain.set_name(name)  
    domain.set_datadir('.')

    #------------------
    # Define topography
    #------------------
    scale_me=1.0

    def topography(x,y):
        return (-x/2.0 +0.05*num.sin((x+y)*50.0))*scale_me

    def stagefun(x,y):
        stage=-0.2*scale_me #+0.01*(x>0.9)
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


domain1 = create_domain('domain_original')
domain1.set_multiprocessor_mode(0)

domain2 = create_domain('domain_cuda')
domain2.set_multiprocessor_mode(0) # will change to 2 once burn in

#------------------------------
#Evolve the system through time
#------------------------------
print('Evolve domain1')
for t in domain1.evolve(yieldstep=0.1,finaltime=0.1):
    domain1.print_timestepping_statistics()

print('Evolve domain2')
for t in domain2.evolve(yieldstep=0.1,finaltime=0.1):
    domain2.print_timestepping_statistics()

#----------------------------------------
# Now just run the cuda code on domain2
#----------------------------------------
domain2.set_multiprocessor_mode(4)
timestep = 0.1

domain1.distribute_to_vertices_and_edges()

domain1.compute_fluxes()
timestep1 = domain1.flux_timestep
boundary_flux1 = domain1.boundary_flux_sum[0]

def compute_fluxes_ext_central_kernel(domain,timestep):
    
    local_timestep = num.zeros((1,), dtype=float)     # InOut

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
    minimum_allowed_height =  domain.minimum_allowed_height
    maximum_allowed_speed =  domain.maximum_allowed_speed
    timestep_fluxcalls =  domain.timestep_fluxcalls
    low_froude =  domain.low_froude
    extrapolate_velocity_second_order =  domain.extrapolate_velocity_second_order
    beta_w =  domain.beta_w
    beta_w_dry =  domain.beta_w_dry
    beta_uh =  domain.beta_uh
    beta_uh_dry =  domain.beta_uh_dry
    beta_vh =  domain.beta_vh
    beta_vh_dry =  domain.beta_vh_dry
    max_flux_update_frequency =  domain.max_flux_update_frequency

    neighbours = domain.neighbours
    surrogate_neighbours = domain.surrogate_neighbours
    neighbour_edges = domain.neighbour_edges
    normals = domain.normals
    edgelengths = domain.edgelengths
    radii = domain.radii
    areas = domain.areas
    edge_flux_type = domain.edge_flux_type
    tri_full_flag = domain.tri_full_flag

    already_computed_flux = domain.already_computed_flux

    vertex_coordinates = domain.vertex_coordinates
    edge_coordinates = domain.edge_coordinates
    centroid_coordinates = domain.centroid_coordinates
    max_speed = domain.max_speed
    number_of_boundaries = domain.number_of_boundaries
    flux_update_frequency = domain.flux_update_frequency
    update_next_flux = domain.update_next_flux
    update_extrapolation = domain.update_extrapolation
    allow_timestep_increase = domain.allow_timestep_increase
    edge_timestep = domain.edge_timestep
    edge_flux_work = domain.edge_flux_work
    neigh_work = domain.neigh_work
    pressuregrad_work = domain.pressuregrad_work
    x_centroid_work = domain.x_centroid_work
    y_centroid_work = domain.y_centroid_work
    boundary_flux_sum = domain.boundary_flux_sum
    edge_river_wall_counter = domain.edge_river_wall_counter

    # Quantity structures
    quantities = domain.quantities
    stage = quantities["stage"]
    xmomentum = quantities["xmomentum"]
    ymomentum = quantities["ymomentum"]
    elevation = quantities["elevation"]
    height = quantities["height"]

    stage_edge_values = stage.edge_values
    xmom_edge_values = xmomentum.edge_values
    ymom_edge_values = ymomentum.edge_values
    bed_edge_values = elevation.edge_values
    height_edge_values = height.edge_values
    stage_centroid_values = stage.centroid_values
    xmom_centroid_values = xmomentum.centroid_values
    ymom_centroid_values = ymomentum.centroid_values
    bed_centroid_values = elevation.centroid_values
    height_centroid_values = height.centroid_values
    stage_vertex_values = stage.vertex_values
    xmom_vertex_values = xmomentum.vertex_values
    ymom_vertex_values = ymomentum.vertex_values
    bed_vertex_values = elevation.vertex_values
    height_vertex_values = height.vertex_values
    stage_boundary_values = stage.boundary_values
    xmom_boundary_values = xmomentum.boundary_values
    ymom_boundary_values = ymomentum.boundary_values
    bed_boundary_values = elevation.boundary_values
    stage_explicit_update = stage.explicit_update
    xmom_explicit_update = stage.explicit_update
    ymom_explicit_update = ymomentum.explicit_update
    #------------------------------------------------------
    # Riverwall structures
    #------------------------------------------------------
    riverwallData = domain.riverwallData
    riverwall_elevation = riverwallData.riverwall_elevation
    riverwall_rowIndex = riverwallData.hydraulic_properties_rowIndex
    ncol_riverwall_hydraulic_properties = riverwallData.ncol_hydraulic_properties
    riverwall_hydraulic_properties = riverwallData.hydraulic_properties



    #-------------------------------------
    # glue code
    #-------------------------------------
    timestep_fluxcalls = 1
    call = 1

    timestep_fluxcalls = domain.timestep_fluxcalls
    base_call = call

    substep_count = (call - base_call) % domain.timestep_fluxcalls

    #------------------------------------
    # create cupy arrays
    #------------------------------------
    import cupy as cp
    
    gpu_local_timestep        = cp.array(local_timestep)         #InOut
    gpu_boundary_flux_sum     = cp.array(boundary_flux_sum )     #InOut
    gpu_max_speed             = cp.array(max_speed)              #InOut
    gpu_stage_explicit_update = cp.array(stage_explicit_update)  #InOut
    gpu_xmom_explicit_update  = cp.array(xmom_explicit_update)   #InOut
    gpu_ymom_explicit_update  = cp.array(ymom_explicit_update)   #InOut

    gpu_stage_centroid_values = cp.array(stage_centroid_values) 
    gpu_stage_edge_values     = cp.array(stage_edge_values)
    gpu_xmom_edge_values      = cp.array(xmom_edge_values)
    gpu_ymom_edge_values      = cp.array(ymom_edge_values)
    gpu_bed_edge_values       = cp.array(bed_edge_values)
    gpu_height_edge_values    = cp.array(height_edge_values)
    gpu_height_centroid_values = cp.array(height_centroid_values)
    gpu_bed_centroid_values   = cp.array(bed_centroid_values)
    gpu_stage_boundary_values = cp.array(stage_boundary_values)  
    gpu_xmom_boundary_values  = cp.array(xmom_boundary_values) 
    gpu_ymom_boundary_values  = cp.array(ymom_boundary_values) 
    gpu_areas                 = cp.array(areas)
    gpu_normals               = cp.array(normals)
    gpu_edgelengths           = cp.array(edgelengths)
    gpu_radii                 = cp.array(radii)
    gpu_tri_full_flag         = cp.array(tri_full_flag)
    gpu_neighbours            = cp.array(neighbours)
    gpu_neighbour_edges       = cp.array(neighbour_edges)
    gpu_edge_flux_type        = cp.array(edge_flux_type)  
    gpu_edge_river_wall_counter = cp.array(edge_river_wall_counter)
    gpu_riverwall_elevation   = cp.array(riverwall_elevation)
    gpu_riverwall_rowIndex    = cp.array(riverwall_rowIndex)
    gpu_riverwall_hydraulic_properties = cp.array(riverwall_hydraulic_properties)


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

    kernel = cp.RawKernel(code=code)


    # call the function with a tuple of grid size, a tuple of block size, 
    # and a tuple of all arguments required by the kernel
    # if the kernel requires shared memory, append `shared_mem=n_bytes` to the function call


    import math

    THREADS_PER_BLOCK = 128
    NO_OF_BLOCKS = int(math.ceil(number_of_elements/THREADS_PER_BLOCK))


    kernel( (NO_OF_BLOCKS, 0, 0), 
                                 (THREADS_PER_BLOCK, 0, 0), 
                                  (
                                    gpu_local_timestep, 
                                    gpu_boundary_flux_sum, 
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
                                    ncol_riverwall_hydraulic_properties,
                                    epsilon,
                                    g,
                                    low_froude,
                                    limiting_threshold  
                                 )
                                 )


    #-------------------------------------
    # Recover values from gpu
    #-------------------------------------

    local_timestep[:]        = cp.as_numpy(gpu_local_timestep)          #InOut
    boundary_flux_sum[:]     = cp.as_numpy(gpu_boundary_flux_sum)       #InOut
    max_speed[:]             = cp.as_numpy(gpu_max_speed)               #InOut
    stage_explicit_update[:] = cp.as_numpy(gpu_stage_explicit_update)   #InOut
    xmom_explicit_update[:]  = cp.as_numpy(gpu_xmom_explicit_update)    #InOut
    ymom_explicit_update[:]  = cp.as_numpy(gpu_ymom_explicit_update)    #InOut


    return local_timestep[0]

#-----------------------------------------
# Test the kernel version of compute fluxes
#----------------------------------------
domain2.distribute_to_vertices_and_edges()

timestep = domain2.evolve_max_timestep
flux_timestep = compute_fluxes_ext_central_kernel(domain2, timestep)
domain2.flux_timestep = flux_timestep



domain2.compute_fluxes()
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


quantities2 = domain2.quantities
stage2 = quantities2["stage"]
xmom2 = quantities2["xmomentum"]
ymom2 = quantities2["ymomentum"]


print('timestep error              ', abs(timestep1-timestep2))
print('boundary_flux error         ', abs(boundary_flux1-boundary_flux2))
print('stage explicit update error ', num.linalg.norm(stage1.explicit_update-stage2.explicit_update))
print('xmom  explicit update error ', num.linalg.norm(xmom1.explicit_update-xmom2.explicit_update))
print('ymom  explicit update error ', num.linalg.norm(ymom1.explicit_update-ymom2.explicit_update))
#print('edge timestep error         ', num.linalg.norm(domain1.edge_timestep-domain2.edge_timestep))
#print('pressure work error         ', num.linalg.norm(domain1.pressuregrad_work-domain2.pressuregrad_work))
#print('edge flux work error        ', num.linalg.norm(domain1.edge_flux_work-domain2.edge_flux_work))



#assert num.allclose(timestep1,timestep2)
#assert num.allclose(boundary_flux1,boundary_flux2)
#assert num.allclose(stage1.explicit_update,stage2.explicit_update)
#assert num.allclose(xmom1.explicit_update,xmom2.explicit_update)
#assert num.allclose(ymom1.explicit_update,ymom2.explicit_update)
#assert num.allclose(domain1.edge_timestep,domain2.edge_timestep)
#assert num.allclose(domain1.pressuregrad_work,domain2.pressuregrad_work)
#assert num.allclose(domain1.edge_flux_work,domain2.edge_flux_work)



