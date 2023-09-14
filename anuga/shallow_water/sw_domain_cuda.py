
import numpy as num

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


class Gpu_interface(object):

    def __init__(self, domain):
    
        self.gpu_arrays_allocated = False

        #--------------------------------
        # create alias to domain variables
        #--------------------------------
        self.domain = domain

        self.number_of_elements  =  domain.number_of_elements
        self.boundary_length     =  domain.boundary_length 
        self.number_of_riverwall_edges =  domain.number_of_riverwall_edges
        self.epsilon             =  domain.epsilon
        self.H0                  =  domain.H0
        self.limiting_threshold  =  10.0 * domain.H0
        self.g                   =  domain.g
        self.optimise_dry_cells  =  domain.optimise_dry_cells
        self.evolve_max_timestep =  domain.evolve_max_timestep
        self.timestep_fluxcalls  =  domain.timestep_fluxcalls
        self.low_froude          =  domain.low_froude

        self.minimum_allowed_height =  domain.minimum_allowed_height
        self.maximum_allowed_speed  =  domain.maximum_allowed_speed
        self.extrapolate_velocity_second_order =  domain.extrapolate_velocity_second_order
        self.beta_w                 =  domain.beta_w
        self.beta_w_dry             =  domain.beta_w_dry
        self.beta_uh                =  domain.beta_uh
        self.beta_uh_dry            =  domain.beta_uh_dry
        self.beta_vh                =  domain.beta_vh
        self.beta_vh_dry            =  domain.beta_vh_dry
        self.max_flux_update_frequency =  domain.max_flux_update_frequency

        # Quantity structures
        quantities = domain.quantities
        stage  = quantities["stage"]
        xmom   = quantities["xmomentum"]
        ymom   = quantities["ymomentum"]
        bed    = quantities["elevation"]
        height = quantities["height"]
        riverwallData = domain.riverwallData

        self.riverwall_ncol_hydraulic_properties = riverwallData.ncol_hydraulic_properties

        #----------------------------------------
        # Arrays to collect local timesteps and boundary fluxes
        # and do reduction after kernel call
        # FIXME SR: Maybe we only need a cupy array and use cupy amin 
        # and sum to collect global timestep and boundary flux
        #---------------------------------------
        self.local_boundary_flux_sum = num.zeros(self.number_of_elements, dtype=float) 
        self.timestep_array = num.zeros(self.number_of_elements, dtype=float) 


    #-----------------------------------------------------
    # compile gpu kernels
    #-----------------------------------------------------
    def compile_gpu_kernels(self):

        #----------------------------------------
        # Read in precompiled kernel function
        #----------------------------------------
        # create a Module object in python
        #mod = cp.cuda.function.Module()

        # load the cubin created by comiling ../cuda_anuga.cu 
        # with 
        # nvcc -arch=sm_70 -cubin -o cuda_anuga.cubin cuda_anuga.cu
        #mod.load_file("../cuda_anuga.cubin")

        # FIXME SR: need to ensure we find location of kernel source 
        with open('cuda_anuga.cu') as f:
            code = f.read()

        mod  = cp.RawModule(code=code, options=("--std=c++17",), name_expressions=("_cuda_compute_fluxes_loop",))

        self.flux_kernel = mod.get_function("_cuda_compute_fluxes_loop")



    
    #-----------------------------------------------------
    # Allocate GPU arrays
    #-----------------------------------------------------
    def allocate_gpu_arrays(self):
    
        #------------------------------------
        # create cupy arrays
        #------------------------------------
        import cupy as cp
        
        nvtxRangePush('to gpu')

        domain = self.domain

        quantities = domain.quantities
        stage  = quantities["stage"]
        xmom   = quantities["xmomentum"]
        ymom   = quantities["ymomentum"]
        bed    = quantities["elevation"]
        height = quantities["height"]


        riverwallData = domain.riverwallData

        # FIXME SR: we should probably allocate all these numpy arrays with 
        # pinned memory to speed movement of data from host to device

        self.gpu_timestep_array          = cp.array(timestep_array)           #InOut
        self.gpu_local_boundary_flux_sum = cp.array(local_boundary_flux_sum ) #InOut

        self.gpu_max_speed              = cp.array(domain.max_speed)         #InOut
        self.gpu_stage_explicit_update  = cp.array(stage.explicit_update)    #InOut
        self.gpu_xmom_explicit_update   = cp.array(xmom.explicit_update)     #InOut
        self.gpu_ymom_explicit_update   = cp.array(ymom.explicit_update)     #InOut

        self.gpu_stage_centroid_values  = cp.array(stage.centroid_values) 
        self.gpu_stage_edge_values      = cp.array(stage.edge_values)
        self.gpu_xmom_edge_values       = cp.array(xmom.edge_values)
        self.gpu_ymom_edge_values       = cp.array(ymom.edge_values)
        self.gpu_bed_edge_values        = cp.array(bed.edge_values)
        self.gpu_height_edge_values     = cp.array(height.edge_values)
        self.gpu_height_centroid_values = cp.array(height.centroid_values)
        self.gpu_bed_centroid_values    = cp.array(bed.centroid_values)
        self.gpu_stage_boundary_values  = cp.array(stage.boundary_values)  
        self.gpu_xmom_boundary_values   = cp.array(xmom.boundary_values) 
        self.gpu_ymom_boundary_values   = cp.array(ymom.boundary_values) 
        self.gpu_areas                  = cp.array(domain.areas)
        self.gpu_normals                = cp.array(domain.normals)
        self.gpu_edgelengths            = cp.array(domain.edgelengths)
        self.gpu_radii                  = cp.array(domain.radii)
        self.gpu_tri_full_flag          = cp.array(domain.tri_full_flag)
        self.gpu_neighbours             = cp.array(domain.neighbours)
        self.gpu_neighbour_edges        = cp.array(domain.neighbour_edges)
        self.gpu_edge_flux_type         = cp.array(domain.edge_flux_type)  
        self.gpu_edge_river_wall_counter = cp.array(domain.edge_river_wall_counter)

        self.gpu_riverwall_elevation    = cp.array(riverwallData.riverwall_elevation)
        self.gpu_riverwall_rowIndex     = cp.array(riverwallData.hydraulic_properties_rowIndex)
        self.gpu_riverwall_hydraulic_properties = cp.array(riverwallData.hydraulic_properties)

        nvtxRangePop()

        self.gpu_arrays_allocated = True
    
    #-----------------------------------------------------
    # compute flux
    #-----------------------------------------------------
    def compute_fluxes_ext_central_kernel(self, timestep):
   
        domain = self.domain

        quantities = domain.quantities
        stage  = quantities["stage"]
        xmom   = quantities["xmomentum"]
        ymom   = quantities["ymomentum"]
        bed    = quantities["elevation"]
        height = quantities["height"]

        #-------------------------------------
        # FIXME SR: Need to calc substep_count
        # properly as used to be static in C code
        #-------------------------------------
        timestep_fluxcalls = 1
        call = 1

        timestep_fluxcalls = self.domain.timestep_fluxcalls
        base_call = call

        substep_count = (call - base_call) % self.domain.timestep_fluxcalls

        import math

        THREADS_PER_BLOCK = 128
        NO_OF_BLOCKS = int(math.ceil(self.number_of_elements/THREADS_PER_BLOCK))


        flux_kernel( (NO_OF_BLOCKS, 0, 0), 
                (THREADS_PER_BLOCK, 0, 0), 
                (  
                self.gpu_timestep_array, 
                self.gpu_local_boundary_flux_sum, 

                self.gpu_max_speed, 
                self.gpu_stage_explicit_update,
                self.gpu_xmom_explicit_update,
                self.gpu_ymom_explicit_update,

                self.gpu_stage_centroid_values,
                self.gpu_stage_edge_values,
                self.gpu_xmom_edge_values, 
                self.gpu_ymom_edge_values,
                self.gpu_bed_edge_values,
                self.gpu_height_edge_values,
                self.gpu_height_centroid_values,
                self.gpu_bed_centroid_values,
                self.gpu_stage_boundary_values, 
                self.gpu_xmom_boundary_values, 
                self.gpu_ymom_boundary_values, 
                self.gpu_areas,
                self.gpu_normals,
                self.gpu_edgelengths,
                self.gpu_radii,
                self.gpu_tri_full_flag,
                self.gpu_neighbours,
                self.gpu_neighbour_edges,
                self.gpu_edge_flux_type, 
                self.gpu_edge_river_wall_counter,

                self.gpu_riverwall_elevation,
                self.gpu_riverwall_rowIndex,
                self.gpu_riverwall_hydraulic_properties,

                num.int64(self.number_of_elements),
                num.int64(substep_count),
                num.int64(self.riverwall_ncol_hydraulic_properties),
                num.float64(self.epsilon),
                num.float64(self.g),
                num.int64(self.low_froude),
                num.float64(self.limiting_threshold)
                ) 
                )


        #-------------------------------------
        # Recover values from gpu
        #-------------------------------------

        nvtxRangePush('calculate flux: cupy reductions')
        # FIXME SR: Does gpu_reduce_timestep live on the GPU or CPU?
        gpu_reduce_timestep = self.gpu_timestep_array.min()

        gpu_reduced_local_boundary_flux_sum = self.gpu_local_boundary_flux_sum.sum()
        nvtxRangePop()


        nvtxRangePush('calculate flux: transfer from GPU')
    
        # cp.asnumpy(gpu_timestep_array,          out = timestep_array)          #InOut
        # cp.asnumpy(gpu_local_boundary_flux_sum, out = local_boundary_flux_sum) #InOut
        
        cp.asnumpy(gpu_max_speed, out = domain.max_speed)                    #InOut
        cp.asnumpy(gpu_stage_explicit_update, out = stage.explicit_update)   #InOut
        cp.asnumpy(gpu_xmom_explicit_update, out = xmom.explicit_update)     #InOut
        cp.asnumpy(gpu_ymom_explicit_update, out = ymom.explicit_update)     #InOut

        nvtxRangePop()


        nvtxRangePush('calculate flux: communicate reduced results')
        # do we need to do this?
        if substep_count == 0:
            timestep = cp.asnumpy(gpu_reduce_timestep)

        domain.boundary_flux_sum[substep_count] = cp.asnumpy(gpu_reduced_local_boundary_flux_sum)

        nvtxRangePop()

        return timestep

