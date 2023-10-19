
import numpy as np

import anuga



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


class GPU_interface(object):

    def __init__(self, domain):
    
        self.gpu_arrays_allocated = False

        #--------------------------------
        # create alias to domain variables
        #--------------------------------
        self.domain = domain

        self.cpu_number_of_elements  =  domain.number_of_elements
        self.cpu_number_of_boundaries  =  domain.number_of_boundaries
        self.cpu_boundary_length     =  domain.boundary_length 
        self.cpu_number_of_riverwall_edges =  domain.number_of_riverwall_edges
        self.cpu_epsilon             =  domain.epsilon
        self.cpu_H0                  =  domain.H0

        # FIXME SR: Why is this hard coded
        self.cpu_limiting_threshold  =  10.0 * domain.H0

        self.cpu_g                   =  domain.g
        self.cpu_optimise_dry_cells  =  domain.optimise_dry_cells
        self.cpu_evolve_max_timestep =  domain.evolve_max_timestep
        self.cpu_timestep_fluxcalls  =  domain.timestep_fluxcalls
        self.cpu_low_froude          =  domain.low_froude
        self.cpu_max_speed           =  domain.max_speed

        self.cpu_minimum_allowed_height =  domain.minimum_allowed_height
        self.cpu_maximum_allowed_speed  =  domain.maximum_allowed_speed
        self.cpu_extrapolate_velocity_second_order =  domain.extrapolate_velocity_second_order
        self.cpu_beta_w                 =  domain.beta_w
        self.cpu_beta_w_dry             =  domain.beta_w_dry
        self.cpu_beta_uh                =  domain.beta_uh
        self.cpu_beta_uh_dry            =  domain.beta_uh_dry
        self.cpu_beta_vh                =  domain.beta_vh
        self.cpu_beta_vh_dry            =  domain.beta_vh_dry
        self.cpu_max_flux_update_frequency =  domain.max_flux_update_frequency

        # Quantity structures
        quantities = domain.quantities
        stage  = quantities["stage"]
        xmom   = quantities["xmomentum"]
        ymom   = quantities["ymomentum"]
        bed    = quantities["elevation"]
        height = quantities["height"]

        riverwallData = domain.riverwallData

        self.cpu_riverwall_ncol_hydraulic_properties = riverwallData.ncol_hydraulic_properties

        
        self.cpu_stage_explicit_update  = stage.explicit_update   
        self.cpu_xmom_explicit_update   = xmom.explicit_update   
        self.cpu_ymom_explicit_update   = ymom.explicit_update
       
        self.cpu_stage_centroid_values  = stage.centroid_values
        self.cpu_xmom_centroid_values   = xmom.centroid_values
        self.cpu_ymom_centroid_values   = ymom.centroid_values
        self.cpu_height_centroid_values = height.centroid_values
        self.cpu_bed_centroid_values    = bed.centroid_values
              
        self.cpu_stage_edge_values      = stage.edge_values
        self.cpu_xmom_edge_values       = xmom.edge_values
        self.cpu_ymom_edge_values       = ymom.edge_values
        self.cpu_height_edge_values     = height.edge_values
        self.cpu_bed_edge_values        = bed.edge_values
       
        self.cpu_stage_boundary_values  = stage.boundary_values
        self.cpu_xmom_boundary_values   = xmom.boundary_values 
        self.cpu_ymom_boundary_values   = ymom.boundary_values

        self.cpu_stage_vertex_values    = stage.vertex_values
        self.cpu_height_vertex_values   = height.vertex_values
        self.cpu_xmom_vertex_values     = xmom.vertex_values
        self.cpu_ymom_vertex_values     = ymom.vertex_values
        self.cpu_bed_vertex_values      = bed.vertex_values


        self.cpu_domain_areas           = domain.areas
        self.cpu_domain_normals         = domain.normals
        self.cpu_domain_edgelengths     = domain.edgelengths
        self.cpu_domain_radii           = domain.radii
        self.cpu_domain_tri_full_flag   = domain.tri_full_flag
        self.cpu_domain_neighbours      = domain.neighbours
        self.cpu_domain_neighbour_edges = domain.neighbour_edges
        self.cpu_domain_edge_flux_type  = domain.edge_flux_type
        self.cpu_domain_edge_river_wall_counter = domain.edge_river_wall_counter
        self.cpu_riverwall_elevation    = riverwallData.riverwall_elevation
        self.cpu_riverwall_rowIndex     = riverwallData.hydraulic_properties_rowIndex
        self.cpu_riverwall_hydraulic_properties = riverwallData.hydraulic_properties

        self.cpu_centroid_coordinates   = domain.centroid_coordinates
        self.cpu_edge_coordinates       = domain.edge_coordinates
        self.cpu_surrogate_neighbours   = domain.surrogate_neighbours
        self.cpu_x_centroid_work        = domain.x_centroid_work
        self.cpu_y_centroid_work        = domain.y_centroid_work

        # FIXME SR: check whether these are arrays or scalars   
        self.cpu_beta_w_dry             = domain.beta_w_dry 
        self.cpu_beta_w                 = domain.beta_w
        self.cpu_beta_uh_dry            = domain.beta_uh_dry 
        self.cpu_beta_uh                = domain.beta_uh
        self.cpu_beta_vh_dry            = domain.beta_vh_dry
        self.cpu_beta_vh                = domain.beta_vh

        #----------------------------------------
        # Arrays to collect local timesteps and boundary fluxes
        # and do reduction after kernel call
        # FIXME SR: Maybe we only need a cupy array and use cupy amin 
        # and sum to collect global timestep and boundary flux
        #---------------------------------------
        self.cpu_local_boundary_flux_sum = np.zeros(self.cpu_number_of_elements, dtype=float) 
        self.cpu_timestep_array = np.zeros(self.cpu_number_of_elements, dtype=float) 



    def compile_gpu_kernels(self):
        """ 
        compile gpu kernels
        """

        import cupy as cp

        #----------------------------------------
        # Read in precompiled kernel function
        #----------------------------------------
        # create a Module object in python
        #mod = cp.cuda.function.Module()

        # load the cubin created by comiling ../cuda_anuga.cu 
        # with 
        # nvcc -arch=sm_70 -cubin -o cuda_anuga.cubin cuda_anuga.cu
        #mod.load_file("../cuda_anuga.cubin")

        # Locate cuda_anuga.cu code
        import anuga.shallow_water as sw
        import os        
        sw_dir = os.path.dirname(sw.__file__)
        cu_file = anuga.join(sw_dir,'cuda_anuga.cu')
        #with open('../cuda_anuga.cu') as f:
        with open(cu_file) as f:
            code = f.read()

        self.mod  = cp.RawModule(code=code, options=("--std=c++17",),
                                 name_expressions=("_cuda_compute_fluxes_loop",
                                                   "_cuda_extrapolate_second_order_edge_sw_loop1",
                                                   "_cuda_extrapolate_second_order_edge_sw_loop2",
                                                   "_cuda_extrapolate_second_order_edge_sw_loop3"))

        #FIXME SR: Only flux_kernel defined at present
        #FIXME SR: other kernels should be added to the file cuda_anuga.cu 
        self.flux_kernel = self.mod.get_function("_cuda_compute_fluxes_loop")

        self.extrapolate_kernel1 = self.mod.get_function("_cuda_extrapolate_second_order_edge_sw_loop1")
        self.extrapolate_kernel2 = self.mod.get_function("_cuda_extrapolate_second_order_edge_sw_loop2")
        self.extrapolate_kernel3 = self.mod.get_function("_cuda_extrapolate_second_order_edge_sw_loop3")

    #-----------------------------------------------------
    # Allocate GPU arrays
    #-----------------------------------------------------
    def allocate_gpu_arrays(self):

        import cupy as cp
        
        nvtxRangePush('to gpu')

        # FIXME SR: we should probably allocate all the cpu numpy arrays with 
        # pinned memory to speed movement of data from host to device

        self.gpu_number_of_boundaries = cp.array(self.cpu_number_of_boundaries)

        # these are just used for the reduction operation of the flux calculation
        self.gpu_timestep_array          = cp.array(self.cpu_timestep_array)
        self.gpu_local_boundary_flux_sum = cp.array(self.cpu_local_boundary_flux_sum )

        # these come out of flux calculation
        self.gpu_max_speed              = cp.array(self.cpu_max_speed)                #InOut

        # these come out of flux calculation
        # these go into the update conserved quantities
        self.gpu_stage_explicit_update  = cp.array(self.cpu_stage_explicit_update)    #InOut
        self.gpu_xmom_explicit_update   = cp.array(self.cpu_xmom_explicit_update)     #InOut
        self.gpu_ymom_explicit_update   = cp.array(self.cpu_ymom_explicit_update)     #InOut

        # centroid values go InOut of extrapolation
        # edge values go Out of extrapolation Into update_boundaries
        # vertex values go Out of extrapolation

        self.gpu_stage_centroid_values  = cp.array(self.cpu_stage_centroid_values) 
        self.gpu_stage_edge_values      = cp.array(self.cpu_stage_edge_values)

        self.gpu_xmom_centroid_values   = cp.array(self.cpu_xmom_centroid_values)
        self.gpu_xmom_edge_values       = cp.array(self.cpu_xmom_edge_values)

        self.gpu_ymom_centroid_values   = cp.array(self.cpu_ymom_centroid_values)        
        self.gpu_ymom_edge_values       = cp.array(self.cpu_ymom_edge_values)

        self.gpu_bed_centroid_values    = cp.array(self.cpu_bed_centroid_values)
        self.gpu_bed_edge_values        = cp.array(self.cpu_bed_edge_values)

        self.gpu_height_centroid_values = cp.array(self.cpu_height_centroid_values)
        self.gpu_height_edge_values     = cp.array(self.cpu_height_edge_values)
        
        self.gpu_stage_boundary_values  = cp.array(self.cpu_stage_boundary_values)  
        self.gpu_xmom_boundary_values   = cp.array(self.cpu_xmom_boundary_values) 
        self.gpu_ymom_boundary_values   = cp.array(self.cpu_ymom_boundary_values) 

        self.gpu_stage_vertex_values    = cp.array(self.cpu_stage_vertex_values)
        self.gpu_height_vertex_values   = cp.array(self.cpu_height_vertex_values)
        self.gpu_xmom_vertex_values     = cp.array(self.cpu_xmom_vertex_values)
        self.gpu_ymom_vertex_values     = cp.array(self.cpu_ymom_vertex_values)
        self.gpu_bed_vertex_values      = cp.array(self.cpu_bed_vertex_values)

        # These should not change during evolve
        self.gpu_areas                  = cp.array(self.cpu_domain_areas)
        self.gpu_normals                = cp.array(self.cpu_domain_normals)
        self.gpu_edgelengths            = cp.array(self.cpu_domain_edgelengths)
        self.gpu_radii                  = cp.array(self.cpu_domain_radii)
        self.gpu_tri_full_flag          = cp.array(self.cpu_domain_tri_full_flag)
        self.gpu_neighbours             = cp.array(self.cpu_domain_neighbours)
        self.gpu_neighbour_edges        = cp.array(self.cpu_domain_neighbour_edges)
        self.gpu_edge_flux_type         = cp.array(self.cpu_domain_edge_flux_type)  
        self.gpu_edge_river_wall_counter = cp.array(self.cpu_domain_edge_river_wall_counter)

        self.gpu_riverwall_elevation    = cp.array(self.cpu_riverwall_elevation)
        self.gpu_riverwall_rowIndex     = cp.array(self.cpu_riverwall_rowIndex)
        self.gpu_riverwall_hydraulic_properties = cp.array(self.cpu_riverwall_hydraulic_properties)

        self.gpu_centroid_coordinates   = cp.array(self.cpu_centroid_coordinates)
        self.gpu_edge_coordinates       = cp.array(self.cpu_edge_coordinates)
        self.gpu_surrogate_neighbours   = cp.array(self.cpu_surrogate_neighbours)              

        self.gpu_x_centroid_work        = cp.array(self.cpu_x_centroid_work)
        self.gpu_y_centroid_work        = cp.array(self.cpu_y_centroid_work)

        nvtxRangePop()

        self.gpu_arrays_allocated = True
    

    def cpu_to_gpu_centroid_values(self):
        """
        Move centroid data from cpu to gpu 
        """

        nvtxRangePush('cpu_to_gpu_centroid_values')
        # FIXME SR: Do we need to transfer height and bed centroid values
        self.gpu_stage_centroid_values.set(self.cpu_stage_centroid_values)
        self.gpu_xmom_centroid_values.set(self.cpu_xmom_centroid_values)
        self.gpu_ymom_centroid_values.set(self.cpu_ymom_centroid_values)
        self.gpu_height_centroid_values.set(self.cpu_height_centroid_values)
        self.gpu_bed_centroid_values.set(self.cpu_bed_centroid_values)
        nvtxRangePop()


    def gpu_to_cpu_centroid_values(self):
        """
        Move centroid data from gpu to cpu
        """

        import cupy as cp

        nvtxRangePush('gpu_to_cpu_centroid_values')
        # FIXME SR: Do we need to transfer height and bed centroid values
        cp.asnumpy(self.gpu_stage_centroid_values,  out = self.cpu_stage_centroid_values)
        cp.asnumpy(self.gpu_xmom_centroid_values,  out = self.cpu_xmom_centroid_values)
        cp.asnumpy(self.gpu_ymom_centroid_values,  out = self.cpu_ymom_centroid_values)
        cp.asnumpy(self.gpu_height_centroid_values,  out = self.cpu_height_centroid_values)
        cp.asnumpy(self.gpu_bed_centroid_values,  out = self.cpu_bed_centroid_values)  
        nvtxRangePop()

    def cpu_to_gpu_edge_values(self):
        """
        Move edge data from cpu to gpu 
        """
        
        nvtxRangePush('cpu_to_gpu_edge_values')
        self.gpu_stage_edge_values.set(self.cpu_stage_edge_values)
        self.gpu_xmom_edge_values.set(self.cpu_xmom_edge_values)
        self.gpu_ymom_edge_values.set(self.cpu_ymom_edge_values)
        self.gpu_bed_edge_values.set(self.cpu_bed_edge_values)
        self.gpu_height_edge_values.set(self.cpu_height_edge_values)
        nvtxRangePop()
        


    def gpu_to_cpu_edge_values(self):
        """
        Move centroid data from gpu to cpu
        """

        import cupy as cp

        nvtxRangePush('gpu_to_cpu_edge_values')
        cp.asnumpy(self.gpu_stage_edge_values,  out = self.cpu_stage_edge_values)
        cp.asnumpy(self.gpu_xmom_edge_values,  out = self.cpu_xmom_edge_values)
        cp.asnumpy(self.gpu_ymom_edge_values,  out = self.cpu_ymom_edge_values)
        cp.asnumpy(self.gpu_bed_edge_values,  out = self.cpu_bed_edge_values)
        cp.asnumpy(self.gpu_height_edge_values,  out = self.cpu_height_edge_values)
        nvtxRangePop()       


    def cpu_to_gpu_boundary_values(self):
        """
        Move boundary data from cpu to gpu 
        """
        import cupy as cp

        nvtxRangePush('cpu_to_gpu_boundary_values')
        self.gpu_stage_boundary_values.set(self.cpu_stage_boundary_values)  
        self.gpu_xmom_boundary_values.set(self.cpu_xmom_boundary_values) 
        self.gpu_ymom_boundary_values.set(self.cpu_ymom_boundary_values) 
        nvtxRangePop()

    def gpu_to_cpu_boundary_values(self):
        """
        Move boundary data from gpu to cpu
        """
        import cupy as cp
        nvtxRangePush('gpu_to_cpu_boundary_values')
        cp.asnumpy(self.gpu_stage_boundary_values, out = self.cpu_stage_boundary_values)  
        cp.asnumpy(self.gpu_xmom_boundary_values,  out = self.cpu_xmom_boundary_values) 
        cp.asnumpy(self.gpu_ymom_boundary_values,  out = self.cpu_ymom_boundary_values) 
        nvtxRangePop()

    def gpu_to_cpu_for_compute_fluxes(self):
        """
        Move compute_fluxes updated data from gpu to cpu
        """
        import cupy as cp
        nvtxRangePush('gpu_to_cpu_for_compute_fluxes')
        cp.asnumpy(self.gpu_max_speed, out = self.cpu_max_speed)
        cp.asnumpy(self.gpu_stage_explicit_update, out = self.cpu_stage_explicit_update)
        cp.asnumpy(self.gpu_xmom_explicit_update, out = self.cpu_xmom_explicit_update)
        cp.asnumpy(self.gpu_ymom_explicit_update, out = self.cpu_ymom_explicit_update)
        nvtxRangePop()


    def compute_fluxes_ext_central_kernel(self, timestep,  transfer_from_cpu=True, transfer_gpu_results=True):
        """
        compute flux

        Ensure transient data has been copied to the GPU via cpu_to_gpu_flux
        """
        import cupy as cp

        #------------------------------------------------
        # Transfer transient values from cpu if necessary
        #------------------------------------------------
        if transfer_from_cpu:
            nvtxRangePush('calculate flux: transfer from GPU')
            self.gpu_stage_centroid_values.set(self.cpu_stage_centroid_values)
            self.gpu_xmom_centroid_values.set(self.cpu_xmom_centroid_values)
            self.gpu_ymom_centroid_values.set(self.cpu_ymom_centroid_values)
            self.gpu_height_centroid_values.set(self.cpu_height_centroid_values)
            self.gpu_bed_centroid_values.set(self.cpu_bed_centroid_values)

            self.gpu_stage_edge_values.set(self.cpu_stage_edge_values)
            self.gpu_xmom_edge_values.set(self.cpu_xmom_edge_values)
            self.gpu_ymom_edge_values.set(self.cpu_ymom_edge_values)
            self.gpu_bed_edge_values.set(self.cpu_bed_edge_values)
            self.gpu_height_edge_values.set(self.cpu_height_edge_values)            
            
            #FIXME SR: Want about boundary values!
            nvtxRangePop()



        #-------------------------------------
        # FIXME SR: Need to calc substep_count
        # properly as used to be static in C code
        # probably could set substep_count in the
        # evolve procedure
        # dont need call in the simd version
        # of compute flux
        #-------------------------------------
        timestep_fluxcalls = 1
        call = 1

        timestep_fluxcalls = self.domain.timestep_fluxcalls
        base_call = call

        substep_count = (call - base_call) % self.domain.timestep_fluxcalls

        import math

        THREADS_PER_BLOCK = 128
        NO_OF_BLOCKS = int(math.ceil(self.cpu_number_of_elements/THREADS_PER_BLOCK))


        self.flux_kernel( (NO_OF_BLOCKS, 0, 0), 
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

                np.int64  (self.cpu_number_of_elements),
                np.int64  (substep_count),
                np.int64  (self.cpu_riverwall_ncol_hydraulic_properties),
                np.float64(self.cpu_epsilon),
                np.float64(self.cpu_g),
                np.int64  (self.cpu_low_froude),
                np.float64(self.cpu_limiting_threshold)
                ) 
                )

        nvtxRangePush('calculate flux: cupy reductions')
        # FIXME SR: Does gpu_reduce_timestep live on the GPU or CPU?
        gpu_reduce_timestep = self.gpu_timestep_array.min()
        gpu_reduced_local_boundary_flux_sum = self.gpu_local_boundary_flux_sum.sum()
        nvtxRangePop()

        nvtxRangePush('calculate flux: communicate reduced results')
        # do we need to do this?
        if substep_count == 0:
            timestep = cp.asnumpy(gpu_reduce_timestep)

        self.domain.boundary_flux_sum[substep_count] = cp.asnumpy(gpu_reduced_local_boundary_flux_sum)
        nvtxRangePop()

        #------------------------------------------------
        # Recover transient values from gpu if necessary
        #------------------------------------------------
        if transfer_gpu_results:
            nvtxRangePush('calculate flux: transfer from GPU')
            cp.asnumpy(self.gpu_max_speed, out = self.cpu_max_speed)
            cp.asnumpy(self.gpu_stage_explicit_update, out = self.cpu_stage_explicit_update)
            cp.asnumpy(self.gpu_xmom_explicit_update, out = self.cpu_xmom_explicit_update)
            cp.asnumpy(self.gpu_ymom_explicit_update, out = self.cpu_ymom_explicit_update)
            nvtxRangePop()

        return timestep


    def extrapolate_second_order_edge_sw_kernel(self, transfer_from_cpu=True, transfer_gpu_results=True, verbose=False):
        """
        compute extrapolation

        Ensure transient data has been copied to the GPU via cpu_to_gpu routines
        """
        import cupy as cp

        #------------------------------------------------
        # Transfer transient values from cpu if necessary
        #------------------------------------------------
        if transfer_from_cpu:
            nvtxRangePush('extrapolate kernel: transfer from cpu')
            self.gpu_stage_centroid_values.set(self.cpu_stage_centroid_values)
            self.gpu_xmom_centroid_values.set(self.cpu_xmom_centroid_values)
            self.gpu_ymom_centroid_values.set(self.cpu_ymom_centroid_values)
            self.gpu_height_centroid_values.set(self.cpu_height_centroid_values)
            self.gpu_bed_centroid_values.set(self.cpu_bed_centroid_values)  
            nvtxRangePop()

        import math
        THREADS_PER_BLOCK = 128
        NO_OF_BLOCKS = int(math.ceil(self.cpu_number_of_elements/THREADS_PER_BLOCK))

        nvtxRangePush('extrapolate kernel: run kernel')
        # FIXME SR: Check to see if we need to read in vertex_values arrays

        # self.extrapolate_kernel( (NO_OF_BLOCKS, 0, 0),
        #         (THREADS_PER_BLOCK, 0, 0), 
        #         (  
        #         self.gpu_stage_edge_values, 
        #         self.gpu_xmom_edge_values, 
        #         self.gpu_ymom_edge_values,
        #         self.gpu_height_edge_values, 
        #         self.gpu_bed_edge_values, 

        #         self.gpu_stage_centroid_values, 
        #         self.gpu_xmom_centroid_values,
        #         self.gpu_ymom_centroid_values, 
        #         self.gpu_height_centroid_values,
        #         self.gpu_bed_centroid_values,

        #         self.gpu_x_centroid_work,
        #         self.gpu_y_centroid_work,
                
        #         self.gpu_centroid_coordinates,
        #         self.gpu_edge_coordinates, 
        #         self.gpu_surrogate_neighbours,               
                
        #         self.gpu_beta_w_dry, 
        #         self.gpu_beta_w,
        #         self.gpu_beta_uh_dry, 
        #         self.gpu_beta_uh, 
        #         self.gpu_beta_vh_dry, 
        #         self.gpu_beta_vh,
                
        #         np.float64 (self.cpu_minimum_allowed_height), 
        #         np.int32   (self.cpu_number_of_elements), 
        #         np.int32 (self.cpu_extrapolate_velocity_second_order)
        #         ) 
        #         )

        

        self.extrapolate_kernel1( (NO_OF_BLOCKS, 0, 0),
                (THREADS_PER_BLOCK, 0, 0), 
                (  
                self.gpu_stage_edge_values, 
                self.gpu_xmom_edge_values, 
                self.gpu_ymom_edge_values,
                self.gpu_height_edge_values, 
                self.gpu_bed_edge_values, 

                self.gpu_stage_centroid_values, 
                self.gpu_xmom_centroid_values,
                self.gpu_ymom_centroid_values, 
                self.gpu_height_centroid_values,
                self.gpu_bed_centroid_values,

                self.gpu_x_centroid_work,
                self.gpu_y_centroid_work,
                
                self.gpu_centroid_coordinates,
                self.gpu_edge_coordinates, 
                self.gpu_surrogate_neighbours,               
                
                np.float64 (self.cpu_minimum_allowed_height), 
                np.int64   (self.cpu_number_of_elements), 
                np.int64   (self.cpu_extrapolate_velocity_second_order)
                ) 
                )

        if verbose:
            print('gpu_x_centroid_work after loop 1')
            print(self.gpu_x_centroid_work)
            print('gpu_x_centroid_work after loop 1')
            print(self.gpu_xmom_centroid_values)

        self.extrapolate_kernel2( (NO_OF_BLOCKS, 0, 0),
                (THREADS_PER_BLOCK, 0, 0), 
                (  
                self.gpu_stage_edge_values,
                self.gpu_xmom_edge_values, 
                self.gpu_ymom_edge_values,
                self.gpu_height_edge_values, 
                self.gpu_bed_edge_values, 

                self.gpu_stage_centroid_values, 
                self.gpu_xmom_centroid_values,
                self.gpu_ymom_centroid_values, 
                self.gpu_height_centroid_values,
                self.gpu_bed_centroid_values,

                self.gpu_x_centroid_work,
                self.gpu_y_centroid_work,

                self.gpu_stage_vertex_values,
                self.gpu_height_vertex_values,
                self.gpu_xmom_vertex_values,
                self.gpu_ymom_vertex_values,
                self.gpu_bed_vertex_values,
                
                self.gpu_number_of_boundaries,
                self.gpu_centroid_coordinates,
                self.gpu_edge_coordinates, 
                self.gpu_surrogate_neighbours,               
                
                np.float64 (self.cpu_beta_w_dry), 
                np.float64 (self.cpu_beta_w),
                np.float64 (self.cpu_beta_uh_dry), 
                np.float64 (self.cpu_beta_uh), 
                np.float64 (self.cpu_beta_vh_dry), 
                np.float64 (self.cpu_beta_vh),
                
                np.float64 (self.cpu_minimum_allowed_height), 
                np.int64   (self.cpu_number_of_elements), 
                np.int64   (self.cpu_extrapolate_velocity_second_order)
                ) 
                )

        self.extrapolate_kernel3( (NO_OF_BLOCKS, 0, 0),
                (THREADS_PER_BLOCK, 0, 0), 
                (   
                self.gpu_xmom_centroid_values,
                self.gpu_ymom_centroid_values, 

                self.gpu_x_centroid_work,
                self.gpu_y_centroid_work,
                
                np.int64 (self.cpu_extrapolate_velocity_second_order),
                np.int64 (self.cpu_number_of_elements)
                ) 
                )


        nvtxRangePop()

        #------------------------------------------------
        # Recover transient values from gpu if necessary
        #------------------------------------------------
        if transfer_gpu_results:
            nvtxRangePush('extrapolate kernel: transfer gpu results')
            
            cp.asnumpy(self.gpu_stage_edge_values,  out = self.cpu_stage_edge_values)
            cp.asnumpy(self.gpu_xmom_edge_values,   out = self.cpu_xmom_edge_values)
            cp.asnumpy(self.gpu_ymom_edge_values,   out = self.cpu_ymom_edge_values)
            cp.asnumpy(self.gpu_height_edge_values, out = self.cpu_height_edge_values)
            cp.asnumpy(self.gpu_bed_edge_values,    out = self.cpu_bed_edge_values)

            # FIXME SR: check to see if we need to transfer all these centroid values
            cp.asnumpy(self.gpu_stage_centroid_values,  out = self.cpu_stage_centroid_values)
            cp.asnumpy(self.gpu_xmom_centroid_values,   out = self.cpu_xmom_centroid_values)
            cp.asnumpy(self.gpu_ymom_centroid_values,   out = self.cpu_ymom_centroid_values)
            cp.asnumpy(self.gpu_height_centroid_values, out = self.cpu_height_centroid_values)
            cp.asnumpy(self.gpu_bed_centroid_values,    out = self.cpu_bed_centroid_values)

            nvtxRangePop()






