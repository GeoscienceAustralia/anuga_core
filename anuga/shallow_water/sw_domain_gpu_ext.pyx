#cython: wraparound=False, boundscheck=False, cdivision=True, profile=False, nonecheck=False, overflowcheck=False, cdivision_warnings=False, unraisable_tracebacks=False, warn.unused=False
"""
Cython wrapper for GPU-accelerated shallow water solver with MPI support.

This module bridges Python (mpi4py, NumPy) to the C GPU implementation.
Key responsibilities:
- Extract MPI_Comm handle from mpi4py
- Convert Python domain arrays to C pointers
- Flatten Python ghost dicts to C arrays for GPU efficiency
"""

import cython
from libc.stdint cimport int64_t, int32_t, uint64_t
from libc.stdlib cimport malloc, free
from libc.string cimport memset

import numpy as np
cimport numpy as np

# Import mpi4py's MPI_Comm type when building with MPI support; otherwise fall
# back to a plain-int MPI_Comm matching the C single-process stubs
# (gpu/gpu_mpi_stubs.h). HAVE_MPI4PY is a Cython compile-time env constant set by
# meson based on whether MPI + mpi4py were found at build time.
IF HAVE_MPI4PY:
    from mpi4py cimport MPI
    from mpi4py.libmpi cimport MPI_Comm
ELSE:
    ctypedef int MPI_Comm

# External C declarations (from gpu/ subdirectory)
cdef extern from "gpu_domain.h" nogil:
    # Forward declare the structs
    struct domain:
        int64_t number_of_elements
        int64_t boundary_length
        double epsilon
        double H0
        double g
        double minimum_allowed_height
        double maximum_allowed_speed
        double evolve_max_timestep
        int64_t low_froude
        int64_t reconstruct_edge_bed
        int64_t num_owned_edges
        int64_t* owned_edges
        int64_t extrapolate_velocity_second_order
        # Beta values for gradient limiting
        double beta_w
        double beta_w_dry
        double beta_uh
        double beta_uh_dry
        double beta_vh
        double beta_vh_dry
        # Centroid values
        double* stage_centroid_values
        double* xmom_centroid_values
        double* ymom_centroid_values
        double* bed_centroid_values
        double* height_centroid_values
        # Edge values
        double* stage_edge_values
        double* xmom_edge_values
        double* ymom_edge_values
        double* bed_edge_values
        double* height_edge_values
        # Updates
        double* stage_explicit_update
        double* xmom_explicit_update
        double* ymom_explicit_update
        double* stage_semi_implicit_update
        double* xmom_semi_implicit_update
        double* ymom_semi_implicit_update
        # Boundary values
        double* stage_boundary_values
        double* xmom_boundary_values
        double* ymom_boundary_values
        double* bed_boundary_values
        double* height_boundary_values
        # Backup for RK2
        double* stage_backup_values
        double* xmom_backup_values
        double* ymom_backup_values
        # Mesh connectivity
        int64_t* neighbours
        int64_t* neighbour_edges
        int64_t* surrogate_neighbours
        int64_t* number_of_boundaries
        double* normals
        double* edgelengths
        double* areas
        double* radii
        double* max_speed
        double* centroid_coordinates
        double* edge_coordinates
        # Boundary flux accumulator (written by core_compute_fluxes_central,
        # read by the Python boundary_flux_integral_operator)
        double* boundary_flux_sum
        # Work arrays for extrapolation
        double* x_centroid_work
        double* y_centroid_work
        # Friction
        double* friction_centroid_values
        # Ghost cell flag (1=full/owned, 0=ghost); NULL for single-process domains
        int64_t* tri_full_flag
        # Riverwall arrays
        int64_t number_of_riverwall_edges
        int64_t number_of_tracers
        double beta_tracer
        int64_t n_sediment_classes
        double* sediment_settling_velocity
        double* sediment_d_star
        double sediment_c_max
        double* sediment_diameter
        double* sediment_R
        double* sediment_tau_c_star
        double sediment_gamma0
        int64_t sediment_erosion_mode
        double sediment_K_partheniades
        double sediment_dl_A
        double sediment_dl_alpha
        double sediment_dl_beta
        double sediment_dl_threshold
        double sediment_dl_ks
        double sediment_dl_ks_factor
        int64_t sediment_vegetation_shear
        int64_t sediment_deposition_mode
        double sediment_tau_d
        double sediment_tau_crit
        double sediment_K_e
        double sediment_rho_w
        int64_t sediment_shear_closure
        double sediment_max_slope
        int64_t sediment_slope_frozen
        double* tracer_external_source
        int64_t sediment_d_star_mode
        double* sediment_reference_height
        double sediment_a_h_floor
        int64_t sediment_rouse_beta_mode
        double sediment_rouse_scale
        int64_t sediment_adaptation_mode
        double sediment_adaptation_alpha
        int64_t sediment_nearbed_base
        double sediment_layer_fraction
        double sediment_exchange_factor
        int64_t sediment_velocity_profile
        double sediment_porosity
        double sediment_morphological_factor
        int64_t sediment_bed_evolution
        int64_t sediment_bedload_mode
        double sediment_bedload_K
        double sediment_bedload_m
        double sediment_bedload_tau_c_star
        double sediment_bedload_h_min
        double* sediment_z_base
        double* sediment_shear_factor
        int64_t sediment_has_z_base
        int64_t sediment_bed_composition
        double sediment_active_layer
        double* sediment_bed_active
        double* sediment_bed_substrate
        double* sediment_repose_dz
        double sediment_repose_tan
        double sediment_repose_relax
        int64_t sediment_repose_max_sweeps
        double* sediment_source_limited
        double* sediment_slope_work
        int64_t* sediment_bed_exhausted
        int64_t* sediment_bedload_open
        double* sediment_bedload_supply
        double* sediment_qbx
        double* sediment_qba
        double* sediment_qby
        double sediment_c_pack
        int64_t vegetation_mode
        int64_t use_sloped_mannings
        double vegetation_Cd
        double vegetation_bed_chezy
        double* veg_density_centroid_values
        double* veg_diameter_centroid_values
        double* veg_height_centroid_values
        int64_t sediment_friction_mode
        double sediment_manning_ll
        int64_t sediment_wilson_bed
        double sediment_wilson_D
        double* tracer_centroid_values
        double* tracer_edge_values
        double* tracer_boundary_values
        double* tracer_explicit_update
        double* tracer_conserved_values
        double* tracer_backup_values
        double* tracer_boundary_flux
        double* tracer_speed_factor
        double* tracer_boundary_flux_sum
        int64_t ncol_riverwall_hydraulic_properties
        int64_t nrow_riverwall_hydraulic_properties
        int64_t* edge_flux_type
        int64_t* edge_river_wall_counter
        double* riverwall_elevation
        int64_t* riverwall_rowIndex
        double* riverwall_hydraulic_properties

    struct halo_exchange:
        int num_neighbors
        int* neighbor_ranks
        int* send_counts
        int* recv_counts
        int total_send_size
        int total_recv_size
        int* flat_send_indices
        int* flat_recv_indices
        double* send_buffer
        double* recv_buffer

    struct inlet_operator_info:
        int num_indices
        int* indices
        double* areas
        double total_area
        double* scratch_stages
        double* scratch_bed
        double* scratch_xmom
        double* scratch_ymom
        double* scratch_depths
        int active
        int mapped

    struct inlet_operators:
        inlet_operator_info *ops
        int num_operators
        int capacity

    struct max_quantities_info:
        int n
        double velocity_zero_height
        double *max_stage
        double *max_depth
        double *max_speed
        double *max_uh
        double *max_tracer
        int n_tracers
        int initialized
        int mapped

    struct gpu_domain:
        domain D
        MPI_Comm comm
        int rank
        int nprocs
        int gpu_initialized
        int device_id
        int gpu_aware_mpi
        int verbose
        char last_error[512]
        halo_exchange halo
        inlet_operators inlet_ops
        max_quantities_info max_qty
        double CFL
        double evolve_max_timestep
        double fixed_flux_timestep
        double recorded_flux_timestep
        int use_sloped_mannings
        int use_active_set
        long as_samples

    # Function declarations - initialization and cleanup
    int gpu_domain_init(gpu_domain *GD, MPI_Comm comm, int rank, int nprocs)
    void gpu_domain_finalize(gpu_domain *GD)
    int gpu_halo_init(gpu_domain *GD, int num_neighbors, int *neighbor_ranks,
                      int *send_counts, int *recv_counts,
                      int *flat_send_indices, int *flat_recv_indices)
    void gpu_halo_finalize(gpu_domain *GD)

    # GPU memory management
    size_t gpu_estimate_required_memory(int64_t n, int64_t nb)
    int    gpu_check_device_memory(gpu_domain *GD)
    int    gpu_domain_map_arrays(gpu_domain *GD)
    void gpu_remap_boundary_arrays(gpu_domain *GD)
    void gpu_domain_unmap_arrays(gpu_domain *GD)
    void gpu_domain_sync_to_device(gpu_domain *GD)
    void gpu_domain_sync_from_device(gpu_domain *GD)
    void gpu_domain_sync_all_from_device(gpu_domain *GD)
    void gpu_sync_boundary_values(gpu_domain *GD)
    void gpu_sync_riverwall_to_device(gpu_domain *GD)
    void gpu_sync_tracer_source_to_device(gpu_domain *GD)
    void gpu_sync_tracer_boundary_to_device(gpu_domain *GD)
    void gpu_sync_edge_values_from_device(gpu_domain *GD)
    int gpu_boundary_edge_sync_init(gpu_domain *GD, int num_boundary_cells, int *boundary_cell_ids)
    void gpu_boundary_edge_sync_finalize(gpu_domain *GD)
    void gpu_boundary_edge_sync(gpu_domain *GD)

    # MPI ghost exchange
    void gpu_exchange_ghosts(gpu_domain *GD)
    void gpu_exchange_bed_ghosts(gpu_domain *GD)

    # Reflective boundary
    int gpu_reflective_init(gpu_domain *GD, int num_edges,
                            int *boundary_indices, int *vol_ids, int *edge_ids)
    void gpu_reflective_finalize(gpu_domain *GD)
    void gpu_evaluate_reflective_boundary(gpu_domain *GD)

    # Dirichlet boundary
    int gpu_dirichlet_init(gpu_domain *GD, int num_edges,
                           int *boundary_indices, int *vol_ids, int *edge_ids,
                           double *stage_values, double *xmom_values, double *ymom_values)
    void gpu_dirichlet_finalize(gpu_domain *GD)
    void gpu_evaluate_dirichlet_boundary(gpu_domain *GD)

    # Transmissive boundary
    int gpu_transmissive_init(gpu_domain *GD, int num_edges,
                              int *boundary_indices, int *vol_ids, int *edge_ids,
                              int use_centroid)
    void gpu_transmissive_finalize(gpu_domain *GD)
    void gpu_evaluate_transmissive_boundary(gpu_domain *GD)

    # Transmissive_n_momentum_zero_t_momentum_set_stage boundary
    int gpu_transmissive_n_zero_t_init(gpu_domain *GD, int num_edges,
                                       int *boundary_indices, int *vol_ids, int *edge_ids)
    void gpu_transmissive_n_zero_t_finalize(gpu_domain *GD)
    void gpu_transmissive_n_zero_t_set_stage(gpu_domain *GD, double stage_value)
    void gpu_evaluate_transmissive_n_zero_t_boundary(gpu_domain *GD)

    # File_boundary / Field_boundary - spatially varying time-dependent values
    int  gpu_file_boundary_init(gpu_domain *GD, int num_edges,
                                int *boundary_indices, int *vol_ids, int *edge_ids)
    void gpu_file_boundary_finalize(gpu_domain *GD)
    void gpu_file_boundary_set_values(gpu_domain *GD,
                                      double *stage, double *xmom, double *ymom)
    void gpu_evaluate_file_boundary(gpu_domain *GD)

    # Time_boundary - time-dependent Dirichlet values
    int gpu_time_boundary_init(gpu_domain *GD, int num_edges,
                               int *boundary_indices, int *vol_ids, int *edge_ids)
    void gpu_time_boundary_finalize(gpu_domain *GD)
    void gpu_time_boundary_set_values(gpu_domain *GD, double *stage, double *xmom, double *ymom)
    void gpu_evaluate_time_boundary(gpu_domain *GD)

    # Absorbing_wave_boundary
    int  gpu_absorbing_wave_init(gpu_domain *GD, int num_edges,
                                 int *boundary_indices, int *vol_ids, int *edge_ids)
    void gpu_absorbing_wave_finalize(gpu_domain *GD)
    void gpu_absorbing_wave_set_value(gpu_domain *GD, double wave_value)
    void gpu_evaluate_absorbing_wave_boundary(gpu_domain *GD)

    # Characteristic_wave_boundary
    int  gpu_characteristic_wave_init(gpu_domain *GD, int num_edges,
                                      int *boundary_indices, int *vol_ids, int *edge_ids,
                                      double background_stage)
    void gpu_characteristic_wave_finalize(gpu_domain *GD)
    void gpu_characteristic_wave_set_value(gpu_domain *GD, double wave_value)
    void gpu_evaluate_characteristic_wave_boundary(gpu_domain *GD)

    # Flather_external_stage_zero_velocity_boundary
    int  gpu_flather_init(gpu_domain *GD, int num_edges,
                          int *boundary_indices, int *vol_ids, int *edge_ids)
    void gpu_flather_finalize(gpu_domain *GD)
    void gpu_flather_set_value(gpu_domain *GD, double stage_outside)
    void gpu_evaluate_flather_boundary(gpu_domain *GD)

    # GPU kernels
    void gpu_extrapolate_second_order(gpu_domain *GD)
    double gpu_compute_fluxes(gpu_domain *GD, int substep_count, int timestep_fluxcalls)
    void gpu_update_conserved_quantities(gpu_domain *GD, double timestep)
    void gpu_backup_conserved_quantities(gpu_domain *GD)
    void gpu_saxpy_conserved_quantities(gpu_domain *GD, double a, double b)
    void gpu_saxpy3_conserved_quantities(gpu_domain *GD, double a, double b, double c)
    double gpu_protect(gpu_domain *GD)
    double gpu_compute_water_volume(gpu_domain *GD)
    double gpu_active_set_mean_fraction(gpu_domain *GD)
    void gpu_manning_friction(gpu_domain *GD)

    # ADER-2 Cauchy-Kovalewski predictor — centroid variant
    void gpu_ader_ck_predictor(gpu_domain *GD, double dt)

    # ADER-2 fused edge predictor (shifts edges to Q^{n+1/2}, centroids unchanged)
    void gpu_ader_ck_predictor_edge(gpu_domain *GD, double dt)

    # Full ADER-2 step (single flux call with fused edge predictor, prev_dt from prior step)
    double gpu_evolve_one_ader2_step(gpu_domain *GD, double max_timestep, int apply_forcing, double prev_dt)

    # Full Euler step
    double gpu_evolve_one_euler_step(gpu_domain *GD, double max_timestep, int apply_forcing)

    # Full RK2 step
    double gpu_evolve_one_rk2_step(gpu_domain *GD, double max_timestep, int apply_forcing)

    # Full SSP-RK3 step
    double gpu_evolve_one_rk3_step(gpu_domain *GD, double max_timestep, int apply_forcing)

    void print_gpu_domain_info(gpu_domain *GD)
    int detect_gpu_aware_mpi()
    int gpu_is_available()
    int gpu_get_num_devices()
    int gpu_get_initial_device()
    int gpu_get_default_device()
    void gpu_set_default_device(int device_id)
    void gpu_set_offload_enabled(int enabled)
    int gpu_get_offload_enabled()

    # Rate operators (rain, extraction, etc.)
    int gpu_rate_operator_init(gpu_domain *GD, int num_indices, int *indices,
                               double *areas, int *full_indices, int num_full)
    void gpu_rate_operator_finalize(gpu_domain *GD, int op_id)
    void gpu_rate_operators_finalize_all(gpu_domain *GD)
    int gpu_rate_operator_tracers_capture(gpu_domain *GD, int op_id)
    void gpu_rate_operator_tracers_apply(gpu_domain *GD, int op_id,
                                         const double *c_in, const int *carry, int ns)
    double gpu_rate_operator_apply(gpu_domain *GD, int op_id,
                                   double rate, double factor, double timestep)
    double gpu_rate_operator_apply_array(gpu_domain *GD, int op_id,
                                         double *rate_array, int rate_array_size,
                                         int use_indices_into_rate,
                                         int rate_changed,
                                         double factor, double timestep)

    # Inlet operators (GPU-accelerated inlet flow)
    int gpu_inlet_operator_init(gpu_domain *GD, int num_indices, int *indices, double *areas)
    void gpu_inlet_operator_finalize(gpu_domain *GD, int op_id)
    void gpu_inlet_operators_finalize_all(gpu_domain *GD)
    double gpu_inlet_get_volume(gpu_domain *GD, int op_id)
    void gpu_inlet_get_velocities(gpu_domain *GD, int op_id, double *u_out, double *v_out)
    void gpu_inlet_set_depths(gpu_domain *GD, int op_id, double depth)
    void gpu_inlet_set_xmoms(gpu_domain *GD, int op_id, double value)
    void gpu_inlet_set_ymoms(gpu_domain *GD, int op_id, double value)
    void gpu_inlet_set_xmoms_array(gpu_domain *GD, int op_id, double *values, int n)
    void gpu_inlet_set_ymoms_array(gpu_domain *GD, int op_id, double *values, int n)
    void gpu_inlet_set_stages_evenly(gpu_domain *GD, int op_id, double volume)
    double gpu_inlet_apply(gpu_domain *GD, int op_id, double volume,
                           double current_volume, double total_area,
                           double *vel_u, double *vel_v, int num_vel,
                           int has_velocity, double ext_vel_u, double ext_vel_v,
                           int zero_velocity,
                           const double *c_in, int n_cin, double *dmass_out)

    # Culvert operators (Boyd box/pipe/weir_trapezoid - batched GPU gather/scatter)
    struct culvert_params:
        int type
        double g
        double width
        double height
        double diameter
        double z1
        double z2
        double length
        double manning
        double sum_loss
        double blockage
        double barrels
        int use_velocity_head
        int use_momentum_jet
        int use_old_momentum_method
        int always_use_Q_wetdry_adjustment
        double max_velocity
        double smoothing_timescale
        double outward_vector_0[2]
        double outward_vector_1[2]
        double invert_elevation_0
        double invert_elevation_1
        int has_invert_elevation_0
        int has_invert_elevation_1

    int gpu_culvert_init(gpu_domain *GD, const culvert_params *params,
                         int enquiry_index_0, int enquiry_index_1,
                         int inlet0_num, int *inlet0_indices, double *inlet0_areas,
                         int inlet1_num, int *inlet1_indices, double *inlet1_areas,
                         int master_proc, int enquiry_proc_0, int enquiry_proc_1,
                         int inlet_master_proc_0, int inlet_master_proc_1,
                         int is_local, int mpi_tag_base,
                         double init_smooth_Q, double init_smooth_delta_total_energy)
    void gpu_culvert_finalize(gpu_domain *GD, int culvert_id)
    void gpu_culverts_finalize_all(gpu_domain *GD)
    int gpu_culverts_map(gpu_domain *GD)
    void gpu_culverts_apply_all(gpu_domain *GD, double timestep)
    int gpu_culverts_get_report(gpu_domain *GD, int culvert_id, double *out)

    # Max-quantities operator
    int  gpu_max_quantities_init(gpu_domain *GD, int n, double velocity_zero_height)
    void gpu_max_quantities_update(gpu_domain *GD)
    void gpu_max_quantities_get(gpu_domain *GD,
                                double *out_stage, double *out_depth,
                                double *out_speed, double *out_uh)
    void gpu_max_tracers_get(gpu_domain *GD, double *out_tracer)
    void gpu_max_quantities_finalize(gpu_domain *GD)

    # FLOP counters (Gordon Bell performance profiling)
    void gpu_flop_counters_init(gpu_domain *GD)
    void gpu_flop_counters_reset(gpu_domain *GD)
    void gpu_flop_counters_enable(gpu_domain *GD, int enable)
    void gpu_flop_counters_start_timer(gpu_domain *GD)
    void gpu_flop_counters_stop_timer(gpu_domain *GD)
    uint64_t gpu_flop_counters_get_total(gpu_domain *GD)
    double gpu_flop_counters_get_flops(gpu_domain *GD)
    void gpu_flop_counters_print(gpu_domain *GD)
    uint64_t gpu_flop_counters_get_extrapolate(gpu_domain *GD)
    uint64_t gpu_flop_counters_get_compute_fluxes(gpu_domain *GD)
    uint64_t gpu_flop_counters_get_update(gpu_domain *GD)
    uint64_t gpu_flop_counters_get_protect(gpu_domain *GD)
    uint64_t gpu_flop_counters_get_manning(gpu_domain *GD)
    uint64_t gpu_flop_counters_get_backup(gpu_domain *GD)
    uint64_t gpu_flop_counters_get_saxpy(gpu_domain *GD)
    uint64_t gpu_flop_counters_get_rate_operator(gpu_domain *GD)
    uint64_t gpu_flop_counters_get_ghost_exchange(gpu_domain *GD)
    # MPI reduction for multi-GPU (Gordon Bell)
    uint64_t gpu_flop_counters_get_global_total(gpu_domain *GD)
    double gpu_flop_counters_get_global_flops(gpu_domain *GD)
    void gpu_flop_counters_print_global(gpu_domain *GD)


cdef extern from "core_kernels.h" nogil:
    void core_apply_sediment_source(domain* D, double timestep)
    void core_apply_bedload(domain* D, double timestep)
    int64_t core_apply_repose(domain* D)


cdef extern from "gpu_culvert_operator.h" nogil:
    void culvert_apply_one_host(
        int type, double g, double width, double height, double diameter,
        double z1, double z2, double length, double manning, double sum_loss,
        double blockage, double barrels,
        int use_velocity_head, int use_momentum_jet, int use_old_momentum_method,
        int always_use_Q_wetdry_adjustment, double max_velocity,
        double smoothing_timescale,
        double ov0x, double ov0y, double ov1x, double ov1y,
        double invert0, double invert1, int has_invert0, int has_invert1,
        double *smooth_delta_total_energy, double *smooth_Q,
        double timestep,
        double e0_stage, double e0_xmom, double e0_ymom, double e0_elev,
        double a0_stage, double a0_depth, double a0_xmom, double a0_ymom, double a0_area,
        double e1_stage, double e1_xmom, double e1_ymom, double e1_elev,
        double a1_stage, double a1_depth, double a1_xmom, double a1_ymom, double a1_area,
        int *inflow_idx,
        double *new_inflow_depth, double *new_inflow_xmom, double *new_inflow_ymom,
        double *new_outflow_depth, double *new_outflow_xmom, double *new_outflow_ymom,
        double *report_gain, double *report_discharge, double *report_velocity,
        double *report_driving_energy, double *report_delta_total_energy,
        double *outlet_culvert_depth)

    void culvert_gather_inlet_host(
        int n, const int *indices, const double *areas,
        const double *stage_c, const double *xmom_c,
        const double *ymom_c, const double *bed_c,
        double total_area,
        double *avg_stage, double *avg_depth,
        double *avg_xmom, double *avg_ymom)


def culvert_gather_inlet_host_py(
        int[::1] indices, double[::1] areas,
        double[::1] stage_c, double[::1] xmom_c,
        double[::1] ymom_c, double[::1] bed_c,
        double total_area):
    """Inlet gather matching the mode-2 device gather bit-for-bit.

    Returns (avg_stage, avg_depth, avg_xmom, avg_ymom).
    """
    cdef int n = indices.shape[0]
    cdef double avg_stage = 0.0, avg_depth = 0.0, avg_xmom = 0.0, avg_ymom = 0.0
    culvert_gather_inlet_host(
        n, &indices[0], &areas[0],
        &stage_c[0], &xmom_c[0], &ymom_c[0], &bed_c[0],
        total_area, &avg_stage, &avg_depth, &avg_xmom, &avg_ymom)
    return (avg_stage, avg_depth, avg_xmom, avg_ymom)


def culvert_apply_one_host_py(
        int type, double g, double width, double height, double diameter,
        double z1, double z2, double length, double manning, double sum_loss,
        double blockage, double barrels,
        int use_velocity_head, int use_momentum_jet, int use_old_momentum_method,
        int always_use_Q_wetdry_adjustment, double max_velocity,
        double smoothing_timescale,
        double ov0x, double ov0y, double ov1x, double ov1y,
        double invert0, double invert1, int has_invert0, int has_invert1,
        double smooth_dte, double smooth_Q, double timestep,
        double e0_stage, double e0_xmom, double e0_ymom, double e0_elev,
        double a0_stage, double a0_depth, double a0_xmom, double a0_ymom, double a0_area,
        double e1_stage, double e1_xmom, double e1_ymom, double e1_elev,
        double a1_stage, double a1_depth, double a1_xmom, double a1_ymom, double a1_area):
    """Bit-identical single-culvert update shared with the mode-2 batch.

    Returns (smooth_dte, smooth_Q, inflow_idx,
             new_inflow_depth, new_inflow_xmom, new_inflow_ymom,
             new_outflow_depth, new_outflow_xmom, new_outflow_ymom,
             report_gain, report_discharge, report_velocity,
             report_driving_energy, report_delta_total_energy,
             outlet_culvert_depth).
    """
    cdef double sdte = smooth_dte
    cdef double sQ = smooth_Q
    cdef int inflow_idx = 0
    cdef double nid = 0.0, nix = 0.0, niy = 0.0
    cdef double nod = 0.0, nox = 0.0, noy = 0.0
    cdef double rg = 0.0, rd = 0.0, rv = 0.0, rde = 0.0, rdte = 0.0
    cdef double ocd = 0.0
    culvert_apply_one_host(
        type, g, width, height, diameter, z1, z2, length, manning, sum_loss,
        blockage, barrels, use_velocity_head, use_momentum_jet,
        use_old_momentum_method, always_use_Q_wetdry_adjustment, max_velocity,
        smoothing_timescale, ov0x, ov0y, ov1x, ov1y,
        invert0, invert1, has_invert0, has_invert1,
        &sdte, &sQ, timestep,
        e0_stage, e0_xmom, e0_ymom, e0_elev,
        a0_stage, a0_depth, a0_xmom, a0_ymom, a0_area,
        e1_stage, e1_xmom, e1_ymom, e1_elev,
        a1_stage, a1_depth, a1_xmom, a1_ymom, a1_area,
        &inflow_idx, &nid, &nix, &niy, &nod, &nox, &noy,
        &rg, &rd, &rv, &rde, &rdte, &ocd)
    return (sdte, sQ, inflow_idx, nid, nix, niy, nod, nox, noy,
            rg, rd, rv, rde, rdte, ocd)


# ============================================================================
# GPU Domain Wrapper Class
# ============================================================================

cdef class GPUDomain:
    """
    Wrapper class holding the C gpu_domain struct.
    This provides a Python-accessible handle to the GPU domain state.
    """
    cdef gpu_domain GD
    cdef bint initialized
    cdef object python_domain      # Keep reference to prevent GC
    cdef public object _file_boundary_meta  # List of (B, vol_id, edge_id) per file-boundary edge

    def __cinit__(self):
        self.initialized = False
        self.python_domain = None
        self._file_boundary_meta = None

    def __dealloc__(self):
        if self.initialized:
            gpu_domain_finalize(&self.GD)

    @property
    def rank(self):
        return self.GD.rank

    @property
    def nprocs(self):
        return self.GD.nprocs

    @property
    def gpu_initialized(self):
        return self.GD.gpu_initialized

    @property
    def device_id(self):
        return self.GD.device_id

    @property
    def gpu_aware_mpi(self):
        return self.GD.gpu_aware_mpi

    @property
    def num_neighbors(self):
        return self.GD.halo.num_neighbors

    @property
    def recorded_flux_timestep(self):
        """CFL-constrained step of the most recent evolve_one_*_step, before the
        yieldstep/finaltime cap. Mirrors what legacy update_timestep() records
        into recorded_min/max_timestep."""
        return self.GD.recorded_flux_timestep


# ============================================================================
# MPI Communicator Extraction
# ============================================================================

cdef MPI_Comm get_mpi_comm() noexcept:
    """
    Extract the raw MPI_Comm handle from mpi4py.

    mpi4py exposes the C MPI_Comm via the .ob_mpi attribute on Comm objects.
    This is the key bridge between Python MPI and C MPI calls.

    Note: noexcept is required because MPI_Comm is an int on some platforms
    (e.g., macOS), and Cython's default error return value (NULL) is incompatible.

    In single-process builds (no mpi4py) this returns a stub communicator that is
    never used: the C kernels only touch the communicator when nprocs > 1.
    """
    IF HAVE_MPI4PY:
        import anuga.utilities.parallel_abstraction as pypar
        return (<MPI.Comm>pypar.comm).ob_mpi
    ELSE:
        return <MPI_Comm>0  # MPI_COMM_WORLD stub; unused single-process


cdef int get_mpi_rank():
    """Get current MPI rank."""
    import anuga.utilities.parallel_abstraction as pypar
    # Use function if available (works in both MPI and non-MPI modes)
    if hasattr(pypar, 'rank'):
        return pypar.rank()
    return pypar.myid


cdef int get_mpi_size():
    """Get total number of MPI processes."""
    import anuga.utilities.parallel_abstraction as pypar
    # Use function if available (works in both MPI and non-MPI modes)
    if hasattr(pypar, 'size'):
        return pypar.size()
    return pypar.numprocs


# ============================================================================
# Domain Array Pointer Extraction
# ============================================================================

cdef void get_domain_pointers(gpu_domain *GD, object domain_object):
    """
    Extract NumPy array pointers from Python domain object.

    This mirrors the pattern from sw_domain_openmp_ext.pyx but targets
    the gpu_domain struct.
    """
    cdef domain *D = &GD.D

    # Typed memoryview declarations
    cdef double[::1] stage_cv, xmom_cv, ymom_cv, bed_cv, height_cv
    cdef double[:,::1] stage_ev, xmom_ev, ymom_ev, bed_ev, height_ev
    cdef double[::1] stage_eu, xmom_eu, ymom_eu
    cdef double[::1] stage_siu, xmom_siu, ymom_siu
    cdef double[::1] stage_bv, xmom_bv, ymom_bv
    cdef double[::1] stage_backup, xmom_backup, ymom_backup
    cdef int64_t[:,::1] neighbours, neighbour_edges, surrogate_neighbours
    cdef int64_t[::1] number_of_boundaries
    cdef double[:,::1] normals, edgelengths
    cdef double[::1] areas, radii, max_speed
    cdef double[:,::1] centroid_coords, edge_coords
    cdef double[::1] x_centroid_work, y_centroid_work
    cdef int64_t[::1] tri_full_flag
    cdef double[::1] boundary_flux_sum_arr

    # Get basic parameters
    D.number_of_elements = domain_object.number_of_elements
    D.boundary_length = domain_object.boundary_length
    D.epsilon = domain_object.epsilon
    D.H0 = domain_object.H0
    D.g = domain_object.g
    D.minimum_allowed_height = domain_object.minimum_allowed_height
    D.maximum_allowed_speed = domain_object.maximum_allowed_speed
    D.evolve_max_timestep = domain_object.evolve_max_timestep

    # Copy CFL and evolve_max_timestep to GPU domain struct (used by C RK2 loop)
    GD.CFL = domain_object.CFL
    GD.evolve_max_timestep = domain_object.evolve_max_timestep
    fft = getattr(domain_object, 'fixed_flux_timestep', None)
    GD.fixed_flux_timestep = fft if fft is not None else -1.0
    GD.use_sloped_mannings = 1 if getattr(domain_object, 'use_sloped_mannings', False) else 0
    GD.use_active_set = 1 if getattr(domain_object, 'use_active_set', False) else 0
    D.low_froude = domain_object.low_froude
    # Explicit even though GPUDomain's memory is zero-initialized: bed-edge
    # reconstruction in compute_fluxes is opt-in and stays off for ANUGA.
    D.reconstruct_edge_bed = 0
    D.num_owned_edges = 0
    D.owned_edges = NULL
    D.extrapolate_velocity_second_order = domain_object.extrapolate_velocity_second_order

    # Beta values for gradient limiting
    D.beta_w = domain_object.beta_w
    D.beta_w_dry = domain_object.beta_w_dry
    D.beta_uh = domain_object.beta_uh
    D.beta_uh_dry = domain_object.beta_uh_dry
    D.beta_vh = domain_object.beta_vh
    D.beta_vh_dry = domain_object.beta_vh_dry

    # Get quantities dict
    quantities = domain_object.quantities

    # Stage
    stage = quantities["stage"]
    stage_cv = stage.centroid_values
    D.stage_centroid_values = &stage_cv[0]
    stage_ev = stage.edge_values
    D.stage_edge_values = &stage_ev[0, 0]
    stage_eu = stage.explicit_update
    D.stage_explicit_update = &stage_eu[0]
    stage_siu = stage.semi_implicit_update
    D.stage_semi_implicit_update = &stage_siu[0]
    stage_bv = stage.boundary_values
    D.stage_boundary_values = &stage_bv[0]

    # X-momentum
    xmom = quantities["xmomentum"]
    xmom_cv = xmom.centroid_values
    D.xmom_centroid_values = &xmom_cv[0]
    xmom_ev = xmom.edge_values
    D.xmom_edge_values = &xmom_ev[0, 0]
    xmom_eu = xmom.explicit_update
    D.xmom_explicit_update = &xmom_eu[0]
    xmom_siu = xmom.semi_implicit_update
    D.xmom_semi_implicit_update = &xmom_siu[0]
    xmom_bv = xmom.boundary_values
    D.xmom_boundary_values = &xmom_bv[0]

    # Y-momentum
    ymom = quantities["ymomentum"]
    ymom_cv = ymom.centroid_values
    D.ymom_centroid_values = &ymom_cv[0]
    ymom_ev = ymom.edge_values
    D.ymom_edge_values = &ymom_ev[0, 0]
    ymom_eu = ymom.explicit_update
    D.ymom_explicit_update = &ymom_eu[0]
    ymom_siu = ymom.semi_implicit_update
    D.ymom_semi_implicit_update = &ymom_siu[0]
    ymom_bv = ymom.boundary_values
    D.ymom_boundary_values = &ymom_bv[0]

    # Elevation (bed)
    elev = quantities["elevation"]
    bed_cv = elev.centroid_values
    D.bed_centroid_values = &bed_cv[0]
    bed_ev = elev.edge_values
    D.bed_edge_values = &bed_ev[0, 0]
    cdef double[::1] bed_bv
    bed_bv = elev.boundary_values
    D.bed_boundary_values = &bed_bv[0]

    # Height
    height = quantities["height"]
    height_cv = height.centroid_values
    D.height_centroid_values = &height_cv[0]
    height_ev = height.edge_values
    D.height_edge_values = &height_ev[0, 0]
    cdef double[::1] height_bv
    height_bv = height.boundary_values
    D.height_boundary_values = &height_bv[0]

    # Mesh connectivity
    neighbours = domain_object.neighbours
    D.neighbours = &neighbours[0, 0]

    neighbour_edges = domain_object.neighbour_edges
    D.neighbour_edges = &neighbour_edges[0, 0]

    surrogate_neighbours = domain_object.surrogate_neighbours
    D.surrogate_neighbours = &surrogate_neighbours[0, 0]

    number_of_boundaries = domain_object.number_of_boundaries
    D.number_of_boundaries = &number_of_boundaries[0]

    normals = domain_object.normals
    D.normals = &normals[0, 0]

    edgelengths = domain_object.edgelengths
    D.edgelengths = &edgelengths[0, 0]

    areas = domain_object.areas
    D.areas = &areas[0]

    radii = domain_object.radii
    D.radii = &radii[0]

    max_speed = domain_object.max_speed
    D.max_speed = &max_speed[0]

    # Boundary flux accumulator. core_compute_fluxes_central writes the per-substep
    # boundary flux here (host-side, after the reduction), and the Python
    # boundary_flux_integral_operator reads domain.boundary_flux_sum.
    boundary_flux_sum_arr = domain_object.boundary_flux_sum
    D.boundary_flux_sum = &boundary_flux_sum_arr[0]

    # tri_full_flag: 1 for owned (full) triangles, 0 for ghost triangles.
    # Required so compute_fluxes excludes ghosts from the local timestep minimum.
    if hasattr(domain_object, 'tri_full_flag'):
        tri_full_flag = domain_object.tri_full_flag
        D.tri_full_flag = &tri_full_flag[0]
    else:
        D.tri_full_flag = NULL

    centroid_coords = domain_object.centroid_coordinates
    D.centroid_coordinates = &centroid_coords[0, 0]

    edge_coords = domain_object.edge_coordinates
    D.edge_coordinates = &edge_coords[0, 0]

    # Work arrays for extrapolation
    x_centroid_work = domain_object.x_centroid_work
    D.x_centroid_work = &x_centroid_work[0]

    y_centroid_work = domain_object.y_centroid_work
    D.y_centroid_work = &y_centroid_work[0]

    # Friction
    cdef double[::1] friction_cv
    friction = quantities["friction"]
    friction_cv = friction.centroid_values
    D.friction_centroid_values = &friction_cv[0]

    # Backup arrays for RK2 - these are on each Quantity object
    cdef double[::1] stage_backup_cv, xmom_backup_cv, ymom_backup_cv
    stage_backup_cv = stage.centroid_backup_values
    xmom_backup_cv = xmom.centroid_backup_values
    ymom_backup_cv = ymom.centroid_backup_values
    D.stage_backup_values = &stage_backup_cv[0]
    D.xmom_backup_values = &xmom_backup_cv[0]
    D.ymom_backup_values = &ymom_backup_cv[0]

    # Riverwall arrays
    cdef int64_t[::1] edge_flux_type
    cdef int64_t[::1] edge_river_wall_counter
    cdef double[::1] riverwall_elevation
    cdef int64_t[::1] riverwall_rowIndex
    cdef double[:,::1] riverwall_hydraulic_properties

    # Get riverwallData object
    riverwallData = domain_object.riverwallData

    # Extract edge_flux_type (NULL when no river walls exist)
    if domain_object.edge_flux_type is not None:
        edge_flux_type = domain_object.edge_flux_type
        D.edge_flux_type = &edge_flux_type[0]
    else:
        D.edge_flux_type = NULL

    # Generic tracers. number_of_tracers MUST be set explicitly -- the mode-2
    # struct is PyMem_Malloc'd and not reliably zeroed, and a garbage count makes
    # the kernel guard fire on the device (CUDA_ERROR_ILLEGAL_ADDRESS).
    D.number_of_tracers = getattr(domain_object, 'number_of_tracers', 0)
    D.beta_tracer = getattr(domain_object, 'beta_tracer', 1.0)
    # Sediment: set explicitly here too -- see the note in
    # sw_domain_openmp_ext.pyx. Missing the mode-2 initialisation of
    # number_of_tracers previously produced CUDA_ERROR_ILLEGAL_ADDRESS.
    D.n_sediment_classes = getattr(domain_object, 'n_sediment_classes', 0)
    D.sediment_c_max = getattr(domain_object, 'sediment_c_max', 0.3)
    D.sediment_gamma0 = getattr(domain_object, 'sediment_gamma0', 0.0024)
    D.sediment_erosion_mode = getattr(domain_object, 'sediment_erosion_mode', 0)
    D.sediment_K_partheniades = getattr(domain_object, 'sediment_K_partheniades', 1.0e-4)
    D.sediment_dl_A = getattr(domain_object, 'sediment_dl_A', 0.000474)
    D.sediment_dl_alpha = getattr(domain_object, 'sediment_dl_alpha', 1.5)
    D.sediment_dl_beta = getattr(domain_object, 'sediment_dl_beta', 1.18)
    D.sediment_dl_threshold = getattr(domain_object, 'sediment_dl_threshold', 0.015)
    D.sediment_dl_ks = getattr(domain_object, 'sediment_dl_ks', 0.0)
    D.sediment_dl_ks_factor = getattr(domain_object, 'sediment_dl_ks_factor', 3.0)
    D.sediment_vegetation_shear = getattr(domain_object, 'sediment_vegetation_shear', 1)
    D.sediment_deposition_mode = getattr(domain_object, 'sediment_deposition_mode', 0)
    D.sediment_tau_d = getattr(domain_object, 'sediment_tau_d', 0.0)
    D.sediment_tau_crit = getattr(domain_object, 'sediment_tau_crit', 0.088)
    D.sediment_K_e = getattr(domain_object, 'sediment_K_e', 0.2e-6/0.088**0.5)
    D.sediment_rho_w = getattr(domain_object, 'sediment_rho_w', 1000.0)
    D.sediment_shear_closure = getattr(domain_object, 'sediment_shear_closure', 0)
    D.sediment_max_slope = getattr(domain_object, 'sediment_max_slope', 0.0)
    D.sediment_slope_frozen = getattr(domain_object, 'sediment_slope_frozen', 0)
    D.sediment_d_star_mode = getattr(domain_object, 'sediment_d_star_mode', 0)
    D.sediment_a_h_floor = getattr(domain_object, 'sediment_a_h_floor', 0.01)
    D.sediment_adaptation_mode = getattr(domain_object, 'sediment_adaptation_mode', 0)
    D.sediment_rouse_beta_mode = getattr(domain_object, 'sediment_rouse_beta_mode', 0)
    D.sediment_rouse_scale = getattr(domain_object, 'sediment_rouse_scale', 1.0)
    D.sediment_adaptation_alpha = getattr(domain_object, 'sediment_adaptation_alpha', 1.0)
    D.sediment_nearbed_base = getattr(domain_object, 'sediment_nearbed_base', -1)
    D.sediment_layer_fraction = getattr(domain_object, 'sediment_layer_fraction', 0.0)
    D.sediment_exchange_factor = getattr(domain_object, 'sediment_exchange_factor', 1.0)
    D.sediment_velocity_profile = getattr(domain_object, 'sediment_velocity_profile', 0)
    D.sediment_c_pack = getattr(domain_object, 'sediment_c_pack', 0.65)
    D.sediment_porosity = getattr(domain_object, 'sediment_porosity', 0.3)
    D.sediment_morphological_factor = getattr(domain_object, 'sediment_morphological_factor', 1.0)
    D.sediment_bed_evolution = 1 if getattr(domain_object, 'sediment_bed_evolution', True) else 0
    D.sediment_bedload_mode = getattr(domain_object, 'sediment_bedload_mode', 0)
    # [L-5]. Set UNCONDITIONALLY, outside the n_sediment_classes guard
    # below: the kernels test this flag before dereferencing
    # sediment_z_base, so leaving it uninitialised makes a NULL
    # pointer look like a configured base.
    D.sediment_has_z_base = getattr(domain_object, 'sediment_has_z_base', 0)
    # Bed composition: flag and thickness unconditionally, like the base;
    # the kernel also requires both arrays before it dereferences them.
    D.sediment_bed_composition = getattr(domain_object, 'sediment_bed_composition', 0)
    D.sediment_active_layer = getattr(domain_object, 'sediment_active_layer', 0.0)
    # spec 7 repose. Scalars again set unconditionally -- the kernel
    # tests sediment_repose_tan before touching sediment_repose_dz.
    D.sediment_repose_tan = getattr(domain_object, 'sediment_repose_tan', 0.0)
    D.sediment_repose_relax = getattr(domain_object, 'sediment_repose_relax', 1.0)
    D.sediment_repose_max_sweeps = getattr(domain_object, 'sediment_repose_max_sweeps', 0)
    D.sediment_bedload_K = getattr(domain_object, 'sediment_bedload_K', 3.97)
    D.sediment_bedload_m = getattr(domain_object, 'sediment_bedload_m', 1.5)
    D.sediment_bedload_tau_c_star = getattr(domain_object, 'sediment_bedload_tau_c_star', 0.0495)
    D.sediment_bedload_h_min = getattr(domain_object, 'sediment_bedload_h_min', 0.0)
    cdef double[::1] veg1
    D.vegetation_mode = getattr(domain_object, 'vegetation_mode', 0)
    D.use_sloped_mannings = 1 if getattr(domain_object, 'use_sloped_mannings', False) else 0
    D.vegetation_Cd = getattr(domain_object, 'vegetation_Cd', 1.68)
    D.vegetation_bed_chezy = getattr(domain_object, 'vegetation_bed_chezy', 65.0)
    if D.vegetation_mode > 0:
        veg1 = domain_object.quantities['veg_density'].centroid_values
        D.veg_density_centroid_values = &veg1[0]
        veg1 = domain_object.quantities['veg_diameter'].centroid_values
        D.veg_diameter_centroid_values = &veg1[0]
        veg1 = domain_object.quantities['veg_height'].centroid_values
        D.veg_height_centroid_values = &veg1[0]
    else:
        D.veg_density_centroid_values = NULL
        D.veg_diameter_centroid_values = NULL
        D.veg_height_centroid_values = NULL
    D.sediment_friction_mode = getattr(domain_object, 'sediment_friction_mode', 0)
    D.sediment_manning_ll = getattr(domain_object, 'sediment_manning_ll', 0.065)
    D.sediment_wilson_bed = getattr(domain_object, 'sediment_wilson_bed', 0)
    D.sediment_wilson_D = getattr(domain_object, 'sediment_wilson_D', 1.0e-3)
    cdef double[::1] sed1
    cdef double[:,::1] sed2
    cdef int64_t[::1] sedi
    if D.n_sediment_classes > 0:
        sed1 = domain_object.sediment_settling_velocity
        D.sediment_settling_velocity = &sed1[0]
        sed1 = domain_object.sediment_d_star
        D.sediment_d_star = &sed1[0]
        sed1 = domain_object.sediment_diameter
        D.sediment_diameter = &sed1[0]
        sed1 = domain_object.sediment_R
        D.sediment_R = &sed1[0]
        sed1 = domain_object.sediment_tau_c_star
        D.sediment_tau_c_star = &sed1[0]
        sed1 = domain_object.sediment_reference_height
        D.sediment_reference_height = &sed1[0]
        sed1 = domain_object.sediment_qbx
        D.sediment_qbx = &sed1[0]
        sed1 = domain_object.sediment_qby
        D.sediment_qby = &sed1[0]
        sed1 = domain_object.sediment_qba
        D.sediment_qba = &sed1[0]
        # [L-5]. See the note in the OpenMP binding: the scratch is
        # dereferenced whenever a class exists, the base only when set.
        sed2 = domain_object.sediment_source_limited
        D.sediment_source_limited = &sed2[0,0]
        sed1 = domain_object.sediment_slope_work
        D.sediment_slope_work = &sed1[0]
        sedi = domain_object.sediment_bed_exhausted
        D.sediment_bed_exhausted = &sedi[0]
        sed1 = domain_object.sediment_repose_dz
        D.sediment_repose_dz = &sed1[0]
        sedi = domain_object.sediment_bedload_open
        if sedi.shape[0] > 0:
            D.sediment_bedload_open = &sedi[0]
        else:
            D.sediment_bedload_open = NULL
            D.sediment_bedload_supply = NULL
        sed1 = domain_object.sediment_bedload_supply
        if sed1.shape[0] > 0:
            D.sediment_bedload_supply = &sed1[0]
        else:
            D.sediment_bedload_supply = NULL
        if getattr(domain_object, 'sediment_bed_composition', 0):
            sed2 = domain_object.sediment_bed_active
            D.sediment_bed_active = &sed2[0,0]
            sed2 = domain_object.sediment_bed_substrate
            D.sediment_bed_substrate = &sed2[0,0]
        else:
            D.sediment_bed_active = NULL
            D.sediment_bed_substrate = NULL
        if domain_object.sediment_has_z_base:
            sed1 = domain_object.sediment_z_base
            D.sediment_z_base = &sed1[0]
        else:
            D.sediment_z_base = NULL
        sfa = getattr(domain_object, 'sediment_shear_factor', None)
        if sfa is not None:
            sed1 = sfa
            D.sediment_shear_factor = &sed1[0]
        else:
            D.sediment_shear_factor = NULL
    else:
        D.sediment_settling_velocity = NULL
        D.sediment_d_star = NULL
        D.sediment_diameter = NULL
        D.sediment_R = NULL
        D.sediment_tau_c_star = NULL
        D.sediment_reference_height = NULL
        D.sediment_qbx = NULL
        D.sediment_qby = NULL
        D.sediment_qba = NULL
        D.sediment_z_base = NULL
        D.sediment_shear_factor = NULL
        D.sediment_source_limited = NULL
        D.sediment_bed_active = NULL
        D.sediment_bed_substrate = NULL
        D.sediment_slope_work = NULL
        D.sediment_bed_exhausted = NULL
        D.sediment_repose_dz = NULL
        D.sediment_bedload_open = NULL
        D.sediment_bedload_supply = NULL
    # Phase 2: wire the tracer arrays for the device. The pointers must be set
    # whenever number_of_tracers > 0 -- the shared kernels guard on that count
    # and dereference all six, so a NULL here is CUDA_ERROR_ILLEGAL_ADDRESS on
    # the device, not a degraded result. gpu_domain_core.c maps them.
    cdef double[:, ::1] tr2
    cdef double[::1] tr1
    if D.number_of_tracers > 0:
        tr2 = domain_object.tracer_centroid_values
        D.tracer_centroid_values = &tr2[0, 0]
        tr2 = domain_object.tracer_edge_values
        D.tracer_edge_values = &tr2[0, 0]
        tr2 = domain_object.tracer_boundary_values
        D.tracer_boundary_values = &tr2[0, 0]
        tr2 = domain_object.tracer_explicit_update
        D.tracer_explicit_update = &tr2[0, 0]
        tr2 = domain_object.tracer_conserved_values
        D.tracer_conserved_values = &tr2[0, 0]
        tr2 = domain_object.tracer_backup_values
        D.tracer_backup_values = &tr2[0, 0]
        if domain_object.tracer_external_source is not None:
            tr2 = domain_object.tracer_external_source
            D.tracer_external_source = &tr2[0, 0]
        else:
            D.tracer_external_source = NULL
        tr2 = domain_object.tracer_boundary_flux
        D.tracer_boundary_flux = &tr2[0, 0]
        tsf = getattr(domain_object, 'tracer_speed_factor', None)
        if tsf is not None:
            tr2 = tsf
            D.tracer_speed_factor = &tr2[0, 0]
        else:
            D.tracer_speed_factor = NULL
        tr1 = domain_object.tracer_boundary_flux_sum
        D.tracer_boundary_flux_sum = &tr1[0]
    else:
        D.tracer_centroid_values = NULL
        D.tracer_edge_values = NULL
        D.tracer_boundary_values = NULL
        D.tracer_explicit_update = NULL
        D.tracer_conserved_values = NULL
        D.tracer_backup_values = NULL
        D.tracer_external_source = NULL
        D.tracer_boundary_flux = NULL
        D.tracer_speed_factor = NULL
        D.tracer_boundary_flux_sum = NULL

    # Extract riverwall arrays (may be empty if no riverwalls)
    D.number_of_riverwall_edges = getattr(domain_object, 'number_of_riverwall_edges', 0)
    D.ncol_riverwall_hydraulic_properties = riverwallData.ncol_hydraulic_properties

    # nrow = number of unique riverwall segments = len(riverwallData.names)
    try:
        D.nrow_riverwall_hydraulic_properties = len(riverwallData.names)
    except:
        D.nrow_riverwall_hydraulic_properties = 0

    # Extract riverwall arrays with try/except for when they don't exist
    try:
        riverwall_elevation = riverwallData.riverwall_elevation
        D.riverwall_elevation = &riverwall_elevation[0]
    except:
        D.riverwall_elevation = NULL

    try:
        riverwall_rowIndex = riverwallData.hydraulic_properties_rowIndex
        D.riverwall_rowIndex = &riverwall_rowIndex[0]
    except:
        D.riverwall_rowIndex = NULL

    try:
        riverwall_hydraulic_properties = riverwallData.hydraulic_properties
        D.riverwall_hydraulic_properties = &riverwall_hydraulic_properties[0, 0]
    except:
        D.riverwall_hydraulic_properties = NULL

    # edge_river_wall_counter is on domain_object, not riverwallData (NULL when no river walls)
    if domain_object.edge_river_wall_counter is not None:
        edge_river_wall_counter = domain_object.edge_river_wall_counter
        D.edge_river_wall_counter = &edge_river_wall_counter[0]
    else:
        D.edge_river_wall_counter = NULL


# ============================================================================
# Halo Structure Building
# ============================================================================

cdef void build_halo_from_dicts(gpu_domain *GD, object domain_object):
    """
    Convert Python ghost dictionaries to flattened C arrays for GPU.

    ANUGA stores ghost exchange info in:
    - domain.full_send_dict[proc] = [indices, global_ids, buffer]
    - domain.ghost_recv_dict[proc] = [indices, global_ids, buffer]

    We flatten these to contiguous arrays for efficient GPU access.
    """
    full_send_dict = domain_object.full_send_dict
    ghost_recv_dict = domain_object.ghost_recv_dict

    # Count neighbors (excluding self)
    my_rank = GD.rank
    neighbors = []
    for proc in full_send_dict:
        if proc != my_rank:
            neighbors.append(proc)

    cdef int num_neighbors = len(neighbors)
    if num_neighbors == 0:
        return

    # Allocate temporary arrays
    cdef np.ndarray[int, ndim=1] neighbor_ranks = np.array(neighbors, dtype=np.int32)
    cdef np.ndarray[int, ndim=1] send_counts = np.zeros(num_neighbors, dtype=np.int32)
    cdef np.ndarray[int, ndim=1] recv_counts = np.zeros(num_neighbors, dtype=np.int32)

    # Count send/recv sizes
    cdef int total_send = 0
    cdef int total_recv = 0
    cdef int ni
    for ni in range(num_neighbors):
        proc = neighbor_ranks[ni]
        send_counts[ni] = len(full_send_dict[proc][0])
        recv_counts[ni] = len(ghost_recv_dict[proc][0])
        total_send += send_counts[ni]
        total_recv += recv_counts[ni]

    # Flatten indices
    cdef np.ndarray[int, ndim=1] flat_send = np.zeros(total_send, dtype=np.int32)
    cdef np.ndarray[int, ndim=1] flat_recv = np.zeros(total_recv, dtype=np.int32)

    cdef int send_offset = 0
    cdef int recv_offset = 0
    cdef int j, idx

    for ni in range(num_neighbors):
        proc = neighbor_ranks[ni]

        # Copy send indices
        send_indices = full_send_dict[proc][0]
        for j in range(send_counts[ni]):
            flat_send[send_offset + j] = send_indices[j]
        send_offset += send_counts[ni]

        # Copy recv indices
        recv_indices = ghost_recv_dict[proc][0]
        for j in range(recv_counts[ni]):
            flat_recv[recv_offset + j] = recv_indices[j]
        recv_offset += recv_counts[ni]

    # Call C function to initialize halo
    gpu_halo_init(GD, num_neighbors,
                  <int*>neighbor_ranks.data,
                  <int*>send_counts.data,
                  <int*>recv_counts.data,
                  <int*>flat_send.data,
                  <int*>flat_recv.data)


# ============================================================================
# Public Python API
# ============================================================================

def get_last_error(GPUDomain gpu_dom):
    """The last error message recorded by the C layer (gpu_set_error).

    Every C init function that returns an error code records the reason
    here first; the wrappers include it in the exception they raise.
    """
    return gpu_dom.GD.last_error.decode('utf-8', 'replace')


def init_gpu_domain(object domain_object, bint verbose=True):
    """
    Initialize GPU domain from Python domain object.

    This is the main entry point called from Python. It:
    1. Extracts MPI_Comm from mpi4py
    2. Creates gpu_domain struct with array pointers
    3. Builds flattened halo exchange structures
    4. Returns a GPUDomain wrapper object

    Parameters
    ----------
    domain_object : anuga.shallow_water.Domain
        The Python domain object (must be already partitioned for MPI)
    verbose : bool, optional
        If True, print GPU initialisation and array-mapping messages.
        Default is False (silent) so pytest runs are not noisy.

    Returns
    -------
    GPUDomain
        Wrapper object holding the GPU domain state
    """
    cdef MPI_Comm comm = get_mpi_comm()
    cdef int rank = get_mpi_rank()
    cdef int nprocs = get_mpi_size()

    # Create wrapper object
    gpu_dom = GPUDomain()

    # Initialize C struct with MPI info
    gpu_domain_init(&gpu_dom.GD, comm, rank, nprocs)

    # Propagate verbose flag to C struct so printf calls respect it
    gpu_dom.GD.verbose = 1 if verbose else 0

    # Ensure lazy work arrays are allocated before the C struct is built.
    # Needed in CPU_ONLY_MODE where GPU kernels call the openmp C extension.
    # edge_flux_type/edge_river_wall_counter are handled separately (may be NULL).
    if hasattr(domain_object, '_ensure_work_arrays'):
        domain_object._ensure_work_arrays()

    # Extract array pointers from Python domain
    get_domain_pointers(&gpu_dom.GD, domain_object)

    # Build halo exchange structures from ghost dicts
    if hasattr(domain_object, 'full_send_dict') and domain_object.full_send_dict:
        build_halo_from_dicts(&gpu_dom.GD, domain_object)

    # Keep reference to Python domain to prevent GC of arrays
    gpu_dom.python_domain = domain_object
    gpu_dom.initialized = True

    if verbose and rank == 0:
        import sys
        sys.stdout.flush()
        print_gpu_domain_info(&gpu_dom.GD)

    return gpu_dom


def estimate_required_memory(int64_t n, int64_t nb):
    """
    Estimate device memory (bytes) needed for a domain with n triangles and nb boundary edges.

    Parameters
    ----------
    n : int
        Number of triangles (number_of_elements).
    nb : int
        Number of boundary edges (boundary_length).

    Returns
    -------
    int
        Estimated bytes that gpu_domain_map_arrays will transfer to the device.
    """
    return gpu_estimate_required_memory(n, nb)


def gpu_available():
    """
    Return True if a real GPU offload target is active for this build.

    Returns False in CPU_ONLY_MODE (standard conda/pip install, CI) and on
    machines where ``omp_get_num_devices()`` reports no devices even though a
    GPU build was compiled.

    Example::

        from anuga.shallow_water.sw_domain_gpu_ext import gpu_available
        if gpu_available():
            domain.set_multiprocessor_mode(2)  # real GPU
        else:
            domain.set_multiprocessor_mode(1)  # CPU OpenMP
    """
    return bool(gpu_is_available())


def get_num_gpu_devices():
    """Return the number of OpenMP offload devices (GPUs) available.

    Returns 0 in CPU_ONLY_MODE or when no GPU offload target is present.
    Use this alongside ``gpu_available()`` to decide whether there are
    enough GPUs for a given number of MPI ranks::

        from anuga.shallow_water.sw_domain_gpu_ext import gpu_available, get_num_gpu_devices
        n = get_num_gpu_devices()
        # safe to run 4-rank GPU test only if 4 GPUs present
        if gpu_available() and n < 4:
            pytest.skip(f'need 4 GPUs, have {n}')
    """
    return gpu_get_num_devices()


def get_initial_device():
    """Return the OpenMP initial (host) device number.

    Pass this to :func:`set_default_device` to route subsequent ``omp target``
    regions back to the host (CPU). Returns 0 in CPU_ONLY_MODE.
    """
    return gpu_get_initial_device()


def get_default_device():
    """Return the current OpenMP default offload device number.

    The unified kernels run on the default device (no ``device()`` clause), so
    this is the device they will use. Returns 0 in CPU_ONLY_MODE.
    """
    return gpu_get_default_device()


def set_default_device(int device_id):
    """Route subsequent ``omp target`` regions to ``device_id`` (process-wide).

    The unified kernels carry no ``device()`` clause, so changing the default
    device is a robust, runtime way to move work between a GPU and the host —
    unlike ``OMP_TARGET_OFFLOAD``, which the runtime only reads at init. Pass
    :func:`get_initial_device` to force the host (CPU). No-op in CPU_ONLY_MODE.
    """
    gpu_set_default_device(device_id)


def set_offload_enabled(enabled):
    """Set the process-global offload flag honoured at GPU-domain init.

    When False, domains built afterwards stay on the host (CPU) even on a GPU
    build, keeping data mapping and kernel execution consistent. Decide before
    building the first 'unified' domain — the device is chosen when arrays are
    mapped. Use :func:`anuga.set_gpu_offload` rather than calling this directly.
    """
    gpu_set_offload_enabled(1 if enabled else 0)


def get_offload_enabled():
    """Return the process-global offload flag (see :func:`set_offload_enabled`)."""
    return bool(gpu_get_offload_enabled())


def gpu_has_mpi():
    """Return True if this extension was compiled with real C MPI support.

    Returns False when built against single-process stubs (mpi.h was not
    found at build time, HAVE_MPI4PY=False).

    Multi-rank GPU halo-exchange tests must be skipped when this returns
    False: the C-level exchange_ghosts() would be a no-op, producing
    incorrect results in parallel runs.
    """
    IF HAVE_MPI4PY:
        return True
    ELSE:
        return False


def is_gpu_aware_mpi():
    """
    Return True if the current MPI library supports GPU-aware communication.

    GPU-aware MPI allows ``MPI_Isend``/``MPI_Irecv`` to operate directly on
    device (GPU) buffers, eliminating the host-side copy in halo exchange.

    Detection priority
    ------------------
    1. Compile-time ``-DGPU_AWARE_MPI`` flag (user explicitly enabled).
    2. Runtime ``MPIX_Query_cuda_support()`` — Open MPI / MVAPICH2 CUDA builds
       (requires ``mpi-ext.h`` at build time).
    3. Runtime ``MPIX_Query_rocm_support()`` — Open MPI / MVAPICH2 ROCm builds.
    4. Returns False if none of the above is available.

    Notes
    -----
    - A True result means the MPI library *claims* GPU-aware support.  The
      buffers must still be allocated with ``omp_target_alloc`` (done
      automatically when the build was compiled with ``-DGPU_AWARE_MPI``).
    - On a standard conda install (CPU_ONLY_MODE) this will return False.
    - To force-enable, rebuild with ``-Dgpu_aware_mpi=true -Dgpu_offload=true``.

    Example::

        from anuga.shallow_water.sw_domain_gpu_ext import is_gpu_aware_mpi
        print("GPU-aware MPI:", is_gpu_aware_mpi())
    """
    return bool(detect_gpu_aware_mpi())


def check_gpu_device_memory(GPUDomain gpu_dom):
    """
    Run the C-level device memory check for this domain.

    Returns 1 if the device has enough memory (or if the device cannot be
    queried, e.g., CPU_ONLY_MODE), 0 if memory is insufficient.
    Also prints a diagnostic line when memory is insufficient.
    """
    return gpu_check_device_memory(&gpu_dom.GD)


def check_device_memory(GPUDomain gpu_dom):
    """
    Check that the target device has enough free memory for this domain.

    Returns (ok, required_mb, free_mb, total_mb).  free_mb and total_mb are 0
    when the device cannot be queried (CPU_ONLY_MODE or unknown vendor).

    Raises RuntimeError with an actionable message if memory is insufficient.
    """
    cdef int64_t n  = gpu_dom.GD.D.number_of_elements
    cdef int64_t nb = gpu_dom.GD.D.boundary_length
    cdef size_t required = gpu_estimate_required_memory(n, nb)
    cdef size_t free_b = 0, total_b = 0

    # Re-use the C query; ignore return value (check_device_memory prints errors)
    ok = gpu_check_device_memory(&gpu_dom.GD)
    if not ok:
        req_mb = required / (1024.0 * 1024.0)
        raise RuntimeError(
            f"Insufficient GPU memory for domain with {n} triangles "
            f"(estimated {req_mb:.0f} MB required). "
            "Use fewer triangles, more MPI ranks, or set_multiprocessor_mode(1)."
        )


def map_to_gpu(GPUDomain gpu_dom):
    """
    Map domain arrays to GPU memory.

    Call this once after init_gpu_domain, before starting the evolve loop.
    Raises RuntimeError if the device does not have enough free memory.
    """
    cdef int ok
    cdef int64_t n = gpu_dom.GD.D.number_of_elements
    cdef int64_t nb = gpu_dom.GD.D.boundary_length
    ok = gpu_domain_map_arrays(&gpu_dom.GD)
    if not ok:
        req_mb = gpu_estimate_required_memory(n, nb) / (1024.0 * 1024.0)
        raise RuntimeError(
            f"GPU memory check failed for domain with {n} triangles "
            f"(estimated {req_mb:.0f} MB required). "
            "Use fewer triangles, more MPI ranks, or set_multiprocessor_mode(1)."
        )


def remap_boundary_arrays(GPUDomain gpu_dom):
    """
    Remap boundary arrays to GPU memory.

    Call this after boundaries are initialized (if they weren't initialized
    before map_to_gpu was called).
    """
    gpu_remap_boundary_arrays(&gpu_dom.GD)


def unmap_from_gpu(GPUDomain gpu_dom):
    """
    Unmap domain arrays from GPU memory.

    Call this when done with GPU computation.
    """
    gpu_domain_unmap_arrays(&gpu_dom.GD)


def sync_to_device(GPUDomain gpu_dom):
    """
    Sync centroid values from host to device.

    Call this after Python modifies arrays (e.g., boundary updates).
    """
    gpu_domain_sync_to_device(&gpu_dom.GD)


def sync_from_device(GPUDomain gpu_dom):
    """
    Sync centroid values from device to host.

    Call this before Python needs to read arrays (e.g., at yieldstep for I/O).
    """
    gpu_domain_sync_from_device(&gpu_dom.GD)


def sync_all_from_device(GPUDomain gpu_dom):
    """
    Sync ALL arrays from device to host (for debugging/testing).

    This syncs centroid_values, edge_values, boundary_values,
    explicit_update, semi_implicit_update, and backup_values.
    Use this when you need to inspect intermediate GPU values.
    """
    gpu_domain_sync_all_from_device(&gpu_dom.GD)


def sync_boundary_values(GPUDomain gpu_dom):
    """
    Sync boundary values from host to device.

    Call this after Python updates boundary conditions.
    """
    gpu_sync_boundary_values(&gpu_dom.GD)


def sync_riverwall_to_device(GPUDomain gpu_dom):
    """
    Sync riverwall crest elevations and hydraulic properties from host to device.

    Call this after Python changes a riverwall at runtime (RiverWall.set_elevation(),
    set_elevation_offset(), set_hydraulic_parameter()). No-op when the domain has no
    riverwalls.
    """
    gpu_sync_riverwall_to_device(&gpu_dom.GD)


def sync_tracer_source_to_device(GPUDomain gpu_dom):
    """Push the tracer external source to the device.

    set_tracer_source writes the host array; a time-varying source rewrites it
    every step, so the device copy has to be refreshed or it goes stale. One
    array, not the whole state. #288
    """
    gpu_sync_tracer_source_to_device(&gpu_dom.GD)


def sync_tracer_boundary_to_device(GPUDomain gpu_dom):
    """Push the tracer boundary concentrations to the device.

    set_tracer_boundary writes the host array, and a callable boundary is
    re-evaluated every step; the per-step boundary push carries only the
    hydrodynamic values, so this one has its own. One small array.
    """
    gpu_sync_tracer_boundary_to_device(&gpu_dom.GD)


def sync_edge_values_from_device(GPUDomain gpu_dom):
    """
    Sync ALL edge values from device to host.

    Call this before CPU boundary evaluation needs to read edge values.
    WARNING: This is expensive - use sync_boundary_edge_values for sparse sync.
    """
    gpu_sync_edge_values_from_device(&gpu_dom.GD)


def init_boundary_edge_sync(GPUDomain gpu_dom, np.ndarray[int, ndim=1, mode="c"] boundary_cell_ids):
    """
    Initialize boundary edge sync buffers - call once after boundaries are set.

    This pre-allocates staging buffers on GPU for efficient sparse sync.

    Parameters
    ----------
    gpu_dom : GPUDomain
        The GPU domain wrapper
    boundary_cell_ids : ndarray
        Unique cell IDs that are adjacent to boundaries
    """
    cdef int num_cells = len(boundary_cell_ids)
    cdef int result = gpu_boundary_edge_sync_init(&gpu_dom.GD, num_cells, &boundary_cell_ids[0])
    if result != 0:
        raise RuntimeError("Failed to initialize boundary edge sync")


def boundary_edge_sync(GPUDomain gpu_dom):
    """
    Sync edge values for boundary-adjacent cells from GPU to host.

    Fast operation - uses pre-allocated buffers. Call init_boundary_edge_sync first.
    Call this before CPU boundary evaluation needs to read edge values.
    """
    gpu_boundary_edge_sync(&gpu_dom.GD)


def exchange_ghosts(GPUDomain gpu_dom):
    """
    Exchange ghost cell data between MPI ranks.

    This uses GPU-aware MPI if available, otherwise does efficient
    D2H/H2D transfers of only the small halo buffers.
    """
    gpu_exchange_ghosts(&gpu_dom.GD)



def exchange_bed_ghosts_gpu(GPUDomain gpu_dom):
    """Give every ghost cell its owner's bed (centroid, edges shifted with it).

    Called after the sediment step in parallel runs whose bed evolves (#424).
    """
    gpu_exchange_bed_ghosts(&gpu_dom.GD)

def init_reflective_boundary(GPUDomain gpu_dom, object domain_object):
    """
    Initialize reflective boundary for GPU evaluation.

    Extracts reflective boundary info from domain and maps to GPU.
    Call this once after domain setup.

    Parameters
    ----------
    gpu_dom : GPUDomain
        The GPU domain wrapper
    domain_object : Domain
        The Python domain object
    """
    cdef int num_edges = 0
    cdef np.ndarray[int, ndim=1, mode="c"] boundary_indices
    cdef np.ndarray[int, ndim=1, mode="c"] vol_ids_arr
    cdef np.ndarray[int, ndim=1, mode="c"] edge_ids_arr

    # Collect all reflective boundary edges (may have multiple tags)
    # boundary_map may not be set up yet if called early
    if domain_object.boundary_map is None:
        return

    all_ids = []
    for tag, boundary in domain_object.boundary_map.items():
        if boundary is not None and boundary.__class__.__name__ == 'Reflective_boundary':
            segment_edges = domain_object.tag_boundary_cells.get(tag, None)
            if segment_edges is not None and len(segment_edges) > 0:
                all_ids.extend(segment_edges)

    if len(all_ids) == 0:
        # No reflective boundary, nothing to do
        return

    # Build arrays
    ids = np.array(all_ids, dtype=np.intc)
    num_edges = len(ids)
    boundary_indices = ids
    vol_ids_arr = np.ascontiguousarray(domain_object.boundary_cells[ids], dtype=np.intc)
    edge_ids_arr = np.ascontiguousarray(domain_object.boundary_edges[ids], dtype=np.intc)

    # Call C init
    if     gpu_reflective_init(&gpu_dom.GD, num_edges,
                           &boundary_indices[0], &vol_ids_arr[0], &edge_ids_arr[0]) != 0:
        raise RuntimeError('init_reflective_boundary: ' + get_last_error(gpu_dom))


def evaluate_reflective_boundary_gpu(GPUDomain gpu_dom):
    """
    Evaluate reflective boundary on GPU.

    This reads edge values and writes boundary values entirely on device.
    No data transfer required.
    """
    gpu_evaluate_reflective_boundary(&gpu_dom.GD)


def init_dirichlet_boundary(GPUDomain gpu_dom, object domain_object):
    """
    Initialize Dirichlet boundary for GPU evaluation.

    Extracts Dirichlet boundary info from domain and prepares for GPU mapping.
    Call this once after domain setup, BEFORE map_to_gpu.
    """
    cdef int num_edges = 0
    cdef np.ndarray[int, ndim=1, mode="c"] boundary_indices
    cdef np.ndarray[int, ndim=1, mode="c"] vol_ids_arr
    cdef np.ndarray[int, ndim=1, mode="c"] edge_ids_arr
    cdef np.ndarray[double, ndim=1, mode="c"] stage_values_arr
    cdef np.ndarray[double, ndim=1, mode="c"] xmom_values_arr
    cdef np.ndarray[double, ndim=1, mode="c"] ymom_values_arr

    if domain_object.boundary_map is None:
        return

    # Collect every Dirichlet boundary edge with that boundary's own values, so
    # multiple Dirichlet boundaries with different values are handled per-edge.
    all_ids = []
    all_stage = []
    all_xmom = []
    all_ymom = []
    for tag, boundary in domain_object.boundary_map.items():
        if boundary is not None and boundary.__class__.__name__ == 'Dirichlet_boundary':
            segment_edges = domain_object.tag_boundary_cells.get(tag, None)
            if segment_edges is not None and len(segment_edges) > 0:
                if hasattr(boundary, 'dirichlet_values') and len(boundary.dirichlet_values) >= 3:
                    s = float(boundary.dirichlet_values[0])
                    x = float(boundary.dirichlet_values[1])
                    y = float(boundary.dirichlet_values[2])
                else:
                    s = x = y = 0.0
                for e in segment_edges:
                    all_ids.append(e)
                    all_stage.append(s)
                    all_xmom.append(x)
                    all_ymom.append(y)

    if len(all_ids) == 0:
        return

    ids = np.array(all_ids, dtype=np.intc)
    num_edges = len(ids)
    boundary_indices = ids
    vol_ids_arr = np.ascontiguousarray(domain_object.boundary_cells[ids], dtype=np.intc)
    edge_ids_arr = np.ascontiguousarray(domain_object.boundary_edges[ids], dtype=np.intc)
    stage_values_arr = np.ascontiguousarray(all_stage, dtype=np.float64)
    xmom_values_arr = np.ascontiguousarray(all_xmom, dtype=np.float64)
    ymom_values_arr = np.ascontiguousarray(all_ymom, dtype=np.float64)

    if     gpu_dirichlet_init(&gpu_dom.GD, num_edges,
                          &boundary_indices[0], &vol_ids_arr[0], &edge_ids_arr[0],
                          &stage_values_arr[0], &xmom_values_arr[0], &ymom_values_arr[0]) != 0:
        raise RuntimeError('init_dirichlet_boundary: ' + get_last_error(gpu_dom))


def evaluate_dirichlet_boundary_gpu(GPUDomain gpu_dom):
    """
    Evaluate Dirichlet boundary on GPU.

    Sets constant values at boundary, entirely on device.
    No data transfer required.
    """
    gpu_evaluate_dirichlet_boundary(&gpu_dom.GD)


def init_transmissive_boundary(GPUDomain gpu_dom, object domain_object):
    """
    Initialize Transmissive boundary for GPU evaluation.

    Extracts Transmissive boundary info from domain and prepares for GPU mapping.
    Call this once after domain setup, BEFORE map_to_gpu.
    """
    cdef int num_edges = 0
    cdef np.ndarray[int, ndim=1, mode="c"] boundary_indices
    cdef np.ndarray[int, ndim=1, mode="c"] vol_ids_arr
    cdef np.ndarray[int, ndim=1, mode="c"] edge_ids_arr
    cdef int use_centroid = 0

    if domain_object.boundary_map is None:
        return

    all_ids = []
    for tag, boundary in domain_object.boundary_map.items():
        if boundary is not None and boundary.__class__.__name__ == 'Transmissive_boundary':
            segment_edges = domain_object.tag_boundary_cells.get(tag, None)
            if segment_edges is not None and len(segment_edges) > 0:
                all_ids.extend(segment_edges)

    if len(all_ids) == 0:
        return

    # Check if domain uses centroid transmissive BC
    if hasattr(domain_object, 'centroid_transmissive_bc'):
        use_centroid = 1 if domain_object.centroid_transmissive_bc else 0

    ids = np.array(all_ids, dtype=np.intc)
    num_edges = len(ids)
    boundary_indices = ids
    vol_ids_arr = np.ascontiguousarray(domain_object.boundary_cells[ids], dtype=np.intc)
    edge_ids_arr = np.ascontiguousarray(domain_object.boundary_edges[ids], dtype=np.intc)

    if     gpu_transmissive_init(&gpu_dom.GD, num_edges,
                             &boundary_indices[0], &vol_ids_arr[0], &edge_ids_arr[0],
                             use_centroid) != 0:
        raise RuntimeError('init_transmissive_boundary: ' + get_last_error(gpu_dom))


def evaluate_transmissive_boundary_gpu(GPUDomain gpu_dom):
    """
    Evaluate Transmissive boundary on GPU.

    Copies interior values to boundary, entirely on device.
    No data transfer required.
    """
    gpu_evaluate_transmissive_boundary(&gpu_dom.GD)


def init_transmissive_n_zero_t_boundary(GPUDomain gpu_dom, object domain_object):
    """
    Initialize Transmissive_n_momentum_zero_t_momentum_set_stage boundary for GPU.

    This boundary type sets stage from a (potentially time-varying) function,
    keeps normal momentum, and zeros tangential momentum.

    Call this once after domain setup, BEFORE map_to_gpu.
    """
    cdef int num_edges = 0
    cdef np.ndarray[int, ndim=1, mode="c"] boundary_indices
    cdef np.ndarray[int, ndim=1, mode="c"] vol_ids_arr
    cdef np.ndarray[int, ndim=1, mode="c"] edge_ids_arr

    if domain_object.boundary_map is None:
        return

    all_ids = []
    for tag, boundary in domain_object.boundary_map.items():
        if boundary is not None and boundary.__class__.__name__ == 'Transmissive_n_momentum_zero_t_momentum_set_stage_boundary':
            segment_edges = domain_object.tag_boundary_cells.get(tag, None)
            if segment_edges is not None and len(segment_edges) > 0:
                all_ids.extend(segment_edges)

    if len(all_ids) == 0:
        return

    ids = np.array(all_ids, dtype=np.intc)
    num_edges = len(ids)
    boundary_indices = ids
    vol_ids_arr = np.ascontiguousarray(domain_object.boundary_cells[ids], dtype=np.intc)
    edge_ids_arr = np.ascontiguousarray(domain_object.boundary_edges[ids], dtype=np.intc)

    if     gpu_transmissive_n_zero_t_init(&gpu_dom.GD, num_edges,
                                      &boundary_indices[0], &vol_ids_arr[0], &edge_ids_arr[0]) != 0:
        raise RuntimeError('init_transmissive_n_zero_t_boundary: ' + get_last_error(gpu_dom))


def set_transmissive_n_zero_t_stage(GPUDomain gpu_dom, double stage_value):
    """
    Update the stage value for Transmissive_n_zero_t boundary.

    Call this each timestep before evaluate_transmissive_n_zero_t_boundary_gpu.
    """
    gpu_transmissive_n_zero_t_set_stage(&gpu_dom.GD, stage_value)


def evaluate_transmissive_n_zero_t_boundary_gpu(GPUDomain gpu_dom):
    """
    Evaluate Transmissive_n_zero_t boundary on GPU.

    Sets stage from external value, keeps normal momentum, zeros tangential.
    Call set_transmissive_n_zero_t_stage first to set the stage value.
    """
    gpu_evaluate_transmissive_n_zero_t_boundary(&gpu_dom.GD)


def init_file_boundary(GPUDomain gpu_dom, object domain_object):
    """
    Initialize File_boundary / Field_boundary for GPU.

    Scans domain.boundary_map for File_boundary and Field_boundary objects,
    collects their edges, and initialises the GPU file_boundary struct.

    Also stores per-edge evaluation metadata on gpu_dom so that
    set_file_boundary_values_from_domain() can update values each timestep.

    Call this once after domain setup, BEFORE map_to_gpu.
    """
    cdef int num_edges = 0
    cdef np.ndarray[int, ndim=1, mode="c"] boundary_indices
    cdef np.ndarray[int, ndim=1, mode="c"] vol_ids_arr
    cdef np.ndarray[int, ndim=1, mode="c"] edge_ids_arr

    FILE_BOUNDARY_TYPES = {'File_boundary', 'Field_boundary'}

    if domain_object.boundary_map is None:
        return

    # Collect (boundary_index, vol_id, edge_id, boundary_object) for all edges
    # belonging to File_boundary or Field_boundary tags.
    all_ids = []
    edge_meta = []  # list of (B, vol_id, edge_id, boundary_index)
    for tag, boundary in domain_object.boundary_map.items():
        if boundary is not None and boundary.__class__.__name__ in FILE_BOUNDARY_TYPES:
            segment_edges = domain_object.tag_boundary_cells.get(tag, None)
            if segment_edges is not None and len(segment_edges) > 0:
                for bid in segment_edges:
                    vol_id  = int(domain_object.boundary_cells[bid])
                    edge_id = int(domain_object.boundary_edges[bid])
                    all_ids.append(bid)
                    edge_meta.append((boundary, vol_id, edge_id, int(bid)))

    if len(all_ids) == 0:
        return

    ids = np.array(all_ids, dtype=np.intc)
    num_edges = len(ids)
    boundary_indices = ids
    vol_ids_arr  = np.ascontiguousarray(domain_object.boundary_cells[ids],  dtype=np.intc)
    edge_ids_arr = np.ascontiguousarray(domain_object.boundary_edges[ids], dtype=np.intc)

    if     gpu_file_boundary_init(&gpu_dom.GD, num_edges,
                              &boundary_indices[0], &vol_ids_arr[0], &edge_ids_arr[0]) != 0:
        raise RuntimeError('init_file_boundary_boundary: ' + get_last_error(gpu_dom))

    # Store edge metadata for per-timestep Python evaluation
    # List of (B, vol_id, edge_id) in the same order as boundary_indices
    gpu_dom._file_boundary_meta = [(B, vid, eid) for (B, vid, eid, bid) in edge_meta]


def set_file_boundary_values_from_domain(GPUDomain gpu_dom, object domain_object):
    """
    Evaluate File_boundary / Field_boundary for the current simulation time
    and push per-edge values to the GPU.

    Call this each timestep before evaluate_file_boundary_gpu().
    """
    meta = getattr(gpu_dom, '_file_boundary_meta', None)
    if meta is None or len(meta) == 0:
        return

    cdef int ne = len(meta)
    cdef np.ndarray[double, ndim=1, mode="c"] stage_arr = np.empty(ne, dtype=np.float64)
    cdef np.ndarray[double, ndim=1, mode="c"] xmom_arr  = np.empty(ne, dtype=np.float64)
    cdef np.ndarray[double, ndim=1, mode="c"] ymom_arr  = np.empty(ne, dtype=np.float64)

    for k, (B, vol_id, edge_id) in enumerate(meta):
        q = B.evaluate(vol_id, edge_id)
        stage_arr[k] = q[0]
        xmom_arr[k]  = q[1]
        ymom_arr[k]  = q[2]

    gpu_file_boundary_set_values(&gpu_dom.GD,
                                 &stage_arr[0], &xmom_arr[0], &ymom_arr[0])


def evaluate_file_boundary_gpu(GPUDomain gpu_dom):
    """
    Evaluate File_boundary / Field_boundary on GPU.

    Writes previously-pushed per-edge values into the global boundary_values
    arrays and fills bed/height from adjacent interior edges.
    Call set_file_boundary_values_from_domain() first.
    """
    gpu_evaluate_file_boundary(&gpu_dom.GD)


def init_time_boundary(GPUDomain gpu_dom, object domain_object):
    """
    Initialize Time_boundary for GPU.

    This boundary type sets conserved quantities from a time-dependent function.
    The function is called from Python each timestep (cheap), and the assignment
    happens on GPU (avoids D2H/H2D transfer overhead).

    Call this once after domain setup, BEFORE map_to_gpu.
    """
    cdef int num_edges = 0
    cdef np.ndarray[int, ndim=1, mode="c"] boundary_indices
    cdef np.ndarray[int, ndim=1, mode="c"] vol_ids_arr
    cdef np.ndarray[int, ndim=1, mode="c"] edge_ids_arr

    if domain_object.boundary_map is None:
        return

    all_ids = []
    for tag, boundary in domain_object.boundary_map.items():
        if boundary is not None and boundary.__class__.__name__ == 'Time_boundary':
            segment_edges = domain_object.tag_boundary_cells.get(tag, None)
            if segment_edges is not None and len(segment_edges) > 0:
                all_ids.extend(segment_edges)

    if len(all_ids) == 0:
        return

    ids = np.array(all_ids, dtype=np.intc)
    num_edges = len(ids)
    boundary_indices = ids
    vol_ids_arr = np.ascontiguousarray(domain_object.boundary_cells[ids], dtype=np.intc)
    edge_ids_arr = np.ascontiguousarray(domain_object.boundary_edges[ids], dtype=np.intc)

    if     gpu_time_boundary_init(&gpu_dom.GD, num_edges,
                              &boundary_indices[0], &vol_ids_arr[0], &edge_ids_arr[0]) != 0:
        raise RuntimeError('init_time_boundary_boundary: ' + get_last_error(gpu_dom))


def set_time_boundary_values(GPUDomain gpu_dom,
                             np.ndarray[double, ndim=1, mode="c"] stage,
                             np.ndarray[double, ndim=1, mode="c"] xmom,
                             np.ndarray[double, ndim=1, mode="c"] ymom):
    """
    Update the per-edge values for Time_boundary.

    Call this each timestep before evaluate_time_boundary_gpu. The arrays hold
    one value per time-boundary edge, concatenated across all Time_boundary tags
    in boundary_map order (matching init_time_boundary), so multiple Time_boundary
    objects with different values do not clobber one another.
    """
    if len(stage) == 0:
        return
    gpu_time_boundary_set_values(&gpu_dom.GD, &stage[0], &xmom[0], &ymom[0])


def evaluate_time_boundary_gpu(GPUDomain gpu_dom):
    """
    Evaluate Time_boundary on GPU.

    Sets stage and momentum from time-dependent values.
    Call set_time_boundary_values first to set the values.
    """
    gpu_evaluate_time_boundary(&gpu_dom.GD)


def init_absorbing_wave_boundary(GPUDomain gpu_dom, object domain_object):
    """
    Initialize Absorbing_wave_boundary for GPU.

    Scans domain.boundary_map for Absorbing_wave_boundary objects, collects
    their edges, and initialises the GPU struct.  Call once after domain setup,
    BEFORE map_to_gpu.
    """
    cdef int num_edges = 0
    cdef np.ndarray[int, ndim=1, mode="c"] boundary_indices
    cdef np.ndarray[int, ndim=1, mode="c"] vol_ids_arr
    cdef np.ndarray[int, ndim=1, mode="c"] edge_ids_arr

    if domain_object.boundary_map is None:
        return 0

    all_ids = []
    for tag, boundary in domain_object.boundary_map.items():
        if boundary is not None and boundary.__class__.__name__ == 'Absorbing_wave_boundary':
            segment_edges = domain_object.tag_boundary_cells.get(tag, None)
            if segment_edges is not None and len(segment_edges) > 0:
                all_ids.extend(segment_edges)

    if len(all_ids) == 0:
        return 0

    ids = np.array(all_ids, dtype=np.intc)
    num_edges = len(ids)
    boundary_indices = ids
    vol_ids_arr  = np.ascontiguousarray(domain_object.boundary_cells[ids],  dtype=np.intc)
    edge_ids_arr = np.ascontiguousarray(domain_object.boundary_edges[ids], dtype=np.intc)

    if     gpu_absorbing_wave_init(&gpu_dom.GD, num_edges,
                               &boundary_indices[0], &vol_ids_arr[0], &edge_ids_arr[0]) != 0:
        raise RuntimeError('init_absorbing_wave_boundary: ' + get_last_error(gpu_dom))
    return num_edges


def set_absorbing_wave_value(GPUDomain gpu_dom, double wave_value):
    """Update the wave stage scalar for Absorbing_wave_boundary (called each timestep)."""
    gpu_absorbing_wave_set_value(&gpu_dom.GD, wave_value)


def evaluate_absorbing_wave_boundary_gpu(GPUDomain gpu_dom):
    """Evaluate Absorbing_wave_boundary on GPU using the current wave_value."""
    gpu_evaluate_absorbing_wave_boundary(&gpu_dom.GD)


def init_characteristic_wave_boundary(GPUDomain gpu_dom, object domain_object):
    """
    Initialize Characteristic_wave_boundary for GPU.

    Scans domain.boundary_map for Characteristic_wave_boundary objects.
    If multiple such boundaries exist with different background_stage values,
    only the first background_stage encountered is used for the GPU struct
    (consistent with the scalar-update pattern used by Time_boundary).
    Call once after domain setup, BEFORE map_to_gpu.
    """
    cdef int num_edges = 0
    cdef np.ndarray[int, ndim=1, mode="c"] boundary_indices
    cdef np.ndarray[int, ndim=1, mode="c"] vol_ids_arr
    cdef np.ndarray[int, ndim=1, mode="c"] edge_ids_arr

    if domain_object.boundary_map is None:
        return 0

    all_ids = []
    bg_stage = 0.0
    found_first = False
    for tag, boundary in domain_object.boundary_map.items():
        if boundary is not None and boundary.__class__.__name__ == 'Characteristic_wave_boundary':
            if not found_first:
                bg_stage = float(boundary.background_stage)
                found_first = True
            segment_edges = domain_object.tag_boundary_cells.get(tag, None)
            if segment_edges is not None and len(segment_edges) > 0:
                all_ids.extend(segment_edges)

    if len(all_ids) == 0:
        return 0

    ids = np.array(all_ids, dtype=np.intc)
    num_edges = len(ids)
    boundary_indices = ids
    vol_ids_arr  = np.ascontiguousarray(domain_object.boundary_cells[ids],  dtype=np.intc)
    edge_ids_arr = np.ascontiguousarray(domain_object.boundary_edges[ids], dtype=np.intc)

    if     gpu_characteristic_wave_init(&gpu_dom.GD, num_edges,
                                    &boundary_indices[0], &vol_ids_arr[0], &edge_ids_arr[0],
                                    bg_stage) != 0:
        raise RuntimeError('init_characteristic_wave_boundary: ' + get_last_error(gpu_dom))
    return num_edges


def set_characteristic_wave_value(GPUDomain gpu_dom, double wave_value):
    """Update the perturbation scalar for Characteristic_wave_boundary (called each timestep)."""
    gpu_characteristic_wave_set_value(&gpu_dom.GD, wave_value)


def evaluate_characteristic_wave_boundary_gpu(GPUDomain gpu_dom):
    """Evaluate Characteristic_wave_boundary on GPU using the current wave_value."""
    gpu_evaluate_characteristic_wave_boundary(&gpu_dom.GD)


def init_flather_boundary(GPUDomain gpu_dom, object domain_object):
    """
    Initialize Flather_external_stage_zero_velocity_boundary for GPU.

    Scans domain.boundary_map for Flather_external_stage_zero_velocity_boundary
    objects, collects their edges, and initialises the GPU struct.
    Call once after domain setup, BEFORE map_to_gpu.

    Returns
    -------
    int
        Number of boundary edges found (0 if none).
    """
    cdef int num_edges = 0
    cdef np.ndarray[int, ndim=1, mode="c"] boundary_indices
    cdef np.ndarray[int, ndim=1, mode="c"] vol_ids_arr
    cdef np.ndarray[int, ndim=1, mode="c"] edge_ids_arr

    if domain_object.boundary_map is None:
        return 0

    all_ids = []
    for tag, boundary in domain_object.boundary_map.items():
        if boundary is not None and \
                boundary.__class__.__name__ == 'Flather_external_stage_zero_velocity_boundary':
            segment_edges = domain_object.tag_boundary_cells.get(tag, None)
            if segment_edges is not None and len(segment_edges) > 0:
                all_ids.extend(segment_edges)

    if len(all_ids) == 0:
        return 0

    ids = np.array(all_ids, dtype=np.intc)
    num_edges = len(ids)
    boundary_indices = ids
    vol_ids_arr  = np.ascontiguousarray(domain_object.boundary_cells[ids],  dtype=np.intc)
    edge_ids_arr = np.ascontiguousarray(domain_object.boundary_edges[ids], dtype=np.intc)

    if     gpu_flather_init(&gpu_dom.GD, num_edges,
                        &boundary_indices[0], &vol_ids_arr[0], &edge_ids_arr[0]) != 0:
        raise RuntimeError('init_flather_boundary: ' + get_last_error(gpu_dom))
    return num_edges


def set_flather_value(GPUDomain gpu_dom, double stage_outside):
    """Update the exterior stage scalar for Flather boundary (called each timestep)."""
    gpu_flather_set_value(&gpu_dom.GD, stage_outside)


def evaluate_flather_boundary_gpu(GPUDomain gpu_dom):
    """Evaluate Flather_external_stage_zero_velocity_boundary on GPU."""
    gpu_evaluate_flather_boundary(&gpu_dom.GD)


def evolve_one_euler_step_gpu(GPUDomain gpu_dom, double max_timestep, int apply_forcing):
    """
    Execute one Euler timestep on GPU.

    Parameters
    ----------
    gpu_dom : GPUDomain
        The GPU domain wrapper
    max_timestep : float
        Maximum allowed timestep (respecting yieldstep/finaltime constraints)
    apply_forcing : int
        Whether to apply GPU-compatible forcing terms

    Returns
    -------
    float
        The timestep used
    """
    return gpu_evolve_one_euler_step(&gpu_dom.GD, max_timestep, apply_forcing)


def evolve_one_rk2_step_gpu(GPUDomain gpu_dom, double max_timestep, int apply_forcing):
    """
    Execute one RK2 timestep on GPU.

    Parameters
    ----------
    gpu_dom : GPUDomain
        The GPU domain wrapper
    max_timestep : float
        Maximum allowed timestep (respecting yieldstep/finaltime constraints)
    apply_forcing : int
        Whether to apply GPU-compatible forcing terms

    Returns
    -------
    float
        The timestep used
    """
    return gpu_evolve_one_rk2_step(&gpu_dom.GD, max_timestep, apply_forcing)


def evolve_one_rk3_step_gpu(GPUDomain gpu_dom, double max_timestep, int apply_forcing):
    """
    Execute one SSP-RK3 timestep on GPU (Shu-Osher 3-stage).

    Parameters
    ----------
    gpu_dom : GPUDomain
        The GPU domain wrapper
    max_timestep : float
        Maximum allowed timestep (respecting yieldstep/finaltime constraints)
    apply_forcing : int
        Whether to apply GPU-compatible forcing terms

    Returns
    -------
    float
        The timestep used
    """
    return gpu_evolve_one_rk3_step(&gpu_dom.GD, max_timestep, apply_forcing)


def ader_ck_predictor_gpu(GPUDomain gpu_dom, double dt):
    """
    Apply the ADER Cauchy-Kovalewski predictor on GPU — centroid variant.

    Advances centroid values from Q^n to Q^{n+1/2} in-place using
    local SWE time derivatives derived from the reconstructed slopes.
    Must be called after extrapolate_second_order_gpu().

    Parameters
    ----------
    gpu_dom : GPUDomain
        The GPU domain wrapper
    dt : float
        Half-timestep (dt/2) to advance centroids by
    """
    gpu_ader_ck_predictor(&gpu_dom.GD, dt)


def ader_ck_predictor_edge_gpu(GPUDomain gpu_dom, double dt):
    """
    Apply the ADER fused edge predictor on GPU.

    Shifts edge values from Q^n to Q^{n+1/2} in-place while leaving centroid
    values unchanged.  Must be called after extrapolate_second_order_gpu().
    No backup/restore of centroid values is needed before/after this call.

    Parameters
    ----------
    gpu_dom : GPUDomain
        The GPU domain wrapper
    dt : float
        Half-timestep (dt/2) to advance edge values by
    """
    gpu_ader_ck_predictor_edge(&gpu_dom.GD, dt)


def evolve_one_ader2_step_gpu(GPUDomain gpu_dom, double max_timestep, int apply_forcing, double prev_dt):
    """
    Execute one ADER-2 timestep on GPU.

    Single-flux-call variant using the fused edge predictor with prev_dt from
    the previous step.  Matches the CPU ADER-2 implementation exactly.
    Cost: 1 extrapolation + 1 edge predictor + 1 flux call.

    Parameters
    ----------
    gpu_dom : GPUDomain
        The GPU domain wrapper
    max_timestep : float
        Maximum allowed timestep (respecting yieldstep/finaltime constraints)
    apply_forcing : int
        Whether to apply GPU-compatible forcing terms
    prev_dt : float
        Timestep from the previous step (pass 0.0 on first call to bootstrap
        with a plain Euler step)

    Returns
    -------
    float
        The timestep used (save as prev_dt for the next call)
    """
    return gpu_evolve_one_ader2_step(&gpu_dom.GD, max_timestep, apply_forcing, prev_dt)


def finalize_gpu_domain(GPUDomain gpu_dom):
    """
    Clean up GPU domain resources.

    This is called automatically when GPUDomain is garbage collected,
    but can be called explicitly for deterministic cleanup.
    """
    if gpu_dom.initialized:
        gpu_domain_finalize(&gpu_dom.GD)
        gpu_dom.initialized = False


# ============================================================================
# GPU Kernel Wrappers
# ============================================================================

def extrapolate_second_order_gpu(GPUDomain gpu_dom):
    """
    Perform second-order edge extrapolation on GPU.

    This extrapolates centroid values to edge values using a limited
    gradient reconstruction. Handles wet-dry fronts with adaptive limiting.
    """
    gpu_extrapolate_second_order(&gpu_dom.GD)


def compute_fluxes_gpu(GPUDomain gpu_dom, int substep_count=0, int timestep_fluxcalls=1):
    """
    Compute fluxes across all edges on GPU.

    Uses the central upwind Kurganov-Noelle-Petrova scheme.

    substep_count / timestep_fluxcalls index domain.boundary_flux_sum so the
    Python boundary_flux_integral_operator gets each RK substep's boundary flux.
    The defaults (0, 1) suit a single flux call (euler / ader2 / a standalone
    compute_fluxes()); multi-substep callers must pass the substep index, e.g.
    rk2 -> (0,2),(1,2) and rk3 -> (0,3),(1,3),(2,3).

    Returns
    -------
    float
        The local minimum timestep (caller should do MPI_Allreduce for global min)
    """
    return gpu_compute_fluxes(&gpu_dom.GD, substep_count, timestep_fluxcalls)


def update_conserved_quantities_gpu(GPUDomain gpu_dom, double timestep):
    """
    Update conserved quantities using explicit and semi-implicit updates.

    Q_new = (Q_old + timestep * explicit_update) / (1 - timestep * semi_implicit_update)
    """
    gpu_update_conserved_quantities(&gpu_dom.GD, timestep)


def backup_conserved_quantities_gpu(GPUDomain gpu_dom):
    """
    Backup centroid values for RK2 timestepping.
    """
    gpu_backup_conserved_quantities(&gpu_dom.GD)


def saxpy_conserved_quantities_gpu(GPUDomain gpu_dom, double a, double b):
    """
    RK2 combination: Q = a*Q_current + b*Q_backup

    Typically called with a=0.5, b=0.5 for standard RK2.
    """
    gpu_saxpy_conserved_quantities(&gpu_dom.GD, a, b)


def saxpy3_conserved_quantities_gpu(GPUDomain gpu_dom, double a, double b, double c):
    """
    RK3 final combination: Q = (a*Q_current + b*Q_backup) / c

    Used for SSP-RK3 final step: saxpy3(2.0, 1.0, 3.0)
    computes Q = (2*Q_current + Q_backup) / 3.
    """
    gpu_saxpy3_conserved_quantities(&gpu_dom.GD, a, b, c)


def apply_sediment_source_gpu(GPUDomain gpu_dom, double timestep):
    """Fractional step on the device: apply E-D to m and the bed change to z.

    Runs the same core kernel as mode 1; on a GPU build its loop is an
    'omp target' region, so the tracer and bed arrays are updated in place on
    the device with no host round trip. That is what keeps this operator
    GPU-safe and off the _gpu_host_writes_suppressed fallback path.
    """
    core_apply_sediment_source(&gpu_dom.GD.D, timestep)


def apply_bedload_gpu(GPUDomain gpu_dom, double timestep):
    """Bedload divergence [G-5] on the device."""
    core_apply_bedload(&gpu_dom.GD.D, timestep)


def apply_repose_gpu(GPUDomain gpu_dom):
    """Angle-of-repose relaxation on the device. Returns sweeps used."""
    return core_apply_repose(&gpu_dom.GD.D)


def protect_gpu(GPUDomain gpu_dom):
    """
    Protect against negative water depths.

    Sets stage = bed where water depth would be negative,
    and zeros momentum where depth is very small.

    Returns
    -------
    float
        Mass error (total volume added to prevent negative depths)
    """
    return gpu_protect(&gpu_dom.GD)


def compute_water_volume_gpu(GPUDomain gpu_dom):
    """
    Compute total water volume on GPU.

    Returns local volume - caller should do MPI_Allreduce for global sum.

    Returns
    -------
    float
        Local water volume (m^3)
    """
    return gpu_compute_water_volume(&gpu_dom.GD)


def set_use_active_set_gpu(GPUDomain gpu_dom, int flag):
    """Push the active-set flag into an already-initialized GPU domain."""
    gpu_dom.GD.use_active_set = flag


def active_set_stats_gpu(GPUDomain gpu_dom):
    """Return (mean_active_fraction, n_rebuilds) for active-set stepping.

    mean_active_fraction is 1.0 when the active set never engaged (disabled,
    MPI fallback, or no rebuild has run yet).
    """
    return (gpu_active_set_mean_fraction(&gpu_dom.GD),
            int(gpu_dom.GD.as_samples))


def manning_friction_gpu(GPUDomain gpu_dom):
    """
    Apply Manning friction on GPU (flat, semi-implicit formulation).

    Adds friction contribution to the semi_implicit_update arrays.
    The friction is then applied during update_conserved_quantities.
    """
    gpu_manning_friction(&gpu_dom.GD)


# ============================================================================
# Rate Operator GPU Support
# ============================================================================

def init_rate_operator(GPUDomain gpu_dom,
                       np.ndarray[int, ndim=1, mode="c"] indices,
                       np.ndarray[double, ndim=1, mode="c"] areas,
                       np.ndarray[int, ndim=1, mode="c"] full_indices=None):
    """
    Initialize a rate operator for GPU execution.

    Parameters
    ----------
    gpu_dom : GPUDomain
        The GPU domain wrapper
    indices : ndarray of int
        Triangle indices where rate is applied
    areas : ndarray of double
        Triangle areas (for mass tracking)
    full_indices : ndarray of int, optional
        Indices of "full" (non-ghost) triangles for mass tracking

    Returns
    -------
    int
        Operator ID (use this to apply rate or finalize)
        Returns -1 on error
    """
    cdef int num_indices = len(indices)
    cdef int num_full = 0
    cdef int *full_ptr = NULL

    if full_indices is not None:
        num_full = len(full_indices)
        full_ptr = &full_indices[0]

    ret = gpu_rate_operator_init(&gpu_dom.GD, num_indices, &indices[0],
                                  &areas[0], full_ptr, num_full)
    if ret < 0:
        raise RuntimeError('init_rate_operator: ' + get_last_error(gpu_dom))
    return ret


def finalize_rate_operator(GPUDomain gpu_dom, int op_id):
    """
    Finalize and free a rate operator.

    Parameters
    ----------
    gpu_dom : GPUDomain
        The GPU domain wrapper
    op_id : int
        Operator ID returned by init_rate_operator
    """
    gpu_rate_operator_finalize(&gpu_dom.GD, op_id)


def finalize_all_rate_operators(GPUDomain gpu_dom):
    """
    Finalize all rate operators.

    Call this during domain cleanup.
    """
    gpu_rate_operators_finalize_all(&gpu_dom.GD)


def rate_operator_tracers_capture_gpu(GPUDomain gpu_dom, int op_id):
    """Record the depths a rate operator's cells have before it is applied."""
    return gpu_rate_operator_tracers_capture(&gpu_dom.GD, op_id)


def rate_operator_tracers_apply_gpu(GPUDomain gpu_dom, int op_id, object c_in, object carry):
    """Move the tracers with the water the rate added or removed.

    `c_in` (ns,) is each tracer's inflow concentration, `carry` (ns,) 1 where
    the tracer leaves with removed water. See gpu_rate_operator.c.
    """
    cdef np.ndarray[double, ndim=1, mode="c"] c_np = np.ascontiguousarray(c_in, dtype=np.float64)
    cdef np.ndarray[int, ndim=1, mode="c"] k_np = np.ascontiguousarray(carry, dtype=np.intc)
    cdef int ns = c_np.shape[0]
    if ns == 0:
        return
    gpu_rate_operator_tracers_apply(&gpu_dom.GD, op_id, &c_np[0], &k_np[0], ns)


def apply_rate_operator_gpu(GPUDomain gpu_dom, int op_id,
                            double rate, double factor, double timestep):
    """
    Apply rate operator on GPU.

    This is the GPU equivalent of Rate_operator.__call__().
    Handles both positive rates (rain) and negative rates (extraction).

    Parameters
    ----------
    gpu_dom : GPUDomain
        The GPU domain wrapper
    op_id : int
        Operator ID returned by init_rate_operator
    rate : double
        Rate value in m/s (scalar - get from Python function if time-dependent)
    factor : double
        Conversion factor
    timestep : double
        Current timestep

    Returns
    -------
    double
        Local mass influx (for mass conservation tracking)
    """
    return gpu_rate_operator_apply(&gpu_dom.GD, op_id, rate, factor, timestep)


def apply_rate_operator_array_gpu(GPUDomain gpu_dom, int op_id,
                                  np.ndarray[double, ndim=1, mode="c"] rate_array,
                                  int use_indices_into_rate,
                                  int rate_changed,
                                  double factor, double timestep):
    """
    Apply rate operator with per-cell rate array on GPU.

    This handles quantity-type rates where each cell has its own rate value.

    Parameters
    ----------
    gpu_dom : GPUDomain
        The GPU domain wrapper
    op_id : int
        Operator ID returned by init_rate_operator
    rate_array : ndarray of double
        Per-cell rate values
    use_indices_into_rate : int
        If 1, rate_array is full domain size (index with indices[k])
        If 0, rate_array matches operator indices size (index with k)
    rate_changed : int
        If 1, transfer new data to GPU; if 0, reuse cached data on GPU
    factor : double
        Conversion factor
    timestep : double
        Current timestep

    Returns
    -------
    double
        Local mass influx (for mass conservation tracking)
    """
    cdef int rate_size = len(rate_array)
    return gpu_rate_operator_apply_array(&gpu_dom.GD, op_id,
                                         &rate_array[0], rate_size,
                                         use_indices_into_rate,
                                         rate_changed,
                                         factor, timestep)


# ============================================================================
# Max-Quantities Operator GPU Support
# ============================================================================

def init_max_quantities_gpu(GPUDomain gpu_dom, int n, double velocity_zero_height):
    """
    Initialise device-resident max-quantities arrays.

    Allocates four arrays of length *n* on the device (stage, depth, speed,
    momentum magnitude), initialised to large-negative / zero as appropriate.
    Call once before the evolve loop starts.

    Returns 0 on success, -1 on allocation failure.
    """
    ret = gpu_max_quantities_init(&gpu_dom.GD, n, velocity_zero_height)
    if ret < 0:
        raise RuntimeError('init_max_quantities_gpu: ' + get_last_error(gpu_dom))
    return ret


def update_max_quantities_gpu(GPUDomain gpu_dom):
    """
    Run one pass of the max-quantities kernel on the GPU.

    Reads stage/bed/xmom/ymom from device memory (already mapped via
    gpu_domain_map_arrays) and updates the four device-resident max arrays.
    No host<->device transfer is performed.
    """
    gpu_max_quantities_update(&gpu_dom.GD)


def get_max_quantities_gpu(GPUDomain gpu_dom,
                           np.ndarray[double, ndim=1, mode="c"] max_stage,
                           np.ndarray[double, ndim=1, mode="c"] max_depth,
                           np.ndarray[double, ndim=1, mode="c"] max_speed,
                           np.ndarray[double, ndim=1, mode="c"] max_uh):
    """
    Sync max arrays from device to host and copy into the supplied numpy arrays.

    The four output arrays must each have length == domain.number_of_elements.
    This is the only call that incurs a device-to-host transfer; call it only
    at yield steps (SWW write) or at end-of-simulation export.
    """
    gpu_max_quantities_get(&gpu_dom.GD,
                           &max_stage[0], &max_depth[0],
                           &max_speed[0], &max_uh[0])


def get_max_tracers_gpu(GPUDomain gpu_dom,
                        np.ndarray[double, ndim=2, mode="c"] max_tracer):
    """
    Sync the per-tracer running maxima from device to host.

    max_tracer must be (n_tracers, number_of_elements) and C-contiguous --
    the same layout the tracer arrays use, which is what the kernel writes.
    Like get_max_quantities_gpu this is a device-to-host transfer; call it
    only at yield steps or on export.
    """
    if max_tracer.shape[0] == 0:
        return
    gpu_max_tracers_get(&gpu_dom.GD, &max_tracer[0, 0])


def finalize_max_quantities_gpu(GPUDomain gpu_dom):
    """
    Unmap device arrays and free host memory for the max-quantities operator.
    Called automatically by gpu_domain_finalize; safe to call earlier.
    """
    gpu_max_quantities_finalize(&gpu_dom.GD)


# ============================================================================
# Inlet Operator GPU Support
# ============================================================================

def init_inlet_operator(GPUDomain gpu_dom,
                        np.ndarray[int, ndim=1, mode="c"] indices,
                        np.ndarray[double, ndim=1, mode="c"] areas):
    """
    Initialize an inlet operator for GPU execution.

    Parameters
    ----------
    gpu_dom : GPUDomain
        The GPU domain wrapper
    indices : ndarray of int
        Triangle indices where inlet operates
    areas : ndarray of double
        Triangle areas

    Returns
    -------
    int
        Operator ID (use this to apply or finalize)
        Returns -1 on error
    """
    cdef int num_indices = len(indices)
    ret = gpu_inlet_operator_init(&gpu_dom.GD, num_indices, &indices[0], &areas[0])
    if ret < 0:
        raise RuntimeError('init_inlet_operator: ' + get_last_error(gpu_dom))
    return ret


def finalize_inlet_operator(GPUDomain gpu_dom, int op_id):
    """Finalize and free an inlet operator."""
    gpu_inlet_operator_finalize(&gpu_dom.GD, op_id)


def finalize_all_inlet_operators(GPUDomain gpu_dom):
    """Finalize all inlet operators."""
    gpu_inlet_operators_finalize_all(&gpu_dom.GD)


def inlet_get_volume_gpu(GPUDomain gpu_dom, int op_id):
    """
    Get total water volume in inlet triangles (GPU reduction).

    Returns scalar volume - no full domain sync needed.
    """
    return gpu_inlet_get_volume(&gpu_dom.GD, op_id)


def inlet_get_velocities_gpu(GPUDomain gpu_dom, int op_id):
    """
    Get area-weighted velocities from inlet triangles (GPU reduction).

    Returns (u_array, v_array) as numpy arrays - small D2H of inlet data only.
    """
    cdef int n = gpu_dom.GD.inlet_ops.ops[op_id].num_indices
    cdef double u_dummy = 0.0
    cdef double v_dummy = 0.0

    gpu_inlet_get_velocities(&gpu_dom.GD, op_id, &u_dummy, &v_dummy)

    # The velocities are stored in scratch buffers - copy to numpy arrays
    cdef double *s_xmom = gpu_dom.GD.inlet_ops.ops[op_id].scratch_xmom
    cdef double *s_ymom = gpu_dom.GD.inlet_ops.ops[op_id].scratch_ymom

    u_arr = np.empty(n, dtype=np.float64)
    v_arr = np.empty(n, dtype=np.float64)
    cdef double[::1] u_view = u_arr
    cdef double[::1] v_view = v_arr
    cdef int k
    for k in range(n):
        u_view[k] = s_xmom[k]
        v_view[k] = s_ymom[k]

    return u_arr, v_arr


def inlet_apply_gpu(GPUDomain gpu_dom, int op_id, double volume,
                    double current_volume, double total_area,
                    object vel_u_arr, object vel_v_arr,
                    int has_velocity, double ext_vel_u, double ext_vel_v,
                    int zero_velocity, object c_in_arr=None, object dmass_arr=None):
    """
    Apply inlet operator on GPU - main entry point.

    Handles all 3 cases (positive volume, negative sustainable, drain), and
    the tracers the water carries: `c_in_arr` is the inflow concentration of
    each tracer (None: all zero); `dmass_arr`, if given, receives the change
    in each tracer's volume in the inlet. Returns actual applied volume.
    """
    cdef np.ndarray[double, ndim=1, mode="c"] cin_np
    cdef np.ndarray[double, ndim=1, mode="c"] dm_np
    cdef const double *cin_ptr = NULL
    cdef double *dm_ptr = NULL
    cdef int n_cin = 0
    if c_in_arr is not None:
        cin_np = np.ascontiguousarray(c_in_arr, dtype=np.float64)
        if cin_np.shape[0] > 0:
            cin_ptr = &cin_np[0]
            n_cin = cin_np.shape[0]
    if dmass_arr is not None:
        dm_np = dmass_arr
        if dm_np.shape[0] > 0:
            dm_ptr = &dm_np[0]
    cdef np.ndarray[double, ndim=1, mode="c"] u_np
    cdef np.ndarray[double, ndim=1, mode="c"] v_np
    cdef double *u_ptr = NULL
    cdef double *v_ptr = NULL
    cdef int n_vel = 0

    if vel_u_arr is not None and vel_v_arr is not None:
        u_np = np.ascontiguousarray(vel_u_arr, dtype=np.float64)
        v_np = np.ascontiguousarray(vel_v_arr, dtype=np.float64)
        u_ptr = &u_np[0]
        v_ptr = &v_np[0]
        n_vel = len(u_np)

    return gpu_inlet_apply(&gpu_dom.GD, op_id, volume,
                           current_volume, total_area,
                           u_ptr, v_ptr, n_vel,
                           has_velocity, ext_vel_u, ext_vel_v,
                           zero_velocity, cin_ptr, n_cin, dm_ptr)


# ============================================================================
# Culvert Operator API (Boyd box/pipe - batched GPU gather/scatter)
# ============================================================================

def init_culvert_operator(GPUDomain gpu_dom,
                          int culvert_type,
                          double width, double height, double diameter,
                          double length, double manning, double sum_loss,
                          double blockage, double barrels,
                          int use_velocity_head, int use_momentum_jet,
                          int use_old_momentum_method,
                          int always_use_Q_wetdry_adjustment,
                          double max_velocity, double smoothing_timescale,
                          outward_vector_0, outward_vector_1,
                          double invert_elevation_0, double invert_elevation_1,
                          int has_invert_elevation_0, int has_invert_elevation_1,
                          int enquiry_index_0, int enquiry_index_1,
                          np.ndarray[int, ndim=1, mode="c"] inlet0_indices,
                          np.ndarray[double, ndim=1, mode="c"] inlet0_areas,
                          np.ndarray[int, ndim=1, mode="c"] inlet1_indices,
                          np.ndarray[double, ndim=1, mode="c"] inlet1_areas,
                          int master_proc=-1,
                          int enquiry_proc_0=-1,
                          int enquiry_proc_1=-1,
                          int inlet_master_proc_0=-1,
                          int inlet_master_proc_1=-1,
                          int is_local=1,
                          int mpi_tag_base=0,
                          double z1=0.0,
                          double z2=0.0,
                          double init_smooth_Q=0.0,
                          double init_smooth_delta_total_energy=0.0):
    """
    Register a culvert (Boyd box, pipe, or weir/orifice trapezoid) with the GPU culvert manager.

    For cross-boundary (parallel) culverts, pass MPI topology:
      master_proc: rank that computes discharge
      enquiry_proc_0/1: ranks owning each enquiry point
      inlet_master_proc_0/1: master rank for each inlet region
      is_local: 0 for cross-boundary, 1 for fully local
      mpi_tag_base: unique MPI tag base for this culvert

    Returns culvert_id (0..63) or -1 on error.
    """
    cdef culvert_params p
    memset(&p, 0, sizeof(culvert_params))

    p.type = culvert_type
    p.g = gpu_dom.GD.D.g
    p.width = width
    p.height = height
    p.diameter = diameter
    p.z1 = z1
    p.z2 = z2
    p.length = length
    p.manning = manning
    p.sum_loss = sum_loss
    p.blockage = blockage
    p.barrels = barrels
    p.use_velocity_head = use_velocity_head
    p.use_momentum_jet = use_momentum_jet
    p.use_old_momentum_method = use_old_momentum_method
    p.always_use_Q_wetdry_adjustment = always_use_Q_wetdry_adjustment
    p.max_velocity = max_velocity
    p.smoothing_timescale = smoothing_timescale
    p.outward_vector_0[0] = outward_vector_0[0]
    p.outward_vector_0[1] = outward_vector_0[1]
    p.outward_vector_1[0] = outward_vector_1[0]
    p.outward_vector_1[1] = outward_vector_1[1]
    p.invert_elevation_0 = invert_elevation_0
    p.invert_elevation_1 = invert_elevation_1
    p.has_invert_elevation_0 = has_invert_elevation_0
    p.has_invert_elevation_1 = has_invert_elevation_1

    cdef int n0 = len(inlet0_indices)
    cdef int n1 = len(inlet1_indices)

    # Default MPI topology: fully local (this rank owns everything)
    if master_proc < 0:
        master_proc = gpu_dom.GD.rank
    if enquiry_proc_0 < 0:
        enquiry_proc_0 = gpu_dom.GD.rank
    if enquiry_proc_1 < 0:
        enquiry_proc_1 = gpu_dom.GD.rank
    if inlet_master_proc_0 < 0:
        inlet_master_proc_0 = gpu_dom.GD.rank
    if inlet_master_proc_1 < 0:
        inlet_master_proc_1 = gpu_dom.GD.rank

    cdef int *p_inlet0 = &inlet0_indices[0] if n0 > 0 else NULL
    cdef double *p_area0 = &inlet0_areas[0] if n0 > 0 else NULL
    cdef int *p_inlet1 = &inlet1_indices[0] if n1 > 0 else NULL
    cdef double *p_area1 = &inlet1_areas[0] if n1 > 0 else NULL

    ret = gpu_culvert_init(&gpu_dom.GD, &p,
                            enquiry_index_0, enquiry_index_1,
                            n0, p_inlet0, p_area0,
                            n1, p_inlet1, p_area1,
                            master_proc, enquiry_proc_0, enquiry_proc_1,
                            inlet_master_proc_0, inlet_master_proc_1,
                            is_local, mpi_tag_base,
                            init_smooth_Q, init_smooth_delta_total_energy)
    if ret < 0:
        raise RuntimeError('init_culvert_operator: ' + get_last_error(gpu_dom))
    return ret


def finalize_culvert_operator(GPUDomain gpu_dom, int culvert_id):
    """Finalize a single culvert operator."""
    gpu_culvert_finalize(&gpu_dom.GD, culvert_id)


def finalize_all_culvert_operators(GPUDomain gpu_dom):
    """Finalize all culvert operators and free scratch buffers."""
    gpu_culverts_finalize_all(&gpu_dom.GD)


def map_culvert_operators(GPUDomain gpu_dom):
    """Map culvert scratch buffers to GPU. Call after all culverts registered."""
    if gpu_culverts_map(&gpu_dom.GD) != 0:
        raise RuntimeError('map_culverts: ' + get_last_error(gpu_dom))


def apply_all_culvert_operators(GPUDomain gpu_dom, double timestep):
    """
    Execute all registered culverts for one timestep.

    Performs batched GPU gather -> CPU discharge calc -> GPU scatter.
    Only 2 GPU sync points regardless of number of culverts.
    """
    gpu_culverts_apply_all(&gpu_dom.GD, timestep)


def get_culvert_report_stats(GPUDomain gpu_dom, int culvert_id):
    """Read back a culvert's per-step reporting stats as a tuple
    ``(gain, discharge, velocity, driving_energy, delta_total_energy)``.

    Values are non-zero only on the proc that computed the discharge (the
    culvert's master proc); other procs get zeros. Used by GPUCulvertManager to
    mirror the C-computed stats onto the Python operator objects so mode-2
    structure ``.log`` files match mode-1.
    """
    cdef double out[5]
    cdef int rc = gpu_culverts_get_report(&gpu_dom.GD, culvert_id, out)
    if rc != 0:
        return (0.0, 0.0, 0.0, 0.0, 0.0)
    return (out[0], out[1], out[2], out[3], out[4])


# ============================================================================
# FLOP Counter API (Gordon Bell Performance Profiling)
# ============================================================================

def flop_counters_reset(GPUDomain gpu_dom):
    """
    Reset all FLOP counters to zero.

    Call this at the start of the profiling period.
    """
    gpu_flop_counters_reset(&gpu_dom.GD)


def flop_counters_enable(GPUDomain gpu_dom, bint enable):
    """
    Enable or disable FLOP counting.

    Parameters
    ----------
    gpu_dom : GPUDomain
        The GPU domain wrapper
    enable : bool
        True to enable counting, False to disable
    """
    gpu_flop_counters_enable(&gpu_dom.GD, 1 if enable else 0)


def flop_counters_start_timer(GPUDomain gpu_dom):
    """
    Start the FLOP counter timer.

    Call this at the start of the profiling period.
    """
    gpu_flop_counters_start_timer(&gpu_dom.GD)


def flop_counters_stop_timer(GPUDomain gpu_dom):
    """
    Stop the FLOP counter timer.

    Call this at the end of the profiling period.
    """
    gpu_flop_counters_stop_timer(&gpu_dom.GD)


def flop_counters_get_total(GPUDomain gpu_dom):
    """
    Get total FLOP count across all kernels.

    Returns
    -------
    int
        Total FLOPs executed since last reset
    """
    return gpu_flop_counters_get_total(&gpu_dom.GD)


def flop_counters_get_flops(GPUDomain gpu_dom):
    """
    Get FLOP rate (FLOPs per second).

    Call flop_counters_stop_timer first to record elapsed time.

    Returns
    -------
    float
        FLOP/s rate
    """
    return gpu_flop_counters_get_flops(&gpu_dom.GD)


def flop_counters_print(GPUDomain gpu_dom):
    """
    Print detailed FLOP counter summary to stdout.

    Shows per-kernel breakdown and total performance in GFLOP/s.
    """
    gpu_flop_counters_print(&gpu_dom.GD)


def flop_counters_get_stats(GPUDomain gpu_dom):
    """
    Get all FLOP counter statistics as a dictionary.

    Returns
    -------
    dict
        Dictionary with per-kernel FLOPs, total, elapsed time, and GFLOP/s
    """
    cdef uint64_t total = gpu_flop_counters_get_total(&gpu_dom.GD)
    cdef double flops = gpu_flop_counters_get_flops(&gpu_dom.GD)

    return {
        'extrapolate': gpu_flop_counters_get_extrapolate(&gpu_dom.GD),
        'compute_fluxes': gpu_flop_counters_get_compute_fluxes(&gpu_dom.GD),
        'update': gpu_flop_counters_get_update(&gpu_dom.GD),
        'protect': gpu_flop_counters_get_protect(&gpu_dom.GD),
        'manning': gpu_flop_counters_get_manning(&gpu_dom.GD),
        'backup': gpu_flop_counters_get_backup(&gpu_dom.GD),
        'saxpy': gpu_flop_counters_get_saxpy(&gpu_dom.GD),
        'rate_operator': gpu_flop_counters_get_rate_operator(&gpu_dom.GD),
        'ghost_exchange': gpu_flop_counters_get_ghost_exchange(&gpu_dom.GD),
        'total_flops': total,
        'flops_per_second': flops,
        'gflops_per_second': flops / 1.0e9,
        'total_gflops': total / 1.0e9,
    }


# ============================================================================
# Global FLOP Counter API (MPI reduction for multi-GPU)
# ============================================================================

def flop_counters_get_global_total(GPUDomain gpu_dom):
    """
    Get global total FLOP count summed across all MPI ranks/GPUs.

    This performs an MPI_Allreduce to sum FLOPs from all ranks.

    Returns
    -------
    int
        Total FLOPs across all GPUs since last reset
    """
    return gpu_flop_counters_get_global_total(&gpu_dom.GD)


def flop_counters_get_global_flops(GPUDomain gpu_dom):
    """
    Get global FLOP rate (FLOPs per second) across all GPUs.

    Call flop_counters_stop_timer first to record elapsed time.

    Returns
    -------
    float
        Global FLOP/s rate (sum of all GPU FLOPs / elapsed time)
    """
    return gpu_flop_counters_get_global_flops(&gpu_dom.GD)


def flop_counters_print_global(GPUDomain gpu_dom):
    """
    Print global FLOP counter summary to stdout (rank 0 only).

    Shows per-kernel breakdown summed across all GPUs, with percentages.
    Includes total GFLOP/s and per-GPU average.
    """
    gpu_flop_counters_print_global(&gpu_dom.GD)


def flop_counters_get_global_stats(GPUDomain gpu_dom):
    """
    Get global FLOP counter statistics as a dictionary (MPI reduced).

    Returns
    -------
    dict
        Dictionary with global total, elapsed time, GFLOP/s, and per-GPU average
    """
    cdef uint64_t global_total = gpu_flop_counters_get_global_total(&gpu_dom.GD)
    cdef double global_flops = gpu_flop_counters_get_global_flops(&gpu_dom.GD)
    cdef int nprocs = gpu_dom.GD.nprocs

    return {
        'global_total_flops': global_total,
        'global_total_gflops': global_total / 1.0e9,
        'global_flops_per_second': global_flops,
        'global_gflops_per_second': global_flops / 1.0e9,
        'num_gpus': nprocs,
        'per_gpu_gflops_per_second': (global_flops / 1.0e9) / nprocs if nprocs > 0 else 0.0,
    }
