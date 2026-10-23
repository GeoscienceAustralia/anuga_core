// GPU-accelerated shallow water solver
// Split from sw_domain_gpu.c for maintainability

#include <stdio.h>
#include <stdarg.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <omp.h>
// MPI (or single-process stubs) come in via gpu_domain.h
#include "gpu_domain.h"
#ifdef HAVE_MPI_EXT_H
#include <mpi-ext.h>
#endif
#include "gpu_omp_macros.h"

// Domain initialization, memory management, and sync functions

// Forward declarations for functions in other files
extern void gpu_flop_counters_init(struct gpu_domain *GD);
extern void gpu_halo_finalize(struct gpu_domain *GD);
extern void gpu_reflective_finalize(struct gpu_domain *GD);
extern void gpu_file_boundary_finalize(struct gpu_domain *GD);
extern void gpu_dirichlet_finalize(struct gpu_domain *GD);
extern void gpu_transmissive_finalize(struct gpu_domain *GD);
extern void gpu_transmissive_n_zero_t_finalize(struct gpu_domain *GD);
extern void gpu_boundary_edge_sync_finalize(struct gpu_domain *GD);

// ============================================================================
// Utility Functions
// ============================================================================

// ---------------------------------------------------------------------------
// Device memory estimation and checking
// ---------------------------------------------------------------------------

// Compute the number of bytes that gpu_domain_map_arrays will transfer to the
// device for a domain with n triangles and nb boundary edges.
//
// Counts only the main domain arrays (double and int64 per-element arrays).
// Boundary and riverwall arrays are small relative to domain size.
size_t gpu_estimate_required_memory(anuga_int n, anuga_int nb) {
    // Per-element double arrays mapped in gpu_domain_map_arrays:
    //   centroid values:       stage, xmom, ymom, bed, height, friction  =  6n
    //   edge values:           stage, xmom, ymom, bed, height             =  5 x 3n = 15n
    //   explicit updates:      stage, xmom, ymom                          =  3n
    //   semi-implicit updates: stage, xmom, ymom                          =  3n
    //   mesh geometry:         normals(6n), edgelengths(3n), areas(n),
    //                          radii(n), max_speed(n), centroid_coords(2n),
    //                          edge_coords(6n), x/y_centroid_work(2n)     = 23n
    //   Total doubles: (6 + 15 + 3 + 3 + 23) = 50n doubles
    //
    // Per-element int64 arrays:
    //   neighbours(3n), neighbour_edges(3n), surrogate_neighbours(3n),
    //   number_of_boundaries(n), tri_full_flag(n)                         = 11n int64
    //
    // Backup values (double): stage, xmom, ymom = 3n  (mapped in map_backup_arrays)
    //   — not mapped here but count as we typically allocate them
    //
    // Boundary arrays (small): 5 x nb doubles
    size_t double_bytes = (size_t)(50) * (size_t)n * sizeof(double);
    size_t int_bytes    = (size_t)(11) * (size_t)n * sizeof(int64_t);
    size_t bndry_bytes  = (size_t)(5)  * (size_t)nb * sizeof(double);
    return double_bytes + int_bytes + bndry_bytes;
}

// Query free and total memory on the current OpenMP target device.
// Returns 1 (OK) if enough memory is available, 0 if not (or if query fails).
// Sets *free_bytes and *total_bytes when available (0 when not queryable).
//
// Vendor-specific query APIs are selected at compile time:
//   -DUSE_CUDA  -> cudaMemGetInfo
//   -DUSE_HIP   -> hipMemGetInfo
// CPU_ONLY_MODE and unknown targets: query is skipped, returns 1.
int gpu_query_device_memory(size_t *free_bytes, size_t *total_bytes) {
    *free_bytes  = 0;
    *total_bytes = 0;

#if defined(USE_CUDA)
    cudaError_t err = cudaMemGetInfo(free_bytes, total_bytes);
    if (err != cudaSuccess) {
        return 1;   // Can't query — don't block
    }
    return 1;  // Caller decides if enough
#elif defined(USE_HIP)
    hipError_t err = hipMemGetInfo(free_bytes, total_bytes);
    if (err != hipSuccess) {
        return 1;
    }
    return 1;
#else
    // No vendor API available (CPU_ONLY_MODE or unknown target).
    // Report unlimited — mapping will silently run on host.
    return 1;
#endif
}

// Check that the device has enough free memory for the domain.
// Prints a clear error message and returns 0 if insufficient.
// Always returns 1 in CPU_ONLY_MODE (no real device to OOM).
int gpu_check_device_memory(struct gpu_domain *GD) {
    anuga_int n  = GD->D.number_of_elements;
    anuga_int nb = GD->D.boundary_length;
    size_t required = gpu_estimate_required_memory(n, nb);
    size_t free_bytes = 0, total_bytes = 0;

    gpu_query_device_memory(&free_bytes, &total_bytes);

    double req_mb   = required    / (1024.0 * 1024.0);
    double free_mb  = free_bytes  / (1024.0 * 1024.0);
    double total_mb = total_bytes / (1024.0 * 1024.0);

    // Print memory info when verbose or when device memory is known
    if (GD->verbose || (total_bytes > 0)) {
        if (GD->rank == 0) {
            if (total_bytes > 0) {
                printf("  GPU memory: %.0f MB required  |  %.0f / %.0f MB free/total\n",
                       req_mb, free_mb, total_mb);
            } else {
                printf("  GPU memory: %.0f MB estimated required (device query not available)\n",
                       req_mb);
            }
            fflush(stdout);
        }
    }

    // Only fail if we have a real memory figure AND it's insufficient
    if (total_bytes > 0 && free_bytes < required) {
        gpu_set_error(GD,
            "\n[ANUGA GPU] ERROR (rank %d): Insufficient GPU memory.\n"
            "  Domain has %" PRId64 " triangles, estimated %.0f MB required.\n"
            "  GPU has %.0f MB free of %.0f MB total.\n"
            "  Options:\n"
            "    1. Reduce domain size (use fewer triangles)\n"
            "    2. Use more MPI ranks to distribute the domain\n"
            "    3. Use set_multiprocessor_mode(1) to run without GPU offloading\n",
            GD->rank, n, req_mb, free_mb, total_mb);
        fflush(stderr);
        return 0;
    }
    return 1;
}

// Detect if MPI library supports GPU-aware communication.
// Priority:
//   1. Compile-time -DGPU_AWARE_MPI flag (user asserts support is present).
//   2. Runtime MPIX_Query_cuda_support()  (Open MPI / MVAPICH2 CUDA builds).
//   3. Runtime MPIX_Query_rocm_support()  (Open MPI / MVAPICH2 ROCm builds).
//   4. Default: 0 (unknown, assume no GPU-aware MPI).
int detect_gpu_aware_mpi(void) {
#ifdef GPU_AWARE_MPI
    return 1;
#elif defined(HAVE_MPIX_CUDA_SUPPORT)
    return MPIX_Query_cuda_support();
#elif defined(HAVE_MPIX_ROCM_SUPPORT)
    return MPIX_Query_rocm_support();
#else
    return 0;
#endif
}

// Return 1 if a real GPU offload target is available (not CPU_ONLY_MODE).
// In CPU_ONLY_MODE (standard install, CI), omp_get_num_devices() may still
// report 0; in GPU builds the pragma is active only when devices > 0.
int gpu_is_available(void) {
#ifdef CPU_ONLY_MODE
    return 0;
#else
    return omp_get_num_devices() > 0 ? 1 : 0;
#endif
}

int gpu_get_num_devices(void) {
#ifdef CPU_ONLY_MODE
    return 0;
#else
    return omp_get_num_devices();
#endif
}

// Default offload device control. The kernels use `#pragma omp target` with no
// device() clause, so they run on the OpenMP default device — redirecting it is
// a robust, runtime way to force a GPU build onto the host (CPU) and back.
// All no-ops in CPU_ONLY_MODE (there is only the host device).

int gpu_get_initial_device(void) {
#ifdef CPU_ONLY_MODE
    return 0;
#else
    return omp_get_initial_device();
#endif
}

int gpu_get_default_device(void) {
#ifdef CPU_ONLY_MODE
    return 0;
#else
    return omp_get_default_device();
#endif
}

// Route subsequent target regions to `device_id`. Pass gpu_get_initial_device()
// to run on the host (CPU). No-op in CPU_ONLY_MODE.
void gpu_set_default_device(int device_id) {
#ifdef CPU_ONLY_MODE
    (void)device_id;
#else
    omp_set_default_device(device_id);
#endif
}

// Process-global offload enable flag (default on). When 0, gpu_domain_init keeps
// domains on the host even on a GPU build, so mode 2 ('unified') runs
// CPU-multicore with data mapping and kernel execution consistently on the host.
// Set via anuga.set_gpu_offload(); decide before the first unified domain is
// built (the device is chosen at init, when arrays are mapped).
static int g_gpu_offload_enabled = 1;

void gpu_set_offload_enabled(int enabled) { g_gpu_offload_enabled = enabled ? 1 : 0; }
int gpu_get_offload_enabled(void) { return g_gpu_offload_enabled; }

// The OpenMP device that target regions for this domain should run on. A real
// GPU id when offloading; the host (initial) device when not (device_id < 0,
// e.g. offload disabled on a GPU build). Use this for omp_set_default_device —
// NEVER pass GD->device_id directly, since -1 is not a valid device number and
// would silently re-route work to the default GPU.
int gpu_compute_device(struct gpu_domain *GD) {
#ifdef CPU_ONLY_MODE
    (void)GD;
    return 0;
#else
    return GD->device_id >= 0 ? GD->device_id : omp_get_initial_device();
#endif
}

void print_gpu_domain_info(struct gpu_domain *GD) {
    if (!GD->verbose) return;
    printf("\n--- GPU Domain (rank %d/%d) ---\n", GD->rank, GD->nprocs);
    printf("  Elements: %" PRId64 "  |  Device: %d  |  GPU-aware MPI: %s",
           GD->D.number_of_elements, GD->device_id,
           GD->gpu_aware_mpi ? "yes" : "no");
    if (GD->halo.num_neighbors > 0) {
        printf("  |  Halo procs: %d (send %d / recv %d)",
               GD->halo.num_neighbors,
               GD->halo.total_send_size, GD->halo.total_recv_size);
    }
    printf("\n");
    fflush(stdout);
}


// ============================================================================
// Initialization and Cleanup
// ============================================================================

void gpu_set_error(struct gpu_domain *GD, const char *fmt, ...) {
    va_list ap;
    va_start(ap, fmt);
    vsnprintf(GD->last_error, GPU_LAST_ERROR_LEN, fmt, ap);
    va_end(ap);
    fprintf(stderr, "[ANUGA GPU rank %d] %s\n", GD->rank, GD->last_error);
    fflush(stderr);
}

const char *gpu_get_last_error(const struct gpu_domain *GD) {
    return GD->last_error;
}

int gpu_domain_init(struct gpu_domain *GD, MPI_Comm comm, int rank, int nprocs) {
    GD->last_error[0] = '\0';
    // Store MPI info
    GD->comm = comm;
    GD->rank = rank;
    GD->nprocs = nprocs;

    // Initialize GPU state
    GD->gpu_initialized = 0;
    GD->backup_arrays_mapped = 0;

    // Select GPU device (round-robin if more ranks than GPUs), unless offload
    // has been disabled process-wide (anuga.set_gpu_offload(False)) — then run
    // on the host so data mapping and kernel execution stay consistent on CPU.
    int num_devices = omp_get_num_devices();
    if (num_devices > 0 && g_gpu_offload_enabled) {
        GD->device_id = rank % num_devices;
        omp_set_default_device(GD->device_id);
    } else {
        GD->device_id = -1;
        omp_set_default_device(omp_get_initial_device());
    }

    // Detect GPU-aware MPI
    GD->gpu_aware_mpi = detect_gpu_aware_mpi();

    // Initialize halo to empty
    GD->halo.num_neighbors = 0;
    GD->halo.neighbor_ranks = NULL;
    GD->halo.send_counts = NULL;
    GD->halo.recv_counts = NULL;
    GD->halo.total_send_size = 0;
    GD->halo.total_recv_size = 0;
    GD->halo.flat_send_indices = NULL;
    GD->halo.flat_recv_indices = NULL;
    GD->halo.send_offsets = NULL;
    GD->halo.recv_offsets = NULL;
    GD->halo.send_buffer = NULL;
    GD->halo.recv_buffer = NULL;
    GD->halo.requests = NULL;

    // Initialize reflective boundary to empty
    GD->reflective.num_edges = 0;
    GD->reflective.boundary_indices = NULL;
    GD->reflective.vol_ids = NULL;
    GD->reflective.edge_ids = NULL;
    GD->reflective.mapped = 0;

    // Initialize dirichlet boundary to empty
    GD->dirichlet.num_edges = 0;
    GD->dirichlet.boundary_indices = NULL;
    GD->dirichlet.vol_ids = NULL;
    GD->dirichlet.edge_ids = NULL;
    GD->dirichlet.stage_values = NULL;
    GD->dirichlet.xmom_values = NULL;
    GD->dirichlet.ymom_values = NULL;
    GD->dirichlet.mapped = 0;

    // Initialize transmissive boundary to empty
    GD->transmissive.num_edges = 0;
    GD->transmissive.boundary_indices = NULL;
    GD->transmissive.vol_ids = NULL;
    GD->transmissive.edge_ids = NULL;
    GD->transmissive.use_centroid = 0;
    GD->transmissive.mapped = 0;

    // Initialize transmissive_n_zero_t boundary to empty
    GD->transmissive_n_zero_t.num_edges = 0;
    GD->transmissive_n_zero_t.boundary_indices = NULL;
    GD->transmissive_n_zero_t.vol_ids = NULL;
    GD->transmissive_n_zero_t.edge_ids = NULL;
    GD->transmissive_n_zero_t.stage_value = 0.0;
    GD->transmissive_n_zero_t.mapped = 0;

    // Initialize time_boundary to empty
    GD->time_bdry.num_edges = 0;
    GD->time_bdry.boundary_indices = NULL;
    GD->time_bdry.vol_ids = NULL;
    GD->time_bdry.edge_ids = NULL;
    GD->time_bdry.stage_values = NULL;
    GD->time_bdry.xmom_values = NULL;
    GD->time_bdry.ymom_values = NULL;
    GD->time_bdry.mapped = 0;

    // Initialize file_boundary to empty
    GD->file_bdry.num_edges        = 0;
    GD->file_bdry.boundary_indices = NULL;
    GD->file_bdry.vol_ids          = NULL;
    GD->file_bdry.edge_ids         = NULL;
    GD->file_bdry.stage_values     = NULL;
    GD->file_bdry.xmom_values      = NULL;
    GD->file_bdry.ymom_values      = NULL;
    GD->file_bdry.mapped           = 0;

    // Initialize absorbing_wave boundary to empty
    GD->absorbing_wave.num_edges        = 0;
    GD->absorbing_wave.boundary_indices = NULL;
    GD->absorbing_wave.vol_ids          = NULL;
    GD->absorbing_wave.edge_ids         = NULL;
    GD->absorbing_wave.wave_value       = 0.0;
    GD->absorbing_wave.mapped           = 0;

    // Initialize characteristic_wave boundary to empty
    GD->characteristic_wave.num_edges        = 0;
    GD->characteristic_wave.boundary_indices = NULL;
    GD->characteristic_wave.vol_ids          = NULL;
    GD->characteristic_wave.edge_ids         = NULL;
    GD->characteristic_wave.wave_value       = 0.0;
    GD->characteristic_wave.background_stage = 0.0;
    GD->characteristic_wave.mapped           = 0;

    // Initialize flather boundary to empty
    GD->flather.num_edges        = 0;
    GD->flather.boundary_indices = NULL;
    GD->flather.vol_ids          = NULL;
    GD->flather.edge_ids         = NULL;
    GD->flather.stage_outside    = 0.0;
    GD->flather.mapped           = 0;

    // Initialize boundary edge sync to empty
    GD->edge_sync.num_boundary_cells = 0;
    GD->edge_sync.cell_ids = NULL;
    GD->edge_sync.buf_size = 0;
    GD->edge_sync.stage_buf = NULL;
    GD->edge_sync.xmom_buf = NULL;
    GD->edge_sync.ymom_buf = NULL;
    GD->edge_sync.bed_buf = NULL;
    GD->edge_sync.height_buf = NULL;
    GD->edge_sync.initialized = 0;

    // Default simulation parameters
    GD->CFL = 1.0;
    GD->evolve_max_timestep = 1.0e10;
    GD->fixed_flux_timestep = -1.0;  // Disabled by default

    // Initialize rate operators to empty (heap-allocated on first use)
    GD->rate_ops.ops = NULL;
    GD->rate_ops.num_operators = 0;
    GD->rate_ops.capacity = 0;
    GD->rate_ops.initialized = 0;

    // Initialize inlet operators to empty (heap-allocated on first use)
    GD->inlet_ops.ops = NULL;
    GD->inlet_ops.num_operators = 0;
    GD->inlet_ops.capacity = 0;

    // Initialize culvert operators to empty (heap-allocated on first use)
    GD->culvert_ops.num_culverts = 0;
    GD->culvert_ops.capacity = 0;
    GD->culvert_ops.params = NULL;
    GD->culvert_ops.indices = NULL;
    GD->culvert_ops.state = NULL;
    GD->culvert_ops.initialized = 0;
    GD->culvert_ops.mapped = 0;
    GD->culvert_ops.total_inlet_triangles = 0;
    GD->culvert_ops.scratch_enquiry_indices = NULL;
    GD->culvert_ops.scratch_stage = NULL;
    GD->culvert_ops.scratch_xmom = NULL;
    GD->culvert_ops.scratch_ymom = NULL;
    GD->culvert_ops.scratch_elev = NULL;
    GD->culvert_ops.scratch_inlet_indices = NULL;
    GD->culvert_ops.scratch_inlet_areas = NULL;
    GD->culvert_ops.scratch_slot_start = NULL;
    GD->culvert_ops.scratch_slot_count = NULL;
    GD->culvert_ops.scratch_avg_stage = NULL;
    GD->culvert_ops.scratch_avg_depth = NULL;
    GD->culvert_ops.scratch_avg_xmom = NULL;
    GD->culvert_ops.scratch_avg_ymom = NULL;
    GD->culvert_ops.scratch_slot_shift = NULL;
    GD->culvert_ops.scratch_slot_xmom = NULL;
    GD->culvert_ops.scratch_slot_ymom = NULL;

    // Initialize FLOP counters (Gordon Bell profiling)
    gpu_flop_counters_init(GD);

    return 0;
}

void gpu_domain_finalize(struct gpu_domain *GD) {
    // Unmap GPU arrays if mapped
    if (GD->gpu_initialized) {
        gpu_domain_unmap_arrays(GD);
    }

    // Free halo structures
    gpu_halo_finalize(GD);

    // Free boundary structures
    gpu_reflective_finalize(GD);
    gpu_dirichlet_finalize(GD);
    gpu_transmissive_finalize(GD);
    gpu_transmissive_n_zero_t_finalize(GD);
    gpu_file_boundary_finalize(GD);
    gpu_absorbing_wave_finalize(GD);
    gpu_characteristic_wave_finalize(GD);
    gpu_flather_finalize(GD);

    // Free boundary edge sync structures
    gpu_boundary_edge_sync_finalize(GD);

    // Free rate and inlet operator structures (heap-allocated ops arrays)
    gpu_rate_operators_finalize_all(GD);
    gpu_inlet_operators_finalize_all(GD);

    // Free culvert operator structures
    gpu_culverts_finalize_all(GD);

    // Free max-quantities operator arrays
    gpu_max_quantities_finalize(GD);

    GD->gpu_initialized = 0;
}


// ============================================================================
// GPU Memory Management
// ============================================================================

int gpu_domain_map_arrays(struct gpu_domain *GD) {
    if (GD->gpu_initialized) return 1;
    // Note: Don't return early if device_id < 0. With OMP_TARGET_OFFLOAD=disabled,
    // the OpenMP target directives run on CPU, so we still need to set up the
    // data structures and flags for the boundary/kernel functions to work.

    // Check device memory before attempting any mapping.
    // Returns 0 and prints a clear error if there isn't enough free memory.
    if (!gpu_check_device_memory(GD)) {
        return 0;
    }

    anuga_int n = GD->D.number_of_elements;
    anuga_int nb = GD->D.boundary_length;
    struct halo_exchange *H = &GD->halo;
    struct reflective_boundary *R = &GD->reflective;

    // Extract pointers to local variables (OpenMP target can't handle struct->member syntax)
    double *stage_cv = GD->D.stage_centroid_values;
    double *xmom_cv = GD->D.xmom_centroid_values;
    double *ymom_cv = GD->D.ymom_centroid_values;
    double *bed_cv = GD->D.bed_centroid_values;
    double *height_cv = GD->D.height_centroid_values;
    double *stage_ev = GD->D.stage_edge_values;
    double *xmom_ev = GD->D.xmom_edge_values;
    double *ymom_ev = GD->D.ymom_edge_values;
    double *bed_ev = GD->D.bed_edge_values;
    double *height_ev = GD->D.height_edge_values;
    double *stage_eu = GD->D.stage_explicit_update;
    double *xmom_eu = GD->D.xmom_explicit_update;
    double *ymom_eu = GD->D.ymom_explicit_update;
    double *stage_siu = GD->D.stage_semi_implicit_update;
    double *xmom_siu = GD->D.xmom_semi_implicit_update;
    double *ymom_siu = GD->D.ymom_semi_implicit_update;
    anuga_int *neighbours = GD->D.neighbours;
    anuga_int *neighbour_edges = GD->D.neighbour_edges;
    anuga_geom_t *normals = GD->D.normals;
    anuga_geom_t *edgelengths = GD->D.edgelengths;
    anuga_geom_t *areas = GD->D.areas;
    anuga_geom_t *radii = GD->D.radii;
    double *max_speed = GD->D.max_speed;
    anuga_geom_t *centroid_coords = GD->D.centroid_coordinates;
    anuga_geom_t *edge_coords = GD->D.edge_coordinates;

    // Additional arrays for extrapolation
    anuga_int *surrogate_neighbours = GD->D.surrogate_neighbours;
    anuga_int *number_of_boundaries = GD->D.number_of_boundaries;
    double *x_centroid_work = GD->D.x_centroid_work;
    double *y_centroid_work = GD->D.y_centroid_work;

    // Friction array
    double *friction_cv = GD->D.friction_centroid_values;

    // tri_full_flag for MPI ghost cell identification
    anuga_int *tri_full_flag = GD->D.tri_full_flag;

    // Map all domain arrays to GPU - persistent for entire simulation
    #pragma omp target enter data map(to: \
        stage_cv[0:n], xmom_cv[0:n], ymom_cv[0:n], \
        bed_cv[0:n], height_cv[0:n], friction_cv[0:n], \
        stage_ev[0:3*n], xmom_ev[0:3*n], ymom_ev[0:3*n], \
        bed_ev[0:3*n], height_ev[0:3*n], \
        stage_eu[0:n], xmom_eu[0:n], ymom_eu[0:n], \
        stage_siu[0:n], xmom_siu[0:n], ymom_siu[0:n], \
        neighbours[0:3*n], neighbour_edges[0:3*n], \
        surrogate_neighbours[0:3*n], number_of_boundaries[0:n], \
        x_centroid_work[0:n], y_centroid_work[0:n], \
        normals[0:6*n], edgelengths[0:3*n], \
        areas[0:n], radii[0:n], max_speed[0:n], \
        centroid_coords[0:2*n], edge_coords[0:6*n])

    // Map tri_full_flag only when present (parallel domains only)
    if (tri_full_flag != NULL) {
        #pragma omp target enter data map(to: tri_full_flag[0:n])
    }

    // Map boundary values if present (including bed and height for reflective boundary)
    if (nb > 0) {
        double *stage_bv = GD->D.stage_boundary_values;
        double *xmom_bv = GD->D.xmom_boundary_values;
        double *ymom_bv = GD->D.ymom_boundary_values;
        double *bed_bv = GD->D.bed_boundary_values;
        double *height_bv = GD->D.height_boundary_values;
        #pragma omp target enter data map(to: \
            stage_bv[0:nb], xmom_bv[0:nb], ymom_bv[0:nb], \
            bed_bv[0:nb], height_bv[0:nb])
    }

    // Map reflective boundary arrays if initialized (must be done BEFORE this function is called)
    if (R->num_edges > 0 && R->boundary_indices != NULL) {
        int ne = R->num_edges;
        int *b_idx = R->boundary_indices;
        int *v_ids = R->vol_ids;
        int *e_ids = R->edge_ids;

        #pragma omp target enter data map(to: b_idx[0:ne], v_ids[0:ne], e_ids[0:ne])
        R->mapped = 1;

        if (GD->rank == 0 && GD->verbose) {
            printf("  Reflective boundary: %d edges\n", ne);
            fflush(stdout);
        }
    }

    // Map dirichlet boundary arrays if initialized
    struct dirichlet_boundary *Dir = &GD->dirichlet;
    if (Dir->num_edges > 0 && Dir->boundary_indices != NULL) {
        int ne = Dir->num_edges;
        int *b_idx = Dir->boundary_indices;
        int *v_ids = Dir->vol_ids;
        int *e_ids = Dir->edge_ids;
        double *s_val = Dir->stage_values;
        double *x_val = Dir->xmom_values;
        double *y_val = Dir->ymom_values;

        #pragma omp target enter data map(to: b_idx[0:ne], v_ids[0:ne], e_ids[0:ne], \
                                              s_val[0:ne], x_val[0:ne], y_val[0:ne])
        Dir->mapped = 1;

        if (GD->rank == 0 && GD->verbose) {
            printf("  Dirichlet boundary:  %d edges\n", ne);
            fflush(stdout);
        }
    }

    // Map transmissive boundary arrays if initialized
    struct transmissive_boundary *T = &GD->transmissive;
    if (T->num_edges > 0 && T->boundary_indices != NULL) {
        int ne = T->num_edges;
        int *b_idx = T->boundary_indices;
        int *v_ids = T->vol_ids;
        int *e_ids = T->edge_ids;

        #pragma omp target enter data map(to: b_idx[0:ne], v_ids[0:ne], e_ids[0:ne])
        T->mapped = 1;

        if (GD->rank == 0 && GD->verbose) {
            printf("  Transmissive boundary: %d edges\n", ne);
            fflush(stdout);
        }
    }

    // Map transmissive_n_zero_t boundary arrays if initialized
    struct transmissive_n_zero_t_boundary *Tnzt = &GD->transmissive_n_zero_t;
    if (Tnzt->num_edges > 0 && Tnzt->boundary_indices != NULL) {
        int ne = Tnzt->num_edges;
        int *b_idx = Tnzt->boundary_indices;
        int *v_ids = Tnzt->vol_ids;
        int *e_ids = Tnzt->edge_ids;

        #pragma omp target enter data map(to: b_idx[0:ne], v_ids[0:ne], e_ids[0:ne])
        Tnzt->mapped = 1;

        if (GD->rank == 0 && GD->verbose) {
            printf("  Transmissive_n_zero_t boundary: %d edges\n", ne);
            fflush(stdout);
        }
    }

    // Map time_boundary arrays if initialized
    struct time_boundary *TB = &GD->time_bdry;
    if (TB->num_edges > 0 && TB->boundary_indices != NULL) {
        int ne = TB->num_edges;
        int *b_idx = TB->boundary_indices;
        int *v_ids = TB->vol_ids;
        int *e_ids = TB->edge_ids;
        double *s_val = TB->stage_values;
        double *x_val = TB->xmom_values;
        double *y_val = TB->ymom_values;

        #pragma omp target enter data map(to: b_idx[0:ne], v_ids[0:ne], e_ids[0:ne], \
                                              s_val[0:ne], x_val[0:ne], y_val[0:ne])
        TB->mapped = 1;

        if (GD->rank == 0 && GD->verbose) {
            printf("  Time boundary: %d edges\n", ne);
            fflush(stdout);
        }
    }

    // Map file_boundary arrays if initialized
    struct file_boundary *FB = &GD->file_bdry;
    if (FB->num_edges > 0 && FB->boundary_indices != NULL) {
        int ne = FB->num_edges;
        int    *b_idx   = FB->boundary_indices;
        int    *v_ids   = FB->vol_ids;
        int    *e_ids   = FB->edge_ids;
        double *stage_v = FB->stage_values;
        double *xmom_v  = FB->xmom_values;
        double *ymom_v  = FB->ymom_values;

        #pragma omp target enter data map(to: b_idx[0:ne], v_ids[0:ne], e_ids[0:ne], \
                                              stage_v[0:ne], xmom_v[0:ne], ymom_v[0:ne])
        FB->mapped = 1;

        if (GD->rank == 0 && GD->verbose) {
            printf("  File boundary: %d edges\n", ne);
            fflush(stdout);
        }
    }

    // Map absorbing_wave boundary arrays if initialized
    struct absorbing_wave_boundary *AW = &GD->absorbing_wave;
    if (AW->num_edges > 0 && AW->boundary_indices != NULL) {
        int ne = AW->num_edges;
        int *b_idx = AW->boundary_indices;
        int *v_ids = AW->vol_ids;
        int *e_ids = AW->edge_ids;

        #pragma omp target enter data map(to: b_idx[0:ne], v_ids[0:ne], e_ids[0:ne])
        AW->mapped = 1;

        if (GD->rank == 0 && GD->verbose) {
            printf("  Absorbing wave boundary: %d edges\n", ne);
            fflush(stdout);
        }
    }

    // Map characteristic_wave boundary arrays if initialized
    struct characteristic_wave_boundary *CW = &GD->characteristic_wave;
    if (CW->num_edges > 0 && CW->boundary_indices != NULL) {
        int ne = CW->num_edges;
        int *b_idx = CW->boundary_indices;
        int *v_ids = CW->vol_ids;
        int *e_ids = CW->edge_ids;

        #pragma omp target enter data map(to: b_idx[0:ne], v_ids[0:ne], e_ids[0:ne])
        CW->mapped = 1;

        if (GD->rank == 0 && GD->verbose) {
            printf("  Characteristic wave boundary: %d edges\n", ne);
            fflush(stdout);
        }
    }

    // Map flather boundary arrays if initialized
    struct flather_boundary *FL = &GD->flather;
    if (FL->num_edges > 0 && FL->boundary_indices != NULL) {
        int ne = FL->num_edges;
        int *b_idx = FL->boundary_indices;
        int *v_ids = FL->vol_ids;
        int *e_ids = FL->edge_ids;

        #pragma omp target enter data map(to: b_idx[0:ne], v_ids[0:ne], e_ids[0:ne])
        FL->mapped = 1;

        if (GD->rank == 0 && GD->verbose) {
            printf("  Flather boundary: %d edges\n", ne);
            fflush(stdout);
        }
    }

    // Map riverwall arrays if present
    // edge_flux_type is always mapped (size 3*n); other arrays only if riverwalls exist
    anuga_int *edge_flux_type = GD->D.edge_flux_type;
    if (edge_flux_type != NULL) {
        #pragma omp target enter data map(to: edge_flux_type[0:3*n])
    }

    anuga_int n_rw_edges = GD->D.number_of_riverwall_edges;
    if (n_rw_edges > 0) {
        anuga_int *edge_river_wall_counter = GD->D.edge_river_wall_counter;
        double *riverwall_elevation = GD->D.riverwall_elevation;
        anuga_int *riverwall_rowIndex = GD->D.riverwall_rowIndex;
        double *riverwall_hydraulic_properties = GD->D.riverwall_hydraulic_properties;
        anuga_int ncol_hp = GD->D.ncol_riverwall_hydraulic_properties;
        anuga_int nrow_hp = GD->D.nrow_riverwall_hydraulic_properties;

        // edge_river_wall_counter has same size as edge_flux_type (3*n)
        if (edge_river_wall_counter != NULL) {
            #pragma omp target enter data map(to: edge_river_wall_counter[0:3*n])
        }
        // riverwall_elevation and riverwall_rowIndex have size n_rw_edges
        if (riverwall_elevation != NULL) {
            #pragma omp target enter data map(to: riverwall_elevation[0:n_rw_edges])
        }
        if (riverwall_rowIndex != NULL) {
            #pragma omp target enter data map(to: riverwall_rowIndex[0:n_rw_edges])
        }
        // riverwall_hydraulic_properties is (nrow_hp x ncol_hp)
        // nrow_hp = number of unique riverwall segments
        if (riverwall_hydraulic_properties != NULL && ncol_hp > 0 && nrow_hp > 0) {
            #pragma omp target enter data map(to: riverwall_hydraulic_properties[0:nrow_hp*ncol_hp])
        }

        if (GD->rank == 0 && GD->verbose) {
            printf("  Riverwall: %" PRId64 " edges, %" PRId64 " segments\n", n_rw_edges, nrow_hp);
            fflush(stdout);
        }
    }

    // Map backup arrays for RK2 if present
    if (GD->D.stage_backup_values != NULL) {
        double *stage_backup = GD->D.stage_backup_values;
        double *xmom_backup = GD->D.xmom_backup_values;
        double *ymom_backup = GD->D.ymom_backup_values;
        #pragma omp target enter data map(to: \
            stage_backup[0:n], xmom_backup[0:n], ymom_backup[0:n])
        GD->backup_arrays_mapped = 1;
    }

    // Map generic tracer arrays if any tracers are registered.
    //
    // Layout is tracer-major and C-contiguous, matching the indexing the shared
    // kernels use: centroid[s*n + k], edge[s*3n + 3k + i], boundary[s*nb + m].
    // All six must be mapped together: the kernels' n_tracers > 0 guard reads
    // all of them, so a partial mapping is an illegal device address, not a
    // degraded result.
    if (GD->D.number_of_tracers > 0) {
        anuga_int ns = GD->D.number_of_tracers;
        double *tr_cv = GD->D.tracer_centroid_values;
        double *tr_ev = GD->D.tracer_edge_values;
        double *tr_eu = GD->D.tracer_explicit_update;
        double *tr_qv = GD->D.tracer_conserved_values;
        double *tr_bk = GD->D.tracer_backup_values;
        // tracer_boundary_flux is per-cell conservation scratch: the kernel
        // accumulates into it and totals it per substep, so it belongs with the
        // rest of the block for the same reason -- the n_tracers > 0 guard
        // reads it, and a partial mapping is an illegal device address.
        double *tr_bf = GD->D.tracer_boundary_flux;
        #pragma omp target enter data map(to: \
            tr_cv[0:ns*n], tr_ev[0:ns*3*n], \
            tr_eu[0:ns*n], tr_qv[0:ns*n], tr_bk[0:ns*n], \
            tr_bf[0:ns*n])
        // [D-5v] the per-tracer speed factor, only when a profile is on.
        double *tr_sf = GD->D.tracer_speed_factor;
        if (tr_sf != NULL) {
            #pragma omp target enter data map(to: tr_sf[0:ns*n])
        }

        // The external source [G-3]. Allocated with the rest of the block by
        // add_tracer, so it is always present when n_tracers > 0 and is mapped
        // unconditionally like the others. It is read by the source kernel, so
        // an unmapped one is an unmapped host address on the device and the
        // source silently contributes nothing (#288).
        double *tr_es = GD->D.tracer_external_source;
        if (tr_es != NULL) {
            #pragma omp target enter data map(to: tr_es[0:ns*n])
        }

        // Boundary values are sized by boundary_length, which is zero on a
        // domain with no boundary edges -- mapping a zero-length array is not
        // meaningful, so follow the same nb > 0 guard the other bv arrays use.
        if (nb > 0) {
            double *tr_bv = GD->D.tracer_boundary_values;
            #pragma omp target enter data map(to: tr_bv[0:ns*nb])
        }
        // Deliberately NOT recording a 'tracer_arrays_mapped' flag in
        // struct gpu_domain: that struct embeds struct domain D, so adding a
        // member risks the silent offset-aliasing failure documented in
        // HANDOVER.md 2.1. number_of_tracers cannot change over the lifetime
        // of a mapped domain (add_tracer tears the GPU interface down), so the
        // unmap path re-tests the same condition instead.
    }

    // Map the per-class sediment parameter arrays (Phase 3). Tiny -- one
    // double per class -- but they are read inside a target region, so they
    // must be present on the device like everything else.
    if (GD->D.n_sediment_classes > 0) {
        anuga_int ncl = GD->D.n_sediment_classes;
        double *sed_vs = GD->D.sediment_settling_velocity;
        double *sed_ds = GD->D.sediment_d_star;
        double *sed_dm = GD->D.sediment_diameter;
        double *sed_rr = GD->D.sediment_R;
        double *sed_tc = GD->D.sediment_tau_c_star;
        double *sed_ar = GD->D.sediment_reference_height;
        double *sed_qx = GD->D.sediment_qbx;
        double *sed_qy = GD->D.sediment_qby;
        double *sed_qa = GD->D.sediment_qba;
        #pragma omp target enter data map(to: sed_vs[0:ncl], sed_ds[0:ncl], \
            sed_dm[0:ncl], sed_rr[0:ncl], sed_tc[0:ncl], sed_ar[0:ncl]) \
            map(alloc: sed_qx[0:n], sed_qy[0:n], sed_qa[0:n])

        // [L-5]. The scratch and the exhaustion snapshot are pure device
        // workspace, so alloc; the base itself is input and must be copied.
        // Mapped inside the same n_sediment_classes guard as everything
        // above, which is what allocates them.
        double *sed_sl = GD->D.sediment_source_limited;
        double *sed_sw = GD->D.sediment_slope_work;
        anuga_int *sed_ex = GD->D.sediment_bed_exhausted;
        double *sed_rdz = GD->D.sediment_repose_dz;
        // sed_sw is uploaded, not just allocated: with a frozen slope it
        // holds the bed slope computed on the host at setup.
        #pragma omp target enter data map(alloc: sed_sl[0:ncl*n], \
            sed_ex[0:n], sed_rdz[0:n]) map(to: sed_sw[0:n])
        if (GD->D.sediment_has_z_base) {
            double *sed_zb = GD->D.sediment_z_base;
            #pragma omp target enter data map(to: sed_zb[0:n])
        }
        // Bed composition: state the device owns from here on, so copied in
        // once and brought back by gpu_domain_sync_from_device.
        if (GD->D.sediment_bed_composition && GD->D.sediment_bed_active != NULL
                && GD->D.sediment_bed_substrate != NULL) {
            double *sed_ba = GD->D.sediment_bed_active;
            double *sed_bb = GD->D.sediment_bed_substrate;
            #pragma omp target enter data map(to: sed_ba[0:ncl*n], sed_bb[0:ncl*n])
        }
        // [T-16] the shear factor, input, only when set.
        double *sed_sf = GD->D.sediment_shear_factor;
        if (sed_sf != NULL) {
            #pragma omp target enter data map(to: sed_sf[0:n])
        }
        // Which boundary edges pass bedload: input, set once from Python.
        anuga_int *sed_bo = GD->D.sediment_bedload_open;
        if (sed_bo != NULL && nb > 0) {
            #pragma omp target enter data map(to: sed_bo[0:nb])
        }
        double *sed_bs = GD->D.sediment_bedload_supply;
        if (sed_bs != NULL && nb > 0) {
            #pragma omp target enter data map(to: sed_bs[0:nb])
        }
    }

    // Vegetation drag fields (spec 8): input, set once from Python.
    if (GD->D.vegetation_mode > 0) {
        double *veg_m = GD->D.veg_density_centroid_values;
        double *veg_d = GD->D.veg_diameter_centroid_values;
        double *veg_h = GD->D.veg_height_centroid_values;
        #pragma omp target enter data map(to: veg_m[0:n], veg_d[0:n], veg_h[0:n])
    }

    // Map halo exchange arrays if we have neighbors
    if (H->num_neighbors > 0) {
        int send_size = H->total_send_size;
        int recv_size = H->total_recv_size;
        int *flat_send = H->flat_send_indices;
        int *flat_recv = H->flat_recv_indices;
        double *send_buf = H->send_buffer;
        double *recv_buf = H->recv_buffer;
        // Must match gpu_halo_stride() in gpu_halo.c: 3 conserved hydro
        // quantities plus one m slot per tracer. A mismatch here maps a buffer
        // shorter than the exchange writes.
        const int halo_stride = 3 + (int)GD->D.number_of_tracers;

        #pragma omp target enter data map(to: flat_send[0:send_size], flat_recv[0:recv_size]) \
            map(alloc: send_buf[0:halo_stride*send_size], recv_buf[0:halo_stride*recv_size])
    }

    GD->gpu_initialized = 1;

    if (GD->rank == 0 && GD->verbose) {
        printf("  Arrays mapped to device %d\n", GD->device_id);
        fflush(stdout);
    }
    return 1;
}

void gpu_remap_boundary_arrays(struct gpu_domain *GD) {
    // Remap boundary arrays that were initialized after the initial map_to_gpu call.
    // This is needed when set_boundary() is called AFTER set_multiprocessor_mode().

    // Map reflective boundary arrays if initialized but not yet mapped
    struct reflective_boundary *R = &GD->reflective;
    if (R->num_edges > 0 && R->boundary_indices != NULL && !R->mapped) {
        int ne = R->num_edges;
        int *b_idx = R->boundary_indices;
        int *v_ids = R->vol_ids;
        int *e_ids = R->edge_ids;

        #pragma omp target enter data map(to: b_idx[0:ne], v_ids[0:ne], e_ids[0:ne])
        R->mapped = 1;

        if (GD->rank == 0 && GD->verbose) {
            printf("Reflective boundary arrays mapped to GPU (late): %d edges\n", ne);
            fflush(stdout);
        }
    }

    // Map dirichlet boundary arrays if initialized but not yet mapped
    struct dirichlet_boundary *Dir = &GD->dirichlet;
    if (Dir->num_edges > 0 && Dir->boundary_indices != NULL && !Dir->mapped) {
        int ne = Dir->num_edges;
        int *b_idx = Dir->boundary_indices;
        int *v_ids = Dir->vol_ids;
        int *e_ids = Dir->edge_ids;

        #pragma omp target enter data map(to: b_idx[0:ne], v_ids[0:ne], e_ids[0:ne])
        Dir->mapped = 1;

        if (GD->rank == 0 && GD->verbose) {
            printf("Dirichlet boundary arrays mapped to GPU (late): %d edges\n", ne);
            fflush(stdout);
        }
    }

    // Map transmissive boundary arrays if initialized but not yet mapped
    struct transmissive_boundary *T = &GD->transmissive;
    if (T->num_edges > 0 && T->boundary_indices != NULL && !T->mapped) {
        int ne = T->num_edges;
        int *b_idx = T->boundary_indices;
        int *v_ids = T->vol_ids;
        int *e_ids = T->edge_ids;

        #pragma omp target enter data map(to: b_idx[0:ne], v_ids[0:ne], e_ids[0:ne])
        T->mapped = 1;

        if (GD->rank == 0 && GD->verbose) {
            printf("Transmissive boundary arrays mapped to GPU (late): %d edges\n", ne);
            fflush(stdout);
        }
    }

    // Map transmissive_n_zero_t boundary arrays if initialized but not yet mapped
    struct transmissive_n_zero_t_boundary *Tnzt = &GD->transmissive_n_zero_t;
    if (Tnzt->num_edges > 0 && Tnzt->boundary_indices != NULL && !Tnzt->mapped) {
        int ne = Tnzt->num_edges;
        int *b_idx = Tnzt->boundary_indices;
        int *v_ids = Tnzt->vol_ids;
        int *e_ids = Tnzt->edge_ids;

        #pragma omp target enter data map(to: b_idx[0:ne], v_ids[0:ne], e_ids[0:ne])
        Tnzt->mapped = 1;

        if (GD->rank == 0 && GD->verbose) {
            printf("Transmissive_n_zero_t boundary arrays mapped to GPU (late): %d edges\n", ne);
            fflush(stdout);
        }
    }

    // Map time_boundary arrays if initialized but not yet mapped
    struct time_boundary *TB = &GD->time_bdry;
    if (TB->num_edges > 0 && TB->boundary_indices != NULL && !TB->mapped) {
        int ne = TB->num_edges;
        int *b_idx = TB->boundary_indices;
        int *v_ids = TB->vol_ids;
        int *e_ids = TB->edge_ids;
        double *s_val = TB->stage_values;
        double *x_val = TB->xmom_values;
        double *y_val = TB->ymom_values;

        #pragma omp target enter data map(to: b_idx[0:ne], v_ids[0:ne], e_ids[0:ne], \
                                              s_val[0:ne], x_val[0:ne], y_val[0:ne])
        TB->mapped = 1;

        if (GD->rank == 0 && GD->verbose) {
            printf("Time_boundary arrays mapped to GPU (late): %d edges\n", ne);
            fflush(stdout);
        }
    }

    // Map absorbing_wave boundary arrays if initialized but not yet mapped
    struct absorbing_wave_boundary *AW = &GD->absorbing_wave;
    if (AW->num_edges > 0 && AW->boundary_indices != NULL && !AW->mapped) {
        int ne = AW->num_edges;
        int *b_idx = AW->boundary_indices;
        int *v_ids = AW->vol_ids;
        int *e_ids = AW->edge_ids;

        #pragma omp target enter data map(to: b_idx[0:ne], v_ids[0:ne], e_ids[0:ne])
        AW->mapped = 1;

        if (GD->rank == 0 && GD->verbose) {
            printf("Absorbing_wave boundary arrays mapped to GPU (late): %d edges\n", ne);
            fflush(stdout);
        }
    }

    // Map characteristic_wave boundary arrays if initialized but not yet mapped
    struct characteristic_wave_boundary *CW = &GD->characteristic_wave;
    if (CW->num_edges > 0 && CW->boundary_indices != NULL && !CW->mapped) {
        int ne = CW->num_edges;
        int *b_idx = CW->boundary_indices;
        int *v_ids = CW->vol_ids;
        int *e_ids = CW->edge_ids;

        #pragma omp target enter data map(to: b_idx[0:ne], v_ids[0:ne], e_ids[0:ne])
        CW->mapped = 1;

        if (GD->rank == 0 && GD->verbose) {
            printf("Characteristic_wave boundary arrays mapped to GPU (late): %d edges\n", ne);
            fflush(stdout);
        }
    }

    // Map flather boundary arrays if initialized but not yet mapped
    struct flather_boundary *FL = &GD->flather;
    if (FL->num_edges > 0 && FL->boundary_indices != NULL && !FL->mapped) {
        int ne = FL->num_edges;
        int *b_idx = FL->boundary_indices;
        int *v_ids = FL->vol_ids;
        int *e_ids = FL->edge_ids;

        #pragma omp target enter data map(to: b_idx[0:ne], v_ids[0:ne], e_ids[0:ne])
        FL->mapped = 1;

        if (GD->rank == 0 && GD->verbose) {
            printf("Flather boundary arrays mapped to GPU (late): %d edges\n", ne);
            fflush(stdout);
        }
    }
}

void gpu_domain_unmap_arrays(struct gpu_domain *GD) {
    if (!GD->gpu_initialized) return;

    anuga_int n = GD->D.number_of_elements;
    anuga_int nb = GD->D.boundary_length;
    struct halo_exchange *H = &GD->halo;

    // Extract pointers to local variables
    double *stage_cv = GD->D.stage_centroid_values;
    double *xmom_cv = GD->D.xmom_centroid_values;
    double *ymom_cv = GD->D.ymom_centroid_values;
    double *bed_cv = GD->D.bed_centroid_values;
    double *height_cv = GD->D.height_centroid_values;
    double *stage_ev = GD->D.stage_edge_values;
    double *xmom_ev = GD->D.xmom_edge_values;
    double *ymom_ev = GD->D.ymom_edge_values;
    double *bed_ev = GD->D.bed_edge_values;
    double *height_ev = GD->D.height_edge_values;
    double *stage_eu = GD->D.stage_explicit_update;
    double *xmom_eu = GD->D.xmom_explicit_update;
    double *ymom_eu = GD->D.ymom_explicit_update;
    double *stage_siu = GD->D.stage_semi_implicit_update;
    double *xmom_siu = GD->D.xmom_semi_implicit_update;
    double *ymom_siu = GD->D.ymom_semi_implicit_update;
    anuga_int *neighbours = GD->D.neighbours;
    anuga_int *neighbour_edges = GD->D.neighbour_edges;
    anuga_geom_t *normals = GD->D.normals;
    anuga_geom_t *edgelengths = GD->D.edgelengths;
    anuga_geom_t *areas = GD->D.areas;
    anuga_geom_t *radii = GD->D.radii;
    double *max_speed = GD->D.max_speed;
    anuga_geom_t *centroid_coords = GD->D.centroid_coordinates;
    anuga_geom_t *edge_coords = GD->D.edge_coordinates;

    // Additional arrays for extrapolation
    anuga_int *surrogate_neighbours = GD->D.surrogate_neighbours;
    anuga_int *number_of_boundaries = GD->D.number_of_boundaries;
    double *x_centroid_work = GD->D.x_centroid_work;
    double *y_centroid_work = GD->D.y_centroid_work;

    // Friction array
    double *friction_cv = GD->D.friction_centroid_values;

    // tri_full_flag for MPI ghost cell identification
    anuga_int *tri_full_flag = GD->D.tri_full_flag;

    // Unmap domain arrays
    #pragma omp target exit data map(delete: \
        stage_cv[0:n], xmom_cv[0:n], ymom_cv[0:n], \
        bed_cv[0:n], height_cv[0:n], friction_cv[0:n], \
        stage_ev[0:3*n], xmom_ev[0:3*n], ymom_ev[0:3*n], \
        bed_ev[0:3*n], height_ev[0:3*n], \
        stage_eu[0:n], xmom_eu[0:n], ymom_eu[0:n], \
        stage_siu[0:n], xmom_siu[0:n], ymom_siu[0:n], \
        neighbours[0:3*n], neighbour_edges[0:3*n], \
        surrogate_neighbours[0:3*n], number_of_boundaries[0:n], \
        x_centroid_work[0:n], y_centroid_work[0:n], \
        normals[0:6*n], edgelengths[0:3*n], \
        areas[0:n], radii[0:n], max_speed[0:n], \
        centroid_coords[0:2*n], edge_coords[0:6*n])

    // Unmap tri_full_flag only when it was mapped (parallel domains only)
    if (tri_full_flag != NULL) {
        #pragma omp target exit data map(delete: tri_full_flag[0:n])
    }

    if (nb > 0) {
        double *stage_bv = GD->D.stage_boundary_values;
        double *xmom_bv = GD->D.xmom_boundary_values;
        double *ymom_bv = GD->D.ymom_boundary_values;
        double *bed_bv = GD->D.bed_boundary_values;
        double *height_bv = GD->D.height_boundary_values;
        #pragma omp target exit data map(delete: \
            stage_bv[0:nb], xmom_bv[0:nb], ymom_bv[0:nb], \
            bed_bv[0:nb], height_bv[0:nb])
    }

    // Unmap reflective boundary arrays if mapped
    struct reflective_boundary *R = &GD->reflective;
    if (R->mapped && R->num_edges > 0) {
        int ne = R->num_edges;
        int *b_idx = R->boundary_indices;
        int *v_ids = R->vol_ids;
        int *e_ids = R->edge_ids;
        #pragma omp target exit data map(delete: b_idx[0:ne], v_ids[0:ne], e_ids[0:ne])
        R->mapped = 0;
    }

    // Unmap dirichlet boundary arrays if mapped
    struct dirichlet_boundary *Dir = &GD->dirichlet;
    if (Dir->mapped && Dir->num_edges > 0) {
        int ne = Dir->num_edges;
        int *b_idx = Dir->boundary_indices;
        int *v_ids = Dir->vol_ids;
        int *e_ids = Dir->edge_ids;
        #pragma omp target exit data map(delete: b_idx[0:ne], v_ids[0:ne], e_ids[0:ne])
        Dir->mapped = 0;
    }

    // Unmap transmissive boundary arrays if mapped
    struct transmissive_boundary *T = &GD->transmissive;
    if (T->mapped && T->num_edges > 0) {
        int ne = T->num_edges;
        int *b_idx = T->boundary_indices;
        int *v_ids = T->vol_ids;
        int *e_ids = T->edge_ids;
        #pragma omp target exit data map(delete: b_idx[0:ne], v_ids[0:ne], e_ids[0:ne])
        T->mapped = 0;
    }

    // Unmap transmissive_n_zero_t boundary arrays if mapped
    struct transmissive_n_zero_t_boundary *Tnzt = &GD->transmissive_n_zero_t;
    if (Tnzt->mapped && Tnzt->num_edges > 0) {
        int ne = Tnzt->num_edges;
        int *b_idx = Tnzt->boundary_indices;
        int *v_ids = Tnzt->vol_ids;
        int *e_ids = Tnzt->edge_ids;
        #pragma omp target exit data map(delete: b_idx[0:ne], v_ids[0:ne], e_ids[0:ne])
        Tnzt->mapped = 0;
    }

    // Unmap time_boundary arrays if mapped
    struct time_boundary *TB = &GD->time_bdry;
    if (TB->mapped && TB->num_edges > 0) {
        int ne = TB->num_edges;
        int *b_idx = TB->boundary_indices;
        int *v_ids = TB->vol_ids;
        int *e_ids = TB->edge_ids;
        double *s_val = TB->stage_values;
        double *x_val = TB->xmom_values;
        double *y_val = TB->ymom_values;
        #pragma omp target exit data map(delete: b_idx[0:ne], v_ids[0:ne], e_ids[0:ne], \
                                                 s_val[0:ne], x_val[0:ne], y_val[0:ne])
        TB->mapped = 0;
    }

    // Unmap file_boundary arrays if mapped
    struct file_boundary *FB = &GD->file_bdry;
    if (FB->mapped && FB->num_edges > 0) {
        int ne = FB->num_edges;
        int    *b_idx   = FB->boundary_indices;
        int    *v_ids   = FB->vol_ids;
        int    *e_ids   = FB->edge_ids;
        double *stage_v = FB->stage_values;
        double *xmom_v  = FB->xmom_values;
        double *ymom_v  = FB->ymom_values;
        #pragma omp target exit data map(delete: b_idx[0:ne], v_ids[0:ne], e_ids[0:ne], \
                                                 stage_v[0:ne], xmom_v[0:ne], ymom_v[0:ne])
        FB->mapped = 0;
    }

    // Unmap absorbing_wave boundary arrays if mapped
    struct absorbing_wave_boundary *AW = &GD->absorbing_wave;
    if (AW->mapped && AW->num_edges > 0) {
        int ne = AW->num_edges;
        int *b_idx = AW->boundary_indices;
        int *v_ids = AW->vol_ids;
        int *e_ids = AW->edge_ids;
        #pragma omp target exit data map(delete: b_idx[0:ne], v_ids[0:ne], e_ids[0:ne])
        AW->mapped = 0;
    }

    // Unmap characteristic_wave boundary arrays if mapped
    struct characteristic_wave_boundary *CW = &GD->characteristic_wave;
    if (CW->mapped && CW->num_edges > 0) {
        int ne = CW->num_edges;
        int *b_idx = CW->boundary_indices;
        int *v_ids = CW->vol_ids;
        int *e_ids = CW->edge_ids;
        #pragma omp target exit data map(delete: b_idx[0:ne], v_ids[0:ne], e_ids[0:ne])
        CW->mapped = 0;
    }

    // Unmap flather boundary arrays if mapped
    struct flather_boundary *FL = &GD->flather;
    if (FL->mapped && FL->num_edges > 0) {
        int ne = FL->num_edges;
        int *b_idx = FL->boundary_indices;
        int *v_ids = FL->vol_ids;
        int *e_ids = FL->edge_ids;
        #pragma omp target exit data map(delete: b_idx[0:ne], v_ids[0:ne], e_ids[0:ne])
        FL->mapped = 0;
    }

    // Unmap riverwall arrays
    anuga_int *edge_flux_type = GD->D.edge_flux_type;
    if (edge_flux_type != NULL) {
        #pragma omp target exit data map(delete: edge_flux_type[0:3*n])
    }

    anuga_int n_rw_edges = GD->D.number_of_riverwall_edges;
    if (n_rw_edges > 0) {
        anuga_int *edge_river_wall_counter = GD->D.edge_river_wall_counter;
        double *riverwall_elevation = GD->D.riverwall_elevation;
        anuga_int *riverwall_rowIndex = GD->D.riverwall_rowIndex;
        double *riverwall_hydraulic_properties = GD->D.riverwall_hydraulic_properties;
        anuga_int ncol_hp = GD->D.ncol_riverwall_hydraulic_properties;
        anuga_int nrow_hp = GD->D.nrow_riverwall_hydraulic_properties;

        if (edge_river_wall_counter != NULL) {
            #pragma omp target exit data map(delete: edge_river_wall_counter[0:3*n])
        }
        if (riverwall_elevation != NULL) {
            #pragma omp target exit data map(delete: riverwall_elevation[0:n_rw_edges])
        }
        if (riverwall_rowIndex != NULL) {
            #pragma omp target exit data map(delete: riverwall_rowIndex[0:n_rw_edges])
        }
        if (riverwall_hydraulic_properties != NULL && ncol_hp > 0 && nrow_hp > 0) {
            #pragma omp target exit data map(delete: riverwall_hydraulic_properties[0:nrow_hp*ncol_hp])
        }
    }

    if (GD->backup_arrays_mapped) {
        double *stage_backup = GD->D.stage_backup_values;
        double *xmom_backup = GD->D.xmom_backup_values;
        double *ymom_backup = GD->D.ymom_backup_values;
        #pragma omp target exit data map(delete: \
            stage_backup[0:n], xmom_backup[0:n], ymom_backup[0:n])
    }

    // Unmap tracer arrays. Same condition as the map path above -- see the note
    // there on why this is re-tested rather than recorded in a struct flag.
    if (GD->D.number_of_tracers > 0) {
        anuga_int ns = GD->D.number_of_tracers;
        double *tr_cv = GD->D.tracer_centroid_values;
        double *tr_ev = GD->D.tracer_edge_values;
        double *tr_eu = GD->D.tracer_explicit_update;
        double *tr_qv = GD->D.tracer_conserved_values;
        double *tr_bk = GD->D.tracer_backup_values;
        double *tr_bf = GD->D.tracer_boundary_flux;
        #pragma omp target exit data map(delete: \
            tr_cv[0:ns*n], tr_ev[0:ns*3*n], \
            tr_eu[0:ns*n], tr_qv[0:ns*n], tr_bk[0:ns*n], \
            tr_bf[0:ns*n])
        double *tr_sf = GD->D.tracer_speed_factor;
        if (tr_sf != NULL) {
            #pragma omp target exit data map(delete: tr_sf[0:ns*n])
        }
        double *tr_es = GD->D.tracer_external_source;
        if (tr_es != NULL) {
            #pragma omp target exit data map(delete: tr_es[0:ns*n])
        }
        if (nb > 0) {
            double *tr_bv = GD->D.tracer_boundary_values;
            #pragma omp target exit data map(delete: tr_bv[0:ns*nb])
        }
    }

    // Unmap the sediment parameter arrays. Same condition as the map path.
    if (GD->D.n_sediment_classes > 0) {
        anuga_int ncl = GD->D.n_sediment_classes;
        double *sed_vs = GD->D.sediment_settling_velocity;
        double *sed_ds = GD->D.sediment_d_star;
        double *sed_dm = GD->D.sediment_diameter;
        double *sed_rr = GD->D.sediment_R;
        double *sed_tc = GD->D.sediment_tau_c_star;
        double *sed_ar = GD->D.sediment_reference_height;
        double *sed_qx = GD->D.sediment_qbx;
        double *sed_qy = GD->D.sediment_qby;
        double *sed_qa = GD->D.sediment_qba;
        #pragma omp target exit data map(delete: sed_vs[0:ncl], sed_ds[0:ncl], \
            sed_dm[0:ncl], sed_rr[0:ncl], sed_tc[0:ncl], sed_ar[0:ncl], \
            sed_qx[0:n], sed_qy[0:n], sed_qa[0:n])

        // [L-5]. Mirrors the enter-data above, on the same guard.
        double *sed_sl = GD->D.sediment_source_limited;
        double *sed_sw = GD->D.sediment_slope_work;
        anuga_int *sed_ex = GD->D.sediment_bed_exhausted;
        double *sed_rdz = GD->D.sediment_repose_dz;
        #pragma omp target exit data map(delete: sed_sl[0:ncl*n], sed_sw[0:n], \
            sed_ex[0:n], sed_rdz[0:n])
        if (GD->D.sediment_has_z_base) {
            double *sed_zb = GD->D.sediment_z_base;
            #pragma omp target exit data map(delete: sed_zb[0:n])
        }
        if (GD->D.sediment_bed_composition && GD->D.sediment_bed_active != NULL
                && GD->D.sediment_bed_substrate != NULL) {
            double *sed_ba = GD->D.sediment_bed_active;
            double *sed_bb = GD->D.sediment_bed_substrate;
            #pragma omp target exit data map(delete: sed_ba[0:ncl*n], sed_bb[0:ncl*n])
        }
        double *sed_sf = GD->D.sediment_shear_factor;
        if (sed_sf != NULL) {
            #pragma omp target exit data map(delete: sed_sf[0:n])
        }
        anuga_int *sed_bo = GD->D.sediment_bedload_open;
        if (sed_bo != NULL && nb > 0) {
            #pragma omp target exit data map(delete: sed_bo[0:nb])
        }
        double *sed_bs = GD->D.sediment_bedload_supply;
        if (sed_bs != NULL && nb > 0) {
            #pragma omp target exit data map(delete: sed_bs[0:nb])
        }
    }

    if (GD->D.vegetation_mode > 0) {
        double *veg_m = GD->D.veg_density_centroid_values;
        double *veg_d = GD->D.veg_diameter_centroid_values;
        double *veg_h = GD->D.veg_height_centroid_values;
        #pragma omp target exit data map(delete: veg_m[0:n], veg_d[0:n], veg_h[0:n])
    }

    // Unmap halo arrays
    if (H->num_neighbors > 0) {
        int send_size = H->total_send_size;
        int recv_size = H->total_recv_size;
        int *flat_send = H->flat_send_indices;
        int *flat_recv = H->flat_recv_indices;
        double *send_buf = H->send_buffer;
        double *recv_buf = H->recv_buffer;
        // Must match gpu_halo_stride() in gpu_halo.c: 3 conserved hydro
        // quantities plus one m slot per tracer. A mismatch here maps a buffer
        // shorter than the exchange writes.
        const int halo_stride = 3 + (int)GD->D.number_of_tracers;

        #pragma omp target exit data map(delete: \
            flat_send[0:send_size], flat_recv[0:recv_size], \
            send_buf[0:halo_stride*send_size], recv_buf[0:halo_stride*recv_size])
    }

    GD->gpu_initialized = 0;
}

void gpu_domain_sync_to_device(struct gpu_domain *GD) {
    // Sync all centroid values to GPU (use at start of GPU computation)
    if (!GD->gpu_initialized) return;

    anuga_int n = GD->D.number_of_elements;
    double *stage_cv = GD->D.stage_centroid_values;
    double *xmom_cv = GD->D.xmom_centroid_values;
    double *ymom_cv = GD->D.ymom_centroid_values;
    double *height_cv = GD->D.height_centroid_values;

    #pragma omp target update to(stage_cv[0:n], xmom_cv[0:n], ymom_cv[0:n], height_cv[0:n])

    // Tracers: the conserved variable m = h*c is what the device integrates,
    // and c is what Python seeds via set_tracer, so both have to go across.
    // Boundary values too -- they are set on the host by boundary conditions.
    if (GD->D.number_of_tracers > 0) {
        anuga_int ns = GD->D.number_of_tracers;
        anuga_int nb = GD->D.boundary_length;
        double *tr_cv = GD->D.tracer_centroid_values;
        double *tr_qv = GD->D.tracer_conserved_values;
        #pragma omp target update to(tr_cv[0:ns*n], tr_qv[0:ns*n])
        if (nb > 0) {
            double *tr_bv = GD->D.tracer_boundary_values;
            #pragma omp target update to(tr_bv[0:ns*nb])
        }
    }
}

// Push the external source [G-3] to the device on its own.
//
// set_tracer_source writes it on the HOST, and a time-varying supply -- a
// manufactured solution, a tributary hydrograph -- rewrites it every step. The
// device copy would otherwise hold whatever it had when the arrays were mapped
// and the source would be silently stale rather than absent, which is worse.
// Separate from gpu_sync_to_device so a per-step push costs one array, not all
// of them. #288
void gpu_sync_tracer_source_to_device(struct gpu_domain *GD)
{
    if (!GD->gpu_initialized) return;
    if (GD->D.number_of_tracers <= 0) return;

    double *tr_es = GD->D.tracer_external_source;
    if (tr_es == NULL) return;

    anuga_int n = GD->D.number_of_elements;
    anuga_int ns = GD->D.number_of_tracers;
    #pragma omp target update to(tr_es[0:ns*n])
}

// Push the tracer boundary concentrations to the device on their own.
//
// set_tracer_boundary() and update_tracer_boundary_values() write the HOST
// array.  The device copy is made once, when the arrays are mapped, and
// gpu_sync_boundary_values() (the per-step push the boundary evaluators use)
// carries only the hydrodynamic boundary values -- so without this a value
// set after the interface exists, or re-evaluated for a later time, is never
// seen by the flux kernel and every inflow edge keeps injecting the mapped
// (t = 0, or zero-filled) concentration.  One small array, like the source.
void gpu_sync_tracer_boundary_to_device(struct gpu_domain *GD)
{
    if (!GD->gpu_initialized) return;
    if (GD->D.number_of_tracers <= 0) return;

    anuga_int nb = GD->D.boundary_length;
    if (nb <= 0) return;   /* not mapped either; see gpu_domain_map_arrays */

    double *tr_bv = GD->D.tracer_boundary_values;
    if (tr_bv == NULL) return;

    anuga_int ns = GD->D.number_of_tracers;
    #pragma omp target update to(tr_bv[0:ns*nb])
}

void gpu_domain_sync_from_device(struct gpu_domain *GD) {
    // Sync centroid values from GPU (use at yieldstep for Python I/O)
    if (!GD->gpu_initialized) return;

    anuga_int n = GD->D.number_of_elements;
    double *stage_cv = GD->D.stage_centroid_values;
    double *xmom_cv = GD->D.xmom_centroid_values;
    double *ymom_cv = GD->D.ymom_centroid_values;
    double *height_cv = GD->D.height_centroid_values;

    #pragma omp target update from(stage_cv[0:n], xmom_cv[0:n], ymom_cv[0:n], height_cv[0:n])

    // Tracers: bring back both the conserved m and the derived c. Python reads
    // c via get_tracer and computes mass from m, so returning only one of them
    // would hand back a domain whose c and m disagree.
    if (GD->D.number_of_tracers > 0) {
        anuga_int ns = GD->D.number_of_tracers;
        double *tr_cv = GD->D.tracer_centroid_values;
        double *tr_qv = GD->D.tracer_conserved_values;
        #pragma omp target update from(tr_cv[0:ns*n], tr_qv[0:ns*n])
    }

    // Bed elevation, once sediment can move it. Normally z is constant and
    // never needs syncing back, which is why it is not in the list above; with
    // the Exner update of [G-4] the device owns it and the host copy would
    // otherwise stay at its initial value for ever. Centroid AND edge, because
    // the DE algorithms use discontinuous elevation and the kernel shifts both.
    if (GD->D.n_sediment_classes > 0) {
        double *bed_cv_s = GD->D.bed_centroid_values;
        double *bed_ev_s = GD->D.bed_edge_values;
        #pragma omp target update from(bed_cv_s[0:n], bed_ev_s[0:3*n])
        if (GD->D.sediment_bed_composition && GD->D.sediment_bed_active != NULL
                && GD->D.sediment_bed_substrate != NULL) {
            anuga_int ncl = GD->D.n_sediment_classes;
            double *sed_ba = GD->D.sediment_bed_active;
            double *sed_bb = GD->D.sediment_bed_substrate;
            #pragma omp target update from(sed_ba[0:ncl*n], sed_bb[0:ncl*n])
        }
    }
}

void gpu_domain_sync_all_from_device(struct gpu_domain *GD) {
    // Sync ALL arrays from GPU (for debugging/testing intermediate values)
    if (!GD->gpu_initialized) return;

    anuga_int n = GD->D.number_of_elements;
    anuga_int nb = GD->D.boundary_length;

    // Centroid values
    double *stage_cv = GD->D.stage_centroid_values;
    double *xmom_cv = GD->D.xmom_centroid_values;
    double *ymom_cv = GD->D.ymom_centroid_values;
    double *height_cv = GD->D.height_centroid_values;

    // Edge values
    double *stage_ev = GD->D.stage_edge_values;
    double *xmom_ev = GD->D.xmom_edge_values;
    double *ymom_ev = GD->D.ymom_edge_values;
    double *height_ev = GD->D.height_edge_values;
    double *bed_ev = GD->D.bed_edge_values;

    // Explicit and semi-implicit updates
    double *stage_eu = GD->D.stage_explicit_update;
    double *xmom_eu = GD->D.xmom_explicit_update;
    double *ymom_eu = GD->D.ymom_explicit_update;
    double *stage_siu = GD->D.stage_semi_implicit_update;
    double *xmom_siu = GD->D.xmom_semi_implicit_update;
    double *ymom_siu = GD->D.ymom_semi_implicit_update;

    // Sync centroid values
    #pragma omp target update from(stage_cv[0:n], xmom_cv[0:n], ymom_cv[0:n], height_cv[0:n])

    // Sync edge values
    #pragma omp target update from(stage_ev[0:3*n], xmom_ev[0:3*n], ymom_ev[0:3*n], \
                                   height_ev[0:3*n], bed_ev[0:3*n])

    // Sync explicit and semi-implicit updates
    #pragma omp target update from(stage_eu[0:n], xmom_eu[0:n], ymom_eu[0:n], \
                                   stage_siu[0:n], xmom_siu[0:n], ymom_siu[0:n])

    // Sync boundary values if present
    if (nb > 0) {
        double *stage_bv = GD->D.stage_boundary_values;
        double *xmom_bv = GD->D.xmom_boundary_values;
        double *ymom_bv = GD->D.ymom_boundary_values;
        double *height_bv = GD->D.height_boundary_values;
        double *bed_bv = GD->D.bed_boundary_values;

        #pragma omp target update from(stage_bv[0:nb], xmom_bv[0:nb], ymom_bv[0:nb], \
                                       height_bv[0:nb], bed_bv[0:nb])
    }

    // Sync backup values if mapped
    if (GD->backup_arrays_mapped) {
        double *stage_backup = GD->D.stage_backup_values;
        double *xmom_backup = GD->D.xmom_backup_values;
        double *ymom_backup = GD->D.ymom_backup_values;
        #pragma omp target update from(stage_backup[0:n], xmom_backup[0:n], ymom_backup[0:n])
    }
}

void gpu_sync_boundary_values(struct gpu_domain *GD) {
    // Sync boundary values TO GPU (after CPU boundary evaluation)
    if (!GD->gpu_initialized) return;

    anuga_int nb = GD->D.boundary_length;
    if (nb == 0) return;

    double *stage_bv = GD->D.stage_boundary_values;
    double *xmom_bv = GD->D.xmom_boundary_values;
    double *ymom_bv = GD->D.ymom_boundary_values;
    double *bed_bv = GD->D.bed_boundary_values;
    double *height_bv = GD->D.height_boundary_values;

    #pragma omp target update to(stage_bv[0:nb], xmom_bv[0:nb], ymom_bv[0:nb], \
                                 bed_bv[0:nb], height_bv[0:nb])
}

void gpu_sync_riverwall_to_device(struct gpu_domain *GD) {
    // Sync riverwall crest elevations and hydraulic properties TO GPU.
    //
    // These arrays are mapped once at setup and are never written on the device,
    // so a host-side change — RiverWall.set_elevation() / set_elevation_offset()
    // / set_hydraulic_parameter(), e.g. operating a gate mid-run — is invisible
    // to the kernels until it is pushed back across. Cheap: sized by the number
    // of riverwall edges, not by the mesh.
    if (!GD->gpu_initialized) return;

    anuga_int n_rw_edges = GD->D.number_of_riverwall_edges;
    if (n_rw_edges <= 0) return;

    double *riverwall_elevation = GD->D.riverwall_elevation;
    if (riverwall_elevation != NULL) {
        #pragma omp target update to(riverwall_elevation[0:n_rw_edges])
    }

    double *riverwall_hydraulic_properties = GD->D.riverwall_hydraulic_properties;
    anuga_int ncol_hp = GD->D.ncol_riverwall_hydraulic_properties;
    anuga_int nrow_hp = GD->D.nrow_riverwall_hydraulic_properties;
    if (riverwall_hydraulic_properties != NULL && ncol_hp > 0 && nrow_hp > 0) {
        #pragma omp target update to(riverwall_hydraulic_properties[0:nrow_hp*ncol_hp])
    }
}

void gpu_sync_edge_values_from_device(struct gpu_domain *GD) {
    // Sync ALL edge values FROM GPU - expensive, use sparse version if possible
    if (!GD->gpu_initialized) return;

    anuga_int n = GD->D.number_of_elements;

    double *stage_ev = GD->D.stage_edge_values;
    double *xmom_ev = GD->D.xmom_edge_values;
    double *ymom_ev = GD->D.ymom_edge_values;
    double *bed_ev = GD->D.bed_edge_values;
    double *height_ev = GD->D.height_edge_values;

    #pragma omp target update from(stage_ev[0:3*n], xmom_ev[0:3*n], ymom_ev[0:3*n], \
                                   bed_ev[0:3*n], height_ev[0:3*n])
}

// Initialize boundary edge sync buffers - call once after boundaries are set
int gpu_boundary_edge_sync_init(struct gpu_domain *GD,
                                int num_boundary_cells,
                                int *boundary_cell_ids) {
    if (!GD->gpu_initialized) return -1;

    struct boundary_edge_sync *S = &GD->edge_sync;

    if (S->initialized) {
        // Already initialized - clean up first
        gpu_boundary_edge_sync_finalize(GD);
    }

    S->num_boundary_cells = num_boundary_cells;
    S->buf_size = num_boundary_cells * 3;  // 3 edges per cell

    if (num_boundary_cells == 0) {
        S->initialized = 1;
        return 0;  // Nothing to do
    }

    // Allocate cell IDs array (copy from Python's array) and staging buffers
    S->cell_ids = (int*)malloc(num_boundary_cells * sizeof(int));
    S->stage_buf = (double*)malloc(S->buf_size * sizeof(double));
    S->xmom_buf = (double*)malloc(S->buf_size * sizeof(double));
    S->ymom_buf = (double*)malloc(S->buf_size * sizeof(double));
    S->bed_buf = (double*)malloc(S->buf_size * sizeof(double));
    S->height_buf = (double*)malloc(S->buf_size * sizeof(double));
    if (!S->cell_ids || !S->stage_buf || !S->xmom_buf || !S->ymom_buf ||
        !S->bed_buf || !S->height_buf) {
        gpu_set_error(GD, "boundary edge sync: could not allocate staging for %d cells", num_boundary_cells);
        free(S->cell_ids);  free(S->stage_buf); free(S->xmom_buf);
        free(S->ymom_buf);  free(S->bed_buf);   free(S->height_buf);
        S->cell_ids = NULL; S->stage_buf = S->xmom_buf = S->ymom_buf = S->bed_buf = S->height_buf = NULL;
        S->num_boundary_cells = 0;
        S->buf_size = 0;
        return -1;
    }
    memcpy(S->cell_ids, boundary_cell_ids, num_boundary_cells * sizeof(int));

    // Map all buffers to GPU once
    int nc = num_boundary_cells;
    int bs = S->buf_size;
    int *cell_ids_ptr = S->cell_ids;
    double *stage_buf = S->stage_buf;
    double *xmom_buf = S->xmom_buf;
    double *ymom_buf = S->ymom_buf;
    double *bed_buf = S->bed_buf;
    double *height_buf = S->height_buf;

    #pragma omp target enter data map(to: cell_ids_ptr[0:nc]) \
        map(alloc: stage_buf[0:bs], xmom_buf[0:bs], ymom_buf[0:bs], \
                   bed_buf[0:bs], height_buf[0:bs])

    S->initialized = 1;
    if (GD->verbose) {
        printf("Rank %d: Boundary edge sync initialized for %d cells (%d edge values)\n",
               GD->rank, num_boundary_cells, S->buf_size);
        fflush(stdout);
    }

    return 0;
}

// Finalize boundary edge sync buffers
void gpu_boundary_edge_sync_finalize(struct gpu_domain *GD) {
    struct boundary_edge_sync *S = &GD->edge_sync;

    if (!S->initialized) return;

    if (S->num_boundary_cells > 0) {
        // Unmap from GPU
        int nc = S->num_boundary_cells;
        int bs = S->buf_size;
        int *cell_ids_ptr = S->cell_ids;
        double *stage_buf = S->stage_buf;
        double *xmom_buf = S->xmom_buf;
        double *ymom_buf = S->ymom_buf;
        double *bed_buf = S->bed_buf;
        double *height_buf = S->height_buf;

        #pragma omp target exit data map(delete: cell_ids_ptr[0:nc], \
            stage_buf[0:bs], xmom_buf[0:bs], ymom_buf[0:bs], \
            bed_buf[0:bs], height_buf[0:bs])

        // Free host memory
        free(S->cell_ids);
        free(S->stage_buf);
        free(S->xmom_buf);
        free(S->ymom_buf);
        free(S->bed_buf);
        free(S->height_buf);
    }

    S->cell_ids = NULL;
    S->stage_buf = NULL;
    S->xmom_buf = NULL;
    S->ymom_buf = NULL;
    S->bed_buf = NULL;
    S->height_buf = NULL;
    S->num_boundary_cells = 0;
    S->buf_size = 0;
    S->initialized = 0;
}

// Sync boundary edge values from GPU - fast version using pre-allocated buffers
void gpu_boundary_edge_sync(struct gpu_domain *GD) {
    struct boundary_edge_sync *S = &GD->edge_sync;

    if (!S->initialized || S->num_boundary_cells == 0) return;

    int nc = S->num_boundary_cells;
    int bs = S->buf_size;
    int *cell_ids = S->cell_ids;
    double *stage_buf = S->stage_buf;
    double *xmom_buf = S->xmom_buf;
    double *ymom_buf = S->ymom_buf;
    double *bed_buf = S->bed_buf;
    double *height_buf = S->height_buf;

    double *stage_ev = GD->D.stage_edge_values;
    double *xmom_ev = GD->D.xmom_edge_values;
    double *ymom_ev = GD->D.ymom_edge_values;
    double *bed_ev = GD->D.bed_edge_values;
    double *height_ev = GD->D.height_edge_values;

    // Gather on GPU into contiguous staging buffers
    OMP_PARALLEL_LOOP
    for (int i = 0; i < nc; i++) {
        int vid = cell_ids[i];
        for (int e = 0; e < 3; e++) {
            int buf_idx = i * 3 + e;
            int src_idx = vid * 3 + e;
            stage_buf[buf_idx] = stage_ev[src_idx];
            xmom_buf[buf_idx] = xmom_ev[src_idx];
            ymom_buf[buf_idx] = ymom_ev[src_idx];
            bed_buf[buf_idx] = bed_ev[src_idx];
            height_buf[buf_idx] = height_ev[src_idx];
        }
    }

    // Sync staging buffers from GPU to host
    #pragma omp target update from(stage_buf[0:bs], xmom_buf[0:bs], ymom_buf[0:bs], \
                                   bed_buf[0:bs], height_buf[0:bs])

    // Scatter on CPU to the actual edge value arrays (host copies)
    for (int i = 0; i < nc; i++) {
        int vid = S->cell_ids[i];  // Use host copy of cell_ids
        for (int e = 0; e < 3; e++) {
            int buf_idx = i * 3 + e;
            int dst_idx = vid * 3 + e;
            stage_ev[dst_idx] = stage_buf[buf_idx];
            xmom_ev[dst_idx] = xmom_buf[buf_idx];
            ymom_ev[dst_idx] = ymom_buf[buf_idx];
            bed_ev[dst_idx] = bed_buf[buf_idx];
            height_ev[dst_idx] = height_buf[buf_idx];
        }
    }
}

