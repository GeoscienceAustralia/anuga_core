// GPU-accelerated domain struct with MPI halo exchange support
//
// This extends the base ANUGA domain struct with:
// - MPI communicator and process info
// - Flattened halo exchange structures for GPU efficiency
// - GPU state tracking
//
// Based on miniapp_mpi.c patterns

#ifndef GPU_DOMAIN_H
#define GPU_DOMAIN_H

#define GPU_LAST_ERROR_LEN 512

#include <stdint.h>
#include <stdio.h>
#include <stdbool.h>
#ifdef HAVE_MPI
#include <mpi.h>
#else
#include "gpu_mpi_stubs.h"
#endif
#include "anuga_typedefs.h"
#include "sw_domain.h"

// Halo exchange structure - flattened for GPU efficiency
// Mirrors the pattern from miniapp_mpi.c halo_info struct
struct halo_exchange {
    int num_neighbors;           // Number of MPI neighbors we communicate with
    int *neighbor_ranks;         // Ranks we communicate with [num_neighbors]
    int *send_counts;            // Elements to send to each neighbor [num_neighbors]
    int *recv_counts;            // Elements to receive from each neighbor [num_neighbors]

    // Total sizes for buffer allocation
    int total_send_size;         // Sum of send_counts
    int total_recv_size;         // Sum of recv_counts

    // Flattened index arrays for GPU (contiguous memory)
    // These replace the per-neighbor arrays from Python dicts
    int *flat_send_indices;      // Local element indices to pack for sending [total_send_size]
    int *flat_recv_indices;      // Local element indices where received data goes [total_recv_size]

    // Offset arrays for navigating flattened indices per neighbor
    int *send_offsets;           // Start offset in flat_send_indices for each neighbor [num_neighbors+1]
    int *recv_offsets;           // Start offset in flat_recv_indices for each neighbor [num_neighbors+1]

    // Communication buffers: device buffers for pack/unpack (GPU_AWARE_MPI: omp_target_alloc,
    // else: host malloc used directly for MPI too)
    double *send_buffer;         // [3 * total_send_size]
    double *recv_buffer;         // [3 * total_recv_size]

    // Host staging buffers used by GPU_AWARE_MPI path to work around UCX transports
    // (e.g. uct_mm shared-memory) that cannot handle omp_target_alloc device pointers.
    // NULL in the non-GPU_AWARE_MPI path (send_buffer/recv_buffer are already host).
    double *host_send_buffer;    // [3 * total_send_size]  (GPU_AWARE_MPI only)
    double *host_recv_buffer;    // [3 * total_recv_size]  (GPU_AWARE_MPI only)

    // MPI request arrays for non-blocking communication
    MPI_Request *requests;       // [2 * num_neighbors] for Isend/Irecv pairs
};

// Reflective boundary info - stored on GPU for efficient evaluation
struct reflective_boundary {
    int num_edges;               // Number of reflective boundary edges
    int *boundary_indices;       // Where to write in boundary_values arrays [num_edges]
    int *vol_ids;                // Interior cell IDs [num_edges]
    int *edge_ids;               // Which edge (0, 1, or 2) [num_edges]
    int mapped;                  // Whether arrays are mapped to GPU
};

// Dirichlet boundary info - constant values at boundary
struct dirichlet_boundary {
    int num_edges;               // Number of Dirichlet boundary edges
    int *boundary_indices;       // Where to write in boundary_values arrays [num_edges]
    int *vol_ids;                // Interior cell IDs [num_edges]
    int *edge_ids;               // Which edge (0, 1, or 2) [num_edges]
    double *stage_values;        // Per-edge constant stage value [num_edges]
    double *xmom_values;         // Per-edge constant xmom value  [num_edges]
    double *ymom_values;         // Per-edge constant ymom value  [num_edges]
    int mapped;                  // Whether arrays are mapped to GPU
};

// Transmissive boundary info - copies interior values to boundary
struct transmissive_boundary {
    int num_edges;               // Number of transmissive boundary edges
    int *boundary_indices;       // Where to write in boundary_values arrays [num_edges]
    int *vol_ids;                // Interior cell IDs [num_edges]
    int *edge_ids;               // Which edge (0, 1, or 2) [num_edges]
    int use_centroid;            // 1 = use centroid values, 0 = use edge values
    int mapped;                  // Whether arrays are mapped to GPU
};

// Transmissive_n_momentum_zero_t_momentum_set_stage boundary
// Sets stage from external value, keeps normal momentum, zeros tangential
struct transmissive_n_zero_t_boundary {
    int num_edges;               // Number of boundary edges
    int *boundary_indices;       // Where to write in boundary_values arrays [num_edges]
    int *vol_ids;                // Interior cell IDs [num_edges]
    int *edge_ids;               // Which edge (0, 1, or 2) [num_edges]
    double stage_value;          // Current stage value (updated each timestep from Python)
    int mapped;                  // Whether arrays are mapped to GPU
};

// Time_boundary - time-dependent Dirichlet values (function called from Python each timestep)
struct time_boundary {
    int num_edges;               // Number of time boundary edges
    int *boundary_indices;       // Where to write in boundary_values arrays [num_edges]
    int *vol_ids;                // Interior cell IDs [num_edges]
    int *edge_ids;               // Which edge (0, 1, or 2) [num_edges]
    double *stage_values;        // Per-edge stage (updated each timestep from Python) [num_edges]
    double *xmom_values;         // Per-edge xmom  [num_edges]
    double *ymom_values;         // Per-edge ymom  [num_edges]
    int mapped;                  // Whether arrays are mapped to GPU
};

// File_boundary / Field_boundary - spatially varying time-dependent values
// The Python side evaluates F(t, point_id=i) for each edge each timestep and
// pushes the resulting per-edge arrays to the device via set_values.
struct file_boundary {
    int num_edges;               // Number of file boundary edges
    int *boundary_indices;       // Where to write in boundary_values arrays [num_edges]
    int *vol_ids;                // Interior cell IDs [num_edges]
    int *edge_ids;               // Which edge (0, 1, or 2) [num_edges]
    double *stage_values;        // Per-edge stage  (updated each timestep from Python) [num_edges]
    double *xmom_values;         // Per-edge xmom   [num_edges]
    double *ymom_values;         // Per-edge ymom   [num_edges]
    int mapped;                  // Whether arrays are mapped to GPU
};

// Absorbing_wave_boundary - active-absorption open boundary.
// Ghost stage = 2*wave(t) - stage_interior; velocity-preserving momentum.
// wave_value (scalar) is updated each timestep from Python.
struct absorbing_wave_boundary {
    int num_edges;
    int *boundary_indices;
    int *vol_ids;
    int *edge_ids;
    double wave_value;           // Current wave stage (updated each timestep from Python)
    int mapped;
};

// Characteristic_wave_boundary - nonlinear Riemann-invariant open boundary.
// background_stage is fixed at init; wave_value (perturbation) updated each timestep.
struct characteristic_wave_boundary {
    int num_edges;
    int *boundary_indices;
    int *vol_ids;
    int *edge_ids;
    double wave_value;           // Stage perturbation (updated each timestep from Python)
    double background_stage;     // Still-water stage (fixed)
    int mapped;
};

// Flather_external_stage_zero_velocity_boundary - weakly-reflecting Flather-type BC.
// stage_outside (exterior stage) updated each timestep from Python; interior state
// read on device to form the characteristic-variable ghost state.
struct flather_boundary {
    int num_edges;
    int *boundary_indices;
    int *vol_ids;
    int *edge_ids;
    double stage_outside;        // Exterior stage (updated each timestep from Python)
    int mapped;
};

// Boundary edge sync buffers - pre-allocated for efficient sparse sync
// Allocated once during setup, reused every timestep
struct boundary_edge_sync {
    int num_boundary_cells;      // Number of unique boundary-adjacent cells
    int *cell_ids;               // Cell IDs to sync [num_boundary_cells] - mapped to GPU

    // Staging buffers for gather/scatter (mapped to GPU once)
    // Size = num_boundary_cells * 3 (3 edges per cell)
    int buf_size;
    double *stage_buf;
    double *xmom_buf;
    double *ymom_buf;
    double *bed_buf;
    double *height_buf;

    int initialized;
};

// Rate operator info - for GPU-accelerated rate application (rain, etc.)
// Supports both positive rates (inflow) and negative rates (extraction)
//
// MAX_RATE_OPERATORS / MAX_INLET_OPERATORS / MAX_CULVERTS are the *initial*
// heap allocation sizes.  The arrays grow automatically (doubling) so there
// is no hard limit on the number of operators.
#define MAX_RATE_OPERATORS  64
#define MAX_INLET_OPERATORS 32

struct rate_operator_info {
    int num_indices;             // Number of triangles this operator applies to
    int *indices;                // Triangle indices [num_indices] - mapped to GPU
    // Triangle area used for mass tracking [num_indices], with the ghost-cell mask
    // ALREADY BAKED IN: mass_areas[k] is 0.0 for a ghost triangle, so ghosts add
    // nothing to the influx reduction. Under MPI a rainfall polygon straddling a
    // partition boundary appears on several ranks; only the rank that OWNS a
    // triangle may count it, or the reported influx is inflated. (The *stage*
    // update is still applied to ghosts, as on the CPU path — the halo exchange
    // overwrites them.) Mirrors Rate_operator.full_indices on the Python side.
    double *mass_areas;
    int active;                  // Whether this operator slot is in use
    int mapped;                  // Whether arrays are mapped to GPU
    // Rate array caching (avoids H2D transfer every call)
    double *rate_array_cache;    // GPU-resident copy of rate array
    int rate_array_size;         // Size of cached rate array
    int rate_array_mapped;       // Whether rate array is mapped to GPU
    // Tracers carried by the rate (gpu_rate_operator_tracers_*). Appended so
    // no earlier member moves; allocated on first use.
    double *scratch_hold;        // depth before the rate is applied [num_indices]
    int hold_mapped;
};

struct rate_operators {
    struct rate_operator_info *ops; // heap-allocated, capacity entries
    int num_operators;              // Number of active operators
    int capacity;                   // Allocated size of ops array (0 = uninitialised)
    int initialized;
};

// Inlet operator info - for GPU-accelerated inlet operations
// The inlet only touches 10-100 triangles, so we do small D2H/H2D
// of just those values (~6KB) instead of full domain sync (~235MB)

struct inlet_operator_info {
    int num_indices;             // Number of inlet triangles
    int *indices;                // Triangle indices [num_indices] - mapped to GPU
    double *areas;               // Triangle areas [num_indices] - mapped to GPU
    double total_area;           // Precomputed sum of areas
    double *scratch_stages;      // Host scratch buffer [num_indices]
    double *scratch_bed;         // Host scratch buffer [num_indices]
    double *scratch_xmom;        // Host scratch buffer [num_indices]
    double *scratch_ymom;        // Host scratch buffer [num_indices]
    double *scratch_depths;      // Host scratch buffer [num_indices]
    int active;
    int mapped;
    // Tracers carried by the inlet (see gpu_inlet_apply). Appended so no
    // earlier member moves. scratch_hold is sized with the others; the tracer
    // mass buffer is sized on first use, since tracers can be registered after
    // the operator is created.
    double *scratch_hold;        // depths before the water moves [num_indices]
    double *scratch_tracer;      // tracer mass [ns * num_indices]
    int scratch_tracer_len;      // allocated (and mapped) length of scratch_tracer
};

struct inlet_operators {
    struct inlet_operator_info *ops; // heap-allocated, capacity entries
    int num_operators;
    int capacity;                    // Allocated size of ops array
};

// Culvert operator types
// MAX_CULVERTS, like MAX_RATE_OPERATORS / MAX_INLET_OPERATORS above, is the
// *initial* heap allocation size; the array grows by doubling, so there is no
// hard limit on the number of culverts. There is likewise no limit on the
// number of triangles in an inlet — see struct culvert_indices.
#define MAX_CULVERTS 64
#define CULVERT_TYPE_BOX              0
#define CULVERT_TYPE_PIPE             1
#define CULVERT_TYPE_WEIR_TRAPEZOID   2

// Static geometry parameters for one culvert
struct culvert_params {
    int type;                    // CULVERT_TYPE_BOX / PIPE / WEIR_TRAPEZOID
    double g;                    // Gravity [m/s^2] (from domain)
    double width;                // Box/trapezoid bottom width [m]
    double height;               // Box/trapezoid height [m]
    double diameter;             // Pipe diameter [m]
    double z1;                   // Trapezoid left side slope (horiz/vert); 0 for box/pipe
    double z2;                   // Trapezoid right side slope (horiz/vert); 0 for box/pipe
    double length;               // Culvert length [m]
    double manning;              // Manning's n for culvert
    double sum_loss;             // Sum of loss coefficients
    double blockage;             // Blockage fraction [0,1]
    double barrels;              // Number of barrels
    int use_velocity_head;
    int use_momentum_jet;
    int use_old_momentum_method;
    int always_use_Q_wetdry_adjustment;
    double max_velocity;
    double smoothing_timescale;
    double outward_vector_0[2];
    double outward_vector_1[2];
    double invert_elevation_0;
    double invert_elevation_1;
    int has_invert_elevation_0;
    int has_invert_elevation_1;
};

// Per-culvert indexing into domain arrays
struct culvert_indices {
    int enquiry_index_0;         // -1 if not on this rank
    int enquiry_index_1;         // -1 if not on this rank
    // Host-side staging for the inlet triangles of each end, between
    // registration (gpu_culvert_init) and the flattening done by
    // gpu_culverts_map(), which is what the device actually sees. Heap arrays
    // sized by the actual triangle count — as in rate_operator_info and
    // inlet_operator_info — so an inlet region may cover any number of
    // triangles. NULL when the corresponding *_num is 0. Owned here and freed
    // by gpu_culverts_finalize_all().
    int inlet0_num;              // 0 if no local triangles
    int *inlet0_indices;         // [inlet0_num]
    double *inlet0_areas;        // [inlet0_num]
    double inlet0_total_area;    // LOCAL area (partial if cross-boundary)
    int inlet1_num;
    int *inlet1_indices;         // [inlet1_num]
    double *inlet1_areas;        // [inlet1_num]
    double inlet1_total_area;

    // MPI topology (for cross-boundary culverts)
    int master_proc;             // Rank that computes discharge
    int enquiry_proc[2];         // Rank owning each enquiry point
    int inlet_master_proc[2];    // Master rank for each inlet region
    int is_local;                // 1 = fully local (no MPI), 0 = cross-boundary
    int mpi_tag_base;            // Base MPI tag for this culvert's messages
};

// Per-culvert dynamic state (smoothing)
struct culvert_state {
    double smooth_delta_total_energy;
    double smooth_Q;

    // Reporting stats (host-side), refreshed every apply_all on the proc that
    // computes the discharge (the master proc; zero on non-master / dry /
    // closed). Read back by the Python structure logger so mode-2 (unified /
    // GPU) .log files carry the same columns as mode-1. See
    // gpu_culverts_get_report().
    double report_gain;               // Q * timestep_star this step   [m^3]
    double report_discharge;          // instantaneous discharge       [m^3/s]
    double report_velocity;           // barrel velocity               [m/s]
    double report_driving_energy;     // inflow driving energy         [m]
    double report_delta_total_energy; // |smoothed delta total energy| [m]
};

// Culvert manager (lives inside gpu_domain)
struct culvert_operators {
    int num_culverts;
    int capacity;                    // Allocated size of params/indices/state arrays
    struct culvert_params *params;   // heap-allocated, capacity entries
    struct culvert_indices *indices; // heap-allocated, capacity entries
    struct culvert_state *state;     // heap-allocated, capacity entries
    int initialized;

    // ------------------------------------------------------------------
    // Host-side per-step working buffers, grown to hold num_culverts.
    // These were once fixed-size stack arrays of MAX_CULVERTS, which
    // silently overflowed for any model with more culverts than that
    // (MAX_CULVERTS is only the *initial* capacity — see the note above).
    // Allocated once and reused, so the hot path stays malloc-free.
    // ------------------------------------------------------------------
    int host_scratch_capacity;             // entries allocated below (0 = none)
    struct inlet_data *host_data0;         // [capacity] inlet 0 gathered data
    struct inlet_data *host_data1;         // [capacity] inlet 1 gathered data
    struct culvert_result *host_results;   // [capacity] per-culvert discharge
    struct culvert_transfer *host_transfers; // [capacity] per-culvert transfer
    struct culvert_mpi_bufs *host_mpi_bufs;  // MPI exchange buffers (opaque)

    // ------------------------------------------------------------------
    // Device-resident scratch. Everything here is mapped ONCE in
    // gpu_culverts_map() and torn down in gpu_culverts_finalize_all().
    // Constant buffers use map(to:); per-step buffers use map(alloc:) and
    // are refreshed on-device each timestep (no per-step map/alloc/free).
    // ------------------------------------------------------------------

    // Enquiry points (2 per culvert). Indices are constant → map(to:) once.
    int *scratch_enquiry_indices;   // [2*nc] centroid index of each enquiry pt
    double *scratch_stage;          // [2*nc] gathered enquiry values (D2H each step)
    double *scratch_xmom;
    double *scratch_ymom;
    double *scratch_elev;

    // Inlet triangles, flattened across all culverts. Constant metadata is
    // mapped map(to:) once; per-triangle values are read straight from the
    // domain centroid arrays on-device, so no per-triangle value buffers.
    int total_inlet_triangles;      // nt
    int *scratch_inlet_indices;     // [nt] centroid index of each triangle
    double *scratch_inlet_areas;    // [nt] triangle area (reduction weight)
    // Per-inlet contiguous range into the flattened triangle arrays. Slot
    // 2*c is inlet 0, 2*c+1 is inlet 1. Lets gather/scatter parallelise over
    // inlets (ne teams) with a sequential inner sum — no atomics, and the
    // summation order matches the old host code exactly.
    int *scratch_slot_start;        // [2*nc] first flattened triangle of the inlet
    int *scratch_slot_count;        // [2*nc] triangle count for the inlet

    // Per-inlet area-weighted sums, accumulated on-device (D2H each step).
    double *scratch_avg_stage;      // [2*nc]
    double *scratch_avg_depth;
    double *scratch_avg_xmom;
    double *scratch_avg_ymom;

    // Per-inlet scatter values, pushed H2D each step then written on-device.
    double *scratch_slot_shift;     // [2*nc] delta average depth per inlet (new - old)
    double *scratch_slot_xmom;      // [2*nc]
    double *scratch_slot_ymom;      // [2*nc]

    int mapped;

    // Tracers carried by the transfer (see gpu_culvert_scatter). Appended so
    // no earlier member moves; sized on first use for the domain's tracers.
    int tracer_scratch_ns;          // tracers the buffers are sized for (0 = none)
    double *scratch_tr_hold;        // [nt] depth before the transfer
    double *scratch_tr_T;           // [2*nc*ns] tracer given up, per inlet
    double *scratch_tr_G;           // [2*nc] water gained, per inlet
    double *scratch_tr_c;           // [nc*ns] concentration each culvert carries
};

// Max-quantities operator — device-resident running maxima of stage, depth,
// speed, and momentum magnitude collected each timestep on the GPU.
struct max_quantities_info {
    int n;                       // Number of triangles
    double velocity_zero_height; // Speed zeroed below this depth
    double *max_stage;           // [n] device-resident running max stage
    double *max_depth;           // [n] device-resident running max depth
    double *max_speed;           // [n] device-resident running max speed
    double *max_uh;              // [n] device-resident running max ||(uh,vh)||
    // Running max of each tracer's concentration, [n_tracers][n] flattened as
    // max_tracer[s*n + i], the same layout the tracer arrays themselves use.
    // NULL and n_tracers == 0 when the domain has no tracers, which is the
    // common case and costs nothing.
    double *max_tracer;
    int n_tracers;
    int initialized;
    int mapped;
};

// FLOP counter structure for performance profiling (Gordon Bell)
// Counts floating-point operations per kernel for FLOPS reporting
struct flop_counters {
    // Per-kernel FLOP counts (accumulated across all calls)
    uint64_t extrapolate_flops;      // gpu_extrapolate_second_order
    uint64_t compute_fluxes_flops;   // gpu_compute_fluxes
    uint64_t update_flops;           // gpu_update_conserved_quantities
    uint64_t protect_flops;          // gpu_protect
    uint64_t manning_flops;          // gpu_manning_friction
    uint64_t backup_flops;           // gpu_backup_conserved_quantities
    uint64_t saxpy_flops;            // gpu_saxpy_conserved_quantities
    uint64_t rate_operator_flops;    // gpu_rate_operator_apply*
    uint64_t ghost_exchange_flops;   // gpu_exchange_ghosts (pack/unpack)

    // Total FLOP count
    uint64_t total_flops;

    // Call counts for verification
    uint64_t extrapolate_calls;
    uint64_t compute_fluxes_calls;
    uint64_t update_calls;
    uint64_t protect_calls;
    uint64_t manning_calls;
    uint64_t backup_calls;
    uint64_t saxpy_calls;
    uint64_t rate_operator_calls;
    uint64_t ghost_exchange_calls;

    // Timing for FLOPS calculation (optional, can use external timer)
    double start_time;
    double elapsed_time;

    // Enable/disable flag
    int enabled;
};

// GPU domain struct - extends base domain with MPI/GPU state
struct gpu_domain {
    // Base domain struct (contains all ANUGA arrays)
    struct domain D;

    // MPI state (passed from Python via mpi4py)
    MPI_Comm comm;
    int rank;
    int nprocs;

    // GPU state
    int gpu_initialized;
    int device_id;
    int gpu_aware_mpi;           // Runtime flag: 1 if GPU-aware MPI available
    int verbose;                 // 0 = silent (default), 1 = print init/mapping messages
    // Last error message set by gpu_set_error(). Read back through
    // get_last_error() in the Cython layer so a failed init raises with the
    // reason, rather than the reason going only to a (possibly discarded)
    // per-rank stderr.
    char last_error[GPU_LAST_ERROR_LEN];

    // Active-set stepping: opt-in wet/dry fast path (domain.set_use_active_set()).
    // Prepared lazily on the first evolve step: builds the compacted
    // owned-edge list, switches the flux kernel to scatter mode, and maps the
    // classification scratch.  Serial (nprocs == 1) only -- silently disabled
    // with a warning under MPI.  See core_build_active_sets() for the
    // 2-ring-halo exactness argument.
    int use_active_set;          // requested from Python
    int as_prepared;             // lazy init has run
    int as_ready;                // >=1 full step since mapping (first step runs full)
    anuga_int *as_wet;           // scratch: wet flags         [n]
    anuga_int *as_ring1;         // scratch: 1-ring flags      [n]
    anuga_int *as_cells;         // compacted active cells     [n]
    anuga_int *as_edges;         // compacted active edges     [num_owned_edges]
    anuga_int as_counts[2];      // {n_active_cells, n_active_edges} per rebuild
    double as_cellfrac_sum;      // stats: sum of active fractions
    long   as_samples;           // stats: number of rebuilds

    // Halo exchange info
    struct halo_exchange halo;

    // Boundary conditions
    struct reflective_boundary reflective;
    struct dirichlet_boundary dirichlet;
    struct transmissive_boundary transmissive;
    struct transmissive_n_zero_t_boundary transmissive_n_zero_t;
    struct time_boundary time_bdry;
    struct file_boundary file_bdry;
    struct absorbing_wave_boundary absorbing_wave;
    struct characteristic_wave_boundary characteristic_wave;
    struct flather_boundary flather;

    // Boundary edge sync (for sparse edge value sync)
    struct boundary_edge_sync edge_sync;

    // Rate operators (rain, extraction, etc.)
    struct rate_operators rate_ops;

    // Inlet operators (GPU-accelerated inlet flow)
    struct inlet_operators inlet_ops;

    // Culvert operators (Boyd box/pipe - batched GPU gather/scatter)
    struct culvert_operators culvert_ops;

    // Max-quantities operator (device-resident running maxima)
    struct max_quantities_info max_qty;

    // Simulation parameters for GPU kernels
    double CFL;
    double evolve_max_timestep;
    double fixed_flux_timestep;  // <= 0 means disabled (use CFL-based timestep)
    // CFL-constrained step from the most recent evolve_one_*_step, BEFORE the
    // yieldstep/finaltime cap is applied. Mirrors the value legacy's
    // update_timestep() feeds into recorded_min/max_timestep, so mode-2 reports
    // the mathematical (CFL) constraint rather than the yield-limited step taken.
    double recorded_flux_timestep;
    // Manning friction variant: 1 => sloped (edge-based) like legacy when
    // domain.use_sloped_mannings is set, 0 => flat. Set from the domain in
    // get_domain_pointers(); selected in gpu_manning_friction().
    int use_sloped_mannings;

    // RK2 backup arrays (allocated on GPU)
    // These may already exist in base domain, but we track GPU copies here
    int backup_arrays_mapped;

    // FLOP counters for performance profiling (Gordon Bell)
    struct flop_counters flops;
};

// ============================================================================
// Function declarations
// ============================================================================

// Initialization and cleanup
int gpu_domain_init(struct gpu_domain *GD, MPI_Comm comm, int rank, int nprocs);
// Lazy first-step setup for active-set stepping (no-op unless
// GD->use_active_set); safe to call every step.
void gpu_active_set_prepare(struct gpu_domain *GD);
// Mean active cell fraction since mapping (1.0 if never engaged).
double gpu_active_set_mean_fraction(struct gpu_domain *GD);
void gpu_domain_finalize(struct gpu_domain *GD);

// Halo exchange setup - called once after Python domain is partitioned
int gpu_halo_init(struct gpu_domain *GD,
                  int num_neighbors,
                  int *neighbor_ranks,
                  int *send_counts,
                  int *recv_counts,
                  int *flat_send_indices,
                  int *flat_recv_indices);
void gpu_halo_finalize(struct gpu_domain *GD);

// GPU memory management
size_t gpu_estimate_required_memory(anuga_int n, anuga_int nb);
int    gpu_query_device_memory(size_t *free_bytes, size_t *total_bytes);
int    gpu_check_device_memory(struct gpu_domain *GD);
int    gpu_domain_map_arrays(struct gpu_domain *GD);
void gpu_remap_boundary_arrays(struct gpu_domain *GD);
void gpu_domain_unmap_arrays(struct gpu_domain *GD);
void gpu_domain_sync_to_device(struct gpu_domain *GD);
void gpu_domain_sync_from_device(struct gpu_domain *GD);
void gpu_domain_sync_all_from_device(struct gpu_domain *GD);  // Debug: sync ALL arrays

// Sync boundary values TO GPU (after CPU boundary evaluation)
void gpu_sync_boundary_values(struct gpu_domain *GD);

// Sync riverwall crest elevations / hydraulic properties TO GPU (after a
// host-side RiverWall.set_elevation() and friends)
void gpu_sync_riverwall_to_device(struct gpu_domain *GD);
void gpu_sync_tracer_source_to_device(struct gpu_domain *GD);
void gpu_sync_tracer_boundary_to_device(struct gpu_domain *GD);

// Sync edge values FROM GPU (before CPU boundary evaluation)
void gpu_sync_edge_values_from_device(struct gpu_domain *GD);

// Boundary edge sync - init once, sync every timestep
// Call init after boundaries are known (after set_boundary)
int gpu_boundary_edge_sync_init(struct gpu_domain *GD,
                                int num_boundary_cells,
                                int *boundary_cell_ids);
void gpu_boundary_edge_sync_finalize(struct gpu_domain *GD);
void gpu_boundary_edge_sync(struct gpu_domain *GD);  // Fast sync using pre-allocated buffers

// Reflective boundary - setup and evaluation on GPU
int gpu_reflective_init(struct gpu_domain *GD, int num_edges,
                        int *boundary_indices, int *vol_ids, int *edge_ids);
void gpu_reflective_finalize(struct gpu_domain *GD);
void gpu_evaluate_reflective_boundary(struct gpu_domain *GD);

// Dirichlet boundary - constant values at boundary
int gpu_dirichlet_init(struct gpu_domain *GD, int num_edges,
                       int *boundary_indices, int *vol_ids, int *edge_ids,
                       double *stage_values, double *xmom_values, double *ymom_values);
void gpu_dirichlet_finalize(struct gpu_domain *GD);
void gpu_evaluate_dirichlet_boundary(struct gpu_domain *GD);

// Transmissive boundary - copies interior values to boundary
int gpu_transmissive_init(struct gpu_domain *GD, int num_edges,
                          int *boundary_indices, int *vol_ids, int *edge_ids,
                          int use_centroid);
void gpu_transmissive_finalize(struct gpu_domain *GD);
void gpu_evaluate_transmissive_boundary(struct gpu_domain *GD);

// Transmissive_n_momentum_zero_t_momentum_set_stage boundary
int gpu_transmissive_n_zero_t_init(struct gpu_domain *GD, int num_edges,
                                   int *boundary_indices, int *vol_ids, int *edge_ids);
void gpu_transmissive_n_zero_t_finalize(struct gpu_domain *GD);
void gpu_transmissive_n_zero_t_set_stage(struct gpu_domain *GD, double stage_value);
void gpu_evaluate_transmissive_n_zero_t_boundary(struct gpu_domain *GD);

// File_boundary / Field_boundary - spatially varying time-dependent values
int  gpu_file_boundary_init(struct gpu_domain *GD, int num_edges,
                             int *boundary_indices, int *vol_ids, int *edge_ids);
void gpu_file_boundary_finalize(struct gpu_domain *GD);
void gpu_file_boundary_set_values(struct gpu_domain *GD,
                                   double *stage, double *xmom, double *ymom);
void gpu_evaluate_file_boundary(struct gpu_domain *GD);

// Absorbing_wave_boundary - active-absorption open boundary
int  gpu_absorbing_wave_init(struct gpu_domain *GD, int num_edges,
                              int *boundary_indices, int *vol_ids, int *edge_ids);
void gpu_absorbing_wave_finalize(struct gpu_domain *GD);
void gpu_absorbing_wave_set_value(struct gpu_domain *GD, double wave_value);
void gpu_evaluate_absorbing_wave_boundary(struct gpu_domain *GD);

// Characteristic_wave_boundary - nonlinear Riemann-invariant open boundary
int  gpu_characteristic_wave_init(struct gpu_domain *GD, int num_edges,
                                   int *boundary_indices, int *vol_ids, int *edge_ids,
                                   double background_stage);
void gpu_characteristic_wave_finalize(struct gpu_domain *GD);
void gpu_characteristic_wave_set_value(struct gpu_domain *GD, double wave_value);
void gpu_evaluate_characteristic_wave_boundary(struct gpu_domain *GD);

// Flather_external_stage_zero_velocity_boundary - weakly-reflecting characteristic BC
int  gpu_flather_init(struct gpu_domain *GD, int num_edges,
                      int *boundary_indices, int *vol_ids, int *edge_ids);
void gpu_flather_finalize(struct gpu_domain *GD);
void gpu_flather_set_value(struct gpu_domain *GD, double stage_outside);
void gpu_evaluate_flather_boundary(struct gpu_domain *GD);

// Time_boundary - time-dependent Dirichlet values
int gpu_time_boundary_init(struct gpu_domain *GD, int num_edges,
                           int *boundary_indices, int *vol_ids, int *edge_ids);
void gpu_time_boundary_finalize(struct gpu_domain *GD);
void gpu_time_boundary_set_values(struct gpu_domain *GD, double *stage, double *xmom, double *ymom);
void gpu_evaluate_time_boundary(struct gpu_domain *GD);

// Rate operators - rain, extraction, etc.
// Returns operator ID (>= 0) or -1 on error
int gpu_rate_operator_init(struct gpu_domain *GD, int num_indices, int *indices,
                           double *areas, int *full_indices, int num_full);
void gpu_rate_operator_finalize(struct gpu_domain *GD, int op_id);
void gpu_rate_operators_finalize_all(struct gpu_domain *GD);

// Apply rate operator on GPU - returns local_influx (mass added to full cells)
// rate: the rate value in m/s (scalar - computed from Python function if time-dependent)
// factor: conversion factor
// timestep: current timestep
// For negative rates, also scales momentum appropriately
// Tracers carried by a rate operator: capture depths before the rate is
// applied, then add c_in * dh where water was added and, for tracers flagged
// carry, keep the concentration where water was removed.
int gpu_rate_operator_tracers_capture(struct gpu_domain *GD, int op_id);
void gpu_rate_operator_tracers_apply(struct gpu_domain *GD, int op_id,
                                     const double *c_in, const int *carry, int ns);

double gpu_rate_operator_apply(struct gpu_domain *GD, int op_id,
                               double rate, double factor, double timestep);

// Apply rate operator with per-cell rate array (for quantity-type rates)
// rate_array: array of rate values, one per cell in indices
// rate_array_size: size of rate_array (must match num_indices or be full domain size)
// use_indices_into_rate: if 1, rate_array is full domain size, index with indices[k]
//                        if 0, rate_array matches indices size, index with k
// rate_changed: if 1, transfer new data to GPU; if 0, reuse cached data on GPU
double gpu_rate_operator_apply_array(struct gpu_domain *GD, int op_id,
                                     double *rate_array, int rate_array_size,
                                     int use_indices_into_rate,
                                     int rate_changed,
                                     double factor, double timestep);

// Max-quantities operator — GPU-resident running maxima
int  gpu_max_quantities_init(struct gpu_domain *GD, int n, double velocity_zero_height);
void gpu_max_quantities_update(struct gpu_domain *GD);
void gpu_max_quantities_get(struct gpu_domain *GD,
                            double *out_stage, double *out_depth,
                            double *out_speed, double *out_uh);
void gpu_max_tracers_get(struct gpu_domain *GD, double *out_tracer);
void gpu_max_quantities_finalize(struct gpu_domain *GD);

// Ghost exchange - the key MPI function
// Uses GPU-aware MPI if available, otherwise does D2H/H2D for small halo buffers
void gpu_exchange_ghosts(struct gpu_domain *GD);
void gpu_exchange_bed_ghosts(struct gpu_domain *GD);

// GPU kernels (stubs - will be implemented in sw_domain_gpu.c)
void gpu_extrapolate_second_order(struct gpu_domain *GD);

// Fused RK2-backup (optional) + protect + extrapolate centroid pass in one
// launch; returns the protect mass error.  Pair with gpu_extrapolate_edges.
double gpu_prepare_step(struct gpu_domain *GD, int do_backup, int zero_eu);
void gpu_extrapolate_edges(struct gpu_domain *GD, double predictor_dt);
double gpu_compute_fluxes(struct gpu_domain *GD, int substep_count, int timestep_fluxcalls);
void gpu_update_conserved_quantities(struct gpu_domain *GD, double timestep);
void gpu_backup_conserved_quantities(struct gpu_domain *GD);
void gpu_saxpy_conserved_quantities(struct gpu_domain *GD, double a, double b);
void gpu_saxpy3_conserved_quantities(struct gpu_domain *GD, double a, double b, double c);
double gpu_protect(struct gpu_domain *GD);
double gpu_compute_water_volume(struct gpu_domain *GD);
void gpu_manning_friction(struct gpu_domain *GD);

// Fused Manning + update (+ RK2 average when do_saxpy) in a single launch.
// timestep must already be known. Falls back to the separate kernels under
// sloped Manning.
void gpu_forcing_and_update(struct gpu_domain *GD, double timestep,
                            int apply_forcing, int do_saxpy,
                            double a, double b);

// Flux + apply phases: select the edge-based kernel pair automatically when
// D->edge_flux_work is allocated (and no riverwalls); cell-based otherwise.
double gpu_flux_phase(struct gpu_domain *GD, int substep_count, int timestep_fluxcalls);
int gpu_prepare_should_zero_eu(struct gpu_domain *GD);
void gpu_apply_phase(struct gpu_domain *GD, double timestep, int apply_forcing,
                     int do_saxpy, double a, double b, int substep_count);

// Full Euler step on GPU
double gpu_evolve_one_euler_step(struct gpu_domain *GD, double max_timestep, int apply_forcing);

// Full RK2 step on GPU (calls all the above in sequence)
// max_timestep: Maximum allowed timestep (respecting yieldstep/finaltime constraints)
double gpu_evolve_one_rk2_step(struct gpu_domain *GD, double max_timestep, int apply_forcing);

// Full SSP-RK3 step on GPU (Shu-Osher 3-stage)
double gpu_evolve_one_rk3_step(struct gpu_domain *GD, double max_timestep, int apply_forcing);

// ADER-2 Cauchy-Kovalewski predictor — centroid variant (advance centroids by dt)
void gpu_ader_ck_predictor(struct gpu_domain *GD, double dt);

// ADER-2 fused edge predictor — shifts edge values to Q^{n+1/2} in-place,
// leaving centroid values unchanged (no backup/restore needed).
void gpu_ader_ck_predictor_edge(struct gpu_domain *GD, double dt);

// Full ADER-2 step: single flux call with fused edge predictor.
// prev_dt is the timestep from the previous step (0 on bootstrap → Euler).
double gpu_evolve_one_ader2_step(struct gpu_domain *GD, double max_timestep, int apply_forcing, double prev_dt);

// Utility functions
int detect_gpu_aware_mpi(void);
int gpu_is_available(void);
int gpu_get_num_devices(void);
int gpu_get_initial_device(void);
int gpu_get_default_device(void);
void gpu_set_default_device(int device_id);
void gpu_set_offload_enabled(int enabled);
int gpu_get_offload_enabled(void);
int gpu_compute_device(struct gpu_domain *GD);
void print_gpu_domain_info(struct gpu_domain *GD);

// FLOP counter functions (Gordon Bell performance profiling)
void gpu_flop_counters_init(struct gpu_domain *GD);
void gpu_flop_counters_reset(struct gpu_domain *GD);
void gpu_flop_counters_enable(struct gpu_domain *GD, int enable);
void gpu_flop_counters_start_timer(struct gpu_domain *GD);
void gpu_flop_counters_stop_timer(struct gpu_domain *GD);
uint64_t gpu_flop_counters_get_total(struct gpu_domain *GD);
double gpu_flop_counters_get_flops(struct gpu_domain *GD);  // FLOP/s
void gpu_flop_counters_print(struct gpu_domain *GD);

// Per-kernel FLOP getters
uint64_t gpu_flop_counters_get_extrapolate(struct gpu_domain *GD);
uint64_t gpu_flop_counters_get_compute_fluxes(struct gpu_domain *GD);
uint64_t gpu_flop_counters_get_update(struct gpu_domain *GD);
uint64_t gpu_flop_counters_get_protect(struct gpu_domain *GD);
uint64_t gpu_flop_counters_get_manning(struct gpu_domain *GD);
uint64_t gpu_flop_counters_get_backup(struct gpu_domain *GD);
uint64_t gpu_flop_counters_get_saxpy(struct gpu_domain *GD);
uint64_t gpu_flop_counters_get_rate_operator(struct gpu_domain *GD);
uint64_t gpu_flop_counters_get_ghost_exchange(struct gpu_domain *GD);

// MPI reduction for multi-GPU (Gordon Bell)
// Returns global totals across all ranks
uint64_t gpu_flop_counters_get_global_total(struct gpu_domain *GD);
double gpu_flop_counters_get_global_flops(struct gpu_domain *GD);
void gpu_flop_counters_print_global(struct gpu_domain *GD);

// Inlet operators - GPU-accelerated inlet flow
// Returns operator ID (>= 0) or -1 on error
int gpu_inlet_operator_init(struct gpu_domain *GD, int num_indices, int *indices, double *areas);
void gpu_inlet_operator_finalize(struct gpu_domain *GD, int op_id);
void gpu_inlet_operators_finalize_all(struct gpu_domain *GD);

// Inlet operator queries (small GPU reductions returning scalars)
double gpu_inlet_get_volume(struct gpu_domain *GD, int op_id);
void gpu_inlet_get_velocities(struct gpu_domain *GD, int op_id, double *u_out, double *v_out);

// Inlet operator modifications (small GPU kernels on 10-100 triangles)
void gpu_inlet_set_depths(struct gpu_domain *GD, int op_id, double depth);
void gpu_inlet_set_xmoms(struct gpu_domain *GD, int op_id, double value);
void gpu_inlet_set_ymoms(struct gpu_domain *GD, int op_id, double value);
void gpu_inlet_set_xmoms_array(struct gpu_domain *GD, int op_id, double *values, int n);
void gpu_inlet_set_ymoms_array(struct gpu_domain *GD, int op_id, double *values, int n);
void gpu_inlet_set_stages_evenly(struct gpu_domain *GD, int op_id, double volume);

// Main entry point combining all 3 cases of inlet application
// Returns actual applied volume
double gpu_inlet_apply(struct gpu_domain *GD, int op_id, double volume,
                       double current_volume, double total_area,
                       double *vel_u, double *vel_v, int num_vel,
                       int has_velocity, double ext_vel_u, double ext_vel_v,
                       int zero_velocity,
                       const double *c_in, int n_cin, double *dmass_out);

// Culvert operators (Boyd box/pipe - batched GPU gather/scatter with MPI)
int gpu_culvert_init(struct gpu_domain *GD,
                     const struct culvert_params *params,
                     int enquiry_index_0, int enquiry_index_1,
                     int inlet0_num, int *inlet0_indices, double *inlet0_areas,
                     int inlet1_num, int *inlet1_indices, double *inlet1_areas,
                     int master_proc, int enquiry_proc_0, int enquiry_proc_1,
                     int inlet_master_proc_0, int inlet_master_proc_1,
                     int is_local, int mpi_tag_base,
                     double init_smooth_Q, double init_smooth_delta_total_energy);
void gpu_culvert_finalize(struct gpu_domain *GD, int culvert_id);
void gpu_culverts_finalize_all(struct gpu_domain *GD);
int gpu_culverts_map(struct gpu_domain *GD);   // 0, or -1 if scratch could not be allocated
void gpu_culverts_apply_all(struct gpu_domain *GD, double timestep);
// Read back a culvert's per-step reporting stats into out[5] =
// {gain, discharge, velocity, driving_energy, delta_total_energy}. Returns 0 on
// success, -1 if culvert_id is out of range. Values are non-zero only on the
// proc that computed the discharge (the master proc).
int gpu_culverts_get_report(struct gpu_domain *GD, int culvert_id, double *out);

// Record an error message on the domain (printf-style) and echo it to stderr.
// Functions that return an error code call this first so the Cython layer can
// raise with the message (get_last_error).
void gpu_set_error(struct gpu_domain *GD, const char *fmt, ...);
const char *gpu_get_last_error(const struct gpu_domain *GD);

#endif // GPU_DOMAIN_H
