// Core kernels for CPU/GPU execution
//
// These functions operate on struct domain* directly (no MPI dependency)
// and use OpenMP for parallelization. When compiled with -DCPU_ONLY_MODE,
// they run on CPU multicore. Otherwise, they can target GPU via OpenMP offload.
//
// Both sw_domain_openmp_ext and sw_domain_gpu_ext use these same kernels.

#ifndef CORE_KERNELS_H
#define CORE_KERNELS_H

#include "sw_domain.h"

// Extrapolate centroid values to edge values (second-order reconstruction).
// Runs the two passes below in sequence.
void core_extrapolate_second_order_edge(struct domain *D);

// The two passes individually: the centroid pass is cell-local (fusable into
// neighbouring cell-local kernels); the edge pass reads neighbour centroids
// and must be its own launch.
void core_extrapolate_centroid_pass(struct domain *D);
// predictor_dt != 0 fuses the ADER-2 C-K edge predictor into the tail
// (identical arithmetic to core_ader_ck_predictor_edge, one launch fewer).
void core_extrapolate_edge_pass(struct domain *D, double predictor_dt);
void core_extrapolate_edge_pass_on(struct domain *D, double predictor_dt,
                                   const anuga_int *iter, anuga_int iter_n);

// Fused RK2-backup (optional) + protect + extrapolate centroid pass in one
// launch.  Returns the protect mass error.  Follow with
// core_extrapolate_edge_pass() to complete the reconstruction.
double core_prepare_step(struct domain *D, int do_backup, int zero_eu);
// _on variants iterate an explicit cell/edge list (NULL = all) -- the
// active-set fast path for wet/dry flood domains; see core_build_active_sets.
double core_prepare_step_on(struct domain *D, int do_backup, int zero_eu,
                            const anuga_int *iter, anuga_int iter_n);

// Distribute edge values to vertices
void core_distribute_edges_to_vertices(struct domain *D);

// Update conserved quantities with explicit/semi-implicit updates
void core_apply_sediment_source(struct domain *D, double timestep);
void core_apply_bedload(struct domain *D, double timestep);
anuga_int core_apply_repose(struct domain *D);
void core_update_conserved_quantities(struct domain *D, double timestep);

// Backup conserved quantities for RK2 timestepping
void core_backup_conserved_quantities(struct domain *D);

// RK2 combination: Q = a*Q_current + b*Q_backup + c*(something)
void core_saxpy_conserved_quantities(struct domain *D, double a, double b, double c);

// Protect against negative water depths, returns mass error
double core_protect(struct domain *D);

// Fix negative cells (after update)
int core_fix_negative_cells(struct domain *D);
double core_negative_cells_volume(struct domain *D);

// Manning friction (flat, semi-implicit)
void core_manning_friction_flat_semi_implicit(struct domain *D);
void core_vegetation_friction_semi_implicit(struct domain *D);

// Manning friction (sloped, semi-implicit)
void core_manning_friction_sloped_semi_implicit(struct domain *D);

// Manning friction (sloped, semi-implicit, edge-based)
// Active per-timestep path when domain.use_sloped_mannings=True
void core_manning_friction_sloped_semi_implicit_edge_based(struct domain *D);

// Gravity term
int core_gravity(struct domain *D);

// Gravity term (well-balanced)
int core_gravity_wb(struct domain *D);

// Fused Manning friction + conserved-quantity update + optional RK2 average.
// All three are cell-local, so they run in a single kernel launch.  Only the
// flat Manning variant is inlined; sloped callers must not use this.
// timestep must already be known (RK2 substep 2 reuses substep 1's dt).
void core_forcing_and_update(struct domain *D, double timestep,
                             int apply_manning, int do_saxpy,
                             double a, double b);
void core_forcing_and_update_on(struct domain *D, double timestep,
                                int apply_manning, int do_saxpy,
                                double a, double b,
                                const anuga_int *iter, anuga_int iter_n);

// Compute fluxes using central upwind scheme
// Returns minimum timestep, stores boundary flux sum in boundary_flux_sum[substep_count]
// substep_count: which substep of RK timestepping (0 = first, only update timestep on first)
// timestep_fluxcalls: total number of flux calls per timestep (for boundary flux array indexing)
double core_compute_fluxes_central(struct domain *D, int substep_count, int timestep_fluxcalls);

// Edge-based flux computation (opt-in; active when the driver allocates
// D->edge_flux_work with 6*3n doubles).  Kernel A solves each unique edge's
// Riemann problem once (half the solves of the cell-based kernel, exactly
// antisymmetric exchange) and performs the dt / boundary-flux reductions;
// kernel B is cell-local and finishes the entire step in one launch: flux
// gather + pressure gradients + Manning + update (+ optional RK2 average).
// Not valid with riverwalls; assumes fluxes follow an extrapolate.
double core_compute_fluxes_edge_based(struct domain *D, int substep_count,
                                      int timestep_fluxcalls);
void core_flux_apply_and_update(struct domain *D, double timestep,
                                int apply_manning, int do_saxpy,
                                double a, double b, int substep_count);

// Scatter-mode fluxes: single Riemann solve per edge, both sides' scaled
// contributions accumulated straight into the (pre-zeroed) explicit updates
// with omp atomics -- no intermediate storage.  Same restrictions as the
// slot variant; max_speed_array is not maintained.
double core_compute_fluxes_scatter(struct domain *D, int substep_count,
                                   int timestep_fluxcalls);
double core_compute_fluxes_scatter_on(struct domain *D, int substep_count,
                                      int timestep_fluxcalls,
                                      const anuga_int *edges, anuga_int nedges);

// Active-set construction for the _on variants (opt-in; see the .c comment
// for the exactness argument and the driver's obligations).
void core_build_active_sets(struct domain *D,
                            anuga_int *wet_flag,
                            anuga_int *ring1_flag,
                            anuga_int *active_cells,
                            anuga_int *active_edges,
                            const anuga_int *owned_edges,
                            anuga_int num_owned_edges,
                            anuga_int *counts_out);

// ADER Cauchy-Kovalewski predictor: advance centroid values forward by dt.
// Must be called after core_extrapolate_second_order_edge().
// Recovers cell slopes from edge values, evaluates SWE time derivatives locally,
// and updates stage/xmom/ymom/height centroid values in-place.
void core_ader_ck_predictor(struct domain *D, double dt);

// Fused ADER-2 predictor: advances edge values to Q^{n+1/2}, leaving
// centroid values unchanged.  Eliminates the second extrapolation pass.
// Call after core_extrapolate_second_order_edge() + boundary update.
void core_ader_ck_predictor_edge(struct domain *D, double dt);

#endif // CORE_KERNELS_H
