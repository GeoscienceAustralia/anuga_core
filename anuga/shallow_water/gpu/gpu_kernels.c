// GPU-accelerated shallow water solver
// Split from sw_domain_gpu.c for maintainability

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <omp.h>
// MPI (or single-process stubs) come in via gpu_domain.h
#include "gpu_domain.h"

#include "gpu_device_helpers.h"

// Include pragma macros for GPU vs CPU execution
#include "gpu_omp_macros.h"

// Core kernels (shared with sw_domain_openmp_ext)
#include "core_kernels.h"

// NVTX profiling hooks (no-ops unless -DNVTX_ENABLED)
#include "gpu_nvtx.h"

// GPU compute kernels: extrapolate, flux, protect, update, etc.

void gpu_extrapolate_second_order(struct gpu_domain *GD) {
    NVTX_PUSH("gpu_extrapolate_second_order");
    // Delegate to core kernel (shared with CPU implementation)
    core_extrapolate_second_order_edge(&GD->D);

    // Count FLOPs
    if (GD->flops.enabled) {
        GD->flops.extrapolate_flops += (uint64_t)GD->D.number_of_elements * FLOPS_EXTRAPOLATE;
        GD->flops.extrapolate_calls++;
    }
    NVTX_POP();
}

// Fused RK2-backup + protect + extrapolate centroid pass: one launch for the
// three cell-local step-opening kernels (and it retires protect's separate
// height-refresh launch, which the centroid pass subsumes).
double gpu_prepare_step(struct gpu_domain *GD, int do_backup, int zero_eu) {
    NVTX_PUSH("gpu_prepare_step");
    double mass_error = core_prepare_step(&GD->D, do_backup, zero_eu);

    if (GD->flops.enabled) {
        anuga_int n = GD->D.number_of_elements;
        GD->flops.protect_flops += (uint64_t)n * FLOPS_PROTECT;
        GD->flops.protect_calls++;
        if (do_backup) {
            GD->flops.backup_flops += (uint64_t)n * FLOPS_BACKUP;
            GD->flops.backup_calls++;
        }
    }
    NVTX_POP();
    return mass_error;
}

// The extrapolation edge pass alone -- pairs with gpu_prepare_step(), which
// already ran the centroid pass.  predictor_dt != 0 fuses the ADER-2 C-K
// edge predictor into the same launch (pass 0.0 everywhere else).
void gpu_extrapolate_edges(struct gpu_domain *GD, double predictor_dt) {
    NVTX_PUSH("gpu_extrapolate_edges");
    core_extrapolate_edge_pass(&GD->D, predictor_dt);

    if (GD->flops.enabled) {
        GD->flops.extrapolate_flops += (uint64_t)GD->D.number_of_elements * FLOPS_EXTRAPOLATE;
        GD->flops.extrapolate_calls++;
    }
    NVTX_POP();
}

double gpu_compute_fluxes(struct gpu_domain *GD, int substep_count, int timestep_fluxcalls) {
    NVTX_PUSH("gpu_compute_fluxes");
    // Unified: calls core_compute_fluxes_central from core_kernels.c.
    // substep_count / timestep_fluxcalls index D.boundary_flux_sum so the Python
    // boundary_flux_integral_operator gets each RK substep's boundary flux
    // (euler/ader2 -> (0,1); rk2 -> (0,2),(1,2); rk3 -> (0,3),(1,3),(2,3)).

    double local_timestep = core_compute_fluxes_central(&GD->D, substep_count, timestep_fluxcalls);

    // Count FLOPs: 380 FLOPs per element (3 edges × flux function)
    if (GD->flops.enabled) {
        GD->flops.compute_fluxes_flops += (uint64_t)GD->D.number_of_elements * FLOPS_COMPUTE_FLUXES;
        GD->flops.compute_fluxes_calls++;
    }

    NVTX_POP();
    return local_timestep;
}

void gpu_update_conserved_quantities(struct gpu_domain *GD, double timestep) {
    NVTX_PUSH("gpu_update_conserved_quantities");
    // Delegate to core kernel
    core_update_conserved_quantities(&GD->D, timestep);

    // Count FLOPs: 21 FLOPs per element (explicit + semi-implicit update)
    if (GD->flops.enabled) {
        GD->flops.update_flops += (uint64_t)GD->D.number_of_elements * FLOPS_UPDATE;
        GD->flops.update_calls++;
    }
    NVTX_POP();
}

void gpu_backup_conserved_quantities(struct gpu_domain *GD) {
    NVTX_PUSH("gpu_backup_conserved_quantities");
    // Delegate to core kernel
    core_backup_conserved_quantities(&GD->D);

    // Count FLOPs: 0 FLOPs per element (memory copy only)
    if (GD->flops.enabled) {
        GD->flops.backup_flops += (uint64_t)GD->D.number_of_elements * FLOPS_BACKUP;
        GD->flops.backup_calls++;
    }
    NVTX_POP();
}

void gpu_saxpy_conserved_quantities(struct gpu_domain *GD, double a, double b) {
    NVTX_PUSH("gpu_saxpy_conserved_quantities");
    // Delegate to core kernel (c=0.0 means "skip division", used for RK2)
    core_saxpy_conserved_quantities(&GD->D, a, b, 0.0);

    // Also update height to match the new stage (needed for volume calculation)
    anuga_int n = GD->D.number_of_elements;
    double * restrict stage_cv = GD->D.stage_centroid_values;
    double * restrict height_cv = GD->D.height_centroid_values;
    double * restrict bed_cv = GD->D.bed_centroid_values;

    OMP_PARALLEL_LOOP
    for (anuga_int k = 0; k < n; k++) {
        height_cv[k] = fmax(stage_cv[k] - bed_cv[k], 0.0);
    }

    // Count FLOPs: 9 FLOPs per element (3 quantities × (2 mul + 1 add) + height calc)
    if (GD->flops.enabled) {
        GD->flops.saxpy_flops += (uint64_t)n * FLOPS_SAXPY;
        GD->flops.saxpy_calls++;
    }
    NVTX_POP();
}

void gpu_saxpy3_conserved_quantities(struct gpu_domain *GD, double a, double b, double c) {
    NVTX_PUSH("gpu_saxpy3_conserved_quantities");
    // Divide-by-c variant used for the final RK3 combination:
    //   Q = (a*Q_current + b*Q_backup) / c
    // Calling core with c != 0 and c != 1 triggers the division pass.
    core_saxpy_conserved_quantities(&GD->D, a, b, c);

    // Update height to match the new stage values
    anuga_int n = GD->D.number_of_elements;
    double * restrict stage_cv = GD->D.stage_centroid_values;
    double * restrict height_cv = GD->D.height_centroid_values;
    double * restrict bed_cv = GD->D.bed_centroid_values;

    OMP_PARALLEL_LOOP
    for (anuga_int k = 0; k < n; k++) {
        height_cv[k] = fmax(stage_cv[k] - bed_cv[k], 0.0);
    }

    if (GD->flops.enabled) {
        GD->flops.saxpy_flops += (uint64_t)n * FLOPS_SAXPY;
        GD->flops.saxpy_calls++;
    }
    NVTX_POP();
}

double gpu_protect(struct gpu_domain *GD) {
    NVTX_PUSH("gpu_protect");
    // Delegate to core kernel
    double mass_error = core_protect(&GD->D);

    // Also update height quantity (core_protect doesn't do this)
    anuga_int n = GD->D.number_of_elements;
    double * restrict stage_cv = GD->D.stage_centroid_values;
    double * restrict bed_cv = GD->D.bed_centroid_values;
    double * restrict height_cv = GD->D.height_centroid_values;

    OMP_PARALLEL_LOOP
    for (anuga_int k = 0; k < n; k++) {
        height_cv[k] = fmax(stage_cv[k] - bed_cv[k], 0.0);
    }

    // Count FLOPs: 5 FLOPs per element (depth check, mass error)
    if (GD->flops.enabled) {
        GD->flops.protect_flops += (uint64_t)GD->D.number_of_elements * FLOPS_PROTECT;
        GD->flops.protect_calls++;
    }

    NVTX_POP();
    return mass_error;
}

double gpu_compute_water_volume(struct gpu_domain *GD) {
    // Compute total water volume on GPU
    // Returns local volume (caller should do MPI_Allreduce for global sum)
    //
    // Volume = sum((stage - elevation) * area) for all elements

    anuga_int n = GD->D.number_of_elements;
    double volume = 0.0;

    double * restrict stage_cv = GD->D.stage_centroid_values;
    double * restrict bed_cv = GD->D.bed_centroid_values;
    anuga_geom_t * restrict areas = GD->D.areas;

    OMP_PARALLEL_LOOP_REDUCTION_PLUS(volume)
    for (anuga_int k = 0; k < n; k++) {
        double h = stage_cv[k] - bed_cv[k];
        if (h > 0.0) {
            volume += h * areas[k];
        }
    }

    return volume;
}

void gpu_manning_friction(struct gpu_domain *GD) {
    NVTX_PUSH("gpu_manning_friction");
    // Delegate to core kernel — sloped (edge-based) or flat, matching legacy's
    // friction.py dispatch on domain.use_sloped_mannings.
    if (GD->use_sloped_mannings) {
        core_manning_friction_sloped_semi_implicit_edge_based(&GD->D);
    } else {
        core_manning_friction_flat_semi_implicit(&GD->D);
    }
    // Vegetation drag (spec 8) on top of Manning; returns at once when off.
    core_vegetation_friction_semi_implicit(&GD->D);

    // Count FLOPs: 15 FLOPs per element (sqrt, pow, semi-implicit)
    if (GD->flops.enabled) {
        GD->flops.manning_flops += (uint64_t)GD->D.number_of_elements * FLOPS_MANNING;
        GD->flops.manning_calls++;
    }
    NVTX_POP();
}

// Fused Manning friction + conserved-quantity update + optional RK2 average:
// one kernel launch instead of two or three.  All three are cell-local, which
// is what makes the fusion legal -- see the note above core_forcing_and_update()
// for why compute_fluxes and extrapolate cannot join them.
//
// Falls back to the separate kernels when sloped Manning is selected, since
// that variant reads vertex values and is not inlined in the fused kernel,
// and when vegetation drag is on, which only the separate path applies.
void gpu_forcing_and_update(struct gpu_domain *GD, double timestep,
                            int apply_forcing, int do_saxpy,
                            double a, double b) {
    if (apply_forcing && (GD->use_sloped_mannings || GD->D.vegetation_mode)) {
        gpu_manning_friction(GD);
        gpu_update_conserved_quantities(GD, timestep);
        if (do_saxpy) gpu_saxpy_conserved_quantities(GD, a, b);
        return;
    }

    NVTX_PUSH("gpu_forcing_and_update");
    core_forcing_and_update(&GD->D, timestep, apply_forcing, do_saxpy, a, b);

    if (GD->flops.enabled) {
        anuga_int n = GD->D.number_of_elements;
        if (apply_forcing) {
            GD->flops.manning_flops += (uint64_t)n * FLOPS_MANNING;
            GD->flops.manning_calls++;
        }
        GD->flops.update_flops += (uint64_t)n * FLOPS_UPDATE;
        GD->flops.update_calls++;
        if (do_saxpy) {
            GD->flops.saxpy_flops += (uint64_t)n * FLOPS_SAXPY;
            GD->flops.saxpy_calls++;
        }
    }
    NVTX_POP();
}

// ----------------------------------------------------------------------------
// Flux + apply phases with automatic edge-based/cell-based selection.
//
// Edge-based (core_compute_fluxes_edge_based + core_flux_apply_and_update)
// runs when the driver has allocated D->edge_flux_work and there are no
// riverwalls (the weir corrections are one-sided and not supported there).
// Otherwise the classic cell-based flux kernel runs, followed by the fused
// forcing+update.  Both paths implement the same discretization; edge-based
// halves the Riemann solves and never materializes the explicit updates.
// ----------------------------------------------------------------------------

// Flux path selection (all phases of a step must agree):
//   FLUX_CELL    -- classic cell-based kernel (always valid)
//   FLUX_SLOT    -- edge-based pair via slot records (D->edge_flux_work set)
//   FLUX_SCATTER -- single-solve atomic scatter (D->neigh_work set, no slots)
// Riverwalls, sloped Manning and passive tracers force the cell-based path:
// the weir corrections are one-sided, sloped Manning needs the separate
// kernels that consume array-resident explicit updates, and only the
// cell-based kernel advects tracers (the single-solve kernels have no upwind
// tracer flux and never touch tracer_explicit_update).
enum { FLUX_CELL = 0, FLUX_SLOT, FLUX_SCATTER };

static int gpu_flux_mode(struct gpu_domain *GD) {
    if (GD->D.number_of_riverwall_edges != 0 || GD->use_sloped_mannings
            || GD->D.number_of_tracers > 0)
        return FLUX_CELL;
    if (GD->D.edge_flux_work != NULL) return FLUX_SLOT;
    if (GD->D.reconstruct_edge_bed == 2 && GD->D.owned_edges != NULL)
        return FLUX_SCATTER;
    return FLUX_CELL;
}

// Scatter mode needs the explicit updates zeroed before the flux kernel;
// the step's prepare launch does it for free.
int gpu_prepare_should_zero_eu(struct gpu_domain *GD) {
    return gpu_flux_mode(GD) == FLUX_SCATTER;
}

double gpu_flux_phase(struct gpu_domain *GD, int substep_count, int timestep_fluxcalls) {
    const int mode = gpu_flux_mode(GD);
    if (GD->verbose) {
        static int printed = 0;
        if (!printed) {
            printf("  flux path : %s\n",
                   mode == FLUX_SLOT ? "slot" : mode == FLUX_SCATTER ? "scatter" : "cell");
            fflush(stdout);
            printed = 1;
        }
    }
    if (mode == FLUX_CELL)
        return gpu_compute_fluxes(GD, substep_count, timestep_fluxcalls);

    NVTX_PUSH(mode == FLUX_SLOT ? "gpu_flux_edge_based" : "gpu_flux_scatter");
    double local_timestep = (mode == FLUX_SLOT)
        ? core_compute_fluxes_edge_based(&GD->D, substep_count, timestep_fluxcalls)
        : core_compute_fluxes_scatter(&GD->D, substep_count, timestep_fluxcalls);
    if (GD->flops.enabled) {
        // Half the Riemann solves of the cell-based kernel (owner side only)
        GD->flops.compute_fluxes_flops +=
            (uint64_t)GD->D.number_of_elements * FLOPS_COMPUTE_FLUXES / 2;
        GD->flops.compute_fluxes_calls++;
    }
    NVTX_POP();
    return local_timestep;
}

void gpu_apply_phase(struct gpu_domain *GD, double timestep, int apply_forcing,
                     int do_saxpy, double a, double b, int substep_count) {
    if (gpu_flux_mode(GD) == FLUX_SLOT) {
        NVTX_PUSH("gpu_flux_apply_and_update");
        core_flux_apply_and_update(&GD->D, timestep, apply_forcing, do_saxpy,
                                   a, b, substep_count);
        if (GD->flops.enabled) {
            anuga_int n = GD->D.number_of_elements;
            if (apply_forcing) {
                GD->flops.manning_flops += (uint64_t)n * FLOPS_MANNING;
                GD->flops.manning_calls++;
            }
            GD->flops.update_flops += (uint64_t)n * FLOPS_UPDATE;
            GD->flops.update_calls++;
            if (do_saxpy) {
                GD->flops.saxpy_flops += (uint64_t)n * FLOPS_SAXPY;
                GD->flops.saxpy_calls++;
            }
        }
        NVTX_POP();
        return;
    }
    gpu_forcing_and_update(GD, timestep, apply_forcing, do_saxpy, a, b);
}

void gpu_ader_ck_predictor(struct gpu_domain *GD, double dt) {
    NVTX_PUSH("gpu_ader_ck_predictor");
    core_ader_ck_predictor(&GD->D, dt);
    if (GD->flops.enabled) {
        GD->flops.extrapolate_flops += (uint64_t)GD->D.number_of_elements * FLOPS_ADER_PREDICTOR;
        GD->flops.extrapolate_calls++;
    }
    NVTX_POP();
}

void gpu_ader_ck_predictor_edge(struct gpu_domain *GD, double dt) {
    NVTX_PUSH("gpu_ader_ck_predictor_edge");
    core_ader_ck_predictor_edge(&GD->D, dt);
    if (GD->flops.enabled) {
        GD->flops.extrapolate_flops += (uint64_t)GD->D.number_of_elements * FLOPS_ADER_PREDICTOR;
        GD->flops.extrapolate_calls++;
    }
    NVTX_POP();
}

// ============================================================================
// Full ADER-2 Step
// ============================================================================

// ============================================================================
// Active-set stepping (opt-in wet/dry fast path)
// ============================================================================

void gpu_active_set_prepare(struct gpu_domain *GD) {
    if (!GD->use_active_set || GD->as_prepared) return;
    struct domain *D = &GD->D;

    if (GD->nprocs > 1) {
        if (GD->rank == 0) {
            fprintf(stderr, "gpu: active-set stepping is serial-only for now; "
                            "disabled under MPI (np=%d)\n", GD->nprocs);
        }
        GD->use_active_set = 0;
        GD->as_prepared = 1;
        return;
    }

    // The active path drives the scatter flux kernel DIRECTLY (it is the only
    // flux kernel with an edge-list variant), bypassing gpu_flux_mode.  So
    // everything gpu_flux_mode would route to the cell-based kernel has to
    // disable the mode here instead, or the active steps would silently use a
    // kernel that cannot represent the physics: riverwalls (one-sided weir
    // corrections), sloped Manning (separate kernels), and passive tracers --
    // sediment classes are tracers -- which only the cell-based kernel
    // advects.  Full stepping stays correct; only the speedup is lost.
    // number_of_tracers cannot change while the interface is mapped
    // (add_tracer tears it down), so a one-time check here is sufficient.
    const char *why = NULL;
    if (D->number_of_tracers > 0)
        why = "passive tracers / sediment classes are registered (only the "
              "cell-based flux kernel advects them)";
    else if (D->number_of_riverwall_edges != 0)
        why = "riverwalls are present (their weir corrections are one-sided)";
    else if (GD->use_sloped_mannings)
        why = "sloped Manning friction is selected";
    if (why != NULL) {
        if (GD->rank == 0) {
            fprintf(stderr, "gpu: active-set stepping disabled: %s; "
                            "running full steps\n", why);
        }
        GD->use_active_set = 0;
        GD->as_prepared = 1;
        return;
    }

    const anuga_int n = D->number_of_elements;

    // Compacted owned-slot list: every boundary slot + the larger-index side
    // of each interior edge -- one scatter thread per physical edge.  Also
    // switches the flux kernel to scatter mode (reconstruct_edge_bed = 2),
    // which the full-set steps use too, so active and full stepping share one
    // discretization.
    if (D->owned_edges == NULL) {
        anuga_int *owned = (anuga_int *)malloc((size_t)(3 * n) * sizeof(anuga_int));
        if (owned == NULL) {
            // Cannot build the compacted edge list: run full steps instead.
            if (GD->rank == 0) {
                fprintf(stderr, "gpu: active-set stepping disabled: could not "
                                "allocate the owned-edge list; running full steps\n");
            }
            GD->use_active_set = 0;
            GD->as_prepared = 1;
            return;
        }
        anuga_int ne = 0;
        for (anuga_int p2 = 0; p2 < 3 * n; p2++) {
            const anuga_int nbr = D->neighbours[p2];
            if (nbr < 0 || nbr > p2 / 3) owned[ne++] = p2;
        }
        D->owned_edges = owned;
        D->num_owned_edges = ne;
        #pragma omp target enter data map(to: owned[0:ne])
    }
    D->reconstruct_edge_bed = 2;

    GD->as_wet   = (anuga_int *)calloc((size_t)n, sizeof(anuga_int));
    GD->as_ring1 = (anuga_int *)calloc((size_t)n, sizeof(anuga_int));
    GD->as_cells = (anuga_int *)calloc((size_t)n, sizeof(anuga_int));
    GD->as_edges = (anuga_int *)calloc((size_t)D->num_owned_edges, sizeof(anuga_int));
    {
        anuga_int *w = GD->as_wet, *r1 = GD->as_ring1;
        anuga_int *c = GD->as_cells, *e = GD->as_edges;
        const anuga_int ne = D->num_owned_edges;
        #pragma omp target enter data map(alloc: w[0:n], r1[0:n], c[0:n], e[0:ne])
    }
    GD->as_ready = 0;      // first step runs the full set
    GD->as_prepared = 1;
    if (GD->verbose) {
        printf("gpu: active-set stepping enabled (%lld cells, %lld owned edges)\n",
               (long long)n, (long long)D->num_owned_edges);
    }
}

// Per-step (re)build.  Returns full-set iteration (NULL lists) until the
// first full step has run -- classification reads centroid stage, which is
// only well-defined device-side after a complete step in every compute path.
static void gpu_active_lists(struct gpu_domain *GD,
                             const anuga_int **cells, anuga_int *ncells,
                             const anuga_int **edges, anuga_int *nedges) {
    gpu_active_set_prepare(GD);
    if (!GD->use_active_set || !GD->as_ready) {
        *cells = NULL; *ncells = 0;
        *edges = GD->D.owned_edges; *nedges = GD->D.num_owned_edges;
        return;
    }
    core_build_active_sets(&GD->D, GD->as_wet, GD->as_ring1,
                           GD->as_cells, GD->as_edges,
                           GD->D.owned_edges, GD->D.num_owned_edges,
                           GD->as_counts);
    *cells = GD->as_cells; *ncells = GD->as_counts[0];
    *edges = GD->as_edges; *nedges = GD->as_counts[1];
    GD->as_cellfrac_sum +=
        (double)GD->as_counts[0] / (double)GD->D.number_of_elements;
    GD->as_samples++;
}

double gpu_active_set_mean_fraction(struct gpu_domain *GD) {
    if (GD->as_samples == 0) return 1.0;
    return GD->as_cellfrac_sum / (double)GD->as_samples;
}

double gpu_evolve_one_ader2_step(struct gpu_domain *GD, double max_timestep, int apply_forcing, double prev_dt) {
    NVTX_PUSH("gpu_evolve_one_ader2_step");
    // ADER-2 step: extrapolate Q^n → fused edge C-K predictor(prev_dt/2) →
    //              single flux call from Q^{n+1/2} → Allreduce → update Q^n
    //
    // Single-flux-call variant matching the CPU ADER-2 implementation.
    // prev_dt is the timestep from the previous step; pass 0.0 on the first
    // call to bootstrap with a plain Euler step.
    //
    // The fused edge predictor shifts edge values to Q^{n+1/2} in-place while
    // centroid values remain at Q^n, so no backup/restore is needed.

    double local_timestep, global_timestep, timestep;

    const anuga_int *as_c, *as_e; anuga_int as_nc, as_ne;
    gpu_active_lists(GD, &as_c, &as_nc, &as_e, &as_ne);

    // ========================================
    // Step 1: protect + extrapolate Q^n → edges + evaluate boundaries
    // ========================================

    // Fused protect + extrapolate centroid pass (no RK2 backup needed: the
    // C-K predictor shifts edge values in place, centroids stay at Q^n).  The
    // edge pass carries the C-K predictor in its tail (predictor_dt = half the
    // previous step's dt; 0.0 on the bootstrap step = plain Euler), so the
    // reconstruction and the shift to Q^{n+1/2} are ONE launch -- and the
    // boundaries are evaluated ONCE, from the shifted edges.  The old
    // sequence's first boundary evaluation (before the standalone predictor
    // kernel) was provably dead: the predictor never reads boundary values,
    // and the second evaluation overwrote every value the first produced.
    if (as_c) {
        const double pdt = prev_dt > 0.0 ? prev_dt * 0.5 : 0.0;
        core_prepare_step_on(&GD->D, 0, gpu_prepare_should_zero_eu(GD), as_c, as_nc);
        core_extrapolate_edge_pass_on(&GD->D, pdt, as_c, as_nc);
    } else {
        gpu_prepare_step(GD, 0, gpu_prepare_should_zero_eu(GD));
        gpu_extrapolate_edges(GD, prev_dt > 0.0 ? prev_dt * 0.5 : 0.0);
    }

    gpu_evaluate_reflective_boundary(GD);
    gpu_evaluate_dirichlet_boundary(GD);
    gpu_evaluate_transmissive_boundary(GD);
    gpu_evaluate_transmissive_n_zero_t_boundary(GD);
    gpu_evaluate_time_boundary(GD);
    gpu_evaluate_file_boundary(GD);
    gpu_evaluate_absorbing_wave_boundary(GD);
    gpu_evaluate_characteristic_wave_boundary(GD);
    gpu_evaluate_flather_boundary(GD);

    // ========================================
    // Step 3: single flux call from Q^{n+1/2} edges (or Q^n on bootstrap step)
    // ========================================

    local_timestep = as_c ? core_compute_fluxes_scatter_on(&GD->D, 0, 1, as_e, as_ne)
                          : gpu_flux_phase(GD, 0, 1);

    // ========================================
    // Step 4: Allreduce for global min CFL timestep + clip to max_timestep
    // ========================================

    static int fixed_ts_printed_ader2 = 0;
    if (GD->fixed_flux_timestep > 0.0) {
        if (GD->rank == 0 && !fixed_ts_printed_ader2) {
            printf("ADER2: Using a fixed timestep! (dt = %e)\n", GD->fixed_flux_timestep);
            fflush(stdout);
            fixed_ts_printed_ader2 = 1;
        }
        timestep = GD->fixed_flux_timestep;
        GD->recorded_flux_timestep = GD->fixed_flux_timestep;
        if (timestep > max_timestep) timestep = max_timestep;
    } else {
        if (GD->nprocs > 1) {
            MPI_Allreduce(&local_timestep, &global_timestep, 1, MPI_DOUBLE, MPI_MIN, GD->comm);
        } else {
            global_timestep = local_timestep;
        }
        timestep = GD->CFL * global_timestep;
        // CFL constraint before the yieldstep/finaltime cap (for recorded stats)
        GD->recorded_flux_timestep =
            (timestep < GD->evolve_max_timestep) ? timestep : GD->evolve_max_timestep;
        if (timestep > max_timestep) timestep = max_timestep;
    }

    // ========================================
    // Step 5: update Q^{n+1} = Q^n + timestep * R(Q^{n+1/2})
    // (Q^n centroids are unchanged — no restore needed)
    // ========================================

    // Manning + update fused into one launch (Manning has no dt dependence,
    // so evaluating it here, after the reduction, is bit-identical to the old
    // pre-reduction call -- it only accumulates into the semi-implicit terms
    // that the update consumes).
    if (as_c) core_forcing_and_update_on(&GD->D, timestep, apply_forcing, 0, 0.0, 0.0, as_c, as_nc);
    else      gpu_apply_phase(GD, timestep, apply_forcing, 0, 0.0, 0.0, 0);
    GD->as_ready = 1;

    NVTX_POP();  // gpu_evolve_one_ader2_step
    return timestep;
}

// ============================================================================
// Full Euler Step - single-step C orchestration
// ============================================================================

double gpu_evolve_one_euler_step(struct gpu_domain *GD, double max_timestep, int apply_forcing) {
    NVTX_PUSH("gpu_evolve_one_euler_step");

    double local_timestep, global_timestep, timestep;

    const anuga_int *as_c, *as_e; anuga_int as_nc, as_ne;
    gpu_active_lists(GD, &as_c, &as_nc, &as_e, &as_ne);

    // Fused protect + extrapolate centroid pass, then the edge pass --
    // same launch structure as the fused RK2/ADER2 steps.
    if (as_c) {
        core_prepare_step_on(&GD->D, 0, gpu_prepare_should_zero_eu(GD), as_c, as_nc);
        core_extrapolate_edge_pass_on(&GD->D, 0.0, as_c, as_nc);
    } else {
        gpu_prepare_step(GD, 0, gpu_prepare_should_zero_eu(GD));
        gpu_extrapolate_edges(GD, 0.0);
    }

    gpu_evaluate_reflective_boundary(GD);
    gpu_evaluate_dirichlet_boundary(GD);
    gpu_evaluate_transmissive_boundary(GD);
    gpu_evaluate_transmissive_n_zero_t_boundary(GD);
    gpu_evaluate_time_boundary(GD);
    gpu_evaluate_file_boundary(GD);
    gpu_evaluate_absorbing_wave_boundary(GD);
    gpu_evaluate_characteristic_wave_boundary(GD);
    gpu_evaluate_flather_boundary(GD);

    local_timestep = as_c ? core_compute_fluxes_scatter_on(&GD->D, 0, 1, as_e, as_ne)
                          : gpu_flux_phase(GD, 0, 1);

    static int fixed_ts_printed_euler = 0;
    if (GD->fixed_flux_timestep > 0.0) {
        if (GD->rank == 0 && !fixed_ts_printed_euler) {
            printf("Using a fixed timestep! (dt = %e)\n", GD->fixed_flux_timestep);
            fflush(stdout);
            fixed_ts_printed_euler = 1;
        }
        timestep = GD->fixed_flux_timestep;
        GD->recorded_flux_timestep = GD->fixed_flux_timestep;
        if (timestep > max_timestep) timestep = max_timestep;
    } else {
        if (GD->nprocs > 1) {
            MPI_Allreduce(&local_timestep, &global_timestep, 1, MPI_DOUBLE, MPI_MIN, GD->comm);
        } else {
            global_timestep = local_timestep;
        }
        timestep = GD->CFL * global_timestep;
        // CFL constraint before the yieldstep/finaltime cap (for recorded stats)
        GD->recorded_flux_timestep =
            (timestep < GD->evolve_max_timestep) ? timestep : GD->evolve_max_timestep;
        if (timestep > max_timestep) timestep = max_timestep;
    }

    // Manning + update fused (edge-based modes: also the flux gather)
    if (as_c) core_forcing_and_update_on(&GD->D, timestep, apply_forcing, 0, 0.0, 0.0, as_c, as_nc);
    else      gpu_apply_phase(GD, timestep, apply_forcing, 0, 0.0, 0.0, 0);
    GD->as_ready = 1;

    NVTX_POP();  // gpu_evolve_one_euler_step
    return timestep;
}

// ============================================================================
// Full RK2 Step - Orchestrates all GPU operations
// ============================================================================

double gpu_evolve_one_rk2_step(struct gpu_domain *GD, double max_timestep, int apply_forcing) {
    NVTX_PUSH("gpu_evolve_one_rk2_step");
    // Full RK2 step orchestrated entirely in C - eliminates Python round-trip overhead
    //
    // This function performs:
    // 1. Backup conserved quantities
    // 2. First Euler step (protect, extrapolate, boundaries, fluxes, forcing, update, ghost exchange)
    // 3. Second Euler step (same pattern)
    // 4. RK2 averaging (saxpy)
    //
    // Parameters:
    // - max_timestep: Maximum allowed timestep (respecting yieldstep/finaltime constraints)
    // - apply_forcing: Whether to apply forcing terms (Manning friction)
    //
    // Time-dependent boundary values (Time_boundary, Transmissive_n_zero_t) must be set
    // by Python BEFORE calling this function via set_time_boundary_values() and
    // set_transmissive_n_zero_t_stage().

    double local_timestep, global_timestep, timestep;

    const anuga_int *as_c, *as_e; anuga_int as_nc, as_ne;
    gpu_active_lists(GD, &as_c, &as_nc, &as_e, &as_ne);

    // ========================================
    // First Euler step
    // ========================================

    // RK2 backup + protect + extrapolate centroid pass, fused into one
    // cell-local launch; the edge pass (which reads neighbour centroids)
    // follows as its own launch.
    if (as_c) {
        core_prepare_step_on(&GD->D, 1, gpu_prepare_should_zero_eu(GD), as_c, as_nc);
        core_extrapolate_edge_pass_on(&GD->D, 0.0, as_c, as_nc);
    } else {
        gpu_prepare_step(GD, 1, gpu_prepare_should_zero_eu(GD));
        gpu_extrapolate_edges(GD, 0.0);
    }

    // Evaluate all GPU-supported boundary conditions
    gpu_evaluate_reflective_boundary(GD);
    gpu_evaluate_dirichlet_boundary(GD);
    gpu_evaluate_transmissive_boundary(GD);
    gpu_evaluate_transmissive_n_zero_t_boundary(GD);
    gpu_evaluate_time_boundary(GD);
    gpu_evaluate_file_boundary(GD);
    gpu_evaluate_absorbing_wave_boundary(GD);
    gpu_evaluate_characteristic_wave_boundary(GD);
    gpu_evaluate_flather_boundary(GD);

    // Compute fluxes - returns local minimum timestep
    local_timestep = as_c ? core_compute_fluxes_scatter_on(&GD->D, 0, 2, as_e, as_ne)
                          : gpu_flux_phase(GD, 0, 2);

    // Compute global timestep
    static int fixed_ts_printed = 0;
    if (GD->fixed_flux_timestep > 0.0) {
        // Fixed timestep - skip MPI allreduce entirely
        if (GD->rank == 0 && !fixed_ts_printed) {
            printf("Using a fixed timestep! (dt = %e)\n", GD->fixed_flux_timestep);
            fflush(stdout);
            fixed_ts_printed = 1;
        }
        timestep = GD->fixed_flux_timestep;
        GD->recorded_flux_timestep = GD->fixed_flux_timestep;
        if (timestep > max_timestep) {
            timestep = max_timestep;
        }
    } else {
        // MPI reduce to get global minimum timestep
        if (GD->nprocs > 1) {
            MPI_Allreduce(&local_timestep, &global_timestep, 1, MPI_DOUBLE, MPI_MIN, GD->comm);
        } else {
            global_timestep = local_timestep;
        }

        // Apply CFL condition and respect max_timestep from Python
        timestep = GD->CFL * global_timestep;
        // CFL constraint before the yieldstep/finaltime cap (for recorded stats)
        GD->recorded_flux_timestep =
            (timestep < GD->evolve_max_timestep) ? timestep : GD->evolve_max_timestep;
        if (timestep > max_timestep) {
            timestep = max_timestep;
        }
    }

    // Forcing + update, fused into one launch (edge-based: also the flux
    // gather itself -- see gpu_apply_phase)
    if (as_c) core_forcing_and_update_on(&GD->D, timestep, apply_forcing, 0, 0.0, 0.0, as_c, as_nc);
    else      gpu_apply_phase(GD, timestep, apply_forcing, 0, 0.0, 0.0, 0);

    // Ghost exchange (MPI) - sync ghost cells between processes
    if (GD->nprocs > 1) {
        gpu_exchange_ghosts(GD);
    }

    // ========================================
    // Second Euler step
    // ========================================

    if (as_c) {   // no backup on the second substep
        core_prepare_step_on(&GD->D, 0, gpu_prepare_should_zero_eu(GD), as_c, as_nc);
        core_extrapolate_edge_pass_on(&GD->D, 0.0, as_c, as_nc);
    } else {
        gpu_prepare_step(GD, 0, gpu_prepare_should_zero_eu(GD));
        gpu_extrapolate_edges(GD, 0.0);
    }

    // Evaluate boundary conditions (same as first step)
    gpu_evaluate_reflective_boundary(GD);
    gpu_evaluate_dirichlet_boundary(GD);
    gpu_evaluate_transmissive_boundary(GD);
    gpu_evaluate_transmissive_n_zero_t_boundary(GD);
    gpu_evaluate_time_boundary(GD);
    gpu_evaluate_file_boundary(GD);
    gpu_evaluate_absorbing_wave_boundary(GD);
    gpu_evaluate_characteristic_wave_boundary(GD);
    gpu_evaluate_flather_boundary(GD);

    // Compute fluxes (ignore timestep from second step)
    if (as_c) core_compute_fluxes_scatter_on(&GD->D, 1, 2, as_e, as_ne);
    else      gpu_flux_phase(GD, 1, 2);

    // Forcing + update + the RK2 average (Q_final = 0.5*Q_backup + 0.5*Q),
    // all three fused into one launch.  The timestep is the one the first
    // substep already computed, so nothing here waits on a reduction.
    if (as_c) core_forcing_and_update_on(&GD->D, timestep, apply_forcing, 1, 0.5, 0.5, as_c, as_nc);
    else      gpu_apply_phase(GD, timestep, apply_forcing, 1, 0.5, 0.5, 1);
    GD->as_ready = 1;

    NVTX_POP();  // gpu_evolve_one_rk2_step
    return timestep;
}

// ============================================================================
// Full SSP-RK3 Step (Shu-Osher)
// ============================================================================

double gpu_evolve_one_rk3_step(struct gpu_domain *GD, double max_timestep, int apply_forcing) {
    NVTX_PUSH("gpu_evolve_one_rk3_step");
    if (GD->use_active_set) {
        static int rk3_as_notice = 0;
        if (!rk3_as_notice) {
            fprintf(stderr, "gpu: active-set stepping is not implemented for "
                            "rk3 (DE2); running full steps\n");
            rk3_as_notice = 1;
        }
    }
    // Full SSP-RK3 step orchestrated entirely in C.
    //
    // Algorithm (Shu-Osher, 3rd-order strong-stability-preserving):
    //   Stage 1:      Q^(1)   = Q^n + h * L(Q^n)
    //   Intermediate: Q^(1)   = 0.25 * Q^(1) + 0.75 * Q^n     [saxpy a=0.25, b=0.75]
    //   Stage 2:      Q^(2)   = Q^(1)_mid + h * L(Q^(1)_mid)
    //   Final:        Q^{n+1} = (2 * Q^(2) + Q^n) / 3          [saxpy3 a=2, b=1, c=3]
    //
    // Ghost exchanges after Stage 1 and after the intermediate combination.
    // Time-dependent boundary values must be set by Python BEFORE calling this.

    double local_timestep, global_timestep, timestep;

    // Backup Q^n
    gpu_backup_conserved_quantities(GD);

    // ========================================
    // Stage 1: Q^(1) = Q^n + h*L(Q^n)
    // ========================================

    gpu_protect(GD);
    gpu_extrapolate_second_order(GD);

    gpu_evaluate_reflective_boundary(GD);
    gpu_evaluate_dirichlet_boundary(GD);
    gpu_evaluate_transmissive_boundary(GD);
    gpu_evaluate_transmissive_n_zero_t_boundary(GD);
    gpu_evaluate_time_boundary(GD);
    gpu_evaluate_file_boundary(GD);
    gpu_evaluate_absorbing_wave_boundary(GD);
    gpu_evaluate_characteristic_wave_boundary(GD);
    gpu_evaluate_flather_boundary(GD);

    local_timestep = gpu_compute_fluxes(GD, 0, 3);

    // Determine global timestep (same logic as RK2)
    static int fixed_ts_printed_rk3 = 0;
    if (GD->fixed_flux_timestep > 0.0) {
        if (GD->rank == 0 && !fixed_ts_printed_rk3) {
            printf("RK3: Using a fixed timestep! (dt = %e)\n", GD->fixed_flux_timestep);
            fflush(stdout);
            fixed_ts_printed_rk3 = 1;
        }
        timestep = GD->fixed_flux_timestep;
        GD->recorded_flux_timestep = GD->fixed_flux_timestep;
        if (timestep > max_timestep) timestep = max_timestep;
    } else {
        if (GD->nprocs > 1) {
            MPI_Allreduce(&local_timestep, &global_timestep, 1, MPI_DOUBLE, MPI_MIN, GD->comm);
        } else {
            global_timestep = local_timestep;
        }
        timestep = GD->CFL * global_timestep;
        // CFL constraint before the yieldstep/finaltime cap (for recorded stats)
        GD->recorded_flux_timestep =
            (timestep < GD->evolve_max_timestep) ? timestep : GD->evolve_max_timestep;
        if (timestep > max_timestep) timestep = max_timestep;
    }

    if (apply_forcing) gpu_manning_friction(GD);
    gpu_update_conserved_quantities(GD, timestep);

    if (GD->nprocs > 1) gpu_exchange_ghosts(GD);

    // ========================================
    // Stage 2: Q^(2) = Q^(1) + h*L(Q^(1))
    // ========================================

    gpu_protect(GD);
    gpu_extrapolate_second_order(GD);

    gpu_evaluate_reflective_boundary(GD);
    gpu_evaluate_dirichlet_boundary(GD);
    gpu_evaluate_transmissive_boundary(GD);
    gpu_evaluate_transmissive_n_zero_t_boundary(GD);
    gpu_evaluate_time_boundary(GD);
    gpu_evaluate_file_boundary(GD);
    gpu_evaluate_absorbing_wave_boundary(GD);
    gpu_evaluate_characteristic_wave_boundary(GD);
    gpu_evaluate_flather_boundary(GD);

    gpu_compute_fluxes(GD, 1, 3);
    if (apply_forcing) gpu_manning_friction(GD);
    gpu_update_conserved_quantities(GD, timestep);

    // Intermediate: Q = 0.25*Q^(2) + 0.75*Q^n, then sync ghost cells
    gpu_saxpy_conserved_quantities(GD, 0.25, 0.75);
    if (GD->nprocs > 1) gpu_exchange_ghosts(GD);

    // ========================================
    // Stage 3: Q^(3) = Q^(1)_mid + h*L(Q^(1)_mid)
    // ========================================

    gpu_protect(GD);
    gpu_extrapolate_second_order(GD);

    gpu_evaluate_reflective_boundary(GD);
    gpu_evaluate_dirichlet_boundary(GD);
    gpu_evaluate_transmissive_boundary(GD);
    gpu_evaluate_transmissive_n_zero_t_boundary(GD);
    gpu_evaluate_time_boundary(GD);
    gpu_evaluate_file_boundary(GD);
    gpu_evaluate_absorbing_wave_boundary(GD);
    gpu_evaluate_characteristic_wave_boundary(GD);
    gpu_evaluate_flather_boundary(GD);

    gpu_compute_fluxes(GD, 2, 3);
    if (apply_forcing) gpu_manning_friction(GD);
    gpu_update_conserved_quantities(GD, timestep);

    // Final: Q^{n+1} = (2*Q^(3) + Q^n) / 3
    gpu_saxpy3_conserved_quantities(GD, 2.0, 1.0, 3.0);

    NVTX_POP();  // gpu_evolve_one_rk3_step
    return timestep;
}

