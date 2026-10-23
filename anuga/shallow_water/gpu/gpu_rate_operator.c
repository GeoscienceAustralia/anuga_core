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
#include "gpu_omp_macros.h"

// Rate operators (rain, extraction, etc.)

// ============================================================================
// Rate Operators (rain, extraction, etc.)
// ============================================================================

// Grow the ops array to at least one more slot.  Returns 0 on success, -1 on OOM.
static int grow_rate_ops(struct rate_operators *RO) {
    int new_cap = RO->capacity == 0 ? MAX_RATE_OPERATORS : RO->capacity * 2;
    struct rate_operator_info *p = (struct rate_operator_info*)
        realloc(RO->ops, new_cap * sizeof(struct rate_operator_info));
    if (!p) {
        fprintf(stderr, "ERROR: Failed to grow rate_operators to %d slots\n", new_cap);
        return -1;
    }
    memset(p + RO->capacity, 0,
           (new_cap - RO->capacity) * sizeof(struct rate_operator_info));
    RO->ops = p;
    RO->capacity = new_cap;
    return 0;
}

int gpu_rate_operator_init(struct gpu_domain *GD, int num_indices, int *indices,
                           double *areas, int *full_indices, int num_full) {
    struct rate_operators *RO = &GD->rate_ops;

    // Find a free slot, growing heap array if needed
    int op_id = -1;
    for (int i = 0; i < RO->capacity; i++) {
        if (!RO->ops[i].active) { op_id = i; break; }
    }
    if (op_id < 0) {
        if (grow_rate_ops(RO) != 0) {
            gpu_set_error(GD, "rate operator table could not be grown");
            return -1;
        }
        for (int i = 0; i < RO->capacity; i++) {
            if (!RO->ops[i].active) { op_id = i; break; }
        }
    }
    if (op_id < 0) {
        gpu_set_error(GD, "ERROR: No free rate operator slots after grow");
        return -1;
    }

    struct rate_operator_info *op = &RO->ops[op_id];

    op->num_indices = num_indices;
    op->active = 1;
    op->mapped = 0;
    // Initialize rate array cache
    op->scratch_hold = NULL;
    op->hold_mapped = 0;
    op->rate_array_cache = NULL;
    op->rate_array_size = 0;
    op->rate_array_mapped = 0;

    if (num_indices == 0) {
        op->indices = NULL;
        op->mass_areas = NULL;
        RO->num_operators++;
        return op_id;
    }

    // Allocate and copy arrays
    op->indices = (int*)malloc(num_indices * sizeof(int));
    op->mass_areas = (double*)malloc(num_indices * sizeof(double));

    if (!op->indices || !op->mass_areas) {
        gpu_set_error(GD, "Failed to allocate rate_operator arrays");
        if (op->indices) free(op->indices);
        if (op->mass_areas) free(op->mass_areas);
        op->active = 0;
        return -1;
    }

    memcpy(op->indices, indices, num_indices * sizeof(int));

    // Bake the ghost mask into the mass-tracking areas: start at zero and fill in
    // only the triangles this rank OWNS. `full_indices` holds positions WITHIN
    // `indices` (see Rate_operator.set_full_indices), not domain triangle ids.
    //
    // NOTE the deliberate default: full_indices == NULL / num_full == 0 means this
    // rank owns NONE of the operator's triangles, so every mass_area stays 0.0 and
    // the rank contributes no influx. That is correct — the ranks that do own those
    // triangles count them. In serial, full_indices covers every index, so
    // mass_areas == areas and nothing changes.
    memset(op->mass_areas, 0, num_indices * sizeof(double));
    if (full_indices != NULL) {
        for (int f = 0; f < num_full; f++) {
            int k = full_indices[f];
            if (k >= 0 && k < num_indices) {
                op->mass_areas[k] = areas[k];
            }
        }
    }

    // Map to GPU immediately if GPU is already initialized
    if (GD->gpu_initialized) {
        int ni = op->num_indices;
        int *idx = op->indices;
        double *ar = op->mass_areas;
        #pragma omp target enter data map(to: idx[0:ni], ar[0:ni])
        op->mapped = 1;
    }

    RO->num_operators++;

    //printf("[Rank %d] Rate_operator %d initialized: %d indices, %d owned (GPU mapped: %d) "
    //       "indices=%p mass_areas=%p\n",
    //       GD->rank, op_id, num_indices, num_full, op->mapped,
    //       (void*)op->indices, (void*)op->mass_areas);
    //fflush(stdout);

    return op_id;
}

void gpu_rate_operator_finalize(struct gpu_domain *GD, int op_id) {
    if (op_id < 0 || op_id >= GD->rate_ops.capacity) return;

    struct rate_operator_info *op = &GD->rate_ops.ops[op_id];
    if (!op->active) return;

    if (op->mapped && op->num_indices > 0) {
        int ni = op->num_indices;
        int *idx = op->indices;
        double *ar = op->mass_areas;
        #pragma omp target exit data map(delete: idx[0:ni], ar[0:ni])
    }

    // Clean up rate array cache
    if (op->rate_array_mapped && op->rate_array_cache != NULL) {
        double *rac = op->rate_array_cache;
        int ras = op->rate_array_size;
        #pragma omp target exit data map(delete: rac[0:ras])
    }
    if (op->hold_mapped && op->scratch_hold != NULL) {
        double *hd = op->scratch_hold;
        int nh = op->num_indices;
        #pragma omp target exit data map(delete: hd[0:nh])
    }
    if (op->scratch_hold) free(op->scratch_hold);
    op->scratch_hold = NULL;
    op->hold_mapped = 0;
    if (op->rate_array_cache) free(op->rate_array_cache);

    if (op->indices) free(op->indices);
    if (op->mass_areas) free(op->mass_areas);

    op->indices = NULL;
    op->mass_areas = NULL;
    op->rate_array_cache = NULL;
    op->num_indices = 0;
    op->rate_array_size = 0;
    op->active = 0;
    op->mapped = 0;
    op->rate_array_mapped = 0;

    GD->rate_ops.num_operators--;
}

void gpu_rate_operators_finalize_all(struct gpu_domain *GD) {
    struct rate_operators *RO = &GD->rate_ops;
    for (int i = 0; i < RO->capacity; i++) {
        if (RO->ops[i].active) {
            gpu_rate_operator_finalize(GD, i);
        }
    }
    if (RO->ops) { free(RO->ops); RO->ops = NULL; }
    RO->capacity = 0;
    RO->initialized = 0;
}

// ============================================================================
// Tracers carried by the rate (mirrors Rate_operator's host path)
// ============================================================================
//
// Water ADDED (dh > 0) carries the inflow concentration: m += c_in * dh.
// Water REMOVED (dh < 0): a tracer flagged carry leaves with it, keeping the
// cell's concentration (m *= h_new / h_old); otherwise it stays behind, which
// is the old behaviour and right for, say, salt under evaporation. Both are
// exact per cell. Called only when a concentration or a carry flag is set, so
// an ordinary rainfall operator costs nothing extra.

int gpu_rate_operator_tracers_capture(struct gpu_domain *GD, int op_id) {
    if (op_id < 0 || op_id >= GD->rate_ops.capacity) return -1;
    struct rate_operator_info *op = &GD->rate_ops.ops[op_id];
    if (!op->active || op->num_indices == 0) return 0;
    int n = op->num_indices;
    if (!op->mapped) {
        int *idx = op->indices;
        double *ar = op->mass_areas;
        #pragma omp target enter data map(to: idx[0:n], ar[0:n])
        op->mapped = 1;
    }
    if (op->scratch_hold == NULL) {
        op->scratch_hold = (double*)malloc((size_t)n * sizeof(double));
        if (!op->scratch_hold) {
            gpu_set_error(GD, "rate operator %d: could not allocate tracer scratch", op_id);
            return -1;
        }
    }
    if (!op->hold_mapped) {
        double *hd = op->scratch_hold;
        #pragma omp target enter data map(alloc: hd[0:n])
        op->hold_mapped = 1;
    }
    int * restrict indices = op->indices;
    double * restrict hold = op->scratch_hold;
    double * restrict stage_c = GD->D.stage_centroid_values;
    double * restrict bed_c = GD->D.bed_centroid_values;
    OMP_PARALLEL_LOOP
    for (int k = 0; k < n; k++) {
        int i = indices[k];
        double h = stage_c[i] - bed_c[i];
        hold[k] = h > 0.0 ? h : 0.0;
    }
    return 0;
}

void gpu_rate_operator_tracers_apply(struct gpu_domain *GD, int op_id,
                                     const double *c_in, const int *carry, int ns) {
    if (op_id < 0 || op_id >= GD->rate_ops.capacity) return;
    struct rate_operator_info *op = &GD->rate_ops.ops[op_id];
    if (!op->active || op->num_indices == 0 || !op->hold_mapped) return;
    const anuga_int nt = GD->D.number_of_tracers;
    if (nt <= 0 || GD->D.tracer_conserved_values == NULL || ns != nt) return;
    int n = op->num_indices;
    const anuga_int N = GD->D.number_of_elements;
    const double hmin = GD->D.minimum_allowed_height;
    int * restrict indices = op->indices;
    double * restrict hold = op->scratch_hold;
    double * restrict stage_c = GD->D.stage_centroid_values;
    double * restrict bed_c = GD->D.bed_centroid_values;
    double * restrict t_cons = GD->D.tracer_conserved_values;
    double * restrict t_cv = GD->D.tracer_centroid_values;
    // The per-tracer settings are a handful of values: copied per call.
    double *cin = (double*)malloc((size_t)ns * sizeof(double));
    int *cry = (int*)malloc((size_t)ns * sizeof(int));
    if (!cin || !cry) { free(cin); free(cry); return; }
    for (int s = 0; s < ns; s++) {
        cin[s] = c_in ? c_in[s] : 0.0;
        cry[s] = carry ? carry[s] : 0;
    }
    #pragma omp target enter data map(to: cin[0:ns], cry[0:ns])
    OMP_PARALLEL_LOOP
    for (int k = 0; k < n; k++) {
        int i = indices[k];
        double hn = stage_c[i] - bed_c[i];
        if (hn < 0.0) hn = 0.0;
        double ho = hold[k];
        double dh = hn - ho;
        if (dh == 0.0) continue;
        double inv_h = hn > hmin ? 1.0 / hn : 0.0;
        double keep = ho > 0.0 ? hn / ho : 0.0;
        for (anuga_int s = 0; s < ns; s++) {
            double m = t_cons[s * N + i];
            if (dh > 0.0) {
                m += cin[s] * dh;
            } else if (cry[s]) {
                m *= keep;
            }
            t_cons[s * N + i] = m;
            t_cv[s * N + i] = m * inv_h;
        }
    }
    #pragma omp target exit data map(delete: cin[0:ns], cry[0:ns])
    free(cin);
    free(cry);
}

double gpu_rate_operator_apply(struct gpu_domain *GD, int op_id,
                               double rate, double factor, double timestep) {
    if (op_id < 0 || op_id >= GD->rate_ops.capacity) return 0.0;

    struct rate_operator_info *op = &GD->rate_ops.ops[op_id];
    if (!op->active || op->num_indices == 0) return 0.0;

    // Ensure mapped
    if (!op->mapped) {
        int ni = op->num_indices;
        int *idx = op->indices;
        double *ar = op->mass_areas;
        #pragma omp target enter data map(to: idx[0:ni], ar[0:ni])
        op->mapped = 1;
    }

    int num_indices = op->num_indices;
    int * restrict indices = op->indices;
    // Ghost-masked areas: 0.0 for triangles this rank does not own, so they add
    // nothing to the influx reduction below. See struct rate_operator_info.
    double * restrict mass_areas = op->mass_areas;

    // Domain arrays (restrict enables better optimization)
    double * restrict stage_c = GD->D.stage_centroid_values;
    double * restrict xmom_c = GD->D.xmom_centroid_values;
    double * restrict ymom_c = GD->D.ymom_centroid_values;
    double * restrict bed_c = GD->D.bed_centroid_values;

    double local_rate = factor * timestep * rate;
    double local_influx = 0.0;

    if (rate >= 0.0) {
        // Simple positive rate - just add to stage
        // Reduction for mass tracking
        OMP_PARALLEL_LOOP_REDUCTION_PLUS(local_influx)
        for (int k = 0; k < num_indices; k++) {
            int i = indices[k];
            stage_c[i] += local_rate;
            local_influx += local_rate * mass_areas[k];
        }
    } else {
        // Negative rate (extraction) - need to limit and scale momentum
        OMP_PARALLEL_LOOP_REDUCTION_PLUS(local_influx)
        for (int k = 0; k < num_indices; k++) {
            int i = indices[k];

            // Current height
            double height = stage_c[i] - bed_c[i];

            // Can't remove more water than exists
            double actual_rate = (local_rate > -height) ? local_rate : -height;

            // Scaling factor for momentum (when extracting water)
            double scale_factor;
            if (actual_rate < 0.0) {
                scale_factor = (actual_rate + height) / (height + 1.0e-10);
            } else {
                scale_factor = 1.0;
            }

            // Apply updates
            stage_c[i] += actual_rate;
            xmom_c[i] *= scale_factor;
            ymom_c[i] *= scale_factor;

            local_influx += actual_rate * mass_areas[k];
        }
    }

    // Count FLOPs: 8 FLOPs per affected cell
    if (GD->flops.enabled) {
        GD->flops.rate_operator_flops += (uint64_t)op->num_indices * FLOPS_RATE_OPERATOR;
        GD->flops.rate_operator_calls++;
    }

    return local_influx;
}

double gpu_rate_operator_apply_array(struct gpu_domain *GD, int op_id,
                                     double *rate_array, int rate_array_size,
                                     int use_indices_into_rate,
                                     int rate_changed,
                                     double factor, double timestep) {
    if (op_id < 0 || op_id >= GD->rate_ops.capacity) return 0.0;

    struct rate_operator_info *op = &GD->rate_ops.ops[op_id];
    if (!op->active || op->num_indices == 0) return 0.0;

    // Ensure operator arrays are mapped
    if (!op->mapped) {
        int ni = op->num_indices;
        int *idx = op->indices;
        double *ar = op->mass_areas;
        #pragma omp target enter data map(to: idx[0:ni], ar[0:ni])
        op->mapped = 1;
    }

    int num_indices = op->num_indices;
    int * restrict indices = op->indices;
    // Ghost-masked areas: 0.0 for triangles this rank does not own, so they add
    // nothing to the influx reduction below. See struct rate_operator_info.
    double * restrict mass_areas = op->mass_areas;

    // Domain arrays (restrict enables better optimization)
    double * restrict stage_c = GD->D.stage_centroid_values;
    double * restrict xmom_c = GD->D.xmom_centroid_values;
    double * restrict ymom_c = GD->D.ymom_centroid_values;
    double * restrict bed_c = GD->D.bed_centroid_values;

    double local_influx = 0.0;
    double ft = factor * timestep;

    // Use cached rate array on GPU (avoids H2D transfer every call)
    // Only reallocate if size changed
    if (op->rate_array_size != rate_array_size) {
        // Size changed - need to reallocate
        if (op->rate_array_mapped && op->rate_array_cache != NULL) {
            double *old_rac = op->rate_array_cache;
            int old_size = op->rate_array_size;
            #pragma omp target exit data map(delete: old_rac[0:old_size])
        }
        if (op->rate_array_cache) free(op->rate_array_cache);

        op->rate_array_cache = (double*)malloc(rate_array_size * sizeof(double));
        if (op->rate_array_cache == NULL) {
            gpu_set_error(GD, "rate operator %d: could not allocate a %d-element rate cache; rate not applied this step", op_id, rate_array_size);
            op->rate_array_size = 0;
            op->rate_array_mapped = 0;
            return 0.0;
        }
        op->rate_array_size = rate_array_size;
        op->rate_array_mapped = 0;
        rate_changed = 1;  // Force update since we reallocated
    }

    // Only transfer data to GPU if rate has changed
    if (rate_changed || !op->rate_array_mapped) {
        // Copy data to cache
        memcpy(op->rate_array_cache, rate_array, rate_array_size * sizeof(double));

        // Map or update cache on GPU
        double *rac = op->rate_array_cache;
        int ras = rate_array_size;
        if (!op->rate_array_mapped) {
            #pragma omp target enter data map(to: rac[0:ras])
            op->rate_array_mapped = 1;
        } else {
            #pragma omp target update to(rac[0:ras])
        }
    }

    // Use the GPU-resident cache
    double *gpu_rate_array = op->rate_array_cache;

    if (use_indices_into_rate) {
        // gpu_rate_array is full domain size, index with indices[k]
        OMP_PARALLEL_LOOP_REDUCTION_PLUS(local_influx)
        for (int k = 0; k < num_indices; k++) {
            int i = indices[k];
            double rate = gpu_rate_array[i];
            double local_rate = ft * rate;

            if (rate >= 0.0) {
                stage_c[i] += local_rate;
                local_influx += local_rate * mass_areas[k];
            } else {
                // Negative rate - limit and scale momentum
                double height = stage_c[i] - bed_c[i];
                double actual_rate = (local_rate > -height) ? local_rate : -height;
                double scale_factor = (actual_rate < 0.0) ?
                    (actual_rate + height) / (height + 1.0e-10) : 1.0;

                stage_c[i] += actual_rate;
                xmom_c[i] *= scale_factor;
                ymom_c[i] *= scale_factor;
                local_influx += actual_rate * mass_areas[k];
            }
        }
    } else {
        // gpu_rate_array matches indices size, index with k
        OMP_PARALLEL_LOOP_REDUCTION_PLUS(local_influx)
        for (int k = 0; k < num_indices; k++) {
            int i = indices[k];
            double rate = gpu_rate_array[k];
            double local_rate = ft * rate;

            if (rate >= 0.0) {
                stage_c[i] += local_rate;
                local_influx += local_rate * mass_areas[k];
            } else {
                // Negative rate - limit and scale momentum
                double height = stage_c[i] - bed_c[i];
                double actual_rate = (local_rate > -height) ? local_rate : -height;
                double scale_factor = (actual_rate < 0.0) ?
                    (actual_rate + height) / (height + 1.0e-10) : 1.0;

                stage_c[i] += actual_rate;
                xmom_c[i] *= scale_factor;
                ymom_c[i] *= scale_factor;
                local_influx += actual_rate * mass_areas[k];
            }
        }
    }

    // Rate array cache stays mapped on GPU for next call

    // Count FLOPs: 8 FLOPs per affected cell
    if (GD->flops.enabled) {
        GD->flops.rate_operator_flops += (uint64_t)op->num_indices * FLOPS_RATE_OPERATOR;
        GD->flops.rate_operator_calls++;
    }

    return local_influx;
}

