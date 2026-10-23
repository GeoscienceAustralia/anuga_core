// GPU-accelerated max-quantities collection operator.
//
// Tracks the running maximum of stage, depth, speed, and momentum magnitude
// (||(uh, vh)||) over every cell.  All four max arrays are kept device-resident
// between calls so per-timestep cost is one parallel pass with no host<->device
// transfer.  A D2H sync is only needed when the Python operator wants to write
// the current maxima to the SWW file (at yield steps).

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <omp.h>
#include "gpu_domain.h"
#include "gpu_device_helpers.h"
#include "gpu_omp_macros.h"

// Initialise the max-quantities operator: allocates host arrays, fills initial
// values, and maps them to the device if GPU is already initialised.
// Returns 0 on success, -1 on allocation failure.
int gpu_max_quantities_init(struct gpu_domain *GD, int n, double velocity_zero_height)
{
    struct max_quantities_info *MQ = &GD->max_qty;

    if (MQ->initialized) return 0;

    MQ->n                   = n;
    MQ->velocity_zero_height = velocity_zero_height;
    MQ->mapped              = 0;
    // Fixed at init: a tracer added later would need a differently-sized
    // array, and could not get an sww variable anyway (the file's variables
    // are defined when it is created). The Python operator refuses that case.
    MQ->n_tracers           = (int)GD->D.number_of_tracers;
    MQ->max_tracer          = NULL;

    MQ->max_stage = (double*)malloc(n * sizeof(double));
    MQ->max_depth = (double*)malloc(n * sizeof(double));
    MQ->max_speed = (double*)malloc(n * sizeof(double));
    MQ->max_uh    = (double*)malloc(n * sizeof(double));

    if (MQ->n_tracers > 0) {
        MQ->max_tracer = (double*)malloc((size_t)MQ->n_tracers * n * sizeof(double));
    }

    if (!MQ->max_stage || !MQ->max_depth || !MQ->max_speed || !MQ->max_uh
            || (MQ->n_tracers > 0 && !MQ->max_tracer)) {
        gpu_set_error(GD, "ERROR: gpu_max_quantities_init: allocation failed");
        free(MQ->max_stage); free(MQ->max_depth);
        free(MQ->max_speed); free(MQ->max_uh); free(MQ->max_tracer);
        MQ->max_stage = MQ->max_depth = MQ->max_speed = MQ->max_uh = NULL;
        MQ->max_tracer = NULL;
        return -1;
    }

    // stage initialised to large negative; depth/speed/uh initialised to 0
    for (int i = 0; i < n; i++) {
        MQ->max_stage[i] = -1.0e38;
        MQ->max_depth[i] = 0.0;
        MQ->max_speed[i] = 0.0;
        MQ->max_uh[i]    = 0.0;
    }

    // Concentration is non-negative, so 0 is the floor rather than a sentinel.
    if (MQ->max_tracer != NULL) {
        for (int k = 0; k < MQ->n_tracers * n; k++) {
            MQ->max_tracer[k] = 0.0;
        }
    }

    if (GD->gpu_initialized) {
        int ni = n;
        double *ms  = MQ->max_stage;
        double *md  = MQ->max_depth;
        double *msp = MQ->max_speed;
        double *mu  = MQ->max_uh;
        #pragma omp target enter data map(to: ms[0:ni], md[0:ni], msp[0:ni], mu[0:ni])
        if (MQ->n_tracers > 0) {
            double *mt = MQ->max_tracer;
            int nt = MQ->n_tracers * ni;
            #pragma omp target enter data map(to: mt[0:nt])
        }
        MQ->mapped = 1;
    }

    MQ->initialized = 1;
    return 0;
}

// Single-pass kernel: reads stage/bed/xmom/ymom (already device-resident via
// gpu_domain_map_arrays) and updates the four device-resident max arrays.
void gpu_max_quantities_update(struct gpu_domain *GD)
{
    struct max_quantities_info *MQ = &GD->max_qty;
    if (!MQ->initialized) return;

    int n   = MQ->n;
    double vzh = MQ->velocity_zero_height;

    double * restrict stage_c = GD->D.stage_centroid_values;
    double * restrict bed_c   = GD->D.bed_centroid_values;
    double * restrict xmom_c  = GD->D.xmom_centroid_values;
    double * restrict ymom_c  = GD->D.ymom_centroid_values;

    double * restrict max_stage = MQ->max_stage;
    double * restrict max_depth = MQ->max_depth;
    double * restrict max_speed = MQ->max_speed;
    double * restrict max_uh    = MQ->max_uh;

    OMP_PARALLEL_LOOP
    for (int i = 0; i < n; i++) {
        double s  = stage_c[i];
        double b  = bed_c[i];
        double xm = xmom_c[i];
        double ym = ymom_c[i];

        double mn = sqrt(xm * xm + ym * ym);
        double d  = s - b;
        if (d < 0.0) d = 0.0;
        double v  = (d > vzh) ? mn / d : 0.0;

        if (s  > max_stage[i]) max_stage[i] = s;
        if (mn > max_uh[i])    max_uh[i]    = mn;
        if (d  > max_depth[i]) max_depth[i] = d;
        if (v  > max_speed[i]) max_speed[i] = v;
    }

    // Tracers, in their own pass so the loop above is untouched when there
    // are none. Base pointers are taken HERE, at function scope, not inside
    // the target region: GD is not mapped to the device, so a GD->... load
    // inside the region reads a host address and silently does nothing.
    // Same rule as the tracer arrays in core_kernels.c.
    const int ns = MQ->n_tracers;
    if (ns > 0) {
        // DERIVE c = m/h rather than reading tracer_centroid_values: that
        // array is only refreshed during extrapolation, at the START of a
        // step, so when this operator runs it is one step behind the
        // conserved m. Mirrors the CPU path, and the kernel's own rule that
        // a dry cell carries no concentration.
        double * restrict t_cons     = GD->D.tracer_conserved_values;
        double * restrict max_tracer = MQ->max_tracer;
        const double mah = GD->D.minimum_allowed_height;
        if (t_cons != NULL && max_tracer != NULL) {
            OMP_PARALLEL_LOOP
            for (int i = 0; i < n; i++) {
                double d = stage_c[i] - bed_c[i];
                if (d < 0.0) d = 0.0;
                double inv_h = (d > mah) ? (1.0 / d) : 0.0;
                for (int s = 0; s < ns; s++) {
                    double c = t_cons[s * n + i] * inv_h;
                    if (c > max_tracer[s * n + i]) max_tracer[s * n + i] = c;
                }
            }
        }
    }
}

// Sync the four max arrays from device to host, then copy into caller-supplied
// output buffers (which may be the Python operator's numpy array buffers).
void gpu_max_quantities_get(struct gpu_domain *GD,
                            double *out_stage, double *out_depth,
                            double *out_speed, double *out_uh)
{
    struct max_quantities_info *MQ = &GD->max_qty;
    if (!MQ->initialized) return;

    int n = MQ->n;

    if (MQ->mapped) {
        double *ms  = MQ->max_stage;
        double *md  = MQ->max_depth;
        double *msp = MQ->max_speed;
        double *mu  = MQ->max_uh;
        #pragma omp target update from(ms[0:n], md[0:n], msp[0:n], mu[0:n])
    }

    memcpy(out_stage, MQ->max_stage, n * sizeof(double));
    memcpy(out_depth, MQ->max_depth, n * sizeof(double));
    memcpy(out_speed, MQ->max_speed, n * sizeof(double));
    memcpy(out_uh,    MQ->max_uh,    n * sizeof(double));
}

// Sync the per-tracer maxima from device to host and copy them out.
// out_tracer must hold n_tracers * n doubles, laid out as [s*n + i].
// Separate from the call above so a domain with no tracers never pays for it.
void gpu_max_tracers_get(struct gpu_domain *GD, double *out_tracer)
{
    struct max_quantities_info *MQ = &GD->max_qty;
    if (!MQ->initialized || MQ->n_tracers <= 0 || MQ->max_tracer == NULL) return;

    int nt = MQ->n_tracers * MQ->n;

    if (MQ->mapped) {
        double *mt = MQ->max_tracer;
        #pragma omp target update from(mt[0:nt])
    }

    memcpy(out_tracer, MQ->max_tracer, (size_t)nt * sizeof(double));
}

// Unmap device arrays, free host memory.
void gpu_max_quantities_finalize(struct gpu_domain *GD)
{
    struct max_quantities_info *MQ = &GD->max_qty;
    if (!MQ->initialized) return;

    if (MQ->mapped) {
        int n   = MQ->n;
        double *ms  = MQ->max_stage;
        double *md  = MQ->max_depth;
        double *msp = MQ->max_speed;
        double *mu  = MQ->max_uh;
        #pragma omp target exit data map(delete: ms[0:n], md[0:n], msp[0:n], mu[0:n])
        if (MQ->n_tracers > 0 && MQ->max_tracer != NULL) {
            double *mt = MQ->max_tracer;
            int nt = MQ->n_tracers * n;
            #pragma omp target exit data map(delete: mt[0:nt])
        }
        MQ->mapped = 0;
    }

    free(MQ->max_stage);
    free(MQ->max_depth);
    free(MQ->max_speed);
    free(MQ->max_uh);
    free(MQ->max_tracer);

    MQ->max_stage = MQ->max_depth = MQ->max_speed = MQ->max_uh = NULL;
    MQ->max_tracer  = NULL;
    MQ->n_tracers   = 0;
    MQ->n           = 0;
    MQ->initialized = 0;
}
