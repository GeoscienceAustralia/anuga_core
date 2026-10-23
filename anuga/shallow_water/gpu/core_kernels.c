// Core kernels for CPU/GPU execution
//
// These functions use OpenMP parallel loops that compile to:
// - CPU multicore: #pragma omp parallel for simd (when -DCPU_ONLY_MODE)
// - GPU offload: #pragma omp target teams loop (otherwise)
//
// Both sw_domain_openmp_ext and sw_domain_gpu_ext use these same kernels.

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <omp.h>

#include "sw_domain.h"
#include "anuga_constants.h"
#include "core_kernels.h"
#include "gpu_omp_macros.h"
#include "gpu_device_helpers.h"

// ============================================================================
// Extrapolation: centroid values -> edge values (second-order reconstruction)
// ============================================================================

void core_extrapolate_centroid_pass(struct domain *D) {
    anuga_int n = D->number_of_elements;
    double minimum_allowed_height = D->minimum_allowed_height;
    anuga_int extrapolate_velocity_second_order = D->extrapolate_velocity_second_order;

    // Parameters for hfactor computation (wet-dry limiting)
    const anuga_int n_tracers_x = D->number_of_tracers;

    double a_tmp = ANUGA_HFACTOR_A;
    double b_tmp = ANUGA_HFACTOR_B;
    double c_tmp = 1.0 / (a_tmp - b_tmp);
    double d_tmp = 1.0 - (c_tmp * a_tmp);

    // Beta values for gradient limiting
    double beta_w = D->beta_w;
    double beta_w_dry = D->beta_w_dry;
    double beta_uh = D->beta_uh;
    double beta_uh_dry = D->beta_uh_dry;
    double beta_vh = D->beta_vh;
    double beta_vh_dry = D->beta_vh_dry;

    // Extract array pointers
    double * restrict stage_cv = D->stage_centroid_values;
    double * restrict xmom_cv = D->xmom_centroid_values;
    double * restrict ymom_cv = D->ymom_centroid_values;
    double * restrict bed_cv = D->bed_centroid_values;
    double * restrict height_cv = D->height_centroid_values;

    double * restrict stage_ev = D->stage_edge_values;
    double * restrict xmom_ev = D->xmom_edge_values;
    double * restrict ymom_ev = D->ymom_edge_values;
    double * restrict bed_ev = D->bed_edge_values;
    double * restrict height_ev = D->height_edge_values;

    anuga_geom_t * restrict centroid_coords = D->centroid_coordinates;
    anuga_geom_t * restrict edge_coords = D->edge_coordinates;

    anuga_int * restrict surrogate_neighbours = D->surrogate_neighbours;
    anuga_int * restrict number_of_boundaries = D->number_of_boundaries;

    // Generic passive tracers (n_tracers_x above). See the note above
    // core_compute_fluxes_central: on a GPU build these must be loaded at
    // function scope, because D is not mapped inside the 'omp target' regions
    // below; on a CPU build they stay inside the guard.
#ifndef CPU_ONLY_MODE
    double * restrict t_cons = D->tracer_conserved_values;
    double * restrict t_cv   = D->tracer_centroid_values;
#endif
    double * restrict x_centroid_work = D->x_centroid_work;
    double * restrict y_centroid_work = D->y_centroid_work;

    // Step 1: Update centroid values
    //
    // x/y_centroid_work carry the *velocity* the limiter reconstructs from;
    // xmom/ymom_cv keep the momentum, so no restore pass is needed afterwards.
    // (Historically this was the other way round -- the work arrays saved the
    // momentum while _cv held the velocity -- which cost a third kernel launch
    // over every cell just to swap them back. See the note above Step 3.)
    OMP_PARALLEL_LOOP
    for (anuga_int k = 0; k < n; k++) {
        double stage = stage_cv[k];
        double bed = bed_cv[k];
        double xmom = xmom_cv[k];
        double ymom = ymom_cv[k];

        double dk = fmax(stage - bed, 0.0);
        height_cv[k] = dk;

        int is_dry = (dk <= minimum_allowed_height);
        int extrapolate = (extrapolate_velocity_second_order == 1) && (dk > minimum_allowed_height);

        double xmom_out = is_dry ? 0.0 : xmom;
        double ymom_out = is_dry ? 0.0 : ymom;

        double inv_dk = extrapolate ? (1.0 / dk) : 1.0;

        x_centroid_work[k] = xmom_out * inv_dk;
        y_centroid_work[k] = ymom_out * inv_dk;

        xmom_cv[k] = xmom_out;
        ymom_cv[k] = ymom_out;

        // Derive tracer concentration c = m/h from the conserved m, exactly as
        // height is derived from stage above. Dry cells carry no concentration.
        if (n_tracers_x > 0) {
#ifdef CPU_ONLY_MODE
            double * restrict t_cons = D->tracer_conserved_values;
            double * restrict t_cv   = D->tracer_centroid_values;
#endif
            double inv_h = is_dry ? 0.0 : (1.0 / dk);
            for (anuga_int s = 0; s < n_tracers_x; s++) {
                t_cv[s * n + k] = t_cons[s * n + k] * inv_h;
            }
        }
    }

}

// Manning drag factor for the ADER-2 predictor, shared by the fused edge pass
// and the standalone predictor kernel so the two cannot drift.
//
// The predictor's time derivatives carry the flux Jacobian and the gravity
// term g*h*grad(w), but no drag.  On a slope in thin flow the balance IS
// gravity against friction, so without this the half state is pulled downhill
// unresisted, the fluxes are taken from it, and the corrector's friction --
// which sees Q^n, not the half state -- cannot recover the difference.  On
// validation_tests/analytical_exact/rundown_mild_slope_coarse (3.9 mm of
// water on a 1:10 slope) the stock predictor settled at 0.448 of the
// prescribed discharge; with this, at beta 1, it holds 0.9997 -- the same as
// DE0 and DE1, which have no predictor and were always exact there.  At beta
// 0.5 it holds 0.880, which is the limiter's error and not the predictor's:
// DE0 at beta 0.5 gives 0.877 on the same case.
//
// core_ader_ck_predictor(), the centroid variant, is deliberately left alone:
// it is exposed through the extension modules but no evolve path uses it.
//
// Applied as the semi-implicit factor 1/(1 + gamma*dt) rather than folded
// into duh_dt as -gamma*uh: Manning drag is stiff in thin flow, where
// gamma*dt is order one or more, and the explicit form then overshoots
// through zero and drives a cell-to-cell momentum oscillation.  The same
// linearisation the corrector uses is unconditionally sign-preserving.
//
// gamma MUST come from the PRE-shift state, and the callers pass it that way.
// On steady uniform flow -- constant depth h0, linear bed so wx = -S0,
// constant uh = q, vh = 0 -- every derivative term above vanishes except
//     duh_dt = -g*h*wx = +g*h0*S0,
// and the friction balance is g*h0*S0 = gamma*q, so the shift gives
//     uh* = q + dt*g*h0*S0 = q*(1 + gamma*dt)
// and this factor returns exactly q.  The predictor is then stationary on the
// steady state, as the DE0/DE1/DE2 steps are, and the corrector inherits it:
// with Q^{n+1/2} == Q^n there, its own coefficient is the same either way.
//
// Evaluating gamma on the POST-shift momentum instead gives
// gamma' = gamma*(1 + gamma*dt) and leaves
//     uh_pred = q*(1 + x)/(1 + x + x^2),   x = gamma*dt,
// which is 2/3 of q at x = 1: a systematic, stiffness-dependent deficit.
//
// zs is the sloped-Manning factor sqrt(1 + zx^2 + zy^2), matching
// core_manning_friction_sloped_semi_implicit_edge_based().  Both callers get
// the bed gradient free as grad(w) - grad(h), since w = z + h.
#pragma omp declare target
static inline double core_ader_friction_factor(double g, double eta, double zs,
                                               double eps, double seven_thirds,
                                               double h_e, double uh_e,
                                               double vh_e, double dt) {
    if (eta <= 1.0e-16 || dt == 0.0 || h_e < eps) {
        return 1.0;
    }
    double speed = sqrt(uh_e * uh_e + vh_e * vh_e);
    if (speed == 0.0) {
        return 1.0;
    }
    double gamma = g * eta * eta * zs * speed / pow(h_e, seven_thirds);
    return 1.0 / (1.0 + gamma * dt);
}
#pragma omp end declare target

// The edge pass: the second-order reconstruction proper.  Reads the
// neighbours' centroid values (via surrogate_neighbours), so it MUST be a
// separate kernel launch from anything that writes centroid values -- there is
// no device-wide barrier inside an `omp target teams loop`.
//
// predictor_dt: 0.0 for a plain reconstruction (RK2 and the ADER-2 bootstrap
// step).  Non-zero fuses the ADER-2 Cauchy-Kovalewski edge predictor into the
// tail of the cell body: the just-reconstructed edge values are shifted to
// Q^{n + predictor_dt} in the same launch, reusing the dxv/dyv edge offsets
// already in registers.  The standalone core_ader_ck_predictor_edge() kernel
// computes the same arithmetic from arrays; fusing it here removes that
// kernel's full read+write sweep over the edge arrays (~200 B/cell), and --
// because the predictor never reads boundary_values -- lets the ADER-2 step
// evaluate boundaries ONCE, after the shift, instead of before and after.
void core_extrapolate_edge_pass_on(struct domain *D, double predictor_dt,
                                   const anuga_int * restrict iter, anuga_int iter_n) {
    anuga_int n = D->number_of_elements;
    double minimum_allowed_height = D->minimum_allowed_height;
    anuga_int extrapolate_velocity_second_order = D->extrapolate_velocity_second_order;
    double g_pred = D->g;                       // used by the predictor tail only

    // Parameters for hfactor computation (wet-dry limiting)
    double a_tmp = ANUGA_HFACTOR_A;
    double b_tmp = ANUGA_HFACTOR_B;
    double c_tmp = 1.0 / (a_tmp - b_tmp);
    double d_tmp = 1.0 - (c_tmp * a_tmp);

    // Beta values for gradient limiting
    double beta_w = D->beta_w;
    double beta_w_dry = D->beta_w_dry;
    double beta_uh = D->beta_uh;
    double beta_uh_dry = D->beta_uh_dry;
    double beta_vh = D->beta_vh;
    double beta_vh_dry = D->beta_vh_dry;

    // Extract array pointers
    double * restrict stage_cv = D->stage_centroid_values;
    double * restrict xmom_cv = D->xmom_centroid_values;
    double * restrict ymom_cv = D->ymom_centroid_values;
    double * restrict bed_cv = D->bed_centroid_values;
    double * restrict height_cv = D->height_centroid_values;

    double * restrict stage_ev = D->stage_edge_values;
    double * restrict xmom_ev = D->xmom_edge_values;
    double * restrict ymom_ev = D->ymom_edge_values;
    double * restrict bed_ev = D->bed_edge_values;
    double * restrict height_ev = D->height_edge_values;
    double * restrict friction_cv = D->friction_centroid_values;   // predictor drag
    const int sloped_mannings = (int)D->use_sloped_mannings;

    anuga_geom_t * restrict centroid_coords = D->centroid_coordinates;
    anuga_geom_t * restrict edge_coords = D->edge_coordinates;

    anuga_int * restrict surrogate_neighbours = D->surrogate_neighbours;
    anuga_int * restrict number_of_boundaries = D->number_of_boundaries;
    double * restrict x_centroid_work = D->x_centroid_work;
    double * restrict y_centroid_work = D->y_centroid_work;

    // Generic passive tracers. n_tracers == 0 in every ordinary run; keep only
    // the loop-invariant count live in the common path.
    const anuga_int n_tracers = D->number_of_tracers;
    const double beta_tracer = D->beta_tracer;
    // See the note above core_compute_fluxes_central: on a GPU build these must
    // be loaded at function scope, because D is not mapped inside the 'omp
    // target' regions below; on a CPU build they stay inside the guard.
#ifndef CPU_ONLY_MODE
    double * restrict t_cv   = D->tracer_centroid_values;
    double * restrict t_ev   = D->tracer_edge_values;
#endif

    // Step 2: Main extrapolation loop
    const anuga_int loop_n = iter ? iter_n : n;
    OMP_PARALLEL_LOOP
    for (anuga_int q = 0; q < loop_n; q++) {
        const anuga_int k = iter ? iter[q] : q;
        anuga_int k2 = k * 2;
        anuga_int k3 = k * 3;
        anuga_int k6 = k * 6;

        double xv0 = edge_coords[k6 + 0];
        double yv0 = edge_coords[k6 + 1];
        double xv1 = edge_coords[k6 + 2];
        double yv1 = edge_coords[k6 + 3];
        double xv2 = edge_coords[k6 + 4];
        double yv2 = edge_coords[k6 + 5];

        double x = centroid_coords[k2 + 0];
        double y = centroid_coords[k2 + 1];

        double dxv0 = xv0 - x;
        double dxv1 = xv1 - x;
        double dxv2 = xv2 - x;
        double dyv0 = yv0 - y;
        double dyv1 = yv1 - y;
        double dyv2 = yv2 - y;

        anuga_int k0 = surrogate_neighbours[k3 + 0];
        anuga_int k1 = surrogate_neighbours[k3 + 1];
        anuga_int sn2 = surrogate_neighbours[k3 + 2];

        double x0 = centroid_coords[2 * k0 + 0];
        double y0 = centroid_coords[2 * k0 + 1];
        double x1 = centroid_coords[2 * k1 + 0];
        double y1 = centroid_coords[2 * k1 + 1];
        double x2 = centroid_coords[2 * sn2 + 0];
        double y2 = centroid_coords[2 * sn2 + 1];

        double dx1 = x1 - x0;
        double dx2 = x2 - x0;
        double dy1 = y1 - y0;
        double dy2 = y2 - y0;

        double area2 = dy2 * dx1 - dy1 * dx2;

        int dry = ((height_cv[k0] < minimum_allowed_height) || (k0 == k)) &&
                  ((height_cv[k1] < minimum_allowed_height) || (k1 == k)) &&
                  ((height_cv[sn2] < minimum_allowed_height) || (sn2 == k));

        if (dry) {
            x_centroid_work[k] = 0.0;
            xmom_cv[k] = 0.0;
            y_centroid_work[k] = 0.0;
            ymom_cv[k] = 0.0;
        }

        int num_boundaries = number_of_boundaries[k];

        if (num_boundaries == 3) {
            double stage_c = stage_cv[k];
            double xmom_c = x_centroid_work[k];
            double ymom_c = y_centroid_work[k];
            double height_c = height_cv[k];
            double bed_c = bed_cv[k];

            for (int i = 0; i < 3; i++) {
                stage_ev[k3 + i] = stage_c;
                xmom_ev[k3 + i] = xmom_c;
                ymom_ev[k3 + i] = ymom_c;
                height_ev[k3 + i] = height_c;
                bed_ev[k3 + i] = bed_c;
            }

            if (n_tracers > 0) {
#ifdef CPU_ONLY_MODE
                double * restrict t_cv = D->tracer_centroid_values;
                double * restrict t_ev = D->tracer_edge_values;
#endif
                for (anuga_int sidx = 0; sidx < n_tracers; sidx++) {
                    double tc = t_cv[sidx * n + k];
                    t_ev[sidx * 3 * n + k3 + 0] = tc;
                    t_ev[sidx * 3 * n + k3 + 1] = tc;
                    t_ev[sidx * 3 * n + k3 + 2] = tc;
                }
            }

        } else if (num_boundaries <= 1) {
            double hc = height_cv[k];
            double h0 = height_cv[k0];
            double h1 = height_cv[k1];
            double h2 = height_cv[sn2];

            double hmin = fmin(fmin(h0, fmin(h1, h2)), hc);
            double hmax = fmax(fmax(h0, fmax(h1, h2)), hc);

            double tmp1 = c_tmp * fmax(hmin, 0.0) / fmax(hc, 1.0e-06) + d_tmp;
            double tmp2 = c_tmp * fmax(hc, 0.0) / fmax(hmax, 1.0e-06) + d_tmp;
            double hfactor = fmax(0.0, fmin(tmp1, fmin(tmp2, 1.0)));

            hfactor = fmin(1.2 * fmax(hmin - minimum_allowed_height, 0.0) /
                           (fmax(hmin, 0.0) + minimum_allowed_height), hfactor);

            double inv_area2 = 1.0 / area2;
            double edge_vals[3];

            // Stage
            double beta_stage = beta_w_dry + (beta_w - beta_w_dry) * hfactor;
            if (beta_stage > 0.0) {
                gpu_calc_edge_values_with_gradient(
                    stage_cv[k], stage_cv[k0], stage_cv[k1], stage_cv[sn2],
                    dxv0, dxv1, dxv2, dyv0, dyv1, dyv2,
                    dx1, dx2, dy1, dy2, inv_area2, beta_stage, edge_vals);
            } else {
                gpu_set_constant_edge_values(stage_cv[k], edge_vals);
            }
            stage_ev[k3 + 0] = edge_vals[0];
            stage_ev[k3 + 1] = edge_vals[1];
            stage_ev[k3 + 2] = edge_vals[2];

            // Height (same beta as stage)
            if (beta_stage > 0.0) {
                gpu_calc_edge_values_with_gradient(
                    height_cv[k], height_cv[k0], height_cv[k1], height_cv[sn2],
                    dxv0, dxv1, dxv2, dyv0, dyv1, dyv2,
                    dx1, dx2, dy1, dy2, inv_area2, beta_stage, edge_vals);
            } else {
                gpu_set_constant_edge_values(height_cv[k], edge_vals);
            }
            height_ev[k3 + 0] = edge_vals[0];
            height_ev[k3 + 1] = edge_vals[1];
            height_ev[k3 + 2] = edge_vals[2];

            // X-momentum
            double beta_xmom = beta_uh_dry + (beta_uh - beta_uh_dry) * hfactor;
            if (beta_xmom > 0.0) {
                gpu_calc_edge_values_with_gradient(
                    x_centroid_work[k], x_centroid_work[k0], x_centroid_work[k1], x_centroid_work[sn2],
                    dxv0, dxv1, dxv2, dyv0, dyv1, dyv2,
                    dx1, dx2, dy1, dy2, inv_area2, beta_xmom, edge_vals);
            } else {
                gpu_set_constant_edge_values(x_centroid_work[k], edge_vals);
            }
            xmom_ev[k3 + 0] = edge_vals[0];
            xmom_ev[k3 + 1] = edge_vals[1];
            xmom_ev[k3 + 2] = edge_vals[2];

            // Y-momentum
            double beta_ymom = beta_vh_dry + (beta_vh - beta_vh_dry) * hfactor;
            if (beta_ymom > 0.0) {
                gpu_calc_edge_values_with_gradient(
                    y_centroid_work[k], y_centroid_work[k0], y_centroid_work[k1], y_centroid_work[sn2],
                    dxv0, dxv1, dxv2, dyv0, dyv1, dyv2,
                    dx1, dx2, dy1, dy2, inv_area2, beta_ymom, edge_vals);
            } else {
                gpu_set_constant_edge_values(y_centroid_work[k], edge_vals);
            }
            ymom_ev[k3 + 0] = edge_vals[0];
            ymom_ev[k3 + 1] = edge_vals[1];
            ymom_ev[k3 + 2] = edge_vals[2];

            // Tracers. Reconstruct c (the intensive variable) rather than the
            // conserved h*c: the limiter then bounds each edge value by the
            // cell-and-neighbour range of c, which is what preserves positivity
            // and prevents spurious extrema where h varies sharply.
            if (n_tracers > 0) {
#ifdef CPU_ONLY_MODE
                double * restrict t_cv = D->tracer_centroid_values;
                double * restrict t_ev = D->tracer_edge_values;
#endif
                double beta_t = beta_tracer * hfactor;
                for (anuga_int sidx = 0; sidx < n_tracers; sidx++) {
                    anuga_int off = sidx * n;
                    if (beta_t > 0.0) {
                        gpu_calc_edge_values_with_gradient(
                            t_cv[off + k], t_cv[off + k0], t_cv[off + k1], t_cv[off + sn2],
                            dxv0, dxv1, dxv2, dyv0, dyv1, dyv2,
                            dx1, dx2, dy1, dy2, inv_area2, beta_t, edge_vals);
                    } else {
                        gpu_set_constant_edge_values(t_cv[off + k], edge_vals);
                    }
                    t_ev[sidx * 3 * n + k3 + 0] = edge_vals[0];
                    t_ev[sidx * 3 * n + k3 + 1] = edge_vals[1];
                    t_ev[sidx * 3 * n + k3 + 2] = edge_vals[2];
                }
            }

        } else {
            // Number of boundaries == 2
            // One internal neighbour, gradient is in direction of neighbour's centroid
            // Find the only internal neighbour
            anuga_int kn = k;
            for (int i = 0; i < 3; i++) {
                anuga_int sn = surrogate_neighbours[k3 + i];
                if (sn != k) {
                    kn = sn;
                    break;
                }
            }

            // Compute gradient projection between centroids
            double xn = centroid_coords[2 * kn + 0];
            double yn = centroid_coords[2 * kn + 1];
            double dx = xn - x;
            double dy = yn - y;
            double dist2 = dx * dx + dy * dy;

            double grad_dx2 = (dist2 > 0.0) ? dx / dist2 : 0.0;
            double grad_dy2 = (dist2 > 0.0) ? dy / dist2 : 0.0;

            double dqv[3], qmin, qmax, dq1;

            // Stage
            dq1 = stage_cv[kn] - stage_cv[k];
            gpu_compute_dqv_from_gradient(dq1, grad_dx2, grad_dy2,
                                          dxv0, dxv1, dxv2, dyv0, dyv1, dyv2, dqv);
            gpu_compute_qmin_qmax_from_dq1(dq1, &qmin, &qmax);
            gpu_limit_gradient(dqv, qmin, qmax, beta_w);
            stage_ev[k3 + 0] = stage_cv[k] + dqv[0];
            stage_ev[k3 + 1] = stage_cv[k] + dqv[1];
            stage_ev[k3 + 2] = stage_cv[k] + dqv[2];

            // Height
            dq1 = height_cv[kn] - height_cv[k];
            gpu_compute_dqv_from_gradient(dq1, grad_dx2, grad_dy2,
                                          dxv0, dxv1, dxv2, dyv0, dyv1, dyv2, dqv);
            gpu_compute_qmin_qmax_from_dq1(dq1, &qmin, &qmax);
            gpu_limit_gradient(dqv, qmin, qmax, beta_w);
            height_ev[k3 + 0] = height_cv[k] + dqv[0];
            height_ev[k3 + 1] = height_cv[k] + dqv[1];
            height_ev[k3 + 2] = height_cv[k] + dqv[2];

            // X-momentum
            dq1 = x_centroid_work[kn] - x_centroid_work[k];
            gpu_compute_dqv_from_gradient(dq1, grad_dx2, grad_dy2,
                                          dxv0, dxv1, dxv2, dyv0, dyv1, dyv2, dqv);
            gpu_compute_qmin_qmax_from_dq1(dq1, &qmin, &qmax);
            gpu_limit_gradient(dqv, qmin, qmax, beta_w);
            xmom_ev[k3 + 0] = x_centroid_work[k] + dqv[0];
            xmom_ev[k3 + 1] = x_centroid_work[k] + dqv[1];
            xmom_ev[k3 + 2] = x_centroid_work[k] + dqv[2];

            // Y-momentum
            dq1 = y_centroid_work[kn] - y_centroid_work[k];
            gpu_compute_dqv_from_gradient(dq1, grad_dx2, grad_dy2,
                                          dxv0, dxv1, dxv2, dyv0, dyv1, dyv2, dqv);
            gpu_compute_qmin_qmax_from_dq1(dq1, &qmin, &qmax);
            gpu_limit_gradient(dqv, qmin, qmax, beta_w);
            ymom_ev[k3 + 0] = y_centroid_work[k] + dqv[0];
            ymom_ev[k3 + 1] = y_centroid_work[k] + dqv[1];
            ymom_ev[k3 + 2] = y_centroid_work[k] + dqv[2];

            // Tracers, 1D gradient toward the one internal neighbour
            if (n_tracers > 0) {
#ifdef CPU_ONLY_MODE
                double * restrict t_cv = D->tracer_centroid_values;
                double * restrict t_ev = D->tracer_edge_values;
#endif
                for (anuga_int sidx = 0; sidx < n_tracers; sidx++) {
                    anuga_int off = sidx * n;
                    double tk = t_cv[off + k];
                    if (beta_tracer > 0.0) {
                        dq1 = t_cv[off + kn] - tk;
                        gpu_compute_dqv_from_gradient(dq1, grad_dx2, grad_dy2,
                                                      dxv0, dxv1, dxv2, dyv0, dyv1, dyv2, dqv);
                        gpu_compute_qmin_qmax_from_dq1(dq1, &qmin, &qmax);
                        gpu_limit_gradient(dqv, qmin, qmax, beta_tracer);
                    } else {
                        dqv[0] = 0.0; dqv[1] = 0.0; dqv[2] = 0.0;
                    }
                    t_ev[sidx * 3 * n + k3 + 0] = tk + dqv[0];
                    t_ev[sidx * 3 * n + k3 + 1] = tk + dqv[1];
                    t_ev[sidx * 3 * n + k3 + 2] = tk + dqv[2];
                }
            }
        }

        // Convert velocity edge values back to momentum if needed
        if (extrapolate_velocity_second_order == 1) {
            for (int i = 0; i < 3; i++) {
                double dk = height_ev[k3 + i];
                xmom_ev[k3 + i] *= dk;
                ymom_ev[k3 + i] *= dk;
            }
        }

        // Compute bed edge values from stage - height
        for (int i = 0; i < 3; i++) {
            bed_ev[k3 + i] = stage_ev[k3 + i] - height_ev[k3 + i];
        }

        // ---- fused ADER-2 C-K edge predictor (see the function comment).
        // Identical arithmetic to core_ader_ck_predictor_edge(), with the
        // edge-offset vectors dxv0/dyv0/dxv1/dyv1 reused from the limiter
        // geometry above and the just-written edge values re-read from this
        // thread's own stores (register/L1-resident, never a remote gather).
        if (predictor_dt != 0.0) {
            double det_p = dxv0 * dyv1 - dxv1 * dyv0;
            if (fabs(det_p) >= 1.0e-20) {
                double inv_det = 1.0 / det_p;

                double w_c  = stage_cv[k];
                double h_c  = fmax(w_c - bed_cv[k], 0.0);
                double uh_c = xmom_cv[k];
                double vh_c = ymom_cv[k];

                double inv_h_c = (h_c > minimum_allowed_height) ? 1.0 / h_c : 0.0;
                double u_c = uh_c * inv_h_c;
                double v_c = vh_c * inv_h_c;

                double dw0 = stage_ev[k3 + 0] - w_c;
                double dw1 = stage_ev[k3 + 1] - w_c;
                double wx  = inv_det * (dyv1 * dw0 - dyv0 * dw1);
                double wy  = inv_det * (dxv0 * dw1 - dxv1 * dw0);

                double dh0 = height_ev[k3 + 0] - h_c;
                double dh1 = height_ev[k3 + 1] - h_c;
                double hx  = inv_det * (dyv1 * dh0 - dyv0 * dh1);
                double hy  = inv_det * (dxv0 * dh1 - dxv1 * dh0);

                double h_e0     = height_ev[k3 + 0];
                double h_e1     = height_ev[k3 + 1];
                double inv_h_e0 = (h_e0 > minimum_allowed_height) ? 1.0 / h_e0 : 0.0;
                double inv_h_e1 = (h_e1 > minimum_allowed_height) ? 1.0 / h_e1 : 0.0;
                double u_e0 = xmom_ev[k3 + 0] * inv_h_e0;
                double u_e1 = xmom_ev[k3 + 1] * inv_h_e1;
                double v_e0 = ymom_ev[k3 + 0] * inv_h_e0;
                double v_e1 = ymom_ev[k3 + 1] * inv_h_e1;

                double du0 = u_e0 - u_c;
                double du1 = u_e1 - u_c;
                double dv0 = v_e0 - v_c;
                double dv1 = v_e1 - v_c;
                double ux  = inv_det * (dyv1 * du0 - dyv0 * du1);
                double uy  = inv_det * (dxv0 * du1 - dxv1 * du0);
                double vx  = inv_det * (dyv1 * dv0 - dyv0 * dv1);
                double vy  = inv_det * (dxv0 * dv1 - dxv1 * dv0);

                double g_h = g_pred * h_c;
                double dw_dt  = -(u_c * hx + h_c * ux + v_c * hy + h_c * vy);
                double duh_dt = -(2.0*u_c*h_c*ux + u_c*u_c*hx + u_c*v_c*hy
                                 + v_c*h_c*uy + u_c*h_c*vy + g_h * wx);
                double dvh_dt = -(v_c*h_c*ux + u_c*h_c*vx + u_c*v_c*hx
                                 + 2.0*v_c*h_c*vy + v_c*v_c*hy + g_h * wy);

                // NOTE: bed_ev is NOT refreshed after the shift, exactly like
                // the standalone predictor: stage and height shift by the same
                // dw_dt, so stage - height still equals the true bed everywhere
                // except clamped near-dry edges -- and the pre-shift bed_ev
                // (the true bed) is what the boundary kernels should read.
                // Manning drag on the half-advanced state; see
                // core_ader_friction_factor().  The bed gradient is
                // grad(w) - grad(h), already in registers.
                double zs_p = 1.0;
                if (sloped_mannings) {
                    double zx_p = wx - hx;               // w = z + h
                    double zy_p = wy - hy;
                    zs_p = sqrt(1.0 + zx_p * zx_p + zy_p * zy_p);
                }
                double eta_p = friction_cv[k];

                for (int i = 0; i < 3; i++) {
                    // gamma from the PRE-shift state: that is what makes the
                    // step exact on steady uniform flow.  See the derivation
                    // above core_ader_friction_factor().
                    double fac = core_ader_friction_factor(
                        g_pred, eta_p, zs_p, minimum_allowed_height,
                        ANUGA_SEVEN_THIRDS, height_ev[k3 + i],
                        xmom_ev[k3 + i], ymom_ev[k3 + i], predictor_dt);

                    stage_ev[k3 + i] += predictor_dt * dw_dt;
                    double uh_e = xmom_ev[k3 + i] + predictor_dt * duh_dt;
                    double vh_e = ymom_ev[k3 + i] + predictor_dt * dvh_dt;
                    height_ev[k3 + i] = fmax(height_ev[k3 + i]
                                             + predictor_dt * dw_dt, 0.0);

                    xmom_ev[k3 + i] = uh_e * fac;
                    ymom_ev[k3 + i] = vh_e * fac;
                }
            }
        }

    }

}

void core_extrapolate_edge_pass(struct domain *D, double predictor_dt) {
    core_extrapolate_edge_pass_on(D, predictor_dt, NULL, 0);
}

void core_extrapolate_second_order_edge(struct domain *D) {
    // Kept as the two passes below so callers that can fuse the (cell-local)
    // centroid pass into a neighbouring kernel may call the edge pass alone.
    core_extrapolate_centroid_pass(D);
    core_extrapolate_edge_pass(D, 0.0);
}

// ============================================================================
// Distribute edge values to vertices
// ============================================================================

void core_distribute_edges_to_vertices(struct domain *D) {
    anuga_int n = D->number_of_elements;

    double * restrict stage_ev = D->stage_edge_values;
    double * restrict xmom_ev = D->xmom_edge_values;
    double * restrict ymom_ev = D->ymom_edge_values;
    double * restrict bed_ev = D->bed_edge_values;
    double * restrict height_ev = D->height_edge_values;

    double * restrict stage_vv = D->stage_vertex_values;
    double * restrict xmom_vv = D->xmom_vertex_values;
    double * restrict ymom_vv = D->ymom_vertex_values;
    double * restrict bed_vv = D->bed_vertex_values;
    double * restrict height_vv = D->height_vertex_values;

    OMP_PARALLEL_LOOP
    for (anuga_int k = 0; k < n; k++) {
        anuga_int k3 = k * 3;

        // Reconstruct vertex values from edge values
        // vertex[i] = edge[i+1] + edge[i+2] - edge[i]
        stage_vv[k3 + 0] = stage_ev[k3 + 1] + stage_ev[k3 + 2] - stage_ev[k3 + 0];
        stage_vv[k3 + 1] = stage_ev[k3 + 2] + stage_ev[k3 + 0] - stage_ev[k3 + 1];
        stage_vv[k3 + 2] = stage_ev[k3 + 0] + stage_ev[k3 + 1] - stage_ev[k3 + 2];

        xmom_vv[k3 + 0] = xmom_ev[k3 + 1] + xmom_ev[k3 + 2] - xmom_ev[k3 + 0];
        xmom_vv[k3 + 1] = xmom_ev[k3 + 2] + xmom_ev[k3 + 0] - xmom_ev[k3 + 1];
        xmom_vv[k3 + 2] = xmom_ev[k3 + 0] + xmom_ev[k3 + 1] - xmom_ev[k3 + 2];

        ymom_vv[k3 + 0] = ymom_ev[k3 + 1] + ymom_ev[k3 + 2] - ymom_ev[k3 + 0];
        ymom_vv[k3 + 1] = ymom_ev[k3 + 2] + ymom_ev[k3 + 0] - ymom_ev[k3 + 1];
        ymom_vv[k3 + 2] = ymom_ev[k3 + 0] + ymom_ev[k3 + 1] - ymom_ev[k3 + 2];

        bed_vv[k3 + 0] = bed_ev[k3 + 1] + bed_ev[k3 + 2] - bed_ev[k3 + 0];
        bed_vv[k3 + 1] = bed_ev[k3 + 2] + bed_ev[k3 + 0] - bed_ev[k3 + 1];
        bed_vv[k3 + 2] = bed_ev[k3 + 0] + bed_ev[k3 + 1] - bed_ev[k3 + 2];

        height_vv[k3 + 0] = height_ev[k3 + 1] + height_ev[k3 + 2] - height_ev[k3 + 0];
        height_vv[k3 + 1] = height_ev[k3 + 2] + height_ev[k3 + 0] - height_ev[k3 + 1];
        height_vv[k3 + 2] = height_ev[k3 + 0] + height_ev[k3 + 1] - height_ev[k3 + 2];
    }
}



// ============================================================================
// Near-bed concentration ratio d*(Z)   -- spec 4.3, open item S1a
// ============================================================================
//
// [S-4] is a 1-D quadrature per cell, far too expensive to run inside a kernel
// every step, so the spec calls for a fitted form. This is that fit.
//
// CORRECTION TO [S-4] -- THE TYPO IS IN DL09 AS PUBLISHED, not in PHYSICS_SPEC,
// which transcribed them faithfully. DL09 Eq 19 gives the Rouse profile factor
// as ((z-a)/(h-a) . a/z)^Z, which is ZERO at the reference height z = a.
//
// Their own paper disproves it: immediately above Eq 19 they write the flux as
// q_S = c_S(a) * integral( [...]^Z u(z) dz ). Factoring c_S(a) out REQUIRES the
// bracket to be 1 at z = a; the printed factor is 0 there, giving c_s(a) = 0.
//
// The Rouse-Vanoni profile [(h-z)/z . a/(h-a)]^Z rearranges exactly to
// ((h-z)/(h-a) . a/z)^Z, which is 1 at z = a and 0 at z = h. One glyph: (z-a)
// is a slip for (h-z). At Z = 2, a/h = 0.005 the printed form gives d* = 41227
// against 356 corrected -- not a value that appears on their Figure 4, so that
// figure was evidently computed with the correct profile.
//
// The fit below is to the CORRECTED integral. See PHYSICS_SPEC 4.3, Draft 5.
//
// FITTED FORM. Near the bed the Rouse profile behaves like z^(-Z), so d*
// diverges roughly as (a/h)^(-Z). Factoring that out first,
//
//     ln d* = -Z ln(a/h) + P(Z, ln(a/h))
//
// leaves a mild remainder a low-order polynomial captures well: 28 terms for
// 0.82 percent max / 0.05 percent mean error. A direct polynomial in
// (ln Z, ln(a/h)) needs far more terms for far worse accuracy. The structure is
// what buys it.
//
// RANGE. Fitted for Z in [0.01, 2.5] and a/h in [1e-3, 0.15]. Beyond Z ~ 2.5
// transport is essentially bedload and this ratio is not the right model. The
// a/h range reaches down to 1e-3 deliberately: the shipped anugaSed operator
// implies a/h ~ 9.3e-4 (spec 12, D4b), so its regime is reachable when the
// a >= floor*h floor is relaxed. (9.3e-4 itself clamps to 1e-3, ~7 percent in
// a/h and so ~7 percent in d* at Z = 1.)
//
// OUT-OF-RANGE INPUTS ARE CLAMPED, NOT EXTRAPOLATED -- the spec asks for that
// explicitly, because a polynomial taken outside its fitted range goes wrong
// quietly. anugaSed's own 8th-degree fit is extrapolated freely and reaches
// p(6) = 282088 (spec 12, D4c); this one cannot.
#define ANUGA_ROUSE_Z_LO   0.01
#define ANUGA_ROUSE_Z_HI   2.5
#define ANUGA_ROUSE_AH_LO  1e-3
#define ANUGA_ROUSE_AH_HI  0.15

static inline double core_rouse_d_star(double Z, double a_h) {
    /* P(Z, L) = sum_i Z^i (c_i0 + c_i1 L + c_i2 L^2 + c_i3 L^3), L = ln(a/h) */
    const double C[7][4] = {
    {+1.097192252266e-03, +9.816426876103e-04, +2.816550608693e-04, +2.216981593577e-05},
    {+8.152552738643e-01, +2.984288438662e-01, +4.717126357513e-02, +2.488390592718e-03},
    {-3.858022145865e-02, +7.497943739332e-01, +1.494530687071e-01, +1.016829763599e-02},
    {-1.416163484237e-01, -6.145585548869e-01, -2.181478641118e-01, -1.989085090511e-02},
    {+2.441798567588e-02, +2.477861105262e-01, +1.172562489413e-01, +1.380260943049e-02},
    {+1.714535144604e-02, -4.453006825886e-02, -2.922351673152e-02, -4.300807416335e-03},
    {-4.557991043496e-03, +2.511784955218e-03, +2.810236330103e-03, +5.048472261972e-04},
    };

    if (Z < ANUGA_ROUSE_Z_LO) Z = ANUGA_ROUSE_Z_LO;
    else if (Z > ANUGA_ROUSE_Z_HI) Z = ANUGA_ROUSE_Z_HI;
    if (a_h < ANUGA_ROUSE_AH_LO) a_h = ANUGA_ROUSE_AH_LO;
    else if (a_h > ANUGA_ROUSE_AH_HI) a_h = ANUGA_ROUSE_AH_HI;

    const double L = log(a_h);
    const double L2 = L * L;

    /* Horner in Z over coefficients that are cubics in L. */
    const double L3 = L2 * L;
    double P = 0.0;
    for (int i = 6; i >= 0; i--) {
        P = P * Z + (C[i][0] + C[i][1] * L + C[i][2] * L2 + C[i][3] * L3);
    }

    const double d = exp(-Z * L + P);
    /* DL09: d* is always > 1. Guard the fit against dipping below it. */
    return (d < 1.0) ? 1.0 : d;
}

// Centroid height of the Rouse profile, for the settling time of [D-4].
//
// (z_c - a)/h with z_c = int z c dz / int c dz over [a, h] of the corrected
// Rouse-Vanoni profile above. A collapsing profile deposits its sediment
// from about its centroid, so the settling time is (z_c - a)/w_s: h/2 over
// w_s when well mixed (Z -> 0, z_c -> h/2), vanishing when the load already
// hugs the bed. Same fitted form and clamps as d*: ln((z_c - a)/h) =
// P(Z, ln(a/h)), 28 terms, 1.0 percent max / 0.27 percent mean error over
// Z in [0.01, 2.5], a/h in [1e-3, 0.15] (least squares to quadrature,
// 2026-09-24).
static inline double core_rouse_z_c(double Z, double a_h) {
    const double C[7][4] = {
    {-1.236620662919e+00, -2.966806356338e-01, -5.451826597743e-02, -3.340785646720e-03},
    {-4.934470965259e-01, +1.986068502121e-01, +2.681666552524e-02, +1.430964808434e-03},
    {+5.693206124955e-01, +4.438916973539e-01, +7.169076944620e-02, +3.074161648449e-03},
    {-8.317934945723e-01, -8.674254433522e-01, -2.464308678880e-01, -1.711508840976e-02},
    {+5.466837296437e-01, +6.103587559357e-01, +2.038814693879e-01, +1.776942243679e-02},
    {-1.671607038599e-01, -1.948379068715e-01, -7.092485415218e-02, -7.151638532153e-03},
    {+1.963501463003e-02, +2.364362131184e-02, +9.087924042213e-03, +1.006461512046e-03},
    };
    if (Z < ANUGA_ROUSE_Z_LO) Z = ANUGA_ROUSE_Z_LO;
    else if (Z > ANUGA_ROUSE_Z_HI) Z = ANUGA_ROUSE_Z_HI;
    if (a_h < ANUGA_ROUSE_AH_LO) a_h = ANUGA_ROUSE_AH_LO;
    else if (a_h > ANUGA_ROUSE_AH_HI) a_h = ANUGA_ROUSE_AH_HI;
    const double L = log(a_h), L2 = L * L, L3 = L2 * L;
    double P = 0.0;
    for (int i = 6; i >= 0; i--) {
        P = P * Z + (C[i][0] + C[i][1] * L + C[i][2] * L2 + C[i][3] * L3);
    }
    return exp(P);
}

/* Slope magnitude of a centroid field at cell k for [T-7]/[T-7e]: the
 * least-squares gradient over the cell and its neighbours, one-sided at
 * boundaries (the gradient along the line when the neighbours are
 * collinear or there is only one). A plane gives its slope in every cell,
 * walls and corners included.
 *
 * Centroid values, not edge values. The DE algorithms rebuild the bed and
 * stage edge values every step as (limited stage) - (limited height), so a
 * slope read from the edges is whatever the hydrodynamic limiter left: zero
 * in the cells along a reflective wall (the ghost mirrors the interior, so
 * the reconstruction sees no gradient toward it) and reduced wherever the
 * limiter engages -- kinks, steps, scarps, exactly where erosion matters.
 * The wall cells of a uniformly sloping bed then never eroded while their
 * neighbours did. The centroid gradient is the actual slope of the field
 * the solver carries, tracks an evolving bed, and at a wall the internal
 * neighbours still determine it.
 *
 * READS ITS NEIGHBOURS, so it must not run inside a loop that writes the
 * field: see sediment_slope_work in sw_domain.h. */
// Friction factor the sediment kernels see in a VEGETATED cell (spec 8 with
// spec 3.3). Baptist et al. (2007) split the resistance: the stems take the
// drag, and the bed feels the flow inside the canopy on the bed roughness Cb,
//
//   u_v = U Cv_r / Cv,   tau_b / rho = g u_v^2 / Cb^2        (veg_shear 1, bed)
//
// with Cv_r = (Cb^-2 + Cd m D min(h, hv) / (2 g))^-1/2 the canopy part of the
// vegetated Chezy coefficient Cv (Cv = Cv_r when the stems are emergent). The
// alternative is the whole vegetated resistance, f_c = g / Cv^2 (veg_shear 2,
// total), which is what the Delta-X Wax Lake sediment model uses. Returns -1
// for a cell without vegetation, or with veg_shear 0, so the caller keeps
// its own closure there. Expressed as f_c in tau_b / rho = f_c U^2.
static inline double core_vegetation_sediment_fc(anuga_int veg_shear, double h,
                                                 double m, double Dstem, double hv,
                                                 double Cd, double Cb, double g) {
    if (veg_shear == 0 || !(m > 0.0) || !(hv > 0.0) || !(h > 0.0)) return -1.0;
    const double Cv_r = 1.0 / sqrt(1.0 / (Cb * Cb) + Dstem * m * Cd / (2.0 * g) * fmin(h, hv));
    double Cv = Cv_r;
    if (h > hv) Cv += sqrt(g) / 0.4 * log(h / hv);
    if (veg_shear == 2) return g / (Cv * Cv);
    const double r = Cv_r / Cv;
    return g * r * r / (Cb * Cb);
}

static inline double core_centroid_slope(const double * restrict f,
                                         const anuga_int * restrict neighbours,
                                         const anuga_geom_t * restrict cc,
                                         anuga_int k) {
    const double fk = f[k];
    const double xk = cc[2 * k];
    const double yk = cc[2 * k + 1];
    double sxx = 0.0, sxy = 0.0, syy = 0.0, sxf = 0.0, syf = 0.0;
    for (anuga_int i = 0; i < 3; i++) {
        const anuga_int n = neighbours[3 * k + i];
        if (n < 0) {
            continue;                            /* boundary edge: one-sided */
        }
        const double dx = cc[2 * n] - xk;
        const double dy = cc[2 * n + 1] - yk;
        const double df = f[n] - fk;
        sxx += dx * dx;
        sxy += dx * dy;
        syy += dy * dy;
        sxf += dx * df;
        syf += dy * df;
    }
    double gx, gy;
    const double det = sxx * syy - sxy * sxy;
    const double trace = sxx + syy;
    if (det > 1.0e-12 * trace * trace) {
        gx = (syy * sxf - sxy * syf) / det;      /* the least-squares plane */
        gy = (sxx * syf - sxy * sxf) / det;
    } else if (trace > 0.0) {
        gx = sxf / trace;                        /* along the one line the */
        gy = syf / trace;                        /* data determine */
    } else {
        gx = 0.0;                                /* isolated cell */
        gy = 0.0;
    }
    return sqrt(gx * gx + gy * gy);
}

/* tau_b / rho for a cell, under either shear closure (spec 3.1 / 3.4).
 *
 *   [T-1]  tau_b/rho = f_c |v|^2          quadratic drag, no equilibrium
 *                                         assumption
 *   [T-7]  tau_b/rho = g h S              depth-slope, aSM16 Eqs 6-7
 *
 * Returning tau_b/rho rather than tau_b keeps the two interchangeable
 * everywhere downstream: the Shields stress is tau* = (tau_b/rho)/(R g d), in
 * which the density cancels, and the dimensional stress the cohesive route
 * needs is simply rho_w times this.
 *
 * S is core_centroid_slope of the bed for [T-7] and of the free surface for
 * [T-7e], computed by the caller beforehand. The two agree only in steady
 * uniform flow, which is [T-7]'s stated assumption; where the flow is not in
 * equilibrium the free-surface slope is the better estimate of what drives
 * it, and it is what the older Bed_shear_erosion_operator used. */
static inline double core_tau_b_over_rho(anuga_int closure, double f_c,
                                         double vel2, double grav, double h,
                                         double S) {
    if (closure != 1 && closure != 2) {
        return f_c * vel2;                       /* [T-1] */
    }
    return grav * h * S;                         /* [T-7] / [T-7e] */
}

// ============================================================================
// Bedload transport and its bed evolution   [G-5]/[K-3], spec 6
// ============================================================================
//
//   [K-1]  q_b* = K tau_x^m                 tau_x = tau* - tau_c*  [T-4]
//   [K-5]  q_b* = 0.05 tau*^2.5 / f_c       Engelund-Hansen, no threshold
//   [K-2]  q_b  = q_b* sqrt(R g) d^1.5
//   [K-6]  q_b  = A_g |u|^m                 Grass, total load, no threshold;
//                                           dimensional as it stands, so [K-2]
//                                           does not apply and the classes do
//                                           not enter
//   [K-4]  q_b is parallel to the bed shear stress, hence to (u, v)
//   [K-3]  dz/dt = -(1/(1-lambda)) div q_b
//
// Unlike the suspended exchange, this is a DIVERGENCE: it moves sediment along
// the bed rather than between bed and water column, so in a closed domain it
// redistributes bed material and conserves total bed volume exactly. That is
// the property to test it with.
//
// Three passes, because the divergence at cell k needs its neighbours' q_b
// and, through the Rusanov term, their bed elevation: pass 1 fills the
// per-cell transport vector, pass 2 takes the divergence into a scratch,
// pass 3 applies it. The bed is not written until every cell has read it,
// so the two sides of an edge always see the same state and the flux stays
// antisymmetric. All are ordinary cell loops, so all offload.
//
// Edge flux is CENTRED WITH RUSANOV DISSIPATION; see the note at the flux
// itself for why plain upwinding was tried and rejected, and why the bare
// centred flux was not enough. Boundary edges carry zero bedload flux, which is what
// makes the closed-domain conservation test exact -- unless the edge has been
// declared OPEN (sediment_bedload_open, set per boundary tag from Python), in
// which case the flux through it is the cell's own q_b.n: zero gradient, so
// bedload leaves through an outflow at the rate it arrives and enters through
// an inflow at the rate the first cell carries it away (equilibrium supply).
// Without that an inflow cell exports and never imports, and digs a hole that
// travels downstream at the bed-wave speed.
void core_apply_bedload(struct domain *D, double timestep) {
    const anuga_int mode = D->sediment_bedload_mode;
    const anuga_int n_classes = D->n_sediment_classes;
    if (mode == 0 || n_classes <= 0 || timestep <= 0.0) {
        return;
    }
    if (!D->sediment_bed_evolution) {
        return;   /* fixed bed: bedload would have nowhere to go */
    }

    const anuga_int n = D->number_of_elements;
    const double grav = D->g;
    const double h_eps = D->epsilon;
    const double minimum_allowed_height = D->minimum_allowed_height;
    const double one_minus_lambda = 1.0 - D->sediment_porosity;
    const double morfac = D->sediment_morphological_factor;
    const double K = D->sediment_bedload_K;
    const double mexp = D->sediment_bedload_m;
    const double tau_c_b = D->sediment_bedload_tau_c_star;
    const double bl_h_min = D->sediment_bedload_h_min;
    const anuga_int fric_mode = D->sediment_friction_mode;
    const double n_ll = D->sediment_manning_ll;
    const anuga_int wbed = D->sediment_wilson_bed;
    const double wD = D->sediment_wilson_D;
    // Vegetation seen by the sediment shear (core_vegetation_sediment_fc).
    const anuga_int veg_shear = (D->vegetation_mode > 0) ? D->sediment_vegetation_shear : 0;
    double * restrict veg_m = D->veg_density_centroid_values;
    double * restrict shear_fac = D->sediment_shear_factor;   /* [T-16], or NULL */
    double * restrict veg_d = D->veg_diameter_centroid_values;
    double * restrict veg_h = D->veg_height_centroid_values;
    const double veg_Cd = D->vegetation_Cd;
    const double veg_Cb = D->vegetation_bed_chezy;
    const anuga_int shear_closure = D->sediment_shear_closure;
    const double max_slope = D->sediment_max_slope;
    const anuga_int slope_frozen = D->sediment_slope_frozen;

    double * restrict stage_cv = D->stage_centroid_values;
    double * restrict bed_cv = D->bed_centroid_values;
    double * restrict bed_ev = D->bed_edge_values;
    double * restrict xmom_cv = D->xmom_centroid_values;
    double * restrict ymom_cv = D->ymom_centroid_values;
    double * restrict friction_cv = D->friction_centroid_values;
    double * restrict diam = D->sediment_diameter;
    double * restrict sedR = D->sediment_R;
    double * restrict qbx = D->sediment_qbx;
    double * restrict qby = D->sediment_qby;
    /* Rusanov coefficient per cell: |d q_b / d z| / |q_b| = gamma / h, so
     * that gamma/h * |q_b . n| bounds the bed-wave speed across the edge. */
    double * restrict qba = D->sediment_qba;
    /* The bed change of this step, held until every cell has read the bed.
     * The repose kernel's Jacobi scratch, free between its calls. */
    double * restrict dzs = D->sediment_repose_dz;
    anuga_int * restrict neighbours = D->neighbours;
    anuga_geom_t * restrict cc = D->centroid_coordinates;
    double * restrict slope_w = D->sediment_slope_work;
    anuga_geom_t * restrict normals = D->normals;
    anuga_geom_t * restrict edgelengths = D->edgelengths;
    anuga_geom_t * restrict areas = D->areas;
    /* [L-5]. Both the flag and the arrays are required; see the note in
     * core_apply_sediment_source. */
    double * restrict z_base = D->sediment_z_base;
    anuga_int * restrict exhausted = D->sediment_bed_exhausted;
    const anuga_int has_z_base = (D->sediment_has_z_base
                                  && z_base != NULL && exhausted != NULL);
    /* Per boundary edge, indexed by -neighbour - 1: 1 where bedload passes. */
    anuga_int * restrict bopen = D->sediment_bedload_open;
    /* Per boundary edge: a prescribed bedload inflow (m^2/s), or < 0. */
    double * restrict bsup = D->sediment_bedload_supply;

    if (one_minus_lambda <= 0.0) {
        return;
    }

    /* ---- pass 1: the transport vector, summed over classes ---- */
    OMP_PARALLEL_LOOP
    for (anuga_int k = 0; k < n; k++) {
        qbx[k] = 0.0;
        qby[k] = 0.0;
        qba[k] = 0.0;

        const double h = fmax(stage_cv[k] - bed_cv[k], 0.0);
        if (h <= minimum_allowed_height) {
            continue;
        }
        /* [K-7] shallow-film ramp: none below h_min, full above 2 h_min. */
        double h_ramp = 1.0;
        if (bl_h_min > 0.0) {
            h_ramp = (h - bl_h_min) / bl_h_min;
            if (h_ramp <= 0.0) continue;
            if (h_ramp > 1.0) h_ramp = 1.0;
        }

        const double denom = h * h + h_eps * h_eps;
        const double u = (denom > 0.0) ? (xmom_cv[k] * h / denom) : 0.0;
        const double v = (denom > 0.0) ? (ymom_cv[k] * h / denom) : 0.0;
        const double vel2 = u * u + v * v;
        if (vel2 <= 0.0) {
            continue;
        }
        const double speed = sqrt(vel2);

        double f_c;
        if (fric_mode == 2) {
            double rel = h / wD;
            if (!(rel > 1.0)) rel = 1.0;
            double X;
            if (wbed == 0)      X = 8.46 * pow(rel, 0.1005);
            else if (wbed == 1) X = 5.75 * log10(rel) + 3.514;
            else                X = 5.62 * log10(rel) + 4.0;
            f_c = 1.0 / (X * X);
        } else {
            const double nman = (fric_mode == 1) ? n_ll : friction_cv[k];
            f_c = grav * nman * nman / cbrt(h);
        }
        if (veg_shear) {
            const double f_v = core_vegetation_sediment_fc(veg_shear, h, veg_m[k], veg_d[k],
                                                           veg_h[k], veg_Cd, veg_Cb, grav);
            if (f_v >= 0.0) f_c = f_v;
        }
        if (shear_fac != NULL) f_c *= shear_fac[k];           /* [T-16] */

        /* Same closure as the suspended source: [T-1], [T-7] or [T-7e].
         * This pass only reads the bed, so the slope can be taken here. */
        double S = (shear_closure == 1)
                     ? (slope_frozen ? slope_w[k]
                                     : core_centroid_slope(bed_cv, neighbours, cc, k))
                 : (shear_closure == 2) ? core_centroid_slope(stage_cv, neighbours, cc, k)
                 : 0.0;
        if (max_slope > 0.0 && S > max_slope) S = max_slope;
        const double tbr = core_tau_b_over_rho(shear_closure, f_c, vel2, grav,
                                               h, S);

        /* q_b_total is the magnitude; dq_total is h |dq_b/dz| at fixed q
         * and w, which the Rusanov term below needs. With u = q/h: Grass
         * q_b ~ u^m gives m q_b; the power law in the Shields stress, which
         * goes as f_c u^2 ~ h^(-7/3) under Manning, gives m (7/3) (tau_star
         * over tau_x) q_b per class; Engelund-Hansen, tau_star^2.5 over f_c
         * ~ h^(-5.5), gives 5.5 q_b. Wilson's f_c varies less with h than
         * Manning's, so 7/3 over-estimates there, which only adds
         * dissipation. */
        double q_b_total = 0.0;
        double dq_total = 0.0;
        if (mode == 3) {
            /* [K-6] Grass. A_g (in K) and the exponent (in m) are the whole
             * law: no stress, no threshold, no grain size. */
            q_b_total = K * pow(speed, mexp);
            dq_total = mexp * q_b_total;
        }
        for (anuga_int s = 0; s < n_classes && mode != 3; s++) {
            const double Rgd = sedR[s] * grav * diam[s];
            if (!(Rgd > 0.0)) continue;
            const double tau_star = tbr / Rgd;

            double q_star;
            double dfac;
            if (mode == 2) {
                /* [K-5] Engelund-Hansen, total load, NO threshold. Subtracting
                 * tau_c* here would silently make it a different model. */
                q_star = (f_c > 0.0) ? 0.05 * pow(tau_star, 2.5) / f_c : 0.0;
                dfac = 5.5;
            } else {
                const double tau_x = tau_star - tau_c_b;
                q_star = (tau_x > 0.0) ? K * pow(tau_x, mexp) : 0.0;
                dfac = (tau_x > 0.0) ? mexp * (7.0 / 3.0) * tau_star / tau_x
                                     : 0.0;
            }
            if (q_star <= 0.0) continue;

            /* [K-2] */
            const double q_bs = q_star * sqrt(sedR[s] * grav) * pow(diam[s], 1.5);
            q_b_total += q_bs;
            dq_total += dfac * q_bs;
        }

        q_b_total *= h_ramp;                     /* [K-7]; dq scales with it */
        dq_total *= h_ramp;
        if (q_b_total > 0.0) {
            /* The Rusanov coefficient, gamma / h. Formed before the [L-5]
             * cap below, which scales q_b and dq alike and can take both to
             * zero on an exhausted cell. */
            qba[k] = dq_total / (q_b_total * h);

            /* [L-5]. A cell cannot export bed material it does not have.
             * The limit is applied to the TRANSPORT VECTOR, not to the
             * divergence: both cells sharing an edge then form their flux
             * from the same limited q, so the flux stays antisymmetric and
             * the scheme stays exactly conservative. Clamping the divergence
             * instead would let one side remove what the other never
             * received, which creates bed material. */
            if (has_z_base) {
                const double avail = bed_cv[k] - z_base[k];
                const double thickness = (avail > 0.0) ? avail : 0.0;
                /* The cell's own contribution to its outflow: its half of
                 * every edge's centred flux, counting only the outgoing
                 * ones. */
                double own_out = 0.0;
                const double ex = q_b_total * u / speed;
                const double ey = q_b_total * v / speed;
                for (anuga_int i = 0; i < 3; i++) {
                    const anuga_int ki = 3 * k + i;
                    const anuga_int nbk = neighbours[ki];
                    double qn = ex * normals[6 * k + 2 * i]
                              + ey * normals[6 * k + 2 * i + 1];
                    if (nbk < 0) {
                        /* An open boundary edge carries the whole of the
                         * cell's own flux; a closed one carries none; one
                         * with a prescribed supply is an inflow. */
                        if (bopen == NULL || !bopen[-nbk - 1]) continue;
                        if (bsup != NULL && bsup[-nbk - 1] >= 0.0) continue;
                    } else {
                        qn *= 0.5;
                    }
                    if (qn > 0.0) own_out += qn * edgelengths[ki];
                }
                if (own_out > 0.0) {
                    /* over the step the bed moves morfac times the
                     * physical change, so the cap shrinks by morfac */
                    const double cap = thickness * one_minus_lambda
                                     * areas[k] / (timestep * morfac);
                    if (own_out > cap) {
                        q_b_total *= cap / own_out;
                    }
                }
            }

            /* [K-4] parallel to (u, v) */
            qbx[k] = q_b_total * u / speed;
            qby[k] = q_b_total * v / speed;
        }
    }

    /* ---- [L-5] pass 1.5: which cells cannot afford what is about to be
     * taken from them ----------------------------------------------------
     *
     * Flagging cells that have ALREADY reached the base is not enough. A cell
     * with a millimetre left can be asked for two in a single step and is
     * only found empty afterwards, which is how the first version of this
     * overshot its base by 6.1e-6 m on a 1e-2 m layer. So the test is
     * predictive: form the divergence this step WILL produce and flag the
     * cell if it cannot pay for it.
     *
     * Separate sweep, not folded into pass 2, because pass 2 writes bed
     * elevation while its neighbours are still reading it -- see the note in
     * sw_domain.h. This one only reads, so every cell sees the same state.
     */
    if (has_z_base) {
        OMP_PARALLEL_LOOP
        for (anuga_int k = 0; k < n; k++) {
            const double avail = bed_cv[k] - z_base[k];
            if (avail <= 0.0) {
                exhausted[k] = 1;      /* nothing left to give */
                continue;
            }
            double outflux = 0.0;
            for (anuga_int i = 0; i < 3; i++) {
                const anuga_int ki = 3 * k + i;
                const anuga_int nb = neighbours[ki];
                double qn;
                if (nb < 0) {
                    if (bopen == NULL || !bopen[-nb - 1]) continue;
                    qn = qbx[k] * normals[6 * k + 2 * i]
                       + qby[k] * normals[6 * k + 2 * i + 1];
                } else {
                    const double nx = normals[6 * k + 2 * i];
                    const double ny = normals[6 * k + 2 * i + 1];
                    const double qk = qbx[k] * nx + qby[k] * ny;
                    const double qb = qbx[nb] * nx + qby[nb] * ny;
                    const double a = fmax(qba[k] * fabs(qk), qba[nb] * fabs(qb));
                    qn = 0.5 * (qk + qb) - 0.5 * a * (bed_cv[nb] - bed_cv[k]);
                }
                outflux += qn * edgelengths[ki];
            }
            /* dz this step, if nothing were blocked. */
            const double drop = (timestep * outflux / areas[k])
                              / one_minus_lambda;
            exhausted[k] = (drop > avail) ? 1 : 0;
        }
    }

    /* ---- pass 2: divergence, into the scratch ---- */
    OMP_PARALLEL_LOOP
    for (anuga_int k = 0; k < n; k++) {
        const double qx_k = qbx[k];
        const double qy_k = qby[k];
        double outflux = 0.0;

        for (anuga_int i = 0; i < 3; i++) {
            const anuga_int ki = 3 * k + i;
            const anuga_int nb = neighbours[ki];
            const double nx = normals[6 * k + 2 * i];
            const double ny = normals[6 * k + 2 * i + 1];
            if (nb < 0) {
                /* Boundary: no bedload across it unless declared open, in
                 * which case the flux is the cell's own q_b.n (zero
                 * gradient). Only k's own exhaustion can close it. */
                if (bopen == NULL || !bopen[-nb - 1]) continue;
                double qn;
                if (bsup != NULL && bsup[-nb - 1] >= 0.0) {
                    /* Prescribed supply: bedload enters at this rate
                     * whatever the cell's own transport, so the import
                     * cannot feed back on the cell's state. */
                    qn = -bsup[-nb - 1];
                } else {
                    qn = qx_k * nx + qy_k * ny;
                    if (has_z_base && qn > 0.0 && exhausted[k]) qn = 0.0;
                }
                outflux += qn * edgelengths[ki];
                continue;
            }

            /* CENTRED edge flux: q_edge = (q_k + q_nb)/2.
             *
             * Exactly conservative -- cell nb forms the same average against
             * the opposite normal and so removes precisely what k gains -- and,
             * unlike an upwind donor choice, CONTINUOUS.
             *
             * Upwinding was tried first and rejected twice over. Deciding the
             * donor from each cell's own q.n is not antisymmetric and creates
             * bed material (measured 1.05e-5 of bed volume in 60 s). Fixing
             * that by deciding from the averaged vector is conservative but
             * DISCONTINUOUS: where q_k.n = -q_nb.n, which is exactly a
             * converging-flow edge, the average passes through zero while the
             * two candidate fluxes differ by a finite amount, so the donor
             * flips on roundoff. That put mode 1 and mode 2 1.3e-4 apart.
             *
             * The centred flux is also the physically right answer at such an
             * edge: bedload converging from both sides should deposit there,
             * not be attributed to one arbitrary donor.
             *
             * They did appear, in the migrating-hump validation case
             * (sediment_bed_hump): under the bare centred flux the hump grew
             * a scour hole at its upstream toe and an overshoot at its crest
             * that fed back into the flow. So the flux is Rusanov: centred
             * plus a dissipation term in (z_nb - z_k) scaled by a bound on
             * the bed-wave speed, |d(q_b . n)/dz| = (gamma/h) |q_b . n| from
             * pass 1, the larger of the two sides. It is still antisymmetric
             * (both sides form the same a and the same difference with the
             * opposite sign), so still exactly conservative, and continuous
             * in the state. The dissipation vanishes where the bed is flat
             * and where nothing moves, so the closed-domain and open-uniform
             * tests are unchanged by it. */
            const double qn_k = qx_k * nx + qy_k * ny;
            const double qn_nb = qbx[nb] * nx + qby[nb] * ny;
            const double a = fmax(qba[k] * fabs(qn_k), qba[nb] * fabs(qn_nb));
            double qn = 0.5 * (qn_k + qn_nb) - 0.5 * a * (bed_cv[nb] - bed_cv[k]);

            /* [L-5]. A cell that cannot pay for this step's removal (pass
             * 1.5) may gain material but must not lose any, so close every
             * edge that would take material OUT of it. qn > 0 is outflow from
             * k, qn < 0 is outflow from nb.
             *
             * This is SYMMETRIC, which is the whole point: cell nb reaches
             * this edge with the opposite normal, hence -qn, and the same
             * two tests in the same order, so both sides close the same edge
             * and neither can remove what the other did not give up. That is
             * what keeps bedload exactly conservative with a base present.
             * It works on a snapshot taken in pass 1.5 so that both sides
             * decide from the same state. */
            if (has_z_base) {
                if (qn > 0.0 && exhausted[k])  qn = 0.0;
                if (qn < 0.0 && exhausted[nb]) qn = 0.0;
            }
            outflux += qn * edgelengths[ki];
        }

        /* [K-3]: dz/dt = -(1/(1-lambda)) div q_b, div q_b = outflux/area;
         * times the morphological factor */
        dzs[k] = -(morfac * timestep * outflux / areas[k]) / one_minus_lambda;
    }

    /* ---- pass 3: the bed update ---- */
    OMP_PARALLEL_LOOP
    for (anuga_int k = 0; k < n; k++) {
        const double dz = dzs[k];
        if (dz != 0.0) {
            bed_cv[k] += dz;
            const anuga_int k3 = 3 * k;
            bed_ev[k3 + 0] += dz;
            bed_ev[k3 + 1] += dz;
            bed_ev[k3 + 2] += dz;
        }
    }
}

// ============================================================================
// Angle-of-repose relaxation  (spec 7)
// ============================================================================

// Diffuse bed material downslope wherever the centroid-to-centroid slope
// exceeds the critical angle, until it does not.
//
// FG21 §2.2.4, who are explicit that this is a NUMERICAL HEURISTIC and not
// physics: real bed slope failures are advective. It exists to stop the rest of
// the model breaking on over-steep slopes -- and it has a side effect worth
// remembering, that it limits the steepness of canyon walls and knickpoints and
// so suppresses knickpoint retreat that may be real.
//
// THE ONLY NON-CELL-LOCAL SEDIMENT KERNEL. A cell's update depends on its
// neighbours' elevation, which forces two things:
//
//   * Jacobi, not Gauss-Seidel. Each sweep reads elevation and writes
//     increments to a separate array, so no cell sees a neighbour that has
//     already moved this sweep. Reading live elevation would make the answer
//     depend on which thread got there first, and mode 1 and mode 2 would
//     diverge.
//   * a hard sweep cap, reported to the caller. Relaxation is iterative and a
//     pathological bed could otherwise spin.
//
// MASS IS CONSERVED, which is what makes this different from the older
// sanddune operator: material removed from an over-steep cell is DEPOSITED ON
// ITS NEIGHBOUR, never discarded. The mechanism is the same one bedload uses --
// the transfer is computed per EDGE from data both cells share, so both compute
// the identical volume and the pair balances exactly:
//
//     V = relax (|dz| - tan(theta) d) / (1/A_k + 1/A_nb)
//
// which is the volume that brings the pair exactly to the threshold slope when
// relax = 1: moving V lowers the donor by V/A_donor and raises the receiver by
// V/A_receiver, closing the excess by V(1/A_k + 1/A_nb).
//
// Interacts with [L-5]: material that cannot be scoured cannot slump either, so
// a transfer is capped by the DONOR's erodible thickness. The cap is a third of
// it per edge, because a cell has three edges and may be the donor on all of
// them; a full-thickness cap on each could lower it to three times its
// available depth in one sweep. The remainder is simply moved on later sweeps.
//
// Returns the number of sweeps used. Equal to max_sweeps means the cap was hit
// and the bed may still be over-steep -- the caller reports that rather than
// letting it pass silently.
/* Relative tolerance on the threshold slope; see the note at its use. */
#define REPOSE_TOL 1.0e-3

anuga_int core_apply_repose(struct domain *D) {
    const double tan_c = D->sediment_repose_tan;
    if (!(tan_c > 0.0)) {
        return 0;                      /* disabled, the default */
    }
    if (!D->sediment_bed_evolution) {
        return 0;                      /* a fixed bed cannot slump */
    }

    const anuga_int n = D->number_of_elements;
    const anuga_int max_sweeps = D->sediment_repose_max_sweeps;
    const double relax = D->sediment_repose_relax;

    double * restrict bed_cv = D->bed_centroid_values;
    double * restrict bed_ev = D->bed_edge_values;
    double * restrict dz = D->sediment_repose_dz;
    anuga_geom_t * restrict areas = D->areas;
    anuga_geom_t * restrict cc = D->centroid_coordinates;
    anuga_int * restrict neighbours = D->neighbours;
    double * restrict z_base = D->sediment_z_base;
    const anuga_int has_z_base = (D->sediment_has_z_base && z_base != NULL);

    if (dz == NULL || max_sweeps <= 0) {
        return 0;
    }

    anuga_int sweeps = 0;
    for (anuga_int it = 0; it < max_sweeps; it++) {
        anuga_int n_steep = 0;

        OMP_PARALLEL_LOOP_REDUCTION_PLUS(n_steep)
        for (anuga_int k = 0; k < n; k++) {
            double acc = 0.0;
            const double z_k = bed_cv[k];
            const double A_k = areas[k];
            const double xk = cc[2 * k];
            const double yk = cc[2 * k + 1];

            for (anuga_int i = 0; i < 3; i++) {
                const anuga_int nb = neighbours[3 * k + i];
                if (nb < 0 || nb == k) {
                    continue;          /* boundary, or a ghost self-reference */
                }
                const double dx = xk - cc[2 * nb];
                const double dy = yk - cc[2 * nb + 1];
                const double d = sqrt(dx * dx + dy * dy);
                if (!(d > 0.0)) {
                    continue;
                }

                const double z_nb = bed_cv[nb];
                const double diff = z_k - z_nb;
                const double adiff = fabs(diff);
                const double thresh = tan_c * d;

                /* Convergence is ASYMPTOTIC: each sweep removes a fraction of
                 * the excess, so a strict `> thresh` test is never satisfied
                 * and the loop runs to its cap every timestep, reporting a
                 * failure that has not happened. Measured on an over-steep
                 * cone: 36.84 -> 30.09 degrees against a 30 degree limit in
                 * 400 sweeps, still "not converged".
                 *
                 * So converged means within REPOSE_TOL of the threshold
                 * slope, which bounds the final angle: at 30 degrees a
                 * tolerance of 1e-3 in tan leaves at most 30.03 degrees.
                 * Deliberately not exposed -- it is the kernel's own
                 * convergence criterion, not a physical parameter, and the
                 * physical knob (the angle) is already there. */
                if (!(adiff > thresh * (1.0 + REPOSE_TOL))) {
                    continue;
                }
                n_steep++;

                const double A_nb = areas[nb];
                const double inv_sum = 1.0 / A_k + 1.0 / A_nb;

                /* The /3 is a STABILITY limit, not a fudge. relax = 1 with no
                 * divisor brings a single over-steep PAIR exactly to the
                 * threshold in one sweep -- but a cell has three edges and can
                 * be the donor on all of them at once, so its total change is
                 * up to three times what any one edge intended. That is an
                 * explicit diffusion step past its stability limit: it
                 * overshoots, creates fresh over-steep edges on the far side,
                 * and oscillates instead of converging. Observed with relax =
                 * 0.5 and no divisor: an over-steep cone stalled at 30.09
                 * degrees against a 30 degree limit and burned all 400 sweeps
                 * every timestep, so the cap was reported hit on a problem
                 * that was simply never going to converge.
                 *
                 * Dividing by the edge count bounds a cell's total movement by
                 * relax times its worst excess, which is stable for any
                 * relax <= 1, and keeps relax meaning what the docstring says
                 * it means. Both cells divide by the same 3, so the transfer
                 * stays symmetric and B1 conservation is untouched. */
                double V = relax * (adiff - thresh) / inv_sum / 3.0;

                /* [L-5]: the donor cannot give up what it may not lose. Both
                 * cells identify the same donor -- the higher one -- and
                 * compute the same cap, so the transfer stays symmetric. */
                if (has_z_base) {
                    const anuga_int donor = (diff > 0.0) ? k : nb;
                    const double avail = (bed_cv[donor] - z_base[donor])
                                       * areas[donor] / 3.0;
                    if (V > avail) {
                        V = (avail > 0.0) ? avail : 0.0;
                    }
                }

                /* The higher cell gives, the lower receives. */
                acc += (diff > 0.0) ? (-V / A_k) : (V / A_k);
            }
            dz[k] = acc;
        }

        sweeps = it + 1;
        if (n_steep == 0) {
            /* Nothing was over-steep, so nothing was written; stop before
             * applying a sweep of zeros. */
            sweeps = it;
            break;
        }

        OMP_PARALLEL_LOOP
        for (anuga_int k = 0; k < n; k++) {
            const double d = dz[k];
            if (d != 0.0) {
                bed_cv[k] += d;
                const anuga_int k3 = 3 * k;
                bed_ev[k3 + 0] += d;
                bed_ev[k3 + 1] += d;
                bed_ev[k3 + 2] += d;
            }
        }
    }

    return sweeps;
}


// ============================================================================
// Suspended sediment source terms  (Phase 3)
// ============================================================================

// Apply the sediment exchange term E_s - D_s of [G-3], and the resulting bed
// change [G-4], for every registered sediment class.
//
// THIS IS A FRACTIONAL STEP. It is called ONCE PER TIMESTEP with the full dt,
// after the hydrodynamic step, not inside the RK substep loop. So it updates
// the STATE directly (m and z) rather than contributing to a tendency:
//
//     m_s  <-  m_s + dt (E_s - D_s)                       [G-3] source part
//     z    <-  z   + dt (D - E)/(1 - lambda)              [G-4]
//
// The two use the SAME limited source, so the sediment volume leaving
// suspension equals (1 - lambda) dz exactly and the budget closes by
// construction, whatever the timestepping method.
//
// Stage is deliberately NOT adjusted here: holding w while z rises makes the
// depth h = w - z fall by exactly dz, which is the quiescent-water behaviour
// LM15 Example 2 requires (their free surface stays flat while the bed
// aggrades). The pore space in the new bed is filled from the water column.
//
// Phase 3 is the FIXED-BED stage of spec 2.4: the bed does not evolve, there is
// no bed -> flow feedback and no sediment -> momentum feedback. Deposited mass
// simply leaves suspension.
//
// Deposition is [D-1]:      D_s = d*(Z) . c_s . v_s
// with d* the near-bed to depth-averaged concentration ratio (1.0 = well
// mixed). Erosion E_s is Phase 3b and is zero here.
//
// TWO LIMITERS, both from spec 4.5, and both applied to the SOURCE TERM rather
// than by clamping the state -- clamping m would break the exact conservation
// the transport scheme provides:
//
//   [L-1] positivity, mandatory:  F_s^net >= -m_s / dt
//         deposition can never remove more sediment than is present.
//   [L-2] concentration ceiling:  c_s <= c_max
//         applied as a cap on how much a cell may GAIN this step.
//
// Called after the flux kernel has filled tracer_explicit_update and before the
// time integration consumes it, so the source lands in the same dm/dt the
// integrator already applies.
void core_apply_sediment_source(struct domain *D, double timestep) {
    const anuga_int n_classes = D->n_sediment_classes;
    if (n_classes <= 0 || timestep <= 0.0) {
        return;
    }

    const anuga_int n = D->number_of_elements;
    const double c_max = D->sediment_c_max;
    const double minimum_allowed_height = D->minimum_allowed_height;

    const double gamma0 = D->sediment_gamma0;
    const anuga_int erosion_mode = D->sediment_erosion_mode;
    const double tau_crit = D->sediment_tau_crit;
    const double K_e = D->sediment_K_e;
    const double rho_w = D->sediment_rho_w;
    const double K_p = D->sediment_K_partheniades;
    const double dl_A = D->sediment_dl_A;
    const double dl_alpha = D->sediment_dl_alpha;
    const double dl_beta = D->sediment_dl_beta;
    const double dl_b = D->sediment_dl_threshold;
    const double dl_ks = D->sediment_dl_ks;
    const double dl_ks_factor = D->sediment_dl_ks_factor;
    const anuga_int dep_mode = D->sediment_deposition_mode;
    const double tau_d = D->sediment_tau_d;
    const anuga_int shear_closure = D->sediment_shear_closure;
    const double max_slope = D->sediment_max_slope;
    const anuga_int slope_frozen = D->sediment_slope_frozen;
    const double h_eps = D->epsilon;
    const double grav = D->g;

    double * restrict stage_cv = D->stage_centroid_values;
    double * restrict bed_cv = D->bed_centroid_values;
    double * restrict xmom_cv = D->xmom_centroid_values;
    double * restrict ymom_cv = D->ymom_centroid_values;
    double * restrict friction_cv = D->friction_centroid_values;
    double * restrict v_s = D->sediment_settling_velocity;
    double * restrict d_star = D->sediment_d_star;
    double * restrict diam = D->sediment_diameter;
    double * restrict sedR = D->sediment_R;
    double * restrict tau_c_star = D->sediment_tau_c_star;
    double * restrict a_ref = D->sediment_reference_height;
    anuga_int * restrict neighbours_r = D->neighbours;
    anuga_geom_t * restrict cc_r = D->centroid_coordinates;
    double * restrict slope_w = D->sediment_slope_work;
    const anuga_int d_star_mode = D->sediment_d_star_mode;
    const double a_h_floor = D->sediment_a_h_floor;
    const anuga_int rouse_beta_mode = D->sediment_rouse_beta_mode;
    const double rouse_scale = D->sediment_rouse_scale;
    const anuga_int adapt_mode = D->sediment_adaptation_mode;
    const double adapt_alpha = D->sediment_adaptation_alpha;
    const anuga_int nb_base = D->sediment_nearbed_base;
    const double layer_frac = D->sediment_layer_fraction;
    const double exch_fac = D->sediment_exchange_factor;
    const anuga_int vprof = D->sediment_velocity_profile;
    double * restrict t_sf = D->tracer_speed_factor;
    double * restrict t_bv = D->tracer_boundary_values;
    const anuga_int t_bl = D->boundary_length;
    const double c_pack = D->sediment_c_pack;
    const anuga_int fric_mode = D->sediment_friction_mode;
    const double n_ll = D->sediment_manning_ll;
    const anuga_int wbed = D->sediment_wilson_bed;
    const double wD = D->sediment_wilson_D;
    // Vegetation seen by the sediment shear (core_vegetation_sediment_fc).
    const anuga_int veg_shear = (D->vegetation_mode > 0) ? D->sediment_vegetation_shear : 0;
    double * restrict veg_m = D->veg_density_centroid_values;
    double * restrict shear_fac = D->sediment_shear_factor;   /* [T-16], or NULL */
    double * restrict veg_d = D->veg_diameter_centroid_values;
    double * restrict veg_h = D->veg_height_centroid_values;
    const double veg_Cd = D->vegetation_Cd;
    const double veg_Cb = D->vegetation_bed_chezy;

    // Hoisted for the same reason as in the update/backup/saxpy kernels: on a
    // GPU build the loop below is an 'omp target' region and D is NOT mapped to
    // the device, so a D->member load inside it reads a host address and the
    // work silently does not happen. See HANDOVER.md 2.4.
    double * restrict t_cons = D->tracer_conserved_values;
    double * restrict ext_src = D->tracer_external_source;
    double * restrict bed_cv_w = D->bed_centroid_values;
    double * restrict bed_ev_w = D->bed_edge_values;
    const double one_minus_lambda = 1.0 - D->sediment_porosity;
    const double morfac = D->sediment_morphological_factor;
    const anuga_int bed_evolves = D->sediment_bed_evolution;
    /* [L-5]. Hoisted for the same device reason as the tracer pointers.
     *
     * The flag alone does not license the dereference: it and the pointer are
     * set by separate lines of the Cython binding, and when one of them was
     * missed the flag read as uninitialised garbage, tested true, and this
     * kernel dereferenced a NULL base. Require both. */
    double * restrict z_base = D->sediment_z_base;
    const anuga_int has_z_base = (D->sediment_has_z_base && z_base != NULL);
    double * restrict src_lim = D->sediment_source_limited;
    /* Bed composition. Hoisted like the rest; both arrays are required, as
     * for z_base above, so a half-bound composition cannot be dereferenced. */
    double *bed_act = D->sediment_bed_active;
    double *bed_sub = D->sediment_bed_substrate;
    const anuga_int has_comp = (D->sediment_bed_composition && bed_act != NULL
                                && bed_sub != NULL && bed_evolves
                                && one_minus_lambda > 0.0);
    const double active_solid = D->sediment_active_layer * one_minus_lambda;
    /* The GPU loop below maps these with explicit lengths, which a NULL
     * pointer cannot take. Off, they are never touched, so stand them on the
     * scratch array that is always mapped. (Not restrict, for that reason.) */
    if (!has_comp) {
        bed_act = src_lim;
        bed_sub = src_lim;
    }

    /* The slope for [T-7]/[T-7e], in a pass of its own. It reads the
     * neighbours' bed (or stage) centroids, and the loop below WRITES the
     * bed when it evolves, so taken inside that loop it would see whatever a
     * preceding iteration or another thread had already done to a
     * neighbour: loop-order dependent, and mode 1 and mode 2 disagreed. */
    if ((shear_closure == 1 && !slope_frozen) || shear_closure == 2) {
        const double * restrict f = (shear_closure == 2) ? stage_cv : bed_cv;
        OMP_PARALLEL_LOOP
        for (anuga_int k = 0; k < n; k++) {
            slope_w[k] = core_centroid_slope(f, neighbours_r, cc_r, k);
        }
    }
    /* With a frozen depth-slope, slope_w holds the slope of the bed at
     * setup (Domain.set_shear_closure / bed_slope_magnitude). */

    /* On the GPU the per-class arrays are mapped with EXPLICIT lengths. Left
     * implicit, nvc infers a length only while the loop is simple enough to
     * see n_classes bound every access; past that it maps them as zero-length
     * sections, the device read uninitialised memory (compute-sanitizer
     * initcheck at the sedR/diam load), and mode 2 went wrong and
     * nondeterministic with nothing in the physics changed. Present data are
     * reused, so the clauses cost nothing. */
#ifdef CPU_ONLY_MODE
    OMP_PARALLEL_LOOP
#else
    #pragma omp target teams loop \
        map(to: v_s[0:n_classes], d_star[0:n_classes], diam[0:n_classes], \
                sedR[0:n_classes], tau_c_star[0:n_classes], a_ref[0:n_classes]) \
        map(tofrom: src_lim[0:n_classes*n], bed_act[0:n_classes*n], \
                    bed_sub[0:n_classes*n])
#endif
    for (anuga_int k = 0; k < n; k++) {
        const double h = fmax(stage_cv[k] - bed_cv[k], 0.0);

        // Dry cells carry no sediment: spec 9.4 sets c_s = 0 below the wet/dry
        // threshold. Leave the advective tendency untouched and add nothing.
        if (h <= minimum_allowed_height) {
            continue;
        }

        const double inv_h = 1.0 / h;
        const double m_max = c_max * h;
        double dz_cell = 0.0;
        double total_E = 0.0;   /* [L-5]: erosive demand on the bed */
        double total_D = 0.0;   /* [L-5]: what deposition returns to it */

        // Velocity by the ANUGA depth-limiting form [T-5], the same
        // regularisation RDy26 A22-A23 adopt. Never a bare (uh)/h.
        const double denom = h * h + h_eps * h_eps;
        const double u = (denom > 0.0) ? (xmom_cv[k] * h / denom) : 0.0;
        const double v = (denom > 0.0) ? (ymom_cv[k] * h / denom) : 0.0;
        const double vel2 = u * u + v * v;

        // Friction closure, spec 3.3. In EVERY mode f_c varies per cell per
        // timestep -- through h in [T-6] for the Manning-based modes, and
        // through the relative submergence h/D for wilson. Recomputed here
        // rather than cached, which spec 3.3 calls the coupling most easily
        // missed.
        double f_c;
        if (fric_mode == 2) {
            // wilson [T-8]..[T-10]. W04's equations give X = (8/f_W04)^1/2
            // with f_W04 the Darcy-Weisbach f; ours is f/8, so f_c = 1/X^2.
            // See the note in sw_domain.h: taking W04's f_c literally as ours
            // would make tau_b 8x too large.
            //
            // The relations assume the clasts are submerged. Below h = D the
            // logarithms go to zero or negative and the power law leaves its
            // calibration, so relative submergence is floored at 1.
            double rel = h / wD;
            if (!(rel > 1.0)) rel = 1.0;
            double X;
            if (wbed == 0) {
                X = 8.46 * pow(rel, 0.1005);           /* sand,    [T-8]  */
            } else if (wbed == 1) {
                X = 5.75 * log10(rel) + 3.514;         /* gravel,  [T-9]  */
            } else {
                X = 5.62 * log10(rel) + 4.0;           /* boulder, [T-10] */
            }
            f_c = 1.0 / (X * X);
        } else {
            // [T-6] f_c = g n^2 h^(-1/3) == RDy26 A5, with n either the
            // per-cell user field (constant mode) or the uniform Manning-
            // Strickler value of [T-14] (larsen_lamb).
            const double nman = (fric_mode == 1) ? n_ll : friction_cv[k];
            f_c = grav * nman * nman / cbrt(h);
        }
        if (veg_shear) {
            const double f_v = core_vegetation_sediment_fc(veg_shear, h, veg_m[k], veg_d[k],
                                                           veg_h[k], veg_Cd, veg_Cb, grav);
            if (f_v >= 0.0) f_c = f_v;
        }
        if (shear_fac != NULL) f_c *= shear_fac[k];           /* [T-16] */

        /* tau_b/rho under the selected closure: [T-1], [T-7] or [T-7e]. */
        double S = (shear_closure == 1 || shear_closure == 2)
                 ? slope_w[k] : 0.0;
        if (max_slope > 0.0 && S > max_slope) S = max_slope;
        const double tbr = core_tau_b_over_rho(shear_closure, f_c, vel2, grav,
                                               h, S);

        /* Bed composition: the active layer's content, for the class shares
         * F_s = bed_act[s] / act_tot at the START of the step. */
        double act_tot = 0.0;
        if (has_comp) {
            for (anuga_int s = 0; s < n_classes; s++) {
                act_tot += bed_act[s * n + k];
            }
        }

        for (anuga_int s = 0; s < n_classes; s++) {
            const anuga_int idx = s * n + k;

            const double m = t_cons[idx];
            const double c = m * inv_h;

            // Deposition is computed from a NON-NEGATIVE concentration.
            //
            // m can go slightly negative through the ADVECTIVE tendency (the
            // transport scheme guarantees positivity only under CFL, and the
            // source is added to a tendency it does not control). Fed through
            // unguarded, deposition = d* c v_s flips sign and starts ADDING
            // sediment -- and [L-1] below compounds it, because -m/dt is a
            // POSITIVE lower bound when m < 0, which forces the source
            // positive. Together they created 957% of the initial mass in a
            // deposition-only run. Both paths are guarded here.
            const double m_pos = (m > 0.0) ? m : 0.0;
            const double c_pos = m_pos * inv_h;

            // [T-2] u* = |v| sqrt(f_c);  [S-2] Z = v_s / (kappa u*)
            const double ustar = sqrt(f_c * vel2);
            double Z = (ustar > 0.0)
                     ? v_s[s] / (0.41 * ustar)
                     : ANUGA_ROUSE_Z_HI;   /* no shear: fully settled */
            // [S-2b] Rouse-number correction. Van Rijn (1984b) divides Z by
            // beta = 1 + 2 (w_s/u*)^2 (at most 2): grains are mixed more
            // strongly than momentum, so the measured profiles are flatter
            // than the Rouse one. His trench profiles (1986b, Fig. 17) need
            // Z about 0.5 where the plain value is 0.79 and beta gives
            // 0.65; the scale is for that remainder and for testing.
            if (ustar > 0.0) {
                if (rouse_beta_mode == 1) {
                    const double r = v_s[s] / ustar;
                    double beta = 1.0 + 2.0 * r * r;
                    if (beta > 2.0) beta = 2.0;
                    Z /= beta;
                }
                Z *= rouse_scale;
            }
            // a/h with the van Rijn-style floor a >= floor*h. The floor is
            // standard practice and stays on by default, but it is exposed:
            // it is the single largest divergence from anugaSed, which uses
            // no floor and so an ~10x smaller a at h = 1 m, giving roughly
            // 8x more deposition (spec 12, D4b). Set it to 0 to reach that
            // regime; the fit now covers a/h down to 1e-3.
            double a_h = a_ref[s] * inv_h;
            if (a_h < a_h_floor) a_h = a_h_floor;

            // [D-1] deposition. d* is either the constant (P14's d* = 1
            // limiting case) or the Rouse ratio evaluated per cell.
            const double ds = (d_star_mode == 0) ? d_star[s]
                                                 : core_rouse_d_star(Z, a_h);

            // [D-3] ADAPTATION LAG of the near-bed concentration.
            //
            // With E = v_s E* and D = d* v_s c the exchange is
            // d* v_s (c_eq - c), c_eq = E*/d*: the near-bed concentration is
            // assumed to be the equilibrium one for the local flow at every
            // instant. The vertical profile actually takes a time of order
            // h/(alpha w_s) to adjust (Galappatti & Vreugdenhil 1985): grains
            // entrained at the bed must diffuse up the column before the
            // load is carried, and grains high in the column must settle
            // through it before deposition is felt. Their depth-integrated
            // model relaxes the load toward the same equilibrium,
            //
            //     E - D = alpha v_s (c_eq - c),
            //
            // so both terms are scaled by alpha/d*: every equilibrium is
            // untouched, only the rate changes. alpha is Armanini & Di
            // Silvio's (1988) closed form of Galappatti's coefficient,
            //
            //     1/alpha = a/h + (1 - a/h) exp[-1.5 (a/h)^(-1/6) w_s/u*],
            //
            // which is 1 in the well-mixed limit (w_s/u* -> 0: the exchange
            // is the well-mixed one, no faster) and h/a when fully
            // stratified. Van Rijn's pick-up flume reaches equilibrium in
            // ~15 depths without it against >40 measured, and his trench
            // fills ~25% too fast (issue #389). Mode 2 takes alpha as given.
            double f_adapt = 1.0;
            // [D-4] CARRIED NEAR-BED RATIO (mode 3). Instead of scaling the
            // rate, keep the stratification of the suspension as a state:
            // tracer (nb_base + s) carries h r_b, with r_b = c_b / c the
            // ratio of near-bed to depth-averaged concentration the column
            // actually has, advected with the flow by the flux kernel. Its
            // equilibrium is d*. The lag is one-sided. When the flow slows,
            // d* rises above r_b: the profile is collapsing, the grains high
            // in the column settle from about the profile centroid, and r_b
            // approaches d* over T = (z_c - a)/w_s. Deposition uses the
            // carried ratio, so a parcel entering a trench starts depositing
            // at the stratification it brought and works up to the trench's
            // over that time. When d* is at or below r_b (the flow has
            // quickened, or nothing changed) the near-bed concentration is
            // at the bed and responds at once: r_b = d*. A ratio of zero is
            // "not set yet" (initial state, inflow) and takes d*. The slow
            // filling of the UPPER column that delays the depth-integrated
            // load when a bed loads clear water is not a near-bed lag and is
            // not represented here; that needs a layered suspension. Exact
            // exponential relaxation over the step; T -> 0 recovers [D-1].
            double c_b_used = -1.0;
            if (adapt_mode == 3 && nb_base >= 0) {
                const anuga_int nidx = (nb_base + s) * n + k;
                double r_b = t_cons[nidx] * inv_h;
                double r_n = ds;
                if (r_b > 0.0 && ds > r_b && v_s[s] > 0.0) {
                    const double T = core_rouse_z_c(Z, a_h) * h / v_s[s];
                    if (T > 0.0) r_n = ds + (r_b - ds) * exp(-timestep / T);
                }
                t_cons[nidx] = r_n * h;
                c_b_used = r_n * c_pos;
                /* The ratio tracer's inflow value is the local equilibrium
                 * d*: written here on the cell's boundary edges every step
                 * (the array is device-resident in mode 2, as this kernel
                 * is), so inflowing water carries no spurious lag in. */
                for (anuga_int i = 0; i < 3; i++) {
                    const anuga_int nbk = neighbours_r[3 * k + i];
                    if (nbk < 0 && t_bv != NULL) {
                        t_bv[(nb_base + s) * t_bl + (-nbk - 1)] = ds;
                    }
                }
            }
            // [D-5] TWO-LAYER SUSPENSION (mode 4). The column is split into
            // a near-bed layer of thickness h1 = a (the reference height,
            // floor included) and the rest, h2 = h - h1. The fraction's own
            // tracer still carries the TOTAL mass m = h c, with its
            // advection, limiters and bed exchange unchanged; tracer
            // (nb_base + s) carries the upper layer's mass m2 = h2 c2, and
            // the lower layer is what remains, m1 = m - m2. Between the
            // layers: settling v_s c2 down, exchange K (c1 - c2) up. K is
            // set so that the two-layer equilibrium c2/c1 = rho reproduces
            // the fitted Rouse ratio, c1/c = d*: rho = (1/d* - f1)/(1 - f1),
            // K = v_s rho/(1 - rho), so every equilibrium is unchanged and
            // only the transient differs. Deposition is v_s c1, entrainment
            // enters layer 1 (through the total). The partition relaxes
            // toward its equilibrium at the rate lam = K/h1 + (K + v_s)/h2,
            // integrated exactly. So: a parcel entering slower water keeps
            // its upper-layer load and deposits at the near-bed value it
            // brought while the upper layer settles out over ~h2/v_s; and a
            // bed loading clear water fills the lower layer first, so the
            // near-bed concentration and hence deposition lead the depth-
            // averaged load, which grows only as sediment is exchanged up.
            // Both layers are advected with the depth-averaged velocity
            // (no velocity profile yet). m2 < 0 means "not set" and takes
            // the equilibrium partition; the kernel writes the boundary
            // value of the upper-layer tracer as the equilibrium partition
            // of the fraction's own boundary concentration every step.
            if (adapt_mode == 4 && nb_base >= 0 && ds > 0.0) {
                /* ds > 0 is required: the partition is built from 1/d*, and
                 * d* = 0 (deposition switched off through the per-fraction
                 * d_star) has no stratification to carry. */
                const anuga_int nidx = (nb_base + s) * n + k;
                double f1 = (layer_frac > 0.0) ? layer_frac : a_h;
                if (f1 > 0.5) f1 = 0.5;
                const double h1 = f1 * h;
                const double h2 = h - h1;
                double rho = (1.0 / ds - f1) / (1.0 - f1);
                if (rho < 1.0e-6) rho = 1.0e-6;
                if (rho > 1.0) rho = 1.0;
                const double K = v_s[s] * rho / (1.0 - rho + 1.0e-12);
                const double m2_eq = h2 * rho * m_pos / (h1 + h2 * rho);
                double m2 = t_cons[nidx];
                double m2n;
                if (m2 < 0.0) {
                    m2n = m2_eq;
                } else {
                    if (m2 > m_pos) m2 = m_pos;
                    const double lam = exch_fac * (K / h1 + (K + v_s[s]) / h2);
                    m2n = m2_eq + (m2 - m2_eq) * exp(-lam * timestep);
                }
                t_cons[nidx] = m2n;
                double c1 = (m_pos - m2n) / h1;
                if (c1 < 0.0) c1 = 0.0;
                c_b_used = c1;
                // [D-5v] VELOCITY PROFILE. With a log law u(z) = (u*/kappa)
                // ln(z/z0) the depth mean is (u*/kappa) L, L = ln(h/z0) - 1
                // = kappa/sqrt(f_c), and the means over the lower layer
                // (0, f1 h) and the rest are
                //     u1/u = 1 + ln(f1)/L,
                //     u2/u = 1 - f1 ln(f1)/((1 - f1) L),
                // so the near-bed layer, which holds most of the sediment,
                // lags the flow (0.55 u for f1 = 0.1 at k_s = 0.025 m,
                // h = 0.39 m) and the upper layer leads it slightly. The
                // fraction's own tracer (the total) is advected at the
                // mass-weighted mean of the two and the upper-layer tracer
                // at u2; the factors are per cell and go through the flux
                // kernel next step. The total sediment flux is then the
                // integral of u c over the profile, less than u h c.
                if (vprof && t_sf != NULL) {
                    double b1 = 1.0, b2 = 1.0;
                    /* L = ln(h/z0) - 1 = kappa/sqrt(f_c). The log law is only
                     * a profile at all once L exceeds |ln f1|, which is where
                     * the near-bed mean is still forward-going; below that
                     * (a thin or very rough cell, f_c large) the split is
                     * meaningless and BOTH factors are left at 1. Without
                     * this the upper factor runs away as h -> 0 at a wetting
                     * front and the tracer flux collapses the time step. */
                    const double Lp = (f_c > 0.0) ? 0.41 / sqrt(f_c) : 0.0;
                    const double lf1 = log(f1);
                    if (Lp > -lf1) {
                        b1 = 1.0 + lf1 / Lp;
                        if (b1 < 0.05) b1 = 0.05;
                        b2 = 1.0 - f1 * lf1 / ((1.0 - f1) * Lp);
                        if (b2 > 2.0) b2 = 2.0;
                    }
                    const double m1n = m_pos - m2n;
                    t_sf[nidx] = b2;
                    t_sf[idx] = (m_pos > 0.0) ? (m1n * b1 + m2n * b2) / m_pos : 1.0;
                }
                for (anuga_int i = 0; i < 3; i++) {
                    const anuga_int nbk = neighbours_r[3 * k + i];
                    if (nbk < 0 && t_bv != NULL) {
                        const double cbd = t_bv[s * t_bl + (-nbk - 1)];
                        /* the tracer's boundary value is m2/h at the
                         * equilibrium partition of the boundary c */
                        t_bv[(nb_base + s) * t_bl + (-nbk - 1)] =
                            (cbd > 0.0) ? h2 * rho * cbd / (h1 + h2 * rho) : 0.0;
                    }
                }
            }
            if (adapt_mode == 1) {
                double alpha;
                if (ustar > 0.0) {
                    const double r = v_s[s] / ustar;
                    const double inv = a_h + (1.0 - a_h)
                                     * exp(-1.5 * pow(a_h, -1.0 / 6.0) * r);
                    alpha = 1.0 / inv;
                } else {
                    alpha = 1.0 / a_h;          /* no shear: stratified limit */
                }
                f_adapt = alpha / ds;
            } else if (adapt_mode == 2) {
                f_adapt = adapt_alpha / ds;
            }
            // [L-4] NEAR-BED CONCENTRATION IS BOUNDED BY PACKING.
            //
            // [D-1] is D = c_b v_s with c_b = d* c, and nothing in the spec
            // bounds c_b. It needs bounding. d* comes from the EQUILIBRIUM
            // Rouse profile, which is not valid as shear vanishes: at rest
            // u* -> 0, so Z -> infinity and d* -> its clamp (~250 at
            // a/h = 0.01), making the deposition rate enormous. A lake at rest
            // then deposits its entire suspended load in under a second,
            // instead of over the physical h/v_s.
            //
            // c_b is a concentration and cannot exceed maximum packing, the
            // same 0.65 that bounds E* in [E-1]. Capping c_b there keeps the
            // still-water limit sane while leaving the well-mixed and
            // moderate-Z regimes untouched, where d* c is far below packing.
            //
            // Added in PHYSICS_SPEC Draft 5 as [L-4]; it has no counterpart in
            // P14, FG21, RDy26, DL09 or aSM16, being required by combining an
            // equilibrium profile with a transient solver.
            double deposition;
            if (dep_mode == 1) {
                /* [D-2] RDy26's threshold form. tau_d = 0 disables deposition
                 * entirely -- the hook their passive benchmarks rely on. */
                const double tau_b_d = rho_w * tbr;
                deposition = (tau_d > 0.0 && tau_b_d < tau_d)
                           ? v_s[s] * c_pos * (1.0 - tau_b_d / tau_d)
                           : 0.0;
            } else {
                double c_bed = (c_b_used >= 0.0) ? c_b_used : ds * c_pos;
                if (c_bed > c_pack) c_bed = c_pack;
                deposition = c_bed * v_s[s];
            }

            // [E-1]/[E-2] entrainment, non-cohesive (Shields) route.
            //
            //   tau* = f_c |v|^2 / (R g d)   [T-3] -- rho cancels
            //   S    = tau*/tau_c* - 1
            //   E*   = 0.65 gamma0 S / (1 + gamma0 S)   saturating
            //
            // Below threshold (S <= 0) there is no entrainment at all; this is
            // a genuine threshold, not a smooth roll-off.
            double erosion = 0.0;
            if (erosion_mode == 2) {
                /* [E-4] Partheniades. K_p is a MASS flux, so divide by the
                 * class density to get the volume flux the rest of the source
                 * term works in. rho_s = (R + 1) rho_w. */
                const double tau_b = rho_w * tbr;
                if (tau_crit > 0.0 && tau_b > tau_crit) {
                    const double rho_s = (sedR[s] + 1.0) * rho_w;
                    if (rho_s > 0.0) {
                        erosion = (K_p * (tau_b - tau_crit) / tau_crit) / rho_s;
                    }
                }
            } else if (erosion_mode == 1) {
                /* [E-3] cohesive, Hanson & Simon. DIMENSIONAL excess shear:
                 * tau_b = rho f_c |v|^2 [T-1], and E = K_e (tau_b - tau_c),
                 * zero below threshold. Note this is per class only through
                 * the loop -- tau_c and K_e are bed properties, not grain
                 * properties, which is precisely the cohesive premise. */
                const double tau_b = rho_w * tbr;
                const double excess = tau_b - tau_crit;
                if (excess > 0.0) {
                    erosion = K_e * excess;
                }
            } else if (erosion_mode == 3) {
                /* [E-6] de Leeuw et al. (2020), Eq 26a, non-cohesive bed
                 * material load (and, with the Nghiem et al. 2022 constants,
                 * flocculated mud as the Delta-X Wax Lake model uses it):
                 *
                 *   E* = A X^beta / (1 + 3 A X^beta)
                 *   X  = (u*_skin / v_s)^alpha Fr - threshold,  Fr = U / sqrt(g h)
                 *
                 * u*_skin from Manning-Strickler on a roughness k_s (Eq 7):
                 *   U / u*_skin = 8.1 (H_sk / k_s)^(1/6),  u*_skin = sqrt(g H_sk S)
                 * so H_sk = (U k_s^(1/6) / (8.1 sqrt(g S)))^(3/2), capped at h
                 * (skin friction cannot exceed the total). The friction slope
                 * is taken from the active shear closure, S = tau_b/(rho g h),
                 * which under quadratic drag is f_c U^2 / (g h), the Delta-X
                 * model's form. k_s is floored at the viscous sublayer,
                 * nu / (8 u*), as there. E* is the near-bed concentration at
                 * 0.1 h, capped at 1/3 by the denominator. */
                const double U = sqrt(vel2);
                if (U > 0.0 && tbr > 0.0 && v_s[s] > 0.0 && h > 0.0) {
                    const double S_f = tbr / (grav * h);
                    const double ustar = sqrt(tbr);
                    double ks = (dl_ks > 0.0) ? dl_ks : dl_ks_factor * diam[s];
                    const double ks_visc = 1.0e-6 / (8.0 * ustar);
                    if (ks < ks_visc) ks = ks_visc;
                    double H_sk = pow(U * pow(ks, 1.0 / 6.0) / (8.1 * sqrt(grav * S_f)), 1.5);
                    if (H_sk > h) H_sk = h;
                    const double ustar_sk = sqrt(grav * H_sk * S_f);
                    const double Fr = U / sqrt(grav * h);
                    const double X = pow(ustar_sk / v_s[s], dl_alpha) * Fr - dl_b;
                    if (X > 0.0) {
                        const double aX = dl_A * pow(X, dl_beta);
                        erosion = v_s[s] * (aX / (1.0 + 3.0 * aX));
                    }
                }
            } else {
                /* [E-1]/[E-2] non-cohesive, Shields route. */
                const double Rgd = sedR[s] * grav * diam[s];
                if (Rgd > 0.0 && tau_c_star[s] > 0.0) {
                    const double tau_star = tbr / Rgd;
                    const double S = tau_star / tau_c_star[s] - 1.0;
                    if (S > 0.0) {
                        const double gS = gamma0 * S;
                        erosion = v_s[s] * (0.65 * gS / (1.0 + gS));
                    }
                }
            }

            // Bed composition: a class is entrained in proportion to its share
            // of the active layer (Hirano), so a bed with none of it cannot
            // supply it, and an emptied active layer supplies nothing.
            if (has_comp) {
                erosion = (act_tot > 0.0) ? erosion * bed_act[idx] / act_tot : 0.0;
            }

            // Net exchange of [G-3]. Deposition removes, erosion adds; [D-3]
            // scales the rate of both toward the same equilibrium.
            double source = f_adapt * (erosion - deposition);

            // [L-1] positivity. The most this term may remove over the step is
            // exactly the sediment PRESENT, so the state can reach zero but
            // never go below it. Applied to the SOURCE, not to m.
            //
            // m_pos, not m: with m < 0 the bound -m/dt is POSITIVE and would
            // force the source to inject sediment. The limiter must only ever
            // restrain removal, never mandate addition.
            const double min_source = -m_pos / timestep;
            if (source < min_source) {
                source = min_source;
            }

            // [L-2] ceiling. Only ever restrains a GAIN, so it cannot fight
            // [L-1] above: the two act on opposite signs of the source.
            if (source > 0.0) {
                const double max_source = (m_max - m) / timestep;
                if (source > max_source) {
                    source = (max_source > 0.0) ? max_source : 0.0;
                }
            }

            // Held, not applied. [L-5] below limits the classes TOGETHER
            // against the cell's erodible thickness, so no class may be
            // applied until every class's demand on the bed is known.
            //
            // The external supply is deliberately NOT folded in here: it is
            // not a bed exchange, so it must not be scaled by a bed-material
            // limiter, and it is added in the apply loop instead.
            //
            // Bed composition: a class may not take more out of the bed this
            // step than the cell holds of it (active layer and substrate).
            // Per class, so a class that is absent cannot borrow another's
            // share of the erodible thickness the way [L-5] alone allows.
            if (has_comp && source > 0.0) {
                const double held = bed_act[idx] + bed_sub[idx];
                const double cap = (held > 0.0) ? held / (timestep * morfac) : 0.0;
                if (source > cap) source = cap;
            }
            src_lim[idx] = source;
            if (source > 0.0) total_E += source;
            else              total_D += source;
        }

        // ---- [L-5] non-erodible base -----------------------------------
        //
        // The bed may be lowered to sediment_z_base and no further. Erosion
        // is a bed-material budget, so the limit belongs on the SOURCE, like
        // [L-1] and [L-2], and not on z: clamping z after the fact would
        // leave sediment in the water column that no longer came from
        // anywhere, which is exactly how [L-1]'s sign bug created 957% of
        // the initial mass.
        //
        // Only the EROSIVE part is scaled. Deposition is not restrained by a
        // shortage of bed material -- it is what supplies it -- and scaling
        // it down would suppress the very process that reopens the cell.
        //
        // The scale is shared and proportional, so the answer does not depend
        // on the order the classes were registered. There is no bed
        // stratigraphy in this model: the bed is not tracked per class, so
        // no class has a better claim on the last millimetre than another,
        // and proportional is the only choice that does not invent one.
        double scale = 1.0;
        if (has_z_base && bed_evolves && one_minus_lambda > 0.0
                && total_E > 0.0) {
            const double avail = bed_cv[k] - z_base[k];
            const double thickness = (avail > 0.0) ? avail : 0.0;
            // The largest net removal from the bed this step, as a source.
            // per step the bed moves morfac times the physical change
            const double S_max = thickness * one_minus_lambda / (timestep * morfac);
            if (total_E + total_D > S_max) {
                scale = (S_max - total_D) / total_E;
                if (scale < 0.0) scale = 0.0;
                if (scale > 1.0) scale = 1.0;
            }
        }

        // ---- apply ------------------------------------------------------
        for (anuga_int s = 0; s < n_classes; s++) {
            const anuga_int idx = s * n + k;
            double source = src_lim[idx];
            if (source > 0.0) {
                source *= scale;
            }

            // [G-4]. source = E - D, so dz = -source dt/(1-lambda). Taken
            // from the bed exchange ALONE, before the external supply is
            // added: sediment introduced from outside the model does not
            // come out of the bed, so it must not move it.
            if (bed_evolves && one_minus_lambda > 0.0) {
                dz_cell += -(morfac * timestep * source) / one_minus_lambda;
            }
            // The same exchange, as solid volume of this class, through the
            // active layer. It may go negative here (the cap above allows
            // what the substrate holds too); the rebalance below settles it.
            if (has_comp) {
                bed_act[idx] -= morfac * timestep * source;
            }

            // [G-3] S_ms: external supply, added AFTER the limiters. They
            // bound bed exchange by what bed and water column can supply;
            // an external source is neither, and clipping it would also make
            // a manufactured solution impossible to impose exactly.
            if (ext_src != NULL) {
                source += ext_src[idx];
            }

            // Fractional step: update the state directly with the full dt.
            t_cons[idx] += timestep * source;
        }

        // Bed composition: restore the active layer to its thickness. Every
        // transfer below moves an amount of a class from one layer to the
        // other, so each class is conserved exactly.
        if (has_comp) {
            // A class eroded past its active-layer content took the rest
            // from the substrate (the cap guarantees the substrate had it).
            double act = 0.0, sub = 0.0;
            for (anuga_int s = 0; s < n_classes; s++) {
                const anuga_int idx = s * n + k;
                if (bed_act[idx] < 0.0) {
                    bed_sub[idx] += bed_act[idx];
                    bed_act[idx] = 0.0;
                    if (bed_sub[idx] < 0.0) bed_sub[idx] = 0.0;   /* round-off */
                }
                act += bed_act[idx];
                sub += bed_sub[idx];
            }
            if (act > active_solid) {
                // Aggradation: the overflow goes down at the ACTIVE composition.
                const double r = (act - active_solid) / act;
                for (anuga_int s = 0; s < n_classes; s++) {
                    const anuga_int idx = s * n + k;
                    const double m = r * bed_act[idx];
                    bed_act[idx] -= m;
                    bed_sub[idx] += m;
                }
            } else if (act < active_solid && sub > 0.0) {
                // Degradation: refill from the substrate at ITS composition.
                const double need = active_solid - act;
                const double r = (need < sub) ? need / sub : 1.0;
                for (anuga_int s = 0; s < n_classes; s++) {
                    const anuga_int idx = s * n + k;
                    const double m = r * bed_sub[idx];
                    bed_sub[idx] -= m;
                    bed_act[idx] += m;
                }
            }
        }
        // Raise the bed by dz. The DE algorithms use DISCONTINUOUS elevation,
        // so edge values are not re-derived from the centroid and must be
        // shifted too; shifting all three by the same dz preserves the
        // within-cell bed slope, which is what keeps a flat bed flat. Vertex
        // values need no action: extrapolation recomputes them from the edges
        // (bed_vv = bed_ev1 + bed_ev2 - bed_ev0).
        //
        // Stage is left alone, so h = w - z falls by exactly dz -- the
        // quiescent-water behaviour of LM15 Example 2.
        if (dz_cell != 0.0) {
            bed_cv_w[k] += dz_cell;
            const anuga_int k3 = 3 * k;
            bed_ev_w[k3 + 0] += dz_cell;
            bed_ev_w[k3 + 1] += dz_cell;
            bed_ev_w[k3 + 2] += dz_cell;
        }
    }
}

// ============================================================================
// Update conserved quantities
// ============================================================================

void core_update_conserved_quantities(struct domain *D, double timestep) {
    anuga_int n = D->number_of_elements;
    const anuga_int n_tracers = D->number_of_tracers;

    double * restrict stage_cv = D->stage_centroid_values;
    double * restrict xmom_cv = D->xmom_centroid_values;
    double * restrict ymom_cv = D->ymom_centroid_values;

    double * restrict stage_eu = D->stage_explicit_update;
    double * restrict xmom_eu = D->xmom_explicit_update;
    double * restrict ymom_eu = D->ymom_explicit_update;

    double * restrict stage_siu = D->stage_semi_implicit_update;
    double * restrict xmom_siu = D->xmom_semi_implicit_update;
    double * restrict ymom_siu = D->ymom_semi_implicit_update;

    // Tracer pointers are hoisted to FUNCTION SCOPE here, not loaded inside the
    // n_tracers > 0 guard as in the flux kernel. On a GPU build OMP_PARALLEL_LOOP
    // is 'omp target teams loop', and D itself is NOT mapped to the device, so a
    // D->member load inside the loop reads a host address on the device: the
    // tracer update silently does nothing (m never changes, while the flux
    // kernel still fills explicit_update). Hoisting lets the pointer values be
    // captured as firstprivate scalars and address-translated. The flux kernel's
    // in-guard loading is a CPU hot-loop optimisation (HANDOVER 2.4, +2.26%% at
    // Ns=0) and does not apply to these much cheaper elementwise loops.
#ifndef CPU_ONLY_MODE
    double * restrict t_cons = D->tracer_conserved_values;
    double * restrict t_eu   = D->tracer_explicit_update;
#endif

    OMP_PARALLEL_LOOP
    for (anuga_int k = 0; k < n; k++) {
        // Get current centroid values
        double stage_c = stage_cv[k];
        double xmom_c = xmom_cv[k];
        double ymom_c = ymom_cv[k];

        // Apply explicit updates
        double stage_new = stage_c + timestep * stage_eu[k];
        double xmom_new  = xmom_c  + timestep * xmom_eu[k];
        double ymom_new  = ymom_c  + timestep * ymom_eu[k];

        // Apply semi-implicit updates, reformulated to ONE division per quantity.
        // The original did two FP64 divisions per quantity (si = siu/c, then cv/denom
        // with denom = 1 - dt*si); algebraically
        //     cv / (1 - dt*siu/c)  ==  cv*c / (c - dt*siu),
        // so num = c - dt*siu = denom*c, and denom>0  <=>  num*c > 0. Halving the
        // divisions matters on GeForce GPUs, where FP64 is 1/64 rate and ncu shows
        // this kernel FP64-pipe-bound (see issue #199); mathematically identical, so
        // results differ only at floating-point roundoff.
        double num;

        num = stage_c - timestep * stage_siu[k];
        if (stage_c != 0.0 && num * stage_c > 0.0) stage_new = stage_new * stage_c / num;

        num = xmom_c - timestep * xmom_siu[k];
        if (xmom_c != 0.0 && num * xmom_c > 0.0) xmom_new = xmom_new * xmom_c / num;

        num = ymom_c - timestep * ymom_siu[k];
        if (ymom_c != 0.0 && num * ymom_c > 0.0) ymom_new = ymom_new * ymom_c / num;

        stage_cv[k] = stage_new;
        xmom_cv[k] = xmom_new;
        ymom_cv[k] = ymom_new;

        // Reset semi-implicit updates for next timestep
        stage_siu[k] = 0.0;
        xmom_siu[k] = 0.0;
        ymom_siu[k] = 0.0;

        // Tracers: integrate the conserved m = h*c. No semi-implicit term and
        // deliberately NO clamping -- clamping would break exact conservation.
        // Positivity is instead a property of the upwind flux under CFL, and is
        // asserted by the tests rather than enforced here.
        if (n_tracers > 0) {
#ifdef CPU_ONLY_MODE
            double * restrict t_cons = D->tracer_conserved_values;
            double * restrict t_eu   = D->tracer_explicit_update;
#endif
            for (anuga_int s = 0; s < n_tracers; s++) {
                t_cons[s * n + k] += timestep * t_eu[s * n + k];
            }
        }
    }
}

#pragma omp declare target
// One cell's Manning friction + conserved-quantity update + optional RK2
// average, entirely in registers.  Shared by core_forcing_and_update (which
// reads the explicit updates from the eu arrays) and by the edge-based
// core_flux_apply_and_update (which computes them in registers and never
// touches the eu arrays at all).  eu_* are the explicit-update values for
// this cell; the semi-implicit arrays are read, consumed and reset here.
static inline void gpu_cell_forcing_update(
    anuga_int k, double timestep, int apply_manning, int do_saxpy,
    double a, double b, double g, double minimum_allowed_height,
    double seven_thirds,
    double eu_stage, double eu_xmom, double eu_ymom,
    double * restrict stage_cv, double * restrict xmom_cv,
    double * restrict ymom_cv, double * restrict bed_cv,
    double * restrict height_cv, double * restrict friction_cv,
    double * restrict stage_siu, double * restrict xmom_siu,
    double * restrict ymom_siu,
    double * restrict stage_bk, double * restrict xmom_bk,
    double * restrict ymom_bk) {

    double stage_c = stage_cv[k];
    double xmom_c = xmom_cv[k];
    double ymom_c = ymom_cv[k];

    double s_siu = stage_siu[k];
    double x_siu = xmom_siu[k];
    double y_siu = ymom_siu[k];

    if (apply_manning) {
        double S = 0.0;
        double eta = friction_cv[k];
        double abs_mom = sqrt(xmom_c * xmom_c + ymom_c * ymom_c);

        if (eta > 1.0e-15) {  // ETA_SMALL
            double h = stage_c - bed_cv[k];
            if (h >= minimum_allowed_height) {
                S = -g * eta * eta * abs_mom;
                S /= pow(h, seven_thirds);
            }
        }
        x_siu += S * xmom_c;
        y_siu += S * ymom_c;
    }

    // Explicit + semi-implicit update (single-division form; see
    // core_update_conserved_quantities for the derivation)
    double stage_new = stage_c + timestep * eu_stage;
    double xmom_new  = xmom_c  + timestep * eu_xmom;
    double ymom_new  = ymom_c  + timestep * eu_ymom;

    double num;

    num = stage_c - timestep * s_siu;
    if (stage_c != 0.0 && num * stage_c > 0.0) stage_new = stage_new * stage_c / num;

    num = xmom_c - timestep * x_siu;
    if (xmom_c != 0.0 && num * xmom_c > 0.0) xmom_new = xmom_new * xmom_c / num;

    num = ymom_c - timestep * y_siu;
    if (ymom_c != 0.0 && num * ymom_c > 0.0) ymom_new = ymom_new * ymom_c / num;

    stage_siu[k] = 0.0;
    xmom_siu[k] = 0.0;
    ymom_siu[k] = 0.0;

    if (do_saxpy) {
        stage_new = a * stage_new + b * stage_bk[k];
        xmom_new  = a * xmom_new  + b * xmom_bk[k];
        ymom_new  = a * ymom_new  + b * ymom_bk[k];
        height_cv[k] = fmax(stage_new - bed_cv[k], 0.0);
    }

    stage_cv[k] = stage_new;
    xmom_cv[k] = xmom_new;
    ymom_cv[k] = ymom_new;
}
#pragma omp end declare target

// ============================================================================
// Fused forcing + update (+ optional RK2 average)
//
// Manning friction, the conserved-quantity update and the RK2 average are all
// strictly cell-local: each reads and writes only index k.  Running them as
// three separate kernels means three launches and three round trips through
// the semi-implicit and centroid arrays, so they are fused here into one.
//
// This is as far as fusion goes in the DE step.  compute_fluxes cannot join
// them: it is a stencil kernel -- it reads height_cv[neighbour],
// bed_cv[neighbour], stage_cv[neighbour] and the neighbours' edge values -- so
// writing any centroid value from inside it would race against another team
// still reading that value, and an `omp target teams loop` has no device-wide
// barrier to order them.  The same argument rules out fusing extrapolate into
// compute_fluxes.
//
// Only the FLAT Manning variant is inlined here.  Callers with
// use_sloped_mannings must keep calling the sloped kernel separately (it reads
// vertex values, which the GPU path does not map).
//
//   timestep     dt to apply (must already be known -- fine for RK2 substep 2,
//                which reuses substep 1's dt)
//   apply_manning  1 => add the flat Manning friction term
//   do_saxpy     1 => finish with Q = a*Q + b*Q_backup and refresh height_cv
// ============================================================================

void core_forcing_and_update_on(struct domain *D, double timestep,
                                int apply_manning, int do_saxpy,
                                double a, double b,
                                const anuga_int * restrict iter, anuga_int iter_n) {
    anuga_int n = D->number_of_elements;
    double g = D->g;
    double minimum_allowed_height = D->minimum_allowed_height;
    double seven_thirds = ANUGA_SEVEN_THIRDS;

    double * restrict stage_cv = D->stage_centroid_values;
    double * restrict xmom_cv = D->xmom_centroid_values;
    double * restrict ymom_cv = D->ymom_centroid_values;
    double * restrict bed_cv = D->bed_centroid_values;
    double * restrict height_cv = D->height_centroid_values;
    double * restrict friction_cv = D->friction_centroid_values;

    double * restrict stage_eu = D->stage_explicit_update;
    double * restrict xmom_eu = D->xmom_explicit_update;
    double * restrict ymom_eu = D->ymom_explicit_update;

    double * restrict stage_siu = D->stage_semi_implicit_update;
    double * restrict xmom_siu = D->xmom_semi_implicit_update;
    double * restrict ymom_siu = D->ymom_semi_implicit_update;

    double * restrict stage_bk = D->stage_backup_values;
    double * restrict xmom_bk = D->xmom_backup_values;
    double * restrict ymom_bk = D->ymom_backup_values;

    // Generic passive tracers: the conserved m = h*c is integrated exactly as
    // core_update_conserved_quantities does it (no semi-implicit term, no
    // clamping) and RK2-averaged exactly as core_saxpy_conserved_quantities
    // does it, on m rather than c.  Pointers hoisted for the GPU build.
    const anuga_int n_tracers = D->number_of_tracers;
#ifndef CPU_ONLY_MODE
    double * restrict t_cons = D->tracer_conserved_values;
    double * restrict t_eu   = D->tracer_explicit_update;
    double * restrict t_bk   = D->tracer_backup_values;
#endif

    const anuga_int loop_n = iter ? iter_n : n;
    OMP_PARALLEL_LOOP
    for (anuga_int q = 0; q < loop_n; q++) {
        const anuga_int k = iter ? iter[q] : q;
        gpu_cell_forcing_update(k, timestep, apply_manning, do_saxpy, a, b,
                                g, minimum_allowed_height, seven_thirds,
                                stage_eu[k], xmom_eu[k], ymom_eu[k],
                                stage_cv, xmom_cv, ymom_cv, bed_cv, height_cv,
                                friction_cv, stage_siu, xmom_siu, ymom_siu,
                                stage_bk, xmom_bk, ymom_bk);

        if (n_tracers > 0) {
#ifdef CPU_ONLY_MODE
            double * restrict t_cons = D->tracer_conserved_values;
            double * restrict t_eu   = D->tracer_explicit_update;
            double * restrict t_bk   = D->tracer_backup_values;
#endif
            for (anuga_int s = 0; s < n_tracers; s++) {
                const anuga_int idx = s * n + k;
                double m = t_cons[idx] + timestep * t_eu[idx];
                if (do_saxpy) m = a * m + b * t_bk[idx];
                t_cons[idx] = m;
            }
        }
    }
}

void core_forcing_and_update(struct domain *D, double timestep,
                             int apply_manning, int do_saxpy,
                             double a, double b) {
    core_forcing_and_update_on(D, timestep, apply_manning, do_saxpy, a, b, NULL, 0);
}

// ============================================================================
// Backup conserved quantities for RK2
// ============================================================================

void core_backup_conserved_quantities(struct domain *D) {
    anuga_int n = D->number_of_elements;
    const anuga_int n_tracers = D->number_of_tracers;

    double * restrict stage_cv = D->stage_centroid_values;
    double * restrict xmom_cv = D->xmom_centroid_values;
    double * restrict ymom_cv = D->ymom_centroid_values;

    double * restrict stage_bk = D->stage_backup_values;
    double * restrict xmom_bk = D->xmom_backup_values;
    double * restrict ymom_bk = D->ymom_backup_values;

    // Tracer pointers are hoisted to FUNCTION SCOPE here, not loaded inside the
    // n_tracers > 0 guard as in the flux kernel. On a GPU build OMP_PARALLEL_LOOP
    // is 'omp target teams loop', and D itself is NOT mapped to the device, so a
    // D->member load inside the loop reads a host address on the device: the
    // tracer update silently does nothing (m never changes, while the flux
    // kernel still fills explicit_update). Hoisting lets the pointer values be
    // captured as firstprivate scalars and address-translated. The flux kernel's
    // in-guard loading is a CPU hot-loop optimisation (HANDOVER 2.4, +2.26%% at
    // Ns=0) and does not apply to these much cheaper elementwise loops.
#ifndef CPU_ONLY_MODE
    double * restrict t_cons = D->tracer_conserved_values;
    double * restrict t_bk   = D->tracer_backup_values;
#endif

    OMP_PARALLEL_LOOP
    for (anuga_int k = 0; k < n; k++) {
        stage_bk[k] = stage_cv[k];
        xmom_bk[k] = xmom_cv[k];
        ymom_bk[k] = ymom_cv[k];

        if (n_tracers > 0) {
#ifdef CPU_ONLY_MODE
            double * restrict t_cons = D->tracer_conserved_values;
            double * restrict t_bk   = D->tracer_backup_values;
#endif
            for (anuga_int s = 0; s < n_tracers; s++) {
                t_bk[s * n + k] = t_cons[s * n + k];
            }
        }
    }
}

// ============================================================================
// SAXPY for RK2/RK3: Q = (a*Q + b*Q_backup) / c
// ============================================================================

void core_saxpy_conserved_quantities(struct domain *D, double a, double b, double c) {
    anuga_int n = D->number_of_elements;
    const anuga_int n_tracers = D->number_of_tracers;

    double * restrict stage_cv = D->stage_centroid_values;
    double * restrict xmom_cv = D->xmom_centroid_values;
    double * restrict ymom_cv = D->ymom_centroid_values;

    double * restrict stage_bk = D->stage_backup_values;
    double * restrict xmom_bk = D->xmom_backup_values;
    double * restrict ymom_bk = D->ymom_backup_values;

    // Tracer pointers are hoisted to FUNCTION SCOPE here, not loaded inside the
    // n_tracers > 0 guard as in the flux kernel. On a GPU build OMP_PARALLEL_LOOP
    // is 'omp target teams loop', and D itself is NOT mapped to the device, so a
    // D->member load inside the loop reads a host address on the device: the
    // tracer update silently does nothing (m never changes, while the flux
    // kernel still fills explicit_update). Hoisting lets the pointer values be
    // captured as firstprivate scalars and address-translated. The flux kernel's
    // in-guard loading is a CPU hot-loop optimisation (HANDOVER 2.4, +2.26%% at
    // Ns=0) and does not apply to these much cheaper elementwise loops.
#ifndef CPU_ONLY_MODE
    double * restrict t_cons = D->tracer_conserved_values;
    double * restrict t_bk   = D->tracer_backup_values;
#endif

    // Standard SAXPY: Q = a*Q + b*Q_backup
    OMP_PARALLEL_LOOP
    for (anuga_int k = 0; k < n; k++) {
        stage_cv[k] = a * stage_cv[k] + b * stage_bk[k];
        xmom_cv[k] = a * xmom_cv[k] + b * xmom_bk[k];
        ymom_cv[k] = a * ymom_cv[k] + b * ymom_bk[k];

        // SAXPY must act on the CONSERVED m, not on c: h differs between RK
        // stages, so averaging c would not average the transported mass.
        if (n_tracers > 0) {
#ifdef CPU_ONLY_MODE
            double * restrict t_cons = D->tracer_conserved_values;
            double * restrict t_bk   = D->tracer_backup_values;
#endif
            for (anuga_int s = 0; s < n_tracers; s++) {
                t_cons[s * n + k] = a * t_cons[s * n + k] + b * t_bk[s * n + k];
            }
        }
    }

    // Apply c scaling if needed: Q = Q / c
    // Used for numerical stability with RK coefficients like a=1/3, b=2/3
    // Skip if c=0.0 (RK2 passes 0.0) or c=1.0 (no scaling needed)
    if (c != 1.0 && c != 0.0) {
        double c_inv = 1.0 / c;
        OMP_PARALLEL_LOOP
        for (anuga_int k = 0; k < n; k++) {
            stage_cv[k] *= c_inv;
            xmom_cv[k] *= c_inv;
            ymom_cv[k] *= c_inv;
            if (n_tracers > 0) {
#ifdef CPU_ONLY_MODE
                double * restrict t_cons = D->tracer_conserved_values;
#endif
                for (anuga_int s = 0; s < n_tracers; s++) t_cons[s * n + k] *= c_inv;
            }
        }
    }
}

// ============================================================================
// Protect against negative depths
// ============================================================================

double core_protect(struct domain *D) {
    anuga_int n = D->number_of_elements;
    double minimum_allowed_height = D->minimum_allowed_height;

    double * restrict stage_cv = D->stage_centroid_values;
    double * restrict xmom_cv = D->xmom_centroid_values;
    double * restrict ymom_cv = D->ymom_centroid_values;
    double * restrict bed_cv = D->bed_centroid_values;
    anuga_geom_t * restrict areas = D->areas;

    double mass_error = 0.0;

    OMP_PARALLEL_LOOP_REDUCTION_PLUS(mass_error)
    for (anuga_int k = 0; k < n; k++) {
        double h = stage_cv[k] - bed_cv[k];

        if (h < minimum_allowed_height) {
            // Very shallow - zero momentum to prevent instability
            xmom_cv[k] = 0.0;
            ymom_cv[k] = 0.0;
        }

        if (h < 0.0) {
            // Negative depth - track mass error and set stage to bed
            mass_error += (-h) * areas[k];
            stage_cv[k] = bed_cv[k];
        }
    }

    return mass_error;
}

// ============================================================================
// Fused step preparation: RK2 backup + protect + extrapolate centroid pass.
//
// All three touch only index k, so they run as ONE kernel: the centroid values
// are read once into registers, backed up, protected, and converted for the
// edge pass without three separate trips through memory.  This also retires
// the standalone protect's follow-up height refresh -- the centroid pass
// recomputes height_cv from the protected stage anyway.
//
// The sequencing inside the loop body reproduces the original kernel order
// (backup BEFORE protect -- the RK2 average must combine with the unprotected
// state, exactly as gpu_backup_conserved_quantities did) so results are
// bit-identical to the unfused sequence.
//
// Returns the protect mass error (same reduction core_protect performs).
// ============================================================================

double core_prepare_step_on(struct domain *D, int do_backup, int zero_eu,
                            const anuga_int * restrict iter, anuga_int iter_n) {
    anuga_int n = D->number_of_elements;
    double minimum_allowed_height = D->minimum_allowed_height;
    anuga_int extrapolate_velocity_second_order = D->extrapolate_velocity_second_order;

    double * restrict stage_cv = D->stage_centroid_values;
    double * restrict xmom_cv = D->xmom_centroid_values;
    double * restrict ymom_cv = D->ymom_centroid_values;
    double * restrict bed_cv = D->bed_centroid_values;
    double * restrict height_cv = D->height_centroid_values;
    anuga_geom_t * restrict areas = D->areas;
    double * restrict x_centroid_work = D->x_centroid_work;
    double * restrict y_centroid_work = D->y_centroid_work;

    double * restrict stage_bk = D->stage_backup_values;
    double * restrict xmom_bk = D->xmom_backup_values;
    double * restrict ymom_bk = D->ymom_backup_values;

    // Scatter-mode fluxes accumulate into the explicit updates with atomics,
    // so they must start the step at zero; the cell-based flux kernel
    // initializes them itself and passes zero_eu = 0.
    double * restrict stage_eu = D->stage_explicit_update;
    double * restrict xmom_eu = D->xmom_explicit_update;
    double * restrict ymom_eu = D->ymom_explicit_update;

    // Generic passive tracers: the same three cell-local pieces the unfused
    // kernels do (core_backup_conserved_quantities backs up m, the flux
    // kernel zeroes the explicit update, core_extrapolate_centroid_pass
    // derives c = m/h from the PROTECTED height).  Pointers hoisted to
    // function scope for the GPU build -- see core_update_conserved_quantities.
    const anuga_int n_tracers = D->number_of_tracers;
#ifndef CPU_ONLY_MODE
    double * restrict t_cons = D->tracer_conserved_values;
    double * restrict t_cv   = D->tracer_centroid_values;
    double * restrict t_bk   = D->tracer_backup_values;
    double * restrict t_eu   = D->tracer_explicit_update;
#endif

    double mass_error = 0.0;
    const anuga_int loop_n = iter ? iter_n : n;

    OMP_PARALLEL_LOOP_REDUCTION_PLUS(mass_error)
    for (anuga_int q = 0; q < loop_n; q++) {
        const anuga_int k = iter ? iter[q] : q;
        double stage = stage_cv[k];
        double bed = bed_cv[k];
        double xmom = xmom_cv[k];
        double ymom = ymom_cv[k];

        if (zero_eu) {
            stage_eu[k] = 0.0;
            xmom_eu[k] = 0.0;
            ymom_eu[k] = 0.0;
        }

        // RK2 backup of the raw (pre-protect) state
        if (do_backup) {
            stage_bk[k] = stage;
            xmom_bk[k] = xmom;
            ymom_bk[k] = ymom;
        }

        // Protect (core_protect's logic, in registers)
        double h = stage - bed;
        if (h < minimum_allowed_height) {
            xmom = 0.0;
            ymom = 0.0;
        }
        if (h < 0.0) {
            mass_error += (-h) * areas[k];
            stage = bed;
        }
        stage_cv[k] = stage;

        // Extrapolate centroid pass (velocity into the work arrays)
        double dk = fmax(stage - bed, 0.0);
        height_cv[k] = dk;

        int is_dry = (dk <= minimum_allowed_height);
        int extrapolate = (extrapolate_velocity_second_order == 1) && (dk > minimum_allowed_height);

        double xmom_out = is_dry ? 0.0 : xmom;
        double ymom_out = is_dry ? 0.0 : ymom;

        double inv_dk = extrapolate ? (1.0 / dk) : 1.0;

        x_centroid_work[k] = xmom_out * inv_dk;
        y_centroid_work[k] = ymom_out * inv_dk;

        xmom_cv[k] = xmom_out;
        ymom_cv[k] = ymom_out;

        if (n_tracers > 0) {
#ifdef CPU_ONLY_MODE
            double * restrict t_cons = D->tracer_conserved_values;
            double * restrict t_cv   = D->tracer_centroid_values;
            double * restrict t_bk   = D->tracer_backup_values;
            double * restrict t_eu   = D->tracer_explicit_update;
#endif
            const double inv_h = is_dry ? 0.0 : (1.0 / dk);
            for (anuga_int s = 0; s < n_tracers; s++) {
                const anuga_int idx = s * n + k;
                const double m = t_cons[idx];
                if (do_backup) t_bk[idx] = m;
                if (zero_eu)   t_eu[idx] = 0.0;
                t_cv[idx] = m * inv_h;
            }
        }
    }

    return mass_error;
}

double core_prepare_step(struct domain *D, int do_backup, int zero_eu) {
    return core_prepare_step_on(D, do_backup, zero_eu, NULL, 0);
}

// ============================================================================
// Fix negative cells
//
// Matches _openmp_fix_negative_cells (the tested CPU reference):
//   - Only acts on cells where stage - bed < 0  AND  tri_full_flag > 0
//     (ghost cells are skipped, matching the openmp & bitwise-and condition)
//   - Zeros xmom/ymom and resets stage to bed for those cells
//   - Returns count of cells fixed (parallel + reduction)
//
// NOTE: The original core version (before unification) used a different
// threshold (minimum_allowed_height) and ignored tri_full_flag — it has
// been updated here to match the _openmp_ reference behaviour exactly.
// ============================================================================

int core_fix_negative_cells(struct domain *D) {
    anuga_int n = D->number_of_elements;

    double * restrict stage_cv = D->stage_centroid_values;
    double * restrict xmom_cv  = D->xmom_centroid_values;
    double * restrict ymom_cv  = D->ymom_centroid_values;
    double * restrict bed_cv   = D->bed_centroid_values;
    anuga_int * restrict tri_full_flag = D->tri_full_flag;

    int num_negative_cells = 0;

    OMP_PARALLEL_LOOP_REDUCTION_PLUS(num_negative_cells)
    for (anuga_int k = 0; k < n; k++) {
        // Use & (bitwise and) matching the original _openmp_ condition.
        // tri_full_flag is always initialised to ones(N) so the pointer is
        // never NULL when called from Cython; the check avoids UB for the
        // standalone / GPU build path where it could theoretically be NULL.
        int full = (tri_full_flag == NULL) ? 1 : (tri_full_flag[k] > 0);
        if ((stage_cv[k] - bed_cv[k] < 0.0) & full) {
            num_negative_cells = num_negative_cells + 1;
            stage_cv[k] = bed_cv[k];
            xmom_cv[k]  = 0.0;
            ymom_cv[k]  = 0.0;
        }
    }

    return num_negative_cells;
}

// ============================================================================
// Negative-cell volume (read-only)
//
// Measures the water volume that fix_negative_cells will ADD by clamping
// negative-depth cells up to zero depth (stage = bed) — i.e. the conservation
// error the clamp introduces this step. Uses the SAME cell selection as
// core_fix_negative_cells (stage - bed < 0 AND tri_full_flag > 0), so it must
// be called AFTER the flux update but BEFORE core_fix_negative_cells (which
// erases the deficit). Does not modify the domain.
// ============================================================================

double core_negative_cells_volume(struct domain *D) {
    anuga_int n = D->number_of_elements;

    double * restrict stage_cv = D->stage_centroid_values;
    double * restrict bed_cv   = D->bed_centroid_values;
    anuga_geom_t * restrict areas    = D->areas;
    anuga_int * restrict tri_full_flag = D->tri_full_flag;

    double volume = 0.0;

    OMP_PARALLEL_LOOP_REDUCTION_PLUS(volume)
    for (anuga_int k = 0; k < n; k++) {
        int full = (tri_full_flag == NULL) ? 1 : (tri_full_flag[k] > 0);
        if ((stage_cv[k] - bed_cv[k] < 0.0) & full) {
            // bed - stage > 0 here: volume needed to raise the cell to zero depth
            volume = volume + (bed_cv[k] - stage_cv[k]) * areas[k];
        }
    }

    return volume;
}

// ============================================================================
// Manning friction (flat, semi-implicit)
// ============================================================================

void core_manning_friction_flat_semi_implicit(struct domain *D) {
    anuga_int n = D->number_of_elements;
    double g = D->g;
    double minimum_allowed_height = D->minimum_allowed_height;
    double seven_thirds = ANUGA_SEVEN_THIRDS;

    double * restrict stage_cv = D->stage_centroid_values;
    double * restrict bed_cv = D->bed_centroid_values;
    double * restrict xmom_cv = D->xmom_centroid_values;
    double * restrict ymom_cv = D->ymom_centroid_values;
    double * restrict friction_cv = D->friction_centroid_values;

    double * restrict xmom_siu = D->xmom_semi_implicit_update;
    double * restrict ymom_siu = D->ymom_semi_implicit_update;

    OMP_PARALLEL_LOOP
    for (anuga_int k = 0; k < n; k++) {
        double S = 0.0;
        double uh = xmom_cv[k];
        double vh = ymom_cv[k];
        double eta = friction_cv[k];
        double abs_mom = sqrt(uh * uh + vh * vh);

        if (eta > 1.0e-15) {  // ETA_SMALL
            double h = stage_cv[k] - bed_cv[k];
            if (h >= minimum_allowed_height) {
                S = -g * eta * eta * abs_mom;
                S /= pow(h, seven_thirds);
            }
        }
        xmom_siu[k] += S * uh;
        ymom_siu[k] += S * vh;
    }
}

// ============================================================================
// Vegetation drag (Baptist et al. 2007), semi-implicit   spec 8
// ============================================================================
//
// The drag of a field of stems on the flow, as a friction slope
//
//   d(uh)/dt = -g uh |q| / (Cv^2 h^2)
//
// with the vegetated Chezy coefficient of Baptist et al. (2007), J. Hydraul.
// Res. 45, Eq (16):
//
//   Cv = (Cb^-2 + Cd m D min(h, hv) / (2 g))^-1/2
//        + (sqrt(g) / kappa) ln(max(h, hv) / hv)
//
// Cb the bed Chezy coefficient, Cd the stem drag coefficient, m the stem
// density (1/m^2), D the stem diameter (m), hv the stem height (m). Emergent
// stems (h <= hv) give the first term alone; submerged stems add the
// logarithmic overflow term. Same semi-implicit treatment as Manning: the
// slope is accumulated into the semi-implicit update and applied as a
// division in update_conserved_quantities, so it can only ever reduce the
// momentum, never reverse it -- the constraint spec 8.2 asks for.
//
// Cells with no stems (m <= 0) are untouched, so the field lives alongside
// Manning friction: open water carries n, vegetation carries (m, D, hv)
// with n = 0, as in the Delta-X Wax Lake model this was written for.
//
// Cell-local, so it offloads like the Manning kernels.
void core_vegetation_friction_semi_implicit(struct domain *D) {
    if (D->vegetation_mode == 0) {
        return;
    }
    const anuga_int n = D->number_of_elements;
    const double g = D->g;
    const double minimum_allowed_height = D->minimum_allowed_height;
    const double a1 = 1.0 / (D->vegetation_bed_chezy * D->vegetation_bed_chezy);
    const double Cd_over_2g = D->vegetation_Cd / (2.0 * g);
    const double a3 = sqrt(g) / 0.4;      /* kappa = 0.4 */

    double * restrict stage_cv = D->stage_centroid_values;
    double * restrict bed_cv = D->bed_centroid_values;
    double * restrict xmom_cv = D->xmom_centroid_values;
    double * restrict ymom_cv = D->ymom_centroid_values;
    double * restrict m_cv = D->veg_density_centroid_values;
    double * restrict D_cv = D->veg_diameter_centroid_values;
    double * restrict hv_cv = D->veg_height_centroid_values;
    double * restrict xmom_siu = D->xmom_semi_implicit_update;
    double * restrict ymom_siu = D->ymom_semi_implicit_update;

    OMP_PARALLEL_LOOP
    for (anuga_int k = 0; k < n; k++) {
        const double m = m_cv[k];
        if (!(m > 0.0)) continue;
        const double h = stage_cv[k] - bed_cv[k];
        if (h < minimum_allowed_height) continue;
        const double hv = hv_cv[k];
        if (!(hv > 0.0)) continue;
        const double uh = xmom_cv[k];
        const double vh = ymom_cv[k];
        const double abs_mom = sqrt(uh * uh + vh * vh);
        const double a2 = D_cv[k] * m * Cd_over_2g;
        double Cv = 1.0 / sqrt(a1 + a2 * fmin(h, hv));
        if (h > hv) {
            Cv += a3 * log(h / hv);
        }
        const double S = -g * abs_mom / (Cv * Cv * h * h);
        xmom_siu[k] += S * uh;
        ymom_siu[k] += S * vh;
    }
}

// ============================================================================
// Manning friction (sloped, semi-implicit)
// ============================================================================

void core_manning_friction_sloped_semi_implicit(struct domain *D) {
    anuga_int n = D->number_of_elements;
    double g = D->g;
    double minimum_allowed_height = D->minimum_allowed_height;

    double * restrict height_cv = D->height_centroid_values;
    double * restrict xmom_cv = D->xmom_centroid_values;
    double * restrict ymom_cv = D->ymom_centroid_values;
    double * restrict friction_cv = D->friction_centroid_values;
    double * restrict bed_vv = D->bed_vertex_values;
    anuga_geom_t * restrict vertex_coords = D->vertex_coordinates;

    double * restrict xmom_siu = D->xmom_semi_implicit_update;
    double * restrict ymom_siu = D->ymom_semi_implicit_update;

    OMP_PARALLEL_LOOP
    for (anuga_int k = 0; k < n; k++) {
        double h = height_cv[k];

        if (h > minimum_allowed_height) {
            anuga_int k3 = k * 3;
            anuga_int k6 = k * 6;

            // Compute bed slope
            double x0 = vertex_coords[k6 + 0];
            double y0 = vertex_coords[k6 + 1];
            double x1 = vertex_coords[k6 + 2];
            double y1 = vertex_coords[k6 + 3];
            double x2 = vertex_coords[k6 + 4];
            double y2 = vertex_coords[k6 + 5];

            double z0 = bed_vv[k3 + 0];
            double z1 = bed_vv[k3 + 1];
            double z2 = bed_vv[k3 + 2];

            double det = (y2 - y0) * (x1 - x0) - (y1 - y0) * (x2 - x0);
            double dzx = ((y2 - y0) * (z1 - z0) - (y1 - y0) * (z2 - z0)) / det;
            double dzy = ((x1 - x0) * (z2 - z0) - (x2 - x0) * (z1 - z0)) / det;

            double slope = sqrt(1.0 + dzx * dzx + dzy * dzy);

            double eta = friction_cv[k];
            double xmom = xmom_cv[k];
            double ymom = ymom_cv[k];

            double S = -g * eta * eta * sqrt(xmom * xmom + ymom * ymom) * slope;
            S /= pow(h, ANUGA_SEVEN_THIRDS);

            xmom_siu[k] += S;
            ymom_siu[k] += S;
        }
    }
}

// ============================================================================
// Manning friction (sloped, semi-implicit, edge-based)
//
// Like core_manning_friction_sloped_semi_implicit but derives the bed slope
// from edge values (bed_edge_values) instead of vertex values.  This is the
// active per-timestep path when domain.use_sloped_mannings=True
// (friction.py selects manning_friction_sloped_semi_implicit_edge_based).
// ============================================================================

void core_manning_friction_sloped_semi_implicit_edge_based(struct domain *D) {
    anuga_int n = D->number_of_elements;
    double g   = D->g;
    double eps = D->minimum_allowed_height;

    double * restrict stage_cv   = D->stage_centroid_values;
    double * restrict bed_ev     = D->bed_edge_values;
    double * restrict xmom_cv    = D->xmom_centroid_values;
    double * restrict ymom_cv    = D->ymom_centroid_values;
    double * restrict friction_cv = D->friction_centroid_values;
    anuga_geom_t * restrict edge_coords = D->edge_coordinates;

    double * restrict xmom_siu   = D->xmom_semi_implicit_update;
    double * restrict ymom_siu   = D->ymom_semi_implicit_update;

    const double one_third   = ANUGA_ONE_THIRD;
    const double seven_thirds = ANUGA_SEVEN_THIRDS;

    OMP_PARALLEL_LOOP
    for (anuga_int k = 0; k < n; k++) {
        double S = 0.0;
        double eta = friction_cv[k];

        if (eta > 1.0e-16) {
            anuga_int k3 = k * 3;
            anuga_int k6 = k * 6;

            // Bed values at edges
            double z0 = bed_ev[k3 + 0];
            double z1 = bed_ev[k3 + 1];
            double z2 = bed_ev[k3 + 2];

            // Edge midpoint coordinates
            double x0 = edge_coords[k6 + 0];
            double y0 = edge_coords[k6 + 1];
            double x1 = edge_coords[k6 + 2];
            double y1 = edge_coords[k6 + 3];
            double x2 = edge_coords[k6 + 4];
            double y2 = edge_coords[k6 + 5];

            // Bed slope via 2x2 determinant (same as _gradient(), inlined for GPU)
            double det = (y2 - y0) * (x1 - x0) - (y1 - y0) * (x2 - x0);
            double zx  = ((y2 - y0) * (z1 - z0) - (y1 - y0) * (z2 - z0)) / det;
            double zy  = ((x1 - x0) * (z2 - z0) - (x2 - x0) * (z1 - z0)) / det;

            double zs = sqrt(1.0 + zx * zx + zy * zy);
            double z  = (z0 + z1 + z2) * one_third;

            double w  = stage_cv[k];
            double h  = w - z;

            if (h >= eps) {
                double uh = xmom_cv[k];
                double vh = ymom_cv[k];
                S = -g * eta * eta * zs * sqrt(uh * uh + vh * vh);
                S /= pow(h, seven_thirds);
            }
        }

        xmom_siu[k] += S * xmom_cv[k];
        ymom_siu[k] += S * ymom_cv[k];
    }
}

// ============================================================================
// Gravity term
//
// Computes bed-slope gravity source term: duh/dt += -g * avg_h * dz/dx
// Uses stage_centroid - bed_centroid for avg_h (matches the original
// _openmp_gravity which computed this directly, so height need not be
// up-to-date when this function is called).
// ============================================================================

int core_gravity(struct domain *D) {
    anuga_int n = D->number_of_elements;
    double g = D->g;

    double * restrict stage_cv = D->stage_centroid_values;
    double * restrict bed_cv   = D->bed_centroid_values;
    double * restrict bed_vv   = D->bed_vertex_values;

    double * restrict xmom_eu = D->xmom_explicit_update;
    double * restrict ymom_eu = D->ymom_explicit_update;

    anuga_geom_t * restrict vertex_coords = D->vertex_coordinates;

    OMP_PARALLEL_LOOP
    for (anuga_int k = 0; k < n; k++) {
        // Average depth: use live stage - bed (height_cv may be stale)
        double avg_h = stage_cv[k] - bed_cv[k];

        anuga_int k3 = k * 3;
        anuga_int k6 = k * 6;

        double x0 = vertex_coords[k6 + 0];
        double y0 = vertex_coords[k6 + 1];
        double x1 = vertex_coords[k6 + 2];
        double y1 = vertex_coords[k6 + 3];
        double x2 = vertex_coords[k6 + 4];
        double y2 = vertex_coords[k6 + 5];

        double z0 = bed_vv[k3 + 0];
        double z1 = bed_vv[k3 + 1];
        double z2 = bed_vv[k3 + 2];

        // Bed gradient via 2x2 determinant (same as _gradient(), inlined for GPU)
        double det = (y2 - y0) * (x1 - x0) - (y1 - y0) * (x2 - x0);
        double dzx = ((y2 - y0) * (z1 - z0) - (y1 - y0) * (z2 - z0)) / det;
        double dzy = ((x1 - x0) * (z2 - z0) - (x2 - x0) * (z1 - z0)) / det;

        xmom_eu[k] += -g * avg_h * dzx;
        ymom_eu[k] += -g * avg_h * dzy;
    }

    return 0;
}

// ============================================================================
// Gravity term (well-balanced)
//
// Well-balanced formulation after Audusse et al. (2004):
//   du/dt += -g * wx * avg_h                    (stage-gradient term)
//   dv/dt += -g * wy * avg_h
//   PLUS side-pressure correction:
//     sum_i  -0.5 * g * h_i^2 * edgelength_i * n_i / area
// where h_i = stage_edge[i] - bed_edge[i] is the depth at edge i,
// and wx, wy is the gradient of stage (not bed), computed from vertex values.
//
// This formulation is exactly what _openmp_gravity_wb computed.
// Still-water equilibrium (u=v=0, stage=const) is preserved exactly
// because the stage-gradient term and edge-pressure terms cancel.
// ============================================================================

int core_gravity_wb(struct domain *D) {
    anuga_int n = D->number_of_elements;
    double g = D->g;

    double * restrict stage_vv  = D->stage_vertex_values;
    double * restrict stage_cv  = D->stage_centroid_values;
    double * restrict bed_cv    = D->bed_centroid_values;
    double * restrict stage_ev  = D->stage_edge_values;
    double * restrict bed_ev    = D->bed_edge_values;
    anuga_geom_t * restrict normals   = D->normals;
    anuga_geom_t * restrict edgelengths = D->edgelengths;
    anuga_geom_t * restrict areas     = D->areas;
    double * restrict xmom_eu   = D->xmom_explicit_update;
    double * restrict ymom_eu   = D->ymom_explicit_update;
    anuga_geom_t * restrict vertex_coords = D->vertex_coordinates;

    OMP_PARALLEL_LOOP
    for (anuga_int k = 0; k < n; k++) {
        anuga_int k3 = k * 3;
        anuga_int k6 = k * 6;

        // --------------------------------------------------
        // Stage-gradient term: -g * avg_h * (wx, wy)
        // --------------------------------------------------

        // Stage at vertices for gradient calculation
        double w0 = stage_vv[k3 + 0];
        double w1 = stage_vv[k3 + 1];
        double w2 = stage_vv[k3 + 2];

        // Vertex coordinates
        double x0 = vertex_coords[k6 + 0];
        double y0 = vertex_coords[k6 + 1];
        double x1 = vertex_coords[k6 + 2];
        double y1 = vertex_coords[k6 + 3];
        double x2 = vertex_coords[k6 + 4];
        double y2 = vertex_coords[k6 + 5];

        // Compute stage gradient using standard 2x2 determinant formula
        // (identical math to _gradient() in util_ext.h, inlined for GPU compat)
        double det = (y2 - y0) * (x1 - x0) - (y1 - y0) * (x2 - x0);
        double wx  = ((y2 - y0) * (w1 - w0) - (y1 - y0) * (w2 - w0)) / det;
        double wy  = ((x1 - x0) * (w2 - w0) - (x2 - x0) * (w1 - w0)) / det;

        // Centroid depth
        double avg_h = stage_cv[k] - bed_cv[k];

        // Apply stage-gradient term
        xmom_eu[k] += -g * wx * avg_h;
        ymom_eu[k] += -g * wy * avg_h;

        // --------------------------------------------------
        // Edge-pressure (side) correction:
        //   sum_i  -0.5 * g * h_i^2 * edgelength_i * n_i / area
        // --------------------------------------------------
        double sidex = 0.0;
        double sidey = 0.0;
        for (int i = 0; i < 3; i++) {
            double h_edge = stage_ev[k3 + i] - bed_ev[k3 + i];
            double fact   = -0.5 * g * h_edge * h_edge * edgelengths[k3 + i];
            sidex += fact * normals[k6 + 2 * i];
            sidey += fact * normals[k6 + 2 * i + 1];
        }

        double inv_area = 1.0 / areas[k];
        xmom_eu[k] += -sidex * inv_area;
        ymom_eu[k] += -sidey * inv_area;
    }

    return 0;
}

// ============================================================================
// Compute fluxes using central upwind scheme (UNIFIED CPU/GPU)
// ============================================================================

double core_compute_fluxes_central(struct domain *D, int substep_count, int timestep_fluxcalls) {
    anuga_int n = D->number_of_elements;
    double g = D->g;
    double epsilon = D->epsilon;
    anuga_int low_froude = D->low_froude;

    // Extract array pointers
    double * restrict stage_cv = D->stage_centroid_values;
    double * restrict bed_cv = D->bed_centroid_values;
    double * restrict height_cv = D->height_centroid_values;

    double * restrict stage_ev = D->stage_edge_values;
    double * restrict xmom_ev = D->xmom_edge_values;
    double * restrict ymom_ev = D->ymom_edge_values;
    double * restrict height_ev = D->height_edge_values;
    double * restrict bed_ev    = D->bed_edge_values;

    // Opt-in: reconstruct edge bed values as stage - height instead of loading
    // bed_ev.  core_extrapolate_edge_pass computes bed_ev with EXACTLY that
    // expression from exactly these arrays, so whenever fluxes follow an
    // extrapolate (every evolve step) the reconstruction is bit-identical and
    // this memory-bound kernel drops one gather per edge -- on both sides,
    // 6 scattered loads per cell.  It is wrong for callers that set edge
    // values independently and invoke fluxes directly (test_flux does), so it
    // stays off unless the driver guarantees the extrapolate-first contract
    // (D->reconstruct_edge_bed = 1; ANUGA leaves it 0).
    const int reconstruct_z = (D->reconstruct_edge_bed != 0);

    double * restrict stage_bv = D->stage_boundary_values;
    double * restrict xmom_bv = D->xmom_boundary_values;
    double * restrict ymom_bv = D->ymom_boundary_values;

    double * restrict stage_eu = D->stage_explicit_update;
    double * restrict xmom_eu = D->xmom_explicit_update;
    double * restrict ymom_eu = D->ymom_explicit_update;

    anuga_int * restrict neighbours = D->neighbours;
    anuga_int * restrict neighbour_edges = D->neighbour_edges;
    anuga_geom_t * restrict normals = D->normals;
    anuga_geom_t * restrict edgelengths = D->edgelengths;
    anuga_geom_t * restrict radii = D->radii;
    anuga_geom_t * restrict areas = D->areas;
    double * restrict max_speed_array = D->max_speed;
    anuga_int * restrict tri_full_flag = D->tri_full_flag;

    // Riverwall arrays (may be NULL if no riverwalls)
    anuga_int n_riverwall_edges = D->number_of_riverwall_edges;
    anuga_int ncol_riverwall_hp = D->ncol_riverwall_hydraulic_properties;
    anuga_int * restrict edge_flux_type = D->edge_flux_type;
    anuga_int * restrict edge_river_wall_counter = D->edge_river_wall_counter;
    double * restrict riverwall_elevation = D->riverwall_elevation;
    anuga_int * restrict riverwall_rowIndex = D->riverwall_rowIndex;
    double * restrict riverwall_hydraulic_properties = D->riverwall_hydraulic_properties;

    // Generic passive tracers.  n_tracers == 0 in every ordinary run;
    // all tracer work below is guarded on this loop-invariant integer.
    const anuga_int n_tracers = D->number_of_tracers;

// Tracer base pointers: WHERE they are loaded is build-dependent, and both
// choices are load-bearing.
//
//   CPU build  -- load them INSIDE the n_tracers > 0 guard. Hoisting them to
//                 function scope keeps them live across the hot loop and cost
//                 +2.26% on the CPU path at Ns=0 (HANDOVER.md 2.4). Only the
//                 benchmark caught that; every correctness test passed.
//   GPU build  -- the loops below are 'omp target' regions and D is NOT mapped
//                 to the device, so a D->member load inside the region reads a
//                 host address on the device. What happens next depends on
//                 what is loaded. A POINTER member yields a garbage pointer
//                 value, and the tracer work silently does nothing: no crash,
//                 explicit_update stays zero on the device, and m never moves.
//                 A SCALAR member is consumed as data, and the read itself
//                 faults: boundary_length below, loaded in-loop, took the whole
//                 process down with CUDA_ERROR_ILLEGAL_ADDRESS as soon as
//                 n_tracers > 0 (the load ran for every edge, before any
//                 inflow test). Everything read from D must be loaded at
//                 function scope so the VALUES are captured as firstprivate
//                 scalars and, for pointers, address-translated via the
//                 present table.
//
// So: hoisted declarations under #ifndef CPU_ONLY_MODE, in-guard declarations
// under #ifdef CPU_ONLY_MODE -- pointers AND the boundary_length scalar alike,
// so the benchmarked CPU hot loop stays byte-identical. The loop bodies are
// identical either way.
#ifndef CPU_ONLY_MODE
    double * restrict t_eu = D->tracer_explicit_update;
    double * restrict t_ev = D->tracer_edge_values;
    double * restrict t_bv = D->tracer_boundary_values;
    double * restrict t_bf = D->tracer_boundary_flux;
    double * restrict t_sf = D->tracer_speed_factor;
#endif
    // Scalar, so it is hoisted on BOTH builds: a D->member load inside the
    // element loop is a host dereference on the device (CUDA_ERROR_ILLEGAL_ADDRESS
    // the first time a boundary edge carries a tracer), and a loop-invariant
    // integer costs nothing on the CPU path.
    const anuga_int t_bl = D->boundary_length;

    // Reduction variables
    double local_timestep = 1.0e+100;
    double boundary_flux_sum_substep = 0.0;

    // Main flux computation loop with reductions
    #ifdef CPU_ONLY_MODE
    #pragma omp parallel for reduction(min:local_timestep) reduction(+:boundary_flux_sum_substep)
    #else
    #pragma omp target teams distribute parallel for reduction(min:local_timestep) reduction(+:boundary_flux_sum_substep)
    #endif
    for (anuga_int k = 0; k < n; k++) {
        double edgeflux[3];
        double ql[3], qr[3];
        double speed_max_last = 0.0;

        // Zero the explicit updates for this element
        stage_eu[k] = 0.0;
        xmom_eu[k] = 0.0;
        ymom_eu[k] = 0.0;
        if (n_tracers > 0) {
#ifdef CPU_ONLY_MODE
            double * restrict t_eu = D->tracer_explicit_update;
#endif
            for (anuga_int s = 0; s < n_tracers; s++) t_eu[s * n + k] = 0.0;
        }

        // Get centroid values for this element
        double hc = height_cv[k];
        double zc = bed_cv[k];

        // Loop over the 3 edges
        for (int i = 0; i < 3; i++) {
            int ki = 3 * k + i;
            int ki2 = 2 * ki;

            // Left state (this element's edge values); see reconstruct_z above
            ql[0] = stage_ev[ki];
            ql[1] = xmom_ev[ki];
            ql[2] = ymom_ev[ki];
            double hle = height_ev[ki];
            double zl = reconstruct_z ? (ql[0] - hle) : bed_ev[ki];

            // Edge geometry
            double length = edgelengths[ki];
            double n1 = normals[ki2];
            double n2 = normals[ki2 + 1];

            // Get neighbour info
            anuga_int neighbour = neighbours[ki];
            int is_boundary = (neighbour < 0);

            double zr, hre, hc_n, zc_n;

            if (is_boundary) {
                // Boundary edge - get values from boundary arrays
                int m = -neighbour - 1;
                qr[0] = stage_bv[m];
                qr[1] = xmom_bv[m];
                qr[2] = ymom_bv[m];
                zr = zl;
                hre = fmax(qr[0] - zr, 0.0);
                hc_n = hc;
                zc_n = zc;
            } else {
                // Internal edge - get values from neighbour element
                int m = neighbour_edges[ki];
                int nm = neighbour * 3 + m;
                qr[0] = stage_ev[nm];
                qr[1] = xmom_ev[nm];
                qr[2] = ymom_ev[nm];
                hre = height_ev[nm];
                zr = reconstruct_z ? (qr[0] - hre) : bed_ev[nm];
                hc_n = height_cv[neighbour];
                zc_n = bed_cv[neighbour];
            }

            // Compute z_half (max bed elevation at edge)
            double z_half = fmax(zl, zr);

            // Check for riverwall elevation override
            int is_riverwall = 0;
            double zwall = 0.0;
            if (n_riverwall_edges > 0 && edge_flux_type != NULL &&
                edge_river_wall_counter != NULL && riverwall_elevation != NULL &&
                edge_flux_type[ki] == 1) {
                int riverwall_index = edge_river_wall_counter[ki] - 1;
                if (riverwall_index >= 0) {
                    is_riverwall = 1;
                    zwall = riverwall_elevation[riverwall_index];
                    z_half = fmax(zwall, z_half);
                }
            }

            // Compute effective heights at the edge
            double h_left = fmax(hle + zl - z_half, 0.0);
            double h_right = fmax(hre + zr - z_half, 0.0);

            double max_speed_local = 0.0;
            double pressure_flux = 0.0;

            if (h_left == 0.0 && h_right == 0.0) {
                // Both heights zero - no flux
                edgeflux[0] = 0.0;
                edgeflux[1] = 0.0;
                edgeflux[2] = 0.0;
            } else {
                // Compute flux using central scheme
                gpu_flux_function_central(ql, qr,
                                          h_left, h_right,
                                          hle, hre,
                                          n1, n2,
                                          epsilon, z_half, g,
                                          edgeflux, &max_speed_local, &pressure_flux,
                                          low_froude);
            }

            // Apply riverwall weir discharge correction if applicable
            if (is_riverwall && zwall > fmax(zc, zc_n) &&
                riverwall_rowIndex != NULL && riverwall_hydraulic_properties != NULL) {
                // Get hydraulic properties for this riverwall
                anuga_int rw_count = edge_river_wall_counter[ki];
                anuga_int hp_row = riverwall_rowIndex[rw_count - 1];
                anuga_int ii = hp_row * ncol_riverwall_hp;

                double Qfactor = riverwall_hydraulic_properties[ii];
                double s1 = riverwall_hydraulic_properties[ii + 1];
                double s2 = riverwall_hydraulic_properties[ii + 2];
                double h1 = riverwall_hydraulic_properties[ii + 3];
                double h2 = riverwall_hydraulic_properties[ii + 4];
                // Column 5 is Cd_through; guard for old files with only 5 columns
                double Cd_through = (ncol_riverwall_hp > 5)
                    ? riverwall_hydraulic_properties[ii + 5]
                    : 0.0;

                // Weir height above minimum bed elevation
                double weir_height = fmax(zwall - fmin(zl, zr), 0.0);

                // Compute depths above weir using centroid values
                double h_left_weir = fmax(stage_cv[k] - z_half, 0.0);
                double h_right_weir = is_boundary
                    ? fmax(hc_n + zr - z_half, 0.0)
                    : fmax(stage_cv[neighbour] - z_half, 0.0);

                // Apply weir discharge correction (Villemonte overtopping)
                gpu_adjust_edgeflux_with_weir(edgeflux, h_left_weir, h_right_weir,
                                              g, weir_height, Qfactor,
                                              s1, s2, h1, h2, &max_speed_local);

                // Apply throughflow (orifice/seepage through wall body), additive
                double stage_left  = stage_cv[k];
                double stage_right = is_boundary
                    ? (hc_n + zr)
                    : stage_cv[neighbour];
                gpu_adjust_edgeflux_with_throughflow(
                    edgeflux,
                    stage_left, stage_right,
                    zl, zr,
                    zwall, g, Cd_through, &max_speed_local);
            }

            // Multiply flux by edge length (and negate for conservation)
            edgeflux[0] *= -length;
            edgeflux[1] *= -length;
            edgeflux[2] *= -length;

            // Track max speed for this element
            speed_max_last = fmax(speed_max_last, max_speed_local);

            // Accumulate flux contributions
            stage_eu[k] += edgeflux[0];
            xmom_eu[k] += edgeflux[1];
            ymom_eu[k] += edgeflux[2];

            // --- Passive tracer advection -------------------------------
            // edgeflux[0] is the water mass flux through this edge, already
            // multiplied by -length.  Sign convention after that negation:
            //     edgeflux[0] < 0  ->  OUTflow from k, donor is k
            //     edgeflux[0] > 0  ->  INflow  to   k, donor is the neighbour
            // Using the same edgeflux[0] for both cells sharing the edge makes
            // tracer mass conservation structural, independent of cell sizes.
            if (n_tracers > 0) {
#ifdef CPU_ONLY_MODE
                double * restrict t_ev = D->tracer_edge_values;
                double * restrict t_bv = D->tracer_boundary_values;
                double * restrict t_eu = D->tracer_explicit_update;
                double * restrict t_bf = D->tracer_boundary_flux;
                double * restrict t_sf = D->tracer_speed_factor;
#endif
                const double wflux = edgeflux[0];
                const int    inflow = (wflux > 0.0);
                /* [D-5v] the donor cell whose speed factor scales the flux;
                 * a boundary donor takes the receiving cell's. */
                const anuga_int donor = (inflow && !is_boundary) ? neighbour : k;
                /* Conservation accounting: record what crosses a DOMAIN boundary edge,
                 * on the same terms the water balance uses -- a real boundary, and this
                 * cell owned rather than a ghost. A ghost cell's copy of the edge
                 * belongs to its owner, so counting it here would double it. */
                const int count_bdry = (t_bf != NULL) && is_boundary
                    && (tri_full_flag == NULL || tri_full_flag[k] == 1);
                for (anuga_int s = 0; s < n_tracers; s++) {
                    double c_up;
                    if (inflow) {
                        c_up = is_boundary
                             ? t_bv[s * t_bl + (-neighbour - 1)]
                             : t_ev[s * 3 * n + neighbour * 3 + neighbour_edges[ki]];
                    } else {
                        c_up = t_ev[s * 3 * n + ki];
                    }
                    double tflux = wflux * c_up;
                    if (t_sf != NULL) tflux *= t_sf[s * n + donor];
                    t_eu[s * n + k] += tflux;
                    if (count_bdry) {
                        /* Same sign as edgeflux[0]: positive is inflow. */
                        t_bf[s * n + k] += tflux;
                    }
                }
            }
            // -------------------------------------------------------------

            // Boundary flux tracking: if this cell is not a ghost, and the neighbour
            // is a boundary condition OR a ghost cell, add the flux to boundary integral
            if (tri_full_flag != NULL) {
                int is_full = (tri_full_flag[k] == 1);
                int neighbour_is_ghost = (!is_boundary && tri_full_flag[neighbour] == 0);
                if ((is_boundary && is_full) || (is_full && neighbour_is_ghost)) {
                    boundary_flux_sum_substep += edgeflux[0];
                }
            }

            // Pressure gradient (gravity) terms
            double pressuregrad_work = length * (-g * 0.5 * (h_left * h_left - hle * hle
                                       - (hle + hc) * (zl - zc)) + pressure_flux);
            xmom_eu[k] -= normals[ki2] * pressuregrad_work;
            ymom_eu[k] -= normals[ki2 + 1] * pressuregrad_work;

        } // End edge loop

        // Update timestep only on first substep and for non-ghost cells
        if (substep_count == 0) {
            if (tri_full_flag == NULL || tri_full_flag[k] == 1) {
                if (speed_max_last > epsilon) {
                    double cell_timestep = radii[k] / speed_max_last;
                    local_timestep = fmin(local_timestep, cell_timestep);
                }
            }
            max_speed_array[k] = speed_max_last;
        }

        // Normalize by area
        double inv_area = 1.0 / areas[k];
        stage_eu[k] *= inv_area;
        xmom_eu[k] *= inv_area;
        ymom_eu[k] *= inv_area;
        if (n_tracers > 0) {
#ifdef CPU_ONLY_MODE
            double * restrict t_eu = D->tracer_explicit_update;
#endif
            for (anuga_int s = 0; s < n_tracers; s++) t_eu[s * n + k] *= inv_area;
        }

    } // End element loop

    // Store boundary flux sum for this substep
    if (D->boundary_flux_sum != NULL && substep_count < timestep_fluxcalls) {
        D->boundary_flux_sum[substep_count] = boundary_flux_sum_substep;
    }

    /* Same, per tracer. The main loop accumulated per cell so it needed no
     * reduction; total it here with ONE SCALAR REDUCTION PER TRACER, which is
     * a pattern every target supports, instead of the runtime-length array
     * reduction the main loop would have required. n_tracers is small, and the
     * whole block is skipped when there are none. The scratch is zeroed in the
     * same pass, ready for the next substep. */
    if (n_tracers > 0 && D->tracer_boundary_flux != NULL
            && D->tracer_boundary_flux_sum != NULL
            && substep_count < timestep_fluxcalls) {
#ifdef CPU_ONLY_MODE
        /* The GPU build has this hoisted to function scope (see the note on
         * tracer base pointers above); the CPU build loads it in-guard inside
         * the element loop, so it needs its own here. */
        double * restrict t_bf = D->tracer_boundary_flux;
#endif
        for (anuga_int s = 0; s < n_tracers; s++) {
            double tsum = 0.0;
            const anuga_int off = s * n;
#ifdef CPU_ONLY_MODE
            #pragma omp parallel for reduction(+:tsum)
#else
            #pragma omp target teams distribute parallel for reduction(+:tsum)
#endif
            for (anuga_int k = 0; k < n; k++) {
                tsum += t_bf[off + k];
                t_bf[off + k] = 0.0;
            }
            D->tracer_boundary_flux_sum[substep_count * n_tracers + s] = tsum;
        }
    }

    // Return timestep (only meaningful on first substep)
    return local_timestep;
}

// ============================================================================
// Edge-based flux computation (opt-in, two kernels)
//
// The cell-based kernel above solves every interior edge's Riemann problem
// TWICE -- once from each side, with swapped inputs and a flipped normal.
// The central-upwind flux is antisymmetric under that swap and its shared
// scalars (pressure_flux, max wave speed, z_half) are swap-invariant, so a
// single owner-side evaluation serves both cells: the same discretization,
// half the Riemann solves, and an EXACTLY antisymmetric flux exchange (the
// dual evaluation is only antisymmetric to floating-point roundoff).
//
// Kernel A (core_compute_fluxes_edge_based) runs one thread per cell-edge
// slot and computes only the slots it owns (boundary edges, or the side
// whose cell index is larger), storing per-slot
//     [F0, F1, F2, pf_len, z_half, speed]      (stride EDGE_SLOT_STRIDE)
// in D->edge_flux_work, where F* = -length * edgeflux (owner's sign) and
// pf_len = length * pressure_flux.  It also performs the min-dt and
// boundary-flux reductions the cell-based kernel does.
//
// Kernel B (core_flux_apply_and_update) is CELL-LOCAL: it gathers the three
// slot records (own sign for owned slots, negated for the neighbour's),
// assembles the one-sided pressure-gradient terms, normalizes by area, and
// -- because it is cell-local -- finishes the whole step in the same launch
// via gpu_cell_forcing_update (Manning + update + optional RK2 average).
// The explicit-update arrays are never written on this path: the values
// live and die in registers.
//
// Opt-in and restrictions: active only when the driver allocates
// D->edge_flux_work (EDGE_SLOT_STRIDE * 3n doubles; ANUGA leaves it NULL) --
// and, like reconstruct_edge_bed, it assumes fluxes follow an extrapolate,
// reconstructing bed values as stage - height.  Riverwalls are NOT
// supported (their weir corrections are one-sided), and neither are passive
// tracers (only the cell-based kernel advects them); callers must fall back
// to the cell-based kernel when riverwall edges or tracers exist
// (gpu_flux_mode in gpu_kernels.c does).
// ============================================================================

#define EDGE_SLOT_STRIDE 6

double core_compute_fluxes_edge_based(struct domain *D, int substep_count,
                                      int timestep_fluxcalls) {
    anuga_int n = D->number_of_elements;
    double g = D->g;
    double epsilon = D->epsilon;
    anuga_int low_froude = D->low_froude;

    double * restrict stage_ev = D->stage_edge_values;
    double * restrict xmom_ev = D->xmom_edge_values;
    double * restrict ymom_ev = D->ymom_edge_values;
    double * restrict height_ev = D->height_edge_values;

    double * restrict stage_bv = D->stage_boundary_values;
    double * restrict xmom_bv = D->xmom_boundary_values;
    double * restrict ymom_bv = D->ymom_boundary_values;

    anuga_int * restrict neighbours = D->neighbours;
    anuga_int * restrict neighbour_edges = D->neighbour_edges;
    anuga_geom_t * restrict normals = D->normals;
    anuga_geom_t * restrict edgelengths = D->edgelengths;
    anuga_geom_t * restrict radii = D->radii;
    anuga_int * restrict tri_full_flag = D->tri_full_flag;

    double * restrict slots = D->edge_flux_work;

    double local_timestep = 1.0e+100;
    double boundary_flux_sum_substep = 0.0;

    const anuga_int nslots = 3 * n;

    #ifdef CPU_ONLY_MODE
    #pragma omp parallel for reduction(min:local_timestep) reduction(+:boundary_flux_sum_substep)
    #else
    #pragma omp target teams distribute parallel for reduction(min:local_timestep) reduction(+:boundary_flux_sum_substep)
    #endif
    for (anuga_int p = 0; p < nslots; p++) {
        const anuga_int k = p / 3;
        const anuga_int nbr = neighbours[p];
        const int is_boundary = (nbr < 0);

        // Owner side only: boundary slots, or the side with the larger index
        if (!is_boundary && nbr < k) continue;

        double ql[3], qr[3], edgeflux[3];

        ql[0] = stage_ev[p];
        ql[1] = xmom_ev[p];
        ql[2] = ymom_ev[p];
        double hle = height_ev[p];
        double zl = ql[0] - hle;          // == bed_ev (post-extrapolate contract)

        double length = edgelengths[p];
        // Normals are read from the owner's slot; the neighbour's copy of the
        // same physical edge is the exact negation.
        double n1 = normals[2 * p];
        double n2 = normals[2 * p + 1];

        double zr, hre;
        if (is_boundary) {
            const anuga_int m = -nbr - 1;
            qr[0] = stage_bv[m];
            qr[1] = xmom_bv[m];
            qr[2] = ymom_bv[m];
            zr = zl;
            hre = fmax(qr[0] - zr, 0.0);
        } else {
            const anuga_int nm = 3 * nbr + neighbour_edges[p];
            qr[0] = stage_ev[nm];
            qr[1] = xmom_ev[nm];
            qr[2] = ymom_ev[nm];
            hre = height_ev[nm];
            zr = qr[0] - hre;
        }

        const double z_half = fmax(zl, zr);
        const double h_left = fmax(hle + zl - z_half, 0.0);
        const double h_right = fmax(hre + zr - z_half, 0.0);

        double max_speed_local = 0.0;
        double pressure_flux = 0.0;

        if (h_left == 0.0 && h_right == 0.0) {
            edgeflux[0] = 0.0;
            edgeflux[1] = 0.0;
            edgeflux[2] = 0.0;
        } else {
            gpu_flux_function_central(ql, qr, h_left, h_right, hle, hre,
                                      n1, n2, epsilon, z_half, g,
                                      edgeflux, &max_speed_local, &pressure_flux,
                                      low_froude);
        }

        const anuga_int base = EDGE_SLOT_STRIDE * p;
        slots[base + 0] = -length * edgeflux[0];
        slots[base + 1] = -length * edgeflux[1];
        slots[base + 2] = -length * edgeflux[2];
        slots[base + 3] = length * pressure_flux;
        slots[base + 4] = z_half;
        slots[base + 5] = max_speed_local;

        // Timestep reduction: min over (cell, edge) pairs of radii/speed --
        // identical to the cell-based min over cells of radii/max(speed),
        // since both evaluate radii/s at the cell's largest edge speed.
        if (substep_count == 0 && max_speed_local > epsilon) {
            if (tri_full_flag == NULL || tri_full_flag[k] == 1)
                local_timestep = fmin(local_timestep, radii[k] / max_speed_local);
            if (!is_boundary && (tri_full_flag == NULL || tri_full_flag[nbr] == 1))
                local_timestep = fmin(local_timestep, radii[nbr] / max_speed_local);
        }

        // Boundary flux integral: the full cell's own-side mass flux across
        // domain-boundary and full<->ghost edges (matches the cell-based sum)
        if (tri_full_flag != NULL) {
            const int k_full = (tri_full_flag[k] == 1);
            if (is_boundary) {
                if (k_full) boundary_flux_sum_substep += slots[base + 0];
            } else {
                const int n_full = (tri_full_flag[nbr] == 1);
                if (k_full && !n_full) boundary_flux_sum_substep += slots[base + 0];
                else if (!k_full && n_full) boundary_flux_sum_substep -= slots[base + 0];
            }
        }
    }

    if (D->boundary_flux_sum != NULL && substep_count < timestep_fluxcalls) {
        D->boundary_flux_sum[substep_count] = boundary_flux_sum_substep;
    }

    return local_timestep;
}

void core_flux_apply_and_update(struct domain *D, double timestep,
                                int apply_manning, int do_saxpy,
                                double a, double b, int substep_count) {
    anuga_int n = D->number_of_elements;
    double g = D->g;
    double minimum_allowed_height = D->minimum_allowed_height;
    double seven_thirds = ANUGA_SEVEN_THIRDS;

    double * restrict stage_cv = D->stage_centroid_values;
    double * restrict xmom_cv = D->xmom_centroid_values;
    double * restrict ymom_cv = D->ymom_centroid_values;
    double * restrict bed_cv = D->bed_centroid_values;
    double * restrict height_cv = D->height_centroid_values;
    double * restrict friction_cv = D->friction_centroid_values;

    double * restrict stage_ev = D->stage_edge_values;
    double * restrict height_ev = D->height_edge_values;

    double * restrict stage_siu = D->stage_semi_implicit_update;
    double * restrict xmom_siu = D->xmom_semi_implicit_update;
    double * restrict ymom_siu = D->ymom_semi_implicit_update;

    double * restrict stage_bk = D->stage_backup_values;
    double * restrict xmom_bk = D->xmom_backup_values;
    double * restrict ymom_bk = D->ymom_backup_values;

    anuga_int * restrict neighbours = D->neighbours;
    anuga_int * restrict neighbour_edges = D->neighbour_edges;
    anuga_geom_t * restrict normals = D->normals;
    anuga_geom_t * restrict edgelengths = D->edgelengths;
    anuga_geom_t * restrict areas = D->areas;
    double * restrict max_speed_array = D->max_speed;

    double * restrict slots = D->edge_flux_work;

    OMP_PARALLEL_LOOP
    for (anuga_int k = 0; k < n; k++) {
        double eu_stage = 0.0, eu_xmom = 0.0, eu_ymom = 0.0;
        double speed_max_last = 0.0;

        const double hc = height_cv[k];
        const double zc = bed_cv[k];

        for (int i = 0; i < 3; i++) {
            const anuga_int p = 3 * k + i;
            const anuga_int nbr = neighbours[p];
            const int owner = (nbr < 0 || nbr > k);
            const anuga_int slot = owner ? p : 3 * nbr + neighbour_edges[p];
            const anuga_int base = EDGE_SLOT_STRIDE * slot;
            const double sgn = owner ? 1.0 : -1.0;

            eu_stage += sgn * slots[base + 0];
            eu_xmom  += sgn * slots[base + 1];
            eu_ymom  += sgn * slots[base + 2];

            const double pf_len = slots[base + 3];
            const double z_half = slots[base + 4];
            speed_max_last = fmax(speed_max_last, slots[base + 5]);

            // One-sided pressure-gradient term, from this cell's own edge
            // values -- the same expression the cell-based kernel evaluates.
            // (pf_len uses the owner's edge length; both sides of a physical
            // edge share endpoints, so the lengths are bit-identical.)
            const double hle = height_ev[p];
            const double zl = stage_ev[p] - hle;
            const double length = edgelengths[p];
            const double h_side = fmax(hle + zl - z_half, 0.0);

            const double pg = pf_len
                - length * g * 0.5 * (h_side * h_side - hle * hle
                                      - (hle + hc) * (zl - zc));
            eu_xmom -= normals[2 * p] * pg;
            eu_ymom -= normals[2 * p + 1] * pg;
        }

        const double inv_area = 1.0 / areas[k];
        eu_stage *= inv_area;
        eu_xmom  *= inv_area;
        eu_ymom  *= inv_area;

        if (substep_count == 0) max_speed_array[k] = speed_max_last;

        gpu_cell_forcing_update(k, timestep, apply_manning, do_saxpy, a, b,
                                g, minimum_allowed_height, seven_thirds,
                                eu_stage, eu_xmom, eu_ymom,
                                stage_cv, xmom_cv, ymom_cv, bed_cv, height_cv,
                                friction_cv, stage_siu, xmom_siu, ymom_siu,
                                stage_bk, xmom_bk, ymom_bk);
    }
}

// ============================================================================
// Active-set construction (opt-in, for wet/dry flood domains)
//
// A dry cell whose three neighbours are also dry cannot change: its fluxes
// are zero and every forcing term is gated on depth.  On flood domains that
// is most of the mesh most of the time (the 11,000 km^2 spec basin starts
// 99.5% dry), so the step kernels accept an optional iteration list and the
// driver rebuilds, each step:
//   active cells = wet cells  U  neighbours of wet cells  U  boundary cells
//   active edges = owned edges with either side active
// Skipping the rest is BIT-EXACT under two provisos the driver must honour:
// the first step runs full (so every cell's edge values and protect clamps
// are populated once), and boundary-adjacent cells stay in the set (their
// edge values feed the boundary-value kernels; also open boundaries may wet
// them from outside).  A skipped cell's stale edge values are only ever read
// by fluxes on edges whose BOTH sides are inactive -- and those edges have
// zero height on both sides, which the flux kernel already short-circuits.
//
// counts_out[0] = active cells, counts_out[1] = active edges.
// ============================================================================

// Wetness threshold for classification: above pure-roundoff films, far
// below physics.  A strict > 0 test suffers roundoff creep -- update-sum
// cancellation deposits ~1e-15 m "films" on shoreline-adjacent dry cells,
// which then activate their neighbours ring by ring until the whole mesh is
// active.  Films below 1e-12 m produce fluxes ~1e-21 relative -- beneath
// double precision of the stored state -- so skipping them remains
// bit-exact (verified against full-run goldens).
#define ACTIVE_WET_EPS 1.0e-12

void core_build_active_sets(struct domain *D,
                            anuga_int * restrict wet_flag,
                            anuga_int * restrict ring1_flag,
                            anuga_int * restrict active_cells,
                            anuga_int * restrict active_edges,
                            const anuga_int * restrict owned_edges,
                            anuga_int num_owned_edges,
                            anuga_int *counts_out) {
    anuga_int n = D->number_of_elements;
    double * restrict stage_cv = D->stage_centroid_values;
    double * restrict bed_cv = D->bed_centroid_values;
    anuga_int * restrict neighbours = D->neighbours;

    // Pass 1: wetness (no gathers).  Classified from stage - bed rather than
    // height_cv: operators (rain, inflows, culverts) modify STAGE directly,
    // and height_cv is only refreshed by prepare -- which visits listed cells
    // only.  Classifying on height would let rain accumulate invisibly on
    // inactive cells, never flowing; stage - bed is always current.
    OMP_PARALLEL_LOOP
    for (anuga_int k = 0; k < n; k++) {
        wet_flag[k] = (stage_cv[k] - bed_cv[k] > ACTIVE_WET_EPS) ? 1 : 0;
    }

    // Pass 2: ring-1 = wet, neighbour-of-wet, or boundary-adjacent.
    // Ring-1 cells are the ones whose edges can carry flux this step.
    OMP_PARALLEL_LOOP
    for (anuga_int k = 0; k < n; k++) {
        int act = wet_flag[k];
        for (int i = 0; i < 3 && !act; i++) {
            const anuga_int nbr = neighbours[3 * k + i];
            if (nbr < 0 || wet_flag[nbr]) act = 1;
        }
        ring1_flag[k] = act;
    }

    // Pass 3: the CELL list is ring-2 (ring-1 plus its neighbours).  Water
    // advances at most one ring per flux call, and a rebuild covers a whole
    // step (two flux calls under RK2), so the update must reach one ring
    // beyond the cells the fluxes can wet -- exactly the MPI ghost-layer
    // width rule.  With a 1-ring cell list, substep 2 could scatter flux
    // into a cell the update never visits: silent mass loss.
    anuga_int n_cells = 0;
    #ifdef CPU_ONLY_MODE
    #pragma omp parallel for
    #else
    #pragma omp target teams distribute parallel for map(tofrom: n_cells)
    #endif
    for (anuga_int k = 0; k < n; k++) {
        int act = ring1_flag[k];
        for (int i = 0; i < 3 && !act; i++) {
            const anuga_int nbr = neighbours[3 * k + i];
            if (nbr >= 0 && ring1_flag[nbr]) act = 1;
        }
        if (act) {
            anuga_int idx;
            #pragma omp atomic capture
            idx = n_cells++;
            active_cells[idx] = k;
        }
    }

    // Pass 4: active owned edges = either side in ring-1
    anuga_int n_edges = 0;
    #ifdef CPU_ONLY_MODE
    #pragma omp parallel for
    #else
    #pragma omp target teams distribute parallel for map(tofrom: n_edges)
    #endif
    for (anuga_int q = 0; q < num_owned_edges; q++) {
        const anuga_int p = owned_edges[q];
        const anuga_int k = p / 3;
        const anuga_int nbr = neighbours[p];
        if (ring1_flag[k] || (nbr >= 0 && ring1_flag[nbr])) {
            anuga_int idx;
            #pragma omp atomic capture
            idx = n_edges++;
            active_edges[idx] = p;
        }
    }

    counts_out[0] = n_cells;
    counts_out[1] = n_edges;
}

// ============================================================================
// Scatter-mode flux computation (opt-in, single kernel + atomics)
//
// Same single-Riemann-solve-per-edge idea as the slot-based pair above, but
// with NO intermediate storage: the owner thread computes the flux once and
// scatters both sides' full contributions -- flux exchange AND each side's
// one-sided pressure-gradient term, already area-normalized -- directly into
// the explicit-update arrays with `omp atomic update` (portable OpenMP; each
// eu entry receives at most 3 concurrent adds, so contention is negligible).
// The slot-based variant measured SLOWER than the cell-based kernel because
// the 144 B/cell of slot records cost more to move than the duplicate
// Riemann solves saved; this variant keeps the saved solves and moves
// nothing.
//
// Requirements (same contract as the slot variant): the explicit updates
// must be ZERO on entry (core_prepare_step's zero_eu flag), no riverwalls,
// no tracers, and fluxes follow an extrapolate (bed reconstructed as
// stage - height).
// max_speed_array is NOT maintained on this path (per-cell max would need an
// atomic max); the wave speeds live and die in registers, so this mode needs
// NO auxiliary arrays at all -- drivers select it with
// D->reconstruct_edge_bed = 2 and pay zero extra device memory.
// ============================================================================

double core_compute_fluxes_scatter_on(struct domain *D, int substep_count,
                                      int timestep_fluxcalls,
                                      const anuga_int * restrict edges,
                                      anuga_int nedges) {
    anuga_int n = D->number_of_elements;
    double g = D->g;
    double epsilon = D->epsilon;
    anuga_int low_froude = D->low_froude;

    double * restrict stage_ev = D->stage_edge_values;
    double * restrict xmom_ev = D->xmom_edge_values;
    double * restrict ymom_ev = D->ymom_edge_values;
    double * restrict height_ev = D->height_edge_values;

    double * restrict stage_bv = D->stage_boundary_values;
    double * restrict xmom_bv = D->xmom_boundary_values;
    double * restrict ymom_bv = D->ymom_boundary_values;

    double * restrict stage_eu = D->stage_explicit_update;
    double * restrict xmom_eu = D->xmom_explicit_update;
    double * restrict ymom_eu = D->ymom_explicit_update;

    double * restrict height_cv = D->height_centroid_values;
    double * restrict bed_cv = D->bed_centroid_values;

    anuga_int * restrict neighbours = D->neighbours;
    anuga_int * restrict neighbour_edges = D->neighbour_edges;
    anuga_geom_t * restrict normals = D->normals;
    anuga_geom_t * restrict edgelengths = D->edgelengths;
    anuga_geom_t * restrict radii = D->radii;
    anuga_geom_t * restrict areas = D->areas;
    anuga_int * restrict tri_full_flag = D->tri_full_flag;

    double local_timestep = 1.0e+100;
    double boundary_flux_sum_substep = 0.0;

    #ifdef CPU_ONLY_MODE
    #pragma omp parallel for reduction(min:local_timestep) reduction(+:boundary_flux_sum_substep)
    #else
    #pragma omp target teams distribute parallel for reduction(min:local_timestep) reduction(+:boundary_flux_sum_substep)
    #endif
    for (anuga_int q = 0; q < nedges; q++) {
        const anuga_int p = edges[q];
        const anuga_int k = p / 3;
        const anuga_int nbr = neighbours[p];
        const int is_boundary = (nbr < 0);

        double ql[3], qr[3], edgeflux[3];

        ql[0] = stage_ev[p];
        ql[1] = xmom_ev[p];
        ql[2] = ymom_ev[p];
        double hle = height_ev[p];
        double zl = ql[0] - hle;

        double length = edgelengths[p];
        double n1 = normals[2 * p];
        double n2 = normals[2 * p + 1];

        double zr, hre;
        anuga_int nm = 0;
        if (is_boundary) {
            const anuga_int m = -nbr - 1;
            qr[0] = stage_bv[m];
            qr[1] = xmom_bv[m];
            qr[2] = ymom_bv[m];
            zr = zl;
            hre = fmax(qr[0] - zr, 0.0);
        } else {
            nm = 3 * nbr + neighbour_edges[p];
            qr[0] = stage_ev[nm];
            qr[1] = xmom_ev[nm];
            qr[2] = ymom_ev[nm];
            hre = height_ev[nm];
            zr = qr[0] - hre;
        }

        const double z_half = fmax(zl, zr);
        const double h_left = fmax(hle + zl - z_half, 0.0);
        const double h_right = fmax(hre + zr - z_half, 0.0);

        double max_speed_local = 0.0;
        double pressure_flux = 0.0;

        if (h_left == 0.0 && h_right == 0.0) {
            edgeflux[0] = 0.0;
            edgeflux[1] = 0.0;
            edgeflux[2] = 0.0;
        } else {
            gpu_flux_function_central(ql, qr, h_left, h_right, hle, hre,
                                      n1, n2, epsilon, z_half, g,
                                      edgeflux, &max_speed_local, &pressure_flux,
                                      low_froude);
        }

        const double F0 = -length * edgeflux[0];
        const double F1 = -length * edgeflux[1];
        const double F2 = -length * edgeflux[2];
        const double pf_len = length * pressure_flux;

        // ---- owner side: flux + its one-sided pressure gradient, scaled
        {
            const double hc = height_cv[k];
            const double zc = bed_cv[k];
            const double pg = pf_len
                - length * g * 0.5 * (h_left * h_left - hle * hle
                                      - (hle + hc) * (zl - zc));
            const double ia = 1.0 / areas[k];
            const double ds = F0 * ia;
            const double dx = (F1 - n1 * pg) * ia;
            const double dy = (F2 - n2 * pg) * ia;

            #pragma omp atomic update
            stage_eu[k] += ds;
            #pragma omp atomic update
            xmom_eu[k] += dx;
            #pragma omp atomic update
            ymom_eu[k] += dy;
        }

        // ---- neighbour side: negated flux, its own pressure gradient.
        // The neighbour's stored normal for this physical edge is the exact
        // FP negation of the owner's (same endpoints, opposite subtraction),
        // so (-n1, -n2) reproduces its cell-based expression bit-for-bit.
        if (!is_boundary) {
            const double hc_n = height_cv[nbr];
            const double zc_n = bed_cv[nbr];
            const double pg = pf_len
                - length * g * 0.5 * (h_right * h_right - hre * hre
                                      - (hre + hc_n) * (zr - zc_n));
            const double ia = 1.0 / areas[nbr];
            const double ds = -F0 * ia;
            const double dx = (-F1 - (-n1) * pg) * ia;
            const double dy = (-F2 - (-n2) * pg) * ia;

            #pragma omp atomic update
            stage_eu[nbr] += ds;
            #pragma omp atomic update
            xmom_eu[nbr] += dx;
            #pragma omp atomic update
            ymom_eu[nbr] += dy;
        }

        if (substep_count == 0 && max_speed_local > epsilon) {
            if (tri_full_flag == NULL || tri_full_flag[k] == 1)
                local_timestep = fmin(local_timestep, radii[k] / max_speed_local);
            if (!is_boundary && (tri_full_flag == NULL || tri_full_flag[nbr] == 1))
                local_timestep = fmin(local_timestep, radii[nbr] / max_speed_local);
        }

        if (tri_full_flag != NULL) {
            const int k_full = (tri_full_flag[k] == 1);
            if (is_boundary) {
                if (k_full) boundary_flux_sum_substep += F0;
            } else {
                const int n_full = (tri_full_flag[nbr] == 1);
                if (k_full && !n_full) boundary_flux_sum_substep += F0;
                else if (!k_full && n_full) boundary_flux_sum_substep -= F0;
            }
        }
    }

    if (D->boundary_flux_sum != NULL && substep_count < timestep_fluxcalls) {
        D->boundary_flux_sum[substep_count] = boundary_flux_sum_substep;
    }

    return local_timestep;
}

double core_compute_fluxes_scatter(struct domain *D, int substep_count,
                                   int timestep_fluxcalls) {
    // One thread per PHYSICAL edge via the driver-built compacted slot list
    return core_compute_fluxes_scatter_on(D, substep_count, timestep_fluxcalls,
                                          D->owned_edges, D->num_owned_edges);
}

// ============================================================================
// ADER Cauchy-Kovalewski predictor
// ============================================================================

void core_ader_ck_predictor(struct domain *D, double dt) {
    // Advance centroid values by dt using a local Cauchy-Kovalewski predictor.
    // Called after core_extrapolate_second_order_edge() so edge_values hold the
    // reconstructed state.  Slopes are recovered from the 2x2 linear system
    // formed by edges 0 and 1 (no new arrays needed).
    //
    // Well-balanced form: bed slope is dz/dx = dw/dx - dh/dx derived from
    // the reconstruction, so still-water equilibrium is preserved exactly.

    anuga_int n = D->number_of_elements;
    double g   = D->g;
    double eps = D->minimum_allowed_height;

    double * restrict stage_cv  = D->stage_centroid_values;
    double * restrict xmom_cv   = D->xmom_centroid_values;
    double * restrict ymom_cv   = D->ymom_centroid_values;
    double * restrict bed_cv    = D->bed_centroid_values;
    double * restrict height_cv = D->height_centroid_values;

    double * restrict stage_ev  = D->stage_edge_values;
    double * restrict xmom_ev   = D->xmom_edge_values;
    double * restrict ymom_ev   = D->ymom_edge_values;
    double * restrict height_ev = D->height_edge_values;

    anuga_geom_t * restrict edge_coords     = D->edge_coordinates;
    anuga_geom_t * restrict centroid_coords = D->centroid_coordinates;

    OMP_PARALLEL_LOOP
    for (anuga_int k = 0; k < n; k++) {
        anuga_int k3 = k * 3;
        anuga_int k6 = k * 6;
        anuga_int k2 = k * 2;

        // Offsets from centroid to edge midpoints 0 and 1
        double xc   = centroid_coords[k2 + 0];
        double yc   = centroid_coords[k2 + 1];
        double dxv0 = edge_coords[k6 + 0] - xc;
        double dyv0 = edge_coords[k6 + 1] - yc;
        double dxv1 = edge_coords[k6 + 2] - xc;
        double dyv1 = edge_coords[k6 + 3] - yc;

        // Determinant of the 2x2 linear system; skip degenerate cells.
        // Use if-block (not continue) for GPU target-loop compatibility.
        double det = dxv0 * dyv1 - dxv1 * dyv0;
        if (fabs(det) >= 1.0e-20) {
        double inv_det = 1.0 / det;

        // Centroid state
        double w_c  = stage_cv[k];
        double h_c  = fmax(w_c - bed_cv[k], 0.0);
        double uh_c = xmom_cv[k];
        double vh_c = ymom_cv[k];

        // Centroid velocity (guarded)
        double inv_h_c = (h_c > eps) ? 1.0 / h_c : 0.0;
        double u_c = uh_c * inv_h_c;
        double v_c = vh_c * inv_h_c;

        // Recover gradients from edge differences using edges 0 and 1.
        // For any variable q:  grad_x = inv_det*(dyv1*dq0 - dyv0*dq1)
        //                      grad_y = inv_det*(dxv0*dq1 - dxv1*dq0)

        // Stage gradient (∂w/∂x, ∂w/∂y)
        double dw0 = stage_ev[k3 + 0] - w_c;
        double dw1 = stage_ev[k3 + 1] - w_c;
        double wx  = inv_det * (dyv1 * dw0 - dyv0 * dw1);
        double wy  = inv_det * (dxv0 * dw1 - dxv1 * dw0);

        // Height gradient (∂h/∂x, ∂h/∂y)
        double dh0 = height_ev[k3 + 0] - h_c;
        double dh1 = height_ev[k3 + 1] - h_c;
        double hx  = inv_det * (dyv1 * dh0 - dyv0 * dh1);
        double hy  = inv_det * (dxv0 * dh1 - dxv1 * dh0);

        // Edge velocities (recover from edge momentum / edge height)
        double h_e0     = height_ev[k3 + 0];
        double h_e1     = height_ev[k3 + 1];
        double inv_h_e0 = (h_e0 > eps) ? 1.0 / h_e0 : 0.0;
        double inv_h_e1 = (h_e1 > eps) ? 1.0 / h_e1 : 0.0;
        double u_e0 = xmom_ev[k3 + 0] * inv_h_e0;
        double u_e1 = xmom_ev[k3 + 1] * inv_h_e1;
        double v_e0 = ymom_ev[k3 + 0] * inv_h_e0;
        double v_e1 = ymom_ev[k3 + 1] * inv_h_e1;

        // Velocity gradients (∂u/∂x, ∂u/∂y, ∂v/∂x, ∂v/∂y)
        double du0 = u_e0 - u_c;
        double du1 = u_e1 - u_c;
        double dv0 = v_e0 - v_c;
        double dv1 = v_e1 - v_c;
        double ux  = inv_det * (dyv1 * du0 - dyv0 * du1);
        double uy  = inv_det * (dxv0 * du1 - dxv1 * du0);
        double vx  = inv_det * (dyv1 * dv0 - dyv0 * dv1);
        double vy  = inv_det * (dxv0 * dv1 - dxv1 * dv0);

        // Cauchy-Kovalewski time derivatives — well-balanced SWE:
        //   dz/dx = dw/dx - dh/dx  (from reconstruction, not stored centroid z)
        // This ensures cancellation in still water (u=v=0, wx=wy=0).
        double g_h = g * h_c;

        double dw_dt  = -(u_c * hx + h_c * ux + v_c * hy + h_c * vy);
        double duh_dt = -(2.0*u_c*h_c*ux + u_c*u_c*hx + u_c*v_c*hy
                         + v_c*h_c*uy + u_c*h_c*vy + g_h * wx);
        double dvh_dt = -(v_c*h_c*ux + u_c*h_c*vx + u_c*v_c*hx
                         + 2.0*v_c*h_c*vy + v_c*v_c*hy + g_h * wy);

        // Predict forward by dt (caller passes dt/2 for midpoint)
        double w_pred  = w_c  + dt * dw_dt;
        double uh_pred = uh_c + dt * duh_dt;
        double vh_pred = vh_c + dt * dvh_dt;
        double h_pred  = fmax(w_pred - bed_cv[k], 0.0);

        stage_cv[k]  = w_pred;
        xmom_cv[k]   = uh_pred;
        ymom_cv[k]   = vh_pred;
        height_cv[k] = h_pred;
        } // end if (fabs(det) >= 1.0e-20)
    }
}

void core_ader_ck_predictor_edge(struct domain *D, double dt) {
    // Fused ADER-2 predictor: advances edge values to Q^{n+1/2} in-place,
    // leaving centroid values unchanged.  This eliminates the second full
    // extrapolation pass needed by core_ader_ck_predictor (centroid variant).
    //
    // For any quantity q, the reconstructed edge value is:
    //   q_edge[i] = q_c + slope * offset_i
    // Since the predictor adds the same centroid shift dq_c to every edge,
    //   q_edge_pred[i] = q_edge[i] + dt * dq_c/dt
    // The cell slopes are preserved exactly.
    //
    // Well-balanced: same dz/dx = dw/dx - dh/dx derivation as the centroid
    // variant; still-water equilibrium is preserved exactly.

    anuga_int n = D->number_of_elements;
    double g   = D->g;
    double eps = D->minimum_allowed_height;

    double * restrict stage_cv  = D->stage_centroid_values;
    double * restrict xmom_cv   = D->xmom_centroid_values;
    double * restrict ymom_cv   = D->ymom_centroid_values;
    double * restrict bed_cv    = D->bed_centroid_values;

    double * restrict stage_ev  = D->stage_edge_values;
    double * restrict xmom_ev   = D->xmom_edge_values;
    double * restrict ymom_ev   = D->ymom_edge_values;
    double * restrict height_ev = D->height_edge_values;
    double * restrict friction_cv = D->friction_centroid_values;   // predictor drag
    const int sloped_mannings = (int)D->use_sloped_mannings;

    anuga_geom_t * restrict edge_coords     = D->edge_coordinates;
    anuga_geom_t * restrict centroid_coords = D->centroid_coordinates;

    const double seven_thirds = ANUGA_SEVEN_THIRDS;

    OMP_PARALLEL_LOOP
    for (anuga_int k = 0; k < n; k++) {
        anuga_int k3 = k * 3;
        anuga_int k6 = k * 6;
        anuga_int k2 = k * 2;

        double xc   = centroid_coords[k2 + 0];
        double yc   = centroid_coords[k2 + 1];
        double dxv0 = edge_coords[k6 + 0] - xc;
        double dyv0 = edge_coords[k6 + 1] - yc;
        double dxv1 = edge_coords[k6 + 2] - xc;
        double dyv1 = edge_coords[k6 + 3] - yc;

        double det = dxv0 * dyv1 - dxv1 * dyv0;
        if (fabs(det) >= 1.0e-20) {
        double inv_det = 1.0 / det;

        double w_c  = stage_cv[k];
        double h_c  = fmax(w_c - bed_cv[k], 0.0);
        double uh_c = xmom_cv[k];
        double vh_c = ymom_cv[k];

        double inv_h_c = (h_c > eps) ? 1.0 / h_c : 0.0;
        double u_c = uh_c * inv_h_c;
        double v_c = vh_c * inv_h_c;

        double dw0 = stage_ev[k3 + 0] - w_c;
        double dw1 = stage_ev[k3 + 1] - w_c;
        double wx  = inv_det * (dyv1 * dw0 - dyv0 * dw1);
        double wy  = inv_det * (dxv0 * dw1 - dxv1 * dw0);

        double dh0 = height_ev[k3 + 0] - h_c;
        double dh1 = height_ev[k3 + 1] - h_c;
        double hx  = inv_det * (dyv1 * dh0 - dyv0 * dh1);
        double hy  = inv_det * (dxv0 * dh1 - dxv1 * dh0);

        double h_e0     = height_ev[k3 + 0];
        double h_e1     = height_ev[k3 + 1];
        double inv_h_e0 = (h_e0 > eps) ? 1.0 / h_e0 : 0.0;
        double inv_h_e1 = (h_e1 > eps) ? 1.0 / h_e1 : 0.0;
        double u_e0 = xmom_ev[k3 + 0] * inv_h_e0;
        double u_e1 = xmom_ev[k3 + 1] * inv_h_e1;
        double v_e0 = ymom_ev[k3 + 0] * inv_h_e0;
        double v_e1 = ymom_ev[k3 + 1] * inv_h_e1;

        double du0 = u_e0 - u_c;
        double du1 = u_e1 - u_c;
        double dv0 = v_e0 - v_c;
        double dv1 = v_e1 - v_c;
        double ux  = inv_det * (dyv1 * du0 - dyv0 * du1);
        double uy  = inv_det * (dxv0 * du1 - dxv1 * du0);
        double vx  = inv_det * (dyv1 * dv0 - dyv0 * dv1);
        double vy  = inv_det * (dxv0 * dv1 - dxv1 * dv0);

        double g_h = g * h_c;
        double dw_dt  = -(u_c * hx + h_c * ux + v_c * hy + h_c * vy);
        double duh_dt = -(2.0*u_c*h_c*ux + u_c*u_c*hx + u_c*v_c*hy
                         + v_c*h_c*uy + u_c*h_c*vy + g_h * wx);
        double dvh_dt = -(v_c*h_c*ux + u_c*h_c*vx + u_c*v_c*hx
                         + 2.0*v_c*h_c*vy + v_c*v_c*hy + g_h * wy);

        // Manning drag on the half-advanced state.  The derivative terms
        // above carry the flux Jacobian and the gravity term g*h*grad(w) but
        // no drag, so on a slope the predictor is pulled downhill with
        // nothing resisting it: in thin flow, where gravity against friction
        // IS the balance, the fluxes are then taken from an over-driven state
        // that the corrector's friction cannot claw back.  On the coarse mild
        // slope that cost 55% of the discharge.  See
        // core_ader_friction_factor() for why this is applied as a factor
        // rather than folded into duh_dt.
        double zs_p = 1.0;
        if (sloped_mannings) {
            double zx_p = wx - hx;                   // w = z + h
            double zy_p = wy - hy;
            zs_p = sqrt(1.0 + zx_p * zx_p + zy_p * zy_p);
        }
        double eta_p = friction_cv[k];

        // Shift all three edges by the same centroid delta (slopes preserved)
        for (int i = 0; i < 3; i++) {
            // gamma from the PRE-shift state; see core_ader_friction_factor().
            double fac = core_ader_friction_factor(g, eta_p, zs_p, eps,
                                                   seven_thirds,
                                                   height_ev[k3 + i],
                                                   xmom_ev[k3 + i],
                                                   ymom_ev[k3 + i], dt);

            stage_ev[k3 + i] += dt * dw_dt;
            double uh_e = xmom_ev[k3 + i] + dt * duh_dt;
            double vh_e = ymom_ev[k3 + i] + dt * dvh_dt;
            height_ev[k3 + i] = fmax(height_ev[k3 + i] + dt * dw_dt, 0.0);

            xmom_ev[k3 + i] = uh_e * fac;
            ymom_ev[k3 + i] = vh_e * fac;
        }
        } // end if (fabs(det) >= 1.0e-20)
    }
}
