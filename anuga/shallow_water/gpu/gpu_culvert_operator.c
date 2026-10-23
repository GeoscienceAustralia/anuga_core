// GPU-accelerated culvert (Boyd box/pipe) operator
//
// Strategy: Batch all culverts into ONE gather → CPU compute → ONE scatter
// per timestep, reducing GPU↔CPU sync from 2×N_culverts to exactly 2.
//
// The Boyd discharge physics runs on CPU (200-300 FLOPs of branchy serial code).
// Only the data movement (gather/scatter ~2KB) touches the GPU.

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <omp.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif
#include "gpu_domain.h"
#include "gpu_culvert_operator.h"
#include "gpu_omp_macros.h"

// Host-side per-step scratch, sized to the culvert count.  Defined further
// down, next to struct culvert_mpi_bufs, but freed by gpu_culverts_finalize_all
// above it.
struct culvert_operators;
static void culvert_host_scratch_free(struct culvert_operators *CO);
#include "gpu_nvtx.h"

#define VELOCITY_PROTECTION 1.0e-6

// ============================================================================
// Pure computation: Boyd BOX discharge
// Direct translation from boyd_box_operator.py:boyd_box_function()
// ============================================================================

void boyd_box_discharge(const struct culvert_params *p,
                        double driving_energy,
                        double delta_total_energy,
                        double outlet_enquiry_depth,
                        double *Q_out, double *barrel_velocity_out,
                        double *outlet_culvert_depth_out, double *flow_area_out) {

    double width = p->width;
    double depth = p->height;
    double blockage = p->blockage;
    double barrels = p->barrels;
    double length = p->length;
    double sum_loss = p->sum_loss;
    double manning = p->manning;

    double bf = 1.0 - blockage;

    if (blockage >= 1.0) {
        *Q_out = 0.0;
        *barrel_velocity_out = 0.0;
        *outlet_culvert_depth_out = 0.0;
        *flow_area_out = 0.00001;
        return;
    }

    // Inlet control: unsubmerged vs submerged
    double Q_inlet_unsubmerged = 0.544 * sqrt(p->g) * bf * width * barrels * pow(driving_energy, 1.50);
    double Q_inlet_submerged = 0.702 * sqrt(p->g) * bf * width * barrels * pow(depth, 0.89) * pow(driving_energy, 0.61);

    double Q, dcrit, flow_area, perimeter;

    if (Q_inlet_unsubmerged < Q_inlet_submerged) {
        Q = Q_inlet_unsubmerged;
    } else {
        Q = Q_inlet_submerged;
    }

    dcrit = pow(Q * Q / p->g / pow(bf * width * barrels, 2.0), 0.333333);

    if (dcrit > depth) {
        dcrit = depth;
        flow_area = bf * width * dcrit * barrels;
        perimeter = barrels * 2.0 * (bf * width + dcrit);
    } else {
        flow_area = bf * width * barrels * dcrit;
        perimeter = barrels * (bf * width + 2.0 * dcrit);
    }

    double outlet_culvert_depth = dcrit;

    // Recompute dcrit (matches Python exactly)
    dcrit = pow(Q * Q / p->g / pow(bf * width * barrels, 2.0), 0.333333);
    outlet_culvert_depth = dcrit;

    if (outlet_culvert_depth > depth) {
        outlet_culvert_depth = depth;
        flow_area = bf * width * barrels * depth;
        perimeter = barrels * 2.0 * (bf * width + depth);
    } else {
        flow_area = bf * width * barrels * outlet_culvert_depth;
        perimeter = barrels * (bf * width + 2.0 * outlet_culvert_depth);
    }

    double hyd_rad = flow_area / perimeter;
    double culvert_velocity = sqrt(delta_total_energy / ((sum_loss / 2.0 / p->g) +
                                    (manning * manning * length) / pow(hyd_rad, 1.33333)));
    double Q_outlet_tailwater = flow_area * culvert_velocity;

    // Outlet control check
    if (delta_total_energy < driving_energy) {
        if (outlet_enquiry_depth > depth) {
            // Outlet submerged
            outlet_culvert_depth = depth;
            flow_area = bf * width * barrels * depth;
            perimeter = barrels * 2.0 * (bf * width + depth);
        } else {
            dcrit = pow(Q * Q / p->g / pow(bf * width * barrels, 2.0), 0.333333);
            outlet_culvert_depth = dcrit;
            if (outlet_culvert_depth > depth) {
                outlet_culvert_depth = depth;
                flow_area = bf * width * barrels * depth;
                perimeter = barrels * 2.0 * (bf * width + depth);
            } else {
                flow_area = bf * width * barrels * outlet_culvert_depth;
                perimeter = barrels * (bf * width + 2.0 * outlet_culvert_depth);
            }
        }

        hyd_rad = flow_area / perimeter;
        culvert_velocity = sqrt(delta_total_energy / ((sum_loss / 2.0 / p->g) +
                                (manning * manning * length) / pow(hyd_rad, 1.33333)));
        Q_outlet_tailwater = flow_area * culvert_velocity;

        if (Q_outlet_tailwater < Q)
            Q = Q_outlet_tailwater;
    } else {
        /* Head difference at or above the driving energy: the barrel's
         * friction still limits the flow, as for pipes (see the Python). */
        if (Q_outlet_tailwater < Q)
            Q = Q_outlet_tailwater;
    }

    // Barrel velocity with protection
    double barrel_velocity = Q / (flow_area + VELOCITY_PROTECTION / flow_area);

    *Q_out = Q;
    *barrel_velocity_out = barrel_velocity;
    *outlet_culvert_depth_out = outlet_culvert_depth;
    *flow_area_out = flow_area;
}

// ============================================================================
// Critical depth in a circular conduit (per barrel).
// Same algorithm as boyd_pipe_operator.py:circular_critical_depth(): solve
// Q^2/g = A^3/T exactly by fixed-iteration bisection on the central angle
// theta, with A = D^2/8 (theta - sin theta), T = D sin(theta/2). A^3/T is
// monotone from 0 (dry) to infinity (full), so the root is unique and < D.
// ============================================================================

static double circular_critical_depth(double Q, double diameter, double g) {
    if (Q <= 0.0 || diameter <= 0.0) return 0.0;
    double target = Q * Q / g;
    double lo = 0.0, hi = 2.0 * M_PI;
    for (int it = 0; it < 60; it++) {
        double theta = 0.5 * (lo + hi);
        double area = diameter * diameter / 8.0 * (theta - sin(theta));
        double top_width = diameter * sin(0.5 * theta);
        if (area * area * area > target * top_width)
            hi = theta;
        else
            lo = theta;
    }
    double theta = 0.5 * (lo + hi);
    return 0.5 * diameter * (1.0 - cos(0.5 * theta));
}

// ============================================================================
// Pure computation: Boyd PIPE discharge
// Direct translation from boyd_pipe_operator.py:boyd_pipe_function()
// ============================================================================

void boyd_pipe_discharge(const struct culvert_params *p,
                         double inflow_depth,
                         double driving_energy,
                         double delta_total_energy,
                         double outlet_enquiry_depth,
                         double *Q_out, double *barrel_velocity_out,
                         double *outlet_culvert_depth_out, double *flow_area_out) {

    double diameter = p->diameter;
    double blockage = p->blockage;
    double barrels = p->barrels;
    double length = p->length;
    double sum_loss = p->sum_loss;
    double manning = p->manning;

    if (blockage >= 1.0) {
        *Q_out = 0.0;
        *barrel_velocity_out = 0.0;
        *outlet_culvert_depth_out = 0.0;
        *flow_area_out = 0.00001;
        return;
    }

    double bf;
    if (blockage > 0.9) {
        bf = 3.333 - 3.333 * blockage;
    } else {
        bf = 1.0 - 0.4012316798 * blockage - 0.3768350138 * (blockage * blockage);
    }

    // Inlet control
    double Q_inlet_unsubmerged = barrels * (0.421 * sqrt(p->g) * pow(bf * diameter, 0.87) * pow(driving_energy, 1.63));
    double Q_inlet_submerged = barrels * (0.530 * sqrt(p->g) * pow(bf * diameter, 1.87) * pow(driving_energy, 0.63));

    double Q = (Q_inlet_unsubmerged < Q_inlet_submerged) ? Q_inlet_unsubmerged : Q_inlet_submerged;

    // Critical depth of the adopted flow, per barrel (exact solve)
    double outlet_culvert_depth = circular_critical_depth(Q / barrels, bf * diameter, p->g);

    double flow_area, perimeter;
    double alpha;
    double bd = bf * diameter;

    if (outlet_culvert_depth >= bd) {
        // Pipe full
        outlet_culvert_depth = bd;
        flow_area = barrels * (bd / 2.0) * (bd / 2.0) * M_PI;
        perimeter = barrels * bd * M_PI;
    } else {
        // Partial fill
        alpha = acos(1.0 - 2.0 * outlet_culvert_depth / bd) * 2.0;
        flow_area = barrels * bd * bd / 8.0 * (alpha - sin(alpha));
        perimeter = barrels * (alpha * bd / 2.0);
    }

    // Outlet control
    if (delta_total_energy < driving_energy) {
        if (outlet_enquiry_depth > bd) {
            // Outlet submerged - pipe full
            outlet_culvert_depth = bd;
            flow_area = barrels * (bd / 2.0) * (bd / 2.0) * M_PI;
            perimeter = barrels * bd * M_PI;
        } else {
            // Partial flow - recalculate critical depth (per barrel)
            outlet_culvert_depth = circular_critical_depth(Q / barrels, bd, p->g);

            if (outlet_culvert_depth > bd) {
                outlet_culvert_depth = bd;
                flow_area = barrels * (bd / 2.0) * (bd / 2.0) * M_PI;
                perimeter = barrels * bd * M_PI;
            } else {
                alpha = acos(1.0 - 2.0 * outlet_culvert_depth / bd) * 2.0;
                flow_area = barrels * bd * bd / 8.0 * (alpha - sin(alpha));
                perimeter = barrels * alpha * bd / 2.0;
            }
        }
    }

    double hyd_rad = flow_area / perimeter;
    double culvert_velocity = sqrt(delta_total_energy / ((sum_loss / 2.0 / p->g) +
                                    (manning * manning * length) / pow(hyd_rad, 1.33333)));
    double Q_outlet_tailwater = flow_area * culvert_velocity;

    if (Q_outlet_tailwater < Q)
        Q = Q_outlet_tailwater;

    double barrel_velocity = Q / (flow_area + VELOCITY_PROTECTION / flow_area);

    *Q_out = Q;
    *barrel_velocity_out = barrel_velocity;
    *outlet_culvert_depth_out = outlet_culvert_depth;
    *flow_area_out = flow_area;
}

// ============================================================================
// Pure computation: Weir-orifice TRAPEZOID discharge
// Direct translation from weir_orifice_trapezoid_operator.py:
//   weir_orifice_trapezoid_function()
// Cross-section: bottom width `width`, side slopes z1 (left) and z2 (right).
// ============================================================================

// Newton iteration to find critical depth for a trapezoidal section given Q.
// Returns dcrit (clamped to [1e-5, depth] on convergence failure).
static double trapezoid_critical_depth(double Q, double bf_barrels_w,
                                       double z12, double sqrt_z1, double sqrt_z2,
                                       double depth, double g) {
    double dcrit = 1.0e-5;
    for (int ic = 0; ic < 100; ic++) {
        double Tc = bf_barrels_w + z12 * dcrit;
        double Ac = 0.5 * dcrit * (bf_barrels_w + Tc);
        if (Tc < 1.0e-12 || Ac < 1.0e-12) break;
        // Uses the domain gravity g (culvert_params.g <- domain.g). The Python
        // reference weir_orifice_trapezoid_function now derives g the same way
        // (from domain.g), so the two match for any g (e.g. non-Earth).
        double fc  = pow(Ac, 1.5) / sqrt(Tc) - Q / sqrt(g);
        double ffc = -0.5 * pow(Ac, 1.5) * z12 / pow(Tc, 1.5)
                     + 1.5 * sqrt(Ac) * sqrt(Tc);
        if (fabs(ffc) < 1.0e-30) break;
        double dyc = -fc / ffc;
        dcrit += dyc;
        if (dcrit < 1.0e-5) dcrit = 1.0e-5;
        if (fabs(dyc) < 1.0e-5) break;
    }
    if (dcrit > depth) dcrit = depth;
    return dcrit;
}


// Total flow area and wetted perimeter of `barrels` trapezoidal barrels, each
// filled to depth y: bottom + sloping sides, plus the roof (as wide as the top
// of the section) when running full.
static void trapezoid_section(double y, int full, double bottom, double z12,
                              double slant, double barrels,
                              double *area, double *perimeter) {
    double a = bottom * y + 0.5 * z12 * y * y;
    double p = bottom + slant * y + (full ? bottom + z12 * y : 0.0);
    *area = barrels * a;
    *perimeter = barrels * p;
}


void weir_orifice_trapezoid_discharge(const struct culvert_params *p,
                                      double driving_energy,
                                      double delta_total_energy,
                                      double outlet_enquiry_depth,
                                      double *Q_out, double *barrel_velocity_out,
                                      double *outlet_culvert_depth_out,
                                      double *flow_area_out) {
    double width   = p->width;
    double depth   = p->height;
    double blockage = p->blockage;
    double barrels  = p->barrels;
    double z1 = p->z1;
    double z2 = p->z2;
    double z12 = z1 + z2;
    double length  = p->length;
    double sum_loss = p->sum_loss;
    double manning  = p->manning;
    double g = p->g;

    if (blockage >= 1.0) {
        *Q_out = 0.0;
        *barrel_velocity_out = 0.0;
        *outlet_culvert_depth_out = 0.0;
        *flow_area_out = 1.0e-5;
        return;
    }

    double bf = 1.0 - blockage;
    // Every area and perimeter is that of ONE barrel times the number of
    // barrels (each barrel has its own sides and, running full, its own roof);
    // blockage narrows the bottom width. Same as the Python reference.
    double bottom = bf * width;
    double slant = sqrt(z1 * z1 + 1.0) + sqrt(z2 * z2 + 1.0);

    // Inlet control, per barrel times the number of barrels
    // Weir flow (unsubmerged): Q = 1.7 * mean width * driving_energy^1.5
    double top_w  = 2.0 * width + depth * z12;
    double Q_inlet_unsubmerged = 1.7 * bf * barrels * (top_w / 2.0)
                                 * pow(driving_energy, 1.5);
    // Orifice flow (submerged): Q = 0.8 * full area * sqrt(g * driving_energy).
    // The area is one unblocked barrel's: blockage and barrels are applied once.
    double full_area = 0.5 * depth * (2.0 * width + z12 * depth);
    double Q_inlet_submerged = 0.8 * bf * barrels * sqrt(g) * full_area
                               * sqrt(driving_energy);

    double Q = (Q_inlet_unsubmerged < Q_inlet_submerged) ? Q_inlet_unsubmerged
                                                         : Q_inlet_submerged;

    // Critical depth of the inlet-control flow (per barrel) and its section
    double outlet_culvert_depth = trapezoid_critical_depth(Q / barrels, bottom, z12,
                                                           0.0, 0.0, depth, g);
    double flow_area, perimeter;
    if (outlet_culvert_depth >= depth) {
        outlet_culvert_depth = depth;
        trapezoid_section(depth, 1, bottom, z12, slant, barrels, &flow_area, &perimeter);
    } else {
        trapezoid_section(outlet_culvert_depth, 0, bottom, z12, slant, barrels,
                          &flow_area, &perimeter);
    }

    // Outlet-control velocity and Q
    double hyd_rad = flow_area / fmax(perimeter, 1.0e-12);
    double culvert_velocity = sqrt(delta_total_energy
                                   / ((sum_loss / (2.0 * g))
                                      + (manning * manning * length)
                                        / pow(hyd_rad, 1.33333)));
    double Q_outlet_tailwater = flow_area * culvert_velocity;

    if (delta_total_energy < driving_energy) {
        // Outlet control
        if (outlet_enquiry_depth > depth) {
            // Outlet submerged: the barrel runs full
            outlet_culvert_depth = depth;
            trapezoid_section(depth, 1, bottom, z12, slant, barrels, &flow_area, &perimeter);
        } else {
            // Critical depth of the inlet-control flow, as in the box and pipe
            // (see the Python for why Q is no longer cut first).
            outlet_culvert_depth = trapezoid_critical_depth(Q / barrels, bottom, z12,
                                                            0.0, 0.0, depth, g);
            if (outlet_culvert_depth >= depth) {
                outlet_culvert_depth = depth;
                trapezoid_section(depth, 1, bottom, z12, slant, barrels, &flow_area, &perimeter);
            } else {
                trapezoid_section(outlet_culvert_depth, 0, bottom, z12, slant, barrels,
                                  &flow_area, &perimeter);
            }
        }

        hyd_rad = flow_area / fmax(perimeter, 1.0e-12);
        culvert_velocity = sqrt(delta_total_energy
                                / ((sum_loss / (2.0 * g))
                                   + (manning * manning * length)
                                     / pow(hyd_rad, 1.33333)));
        Q_outlet_tailwater = flow_area * culvert_velocity;
        Q = fmin(Q, Q_outlet_tailwater);
    } else {
        /* Head difference at or above the driving energy: the barrel's
         * friction still limits the flow, as for pipes (see the Python). */
        Q = fmin(Q, Q_outlet_tailwater);
    }

    double barrel_velocity = Q / (flow_area + VELOCITY_PROTECTION / fmax(flow_area, 1.0e-12));

    *Q_out = Q;
    *barrel_velocity_out = barrel_velocity;
    *outlet_culvert_depth_out = outlet_culvert_depth;
    *flow_area_out = flow_area;
}

// ============================================================================
// Energy smoothing (from boyd_box_operator.py:total_energy())
// ============================================================================

void culvert_smooth_energy(double *smooth_delta_total_energy,
                           double delta_total_energy,
                           double timestep,
                           double smoothing_timescale,
                           double *ts_out) {
    double ts;
    if (timestep > 0.0) {
        double denom = timestep;
        if (smoothing_timescale > denom) denom = smoothing_timescale;
        if (1.0e-06 > denom) denom = 1.0e-06;
        ts = timestep / denom;
    } else {
        ts = 1.0;
    }

    *smooth_delta_total_energy = *smooth_delta_total_energy +
        ts * (delta_total_energy - *smooth_delta_total_energy);
    *ts_out = ts;
}

// ============================================================================
// Discharge smoothing (from boyd_box_operator.py:smooth_discharge())
// ============================================================================

void culvert_smooth_discharge(double smooth_delta_total_energy,
                              double *smooth_Q,
                              double Q_in,
                              double flow_area,
                              double ts,
                              double *Q_out, double *velocity_out) {
    double Qsign = (smooth_delta_total_energy >= 0.0) ? 1.0 : -1.0;

    *smooth_Q = *smooth_Q + ts * (Q_in * Qsign - *smooth_Q);

    double Q;
    if ((*smooth_Q >= 0.0) != (Qsign >= 0.0)) {
        // Flow direction mismatch - set Q to zero
        Q = 0.0;
    } else {
        double abs_smooth_Q = fabs(*smooth_Q);
        Q = (abs_smooth_Q < Q_in) ? abs_smooth_Q : Q_in;
    }

    double barrel_velocity;
    if (flow_area == 0.0) {
        barrel_velocity = 0.0;
    } else {
        barrel_velocity = Q / flow_area;
    }

    *Q_out = Q;
    *velocity_out = barrel_velocity;
}

// ============================================================================
// Helper: compute enquiry-derived values from raw gathered data
// Mirrors inlet_enquiry.py get_enquiry_* methods
// ============================================================================

static void compute_enquiry_values(const struct inlet_data *data,
                                   const struct culvert_params *p,
                                   int inlet_idx,
                                   double *depth, double *velocity_head,
                                   double *total_energy, double *specific_energy) {
    double invert_elev;
    if (inlet_idx == 0 && p->has_invert_elevation_0)
        invert_elev = p->invert_elevation_0;
    else if (inlet_idx == 1 && p->has_invert_elevation_1)
        invert_elev = p->invert_elevation_1;
    else
        invert_elev = data->enquiry_elevation;

    double d = data->enquiry_stage - invert_elev;
    if (d < 0.0) d = 0.0;
    *depth = d;

    // Velocity head: 0.5 * speed² / g
    double water_depth = data->enquiry_stage - data->enquiry_elevation;
    double denom = water_depth * water_depth + VELOCITY_PROTECTION;
    double u = water_depth * data->enquiry_xmom / denom;
    double v = water_depth * data->enquiry_ymom / denom;
    double speed_sq = u * u + v * v;
    *velocity_head = 0.5 * speed_sq / p->g;

    *total_energy = *velocity_head + data->enquiry_stage;
    *specific_energy = *velocity_head + d;
}

// ============================================================================
// GPU Culvert Manager: Init / Finalize
// ============================================================================

// Free the host-side inlet staging arrays owned by one culvert slot.
static void culvert_free_inlet_staging(struct culvert_indices *ci) {
    if (ci->inlet0_indices) { free(ci->inlet0_indices); ci->inlet0_indices = NULL; }
    if (ci->inlet0_areas)   { free(ci->inlet0_areas);   ci->inlet0_areas   = NULL; }
    if (ci->inlet1_indices) { free(ci->inlet1_indices); ci->inlet1_indices = NULL; }
    if (ci->inlet1_areas)   { free(ci->inlet1_areas);   ci->inlet1_areas   = NULL; }
    ci->inlet0_num = 0;
    ci->inlet1_num = 0;
}

// Copy one inlet's triangle indices/areas into freshly allocated arrays.
// Returns 0 on success, -1 on allocation failure (leaving both outputs NULL).
static int culvert_copy_inlet_staging(int num, const int *indices, const double *areas,
                                      int **out_indices, double **out_areas) {
    *out_indices = NULL;
    *out_areas = NULL;
    if (num <= 0) return 0;

    int *idx = (int*)malloc(num * sizeof(int));
    double *ar = (double*)malloc(num * sizeof(double));
    if (!idx || !ar) {
        free(idx);
        free(ar);
        return -1;
    }
    memcpy(idx, indices, num * sizeof(int));
    memcpy(ar, areas, num * sizeof(double));
    *out_indices = idx;
    *out_areas = ar;
    return 0;
}

int gpu_culvert_init(struct gpu_domain *GD,
                     const struct culvert_params *params,
                     int enquiry_index_0, int enquiry_index_1,
                     int inlet0_num, int *inlet0_indices, double *inlet0_areas,
                     int inlet1_num, int *inlet1_indices, double *inlet1_areas,
                     int master_proc, int enquiry_proc_0, int enquiry_proc_1,
                     int inlet_master_proc_0, int inlet_master_proc_1,
                     int is_local, int mpi_tag_base,
                     double init_smooth_Q, double init_smooth_delta_total_energy) {

    struct culvert_operators *CO = &GD->culvert_ops;

    // Grow params/indices/state arrays if full
    if (CO->num_culverts >= CO->capacity) {
        int new_cap = CO->capacity == 0 ? MAX_CULVERTS : CO->capacity * 2;
        struct culvert_params *np = (struct culvert_params*)
            realloc(CO->params, new_cap * sizeof(struct culvert_params));
        struct culvert_indices *ni = (struct culvert_indices*)
            realloc(CO->indices, new_cap * sizeof(struct culvert_indices));
        struct culvert_state *ns = (struct culvert_state*)
            realloc(CO->state, new_cap * sizeof(struct culvert_state));
        if (!np || !ni || !ns) {
            gpu_set_error(GD, "ERROR: Failed to grow culvert_operators to %d slots", new_cap);
            // The realloc that DID succeed owns the surviving copy of the inlet
            // staging pointers; free them through whichever array that is, so the
            // teardown below does not leak them.
            struct culvert_indices *live = ni ? ni : CO->indices;
            if (live) {
                for (int c = 0; c < CO->num_culverts; c++) culvert_free_inlet_staging(&live[c]);
            }
            if (np) free(np); else if (CO->params) free(CO->params);
            if (ni) free(ni); else if (CO->indices) free(CO->indices);
            if (ns) free(ns); else if (CO->state) free(CO->state);
            CO->params = NULL; CO->indices = NULL; CO->state = NULL;
            return -1;
        }
        // Zero-init new entries so state is clean
        memset(np + CO->capacity, 0,
               (new_cap - CO->capacity) * sizeof(struct culvert_params));
        memset(ni + CO->capacity, 0,
               (new_cap - CO->capacity) * sizeof(struct culvert_indices));
        memset(ns + CO->capacity, 0,
               (new_cap - CO->capacity) * sizeof(struct culvert_state));
        CO->params = np;
        CO->indices = ni;
        CO->state = ns;
        CO->capacity = new_cap;
    }

    int id = CO->num_culverts;

    // Copy params
    CO->params[id] = *params;

    // Copy indices
    struct culvert_indices *ci = &CO->indices[id];
    ci->enquiry_index_0 = enquiry_index_0;
    ci->enquiry_index_1 = enquiry_index_1;

    // Re-registering into a slot that already holds staging (possible only if a
    // caller reuses a slot) would leak, so clear it first.
    culvert_free_inlet_staging(ci);

    if (culvert_copy_inlet_staging(inlet0_num, inlet0_indices, inlet0_areas,
                                   &ci->inlet0_indices, &ci->inlet0_areas) != 0 ||
        culvert_copy_inlet_staging(inlet1_num, inlet1_indices, inlet1_areas,
                                   &ci->inlet1_indices, &ci->inlet1_areas) != 0) {
        gpu_set_error(GD, "ERROR: Failed to allocate inlet staging for %d/%d triangles",
                inlet0_num, inlet1_num);
        culvert_free_inlet_staging(ci);
        return -1;
    }

    ci->inlet0_num = inlet0_num;
    ci->inlet0_total_area = 0.0;
    for (int k = 0; k < inlet0_num; k++) ci->inlet0_total_area += inlet0_areas[k];

    ci->inlet1_num = inlet1_num;
    ci->inlet1_total_area = 0.0;
    for (int k = 0; k < inlet1_num; k++) ci->inlet1_total_area += inlet1_areas[k];

    // MPI topology
    ci->master_proc = master_proc;
    ci->enquiry_proc[0] = enquiry_proc_0;
    ci->enquiry_proc[1] = enquiry_proc_1;
    ci->inlet_master_proc[0] = inlet_master_proc_0;
    ci->inlet_master_proc[1] = inlet_master_proc_1;
    ci->is_local = is_local;
    ci->mpi_tag_base = mpi_tag_base;

    // Initialize smoothing state from CPU operator's pre-seeded values.
    // Python __init__ runs discharge_routine() once to seed these (boyd_box_operator.py:193);
    // matters when smoothing_timescale > 0 so the GPU side starts from the same smoothed state.
    CO->state[id].smooth_delta_total_energy = init_smooth_delta_total_energy;
    CO->state[id].smooth_Q = init_smooth_Q;

    CO->num_culverts++;
    return id;
}

void gpu_culvert_finalize(struct gpu_domain *GD, int culvert_id) {
    // Deliberately a no-op. The slot's inlet staging arrays are freed by
    // gpu_culverts_finalize_all(); releasing them here would leave a slot that
    // is still counted in num_culverts — and still visited by the batched
    // gather/scatter — holding freed pointers.
    (void)GD;
    (void)culvert_id;
}

void gpu_culverts_finalize_all(struct gpu_domain *GD) {
    struct culvert_operators *CO = &GD->culvert_ops;

    if (CO->mapped) {
        int ne = 2 * CO->num_culverts;
        int nt = CO->total_inlet_triangles;

        int *eid = CO->scratch_enquiry_indices;
        double *ss = CO->scratch_stage;
        double *sx = CO->scratch_xmom;
        double *sy = CO->scratch_ymom;
        double *se = CO->scratch_elev;
        double *as = CO->scratch_avg_stage;
        double *ad = CO->scratch_avg_depth;
        double *ax = CO->scratch_avg_xmom;
        double *ay = CO->scratch_avg_ymom;
        double *nd = CO->scratch_slot_shift;
        double *nx = CO->scratch_slot_xmom;
        double *ny = CO->scratch_slot_ymom;
        int *sst = CO->scratch_slot_start;
        int *scn = CO->scratch_slot_count;

        if (ne > 0) {
            #pragma omp target exit data map(delete: eid[0:ne], sst[0:ne], scn[0:ne], \
                ss[0:ne], sx[0:ne], sy[0:ne], se[0:ne], \
                as[0:ne], ad[0:ne], ax[0:ne], ay[0:ne], \
                nd[0:ne], nx[0:ne], ny[0:ne])
        }

        if (nt > 0) {
            int *si = CO->scratch_inlet_indices;
            double *sa = CO->scratch_inlet_areas;
            #pragma omp target exit data map(delete: si[0:nt], sa[0:nt])
        }
        if (CO->tracer_scratch_ns > 0) {
            int ne_t = 2 * CO->num_culverts;
            int nt_t = CO->total_inlet_triangles;
            int nk = CO->tracer_scratch_ns;
            double *th = CO->scratch_tr_hold;
            double *tT = CO->scratch_tr_T;
            double *tG = CO->scratch_tr_G;
            double *tc = CO->scratch_tr_c;
            #pragma omp target exit data map(delete: th[0:nt_t], tT[0:ne_t*nk], \
                tG[0:ne_t], tc[0:CO->num_culverts*nk])
        }
        CO->mapped = 0;
    }

    if (CO->scratch_enquiry_indices) { free(CO->scratch_enquiry_indices); CO->scratch_enquiry_indices = NULL; }
    if (CO->scratch_stage) { free(CO->scratch_stage); CO->scratch_stage = NULL; }
    if (CO->scratch_xmom) { free(CO->scratch_xmom); CO->scratch_xmom = NULL; }
    if (CO->scratch_ymom) { free(CO->scratch_ymom); CO->scratch_ymom = NULL; }
    if (CO->scratch_elev) { free(CO->scratch_elev); CO->scratch_elev = NULL; }
    if (CO->scratch_avg_stage) { free(CO->scratch_avg_stage); CO->scratch_avg_stage = NULL; }
    if (CO->scratch_avg_depth) { free(CO->scratch_avg_depth); CO->scratch_avg_depth = NULL; }
    if (CO->scratch_avg_xmom) { free(CO->scratch_avg_xmom); CO->scratch_avg_xmom = NULL; }
    if (CO->scratch_avg_ymom) { free(CO->scratch_avg_ymom); CO->scratch_avg_ymom = NULL; }
    if (CO->scratch_slot_shift) { free(CO->scratch_slot_shift); CO->scratch_slot_shift = NULL; }
    if (CO->scratch_slot_xmom) { free(CO->scratch_slot_xmom); CO->scratch_slot_xmom = NULL; }
    if (CO->scratch_slot_ymom) { free(CO->scratch_slot_ymom); CO->scratch_slot_ymom = NULL; }
    if (CO->scratch_inlet_indices) { free(CO->scratch_inlet_indices); CO->scratch_inlet_indices = NULL; }
    if (CO->scratch_inlet_areas) { free(CO->scratch_inlet_areas); CO->scratch_inlet_areas = NULL; }
    if (CO->scratch_slot_start) { free(CO->scratch_slot_start); CO->scratch_slot_start = NULL; }
    if (CO->scratch_slot_count) { free(CO->scratch_slot_count); CO->scratch_slot_count = NULL; }
    free(CO->scratch_tr_hold); CO->scratch_tr_hold = NULL;
    free(CO->scratch_tr_T);    CO->scratch_tr_T = NULL;
    free(CO->scratch_tr_G);    CO->scratch_tr_G = NULL;
    free(CO->scratch_tr_c);    CO->scratch_tr_c = NULL;
    CO->tracer_scratch_ns = 0;

    culvert_host_scratch_free(CO);

    if (CO->params)  { free(CO->params);  CO->params  = NULL; }
    if (CO->indices) {
        for (int c = 0; c < CO->num_culverts; c++) {
            culvert_free_inlet_staging(&CO->indices[c]);
        }
        free(CO->indices); CO->indices = NULL;
    }
    if (CO->state)   { free(CO->state);   CO->state   = NULL; }
    CO->num_culverts = 0;
    CO->capacity = 0;
    CO->initialized = 0;
}

// ============================================================================
// GPU Mapping: allocate and map scratch buffers
// Call AFTER all culverts are registered AND GPU domain is initialized
// ============================================================================

int gpu_culverts_map(struct gpu_domain *GD) {
    struct culvert_operators *CO = &GD->culvert_ops;

    if (CO->num_culverts == 0) return 0;
    if (CO->mapped) return 0;

    omp_set_default_device(gpu_compute_device(GD));

    int nc = CO->num_culverts;
    int ne = 2 * nc;  // 2 enquiry/inlet slots per culvert

    // --- Enquiry scratch: constant indices + per-step gathered values ---
    CO->scratch_enquiry_indices = (int*)calloc(ne, sizeof(int));
    CO->scratch_stage = (double*)calloc(ne, sizeof(double));
    CO->scratch_xmom = (double*)calloc(ne, sizeof(double));
    CO->scratch_ymom = (double*)calloc(ne, sizeof(double));
    CO->scratch_elev = (double*)calloc(ne, sizeof(double));
    if (!CO->scratch_enquiry_indices || !CO->scratch_stage || !CO->scratch_xmom ||
        !CO->scratch_ymom || !CO->scratch_elev) {
        gpu_set_error(GD, "gpu_culverts_map: could not allocate enquiry scratch for %d culverts", nc);
        return -1;
    }

    // Enquiry indices are constant for the life of the domain. Remote enquiry
    // points (index < 0) are parked at 0; their gathered values are overwritten
    // by MPI later. Filled once here, mapped map(to:) once below.
    for (int c = 0; c < nc; c++) {
        int ei0 = CO->indices[c].enquiry_index_0;
        int ei1 = CO->indices[c].enquiry_index_1;
        CO->scratch_enquiry_indices[2 * c]     = (ei0 >= 0) ? ei0 : 0;
        CO->scratch_enquiry_indices[2 * c + 1] = (ei1 >= 0) ? ei1 : 0;
    }

    // --- Per-inlet reduction / scatter accumulators (2 per culvert) ---
    CO->scratch_avg_stage = (double*)calloc(ne, sizeof(double));
    CO->scratch_avg_depth = (double*)calloc(ne, sizeof(double));
    CO->scratch_avg_xmom  = (double*)calloc(ne, sizeof(double));
    CO->scratch_avg_ymom  = (double*)calloc(ne, sizeof(double));
    CO->scratch_slot_shift = (double*)calloc(ne, sizeof(double));
    CO->scratch_slot_xmom  = (double*)calloc(ne, sizeof(double));
    CO->scratch_slot_ymom  = (double*)calloc(ne, sizeof(double));
    if (!CO->scratch_avg_stage || !CO->scratch_avg_depth || !CO->scratch_avg_xmom ||
        !CO->scratch_avg_ymom || !CO->scratch_slot_shift || !CO->scratch_slot_xmom ||
        !CO->scratch_slot_ymom) {
        gpu_set_error(GD, "gpu_culverts_map: could not allocate inlet accumulators for %d culverts", nc);
        return -1;
    }

    // --- Flattened inlet-triangle metadata (constant) ---
    CO->total_inlet_triangles = 0;
    for (int c = 0; c < nc; c++) {
        CO->total_inlet_triangles += CO->indices[c].inlet0_num + CO->indices[c].inlet1_num;
    }

    int nt = CO->total_inlet_triangles;
    CO->scratch_inlet_indices = (int*)calloc(nt, sizeof(int));
    CO->scratch_inlet_areas = (double*)calloc(nt, sizeof(double));
    CO->scratch_slot_start = (int*)calloc(ne, sizeof(int));
    CO->scratch_slot_count = (int*)calloc(ne, sizeof(int));
    if (!CO->scratch_inlet_indices || !CO->scratch_inlet_areas ||
        !CO->scratch_slot_start || !CO->scratch_slot_count) {
        gpu_set_error(GD, "gpu_culverts_map: could not allocate inlet index scratch for %d culverts", nc);
        return -1;
    }

    // Flatten inlet indices/areas and record each inlet's contiguous range.
    int offset = 0;
    for (int c = 0; c < nc; c++) {
        struct culvert_indices *ci = &CO->indices[c];

        CO->scratch_slot_start[2 * c] = offset;
        CO->scratch_slot_count[2 * c] = ci->inlet0_num;
        for (int k = 0; k < ci->inlet0_num; k++) {
            CO->scratch_inlet_indices[offset] = ci->inlet0_indices[k];
            CO->scratch_inlet_areas[offset] = ci->inlet0_areas[k];
            offset++;
        }

        CO->scratch_slot_start[2 * c + 1] = offset;
        CO->scratch_slot_count[2 * c + 1] = ci->inlet1_num;
        for (int k = 0; k < ci->inlet1_num; k++) {
            CO->scratch_inlet_indices[offset] = ci->inlet1_indices[k];
            CO->scratch_inlet_areas[offset] = ci->inlet1_areas[k];
            offset++;
        }
    }

    // --- Map everything to the device ONCE ---
    int *eid = CO->scratch_enquiry_indices;
    double *ss = CO->scratch_stage;
    double *sx = CO->scratch_xmom;
    double *sy = CO->scratch_ymom;
    double *se = CO->scratch_elev;
    double *as = CO->scratch_avg_stage;
    double *ad = CO->scratch_avg_depth;
    double *ax = CO->scratch_avg_xmom;
    double *ay = CO->scratch_avg_ymom;
    double *nd = CO->scratch_slot_shift;
    double *nx = CO->scratch_slot_xmom;
    double *ny = CO->scratch_slot_ymom;
    int *sst = CO->scratch_slot_start;
    int *scn = CO->scratch_slot_count;
    #pragma omp target enter data map(to: eid[0:ne], sst[0:ne], scn[0:ne]) \
        map(alloc: ss[0:ne], sx[0:ne], sy[0:ne], se[0:ne], \
                   as[0:ne], ad[0:ne], ax[0:ne], ay[0:ne], \
                   nd[0:ne], nx[0:ne], ny[0:ne])

    if (nt > 0) {
        int *si = CO->scratch_inlet_indices;
        double *sa = CO->scratch_inlet_areas;
        #pragma omp target enter data map(to: si[0:nt], sa[0:nt])
    }

    CO->mapped = 1;
    CO->initialized = 1;
    return 0;
}

// ============================================================================
// Batched Gather: read enquiry + inlet data from GPU in TWO transfers
// ============================================================================

static void gpu_culvert_gather_enquiry(struct gpu_domain *GD,
                                       struct inlet_data *data0,
                                       struct inlet_data *data1) {
    struct culvert_operators *CO = &GD->culvert_ops;
    int nc = CO->num_culverts;
    int ne = 2 * nc;

    double * restrict stage_c = GD->D.stage_centroid_values;
    double * restrict xmom_c = GD->D.xmom_centroid_values;
    double * restrict ymom_c = GD->D.ymom_centroid_values;
    double * restrict bed_c = GD->D.bed_centroid_values;

    double *ss = CO->scratch_stage;
    double *sx = CO->scratch_xmom;
    double *sy = CO->scratch_ymom;
    double *se = CO->scratch_elev;

    // Enquiry indices are persistently mapped (map(to:) in gpu_culverts_map);
    // remote enquiry points were parked at index 0 there and are overwritten by
    // MPI later. Gather straight from the device-resident index buffer — no
    // per-step host allocation, no per-step map.
    int *eid = CO->scratch_enquiry_indices;
    OMP_PARALLEL_LOOP
    for (int k = 0; k < ne; k++) {
        int i = eid[k];
        ss[k] = stage_c[i];
        sx[k] = xmom_c[i];
        sy[k] = ymom_c[i];
        se[k] = bed_c[i];
    }

    // Single D2H transfer (~1KB for 20 culverts)
    #pragma omp target update from(ss[0:ne], sx[0:ne], sy[0:ne], se[0:ne])

    // Unpack into per-culvert inlet_data structs
    for (int c = 0; c < nc; c++) {
        data0[c].enquiry_stage = ss[2 * c];
        data0[c].enquiry_xmom = sx[2 * c];
        data0[c].enquiry_ymom = sy[2 * c];
        data0[c].enquiry_elevation = se[2 * c];

        data1[c].enquiry_stage = ss[2 * c + 1];
        data1[c].enquiry_xmom = sx[2 * c + 1];
        data1[c].enquiry_ymom = sy[2 * c + 1];
        data1[c].enquiry_elevation = se[2 * c + 1];
    }
}

static void gpu_culvert_gather_inlets(struct gpu_domain *GD,
                                      struct inlet_data *data0,
                                      struct inlet_data *data1) {
    struct culvert_operators *CO = &GD->culvert_ops;
    int nc = CO->num_culverts;
    int nt = CO->total_inlet_triangles;

    if (nt == 0) return;

    double * restrict stage_c = GD->D.stage_centroid_values;
    double * restrict xmom_c = GD->D.xmom_centroid_values;
    double * restrict ymom_c = GD->D.ymom_centroid_values;
    double * restrict bed_c = GD->D.bed_centroid_values;

    int ne = 2 * nc;
    int *si = CO->scratch_inlet_indices;
    double *sa = CO->scratch_inlet_areas;
    int *sst = CO->scratch_slot_start;
    int *scn = CO->scratch_slot_count;
    double *as = CO->scratch_avg_stage;
    double *ad = CO->scratch_avg_depth;
    double *ax = CO->scratch_avg_xmom;
    double *ay = CO->scratch_avg_ymom;

    // On-device area-weighted reduction: one team per inlet (ne total), each
    // summing its contiguous triangle range sequentially. No atomics, and the
    // summation order matches the old host loop exactly. Only the per-inlet
    // sums (2*nc doubles ×4 ≈ a few KB) travel back to the host — not every
    // triangle value.
    #pragma omp target teams distribute parallel for
    for (int s = 0; s < ne; s++) {
        double sum_stage = 0.0, sum_depth = 0.0, sum_xmom = 0.0, sum_ymom = 0.0;
        int start = sst[s];
        int cnt = scn[s];
        for (int j = 0; j < cnt; j++) {
            int k = start + j;
            int i = si[k];
            double area = sa[k];
            double depth = stage_c[i] - bed_c[i];
            if (depth < 0.0) depth = 0.0;
            sum_stage += stage_c[i] * area;
            sum_depth += depth * area;
            sum_xmom += xmom_c[i] * area;
            sum_ymom += ymom_c[i] * area;
        }
        as[s] = sum_stage;
        ad[s] = sum_depth;
        ax[s] = sum_xmom;
        ay[s] = sum_ymom;
    }

    // Single D2H transfer of the per-inlet sums.
    #pragma omp target update from(as[0:ne], ad[0:ne], ax[0:ne], ay[0:ne])

    // Divide by (constant, host-side) inlet area to get averages. An inlet
    // with zero local area (e.g. a cross-boundary inlet this rank doesn't own)
    // yields zeros here and is overwritten by the MPI exchange.
    for (int c = 0; c < nc; c++) {
        struct culvert_indices *ci = &CO->indices[c];
        double a0 = ci->inlet0_total_area;
        double a1 = ci->inlet1_total_area;
        int s0 = 2 * c, s1 = 2 * c + 1;

        if (a0 > 0.0) {
            data0[c].avg_stage = as[s0] / a0;
            data0[c].avg_depth = ad[s0] / a0;
            data0[c].avg_xmom  = ax[s0] / a0;
            data0[c].avg_ymom  = ay[s0] / a0;
        } else {
            data0[c].avg_stage = 0.0; data0[c].avg_depth = 0.0;
            data0[c].avg_xmom  = 0.0; data0[c].avg_ymom  = 0.0;
        }
        data0[c].total_area = a0;

        if (a1 > 0.0) {
            data1[c].avg_stage = as[s1] / a1;
            data1[c].avg_depth = ad[s1] / a1;
            data1[c].avg_xmom  = ax[s1] / a1;
            data1[c].avg_ymom  = ay[s1] / a1;
        } else {
            data1[c].avg_stage = 0.0; data1[c].avg_depth = 0.0;
            data1[c].avg_xmom  = 0.0; data1[c].avg_ymom  = 0.0;
        }
        data1[c].total_area = a1;
    }
}

// ============================================================================
// Batched Scatter: write updated depths/momenta back to GPU
// ============================================================================

// Level one inlet to a new average depth — water finds its level (issue #229).
//
// The single implementation of the structure write-back on the device, matching
// level_stages_to_average() + Inlet.set_average_momenta() on the host path.
//
// Stage: adding volume raises the LOWEST stages to a common level; removing
// volume lowers the HIGHEST stages to a common level, clamping each cell at its
// bed (if the inlet holds less than asked for, it is drained dry and no more).
// A zero transfer leaves the stages untouched, which is what makes the update
// well balanced: a lake at rest on a sloping bed is not disturbed. Writing a
// uniform DEPTH (the old behaviour) tilted the surface onto the bed; a shape-
// preserving SHIFT (tried first) destabilised discharge into an initially dry
// inlet. The level is found by bisection on the monotone volume(L) function;
// 100 iterations pins it to the last bit of a double.
//
// Momentum: the physics produces ONE momentum value per inlet. It is written
// depth-weighted — cell i gets m * depth_i / average_depth — so the VELOCITY
// field is uniform and the area-weighted average momentum is exactly m. A
// uniform-momentum write over the now non-uniform depths would give a nearly
// dry cell an enormous velocity and collapse the global timestep.
#pragma omp declare target
static void culvert_level_inlet_surface(const int * restrict idx,
                                        const double * restrict areas,
                                        int ntri, double delta_avg_depth,
                                        double new_xmom, double new_ymom,
                                        double * restrict stage_c,
                                        const double * restrict bed_c,
                                        double * restrict xmom_c,
                                        double * restrict ymom_c) {
    if (ntri <= 0) return;

    double total_area = 0.0;
    for (int j = 0; j < ntri; j++) total_area += areas[j];
    double volume = delta_avg_depth * total_area;

    if (volume > 0.0 && total_area > 0.0) {
        // Fill: find L with  A(L) = sum a_i * max(L - s_i, 0) == volume.
        // A(max_s + volume/total_area) >= volume, so L lies in [lo, hi].
        double lo = stage_c[idx[0]];
        double hi = stage_c[idx[0]];
        for (int j = 1; j < ntri; j++) {
            double sj = stage_c[idx[j]];
            if (sj < lo) lo = sj;
            if (sj > hi) hi = sj;
        }
        hi += volume / total_area;
        for (int it = 0; it < 100; it++) {
            double mid = 0.5 * (lo + hi);
            double added = 0.0;
            for (int j = 0; j < ntri; j++) {
                double a = mid - stage_c[idx[j]];
                if (a > 0.0) added += a * areas[j];
            }
            if (added < volume) lo = mid; else hi = mid;
        }
        for (int j = 0; j < ntri; j++) {
            int i = idx[j];
            if (stage_c[i] < hi) stage_c[i] = hi;
        }
    } else if (volume < 0.0 && total_area > 0.0) {
        // Drawdown: find L with
        //   R(L) = sum a_i * min(depth_i, max(s_i - L, 0)) == -volume,
        // then s_i' = max(bed_i, min(s_i, L)).
        double to_remove = -volume;
        double water = 0.0;
        double lo = bed_c[idx[0]];
        double hi = stage_c[idx[0]];
        for (int j = 0; j < ntri; j++) {
            int i = idx[j];
            double d = stage_c[i] - bed_c[i];
            if (d > 0.0) water += d * areas[j];
            if (bed_c[i] < lo) lo = bed_c[i];
            if (stage_c[i] > hi) hi = stage_c[i];
        }
        if (to_remove >= water) {
            // Asked for more than the inlet holds: drain it dry, no further.
            for (int j = 0; j < ntri; j++) {
                int i = idx[j];
                if (stage_c[i] > bed_c[i]) stage_c[i] = bed_c[i];
            }
        } else {
            for (int it = 0; it < 100; it++) {
                double mid = 0.5 * (lo + hi);
                double removed = 0.0;
                for (int j = 0; j < ntri; j++) {
                    int i = idx[j];
                    double d = stage_c[i] - bed_c[i];
                    if (d < 0.0) d = 0.0;
                    double r = stage_c[i] - mid;
                    if (r < 0.0) r = 0.0;
                    if (r > d) r = d;
                    removed += r * areas[j];
                }
                if (removed > to_remove) lo = mid; else hi = mid;
            }
            for (int j = 0; j < ntri; j++) {
                int i = idx[j];
                double s_new = stage_c[i] < lo ? stage_c[i] : lo;
                if (s_new < bed_c[i]) s_new = bed_c[i];
                stage_c[i] = s_new;
            }
        }
    }
    // volume == 0: stages untouched (the exact lake-at-rest no-op).

    // Depth-weighted momentum over the post-level depths.
    double avg_depth = 0.0;
    for (int j = 0; j < ntri; j++) {
        int i = idx[j];
        double d = stage_c[i] - bed_c[i];
        if (d > 0.0) avg_depth += d * areas[j];
    }
    avg_depth = (total_area > 0.0) ? avg_depth / total_area : 0.0;
    for (int j = 0; j < ntri; j++) {
        int i = idx[j];
        double d = stage_c[i] - bed_c[i];
        double w = (avg_depth > 0.0 && d > 0.0) ? d / avg_depth : 0.0;
        xmom_c[i] = new_xmom * w;
        ymom_c[i] = new_ymom * w;
    }
}
#pragma omp end declare target


// Tracer scratch for the transfer, sized for ns tracers and mapped once.
// Returns 0 on success; on failure the caller carries no tracers this step.
static int culvert_tracer_scratch(struct gpu_domain *GD, int ns) {
    struct culvert_operators *CO = &GD->culvert_ops;
    if (CO->tracer_scratch_ns == ns) return 0;
    int nc = CO->num_culverts, ne = 2 * nc, nt = CO->total_inlet_triangles;
    if (CO->tracer_scratch_ns > 0) {
        int nk = CO->tracer_scratch_ns;
        double *th = CO->scratch_tr_hold, *tT = CO->scratch_tr_T;
        double *tG = CO->scratch_tr_G, *tc = CO->scratch_tr_c;
        #pragma omp target exit data map(delete: th[0:nt], tT[0:ne*nk], tG[0:ne], tc[0:nc*nk])
        free(th); free(tT); free(tG); free(tc);
        CO->tracer_scratch_ns = 0;
    }
    CO->scratch_tr_hold = (double*)calloc(nt > 0 ? nt : 1, sizeof(double));
    CO->scratch_tr_T = (double*)calloc((size_t)ne * ns, sizeof(double));
    CO->scratch_tr_G = (double*)calloc(ne, sizeof(double));
    CO->scratch_tr_c = (double*)calloc((size_t)nc * ns, sizeof(double));
    if (!CO->scratch_tr_hold || !CO->scratch_tr_T || !CO->scratch_tr_G || !CO->scratch_tr_c) {
        gpu_set_error(GD, "culvert tracers: could not allocate scratch for %d tracers", ns);
        free(CO->scratch_tr_hold); free(CO->scratch_tr_T);
        free(CO->scratch_tr_G); free(CO->scratch_tr_c);
        CO->scratch_tr_hold = CO->scratch_tr_T = CO->scratch_tr_G = CO->scratch_tr_c = NULL;
        return -1;
    }
    double *th = CO->scratch_tr_hold, *tT = CO->scratch_tr_T;
    double *tG = CO->scratch_tr_G, *tc = CO->scratch_tr_c;
    int ntm = nt > 0 ? nt : 1;
    #pragma omp target enter data map(alloc: th[0:ntm], tT[0:ne*ns], tG[0:ne], tc[0:nc*ns])
    CO->tracer_scratch_ns = ns;
    return 0;
}

static void gpu_culvert_scatter(struct gpu_domain *GD,
                                struct culvert_transfer *transfers) {
    struct culvert_operators *CO = &GD->culvert_ops;
    int nc = CO->num_culverts;
    int nt = CO->total_inlet_triangles;
    int ne = 2 * nc;

    if (nt == 0) return;

    // Build ONE (delta average depth, xmom, ymom) triple per inlet on the host —
    // the physics already produced these as a single value per inlet region.
    // Inlet inlet_local is the inflow when inlet_local == inflow_idx, else the
    // outflow.
    //
    // The delta is (new average depth - the average depth the physics started
    // from); the kernel levels the inlet to absorb exactly that volume change.
    // Writing the new average depth to every cell instead would flatten the
    // water surface onto the bed and tilt a lake at rest (issue #229).
    double *nd = CO->scratch_slot_shift;
    double *nx = CO->scratch_slot_xmom;
    double *ny = CO->scratch_slot_ymom;
    for (int c = 0; c < nc; c++) {
        struct culvert_transfer *t = &transfers[c];
        for (int inlet = 0; inlet < 2; inlet++) {
            int s = 2 * c + inlet;
            double old_avg_depth = (inlet == 0) ? CO->host_data0[c].avg_depth
                                                : CO->host_data1[c].avg_depth;
            if (inlet == t->inflow_idx) {
                nd[s] = t->new_inflow_depth - old_avg_depth;
                nx[s] = t->new_inflow_xmom;
                ny[s] = t->new_inflow_ymom;
            } else {
                nd[s] = t->new_outflow_depth - old_avg_depth;
                nx[s] = t->new_outflow_xmom;
                ny[s] = t->new_outflow_ymom;
            }
        }
    }

    // Single H2D transfer of the per-inlet values (2*nc doubles ×3).
    #pragma omp target update to(nd[0:ne], nx[0:ne], ny[0:ne])

    // On-device scatter: one team per inlet writes its contiguous triangle
    // range, reading bed elevation straight from the domain array so stage =
    // bed + depth is computed on-device (no gathered bed buffer needed).
    int *si = CO->scratch_inlet_indices;
    int *sst = CO->scratch_slot_start;
    int *scn = CO->scratch_slot_count;
    double * restrict stage_c = GD->D.stage_centroid_values;
    double * restrict xmom_c = GD->D.xmom_centroid_values;
    double * restrict ymom_c = GD->D.ymom_centroid_values;
    double * restrict bed_c = GD->D.bed_centroid_values;

    double *sa = CO->scratch_inlet_areas;

    // Tracers ride with the water (anuga/structures/inlet_tracers.py,
    // StructureTracers): record each inlet cell's depth before it is levelled.
    // Only culverts local to this rank reach here with tracers present; the
    // Python manager leaves cross-boundary ones to the host MPI path.
    const anuga_int ns_tr = GD->D.number_of_tracers;
    const anuga_int N_tr = GD->D.number_of_elements;
    int carry = (ns_tr > 0 && GD->D.tracer_conserved_values != NULL);
    if (carry && culvert_tracer_scratch(GD, (int)ns_tr) != 0) carry = 0;
    double *hold = CO->scratch_tr_hold;
    if (carry) {
        #pragma omp target teams distribute parallel for
        for (int j = 0; j < nt; j++) {
            int i = si[j];
            double dj = stage_c[i] - bed_c[i];
            hold[j] = dj > 0.0 ? dj : 0.0;
        }
    }

    #pragma omp target teams distribute parallel for
    for (int s = 0; s < ne; s++) {
        culvert_level_inlet_surface(si + sst[s], sa + sst[s], scn[s],
                                    nd[s], nx[s], ny[s],
                                    stage_c, bed_c, xmom_c, ymom_c);
    }

    if (carry) {
        // Cells that lost water keep their concentration; sum the tracer they
        // gave up and the water gained, per inlet (one thread per inlet).
        const anuga_int ns = ns_tr, N = N_tr;
        const double hmin = GD->D.minimum_allowed_height;
        double * restrict t_cons = GD->D.tracer_conserved_values;
        double * restrict t_cv = GD->D.tracer_centroid_values;
        double *tT = CO->scratch_tr_T;
        double *tG = CO->scratch_tr_G;
        double *tc = CO->scratch_tr_c;
        #pragma omp target teams distribute parallel for
        for (int s = 0; s < ne; s++) {
            double g = 0.0;
            for (anuga_int k = 0; k < ns; k++) tT[s * ns + k] = 0.0;
            for (int j = sst[s]; j < sst[s] + scn[s]; j++) {
                int i = si[j];
                double hn = stage_c[i] - bed_c[i];
                if (hn < 0.0) hn = 0.0;
                double dh = hn - hold[j];
                if (dh < 0.0) {
                    double keep = hold[j] > 0.0 ? hn / hold[j] : 0.0;
                    double inv_h = hn > hmin ? 1.0 / hn : 0.0;
                    for (anuga_int k = 0; k < ns; k++) {
                        double mo = t_cons[k * N + i];
                        double mn = mo * keep;
                        tT[s * ns + k] += (mo - mn) * sa[j];
                        t_cons[k * N + i] = mn;
                        t_cv[k * N + i] = mn * inv_h;
                    }
                } else if (dh > 0.0) {
                    g += dh * sa[j];
                }
            }
            tG[s] = g;
        }
        #pragma omp target update from(tT[0:ne*ns], tG[0:ne])

        // Each culvert carries (tracer given up) / (water gained).
        for (int c = 0; c < nc; c++) {
            double G = tG[2 * c] + tG[2 * c + 1];
            for (anuga_int k = 0; k < ns; k++) {
                double T = tT[(2 * c) * ns + k] + tT[(2 * c + 1) * ns + k];
                tc[c * ns + k] = (G > 0.0) ? T / G : 0.0;
            }
        }
        #pragma omp target update to(tc[0:nc*ns])

        // Cells that gained water receive it at that concentration.
        #pragma omp target teams distribute parallel for
        for (int s = 0; s < ne; s++) {
            int c = s / 2;
            for (int j = sst[s]; j < sst[s] + scn[s]; j++) {
                int i = si[j];
                double hn = stage_c[i] - bed_c[i];
                if (hn < 0.0) hn = 0.0;
                double dh = hn - hold[j];
                if (dh > 0.0) {
                    double inv_h = hn > hmin ? 1.0 / hn : 0.0;
                    for (anuga_int k = 0; k < ns; k++) {
                        double mn = t_cons[k * N + i] + tc[c * ns + k] * dh;
                        t_cons[k * N + i] = mn;
                        t_cv[k * N + i] = mn * inv_h;
                    }
                }
            }
        }
    }
}

// ============================================================================
// MPI exchange helpers for cross-boundary culverts
// Uses non-blocking MPI to avoid deadlocks when multiple culverts cross
// different rank boundaries.
//
// Message protocol per cross-boundary culvert:
//   enquiry_proc[i] → master: 4 doubles (stage, xmom, ymom, elev)  tag_base+i
//   inlet_master[i] → master: 5 doubles (sum_s, sum_d, sum_xm, sum_ym, area)  tag_base+2+i
//   master → inlet_master[i]: 3 doubles (new_depth, new_xmom, new_ymom)  tag_base+4+i
// ============================================================================

// MPI message buffers for cross-boundary exchange.  Sized to the culvert
// count at allocation time (see culvert_host_scratch_ensure) rather than to
// MAX_CULVERTS, which is only the initial capacity.  The pointer-to-array
// types keep the [culvert][inlet][field] indexing of the original arrays.
struct culvert_mpi_bufs {
    int capacity;                  // culverts these buffers can hold
    int nreq_capacity;             // entries in requests[]
    double (*enquiry_send)[2][4];  // [culvert][inlet][stage,xmom,ymom,elev]
    double (*enquiry_recv)[2][4];
    double (*inlet_send)[2][5];    // [culvert][inlet][sum_s,sum_d,sum_xm,sum_ym,area]
    double (*inlet_recv)[2][5];
    double (*result_send)[2][3];   // [culvert][inlet][new_depth,new_xmom,new_ymom]
    double (*result_recv)[2][3];
    MPI_Request *requests;         // 6 per culvert (2 enquiry + 2 inlet + 2 result)
};

static void culvert_host_scratch_free(struct culvert_operators *CO) {
    if (CO->host_data0) { free(CO->host_data0); CO->host_data0 = NULL; }
    if (CO->host_data1) { free(CO->host_data1); CO->host_data1 = NULL; }
    if (CO->host_results) { free(CO->host_results); CO->host_results = NULL; }
    if (CO->host_transfers) { free(CO->host_transfers); CO->host_transfers = NULL; }
    if (CO->host_mpi_bufs) {
        struct culvert_mpi_bufs *b = CO->host_mpi_bufs;
        free(b->enquiry_send); free(b->enquiry_recv);
        free(b->inlet_send);   free(b->inlet_recv);
        free(b->result_send);  free(b->result_recv);
        free(b->requests);
        free(b);
        CO->host_mpi_bufs = NULL;
    }
    CO->host_scratch_capacity = 0;
}

// Exchange enquiry data: non-blocking sends/recvs, then waitall
static void mpi_exchange_enquiry(struct gpu_domain *GD,
                                  struct inlet_data *data0,
                                  struct inlet_data *data1,
                                  struct culvert_mpi_bufs *bufs) {
    struct culvert_operators *CO = &GD->culvert_ops;
    int nc = CO->num_culverts;
    int myrank = GD->rank;
    MPI_Comm comm = GD->comm;
    MPI_Request *requests = bufs->requests;
    int nreq = 0;

    for (int c = 0; c < nc; c++) {
        struct culvert_indices *ci = &CO->indices[c];
        if (ci->is_local) continue;

        struct inlet_data *data[2] = {&data0[c], &data1[c]};

        for (int inlet = 0; inlet < 2; inlet++) {
            int tag = ci->mpi_tag_base + inlet;

            if (myrank == ci->enquiry_proc[inlet] && myrank != ci->master_proc) {
                // I have this enquiry point, send to master
                bufs->enquiry_send[c][inlet][0] = data[inlet]->enquiry_stage;
                bufs->enquiry_send[c][inlet][1] = data[inlet]->enquiry_xmom;
                bufs->enquiry_send[c][inlet][2] = data[inlet]->enquiry_ymom;
                bufs->enquiry_send[c][inlet][3] = data[inlet]->enquiry_elevation;
                MPI_Isend(bufs->enquiry_send[c][inlet], 4, MPI_DOUBLE,
                          ci->master_proc, tag, comm, &requests[nreq++]);
            }

            if (myrank == ci->master_proc && ci->enquiry_proc[inlet] != myrank) {
                // Master receives from enquiry proc
                MPI_Irecv(bufs->enquiry_recv[c][inlet], 4, MPI_DOUBLE,
                          ci->enquiry_proc[inlet], tag, comm, &requests[nreq++]);
            }
        }
    }

    if (nreq > 0)
        MPI_Waitall(nreq, requests, MPI_STATUSES_IGNORE);

    // Unpack received enquiry data into inlet_data structs
    for (int c = 0; c < nc; c++) {
        struct culvert_indices *ci = &CO->indices[c];
        if (ci->is_local) continue;
        if (myrank != ci->master_proc) continue;

        struct inlet_data *data[2] = {&data0[c], &data1[c]};

        for (int inlet = 0; inlet < 2; inlet++) {
            if (ci->enquiry_proc[inlet] != myrank) {
                data[inlet]->enquiry_stage     = bufs->enquiry_recv[c][inlet][0];
                data[inlet]->enquiry_xmom      = bufs->enquiry_recv[c][inlet][1];
                data[inlet]->enquiry_ymom      = bufs->enquiry_recv[c][inlet][2];
                data[inlet]->enquiry_elevation  = bufs->enquiry_recv[c][inlet][3];
            }
        }
    }
}

// Exchange inlet averages: each inlet_master sends local sums to structure master
static void mpi_exchange_inlet_averages(struct gpu_domain *GD,
                                         struct inlet_data *data0,
                                         struct inlet_data *data1,
                                         struct culvert_mpi_bufs *bufs) {
    struct culvert_operators *CO = &GD->culvert_ops;
    int nc = CO->num_culverts;
    int myrank = GD->rank;
    MPI_Comm comm = GD->comm;
    MPI_Request *requests = bufs->requests;
    int nreq = 0;

    for (int c = 0; c < nc; c++) {
        struct culvert_indices *ci = &CO->indices[c];
        if (ci->is_local) continue;

        struct inlet_data *data[2] = {&data0[c], &data1[c]};

        for (int inlet = 0; inlet < 2; inlet++) {
            int tag = ci->mpi_tag_base + 2 + inlet;

            if (myrank == ci->inlet_master_proc[inlet] && myrank != ci->master_proc) {
                // Send local area-weighted sums to structure master
                double area = data[inlet]->total_area;
                bufs->inlet_send[c][inlet][0] = data[inlet]->avg_stage * area;
                bufs->inlet_send[c][inlet][1] = data[inlet]->avg_depth * area;
                bufs->inlet_send[c][inlet][2] = data[inlet]->avg_xmom * area;
                bufs->inlet_send[c][inlet][3] = data[inlet]->avg_ymom * area;
                bufs->inlet_send[c][inlet][4] = area;
                MPI_Isend(bufs->inlet_send[c][inlet], 5, MPI_DOUBLE,
                          ci->master_proc, tag, comm, &requests[nreq++]);
            }

            if (myrank == ci->master_proc && ci->inlet_master_proc[inlet] != myrank) {
                MPI_Irecv(bufs->inlet_recv[c][inlet], 5, MPI_DOUBLE,
                          ci->inlet_master_proc[inlet], tag, comm, &requests[nreq++]);
            }
        }
    }

    if (nreq > 0)
        MPI_Waitall(nreq, requests, MPI_STATUSES_IGNORE);

    // Master combines local + remote inlet averages
    for (int c = 0; c < nc; c++) {
        struct culvert_indices *ci = &CO->indices[c];
        if (ci->is_local) continue;
        if (myrank != ci->master_proc) continue;

        struct inlet_data *data[2] = {&data0[c], &data1[c]};

        for (int inlet = 0; inlet < 2; inlet++) {
            if (ci->inlet_master_proc[inlet] != myrank) {
                // Replace local (placeholder) with remote data
                double remote_area = bufs->inlet_recv[c][inlet][4];
                if (remote_area > 0.0) {
                    data[inlet]->avg_stage = bufs->inlet_recv[c][inlet][0] / remote_area;
                    data[inlet]->avg_depth = bufs->inlet_recv[c][inlet][1] / remote_area;
                    data[inlet]->avg_xmom  = bufs->inlet_recv[c][inlet][2] / remote_area;
                    data[inlet]->avg_ymom  = bufs->inlet_recv[c][inlet][3] / remote_area;
                    data[inlet]->total_area = remote_area;
                }
            }
            // If master has local inlet data too (inlet_master == master), it already
            // has the correct averages from the GPU gather phase.
        }
    }
}

// Send computed results from master to remote inlet procs for scatter
static void mpi_exchange_results(struct gpu_domain *GD,
                                  struct culvert_transfer *transfers,
                                  struct culvert_mpi_bufs *bufs) {
    struct culvert_operators *CO = &GD->culvert_ops;
    int nc = CO->num_culverts;
    int myrank = GD->rank;
    MPI_Comm comm = GD->comm;
    MPI_Request *requests = bufs->requests;
    int nreq = 0;

    for (int c = 0; c < nc; c++) {
        struct culvert_indices *ci = &CO->indices[c];
        if (ci->is_local) continue;

        struct culvert_transfer *t = &transfers[c];

        for (int inlet = 0; inlet < 2; inlet++) {
            int tag = ci->mpi_tag_base + 4 + inlet;

            if (myrank == ci->master_proc && ci->inlet_master_proc[inlet] != myrank) {
                // Determine if this inlet is inflow or outflow
                if (inlet == t->inflow_idx) {
                    bufs->result_send[c][inlet][0] = t->new_inflow_depth;
                    bufs->result_send[c][inlet][1] = t->new_inflow_xmom;
                    bufs->result_send[c][inlet][2] = t->new_inflow_ymom;
                } else {
                    bufs->result_send[c][inlet][0] = t->new_outflow_depth;
                    bufs->result_send[c][inlet][1] = t->new_outflow_xmom;
                    bufs->result_send[c][inlet][2] = t->new_outflow_ymom;
                }
                MPI_Isend(bufs->result_send[c][inlet], 3, MPI_DOUBLE,
                          ci->inlet_master_proc[inlet], tag, comm, &requests[nreq++]);
            }

            if (myrank == ci->inlet_master_proc[inlet] && ci->master_proc != myrank) {
                MPI_Irecv(bufs->result_recv[c][inlet], 3, MPI_DOUBLE,
                          ci->master_proc, tag, comm, &requests[nreq++]);
            }
        }
    }

    if (nreq > 0)
        MPI_Waitall(nreq, requests, MPI_STATUSES_IGNORE);

    // Non-master inlet procs build transfer structs from received data
    for (int c = 0; c < nc; c++) {
        struct culvert_indices *ci = &CO->indices[c];
        if (ci->is_local) continue;
        if (myrank == ci->master_proc) continue;

        struct culvert_transfer *t = &transfers[c];

        for (int inlet = 0; inlet < 2; inlet++) {
            if (myrank == ci->inlet_master_proc[inlet]) {
                // This rank owns inlet triangles — use received scatter values
                // We don't know the inflow direction, but we know which inlet we own.
                // Set both inflow/outflow fields; scatter uses inflow_idx to pick.
                // Since we have triangles only for this inlet, we just need the values
                // for this inlet's role (inflow or outflow).
                double new_depth = bufs->result_recv[c][inlet][0];
                double new_xmom  = bufs->result_recv[c][inlet][1];
                double new_ymom  = bufs->result_recv[c][inlet][2];

                // Store the values. The scatter function uses inflow_idx to determine
                // which inlet gets inflow vs outflow values. Since the master computed
                // the correct values for each inlet role, we store them in the right slot.
                if (inlet == 0) {
                    t->new_inflow_depth = new_depth;
                    t->new_inflow_xmom  = new_xmom;
                    t->new_inflow_ymom  = new_ymom;
                    t->inflow_idx = 0;  // Treat inlet 0 as "inflow" for scatter
                } else {
                    t->new_outflow_depth = new_depth;
                    t->new_outflow_xmom  = new_xmom;
                    t->new_outflow_ymom  = new_ymom;
                    t->inflow_idx = 0;  // inlet 1 is "outflow" with inflow_idx=0
                }
            }
        }
    }
}

// ============================================================================
// Per-culvert GPU scatter for cross-boundary culverts on non-master ranks
// These ranks only have triangles for one or both inlets (not all data).
// The triangle indices here are per-culvert (not part of the persistently
// mapped flattened array), so this rare MPI path maps the tiny index list
// per call; stage = bed + depth is computed entirely on-device.
// ============================================================================


// ============================================================================
// Main entry point: execute ALL culverts in one batched cycle
// Local culverts use batched gather/scatter (existing fast path).
// Cross-boundary culverts use MPI between gather and scatter phases.
// ============================================================================

// ============================================================================
// Pure per-culvert compute: discharge determination + semi-implicit water
// transfer.  This is THE single implementation of the Boyd/weir update, shared
// by the mode-2 batch (gpu_culverts_apply_all, below) and the mode-1 host path
// (via the Cython wrapper), so the two compute modes agree bit-for-bit.
//
// Assumes a local/master, in-bounds culvert; the caller handles non-master
// parallel skips.  Reads gathered inlet data (data0/data1), reads+updates the
// smoothing state (st), and writes the result (r) and the water transfer (t).
// ============================================================================
void culvert_compute_one(const struct inlet_data *data0,
                         const struct inlet_data *data1,
                         const struct culvert_params *p,
                         struct culvert_state *st,
                         double timestep,
                         struct culvert_result *r,
                         struct culvert_transfer *t) {
    // Reset per-step reporting stats; filled in below where a discharge runs.
    st->report_gain = 0.0;
    st->report_discharge = 0.0;
    st->report_velocity = 0.0;
    st->report_driving_energy = 0.0;
    st->report_delta_total_energy = 0.0;

    double dim = (p->type == CULVERT_TYPE_BOX || p->type == CULVERT_TYPE_WEIR_TRAPEZOID)
                 ? p->height : p->diameter;

    if (dim <= 0.0) {
        // Closed culvert: no discharge. Transfer below still runs (Q=0 => no-op).
        r->Q = 0.0;
        r->barrel_velocity = 0.0;
        r->outlet_culvert_depth = 0.0;
        r->flow_area = 0.00001;
        r->inflow_idx = 0;
    } else {
        // Compute delta_total_energy to determine flow direction
        double delta_total_energy;
        if (p->use_velocity_head) {
            double depth0, vh0, te0, se0;
            double depth1, vh1, te1, se1;
            compute_enquiry_values(data0, p, 0, &depth0, &vh0, &te0, &se0);
            compute_enquiry_values(data1, p, 1, &depth1, &vh1, &te1, &se1);
            delta_total_energy = te0 - te1;
        } else {
            delta_total_energy = data0->enquiry_stage - data1->enquiry_stage;
        }

        // Smooth delta_total_energy
        double ts;
        culvert_smooth_energy(&st->smooth_delta_total_energy,
                              delta_total_energy, timestep,
                              p->smoothing_timescale, &ts);

        // Determine inflow/outflow
        const struct inlet_data *inflow_data, *outflow_data;
        if (st->smooth_delta_total_energy >= 0.0) {
            r->inflow_idx = 0;
            inflow_data = data0;
            outflow_data = data1;
            delta_total_energy = st->smooth_delta_total_energy;
        } else {
            r->inflow_idx = 1;
            inflow_data = data1;
            outflow_data = data0;
            delta_total_energy = -st->smooth_delta_total_energy;
        }

        st->report_delta_total_energy = delta_total_energy;

        double inflow_depth, inflow_vh, inflow_te, inflow_se;
        compute_enquiry_values(inflow_data, p, r->inflow_idx,
                               &inflow_depth, &inflow_vh, &inflow_te, &inflow_se);

        double outflow_depth, outflow_vh, outflow_te, outflow_se;
        compute_enquiry_values(outflow_data, p, 1 - r->inflow_idx,
                               &outflow_depth, &outflow_vh, &outflow_te, &outflow_se);

        if (inflow_depth > 0.01) {
            double driving_energy;
            if (p->use_velocity_head)
                driving_energy = inflow_se;
            else
                driving_energy = inflow_depth;

            double Q, bv, ocd, fa;

            if (p->type == CULVERT_TYPE_BOX) {
                boyd_box_discharge(p, driving_energy, delta_total_energy,
                                   outflow_depth, &Q, &bv, &ocd, &fa);
            } else if (p->type == CULVERT_TYPE_WEIR_TRAPEZOID) {
                weir_orifice_trapezoid_discharge(p, driving_energy,
                                                 delta_total_energy,
                                                 outflow_depth,
                                                 &Q, &bv, &ocd, &fa);
            } else {
                boyd_pipe_discharge(p, inflow_depth, driving_energy,
                                    delta_total_energy, outflow_depth,
                                    &Q, &bv, &ocd, &fa);
            }

            // Apply discharge smoothing
            culvert_smooth_discharge(st->smooth_delta_total_energy,
                                     &st->smooth_Q, Q, fa, ts, &Q, &bv);

            r->Q = Q;
            r->barrel_velocity = bv;
            r->outlet_culvert_depth = ocd;
            r->flow_area = fa;

            // Clamp velocity
            if (r->barrel_velocity > p->max_velocity) {
                r->barrel_velocity = p->max_velocity;
                r->Q = r->flow_area * r->barrel_velocity;
            }

            st->report_driving_energy = driving_energy;
            st->report_velocity = r->barrel_velocity;
        } else {
            r->Q = 0.0;
            r->barrel_velocity = 0.0;
            r->outlet_culvert_depth = 0.0;
            r->flow_area = 0.00001;
        }
    }

    // ---- Semi-implicit water transfer (was PHASE 3) ----
    t->inflow_idx = r->inflow_idx;

    const struct inlet_data *inflow_data = (r->inflow_idx == 0) ? data0 : data1;
    const struct inlet_data *outflow_data = (r->inflow_idx == 0) ? data1 : data0;
    double inflow_area = inflow_data->total_area;
    double outflow_area = outflow_data->total_area;

    double old_inflow_depth = inflow_data->avg_depth;
    double old_inflow_xmom = inflow_data->avg_xmom;
    double old_inflow_ymom = inflow_data->avg_ymom;

    // Semi-implicit factor
    double dt_Q_on_d;
    if (old_inflow_depth > 0.0)
        dt_Q_on_d = timestep * r->Q / old_inflow_depth;
    else
        dt_Q_on_d = 0.0;

    double factor = 1.0 / (1.0 + dt_Q_on_d / inflow_area);

    double new_inflow_depth, timestep_star;
    if (p->always_use_Q_wetdry_adjustment) {
        new_inflow_depth = old_inflow_depth * factor;
        if (old_inflow_depth > 0.0)
            timestep_star = timestep * new_inflow_depth / old_inflow_depth;
        else
            timestep_star = 0.0;
    } else {
        new_inflow_depth = old_inflow_depth - timestep * r->Q / inflow_area;
        timestep_star = timestep;
    }

    st->report_gain = r->Q * timestep_star;
    st->report_discharge = (timestep > 0.0) ? (r->Q * timestep_star / timestep) : 0.0;

    double new_inflow_xmom, new_inflow_ymom;
    if (p->use_old_momentum_method) {
        new_inflow_xmom = old_inflow_xmom * factor;
        new_inflow_ymom = old_inflow_ymom * factor;
    } else {
        double factor2;
        if (old_inflow_depth > 0.0) {
            if (p->always_use_Q_wetdry_adjustment)
                factor2 = 1.0 / (1.0 + dt_Q_on_d * new_inflow_depth / (old_inflow_depth * inflow_area));
            else
                factor2 = 1.0 / (1.0 + timestep * r->Q / (old_inflow_depth * inflow_area));
        } else {
            factor2 = 0.0;
        }
        new_inflow_xmom = old_inflow_xmom * factor2;
        new_inflow_ymom = old_inflow_ymom * factor2;
    }

    t->new_inflow_depth = new_inflow_depth;
    t->new_inflow_xmom = new_inflow_xmom;
    t->new_inflow_ymom = new_inflow_ymom;

    // Outflow
    double outflow_extra_depth = r->Q * timestep_star / outflow_area;
    double new_outflow_depth = outflow_data->avg_depth + outflow_extra_depth;

    const double *outflow_vec = (r->inflow_idx == 0) ? p->outward_vector_1 : p->outward_vector_0;
    double dir0 = -outflow_vec[0];
    double dir1 = -outflow_vec[1];

    double new_outflow_xmom, new_outflow_ymom;
    if (p->use_momentum_jet) {
        new_outflow_xmom = r->barrel_velocity * new_outflow_depth * dir0;
        new_outflow_ymom = r->barrel_velocity * new_outflow_depth * dir1;
    } else {
        new_outflow_xmom = 0.0;
        new_outflow_ymom = 0.0;
    }

    t->new_outflow_depth = new_outflow_depth;
    t->new_outflow_xmom = new_outflow_xmom;
    t->new_outflow_ymom = new_outflow_ymom;
}

// ============================================================================
// Flat host entry point for the mode-1 (Python) path. Marshals scalars into the
// culvert_params / culvert_state / inlet_data structs, calls the shared
// culvert_compute_one(), and returns the transfer + reporting via out-params.
// No gpu_domain, no device memory — operates purely on values gathered by the
// Python operator, so mode-1 and mode-2 run identical culvert arithmetic.
// ============================================================================
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
        double *outlet_culvert_depth) {

    struct culvert_params p;
    memset(&p, 0, sizeof(p));
    p.type = type; p.g = g; p.width = width; p.height = height; p.diameter = diameter;
    p.z1 = z1; p.z2 = z2; p.length = length; p.manning = manning; p.sum_loss = sum_loss;
    p.blockage = blockage; p.barrels = barrels;
    p.use_velocity_head = use_velocity_head;
    p.use_momentum_jet = use_momentum_jet;
    p.use_old_momentum_method = use_old_momentum_method;
    p.always_use_Q_wetdry_adjustment = always_use_Q_wetdry_adjustment;
    p.max_velocity = max_velocity; p.smoothing_timescale = smoothing_timescale;
    p.outward_vector_0[0] = ov0x; p.outward_vector_0[1] = ov0y;
    p.outward_vector_1[0] = ov1x; p.outward_vector_1[1] = ov1y;
    p.invert_elevation_0 = invert0; p.invert_elevation_1 = invert1;
    p.has_invert_elevation_0 = has_invert0; p.has_invert_elevation_1 = has_invert1;

    struct culvert_state st;
    memset(&st, 0, sizeof(st));
    st.smooth_delta_total_energy = *smooth_delta_total_energy;
    st.smooth_Q = *smooth_Q;

    struct inlet_data d0, d1;
    d0.enquiry_stage = e0_stage; d0.enquiry_xmom = e0_xmom;
    d0.enquiry_ymom = e0_ymom; d0.enquiry_elevation = e0_elev;
    d0.avg_stage = a0_stage; d0.avg_depth = a0_depth;
    d0.avg_xmom = a0_xmom; d0.avg_ymom = a0_ymom; d0.total_area = a0_area;
    d1.enquiry_stage = e1_stage; d1.enquiry_xmom = e1_xmom;
    d1.enquiry_ymom = e1_ymom; d1.enquiry_elevation = e1_elev;
    d1.avg_stage = a1_stage; d1.avg_depth = a1_depth;
    d1.avg_xmom = a1_xmom; d1.avg_ymom = a1_ymom; d1.total_area = a1_area;

    struct culvert_result r;
    struct culvert_transfer t;
    memset(&r, 0, sizeof(r));
    memset(&t, 0, sizeof(t));

    culvert_compute_one(&d0, &d1, &p, &st, timestep, &r, &t);

    *smooth_delta_total_energy = st.smooth_delta_total_energy;
    *smooth_Q = st.smooth_Q;
    *inflow_idx = t.inflow_idx;
    *new_inflow_depth = t.new_inflow_depth;
    *new_inflow_xmom = t.new_inflow_xmom;
    *new_inflow_ymom = t.new_inflow_ymom;
    *new_outflow_depth = t.new_outflow_depth;
    *new_outflow_xmom = t.new_outflow_xmom;
    *new_outflow_ymom = t.new_outflow_ymom;
    *report_gain = st.report_gain;
    *report_discharge = st.report_discharge;
    *report_velocity = st.report_velocity;
    *report_driving_energy = st.report_driving_energy;
    *report_delta_total_energy = st.report_delta_total_energy;
    *outlet_culvert_depth = r.outlet_culvert_depth;
}

// Host inlet gather for the mode-1 path. Mirrors the inner reduction of
// gpu_culvert_gather_inlets() exactly (same accumulation order, same depth
// clamp, same divisor) so mode-1's inlet averages match mode-2's bit-for-bit.
void culvert_gather_inlet_host(int n, const int *indices, const double *areas,
                               const double *stage_c, const double *xmom_c,
                               const double *ymom_c, const double *bed_c,
                               double total_area,
                               double *avg_stage, double *avg_depth,
                               double *avg_xmom, double *avg_ymom) {
    double sum_stage = 0.0, sum_depth = 0.0, sum_xmom = 0.0, sum_ymom = 0.0;
    for (int j = 0; j < n; j++) {
        int i = indices[j];
        double area = areas[j];
        double depth = stage_c[i] - bed_c[i];
        if (depth < 0.0) depth = 0.0;
        sum_stage += stage_c[i] * area;
        sum_depth += depth * area;
        sum_xmom += xmom_c[i] * area;
        sum_ymom += ymom_c[i] * area;
    }
    if (total_area > 0.0) {
        *avg_stage = sum_stage / total_area;
        *avg_depth = sum_depth / total_area;
        *avg_xmom = sum_xmom / total_area;
        *avg_ymom = sum_ymom / total_area;
    } else {
        *avg_stage = 0.0; *avg_depth = 0.0; *avg_xmom = 0.0; *avg_ymom = 0.0;
    }
}

// Make sure the host-side working buffers can hold nc culverts, growing them
// if not.  Called once per step; allocation only happens when the culvert
// count first exceeds what is already allocated, so the steady state is
// malloc-free.  Returns 0 on success, -1 if allocation failed.
static int culvert_host_scratch_ensure(struct culvert_operators *CO, int nc) {
    if (CO->host_scratch_capacity >= nc && CO->host_mpi_bufs) {
        return 0;
    }

    struct inlet_data *d0 = (struct inlet_data*)
        realloc(CO->host_data0, nc * sizeof(struct inlet_data));
    if (d0) CO->host_data0 = d0;
    struct inlet_data *d1 = (struct inlet_data*)
        realloc(CO->host_data1, nc * sizeof(struct inlet_data));
    if (d1) CO->host_data1 = d1;
    struct culvert_result *rs = (struct culvert_result*)
        realloc(CO->host_results, nc * sizeof(struct culvert_result));
    if (rs) CO->host_results = rs;
    struct culvert_transfer *tr = (struct culvert_transfer*)
        realloc(CO->host_transfers, nc * sizeof(struct culvert_transfer));
    if (tr) CO->host_transfers = tr;

    struct culvert_mpi_bufs *b = CO->host_mpi_bufs;
    if (!b) {
        b = (struct culvert_mpi_bufs*) calloc(1, sizeof(*b));
        CO->host_mpi_bufs = b;
    }
    if (b) {
        void *es = realloc(b->enquiry_send, nc * sizeof(double[2][4]));
        void *er = realloc(b->enquiry_recv, nc * sizeof(double[2][4]));
        void *is = realloc(b->inlet_send,   nc * sizeof(double[2][5]));
        void *ir = realloc(b->inlet_recv,   nc * sizeof(double[2][5]));
        void *rls = realloc(b->result_send, nc * sizeof(double[2][3]));
        void *rlr = realloc(b->result_recv, nc * sizeof(double[2][3]));
        void *rq = realloc(b->requests, (size_t)nc * 6 * sizeof(MPI_Request));
        if (es) b->enquiry_send = (double(*)[2][4]) es;
        if (er) b->enquiry_recv = (double(*)[2][4]) er;
        if (is) b->inlet_send   = (double(*)[2][5]) is;
        if (ir) b->inlet_recv   = (double(*)[2][5]) ir;
        if (rls) b->result_send = (double(*)[2][3]) rls;
        if (rlr) b->result_recv = (double(*)[2][3]) rlr;
        if (rq) b->requests     = (MPI_Request*) rq;
        if (es && er && is && ir && rls && rlr && rq) {
            b->capacity = nc;
            b->nreq_capacity = nc * 6;
        }
    }

    if (!d0 || !d1 || !rs || !tr || !b || b->capacity < nc) {
        fprintf(stderr, "ERROR: failed to allocate host scratch for %d culverts\n", nc);
        return -1;
    }

    CO->host_scratch_capacity = nc;
    return 0;
}

void gpu_culverts_apply_all(struct gpu_domain *GD, double timestep) {
    NVTX_PUSH("gpu_culverts_apply_all");
    struct culvert_operators *CO = &GD->culvert_ops;
    int nc = CO->num_culverts;

    if (nc == 0 || !CO->initialized) {
        NVTX_POP();
        return;
    }

    omp_set_default_device(gpu_compute_device(GD));
    int myrank = GD->rank;

    // Check if any parallel culverts exist
    int has_parallel = 0;
    for (int c = 0; c < nc; c++) {
        if (!CO->indices[c].is_local) { has_parallel = 1; break; }
    }

    // Per-culvert working data, sized to the actual culvert count.
    if (culvert_host_scratch_ensure(CO, nc) != 0) {
        gpu_set_error(GD, "host scratch for %d culverts could not be allocated; culverts not applied this step", nc);
        NVTX_POP();
        return;
    }
    struct inlet_data *data0 = CO->host_data0;
    struct inlet_data *data1 = CO->host_data1;
    struct culvert_result *results = CO->host_results;
    struct culvert_transfer *transfers = CO->host_transfers;

    // ----------------------------------------------------------------
    // PHASE 1: Batched GPU gather (2 target update from's)
    // Gathers LOCAL data for ALL culverts (local + parallel).
    // Remote enquiry points get placeholder values (overwritten by MPI).
    // ----------------------------------------------------------------
    gpu_culvert_gather_enquiry(GD, data0, data1);
    gpu_culvert_gather_inlets(GD, data0, data1);

    // ----------------------------------------------------------------
    // PHASE 1b: MPI exchange for cross-boundary culverts
    // ----------------------------------------------------------------
    // MPI buffers live with the domain (grown above), not in a function-local
    // static: a static would be shared by every domain in the process and
    // fixed at MAX_CULVERTS entries.
    struct culvert_mpi_bufs *mpi_bufs = CO->host_mpi_bufs;

    if (has_parallel) {
        mpi_exchange_enquiry(GD, data0, data1, mpi_bufs);
        mpi_exchange_inlet_averages(GD, data0, data1, mpi_bufs);
    }

    // ----------------------------------------------------------------
    // PHASE 2+3: per-culvert discharge + semi-implicit water transfer.
    // The actual physics lives in culvert_compute_one() (above), which mode-1
    // also calls via Cython so both compute modes are bit-for-bit identical.
    // Only master_proc computes for cross-boundary culverts.
    // ----------------------------------------------------------------
    for (int c = 0; c < nc; c++) {
        struct culvert_indices *ci = &CO->indices[c];
        struct culvert_params *p = &CO->params[c];
        struct culvert_state *st = &CO->state[c];
        struct culvert_result *r = &results[c];
        struct culvert_transfer *t = &transfers[c];

        // Non-master ranks skip: results/transfer arrive via MPI below.
        if (!ci->is_local && myrank != ci->master_proc) {
            st->report_gain = 0.0;
            st->report_discharge = 0.0;
            st->report_velocity = 0.0;
            st->report_driving_energy = 0.0;
            st->report_delta_total_energy = 0.0;
            r->Q = 0.0;
            r->barrel_velocity = 0.0;
            r->outlet_culvert_depth = 0.0;
            r->flow_area = 0.00001;
            r->inflow_idx = 0;
            memset(t, 0, sizeof(*t));
            continue;
        }

        culvert_compute_one(&data0[c], &data1[c], p, st, timestep, r, t);
    }

    // ----------------------------------------------------------------
    // PHASE 3b: MPI exchange results for cross-boundary culverts
    // Master sends scatter values to remote inlet procs.
    // ----------------------------------------------------------------
    if (has_parallel) {
        mpi_exchange_results(GD, transfers, mpi_bufs);

        // No scatter here: a non-master rank that owns an inlet received the
        // master's values into its transfers[c] above, and PHASE 4's batched
        // scatter writes every slot this rank holds triangles for — including
        // that one. Scattering here as well used to be harmless because the
        // write was an idempotent "stage = bed + depth"; now that the write is
        // a surface SHIFT (issue #229), applying it twice would double it.
    }

    // ----------------------------------------------------------------
    // PHASE 4: Batched GPU scatter for local culverts (and master's
    // local inlets for parallel culverts)
    // ----------------------------------------------------------------
    gpu_culvert_scatter(GD, transfers);
    NVTX_POP();
}

// Read back a culvert's per-step reporting stats. See gpu_domain.h.
int gpu_culverts_get_report(struct gpu_domain *GD, int culvert_id, double *out) {
    struct culvert_operators *CO = &GD->culvert_ops;
    if (culvert_id < 0 || culvert_id >= CO->num_culverts) {
        return -1;
    }
    struct culvert_state *st = &CO->state[culvert_id];
    out[0] = st->report_gain;
    out[1] = st->report_discharge;
    out[2] = st->report_velocity;
    out[3] = st->report_driving_energy;
    out[4] = st->report_delta_total_energy;
    return 0;
}
