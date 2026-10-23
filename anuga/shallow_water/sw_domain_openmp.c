// Python - C extension module for shallow_water.py
//
// To compile (Python2.6):
//  gcc -c swb2_domain_ext.c -I/usr/include/python2.6 -o domain_ext.o -Wall -O
//  gcc -shared swb2_domain_ext.o  -o swb2_domain_ext.so
//
// or use python compile.py
//
// See the module swb_domain.py for more documentation on
// how to use this module
//
//
// Ole Nielsen, GA 2004
// Stephen Roberts, ANU 2009
// Gareth Davies, GA 2011
// Jorge Galvez, NCI 2026 :D

#include "math.h"
#include <math.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <stdint.h>

#include "sw_domain_math.h"
#include "util_ext.h"
#include "sw_domain.h"
#include "anuga_constants.h"

// Core kernels for CPU/GPU (shared with sw_domain_gpu_ext)
#include "core_kernels.h"

// Shared device helper functions (flux, rotation, gradient limiting, etc.)
#include "gpu/gpu_device_helpers.h"

// FIXME: Perhaps use the epsilon used elsewhere.

// Trick to compute n modulo d (n%d in python) when d is a power of 2
anuga_uint __mod_of_power_2(anuga_uint n, anuga_uint d)
{
  return (n & (d - 1));
}

// Wrapper for scalar interface (used by Cython)
// Calls unified gpu_flux_function_central from gpu_device_helpers.h
anuga_int __openmp__flux_function_central(double q_left0, double q_left1, double q_left2,
                                        double q_right0, double q_right1, double q_right2,
                                        double h_left, double h_right,
                                        double hle, double hre,
                                        double n1, double n2,
                                        double epsilon,
                                        double ze,
                                        double g,
                                        double *edgeflux0, double *edgeflux1, double *edgeflux2,
                                        double *max_speed,
                                        double *pressure_flux,
                                        anuga_int low_froude)
{
  double edgeflux[3];
  double q_left[3];
  double q_right[3];

  edgeflux[0] = *edgeflux0;
  edgeflux[1] = *edgeflux1;
  edgeflux[2] = *edgeflux2;

  q_left[0] = q_left0;
  q_left[1] = q_left1;
  q_left[2] = q_left2;

  q_right[0] = q_right0;
  q_right[1] = q_right1;
  q_right[2] = q_right2;

  gpu_flux_function_central(q_left, q_right,
                            h_left, h_right,
                            hle, hre,
                            n1, n2,
                            epsilon, ze, g,
                            edgeflux, max_speed,
                            pressure_flux,
                            low_froude);

  *edgeflux0 = edgeflux[0];
  *edgeflux1 = edgeflux[1];
  *edgeflux2 = edgeflux[2];

  return 0;
}

// Wrapper for Cython - calls gpu_rotate from gpu_device_helpers.h
anuga_int __rotate(double *q, double n1, double n2)
{
  gpu_rotate(q, n1, n2);
  return 0;
}

// Unified: calls core_compute_fluxes_central from core_kernels.c
// Substep tracking lives in the domain struct (flux_call_count,
// flux_base_call, flux_timestep_fluxcalls_seen) so it is per domain rather
// than per module: two domains in one process no longer share counters.
double _openmp_compute_fluxes_central(const struct domain *__restrict D,
                                      double timestep)
{
  struct domain *Dm = (struct domain *)D;

  Dm->flux_call_count++;

  if (D->timestep_fluxcalls != D->flux_timestep_fluxcalls_seen) {
    Dm->flux_timestep_fluxcalls_seen = D->timestep_fluxcalls;
    Dm->flux_base_call = D->flux_call_count;
  }

  // Which substep of the timestepping method are we on?
  int substep_count = (D->flux_call_count - D->flux_base_call) % D->timestep_fluxcalls;

  // Call unified flux computation
  double local_timestep = core_compute_fluxes_central((struct domain *)D, substep_count, D->timestep_fluxcalls);

  // Return timestep only on first substep
  if (substep_count == 0) {
    timestep = local_timestep;
  }

  return timestep;
}

// Protect against the water elevation falling below the triangle bed
// Unified: calls core_protect from core_kernels.c
double _openmp_protect(const struct domain *__restrict D)
{
    return core_protect((struct domain *)D);
}

// Extrapolation helpers now unified in gpu_device_helpers.h:
// - gpu_find_qmin_and_qmax_dq1_dq2
// - gpu_limit_gradient
// - gpu_calc_edge_values_with_gradient
// - gpu_set_constant_edge_values
// - gpu_compute_qmin_qmax_from_dq1

// Extrapolate second order edge values from centroid values
// Unified: calls core_extrapolate_second_order_edge from core_kernels.c
void _openmp_extrapolate_second_order_edge_sw(struct domain *__restrict D)
{
    core_extrapolate_second_order_edge(D);
}

void _openmp_distribute_edges_to_vertices(struct domain *__restrict D)
{
    // Unified: calls core_distribute_edges_to_vertices from core_kernels.c
    core_distribute_edges_to_vertices(D);
}

void _openmp_ader_ck_predictor(struct domain *__restrict D, double dt)
{
    // ADER Cauchy-Kovalewski predictor — calls core_ader_ck_predictor from core_kernels.c
    core_ader_ck_predictor(D, dt);
}

void _openmp_ader_ck_predictor_edge(struct domain *__restrict D, double dt)
{
    // Fused ADER-2 predictor: writes Q^{n+1/2} to edge values, centroids untouched
    core_ader_ck_predictor_edge(D, dt);
}

void _openmp_manning_friction_flat_semi_implicit(const struct domain *__restrict D)
{
    // Unified: calls core_manning_friction_flat_semi_implicit from core_kernels.c
    core_manning_friction_flat_semi_implicit((struct domain *)D);
}




    

void _openmp_vegetation_friction_semi_implicit(const struct domain *__restrict D)
{
    // Unified: calls core_vegetation_friction_semi_implicit from core_kernels.c
    core_vegetation_friction_semi_implicit((struct domain *)D);
}

void _openmp_manning_friction_sloped_semi_implicit(const struct domain *__restrict D)
{
    // Unified: calls core_manning_friction_sloped_semi_implicit from core_kernels.c
    core_manning_friction_sloped_semi_implicit((struct domain *)D);
}

// Unified: calls core_manning_friction_sloped_semi_implicit_edge_based from core_kernels.c
// Edge-based variant of sloped Manning friction — the active per-timestep path
// when domain.use_sloped_mannings=True (friction.py:29 selects this function).
void _openmp_manning_friction_sloped_semi_implicit_edge_based(const struct domain *__restrict D)
{
    core_manning_friction_sloped_semi_implicit_edge_based((struct domain *)D);
}

// Legacy explicit (non semi-implicit) Manning friction variants.
// These are NOT reachable from any non-test production code path:
//   - friction.py only calls the semi-implicit variants above.
//   - The functions below are exposed through sw_domain_openmp_ext.pyx purely
//     for backwards compatibility and unit tests in test_friction.py.
// Do NOT GPU-port these; keep their bodies here for test/back-compat.
void _openmp_manning_friction_flat(const double g, const double eps, const anuga_int N,
                                   double *__restrict w, double *__restrict z_centroid,
                                   double *__restrict uh, double *__restrict vh,
                                   double *__restrict eta, double *__restrict xmom_update, double *__restrict ymom_update)
{

  anuga_int k;
  const double seven_thirds = ANUGA_SEVEN_THIRDS;

#pragma omp parallel for schedule(static) firstprivate(eps, g, seven_thirds)
  for (k = 0; k < N; k++)
  {
    double S, h, z, abs_mom;
    abs_mom = sqrt((uh[k] * uh[k] + vh[k] * vh[k]));
    S = 0.0;

    if (eta[k] > eps)
    {
      z = z_centroid[k];
      h = w[k] - z;
      if (h >= eps)
      {
        S = -g * eta[k] * eta[k] * abs_mom;
        S /= pow(h, seven_thirds); 
      }
    }
    xmom_update[k] += S * uh[k];
    ymom_update[k] += S * vh[k];
  }
}


// Legacy explicit (non semi-implicit) sloped friction — see note above.
void _openmp_manning_friction_sloped(const double g, const double eps, const anuga_int N,
                                     double *__restrict x_vertex, double *__restrict w, double *__restrict z_vertex,
                                     double *__restrict uh, double *__restrict vh,
                                     double *__restrict eta, double *__restrict xmom_update, double *__restrict ymom_update)
{

  const double one_third = ANUGA_ONE_THIRD;
  const double seven_thirds = ANUGA_SEVEN_THIRDS;

#pragma omp parallel for schedule(static) firstprivate(eps, g, one_third, seven_thirds)
  for (anuga_int k = 0; k < N; k++)
  {
    double S = 0.0;
    anuga_int k3 = 3 * k;
    // Get bathymetry
    double z0 = z_vertex[k3 + 0];
    double z1 = z_vertex[k3 + 1];
    double z2 = z_vertex[k3 + 2];

    // Compute bed slope
    anuga_int k6 = 6 * k; // base index

    double x0 = x_vertex[k6 + 0];
    double y0 = x_vertex[k6 + 1];
    double x1 = x_vertex[k6 + 2];
    double y1 = x_vertex[k6 + 3];
    double x2 = x_vertex[k6 + 4];
    double y2 = x_vertex[k6 + 5];

    if (eta[k] > eps)
    {
      double zx, zy, zs, z, h;
      _gradient(x0, y0, x1, y1, x2, y2, z0, z1, z2, &zx, &zy);

      zs = sqrt(1.0 + zx * zx + zy * zy);
      z = (z0 + z1 + z2) * one_third;
      h = w[k] - z;
      if (h >= eps)
      {
        S = -g * eta[k] * eta[k] * zs * sqrt((uh[k] * uh[k] + vh[k] * vh[k]));
        S /= pow(h, seven_thirds); 
      }
    }
    xmom_update[k] += S * uh[k];
    ymom_update[k] += S * vh[k];
  }
}

// Legacy explicit (non semi-implicit) sloped edge-based friction — see note above.
void _openmp_manning_friction_sloped_edge_based(const double g, const double eps, const anuga_int N,
                                     double *__restrict x_edge, double *__restrict w, double *__restrict z_edge,
                                     double *__restrict uh, double *__restrict vh,
                                     double *__restrict eta, double *__restrict xmom_update, double *__restrict ymom_update)
{

  const double one_third = ANUGA_ONE_THIRD;
  const double seven_thirds = ANUGA_SEVEN_THIRDS;

#pragma omp parallel for schedule(static) firstprivate(eps, g, one_third, seven_thirds)
  for (anuga_int k = 0; k < N; k++)
  {
    double S = 0.0;
    anuga_int k3 = 3 * k;
    // Get bathymetry
    double z0 = z_edge[k3 + 0];
    double z1 = z_edge[k3 + 1];
    double z2 = z_edge[k3 + 2];

    // Compute bed slope
    anuga_int k6 = 6 * k; // base index

    double x0 = x_edge[k6 + 0];
    double y0 = x_edge[k6 + 1];
    double x1 = x_edge[k6 + 2];
    double y1 = x_edge[k6 + 3];
    double x2 = x_edge[k6 + 4];
    double y2 = x_edge[k6 + 5];

    if (eta[k] > eps)
    {
      double zx, zy, zs, z, h;
      _gradient(x0, y0, x1, y1, x2, y2, z0, z1, z2, &zx, &zy);

      zs = sqrt(1.0 + zx * zx + zy * zy);
      z = (z0 + z1 + z2) * one_third;
      h = w[k] - z;
      if (h >= eps)
      {
        S = -g * eta[k] * eta[k] * zs * sqrt((uh[k] * uh[k] + vh[k] * vh[k]));
        S /= pow(h, seven_thirds); 
      }
    }
    xmom_update[k] += S * uh[k];
    ymom_update[k] += S * vh[k];
  }
}


// Unified: calls core_fix_negative_cells from core_kernels.c
// core version matches this reference implementation exactly (same & condition,
// same tri_full_flag guard, same reduction pattern).
anuga_int _openmp_fix_negative_cells(const struct domain *__restrict D)
{
    return (anuga_int)core_fix_negative_cells((struct domain *)D);
}


// Read-only: volume that _openmp_fix_negative_cells would add by clamping
// negative-depth cells to zero depth. Call after the flux update but before
// _openmp_fix_negative_cells.
double _openmp_negative_cells_volume(const struct domain *__restrict D)
{
    return core_negative_cells_volume((struct domain *)D);
}


// Unified: calls core_gravity from core_kernels.c
// core_gravity computes avg_h = stage_cv - bed_cv directly (it no longer
// reads height_cv, which could be stale when gravity is used as a forcing
// term), matching this function's original behaviour.
anuga_int _openmp_gravity(const struct domain *__restrict D) {
    return (anuga_int)core_gravity((struct domain *)D);
}

// Unified: calls core_gravity_wb from core_kernels.c
// The well-balanced implementation (stage-gradient + edge-pressure terms)
// was ported there from the serial body that used to live here.
anuga_int _openmp_gravity_wb(const struct domain *__restrict D) {
    return (anuga_int)core_gravity_wb((struct domain *)D);
}



anuga_int _openmp_update_conserved_quantities(const struct domain *__restrict D,
                                              const double timestep)
{
    // Unified: calls core_update_conserved_quantities from core_kernels.c
    core_update_conserved_quantities((struct domain *)D, timestep);
    return 0;
}

anuga_int _openmp_saxpy_conserved_quantities(const struct domain *__restrict D,
                                             const double a,
                                             const double b,
                                             const double c)
{
    // Unified: calls core_saxpy_conserved_quantities from core_kernels.c
    // Computes Q = (a*Q + b*Q_backup) / c
    core_saxpy_conserved_quantities((struct domain *)D, a, b, c);
    return 0;
}

anuga_int _openmp_backup_conserved_quantities(const struct domain *__restrict D)
{
    // Unified: calls core_backup_conserved_quantities from core_kernels.c
    core_backup_conserved_quantities((struct domain *)D);
    return 0;
}

void _openmp_set_omp_num_threads(anuga_int num_threads)
{
  // Set the number of threads for OpenMP
  // This is a global setting and will affect all subsequent OpenMP parallel regions
  omp_set_num_threads(num_threads);
}

// PORTED TO GPU: See gpu_evaluate_reflective_boundary() in gpu_boundaries.c
// TODO: NOT UNIFIED - GPU and CPU use different architectural approaches:
//   - GPU: Pre-collects ALL boundary info during init, evaluates all edges in one kernel call
//   - CPU: Python calls this C function per boundary segment/tag, passing arrays each time
// Unifying would require changing how mode 1 handles boundaries (Python-side changes),
// which could break existing functionality. Keep separate implementations for now.
void _openmp_evaluate_reflective_segment(const struct domain *__restrict D, anuga_int N,
   anuga_int *edge_segment, anuga_int *vol_ids, anuga_int *edge_ids){

  anuga_int number_of_edges = N;

    #pragma omp parallel for schedule(static)
     for(anuga_int k = 0; k < number_of_edges; k++){


      // get vol_ids 
      anuga_int edge_segment_id = edge_segment[k];
      anuga_int vid = vol_ids[k];
      anuga_int edge_id = edge_ids[k];
      double n1 = D->normals[vid * 6 + 2 * edge_id];
      double n2 = D->normals[vid * 6 + 2 * edge_id + 1];

      D->stage_boundary_values[edge_segment_id] = D->stage_edge_values[3 * vid + edge_id];
      // the bed is the elevation
      D->bed_boundary_values[edge_segment_id] = D->bed_edge_values[3 * vid + edge_id];
      D->height_boundary_values[edge_segment_id] = D->height_edge_values[3 * vid + edge_id];

      double q1 = D->xmom_edge_values[3 * vid + edge_id];
      double q2 = D->ymom_edge_values[3 * vid + edge_id];

      double r1 = -q1*n1 - q2*n2;
      double r2 = -q1*n2 + q2*n1;

      double x_mom_boundary_value = n1*r1 - n2*r2;
      double y_mom_boundary_value = n2*r1 + n1*r2;

      D->xmom_boundary_values[edge_segment_id] = x_mom_boundary_value;
      D->ymom_boundary_values[edge_segment_id] = y_mom_boundary_value;

	  // FIXME SR: Check that we really need to work with velocities
      q1 = D->xvelocity_edge_values[3 * vid + edge_id];
      q2 = D->yvelocity_edge_values[3 * vid + edge_id];

      r1 = q1*n1 + q2*n2;
      r2 = q1*n2 - q2*n1;

      double x_vel_boundary_value = n1*r1 - n2*r2;
      double y_vel_boundary_value = n2*r1 + n1*r2;

      D->xvelocity_boundary_values[edge_segment_id] = x_vel_boundary_value;
      D->yvelocity_boundary_values[edge_segment_id] = y_vel_boundary_value;

     }


}
