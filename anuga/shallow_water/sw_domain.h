// C struct for domain and quantities
//
// Stephen Roberts 2012



#ifndef SW_DOMAIN_H
#define SW_DOMAIN_H

#include <stdint.h>
#include <math.h>
#include <stdio.h>
#include <stdbool.h>
#include <inttypes.h>
#include "anuga_typedefs.h"

// Precision of the STATIC GEOMETRY arrays (normals, edgelengths, radii,
// areas, centroid/edge/vertex coordinates).  Default double: bit-compatible
// with every existing consumer (ANUGA never defines the macro).  Building
// with -DANUGA_GEOM_FP32 halves the geometry traffic of the memory-bound
// kernels -- an EXPERIMENT: geometry enters the scheme multiplicatively so
// relative 1e-7 perturbations are benign for convergence, but results are
// not bit-compatible with the double build and the well-balance gate must be
// re-verified.  State (stage/bed/momenta) stays double unconditionally: the
// stage - bed cancellation that well-balancing rests on has ~mm significance
// against ~100 m magnitudes and does not survive FP32.
#ifdef ANUGA_GEOM_FP32
typedef float anuga_geom_t;
#else
typedef double anuga_geom_t;
#endif

// structures
struct domain {
    // Changing these don't change the data in python object
    anuga_int number_of_elements;
    anuga_int boundary_length;
    anuga_int number_of_riverwall_edges;
    anuga_int optimise_dry_cells;
    anuga_int extrapolate_velocity_second_order;
    anuga_int low_froude;
    anuga_int timestep_fluxcalls;
    /* Substep bookkeeping for _openmp_compute_fluxes_central. Per domain,
       so two domains evolving in one process (or one process, two threads)
       do not interleave counters. Initialised once in Domain_C_struct and
       never refreshed from Python. */
    anuga_int flux_call_count;
    anuga_int flux_base_call;
    anuga_int flux_timestep_fluxcalls_seen;
    anuga_int ncol_riverwall_hydraulic_properties;
    anuga_int nrow_riverwall_hydraulic_properties;

    double epsilon;
    double H0;
    double g;
    double evolve_max_timestep;
    double evolve_min_timestep;
    double minimum_allowed_height;
    double maximum_allowed_speed;
    double beta_w;
    double beta_w_dry;
    double beta_uh;
    double beta_uh_dry;
    double beta_vh;
    double beta_vh_dry;

 

    // Changing values in these arrays will change the values in the python object
    anuga_int*   neighbours;
    anuga_int*   neighbour_edges;
    anuga_int*   surrogate_neighbours;
    anuga_geom_t* normals;
    anuga_geom_t* edgelengths;
    anuga_geom_t* radii;
    anuga_geom_t* areas;

    anuga_int*   edge_flux_type;

    anuga_int*   tri_full_flag;
    anuga_int*   already_computed_flux;
    double*    max_speed;

    anuga_geom_t* vertex_coordinates;
    anuga_geom_t* edge_coordinates;
    anuga_geom_t* centroid_coordinates;

    anuga_int* number_of_boundaries;
    double* stage_edge_values;
    double* xmom_edge_values;
    double* ymom_edge_values;
    double* bed_edge_values;
    double* height_edge_values;
    double* xvelocity_edge_values;
    double* yvelocity_edge_values;

    double* stage_centroid_values;
    double* xmom_centroid_values;
    double* ymom_centroid_values;
    double* bed_centroid_values;
    double* height_centroid_values;

    double* stage_vertex_values;
    double* xmom_vertex_values;
    double* ymom_vertex_values;
    double* bed_vertex_values;
    double* height_vertex_values;


    double* stage_boundary_values;
    double* xmom_boundary_values;
    double* ymom_boundary_values;
    double* bed_boundary_values;
    double* height_boundary_values;
    double* xvelocity_boundary_values;
    double* yvelocity_boundary_values;

    double* stage_explicit_update;
    double* xmom_explicit_update;
    double* ymom_explicit_update;

    double* edge_flux_work;
    double* neigh_work;
    double* pressuregrad_work;
    double* x_centroid_work;
    double* y_centroid_work;
    double* boundary_flux_sum;

    anuga_int* edge_river_wall_counter;
    double* riverwall_elevation;
    anuga_int* riverwall_rowIndex;
    double* riverwall_hydraulic_properties;

    double* stage_semi_implicit_update;
    double* xmom_semi_implicit_update;
    double* ymom_semi_implicit_update; 
    
    double* friction_centroid_values;

    /* Vegetation drag, spec 8 (Baptist et al. 2007), alongside Manning.
     *   vegetation_mode 0 = off, 1 = baptist
     * The three fields are the centroid values of the veg_density (stems per
     * m^2), veg_diameter (m) and veg_height (m) quantities; NULL when off. */
    anuga_int vegetation_mode;

    /* Manning friction variant, mirroring friction.py's dispatch:
     *   0 = flat  (the default): no bed-slope factor
     *   1 = sloped: the drag carries zs = sqrt(1 + zx^2 + zy^2)
     * The ADER-2 predictor drag reads this so that it applies exactly the
     * same coefficient the corrector will, which the steady-state balance
     * derived above core_ader_friction_factor() depends on: on a 1:10 bed
     * a zs the corrector does not have is a 0.5% error in gamma. */
    anuga_int use_sloped_mannings;
    double vegetation_Cd;            /* stem drag coefficient, 1.68 */
    double vegetation_bed_chezy;     /* Cb, 65 */
    double* veg_density_centroid_values;
    double* veg_diameter_centroid_values;
    double* veg_height_centroid_values;

    double* stage_backup_values;
    double* xmom_backup_values;
    double* ymom_backup_values;

    /* ------------------------------------------------------------------
     * Generic passive tracers (sediment concentration, salinity, ...).
     *
     * Appended at the END of the struct on purpose so every pre-existing
     * field keeps its offset and the cache layout of the hot arrays is
     * untouched. Adding a member anywhere else silently shifts the offset of
     * everything after it, and the Cython extensions are not rebuilt when
     * this header changes -- so the failure is a wrong binary, not a
     * compile error.
     *
     * number_of_tracers == 0 is the ordinary case and must cost nothing:
     * the flux kernel guards all tracer work behind a single test on this
     * loop-invariant integer.
     *
     * Layout is tracer-major:
     *   centroid[s*n + k]          n = number_of_elements
     *   edge    [s*3n + 3k + i]
     *   boundary[s*boundary_length + m]
     *   explicit_update[s*n + k]
     * ------------------------------------------------------------------ */
    anuga_int number_of_tracers;
    double* tracer_centroid_values;
    double* tracer_edge_values;
    double* tracer_boundary_values;
    double* tracer_explicit_update;
    /* m = h*c, the CONSERVED tracer variable. Integrated by
     * update/backup/saxpy exactly like stage; tracer_centroid_values (c) is
     * DERIVED from it each substep, exactly as height is derived from stage. */
    double* tracer_conserved_values;
    double* tracer_backup_values;
    /* Per-cell accumulation of tracer flux across DOMAIN BOUNDARY edges, for
     * conservation accounting; (n_tracers, N), same layout as the rest.
     *
     * Per-cell rather than a scalar reduction on purpose. The water equivalent
     * uses reduction(+:boundary_flux_sum_substep), but the tracer version would
     * need a runtime-length ARRAY reduction, reduction(+:t[0:ns]), which is
     * unproven on nvc GPU targets. Each cell only ever writes its own edges, so
     * accumulating per cell needs no reduction and no atomics, and the sum is
     * taken on the host at yield time.
     *
     * Sign follows edgeflux[0]: positive is INTO the domain. */
    double* tracer_boundary_flux;
    /* [D-5v] Per-tracer, per-cell factor on the advection speed, (n_tracers,
     * N), or NULL (every tracer moves with the depth-averaged flow). Filled
     * by the sediment source kernel from the log-law velocity of each layer
     * when a velocity profile is on; the flux kernel multiplies the donor
     * cell's factor into the tracer flux, which stays conservative because
     * both cells of an edge see the same product. */
    double* tracer_speed_factor;
    /* Per-substep, per-tracer totals of the above; (timestep_fluxcalls,
     * n_tracers), mirroring boundary_flux_sum for water. Filled by ns SCALAR
     * reductions after the flux loop -- one per tracer, which every target
     * supports -- rather than one array reduction inside it. */
    double* tracer_boundary_flux_sum;
    /* Reconstruction aggressiveness for tracers, analogous to beta_w for stage.
     * 0.0 => first order; >0 => limited second order. Appended at the very end
     * so no previously-existing field offset moves. */
    double beta_tracer;

    /* ------------------------------------------------------------------
     * Phase 3: suspended sediment source terms.
     *
     * Appended at the very END of the struct, after the tracer block, so no
     * previously-existing field offset moves. See HANDOVER.md 2.1: a shifted
     * offset here aliases members silently, with no compile error.
     *
     * All arrays are (n_sediment_classes,) and indexed by TRACER index, so
     * grain size s occupies tracer slot s. Grain sizes are registered through
     * Domain.add_sediment_fraction(), which registers the tracer first.
     *
     * n_sediment_classes == 0 is the ordinary case and must cost nothing: the
     * source kernel returns immediately on a single test, exactly as the flux
     * kernel does for number_of_tracers.
     * ------------------------------------------------------------------ */
    anuga_int n_sediment_classes;
    /* Settling velocity v_s, [S-1] Ferguson & Church (2004), precomputed on
     * the host: it depends only on grain size and fluid properties, so there
     * is no reason to evaluate a square root per cell per step. [m/s] */
    double* sediment_settling_velocity;
    /* d*(Z): ratio of near-bed to depth-averaged concentration, [D-1].
     * 1.0 is the well-mixed limit. Phase 3b replaces this constant with the
     * Rouse profile of spec 4.3. */
    double* sediment_d_star;
    /* [L-2] maximum volumetric concentration. FG21 use 0.30, aS16 0.20. */
    double sediment_c_max;

    /* Phase 3b: entrainment, non-cohesive route [E-1]/[E-2].
     *
     *   tau_star  = f_c |v|^2 / (R g d)      [T-3], f_c = g n^2 h^(-1/3)
     *   S         = tau_star / tau_c_star - 1
     *   E_star    = 0.65 gamma0 S / (1 + gamma0 S)  for S > 0, else 0  [E-1]
     *   E         = v_s E_star                                        [E-2]
     *
     * (Written with _star spelled out: the natural notation tau-star-slash
     * would close this comment.)
     *
     * Water density cancels in tau*, so only the SUBMERGED specific gravity
     * R = rho_s/rho - 1 is needed, not the two densities separately. 0.65 is
     * the maximum packing fraction, so E* saturates rather than growing without
     * bound. */
    double* sediment_diameter;      /* (ncl) grain diameter d_g   [m] */
    double* sediment_R;             /* (ncl) submerged specific gravity */
    double* sediment_tau_c_star;    /* (ncl) critical Shields stress; FG21 0.04 */
    double sediment_gamma0;         /* [E-1] empirical, FG21 0.0024 */

    /* Erosion route, spec 4.1.1. This is a statement about the BED MATERIAL,
     * not a numerical preference: [E-1] and [E-3] describe different sediment
     * and getting it wrong is a physics error, not a tuning error.
     *   0 = NON-COHESIVE, Shields/Smith-McLean [E-1]/[E-2]: sand, gravel
     *   1 = COHESIVE, Hanson & Simon [E-3]/[E-5]: silt, clay, cohesive banks
     *
     *   [E-3]  E   = K_e (tau_b - tau_c)     tau DIMENSIONAL, Pa
     *   [E-5]  K_e = 0.2e-6 / sqrt(tau_c)    [m3 N-1 s-1]
     *
     * Units check (spec open item E1): [m3 N-1 s-1][N m-2] = m/s, the same
     * units as the non-cohesive v_s E*, so it drops into the same slot. */
    anuga_int sediment_erosion_mode;      /* 0 [E-1], 1 [E-3], 2 [E-4] */
    /* [E-4] Partheniades, as RDy26 use it:
     *     E = K_p (tau_b - tau_c)/tau_c   for tau_b > tau_c, else 0
     * K_p is a MASS flux in kg m-2 s-1 (RDy26 give 1e-4), so it is divided by
     * the class density to reach the m/s that E and D use here. */
    double sediment_K_partheniades;
    /* Deposition law:
     *   0 = [D-1]  D = d*(Z) c v_s                          (default)
     *   1 = [D-2]  D = v_s c (1 - tau_b/tau_d),  tau_b < tau_d, else 0
     * [D-2] is RDy26's form. tau_d = 0 disables deposition entirely, which is
     * the hook their passive-transport benchmarks rely on. */
    anuga_int sediment_deposition_mode;
    double sediment_tau_d;                /* critical DEPOSITION stress [Pa] */
    double sediment_tau_crit;       /* tau_c, DIMENSIONAL [Pa] */
    double sediment_K_e;            /* [E-5], [m3 N-1 s-1] */
    double sediment_rho_w;          /* needed for tau_b = rho f_c |v|^2 [T-1] */

    /* Bed shear closure, spec 3.1 and 3.4. Divergence D1.
     *   0 = [T-1] quadratic drag,  tau_b = rho f_c |v|^2   (default)
     *   1 = [T-7] depth-slope,     tau_b = rho g h S       (aSM16 Eqs 6-7)
     *
     * [T-7] is the steady uniform (normal) flow approximation: it assumes the
     * energy slope equals the BED slope and that the flow is locally in
     * equilibrium -- exactly what does not hold in the dam-breach and outburst
     * floods this work targets. It is retained only for reproducing published
     * anugaSed results, as spec 3.4 recommends, and is not the default.
     *
     * S is the magnitude of the bed gradient, taken cell-locally from the
     * divergence theorem, grad z = (1/A) sum_e z_e n_e L_e, using the edge
     * values already to hand. No neighbour access, so it offloads. */
    anuga_int sediment_shear_closure;
    /* Bounds on S for [T-7]/[T-7e]. max_slope caps S (<= 0: no cap).
     * slope_frozen (depth-slope only) keeps the S computed from the bed at
     * setup, held in sediment_slope_work, instead of re-reading the evolving
     * bed each step: the closure is otherwise self-amplifying, erosion
     * roughening the bed it takes its stress from. */
    double sediment_max_slope;
    anuga_int sediment_slope_frozen;

    /* S_ms of [G-3]: an optional EXTERNAL source of suspended sediment --
     * hillslope yield, tributary load, rainfall washoff -- in m/s, i.e. the
     * same units as E and D. (ns, n), or NULL when unused.
     *
     * Applied AFTER the [L-1]/[L-2] limiters, deliberately. Those bound the
     * BED EXCHANGE by what the bed and the water column can supply; an
     * external supply is not bed exchange and must not be clipped by them. It
     * is also what makes a Method of Manufactured Solutions possible: the
     * manufactured source has to reach the equation unmodified. */
    double* tracer_external_source;

    /* Near-bed concentration ratio d*(Z), spec 4.3 / open item S1a.
     *   0 = constant, use sediment_d_star (the P14/P13 d* = 1 limiting case)
     *   1 = Rouse, evaluate the fitted form of [S-4] per cell
     * sediment_reference_height is 'a' in [S-4], the near-bed reference height
     * that aSM16 requires but never states. Exposed rather than hidden: d*
     * varies by up to 13x across plausible a/h at high Z. */
    anuga_int sediment_d_star_mode;
    double* sediment_reference_height;   /* (ncl) a [m] */
    /* van Rijn-style floor a >= sediment_a_h_floor * h. Standard practice and
     * on by default (0.01), but exposed: it is the largest single divergence
     * from anugaSed, which applies no floor. See PHYSICS_SPEC 12, D4b. */
    double sediment_a_h_floor;
    /* [S-2b] Rouse-number correction: mode 1 divides Z by van Rijn's
     * beta = 1 + 2 (w_s/u*)^2 (capped at 2), and the scale multiplies Z
     * afterwards (1 = none). Both flatten or steepen the fitted profile. */
    anuga_int sediment_rouse_beta_mode;
    double sediment_rouse_scale;
    /* [D-3] adaptation lag of the near-bed concentration (Galappatti &
     * Vreugdenhil 1985). 0 = none: the exchange E - D = d* v_s (c_eq - c)
     * responds at once to the local flow. 1 = Armanini & Di Silvio (1988)
     * closed form for alpha(w_s/u*, a/h); 2 = a constant alpha. The
     * exchange becomes alpha v_s (c_eq - c): both E and D are scaled by
     * alpha/d*, so every equilibrium is unchanged and only the transient
     * slows. alpha -> 1 in the well-mixed limit, alpha -> h/a when stratified. */
    anuga_int sediment_adaptation_mode;
    double sediment_adaptation_alpha;
    /* [D-4] carried near-bed concentration (adaptation mode 3): fraction s
     * keeps its near-bed concentration c_b as tracer (nearbed_base + s),
     * advected with the flow and relaxed toward d* c over a settling time
     * (z_c - a)/w_s when the profile is collapsing; instantaneous in the
     * other direction. -1 = not registered. */
    anuga_int sediment_nearbed_base;
    /* [D-5] two-layer tuning, both equilibrium-preserving: the near-bed
     * layer thickness as a fraction of depth (<= 0: the reference height
     * a/h, floor included) and a factor on the rate at which the partition
     * relaxes (1 = the two-box estimate). */
    double sediment_layer_fraction;
    double sediment_exchange_factor;
    /* [D-5v] 1: the two layers are advected at their log-law mean velocities
     * (through tracer_speed_factor), 0: both at the depth-averaged one. */
    anuga_int sediment_velocity_profile;
    /* [L-4] maximum packing fraction bounding the near-bed concentration
     * c_b = d* c. Same constant that bounds E* in [E-1]. */
    double sediment_c_pack;
    /* Bed porosity lambda in [G-4]: dz/dt = (D-E)/(1-lambda). The sediment
     * VOLUME leaving suspension is (1-lambda) dz, the rest being pore space. */
    double sediment_porosity;
    /* Morphological acceleration factor M (Delft3D's MORFAC): the bed change
     * of every step, from the suspended exchange [G-4] and from bedload
     * [G-5], is multiplied by M, and the erodible-base caps [L-5] are
     * scaled to match. The water column is not touched: the suspension
     * adapts in seconds, the bed in hours, so M steps of bed change per
     * hydrodynamic step is the standard way to reach a morphological time
     * at a hydrodynamic cost. 1 = off. Bed and water-column sediment
     * budgets then differ by exactly M, by construction. */
    double sediment_morphological_factor;
    /* Coupling stage of spec 2.4:
     *   0 = FIXED BED (Phase 3): exchange acts on m only, z never moves.
     *   1 = EVOLVING BED (Phase 4): [G-4] Exner update is applied.
     * Both are legitimate published configurations, not a debug switch. */
    anuga_int sediment_bed_evolution;

    /* Bedload, spec 6.  dz/dt = -(1/(1-lambda)) div q_b   [K-3]/[G-5]
     *
     *   0 = off
     *   1 = power law [K-1]/[K-2]:  q_b* = K tau_x^m,
     *       q_b = q_b* sqrt(R g) d^1.5,  tau_x = tau* - tau_c*
     *   2 = Engelund-Hansen [K-5]:  q_b* = 0.05 tau*^2.5 / f_c, NO threshold
     *
     * Mode 2 is a TOTAL LOAD relation with suspension already in it, so the
     * suspended source must not run alongside it (spec 6, "critical usage
     * rule"). The Python API refuses the combination rather than warning.
     *
     * K, m and tau_c* are exposed rather than hard-coded because which Wong &
     * Parker relation FG21 used is still open: Eq 24 gives 3.97/1.5/0.0495,
     * Eq 23 gives 4.93/1.60/0.0470. Resolving it should be a default change. */
    anuga_int sediment_bedload_mode;
    double sediment_bedload_K;
    double sediment_bedload_m;
    double sediment_bedload_tau_c_star;
    /* [K-7] depth below which bedload is ramped to zero: q_b is scaled by
     * clamp((h - h_min)/h_min, 0, 1), so it is off below h_min and full
     * above 2 h_min. 0 (default) = no ramp. Bedload relations assume a flow
     * many grains deep; at a wetting front over an erodible bed the film
     * cells carry the momentum of the front and a large nominal shear, and
     * without the ramp they dig steps that collapse the time step. */
    double sediment_bedload_h_min;
    /* Per-cell bedload transport vector, [K-4]. Scratch: filled and consumed
     * within one call, but device-resident so the divergence pass can read a
     * neighbour's value. */
    double* sediment_qbx;
    double* sediment_qby;
    /* Per-cell Rusanov coefficient for the bedload edge flux, gamma / h with
     * gamma = h |dq_b/dz| / q_b at fixed discharge and stage. Scratch. */
    double* sediment_qba;

    /* Friction closure for the sediment kernel, spec 3.3.
     *   0 = constant     n from friction_centroid_values (default)
     *   1 = larsen_lamb  n uniform, from [T-14]/[T-15]; sediment_manning_ll
     *   2 = wilson       f_c per cell from bed type, [T-8]..[T-10]
     * In ALL modes f_c still varies per cell per timestep, because [T-6]
     * depends on h. Spec 3.3 calls that the coupling most easily missed.
     *
     * NOTE ON THE FACTOR OF 8. W04 write their equations as (8/f_c)^1/2, but
     * THEIR f_c is the Darcy-Weisbach f (their Eq 4: U = [(8gRS)/f_c]^1/2),
     * whereas f_c in this struct is the spec's f_c = f/8. The two collide on
     * the same symbol. Working it through, f_c(ours) = f/8 = 1/X^2 where X is
     * the right-hand side of [T-8]..[T-10] -- the eights cancel. Using W04's
     * value directly would make tau_b 8x too large. */
    anuga_int sediment_friction_mode;
    double sediment_manning_ll;     /* [T-14] uniform n for larsen_lamb */
    anuga_int sediment_wilson_bed;  /* 0 sand [T-8], 1 gravel [T-9], 2 boulder [T-10] */
    double sediment_wilson_D;       /* D50 for sand, D84 for gravel/boulder [m] */

    /* ---- [L-5] non-erodible base -----------------------------------------
     *
     * sediment_z_base is a per-CENTROID bedrock elevation: erosion may lower
     * the bed to it and no further. It is a full (n) array, not a scalar,
     * because bedrock is a surface -- a scoured channel over an outcrop, a
     * lined culvert, a dam apron -- and a single number cannot describe one.
     *
     * The available erodible thickness in cell k is bed_centroid_values[k] -
     * sediment_z_base[k], which the limiter never lets go negative.
     *
     * sediment_has_z_base is the enable flag. With no base configured the
     * kernels take exactly the path they took before this existed, so the
     * default costs nothing and cannot change an existing answer.
     */
    double* sediment_z_base;           /* (n) bedrock centroid elevation [m] */
    /* [T-16] Per-centroid factor on the bed shear seen by the sediment
     * (f_c is multiplied by it in both the suspended source and the bedload
     * kernel), or NULL for 1 everywhere. The depth-averaged shear misses
     * the local amplification at obstacles (the horseshoe vortex at a pier,
     * two to four times the approach shear within a diameter of it); this
     * is where a user or a structure operator puts it. Input, set once. */
    double* sediment_shear_factor;
    anuga_int sediment_has_z_base;     /* 0 = unlimited depth (default) */


    /* Scratch, (ncl x n), tracer-major like the tracer arrays.
     *
     * The suspended source is computed for every class BEFORE any of it is
     * applied, because [L-5] scales the classes against a shared budget: the
     * erodible thickness is a property of the CELL, not of a class, so the
     * classes must be limited together. Applying class 0 in full and letting
     * class 1 take what is left would make the answer depend on registration
     * order, which is not physics. */
    double* sediment_source_limited;   /* (ncl*n) [G-3] bed exchange, per class */

    /* Bed (or free-surface) slope magnitude per cell for the depth-slope
     * closures [T-7]/[T-7e], filled by a pass of its own BEFORE the source
     * loop. The slope is a gradient over the cell's neighbours' centroid
     * values, and the source loop writes those same values (bed evolution),
     * so computing it inside that loop would read a neighbour a preceding
     * iteration or another thread had already moved: the answer then
     * depended on loop order, and mode 1 and mode 2 disagreed. */
    double* sediment_slope_work;       /* (n) S for [T-7]/[T-7e], scratch */

    /* Whether each cell had reached its base at the START of the bedload
     * step, snapshotted because bedload's second pass WRITES bed elevation
     * while reading its neighbours' state. Testing bed_cv against z_base
     * inside that loop would read a value a neighbouring thread may already
     * have updated, so the two cells sharing an edge could disagree about
     * whether it is open -- and a bedload edge flux the two sides disagree
     * about is no longer antisymmetric, which is precisely what makes the
     * scheme conservative. */
    anuga_int* sediment_bed_exhausted; /* (n) */
    /* Per boundary edge (boundary_length), indexed by -neighbour - 1: 1 where
     * bedload passes through the edge with the cell's own q_b.n (an inflow or
     * outflow declared open in set_bedload), 0 where it is closed. NULL means
     * every boundary edge is closed. */
    anuga_int* sediment_bedload_open;
    /* Per boundary edge, indexed like sediment_bedload_open: a PRESCRIBED
     * bedload inflow across that edge, m^2/s volumetric per unit width,
     * or < 0 where none is prescribed (the open edge then carries the
     * cell's own q_b . n, zero gradient). A prescribed supply does not
     * depend on the inflow cell's state, which is what removes the
     * feedback of the zero-gradient import: a cell that aggrades under a
     * fixed inflow stage sees its transport and hence its import rise
     * without bound. */
    double* sediment_bedload_supply;

    /* ---- spec 7, angle-of-repose relaxation ------------------------------
     *
     * Where the centroid-to-centroid bed slope exceeds a threshold angle, bed
     * material is moved downslope until it does not. FG21 are explicit that
     * this is a NUMERICAL HEURISTIC, not physics -- real slope failures are
     * advective -- and that it suppresses knickpoint retreat that may be real.
     * It exists to stop the rest of the model breaking on over-steep slopes.
     *
     * This is the ONLY non-cell-local sediment kernel: a cell's update depends
     * on its neighbours' elevation, so it is swept iteratively, Jacobi style,
     * with sediment_repose_dz holding the sweep's increments. Reading live
     * elevation instead would make the result depend on thread order.
     *
     * sediment_repose_tan = 0 disables it, which is the default.
     */
    double* sediment_repose_dz;            /* (n) per-sweep increment [m] */
    double sediment_repose_tan;            /* tan of the critical angle */
    double sediment_repose_relax;          /* under-relaxation, (0, 1] */
    anuga_int sediment_repose_max_sweeps;  /* hard cap; reported when hit */

    /* ------------------------------------------------------------------
     * Flux-path opt-ins for the fused GPU step (gpu/core_kernels.c).
     *
     * Appended at the END of the struct, like the tracer and sediment blocks
     * above, so no pre-existing field offset moves (see the note at the top
     * of the tracer block).  Both pyx extensions set all three explicitly:
     * the legacy Domain_C_struct comes from PyMem_Malloc, so an unset field
     * is garbage, not zero (claude/KNOWN_ISSUES.md).
     *
     * reconstruct_edge_bed: 0 unless the driver guarantees that fluxes
     * always follow an extrapolate (direct compute_fluxes callers such as
     * test_flux rely on an independently set bed_ev):
     *   1 -- cell-based kernel reconstructs edge bed values as stage - height
     *        (bit-identical to bed_ev after any extrapolate, one less gather)
     *   2 -- additionally selects the scatter flux kernel (single Riemann
     *        solve per edge, atomic accumulation; not valid with riverwalls,
     *        sloped Manning or tracers).  Scatter also requires owned_edges.
     * ------------------------------------------------------------------ */
    anuga_int reconstruct_edge_bed;
    /* Compacted list of the cell-edge slots the scatter flux kernel computes:
     * every boundary slot, plus the side of each interior edge whose
     * neighbour index is larger.  Built once by the driver from `neighbours`
     * (about 1.5 elements per cell + boundary).  Lets the kernel run one
     * thread per PHYSICAL edge instead of one per slot with half of them
     * exiting immediately.  NULL + 0 when scatter mode is unused. */
    anuga_int  num_owned_edges;
    anuga_int* owned_edges;

    /* [E-6] de Leeuw et al. (2020) entrainment, sediment_erosion_mode 3:
     *   E* = A X^beta / (1 + 3 A X^beta),  X = (u*_skin / v_s)^alpha Fr - threshold
     * with the skin-friction shear velocity from Manning-Strickler on a
     * roughness k_s (de Leeuw Eq 7). Appended at the end so no earlier member
     * moves (see the note on the tracer block). */
    double sediment_dl_A;
    double sediment_dl_alpha;
    double sediment_dl_beta;
    double sediment_dl_threshold;
    double sediment_dl_ks;          /* skin roughness height [m]; <= 0: ks_factor * diameter */
    double sediment_dl_ks_factor;   /* k_s = 3 D84 in de Leeuw et al. */

    /* Vegetation seen by the sediment shear, when vegetation_mode > 0:
     *   0 ignore (the cell's own friction closure),
     *   1 bed share of the Baptist resistance (default),
     *   2 the whole vegetated resistance.
     * See core_vegetation_sediment_fc. Appended at the end. */
    anuga_int sediment_vegetation_shear;

    /* ---- bed composition (active layer) -----------------------------------
     *
     * Per class and cell, the solid volume per bed area [m] held in the active
     * layer and in the substrate below it, tracer-major (ncl x n) like the
     * tracers. Entrainment of class s is multiplied by its share of the
     * active layer, erosion of s is capped at what the cell holds of s, and
     * after the exchange the active layer is restored to
     * sediment_active_layer (a bed thickness): overflow goes down at the
     * active composition, a deficit is refilled from the substrate at the
     * substrate composition. Amounts, not fractions, so each class is
     * conserved exactly. sediment_bed_composition == 0 (default) leaves the
     * kernels on the path they took before this existed. */
    anuga_int sediment_bed_composition;
    double sediment_active_layer;      /* active-layer bed thickness [m] */
    double* sediment_bed_active;       /* (ncl*n) */
    double* sediment_bed_substrate;    /* (ncl*n) */
};


struct edge {

    anuga_int cell_id;
    anuga_int edge_id;

    // mid point values
    double w;
    double h;
    double z;
    double uh;
    double vh;
    double u;
    double v;

    // vertex values
    double w1;
    double h1;
    double z1;
    double uh1;
    double vh1;
    double u1;
    double v1;

    double w2;
    double h2;
    double z2;
    double uh2;
    double vh2;
    double u2;
    double v2;

};


static inline void get_edge_data(struct edge *E, struct domain *D, anuga_int k, anuga_int i) {
    // fill edge data (conserved and bed) for ith edge of kth triangle

    anuga_int k3i, k3i1, k3i2;

    k3i = 3 * k + i;
    k3i1 = 3 * k + (i + 1) % 3;
    k3i2 = 3 * k + (i + 2) % 3;

    E->cell_id = k;
    E->edge_id = i;

    E->w = D->stage_edge_values[k3i];
    E->z = D->bed_edge_values[k3i];
    E->h = E->w - E->z;
    E->uh = D->xmom_edge_values[k3i];
    E->vh = D->ymom_edge_values[k3i];

    E->w1 = D->stage_vertex_values[k3i1];
    E->z1 = D->bed_vertex_values[k3i1];
    E->h1 = E->w1 - E->z1;
    E->uh1 = D->xmom_vertex_values[k3i1];
    E->vh1 = D->ymom_vertex_values[k3i1];


    E->w2 = D->stage_vertex_values[k3i2];
    E->z2 = D->bed_vertex_values[k3i2];
    E->h2 = E->w2 - E->z2;
    E->uh2 = D->xmom_vertex_values[k3i2];
    E->vh2 = D->ymom_vertex_values[k3i2];

}

static inline anuga_int print_domain_struct(struct domain *D) {


    printf("D->number_of_elements     %" PRId64 "  \n", D->number_of_elements);
    printf("D->boundary_length        %" PRId64 "  \n", D->boundary_length);
    printf("D->number_of_riverwall_edges %" PRId64 "  \n", D->number_of_riverwall_edges);
    printf("D->epsilon                %g \n", D->epsilon);
    printf("D->H0                     %g \n", D->H0);
    printf("D->g                      %g \n", D->g);
    printf("D->optimise_dry_cells     %" PRId64 " \n", D->optimise_dry_cells);
    printf("D->evolve_max_timestep    %g \n", D->evolve_max_timestep);
    printf("D->evolve_min_timestep    %g \n", D->evolve_min_timestep);
    printf("D->minimum_allowed_height %g \n", D->minimum_allowed_height);
    printf("D->maximum_allowed_speed  %g \n", D->maximum_allowed_speed);
    printf("D->low_froude             %" PRId64 " \n", D->low_froude);
    printf("D->extrapolate_velocity_second_order %" PRId64 " \n", D->extrapolate_velocity_second_order);
    printf("D->beta_w                 %g \n", D->beta_w);
    printf("D->beta_w_dry             %g \n", D->beta_w_dry);
    printf("D->beta_uh                %g \n", D->beta_uh);
    printf("D->beta_uh_dry            %g \n", D->beta_uh_dry);
    printf("D->beta_vh                %g \n", D->beta_vh);
    printf("D->beta_vh_dry            %g \n", D->beta_vh_dry);



    printf("D->neighbours             %p \n", (void *) D->neighbours);
    printf("D->surrogate_neighbours   %p \n", (void *) D->surrogate_neighbours);
    printf("D->neighbour_edges        %p \n", (void *) D->neighbour_edges);
    printf("D->normals                %p \n", (void *) D->normals);
    printf("D->edgelengths            %p \n", (void *) D->edgelengths);
    printf("D->radii                  %p \n", (void *) D->radii);
    printf("D->areas                  %p \n", (void *) D->areas);
    printf("D->tri_full_flag          %p \n", (void *) D->tri_full_flag);
    printf("D->already_computed_flux  %p \n", (void *) D->already_computed_flux);
    printf("D->vertex_coordinates     %p \n", (void *) D->vertex_coordinates);
    printf("D->edge_coordinates       %p \n", (void *) D->edge_coordinates);
    printf("D->centroid_coordinates   %p \n", (void *) D->centroid_coordinates);
    printf("D->max_speed              %p \n", (void *) D->max_speed);
    printf("D->number_of_boundaries   %p \n", (void *) D->number_of_boundaries);
    printf("D->stage_edge_values      %p \n", (void *) D->stage_edge_values);
    printf("D->xmom_edge_values       %p \n", (void *) D->xmom_edge_values);
    printf("D->ymom_edge_values       %p \n", (void *) D->ymom_edge_values);
    printf("D->bed_edge_values        %p \n", (void *) D->bed_edge_values);
    printf("D->stage_centroid_values  %p \n", (void *) D->stage_centroid_values);
    printf("D->xmom_centroid_values   %p \n", (void *) D->xmom_centroid_values);
    printf("D->ymom_centroid_values   %p \n", (void *) D->ymom_centroid_values);
    printf("D->bed_centroid_values    %p \n", (void *) D->bed_centroid_values);
    printf("D->stage_vertex_values    %p \n", (void *) D->stage_vertex_values);
    printf("D->xmom_vertex_values     %p \n", (void *) D->xmom_vertex_values);
    printf("D->ymom_vertex_values     %p \n", (void *) D->ymom_vertex_values);
    printf("D->bed_vertex_values      %p \n", (void *) D->bed_vertex_values);
    printf("D->height_vertex_values      %p \n", (void *) D->height_vertex_values);
    printf("D->stage_boundary_values  %p \n", (void *) D->stage_boundary_values);
    printf("D->xmom_boundary_values   %p \n", (void *) D->xmom_boundary_values);
    printf("D->ymom_boundary_values   %p \n", (void *) D->ymom_boundary_values);
    printf("D->bed_boundary_values    %p \n", (void *) D->bed_boundary_values);
    printf("D->stage_explicit_update  %p \n", (void *) D->stage_explicit_update);
    printf("D->xmom_explicit_update   %p \n", (void *) D->xmom_explicit_update);
    printf("D->ymom_explicit_update   %p \n", (void *) D->ymom_explicit_update);
    printf("D->edge_river_wall_counter   %p \n",   (void *) D->edge_river_wall_counter);
    printf("D->stage_semi_implicit_update  %p \n", (void *) D->stage_semi_implicit_update);
    printf("D->xmom_semi_implicit_update   %p \n", (void *) D->xmom_semi_implicit_update);
    printf("D->ymom_semi_implicit_update   %p \n", (void *) D->ymom_semi_implicit_update);
    printf("D->friction_centroid_values   %p \n", (void *) D->friction_centroid_values);


    return 0;
}


typedef struct {
    double ql[3], qr[3];
    double zl, zr;
    double hle, hre;
    double h_left, h_right;
    double hc, zc, hc_n, zc_n;
    double z_half;
    double normal_x, normal_y;
    double length;
    int n; // neighbour index
    int ki, ki2;
    bool is_boundary;
    bool is_riverwall;
    int riverwall_index;
} EdgeData;

// Extract edge-related data and organize it into EdgeData
static inline void get_edge_data_central_flux(const struct domain * __restrict D, const int k, const int i, EdgeData * __restrict E) {
    E->ki = 3 * k + i;
    E->ki2 = 2 * E->ki;

    E->ql[0] = D->stage_edge_values[E->ki];
    E->ql[1] = D->xmom_edge_values[E->ki];
    E->ql[2] = D->ymom_edge_values[E->ki];
    E->zl = D->bed_edge_values[E->ki];
    E->hle = D->height_edge_values[E->ki];
    E->length = D->edgelengths[E->ki];

    E->n = D->neighbours[E->ki];
    E->is_boundary = (E->n < 0);
    E->normal_x = D->normals[E->ki2];
    E->normal_y = D->normals[E->ki2 + 1];

    E->hc = D->height_centroid_values[k];
    E->zc = D->bed_centroid_values[k];
    E->hc_n=E->hc;
    E->zc_n=D->bed_centroid_values[k];

    if (E->is_boundary) {
        int m = -E->n - 1;
        E->qr[0] = D->stage_boundary_values[m];
        E->qr[1] = D->xmom_boundary_values[m];
        E->qr[2] = D->ymom_boundary_values[m];
        E->zr = E->zl;
        E->hre = fmax(E->qr[0] - E->zr, 0.0);
    } else {
        E->hc_n = D->height_centroid_values[E->n];
        E->zc_n = D->bed_centroid_values[E->n];
        int m = D->neighbour_edges[E->ki];
        int nm = E->n * 3 + m;
        E->qr[0] = D->stage_edge_values[nm];
        E->qr[1] = D->xmom_edge_values[nm];
        E->qr[2] = D->ymom_edge_values[nm];
        E->zr = D->bed_edge_values[nm];
        E->hre = D->height_edge_values[nm];
    }

    E->z_half = fmax(E->zl, E->zr);

    // Check for riverwall elevation override (skip entirely when no riverwalls)
    E->is_riverwall = (D->number_of_riverwall_edges > 0 &&
                       D->edge_flux_type != NULL &&
                       D->edge_flux_type[E->ki] == 1);
    if (E->is_riverwall) {
        E->riverwall_index = D->edge_river_wall_counter[E->ki] - 1;
        double zwall = D->riverwall_elevation[E->riverwall_index];
        E->z_half = fmax(zwall, E->z_half);
    }

    E->h_left = fmax(E->hle + E->zl - E->z_half, 0.0);
    E->h_right = fmax(E->hre + E->zr - E->z_half, 0.0);
}


#endif
