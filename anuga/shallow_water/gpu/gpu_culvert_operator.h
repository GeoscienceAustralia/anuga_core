// GPU-accelerated culvert (Boyd box/pipe) operator
// Batches all culverts into a single gather -> CPU compute -> scatter cycle
//
// Public struct definitions (culvert_params, culvert_indices, etc.) live in
// gpu_domain.h. This header adds internal structs used only within the
// culvert implementation, plus pure-computation function prototypes.

#ifndef GPU_CULVERT_OPERATOR_H
#define GPU_CULVERT_OPERATOR_H

// Forward declare - full definition in gpu_domain.h
struct gpu_domain;
struct culvert_params;

// ============================================================================
// Internal structs (used by gpu_culvert_operator.c only)
// ============================================================================

// Gathered data for one inlet (filled by GPU gather)
struct inlet_data {
    // Enquiry point values
    double enquiry_stage;
    double enquiry_xmom;
    double enquiry_ymom;
    double enquiry_elevation;

    // Area-weighted averages over inlet region
    double avg_stage;
    double avg_depth;
    double avg_xmom;
    double avg_ymom;
    double total_area;
};

// Result from discharge calculation for one culvert
struct culvert_result {
    double Q;                    // Discharge [m^3/s]
    double barrel_velocity;      // Velocity through culvert [m/s]
    double outlet_culvert_depth; // Depth at outlet [m]
    double flow_area;            // Flow cross-section area [m^2]
    int inflow_idx;              // Which inlet (0 or 1) is the inflow
};

// Transfer specification for scattering results back to GPU
struct culvert_transfer {
    // Inflow region updates
    double new_inflow_depth;
    double new_inflow_xmom;
    double new_inflow_ymom;

    // Outflow region updates
    double new_outflow_depth;
    double new_outflow_xmom;
    double new_outflow_ymom;

    int inflow_idx;              // Which inlet (0 or 1) is inflow
};

// ============================================================================
// Pure computation functions (CPU-side, no GPU dependencies)
// ============================================================================

void weir_orifice_trapezoid_discharge(const struct culvert_params *p,
                                      double driving_energy,
                                      double delta_total_energy,
                                      double outlet_enquiry_depth,
                                      double *Q_out, double *barrel_velocity_out,
                                      double *outlet_culvert_depth_out,
                                      double *flow_area_out);

void boyd_box_discharge(const struct culvert_params *p,
                        double driving_energy,
                        double delta_total_energy,
                        double outlet_enquiry_depth,
                        double *Q, double *barrel_velocity,
                        double *outlet_culvert_depth, double *flow_area);

void boyd_pipe_discharge(const struct culvert_params *p,
                         double inflow_depth,
                         double driving_energy,
                         double delta_total_energy,
                         double outlet_enquiry_depth,
                         double *Q, double *barrel_velocity,
                         double *outlet_culvert_depth, double *flow_area);

void culvert_smooth_energy(double *smooth_delta_total_energy,
                           double delta_total_energy,
                           double timestep,
                           double smoothing_timescale,
                           double *ts_out);

void culvert_smooth_discharge(double smooth_delta_total_energy,
                              double *smooth_Q,
                              double Q_in,
                              double flow_area,
                              double ts,
                              double *Q_out, double *velocity_out);

// Single implementation of the per-culvert Boyd/weir update (discharge +
// semi-implicit water transfer), shared by the mode-2 batch and the mode-1
// host path so both compute modes agree bit-for-bit. Assumes a local/master,
// in-bounds culvert; the caller handles non-master parallel skips.
void culvert_compute_one(const struct inlet_data *data0,
                         const struct inlet_data *data1,
                         const struct culvert_params *p,
                         struct culvert_state *st,
                         double timestep,
                         struct culvert_result *r,
                         struct culvert_transfer *t);

// Flat host entry point (scalars in/out) for the mode-1 Python path.
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
        double *outlet_culvert_depth);

// Host inlet gather (mode-1), bit-identical to the mode-2 device gather.
void culvert_gather_inlet_host(int n, const int *indices, const double *areas,
                               const double *stage_c, const double *xmom_c,
                               const double *ymom_c, const double *bed_c,
                               double total_area,
                               double *avg_stage, double *avg_depth,
                               double *avg_xmom, double *avg_ymom);

#endif // GPU_CULVERT_OPERATOR_H
