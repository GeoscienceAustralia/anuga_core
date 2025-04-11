// C struct for domain and quantities
//
// Stephen Roberts 2012



#ifndef SW_DOMAIN_H
#define SW_DOMAIN_H

#include <stdint.h>

// structures
struct domain {
    // Changing these don't change the data in python object
    int64_t    number_of_elements;
    int64_t    boundary_length;
    int64_t    number_of_riverwall_edges;
    double  epsilon;
    double  H0;
    double  g;
    int64_t    optimise_dry_cells;
    double  evolve_max_timestep;
    int64_t    extrapolate_velocity_second_order;
    double  minimum_allowed_height;
    double  maximum_allowed_speed;
    int64_t    low_froude;


    int64_t timestep_fluxcalls;

    double beta_w;
    double beta_w_dry;
    double beta_uh;
    double beta_uh_dry;
    double beta_vh;
    double beta_vh_dry;

    int64_t max_flux_update_frequency;
    int64_t ncol_riverwall_hydraulic_properties;

    // Changing values in these arrays will change the values in the python object
    int64_t*   neighbours;
    int64_t*   neighbour_edges;
    int64_t*   surrogate_neighbours;
    double* normals;
    double* edgelengths;
    double* radii;
    double* areas;

    int64_t* edge_flux_type;

    int64_t*   tri_full_flag;
    int64_t*   already_computed_flux;
    double* max_speed;

    double* vertex_coordinates;
    double* edge_coordinates;
    double* centroid_coordinates;

    int64_t*   number_of_boundaries;
    double* stage_edge_values;
    double* xmom_edge_values;
    double* ymom_edge_values;
    double* bed_edge_values;
    double* height_edge_values;

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

    double* stage_explicit_update;
    double* xmom_explicit_update;
    double* ymom_explicit_update;

    int64_t* flux_update_frequency;
    int64_t* update_next_flux;
    int64_t* update_extrapolation;
    double* edge_timestep;
    double* edge_flux_work;
    double* neigh_work;
    double* pressuregrad_work;
    double* x_centroid_work;
    double* y_centroid_work;
    double* boundary_flux_sum;

    int64_t* allow_timestep_increase;

    long* edge_river_wall_counter;
    double* riverwall_elevation;
    int64_t* riverwall_rowIndex;
    double* riverwall_hydraulic_properties;

    
};


struct edge {

    int64_t cell_id;
    int64_t edge_id;

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


void get_edge_data(struct edge *E, struct domain *D, int64_t k, int64_t i) {
    // fill edge data (conserved and bed) for ith edge of kth triangle

    int64_t k3i, k3i1, k3i2;

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

int64_t print_domain_struct(struct domain *D) {


    printf("D->number_of_elements     %ld  \n", D->number_of_elements);
    printf("D->boundary_length        %ld  \n", D->boundary_length);
    printf("D->number_of_riverwall_edges %ld  \n", D->number_of_riverwall_edges);
    printf("D->epsilon                %g \n", D->epsilon);
    printf("D->H0                     %g \n", D->H0);
    printf("D->g                      %g \n", D->g);
    printf("D->optimise_dry_cells     %ld \n", D->optimise_dry_cells);
    printf("D->evolve_max_timestep    %g \n", D->evolve_max_timestep);
    printf("D->minimum_allowed_height %g \n", D->minimum_allowed_height);
    printf("D->maximum_allowed_speed  %g \n", D->maximum_allowed_speed);
    printf("D->low_froude             %ld \n", D->low_froude);
    printf("D->extrapolate_velocity_second_order %ld \n", D->extrapolate_velocity_second_order);
    printf("D->beta_w                 %g \n", D->beta_w);
    printf("D->beta_w_dry             %g \n", D->beta_w_dry);
    printf("D->beta_uh                %g \n", D->beta_uh);
    printf("D->beta_uh_dry            %g \n", D->beta_uh_dry);
    printf("D->beta_vh                %g \n", D->beta_vh);
    printf("D->beta_vh_dry            %g \n", D->beta_vh_dry);



    printf("D->neighbours             %p \n", D->neighbours);
    printf("D->surrogate_neighbours   %p \n", D->surrogate_neighbours);
    printf("D->neighbour_edges        %p \n", D->neighbour_edges);
    printf("D->normals                %p \n", D->normals);
    printf("D->edgelengths            %p \n", D->edgelengths);
    printf("D->radii                  %p \n", D->radii);
    printf("D->areas                  %p \n", D->areas);
    printf("D->tri_full_flag          %p \n", D->tri_full_flag);
    printf("D->already_computed_flux  %p \n", D->already_computed_flux);
    printf("D->vertex_coordinates     %p \n", D->vertex_coordinates);
    printf("D->edge_coordinates       %p \n", D->edge_coordinates);
    printf("D->centroid_coordinates   %p \n", D->centroid_coordinates);
    printf("D->max_speed              %p \n", D->max_speed);
    printf("D->number_of_boundaries   %p \n", D->number_of_boundaries);
    printf("D->stage_edge_values      %p \n", D->stage_edge_values);
    printf("D->xmom_edge_values       %p \n", D->xmom_edge_values);
    printf("D->ymom_edge_values       %p \n", D->ymom_edge_values);
    printf("D->bed_edge_values        %p \n", D->bed_edge_values);
    printf("D->stage_centroid_values  %p \n", D->stage_centroid_values);
    printf("D->xmom_centroid_values   %p \n", D->xmom_centroid_values);
    printf("D->ymom_centroid_values   %p \n", D->ymom_centroid_values);
    printf("D->bed_centroid_values    %p \n", D->bed_centroid_values);
    printf("D->stage_vertex_values    %p \n", D->stage_vertex_values);
    printf("D->xmom_vertex_values     %p \n", D->xmom_vertex_values);
    printf("D->ymom_vertex_values     %p \n", D->ymom_vertex_values);
    printf("D->bed_vertex_values      %p \n", D->bed_vertex_values);
    printf("D->height_vertex_values      %p \n", D->height_vertex_values);
    printf("D->stage_boundary_values  %p \n", D->stage_boundary_values);
    printf("D->xmom_boundary_values   %p \n", D->xmom_boundary_values);
    printf("D->ymom_boundary_values   %p \n", D->ymom_boundary_values);
    printf("D->bed_boundary_values    %p \n", D->bed_boundary_values);
    printf("D->stage_explicit_update  %p \n", D->stage_explicit_update);
    printf("D->xmom_explicit_update   %p \n", D->xmom_explicit_update);
    printf("D->ymom_explicit_update   %p \n", D->ymom_explicit_update);
    printf("D->edge_river_wall_counter   %p \n", D->edge_river_wall_counter);


    return 0;
}

#endif
