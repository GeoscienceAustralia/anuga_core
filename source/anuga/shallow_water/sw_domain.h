// C struct for domain and quantities
//
// Stephen Roberts 2012



// Shared code snippets
#include "util_ext.h"


// structures
struct domain {
    // Changing these don't change the data in python object
    long    number_of_elements;
    double  epsilon;
    double  H0;
    double  g;
    long    optimise_dry_cells;
    double  evolve_max_timestep;

    // The values in the python object will be changed
    long*   neighbours;
    long*   neighbour_edges;
    double* normals;
    double* edgelengths;
    double* radii;
    double* areas;

    long*   tri_full_flag;
    long*   already_computed_flux;
    double* max_speed;

    double* vertex_coordinates;

    double* stage_edge_values;
    double* xmom_edge_values;
    double* ymom_edge_values;
    double* bed_edge_values;

    double* stage_centroid_values;
    double* xmom_centroid_values;
    double* ymom_centroid_values;
    double* bed_centroid_values;

    double* stage_vertex_values;
    double* xmom_vertex_values;
    double* ymom_vertex_values;
    double* bed_vertex_values;


    double* stage_boundary_values;
    double* xmom_boundary_values;
    double* ymom_boundary_values;
    double* bed_boundary_values;

    double* stage_explicit_update;
    double* xmom_explicit_update;
    double* ymom_explicit_update;
};


struct edge {

    int cell_id;
    int edge_id;

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


void get_edge_data(struct edge *E, struct domain *D, int k, int i) {
    // fill edge data (conserved and bed) for ith edge of kth triangle

    int k3i, k3i1, k3i2;

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


struct domain* get_python_domain(struct domain *D, PyObject *domain) {
    PyArrayObject
            *neighbours,
            *neighbour_edges,
            *normals,
            *edgelengths,
            *radii,
            *areas,
            *tri_full_flag,
            *already_computed_flux,
            *vertex_coordinates,
            *max_speed;

    PyObject *quantities;

    D->number_of_elements = get_python_integer(domain, "number_of_elements");
    D->epsilon = get_python_double(domain, "epsilon");
    D->H0 = get_python_double(domain, "H0");
    D->g = get_python_double(domain, "g");
    D->optimise_dry_cells = get_python_integer(domain, "optimise_dry_cells");
    D->evolve_max_timestep = get_python_double(domain, "evolve_max_timestep");

    neighbours = get_consecutive_array(domain, "neighbours");
    D->neighbours = (long *) neighbours->data;

    neighbour_edges = get_consecutive_array(domain, "neighbour_edges");
    D->neighbour_edges = (long *) neighbour_edges->data;

    normals = get_consecutive_array(domain, "normals");
    D->normals = (double *) normals->data;

    edgelengths = get_consecutive_array(domain, "edgelengths");
    D->edgelengths = (double *) edgelengths->data;

    radii = get_consecutive_array(domain, "radii");
    D->radii = (double *) radii->data;

    areas = get_consecutive_array(domain, "areas");
    D->areas = (double *) areas->data;

    tri_full_flag = get_consecutive_array(domain, "tri_full_flag");
    D->tri_full_flag = (long *) tri_full_flag->data;

    already_computed_flux = get_consecutive_array(domain, "already_computed_flux");
    D->already_computed_flux = (long *) already_computed_flux->data;

    vertex_coordinates = get_consecutive_array(domain, "vertex_coordinates");
    D->vertex_coordinates = (double *) vertex_coordinates->data;
    
    max_speed = get_consecutive_array(domain, "max_speed");
    D->max_speed = (double *) max_speed->data;


    quantities = get_python_object(domain, "quantities");

    D->stage_edge_values     = get_python_array_data_from_dict(quantities, "stage",     "edge_values");
    D->xmom_edge_values      = get_python_array_data_from_dict(quantities, "xmomentum", "edge_values");
    D->ymom_edge_values      = get_python_array_data_from_dict(quantities, "ymomentum", "edge_values");
    D->bed_edge_values       = get_python_array_data_from_dict(quantities, "elevation", "edge_values");

    D->stage_centroid_values     = get_python_array_data_from_dict(quantities, "stage",     "centroid_values");
    D->xmom_centroid_values      = get_python_array_data_from_dict(quantities, "xmomentum", "centroid_values");
    D->ymom_centroid_values      = get_python_array_data_from_dict(quantities, "ymomentum", "centroid_values");
    D->bed_centroid_values       = get_python_array_data_from_dict(quantities, "elevation", "centroid_values");

    D->stage_vertex_values     = get_python_array_data_from_dict(quantities, "stage",     "vertex_values");
    D->xmom_vertex_values      = get_python_array_data_from_dict(quantities, "xmomentum", "vertex_values");
    D->ymom_vertex_values      = get_python_array_data_from_dict(quantities, "ymomentum", "vertex_values");
    D->bed_vertex_values       = get_python_array_data_from_dict(quantities, "elevation", "vertex_values");

    D->stage_boundary_values = get_python_array_data_from_dict(quantities, "stage",     "boundary_values");
    D->xmom_boundary_values  = get_python_array_data_from_dict(quantities, "xmomentum", "boundary_values");
    D->ymom_boundary_values  = get_python_array_data_from_dict(quantities, "ymomentum", "boundary_values");
    D->bed_boundary_values   = get_python_array_data_from_dict(quantities, "elevation", "boundary_values");

    D->stage_explicit_update = get_python_array_data_from_dict(quantities, "stage",     "explicit_update");
    D->xmom_explicit_update  = get_python_array_data_from_dict(quantities, "xmomentum", "explicit_update");
    D->ymom_explicit_update  = get_python_array_data_from_dict(quantities, "ymomentum", "explicit_update");


    return D;
}