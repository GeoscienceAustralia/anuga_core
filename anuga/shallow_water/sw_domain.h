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
    long    extrapolate_velocity_second_order;
    double  minimum_allowed_height;

    long timestep_fluxcalls;

    double beta_w;
    double beta_w_dry;
    double beta_uh;
    double beta_uh_dry;
    double beta_vh;
    double beta_vh_dry;

    long max_flux_update_frequency;
    long ncol_riverwall_hydraulic_properties;

    // Changing values in these arrays will change the values in the python object
    long*   neighbours;
    long*   neighbour_edges;
    long*   surrogate_neighbours;
    double* normals;
    double* edgelengths;
    double* radii;
    double* areas;

    long* edge_flux_type;

    long*   tri_full_flag;
    long*   already_computed_flux;
    double* max_speed;

    double* vertex_coordinates;
    double* edge_coordinates;
    double* centroid_coordinates;

    long*   number_of_boundaries;
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

    long* flux_update_frequency;    
    long* update_next_flux;
    long* update_extrapolation;
    double* edge_timestep;
    double* edge_flux_work;
    double* pressuregrad_work;
    double* x_centroid_work;
    double* y_centroid_work;
    double* boundary_flux_sum;

    long* allow_timestep_increase;

    double* riverwall_elevation;
    long* riverwall_rowIndex;
    double* riverwall_hydraulic_properties;
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
            *edge_flux_type,
            *tri_full_flag,
            *already_computed_flux,
            *vertex_coordinates,
            *edge_coordinates,
            *centroid_coordinates,
            *number_of_boundaries,
            *surrogate_neighbours,
            *max_speed,
            *flux_update_frequency,
            *update_next_flux,
            *update_extrapolation,
            *allow_timestep_increase,
            *edge_timestep,
            *edge_flux_work,
            *pressuregrad_work,
            *x_centroid_work,
            *y_centroid_work,
            *boundary_flux_sum,
            *riverwall_elevation,
            *riverwall_rowIndex,
            *riverwall_hydraulic_properties;

    PyObject *quantities;
    PyObject *riverwallData;

    D->number_of_elements   = get_python_integer(domain, "number_of_elements");
    D->epsilon              = get_python_double(domain, "epsilon");
    D->H0                   = get_python_double(domain, "H0");
    D->g                    = get_python_double(domain, "g");
    D->optimise_dry_cells   = get_python_integer(domain, "optimise_dry_cells");
    D->evolve_max_timestep  = get_python_double(domain, "evolve_max_timestep");
    D->minimum_allowed_height = get_python_double(domain, "minimum_allowed_height");
    D->timestep_fluxcalls = get_python_integer(domain,"timestep_fluxcalls");
    

    D->extrapolate_velocity_second_order  = get_python_integer(domain, "extrapolate_velocity_second_order");

    D->beta_w      = get_python_double(domain, "beta_w");;
    D->beta_w_dry  = get_python_double(domain, "beta_w_dry");
    D->beta_uh     = get_python_double(domain, "beta_uh");
    D->beta_uh_dry = get_python_double(domain, "beta_uh_dry");
    D->beta_vh     = get_python_double(domain, "beta_vh");
    D->beta_vh_dry = get_python_double(domain, "beta_vh_dry");

    D->max_flux_update_frequency = get_python_integer(domain,"max_flux_update_frequency");
    
    neighbours = get_consecutive_array(domain, "neighbours");
    D->neighbours = (long *) neighbours->data;

    surrogate_neighbours = get_consecutive_array(domain, "surrogate_neighbours");
    D->surrogate_neighbours = (long *) surrogate_neighbours->data;

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

    edge_flux_type = get_consecutive_array(domain, "edge_flux_type");
    D->edge_flux_type = (long *) edge_flux_type->data;


    tri_full_flag = get_consecutive_array(domain, "tri_full_flag");
    D->tri_full_flag = (long *) tri_full_flag->data;

    already_computed_flux = get_consecutive_array(domain, "already_computed_flux");
    D->already_computed_flux = (long *) already_computed_flux->data;

    vertex_coordinates = get_consecutive_array(domain, "vertex_coordinates");
    D->vertex_coordinates = (double *) vertex_coordinates->data;

    edge_coordinates = get_consecutive_array(domain, "edge_coordinates");
    D->edge_coordinates = (double *) edge_coordinates->data;


    centroid_coordinates = get_consecutive_array(domain, "centroid_coordinates");
    D->centroid_coordinates = (double *) centroid_coordinates->data;
    
    max_speed = get_consecutive_array(domain, "max_speed");
    D->max_speed = (double *) max_speed->data;

    number_of_boundaries = get_consecutive_array(domain, "number_of_boundaries");
    D->number_of_boundaries = (long *) number_of_boundaries->data;

    flux_update_frequency = get_consecutive_array(domain, "flux_update_frequency");
    D->flux_update_frequency = (long*) flux_update_frequency->data;
    
    update_next_flux = get_consecutive_array(domain, "update_next_flux");
    D->update_next_flux = (long*) update_next_flux->data;
    
    update_extrapolation = get_consecutive_array(domain, "update_extrapolation");
    D->update_extrapolation = (long*) update_extrapolation->data;
    
    allow_timestep_increase = get_consecutive_array(domain, "allow_timestep_increase");
    D->allow_timestep_increase = (long*) allow_timestep_increase->data;

    edge_timestep = get_consecutive_array(domain, "edge_timestep");
    D->edge_timestep = (double*) edge_timestep->data;
    
    edge_flux_work = get_consecutive_array(domain, "edge_flux_work");
    D->edge_flux_work = (double*) edge_flux_work->data;
    
    pressuregrad_work = get_consecutive_array(domain, "pressuregrad_work");
    D->pressuregrad_work = (double*) pressuregrad_work->data;
    
    x_centroid_work = get_consecutive_array(domain, "x_centroid_work");
    D->x_centroid_work = (double*) x_centroid_work->data;

    y_centroid_work = get_consecutive_array(domain, "y_centroid_work");
    D->y_centroid_work = (double*) y_centroid_work->data;
    
    boundary_flux_sum = get_consecutive_array(domain, "boundary_flux_sum");
    D->boundary_flux_sum = (double*) boundary_flux_sum->data;

    quantities = get_python_object(domain, "quantities");

    D->stage_edge_values     = get_python_array_data_from_dict(quantities, "stage",     "edge_values");
    D->xmom_edge_values      = get_python_array_data_from_dict(quantities, "xmomentum", "edge_values");
    D->ymom_edge_values      = get_python_array_data_from_dict(quantities, "ymomentum", "edge_values");
    D->bed_edge_values       = get_python_array_data_from_dict(quantities, "elevation", "edge_values");
    D->height_edge_values    = get_python_array_data_from_dict(quantities, "height", "edge_values");

    D->stage_centroid_values     = get_python_array_data_from_dict(quantities, "stage",     "centroid_values");
    D->xmom_centroid_values      = get_python_array_data_from_dict(quantities, "xmomentum", "centroid_values");
    D->ymom_centroid_values      = get_python_array_data_from_dict(quantities, "ymomentum", "centroid_values");
    D->bed_centroid_values       = get_python_array_data_from_dict(quantities, "elevation", "centroid_values");
    D->height_centroid_values    = get_python_array_data_from_dict(quantities, "height", "centroid_values");

    D->stage_vertex_values     = get_python_array_data_from_dict(quantities, "stage",     "vertex_values");
    D->xmom_vertex_values      = get_python_array_data_from_dict(quantities, "xmomentum", "vertex_values");
    D->ymom_vertex_values      = get_python_array_data_from_dict(quantities, "ymomentum", "vertex_values");
    D->bed_vertex_values       = get_python_array_data_from_dict(quantities, "elevation", "vertex_values");
    D->height_vertex_values       = get_python_array_data_from_dict(quantities, "height", "vertex_values");

    D->stage_boundary_values = get_python_array_data_from_dict(quantities, "stage",     "boundary_values");
    D->xmom_boundary_values  = get_python_array_data_from_dict(quantities, "xmomentum", "boundary_values");
    D->ymom_boundary_values  = get_python_array_data_from_dict(quantities, "ymomentum", "boundary_values");
    D->bed_boundary_values   = get_python_array_data_from_dict(quantities, "elevation", "boundary_values");

    D->stage_explicit_update = get_python_array_data_from_dict(quantities, "stage",     "explicit_update");
    D->xmom_explicit_update  = get_python_array_data_from_dict(quantities, "xmomentum", "explicit_update");
    D->ymom_explicit_update  = get_python_array_data_from_dict(quantities, "ymomentum", "explicit_update");


    riverwallData = get_python_object(domain,"riverwallData");

    riverwall_elevation = get_consecutive_array(riverwallData, "riverwall_elevation");
    D->riverwall_elevation = (double*) riverwall_elevation->data;
    
    riverwall_rowIndex = get_consecutive_array(riverwallData, "hydraulic_properties_rowIndex");
    D->riverwall_rowIndex = (long*) riverwall_rowIndex->data;

    D->ncol_riverwall_hydraulic_properties = get_python_integer(riverwallData, "ncol_hydraulic_properties");

    riverwall_hydraulic_properties = get_consecutive_array(riverwallData, "hydraulic_properties");
    D->riverwall_hydraulic_properties = (double*) riverwall_hydraulic_properties->data;

    Py_DECREF(quantities);
    Py_DECREF(riverwallData);

    Py_DECREF(neighbours);
    Py_DECREF(surrogate_neighbours);
    Py_DECREF(neighbour_edges);
    Py_DECREF(normals);
    Py_DECREF(edgelengths);
    Py_DECREF(radii);
    Py_DECREF(areas);
    Py_DECREF(edge_flux_type);
    Py_DECREF(tri_full_flag);
    Py_DECREF(already_computed_flux);
    Py_DECREF(vertex_coordinates);
    Py_DECREF(edge_coordinates);
    Py_DECREF(centroid_coordinates);
    Py_DECREF(max_speed);
    Py_DECREF(number_of_boundaries);
    Py_DECREF(flux_update_frequency);
    Py_DECREF(update_next_flux);
    Py_DECREF(update_extrapolation);
    Py_DECREF(edge_timestep);
    Py_DECREF(edge_flux_work);
    Py_DECREF(pressuregrad_work);
    Py_DECREF(x_centroid_work);
    Py_DECREF(y_centroid_work);
    Py_DECREF(boundary_flux_sum);
    Py_DECREF(allow_timestep_increase);

    return D;
}



int print_domain_struct(struct domain *D) {


    printf("D->number_of_elements     %ld  \n", D->number_of_elements);
    printf("D->epsilon                %g \n", D->epsilon);
    printf("D->H0                     %g \n", D->H0);
    printf("D->g                      %g \n", D->g);
    printf("D->optimise_dry_cells     %ld \n", D->optimise_dry_cells);
    printf("D->evolve_max_timestep    %g \n", D->evolve_max_timestep);
    printf("D->minimum_allowed_height %g \n", D->minimum_allowed_height);
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


    return 0;
}

