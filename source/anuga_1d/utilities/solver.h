/**********************************************************
Generic params class
**********************************************************/
struct params;

// vtable
struct params_vtable {  
  struct params* (*destroy) (struct params *q);
};

// structure
struct params {
  struct params_vtable *vtable;
  
  double epsilon;
  double cfl;
};

// Public API

// Destruction
struct params* params_destroy(struct params *q);

// method declarations

// method implementations
struct params* params_destroy(struct params *q) {
  return q->vtable->destroy(q);
}
/**********************************************************/


/**********************************************************
Generic quantity class
**********************************************************/
struct quantity;

// vtable
struct quantity_vtable {
  double* (*flux_formula) (struct quantity *q, double normal, struct params *p, double *quantityflux);
  double (*sound_speed) (struct quantity *q, struct params *p);
  double (*get_conserved) (struct quantity *q, int k, double normal);
  struct quantity* (*destroy) (struct quantity *q);
};

// structure
struct quantity {
  struct quantity_vtable *vtable;
};

// Public API
// Destruction
struct quantity* quantity_destroy(struct quantity *q);

// method declarations
double* quantity_flux_formula (struct quantity *q, double normal, struct params *p, double *quantityflux);
double quantity_sound_speed (struct quantity *q, struct params *p);
double quantity_get_conserved (struct quantity *q, int k, double normal);

// method implementations
double* quantity_flux_formula (struct quantity *q, double normal, struct params *p, double *quantityflux) {
  return q->vtable->flux_formula(q, normal, p, quantityflux);
}

double quantity_sound_speed (struct quantity *q, struct params *p) {
  return q->vtable->sound_speed(q, p);
}

double quantity_get_conserved (struct quantity *q, int k, double normal) {
  return q->vtable->get_conserved(q, k, normal);
}

struct quantity* quantity_destroy(struct quantity *q) {
  return q->vtable->destroy(q);
}
/**********************************************************/


/**********************************************************
Generic edge class
**********************************************************/
struct edge;

// vtable
struct edge_vtable{
  double* (*flux_function) (struct edge *e, struct params *p, double *edgeflux, int num_eqns);
  double (*extreme_sound_speeds) (struct edge *e, struct params *p);
  struct edge* (*destroy) (struct edge *e);
};

// structure
struct edge {
  struct edge_vtable *vtable;
  struct quantity *qin, *qout;
  double normal;
  double smin, smax;
};

// Public API
// Destruction 
struct edge* edge_destroy(struct edge *e);

// method declarations
double* edge_flux_function (struct edge *e, struct params *p, double *edgeflux, int num_eqns);
double edge_extreme_sound_speeds (struct edge *e, struct params *p);

// method implementations
double* edge_flux_function (struct edge *e, struct params *p, double *edgeflux, int num_eqns) {
  return e->vtable->flux_function(e, p, edgeflux, num_eqns);
}

double edge_extreme_sound_speeds (struct edge *e, struct params *p) {
  return e->vtable->extreme_sound_speeds(e, p);
}

struct edge* edge_destroy(struct edge *e) {
  return e->vtable->destroy(e);
}
/**********************************************************/


/**********************************************************
Generic cell class
**********************************************************/
struct cell;

// vtable
struct cell_vtable{
  double* (*flux_function) (struct cell *c, struct params *p, double *cellflux, int num_eqns);
  double (*extreme_sound_speeds) (struct cell *c, struct params *p);
  double* (*forcing_terms) (struct cell *c, struct params *p, double *cellforce);
  struct cell* (*destroy) (struct cell *c);
};

// structure
struct cell {
  struct cell_vtable *vtable;
  struct edge **edges;
  int num_edges;
};

// Public API

// Destruction
struct cell* cell_destroy(struct cell *c);

// method declarations
double* cell_flux_function (struct cell *c, struct params *p, double *cellflux, int num_eqns);
double cell_extreme_sound_speeds (struct cell *c, struct params *p);
double* cell_forcing_terms (struct cell *c, struct params *p, double *cellforce);

// method implementations
double* cell_flux_function (struct cell *c, struct params *p, double *cellflux, int num_eqns) {
  return c->vtable->flux_function(c, p, cellflux, num_eqns);
}

double cell_extreme_sound_speeds (struct cell *c, struct params *p) {
  return c->vtable->extreme_sound_speeds(c, p);
}

double* cell_forcing_terms (struct cell *c, struct params *p, double *cellforce) {
  return c->vtable->forcing_terms(c, p, cellforce);
}

struct cell* cell_destroy(struct cell *c) {
  return c->vtable->destroy(c);
}
/**********************************************************/


/**********************************************************
Generic quantities class
**********************************************************/
struct quantities;

// vtable
struct quantities_vtable{
  struct quantity* (*get_quantity) (struct quantities *qs, int i, struct quantity* q);
  void (*update) (struct quantities *qs, double *flux, int k);
  struct quantities* (*destroy) (struct quantities *qs);
};

// structure
struct quantities {
  struct quantities_vtable *vtable;
};

// Public API

// Destruction 
struct quantities* quantities_destroy(struct quantities *qs);

// method declarations
struct quantity* quantities_get_quantity(struct quantities *qs, int i, struct quantity *q); 
void quantities_update (struct quantities *qs, double *flux, int k);

// method implementations
struct quantity* quantities_get_quantity(struct quantities *qs, int i, struct quantity *q) {
  return qs->vtable->get_quantity(qs, i, q);
}

void quantities_update (struct quantities *qs, double *flux, int k) {
  return qs->vtable->update(qs, flux, k);
}

struct quantities* quantities_destroy(struct quantities *qs) {
  return qs->vtable->destroy(qs);
}
/**********************************************************/

/**********************************************************
Generic domain class
**********************************************************/
struct domain;

// vtable
struct domain_vtable {
  struct edge* (*new_edge) (struct domain *D);
  struct cell* (*new_cell) (struct domain *D);
  double (*compute_fluxes) (struct domain *D);
  struct cell* (*get_cell) (struct domain *d, int i, struct cell *c); 
  struct domain* (*destroy) (struct domain *D);
};

// structure
struct domain {
  struct domain_vtable *vtable;

  struct quantities *vertex_values;
  struct quantities *boundary_values;
  struct quantities *explicit_update;
  long* neighbours;
  long* neighbour_vertices;
  double* normals;
  double* areas;
  int number_of_elements;
  int number_of_equations;
  double *max_speed_array;
  double timestep;
  struct params *params;
};

// Public API
// Destruction 
struct domain* domain_destroy(struct domain* D);

// method declarations
struct edge* domain_new_edge(struct domain *qs);
struct cell* domain_new_cell(struct domain *qs);
double domain_compute_fluxes (struct domain *D);
struct cell* domain_get_cell(struct domain *d, int i, struct cell *c); 

// method implementations
struct edge* domain_new_edge(struct domain *qs) {
  return qs->vtable->new_edge(qs);
}

struct cell* domain_new_cell(struct domain *qs) {
  return qs->vtable->new_cell(qs);
}

double domain_compute_fluxes (struct domain *D) {
  return D->vtable->compute_fluxes(D);
}

struct cell* domain_get_cell(struct domain *D, int i, struct cell *c) {
  return D->vtable->get_cell(D, i, c);
}

struct domain* domain_destroy(struct domain* D) {
  return D->vtable->destroy(D);
}
/**********************************************************/


/**********************************************************
Generic methods
**********************************************************/
// Declarations
double domain_compute_fluxes_cells_generic(struct domain *D);
double domain_compute_fluxes_edges_generic(struct domain *D);
struct cell* sqpipe_domain_get_cell(struct domain *D, int i, struct cell *c);
double* edge_flux_function_godunov (struct edge *c, struct params *p, double *edgeflux, int num_eqns);
double* cell_flux_function_generic (struct cell *c, struct params *p, double *cellflux, int num_eqns);

// Implementations

// Computational function for flux computation
// This iterates over cells.
double domain_compute_fluxes_cells_generic (struct domain *D) {
  double flux[D->number_of_equations];
  double max_speed;
  int k, i;
  struct cell *c;
   
  // pre-create a cell object to use in the loop
  // so we need only allocate memory once
  c = domain_new_cell(D);

  for (k=0; k<D->number_of_elements; k++) {
    c = domain_get_cell(D, k, c);
   
    // Compute sound speeds (this updates c) and update timestep
    max_speed = cell_extreme_sound_speeds(c, D->params);
    if (max_speed > D->params->epsilon) {				    
      D->timestep = min(D->timestep, 0.5*D->params->cfl * D->areas[k]/max_speed); 
    }

    // Compute flux
    cell_flux_function(c, D->params, flux, D->number_of_equations);    
    for (i=0; i<D->number_of_equations; i++) {
      flux[i] /= D->areas[k];
    }
    quantities_update(D->explicit_update, flux, k);

    //Keep track of maximal speeds
    D->max_speed_array[k]=max_speed;
  }
  
  cell_destroy(c);

  return D->timestep;	
}

// Computational function for flux computation
// This iterates over edges.
double domain_compute_fluxes_edges_generic (struct domain *D) {
  double flux[D->number_of_equations];
  double max_speed;
  int k, i;
  // pre-create a cell object to use in the loop
  // so we need only allocate memory once
  struct cell *c = domain_new_cell(D);

  for (k=0; k<D->number_of_elements; k++) {
    c = domain_get_cell(D, k, c);
   
    // Compute sound speeds (this updates c) and update timestep
    max_speed = cell_extreme_sound_speeds(c, D->params);
    if (max_speed > D->params->epsilon) {				    
      D->timestep = min(D->timestep, 0.5*D->params->cfl * D->areas[k]/max_speed); 
    }

    // Compute flux and store in explicit_update
    cell_flux_function(c, D->params, flux, D->number_of_equations);    
    for (i=0; i<D->number_of_equations; i++) {
      flux[i] /= D->areas[k];
    }
    quantities_update(D->explicit_update, flux, k);

    //Keep track of maximal speeds
    D->max_speed_array[k]=max_speed;
  }
  
  cell_destroy(c);

  return D->timestep;	
}


struct cell* domain_get_cell_generic(struct domain *D, int k, struct cell *c) {
  int i, ki, n, m, nm;

  for (i=0; i<c->num_edges; i++) {
    ki = k*(c->num_edges) + i;
    c->edges[i]->normal = D->normals[ki];      
    n = D->neighbours[ki];

    c->edges[i]->qin = quantities_get_quantity(D->vertex_values, ki, c->edges[i]->qin);

    if (n<0) {
      m = -n-1;
      c->edges[i]->qout = quantities_get_quantity(D->boundary_values, m, c->edges[i]->qout);
    } else {
      m = D->neighbour_vertices[ki];
      nm = n*2+m;
      c->edges[i]->qout = quantities_get_quantity(D->vertex_values, nm, c->edges[i]->qout);
    }
  } 

  return c;
}


// Generic Godunov flux computation method
double* edge_flux_function_godunov (struct edge *e, struct params *p, double *edgeflux, int num_eqns) {
  int i;
  double flux_in[num_eqns], flux_out[num_eqns];
  double denom;

  // Flux contribution from inside of  edge
  quantity_flux_formula(e->qin, e->normal, p, flux_in);

  // Flux contribution from outside of edge
  quantity_flux_formula(e->qout, e->normal, p, flux_out);

  // Flux computation
  denom = e->smax - e->smin;
  if (denom < p->epsilon) {
    for (i=0; i<num_eqns; i++) edgeflux[i] = 0.0;
  } else {
    for (i=0; i<num_eqns; i++) {
      edgeflux[i] = e->smax*flux_in[i] - e->smin*flux_out[i];
      edgeflux[i] += e->smax*e->smin * (quantity_get_conserved(e->qout, i, e->normal) - quantity_get_conserved(e->qin, i, e->normal));
      edgeflux[i] /= denom;
    }
  }

  return edgeflux;
}

double* cell_flux_function_generic (struct cell *c, struct params *p, double *cellflux, int num_eqns) {
  int i, j;
  double cellforce[num_eqns];
  double flux[num_eqns];

  // Compute forcing terms first
  cell_forcing_terms(c, p, cellforce);
  for (i=0; i<num_eqns; i++) {
    cellflux[i] = -cellforce[i];
  }

  // Compute flux through each edge and append to flux
  for (j=0; j<c->num_edges; j++) {
    edge_flux_function(c->edges[j], p, flux, num_eqns);
    for (i=0; i<num_eqns; i++) {
      cellflux[i] -= flux[i];
    }
  }

  return cellflux;
}
/**********************************************************/
