#include "sqpipe.h"

/**********************************************************
Square pipe default quantity class
**********************************************************/

struct sqpipe_default_quantity {
  struct sqpipe_quantity super;
};

// Construction
void sqpipe_default_quantity_init (struct sqpipe_default_quantity *q);
struct sqpipe_default_quantity *sqpipe_default_quantity_new();

// Destruction

// Methods

// Implementation
void sqpipe_default_quantity_init (struct sqpipe_default_quantity *q) {
  struct sqpipe_quantity *p = (struct sqpipe_quantity *)q;

  sqpipe_quantity_init(p);  
}

struct sqpipe_default_quantity *sqpipe_default_quantity_new() {
  struct sqpipe_default_quantity *p = malloc(sizeof(struct sqpipe_default_quantity));

  sqpipe_default_quantity_init(p);

  return p;
}
/**********************************************************/


/**********************************************************
Square pipe default edge class
**********************************************************/
struct sqpipe_default_edge {
  struct sqpipe_edge super;
};

// Construction
void sqpipe_default_edge_init (struct sqpipe_default_edge *e);
struct sqpipe_default_edge *sqpipe_default_edge_new();

// Destruction

// Methods

// Implementation
void sqpipe_default_edge_init (struct sqpipe_default_edge *q) {
  sqpipe_edge_init((struct sqpipe_edge *)q);
}

struct sqpipe_default_edge *sqpipe_default_edge_new() {
  struct sqpipe_default_edge *e = (struct sqpipe_default_edge *) sqpipe_edge_new();

  sqpipe_default_edge_init(e);

  // Create inner and outer quantities
  ((struct edge *)e)->qin = (struct quantity *)sqpipe_default_quantity_new();
  ((struct edge *)e)->qout = (struct quantity *)sqpipe_default_quantity_new();
  
  return e;
}
/**********************************************************/


/**********************************************************
Square pipe default cell class
**********************************************************/
struct sqpipe_default_cell {
  struct sqpipe_cell super;
};

// Construction
void sqpipe_default_cell_init (struct sqpipe_default_cell *e);
struct sqpipe_default_cell *sqpipe_default_cell_new();

// Destruction

// Methods

// Implementation
void sqpipe_default_cell_init (struct sqpipe_default_cell *c) {
  struct sqpipe_cell *p = (struct sqpipe_cell *) c;

  sqpipe_cell_init(p);
}

struct sqpipe_default_cell *sqpipe_default_cell_new() {
  int i;
  struct sqpipe_default_cell *c = malloc(sizeof(struct sqpipe_default_cell));

  // Allocate memory for a length 2 array of pointers to edges
  ((struct cell *)c)->edges = (struct edge  **) malloc(2 * sizeof(struct sqpipe_edge *));

  sqpipe_default_cell_init(c);

  // Create edges
  for (i=0; i<((struct cell *)c)->num_edges; i++) {
    ((struct cell *)c)->edges[i] = (struct edge *) sqpipe_default_edge_new();
  }

  return c;
}
/**********************************************************/


/**********************************************************
Square pipe default quantities class
**********************************************************/
struct sqpipe_default_quantities {
  struct sqpipe_quantities super;
};

// Construction
void sqpipe_default_quantities_init (struct sqpipe_default_quantities *q);
struct sqpipe_default_quantities *sqpipe_default_quantities_new();

// Destruction

// Methods

// Implementation
void sqpipe_default_quantities_init (struct sqpipe_default_quantities *q) {
  struct sqpipe_quantities *p = (struct sqpipe_quantities *)q;

  sqpipe_quantities_init(p);
}

struct sqpipe_default_quantities *sqpipe_default_quantities_new() {
  struct sqpipe_default_quantities *p = (struct sqpipe_default_quantities *) sqpipe_quantities_new();

  sqpipe_default_quantities_init(p);

  return p;
}
/**********************************************************/


/**********************************************************
Square pipe default domain class
**********************************************************/

struct sqpipe_default_domain {
  struct sqpipe_domain super;
};

// Construction
void sqpipe_default_domain_init (struct sqpipe_default_domain *d);
struct sqpipe_default_domain *sqpipe_default_domain_new();

// Destruction

// Methods
struct edge* sqpipe_default_domain_new_edge(struct domain *d);
struct cell* sqpipe_default_domain_new_cell(struct domain *d);

// Implementation
struct edge* sqpipe_default_domain_new_edge(struct domain *d) {
  struct edge *e = (struct edge*) sqpipe_default_edge_new();

  return e;
}

struct cell* sqpipe_default_domain_new_cell(struct domain *d) {
  struct cell *c = (struct cell*) sqpipe_default_cell_new();
  return c;
}

void sqpipe_default_domain_init (struct sqpipe_default_domain *d) {
  struct sqpipe_domain *p = (struct sqpipe_domain *)d;

  sqpipe_domain_init(p);

  static struct domain_vtable vtable = {
    &sqpipe_default_domain_new_edge,
    &sqpipe_default_domain_new_cell,
    &domain_compute_fluxes_cells_generic,
    &sqpipe_domain_get_cell,
    &sqpipe_domain_destroy
  };
  p->super.vtable = &vtable;
}

struct sqpipe_default_domain *sqpipe_default_domain_new() {
  struct sqpipe_default_domain *d = (struct sqpipe_default_domain *) sqpipe_domain_new();

  sqpipe_default_domain_init(d);

  ((struct domain *) d)->vertex_values = (struct quantities *) sqpipe_default_quantities_new();
  ((struct domain *) d)->boundary_values = (struct quantities *) sqpipe_default_quantities_new();
  ((struct domain *) d)->explicit_update = (struct quantities *) sqpipe_default_quantities_new();

  return d;
}
/**********************************************************/
