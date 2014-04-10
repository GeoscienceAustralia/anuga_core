/**********************************************************
Pipe params class
**********************************************************/

struct sqpipe_params {
  struct params super;

  double g;
  double h0;
  double bulk_modulus;
};

// Construction
void sqpipe_params_init(struct sqpipe_params *q);
struct sqpipe_params* sqpipe_params_new();

// Destruction
struct params* sqpipe_params_destroy(struct params *q);

// Methods

// Implementation
void sqpipe_params_init(struct sqpipe_params *q) {
  static struct params_vtable vtable = {
    &sqpipe_params_destroy
  };
  q->super.vtable = &vtable;
  
  ((struct params *)q)->cfl = 0.0;  
  ((struct params *)q)->epsilon = 0.0;
  q->g = 0.0;
  q->h0 = 0.0;
}

struct sqpipe_params* sqpipe_params_new() {
  struct sqpipe_params *p = malloc(sizeof(struct sqpipe_params));

  sqpipe_params_init(p);

  return p;  
}

struct params* sqpipe_params_destroy(struct params *q) {
  free(q);
  q = NULL;
  
  return q;
}
/**********************************************************/


/**********************************************************
Pipe quantity class
This is a virtual class and should never be created.
The inherting class needs to implement
    init, new, flux_formula0, flux_formula1, sound_speed, 
    get_conserved
**********************************************************/

struct sqpipe_quantity {
  struct quantity super;

  double a;
  double d;
  double w;
  double h;
  double u;
  double z;  
  double b;
  double t;
  long state;
};

// Construction
void sqpipe_quantity_init(struct sqpipe_quantity *q);
struct sqpipe_quantity* sqpipe_quantity_new();

// Destruction
struct quantity* sqpipe_quantity_destroy(struct quantity *q);

// Methods
double* sqpipe_quantity_flux_formula (struct quantity *q, double normal, struct params *p, double *quantityflux);
double sqpipe_quantity_sound_speed (struct quantity *q, struct params *p);
double sqpipe_quantity_get_conserved (struct quantity *q, int k, double normal);

// Helper methods
double* sqpipe_quantity_free_surface_flux_formula (struct quantity *q, double normal, struct params *p, double *quantityflux);
double* sqpipe_quantity_pressurised_flux_formula (struct quantity *q, double normal, struct params *p, double *quantityflux);
double sqpipe_quantity_free_surface_sound_speed (struct quantity *q, struct params *p);
double sqpipe_quantity_pressurised_sound_speed (struct quantity *q, struct params *p);

// Implementation
double* sqpipe_quantity_flux_formula (struct quantity *q, double normal, struct params *p, double *quantityflux) {
  if (((struct sqpipe_quantity *)q)->state == 0) {
    return sqpipe_quantity_free_surface_flux_formula(q, normal, p, quantityflux);
  } else {
    return sqpipe_quantity_pressurised_flux_formula(q, normal, p, quantityflux);
  }
}

double sqpipe_quantity_sound_speed (struct quantity *q, struct params *p) {
  if (((struct sqpipe_quantity *)q)->state == 0) {
    return sqpipe_quantity_free_surface_sound_speed(q, p);
  } else {
    return sqpipe_quantity_pressurised_sound_speed(q, p);
  }
}

double sqpipe_quantity_get_conserved (struct quantity *q, int k, double normal) {
  struct sqpipe_quantity *p = (struct sqpipe_quantity*) q;
  double c;
  
  switch (k) {
  case 0:
    c = p->a;
    break;
  case 1:
    // This should be normal^2 p->d and normal is +/- 1
    c = p->d;
    break;
  default:
    c = 0;
  }
  return c;
}


// Flux formula taking into account the normal
// Note u, d should be normal*u, normal*d but since normal is +/- 1
// (normal*u) * (normal*d) = u*d
double* sqpipe_quantity_free_surface_flux_formula (struct quantity *q, double normal, struct params *p, double *quantityflux) {
  struct sqpipe_quantity *pq = (struct sqpipe_quantity*) q;
  struct sqpipe_params *sp = (struct sqpipe_params *) p;

  //quantityflux[0] = (normal * pq->u) * quantity_get_conserved(q, 0, normal);
  //quantityflux[1] = normal * (pq->u * quantity_get_conserved(q, 1, normal) + 0.5 * sp->g * pq->h *pq->h);
  quantityflux[0] = (normal * pq->u) * pq->a;
  quantityflux[1] = normal * (pq->u * pq->d + 0.5 * sp->g * pq->h *pq->h);

  return quantityflux;
}

double* sqpipe_quantity_pressurised_flux_formula (struct quantity *q, double normal, struct params *p, double *quantityflux) {
  struct sqpipe_quantity *pq = (struct sqpipe_quantity*) q;
  struct sqpipe_params *sp = (struct sqpipe_params *) p;
  
  quantityflux[0] = (normal * pq->u) * pq->a; 
  quantityflux[1] = normal * (pq->u * pq->d + sp->bulk_modulus*(pq->h - pq->t) + 0.5 * sp->g * pq->t * pq->t);

  return quantityflux;
}

double sqpipe_quantity_free_surface_sound_speed (struct quantity *q, struct params *p) {
  struct sqpipe_quantity *sq = (struct sqpipe_quantity*) q;
  struct sqpipe_params *sp = (struct sqpipe_params*) p;  

  return sqrt(sp->g * sq->h);
}

double sqpipe_quantity_pressurised_sound_speed (struct quantity *q, struct params *p) {
  return ((struct sqpipe_params *)p)->bulk_modulus;
}

void sqpipe_quantity_init(struct sqpipe_quantity *q) {
  static struct quantity_vtable vtable = {
    &sqpipe_quantity_flux_formula,
    &sqpipe_quantity_sound_speed,
    &sqpipe_quantity_get_conserved,
    &sqpipe_quantity_destroy
  };
  q->super.vtable = &vtable;
  
  q->a = 0.0;
  q->d = 0.0;
  q->w = 0.0;
  q->h = 0.0;
  q->u = 0.0;
  q->z = 0.0;  
  q->b = 0.0;
  q->t = 0.0;
  q->state = 0;
}

struct sqpipe_quantity* sqpipe_quantity_new() {
  struct sqpipe_quantity *p = malloc(sizeof(struct sqpipe_quantity));

  return p;  
}

struct quantity* sqpipe_quantity_destroy(struct quantity *q) {
  free(q);
  q = NULL;
  
  return q;
}
/**********************************************************/


/**********************************************************
Pipe edge class
This is a virtual class and should never be created.
The inherting class needs to implement
    init, new, flux_function
**********************************************************/
struct sqpipe_edge {
  struct edge super;
};

// Construction
void sqpipe_edge_init(struct sqpipe_edge *e);
struct sqpipe_edge* sqpipe_edge_new();

// Destruction
struct edge* sqpipe_edge_destroy(struct edge *e);

// Methods
double sqpipe_edge_extreme_sound_speeds (struct edge *e, struct params *p);

// Implementation
double sqpipe_edge_extreme_sound_speeds (struct edge *e, struct params *p) {
  double soundspeed_left, soundspeed_right, max_speed;
  struct sqpipe_quantity *qin = (struct sqpipe_quantity*) e->qin;
  struct sqpipe_quantity *qout = (struct sqpipe_quantity*) e->qout;

  soundspeed_left = quantity_sound_speed(e->qin, p);
  soundspeed_right = quantity_sound_speed(e->qout, p);

  e->smax = max(e->normal * qin->u + soundspeed_left, e->normal * qout->u + soundspeed_right);
  if (e->smax < 0.0) e->smax = 0.0;
	
  e->smin = min(e->normal * qin->u - soundspeed_left, e->normal * qout->u - soundspeed_right);
  if (e->smin > 0.0) e->smin = 0.0;

  if (e->smax - e->smin < p->epsilon) {
    max_speed = 0.0;
  } else {
    max_speed = max(fabs(e->smin), fabs(e->smax));
  }

  return max_speed;
}

void sqpipe_edge_init(struct sqpipe_edge *e) {
  static struct edge_vtable vtable = {
    &edge_flux_function_godunov,
    &sqpipe_edge_extreme_sound_speeds,
    &sqpipe_edge_destroy
  };
  e->super.vtable = &vtable;

  ((struct edge *)e)->qin = NULL;
  ((struct edge *)e)->qout = NULL;
  ((struct edge *)e)->normal = 0.0;
  ((struct edge *)e)->smax = 0.0;
  ((struct edge *)e)->smin = 0.0;
}

struct sqpipe_edge* sqpipe_edge_new() {
  struct sqpipe_edge *p = malloc(sizeof(struct sqpipe_edge));

  return p;  
}

struct edge* sqpipe_edge_destroy(struct edge *e) {  
  e->qin = quantity_destroy(e->qin);
  e->qout = quantity_destroy(e->qout);
  
  free(e);
  e = NULL;
  
  return e;
}
/**********************************************************/

/**********************************************************
Pipe cell class
This is a virtual class and should never be created.
The inherting class needs to implement
    init, new, flux_function, forcing_terms
**********************************************************/
struct sqpipe_cell {
  struct cell super;
  
  long state;
};

// Construction
void sqpipe_cell_init(struct sqpipe_cell *e);
struct sqpipe_cell* sqpipe_cell_new();

// Destruction
struct cell* sqpipe_cell_destroy(struct cell *c);

// Methods
double sqpipe_cell_extreme_sound_speeds (struct cell *c, struct params *p);
double* sqpipe_cell_forcing_terms (struct cell *c, struct params *p, double *cellforce);

// Helper methods
double* sqpipe_cell_free_surface_forcing_terms (struct cell *c, struct params *p, double *cellforce);
double* sqpipe_cell_pressurised_forcing_terms (struct cell *c, struct params *p, double *cellforce);

// Implementation
double sqpipe_cell_extreme_sound_speeds (struct cell *c, struct params *p) {
  int i;
  // The compiler doesn't know c->num_edges is positive hence whether the
  // loop will run. This throws the warning:
  // ‘max_speed’ may be used uninitialized in this function
  // and the same for tmp_max_speed
  // We initiallise them to avoid the warning since we know nothing could
  // possibly go wrong...
  double max_speed = 0.0, tmp_max_speed = 0.0;

  for (i=0; i<c->num_edges; i++) {
    tmp_max_speed = edge_extreme_sound_speeds(c->edges[i], p);
    max_speed = max(max_speed, tmp_max_speed);
  }

  return max_speed;
}

double* sqpipe_cell_forcing_terms (struct cell *c, struct params *p, double *cellforce) {
  if (((struct sqpipe_cell *)c)->state == 0) {
    return sqpipe_cell_free_surface_forcing_terms(c, p, cellforce);    
  } else {
    return sqpipe_cell_pressurised_forcing_terms(c, p, cellforce);
  }
}

double* sqpipe_cell_free_surface_forcing_terms (struct cell *c, struct params *p, double *cellforce) {
  struct sqpipe_quantity *qin = (struct sqpipe_quantity *) c->edges[0]->qout;
  struct sqpipe_quantity *qout = (struct sqpipe_quantity *) c->edges[1]->qin;
  struct sqpipe_params *sp = (struct sqpipe_params *) p;

  cellforce[0] = 0.0;  
  cellforce[1] = 0.5 * sp->g * (qout->h + qin->h) * (qout->z - qin->z);
  
  return cellforce;
}

double* sqpipe_cell_pressurised_forcing_terms (struct cell *c, struct params *p, double *cellforce) {
  struct sqpipe_quantity *qin = (struct sqpipe_quantity *) c->edges[0]->qout;
  struct sqpipe_quantity *qout = (struct sqpipe_quantity *) c->edges[1]->qin;
  struct sqpipe_params *sp = (struct sqpipe_params *) p;

  cellforce[0] = 0.0;
  cellforce[1] = sp->bulk_modulus * (sp->g * (0.5*(qout->h + qin->h) - 0.5*(qout->t + qin->t)) * (qout->b - qin->b) + sp->g * (0.5*(qout->h + qin->h) - 0.5*(qout->t + qin->t)) * 0.5*(qout->b + qin->b) * (qout->t - qin->t) * 1/(0.5*(qout->t + qin->t)));

  return cellforce;
}


void sqpipe_cell_init(struct sqpipe_cell *c) {
  static struct cell_vtable vtable = {
    &cell_flux_function_generic, 
    &sqpipe_cell_extreme_sound_speeds,
	&sqpipe_cell_forcing_terms,
    &sqpipe_cell_destroy
  };
  c->super.vtable = &vtable;

  // Can't set edges to NULL here since sqpipe_cell_new() has allocated memory
  // If a cell is created without using cell_new(), then edges will be an
  // unitialised pointer.
  // Can't allocate memory here since init should not allocate memory
  // The moral: Use cell_new() to create cells!

  c->super.num_edges = 2;
  ((struct sqpipe_cell *)c)->state = 0;
}

struct sqpipe_cell* sqpipe_cell_new() {
  struct sqpipe_cell *c = malloc(sizeof(struct sqpipe_cell));

  // Allocate memory for a length 2 array of pointers to edges
  ((struct cell *)c)->edges = (struct edge  **) malloc(2 * sizeof(struct sqpipe_edge *));
  
  return c;
}

struct cell* sqpipe_cell_destroy(struct cell *c) {  
  int i;

  for (i=0; i<c->num_edges; i++) {
    c->edges[i] = edge_destroy(c->edges[i]);
  }  

  free(c->edges);
  c->edges = NULL;
  free(c);
  c = NULL;
  
  return c;
}
/**********************************************************/


/**********************************************************
Pipe quantities class
This is a virtual class and should never be created.
The inherting class needs to implement:
    init, new, get_quantities
**********************************************************/
struct sqpipe_quantities {
  struct quantities super;

  double *a;
  double *d;
  double *w;
  double *h;
  double *u;
  double *z;  
  double *b;
  double *t;
};

// Construction
void sqpipe_quantities_init(struct sqpipe_quantities *q);
struct sqpipe_quantities* sqpipe_quantities_new();

// Destruction
struct quantities* sqpipe_quantities_destroy(struct quantities *qs);

// Methods
struct quantity* sqpipe_quantities_get_quantity(struct quantities *qs, int i, struct quantity *q);
void sqpipe_quantities_update (struct quantities *qs, double *flux, int k);

// Implementation
struct quantity* sqpipe_quantities_get_quantity(struct quantities *qs, int i, struct quantity *q) {  
  struct sqpipe_quantities *ps = (struct sqpipe_quantities *) qs;
  struct sqpipe_quantity *p = (struct sqpipe_quantity *) q;

  p->a = ps->a[i];
  p->d = ps->d[i];
  p->w = ps->w[i];
  p->h = ps->h[i];
  p->u = ps->u[i];  
  p->z = ps->z[i];
  p->b = ps->b[i];
  p->t = ps->t[i];

  return q;
}

void sqpipe_quantities_update (struct quantities *qs, double *flux, int k) {
  struct sqpipe_quantities *ps = (struct sqpipe_quantities *) qs;
  ps->a[k] = flux[0];
  ps->d[k] = flux[1];
}


void sqpipe_quantities_init(struct sqpipe_quantities *q) {
  static struct quantities_vtable vtable = {
    &sqpipe_quantities_get_quantity,
    &sqpipe_quantities_update,
    &sqpipe_quantities_destroy
  };
  q->super.vtable = &vtable;


  q->a = NULL;
  q->d = NULL;
  q->w = NULL;
  q->h = NULL;
  q->u = NULL;
  q->z = NULL;  
  q->b = NULL;
  q->t = NULL;
}

struct sqpipe_quantities* sqpipe_quantities_new() {
  struct sqpipe_quantities *p = malloc(sizeof(struct sqpipe_quantities));

  return p;  
}


struct quantities* sqpipe_quantities_destroy(struct quantities *qs) {
  free(qs);
  qs = NULL;
  
  return qs;
}
/**********************************************************/


/**********************************************************
Pipe domain class
This is a virtual class and should never be created.
The inheriting class needs to implement:
  init, new, compute_fluxes
**********************************************************/
struct sqpipe_domain {
  struct domain super;
  long *state;
};

// Construction
void sqpipe_domain_init(struct sqpipe_domain *d);
struct sqpipe_domain* sqpipe_domain_new();

// Destruction
struct domain* sqpipe_domain_destroy(struct domain *d);

// Methods
struct cell* sqpipe_domain_get_cell(struct domain *D, int k, struct cell *c);

// Implementation
struct cell* sqpipe_domain_get_cell(struct domain *D, int k, struct cell *c) {
  struct sqpipe_cell *sc = (struct sqpipe_cell *) c;
  struct sqpipe_domain *sD = (struct sqpipe_domain *) D;
  int i, ki, n;
  long state_in, state_out, state;
  
  // First get a generic cell
  c = domain_get_cell_generic(D, k, c);

  // Set state
  // The state should be set to state_in for qin and state_out for qout
  // Temporarily, I set the state of the edge, aribtrarily choosing
  // pressurised for transition points
  // This will be fixed when the flux is computed at the edge rather than
  // at the quantity
  sc->state = sD->state[k];
  for (i = 0; i<c->num_edges; i++) {
	// The inner state is that of this cell
	state_in = sc->state;
	// The outer state is that of the corresponding neighbour cell
	ki = k*(c->num_edges) + i;
	n = D->neighbours[ki];
	state_out = sD->state[n];
	
	// If both in and out state are the same, choose that state
	// Otherwise, choose pressurised
	if (state_in == state_out) {
		state = state_in;
	} else {
		state = 1;
	}
    ((struct sqpipe_quantity *)(c->edges[i]->qin))->state = state;	
    ((struct sqpipe_quantity *)(c->edges[i]->qout))->state = state;
  }
  
  return c;
}


void sqpipe_domain_init(struct sqpipe_domain *d) {
  static struct domain_vtable vtable = {
    NULL,
    NULL,
    &domain_compute_fluxes_cells_generic,
    &sqpipe_domain_get_cell,
    &sqpipe_domain_destroy
  };
  d->super.vtable = &vtable;

  ((struct domain *)d)->number_of_equations = 2;

  ((struct domain *)d)->vertex_values = NULL;
  ((struct domain *)d)->boundary_values = NULL;
  ((struct domain *)d)->explicit_update = NULL;
  ((struct domain *)d)->neighbours = NULL;
  ((struct domain *)d)->neighbour_vertices = NULL;
  ((struct domain *)d)->normals = NULL;
  ((struct domain *)d)->areas = NULL;
  ((struct domain *)d)->number_of_elements = 0;  

  ((struct domain *)d)->timestep = 0.0;
  ((struct domain *)d)->max_speed_array = NULL;

  d->state = NULL;
  // Can't set params to NULL here since sqpipe_domain_new() has allocated 
  // memory
  // If a domain is created without using domain_new(), then params will be an
  // unitialised pointer.

  // Can't allocate memory here since init should not allocate memory
  // The moral: Use domain_new() to create domains!  
}

struct sqpipe_domain* sqpipe_domain_new() {
  struct sqpipe_domain *p = malloc(sizeof(struct sqpipe_domain));

  ((struct domain *)p)->params = (struct params *) sqpipe_params_new();

  return p;
}


struct domain* sqpipe_domain_destroy(struct domain *D) {
  D->vertex_values = quantities_destroy(D->vertex_values);
  D->boundary_values = quantities_destroy(D->boundary_values);
  D->explicit_update = quantities_destroy(D->explicit_update);
  D->params = params_destroy(D->params);
  free(D);
  D = NULL;

  return D;
}
/**********************************************************/


