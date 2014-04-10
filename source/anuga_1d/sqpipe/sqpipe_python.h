struct sqpipe_domain* sqpipe_domain_python_get(struct sqpipe_domain *sD, PyObject *domain, double timestep);
struct quantities* sqpipe_quantities_python_get(struct quantities *qs, PyObject* quantities, char *name);

struct sqpipe_domain* sqpipe_domain_python_get(struct sqpipe_domain *sD, PyObject *domain, double timestep) {
  PyObject *quantities;
  PyArrayObject *state;
  struct domain *D = (struct domain *)sD;
  struct sqpipe_params *sp = (struct sqpipe_params *) D->params;
  
  // Get a generic domain
  D = get_python_domain(D, domain, timestep);

  // Get the specific quantities for square pipes
  quantities = get_python_object(domain, "quantities");
  D->vertex_values = sqpipe_quantities_python_get(D->vertex_values, quantities, "vertex_values");
  D->boundary_values = sqpipe_quantities_python_get(D->boundary_values, quantities, "boundary_values");
  D->explicit_update = sqpipe_quantities_python_get(D->explicit_update, quantities, "explicit_update");

  // Get cell level quantities
  state = get_consecutive_array(domain, "state");
  sD->state = (long *) state->data;

  // Get the specific parameters for square pipes
  sp->g = get_python_double(domain,"g");
  sp->h0 = get_python_double(domain,"h0");
  sp->bulk_modulus = get_python_double(domain, "bulk_modulus");

  return sD;
}

struct quantities* sqpipe_quantities_python_get(struct quantities *qs, PyObject* quantities, char *name) {
  struct sqpipe_quantities *ps = (struct sqpipe_quantities *) qs;
  
  ps->a = get_python_array_data_from_dict(quantities, "area", name);
  ps->d = get_python_array_data_from_dict(quantities, "discharge", name);
  ps->w = get_python_array_data_from_dict(quantities, "stage", name);
  ps->h = get_python_array_data_from_dict(quantities, "height", name);
  ps->u = get_python_array_data_from_dict(quantities, "velocity", name);  
  ps->z = get_python_array_data_from_dict(quantities, "elevation", name);
  ps->b = get_python_array_data_from_dict(quantities, "width", name);
  ps->t = get_python_array_data_from_dict(quantities, "top", name);

  return qs;
}
