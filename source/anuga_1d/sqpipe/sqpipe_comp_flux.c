#include "Python.h"
#include "numpy/arrayobject.h"
#include "math.h"
#include <stdio.h>

const double pi = 3.14159265358979;

// Shared code snippets
// Generic routines
#include "util_ext.h"
// Solver
#include "solver.h"
// Python-solver interface
#include "python_solver.h"
// square pipe solver implementation
#include "sqpipe_default.h"
// sqpipe-python glue
#include "sqpipe_python.h"

//=========================================================================
// Python Glue
//=========================================================================

PyObject *compute_fluxes_ext(PyObject *self, PyObject *args) {
  struct domain *D = (struct domain *) sqpipe_default_domain_new();
  double timestep;
  PyObject *domain;

  // Convert Python arguments to C
  if (!PyArg_ParseTuple(args, "dO", &timestep, &domain)) {
    PyErr_SetString(PyExc_RuntimeError, "comp_flux_ext.c: compute_fluxes_ext could not parse input");
    return NULL;
  }

  // Generic sqpipe domain
  D = (struct domain *) sqpipe_domain_python_get((struct sqpipe_domain *)D, domain, timestep);

  // We need the state 
  //state = get_consecutive_array(domain, "state");
  //((struct sqpipe_default_domain *)D)->state = (long *) state->data;

  // Compute the fluxes and return the new timestep
  timestep = domain_compute_fluxes(D);
  domain_destroy(D);
  
  return Py_BuildValue("d", timestep);
}

//-------------------------------
// Method table for python module
//-------------------------------

static struct PyMethodDef MethodTable[] = {
  {"compute_fluxes_ext", compute_fluxes_ext, METH_VARARGS, "Print out"},
  {NULL, NULL}
};

// Module initialisation
void initsqpipe_comp_flux(void){
  Py_InitModule("sqpipe_comp_flux", MethodTable);
  import_array();
}
