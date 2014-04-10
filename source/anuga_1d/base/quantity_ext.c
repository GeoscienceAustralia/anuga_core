#include "Python.h"
#include "numpy/arrayobject.h"
#include "math.h"
#include <stdio.h>
const double pi = 3.14159265358979;


// Shared code snippets
#include "util_ext.h"








//=========================================================================
// Python Glue
//=========================================================================

PyObject *limit_minmod_ext(PyObject *self, PyObject *args) {
  
   PyObject
           *domain,
           *quantity;

    PyArrayObject 
	*qco,
	*qvo,
	*xco,
	*xvo;

    double *qc,
            *qv,
            *xc,
            *xv;
    
  double a, b;
  double phi, dx0, dx1;
  int N, k, k2;

    
  // Convert Python arguments to C
  if (!PyArg_ParseTuple(args, "O", &quantity)) {
      PyErr_SetString(PyExc_RuntimeError, "quantity_ext.c: limit_minmod_ext could not parse input");
      return NULL;
  }
  

  domain = PyObject_GetAttrString(quantity, "domain");

  //printf("B = %p\n",(void*)domain);
  if (!domain) {
    printf("quantity_ext.c: Could not obtain python object");
    fflush(stdout);
    PyErr_SetString(PyExc_RuntimeError, "quantity_ext.c: Could not obtain python object domain");
    return NULL;
  }
  

  
  N  = get_python_integer(quantity,"N");
 
  qco = get_consecutive_array(quantity, "centroid_values");
  qvo = get_consecutive_array(quantity, "vertex_values");
  xco = get_consecutive_array(domain, "centroids");
  xvo = get_consecutive_array(domain, "vertices");
    
  qc = (double *) qco -> data;
  qv = (double *) qvo -> data;
  xc = (double *) xco -> data;
  xv = (double *) xvo -> data;


  for (k=0; k<N; k++) {
      k2 = 2*k;
      if (k == 0) {
        phi = (qc[1]-qc[0])/(xc[1] - xc[0]);
      }
      else if (k==N-1) {
          phi = (qc[N-1] - qc[N-2])/(xc[N-1] - xc[N-2]);
      }
      else {
          a  = (qc[k]-qc[k-1])/(xc[k]-xc[k-1]);
          b  = (qc[k+1]-qc[k])/(xc[k+1]-xc[k]);
          //c  = (qc[K+1]-qc[k-1])/(xc[k+1]-xc[k-1]);

          phi = 0.0;
          if ((fabs(a) < fabs(b)) & (a*b >= 0.0 )) {
              phi = a;
          }
          if ((fabs(b) < fabs(a)) & (a*b >= 0.0 )) {
              phi = b;
          }

      }
      
      dx0 = xv[k2] - xc[k];
      dx1 = xv[k2+1] - xc[k];


      qv[k2]   = qc[k] + phi*dx0;
      qv[k2+1] = qc[k] + phi*dx1;
  }


  Py_DECREF(qco);
  Py_DECREF(qvo);
  Py_DECREF(xco);
  Py_DECREF(xvo);

  // Return updated flux timestep
  return Py_BuildValue("");
}


//====================================================================
PyObject *limit_minmod_kurganov_ext(PyObject *self, PyObject *args) {

   PyObject
           *domain,
           *quantity;

    PyArrayObject
	*qco,
	*qvo,
	*xco,
	*xvo;

    double *qc,
            *qv,
            *xc,
            *xv;

  double a, b, c;
  double phi, dx0, dx1;
  double theta;
  int N, k, k2;


  // Convert Python arguments to C
  if (!PyArg_ParseTuple(args, "O", &quantity)) {
      PyErr_SetString(PyExc_RuntimeError, "quantity_ext.c: limit_minmod_kurganov_ext could not parse input");
      return NULL;
  }


  domain = PyObject_GetAttrString(quantity, "domain");

  //printf("B = %p\n",(void*)domain);
  if (!domain) {
    printf("quantity_ext.c: Could not obtain python object");
    fflush(stdout);
    PyErr_SetString(PyExc_RuntimeError, "quantity_ext.c: Could not obtain python object domain");
    return NULL;
  }



  N  = get_python_integer(quantity,"N");
  theta = get_python_double(quantity,"beta");


  //printf("beta = %f",theta);
  //fflush(stdout);

  qco = get_consecutive_array(quantity, "centroid_values");
  qvo = get_consecutive_array(quantity, "vertex_values");
  xco = get_consecutive_array(domain, "centroids");
  xvo = get_consecutive_array(domain, "vertices");

  qc = (double *) qco -> data;
  qv = (double *) qvo -> data;
  xc = (double *) xco -> data;
  xv = (double *) xvo -> data;


  for (k=0; k<N; k++) {
      k2 = 2*k;
      if (k == 0) {
        phi = (qc[1]-qc[0])/(xc[1] - xc[0]);
      }
      else if (k==N-1) {
          phi = (qc[N-1] - qc[N-2])/(xc[N-1] - xc[N-2]);
      }
      else {
          a  = (qc[k]-qc[k-1])/(xc[k]-xc[k-1]);
          b  = (qc[k+1]-qc[k])/(xc[k+1]-xc[k]);
          c  = (qc[k+1]-qc[k-1])/(xc[k+1]-xc[k-1]);

          phi = 0.0;
          if ((sign(a)*sign(b) > 0.0) & (sign(a)*sign(c) >= 0.0 )) {
              phi = sign(a)*min(theta*min(fabs(a),fabs(b)),fabs(c));
          }


      }

      dx0 = xv[k2] - xc[k];
      dx1 = xv[k2+1] - xc[k];


      qv[k2]   = qc[k] + phi*dx0;
      qv[k2+1] = qc[k] + phi*dx1;
  }


  Py_DECREF(qco);
  Py_DECREF(qvo);
  Py_DECREF(xco);
  Py_DECREF(xvo);

  // Return updated flux timestep
  return Py_BuildValue("");
}

 
//====================================================================
PyObject *limit_vanleer_ext(PyObject *self, PyObject *args) {

   PyObject
           *domain,
           *quantity;

    PyArrayObject
	*qco,
	*qvo,
	*xco,
	*xvo;

    double *qc,
            *qv,
            *xc,
            *xv;

  double a, b;
  double phi, dx0, dx1;
  double theta;
  int N, k, k2;


  // Convert Python arguments to C
  if (!PyArg_ParseTuple(args, "O", &quantity)) {
      PyErr_SetString(PyExc_RuntimeError, "quantity_ext.c: limit_vanleer_ext could not parse input");
      return NULL;
  }


  domain = PyObject_GetAttrString(quantity, "domain");

  //printf("B = %p\n",(void*)domain);
  if (!domain) {
    printf("quantity_ext.c: Could not obtain python object");
    fflush(stdout);
    PyErr_SetString(PyExc_RuntimeError, "quantity_ext.c: Could not obtain python object domain");
    return NULL;
  }



  N  = get_python_integer(quantity,"N");
  theta = get_python_double(quantity,"beta");


  //printf("beta = %f",theta);
  //fflush(stdout);

  qco = get_consecutive_array(quantity, "centroid_values");
  qvo = get_consecutive_array(quantity, "vertex_values");
  xco = get_consecutive_array(domain, "centroids");
  xvo = get_consecutive_array(domain, "vertices");

  qc = (double *) qco -> data;
  qv = (double *) qvo -> data;
  xc = (double *) xco -> data;
  xv = (double *) xvo -> data;


  for (k=0; k<N; k++) {
      k2 = 2*k;
      if (k == 0) {
        phi = (qc[1]-qc[0])/(xc[1] - xc[0]);
      }
      else if (k==N-1) {
          phi = (qc[N-1] - qc[N-2])/(xc[N-1] - xc[N-2]);
      }
      else {
          a  = (qc[k]-qc[k-1])/(xc[k]-xc[k-1]);
          b  = (qc[k+1]-qc[k])/(xc[k+1]-xc[k]);
          //c  = (qc[k+1]-qc[k-1])/(xc[k+1]-xc[k-1]);


          phi = 0.0;
          if ((fabs(a)+fabs(b)) > 1.0e-12) {
              phi = (a*fabs(b)+fabs(a)*b)/(fabs(a)+fabs(b));
          }
          //printf("phi = %f",phi);
          //fflush(stdout);

      }

      dx0 = xv[k2] - xc[k];
      dx1 = xv[k2+1] - xc[k];


      qv[k2]   = qc[k] + phi*dx0;
      qv[k2+1] = qc[k] + phi*dx1;
  }


  Py_DECREF(qco);
  Py_DECREF(qvo);
  Py_DECREF(xco);
  Py_DECREF(xvo);

  // Return updated flux timestep
  return Py_BuildValue("");
}


//====================================================================
PyObject *limit_vanalbada_ext(PyObject *self, PyObject *args) {

   PyObject
           *domain,
           *quantity;

    PyArrayObject
	*qco,
	*qvo,
	*xco,
	*xvo;

    double *qc,
            *qv,
            *xc,
            *xv;

  double a, b;
  //double c;
  double phi, dx0, dx1;
  double theta;
  int N, k, k2;


  // Convert Python arguments to C
  if (!PyArg_ParseTuple(args, "O", &quantity)) {
      PyErr_SetString(PyExc_RuntimeError, "quantity_ext.c: limit_vanalbada_ext could not parse input");
      return NULL;
  }


  domain = PyObject_GetAttrString(quantity, "domain");

  //printf("B = %p\n",(void*)domain);
  if (!domain) {
    printf("quantity_ext.c: Could not obtain python object");
    fflush(stdout);
    PyErr_SetString(PyExc_RuntimeError, "quantity_ext.c: Could not obtain python object domain");
    return NULL;
  }



  N  = get_python_integer(quantity,"N");
  theta = get_python_double(quantity,"beta");


  //printf("beta = %f",theta);
  //fflush(stdout);

  qco = get_consecutive_array(quantity, "centroid_values");
  qvo = get_consecutive_array(quantity, "vertex_values");
  xco = get_consecutive_array(domain, "centroids");
  xvo = get_consecutive_array(domain, "vertices");

  qc = (double *) qco -> data;
  qv = (double *) qvo -> data;
  xc = (double *) xco -> data;
  xv = (double *) xvo -> data;


  for (k=0; k<N; k++) {
      k2 = 2*k;
      if (k == 0) {
        phi = (qc[1]-qc[0])/(xc[1] - xc[0]);
      }
      else if (k==N-1) {
          phi = (qc[N-1] - qc[N-2])/(xc[N-1] - xc[N-2]);
      }
      else {
          a  = (qc[k]-qc[k-1])/(xc[k]-xc[k-1]);
          b  = (qc[k+1]-qc[k])/(xc[k+1]-xc[k]);
          //c  = (qc[k+1]-qc[k-1])/(xc[k+1]-xc[k-1]);


          phi = 0.0;
          if (a*a + b*b >= 1.0e-32) {
              phi = (a*a*b+a*b*b)/(a*a+b*b);
          }


      }

      dx0 = xv[k2] - xc[k];
      dx1 = xv[k2+1] - xc[k];


      qv[k2]   = qc[k] + phi*dx0;
      qv[k2+1] = qc[k] + phi*dx1;
  }


  Py_DECREF(qco);
  Py_DECREF(qvo);
  Py_DECREF(xco);
  Py_DECREF(xvo);

  // Return updated flux timestep
  return Py_BuildValue("");
}



//-------------------------------
// Method table for python module
//-------------------------------

static struct PyMethodDef MethodTable[] = {
  {"limit_minmod_ext",          limit_minmod_ext, METH_VARARGS, "Print out"},
  {"limit_minmod_kurganov_ext", limit_minmod_kurganov_ext, METH_VARARGS, "Print out"},
  {"limit_vanleer_ext",         limit_vanleer_ext, METH_VARARGS, "Print out"},
  {"limit_vanalbada_ext",       limit_vanalbada_ext, METH_VARARGS, "Print out"},
  {NULL, NULL}
};

// Module initialisation
void initquantity_ext(void){
  Py_InitModule("quantity_ext", MethodTable);
  import_array();
}
