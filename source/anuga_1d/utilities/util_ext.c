// Python - C extension for finite_volumes util module.
//
// To compile (Python2.3):
//  gcc -c util_ext.c -I/usr/include/python2.3 -o util_ext.o -Wall -O
//  gcc -shared util_ext.o  -o util_ext.so
//
// See the module util.py
//
//
// Ole Nielsen, GA 2004
//
//NOTE: On 64 bit systems use long* instead of int* for numeric arrays
//this will also work on 32 bit systems

#include <float.h>

#include "Python.h"
#include "numpy/arrayobject.h"
#include "math.h"

//Shared snippets
#include "util_ext.h"




PyObject *gradient(PyObject *self, PyObject *args) {
  //
  // a,b = gradient(x0, y0, x1, y1, x2, y2, q0, q1, q2)
  //

  double x0, y0, x1, y1, x2, y2, q0, q1, q2, a, b;
  PyObject *result;

  // Convert Python arguments to C
  if (!PyArg_ParseTuple(args, "ddddddddd", &x0, &y0, &x1, &y1, &x2, &y2,
			&q0, &q1, &q2)) {
    
    PyErr_SetString(PyExc_RuntimeError, 
		    "gradient could not parse input");    
    return NULL;
  }


  // Call underlying routine
  _gradient(x0, y0, x1, y1, x2, y2,
	    q0, q1, q2, &a, &b);

  // Return values a and b
  result = Py_BuildValue("dd", a, b);
  return result;
}

PyObject *gradient2(PyObject *self, PyObject *args) {
  //
  // a,b = gradient2(x0, y0, x1, y1, q0, q1)
  //

  double x0, y0, x1, y1, q0, q1, a, b;
  PyObject *result;

  // Convert Python arguments to C
  if (!PyArg_ParseTuple(args, "dddddd", &x0, &y0, &x1, &y1, &q0, &q1)) {
    PyErr_SetString(PyExc_RuntimeError, 
		    "gradient2 could not parse input");      
    return NULL;
  }


  // Call underlying routine
  _gradient2(x0, y0, x1, y1, q0, q1, &a, &b);

  // Return values a and b
  result = Py_BuildValue("dd", a, b);
  return result;
}

PyObject * double_precision(PyObject * self, PyObject * args){
  // Get the precision of the double datatype on this system.
  return Py_BuildValue("i", DBL_DIG);
}

// Method table for python module
static struct PyMethodDef MethodTable[] = {
  /* The cast of the function is necessary since PyCFunction values
   * only take two PyObject* parameters, and rotate() takes
   * three.
   */

  //{"rotate", (PyCFunction)rotate, METH_VARARGS | METH_KEYWORDS, "Print out"},
  {"gradient", gradient, METH_VARARGS, "Print out"},
  {"gradient2", gradient2, METH_VARARGS, "Print out"},  
  {"double_precision", double_precision, METH_VARARGS, "Precision of this machine\'s \'double\' type"},
  {NULL, NULL, 0, NULL}   /* sentinel */
};



// Module initialisation
void initutil_ext(void){
  Py_InitModule("util_ext", MethodTable);

  import_array();     //Necessary for handling of NumPY structures
}




