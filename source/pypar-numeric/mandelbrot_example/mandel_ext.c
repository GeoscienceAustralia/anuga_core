/* Python - C extension for fast computation of the Mandelbrot iteration
//
// To compile:
//  python compile.py mandel_ext.c
//
// Example of how to use within Python:
//
// from mandel_ext import calculate_point
//
// kmax = 256
// c = complex(0.5, 0.5)
// k = calculate_point(c, kmax)
//
// Ole Nielsen, SUT 2003 */
	
#include "Python.h"

/* Computational function */
int _calculate_point(Py_complex c, int kmax) {
  int count;
  double temp, lengthsq = 0.0;  
  Py_complex z;
  
  z.real = 0.0; z.imag = 0.0; 
  count = 0;

  while (count < kmax && lengthsq <= 4) {
       temp = z.real * z.real - z.imag * z.imag + c.real;
       z.imag = 2.0 * z.real * z.imag + c.imag;
       z.real = temp;
       lengthsq = z.real * z.real + z.imag * z.imag;
       count++;    
  }

  return count;
}


/* Interface to Python */
PyObject *calculate_point(PyObject *self, PyObject *args) {
  PyComplexObject *C;  /* Python Complex object */
  Py_complex c;        /* C Complex structure */
  int kmax, count;

  /* Convert Python arguments to C  */
  if (!PyArg_ParseTuple(args, "Oi", &C, &kmax))
    return NULL;
  c = PyComplex_AsCComplex((PyObject*) C);   
  
  /* Call underlying routine */
  count = _calculate_point(c, kmax); 

  /* Return results as a Python object */
  return Py_BuildValue("i", count);
}


/* Method table for python module */
static struct PyMethodDef MethodTable[] = {
  {"calculate_point", calculate_point, METH_VARARGS},  
  {NULL, NULL}
};


/* Module initialisation  */
void initmandel_ext(void){
  Py_InitModule("mandel_ext", MethodTable);
}

