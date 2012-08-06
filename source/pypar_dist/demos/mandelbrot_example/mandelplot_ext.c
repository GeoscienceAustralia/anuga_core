// Python - C extension for fast computation of the Mandelbrot plot routine
//
// To compile (Python2.3):
//  gcc -c mandelplot_ext.c -I/usr/include/python2.3 -o mandelplot_ext.o -Wall -O3 
//  gcc -shared mandelplot_ext.o  -o mandelplot_ext.so	
//
// kmax = 256
// c = complex(0.5, 0.5)
// k = calculate_point(c, kmax)
//
// Ole Nielsen, ANU 2004
	
	
#include "Python.h"
#include "numpy/arrayobject.h"

// Computational function
int  _normalise_and_convert(double* A, 
			    int M, 
			    int N,
			    PyObject* L, 
			    int kmax, 
			    int rgbmax, 
			    double exponent) {
			     
  /*			     
  for i in range(A.shape[0]):
    for j in range(A.shape[1]):    
      
      c = A[i,j]/kmax
      if c == 1: c = 0       #Map convergent point (kmax) to black (0)
      c = c**exponent        #Morph slightly
      c = int(c * rgbmax)    #Normalise to 256 levels per channel
 		
      red   = c / 256 / 256
      green = (c / 256) % 256
      blue  = c % 256
      
      L.append( (red, green, blue) ) */
	    

  
  PyObject *tuple;
  int i, j, iN, red, green, blue, ci;
  double c;
  
  for (i=0; i<M; i++) {
    iN = i*N;
    
    for (j=0; j<N; j++) {  
    
      c = A[iN  + j]/kmax;     //Normalise to unit interval	    	    
      c = pow(c, exponent);        //Morph slightly
      ci = (int) (c*rgbmax);       //Normalise to rgbmax levels per channel
      if (ci == rgbmax) ci = 0.0;  //Map convergent point (kmax) to black (0)
      

      //Convert to RGB
      red   = ci / 256 / 256;
      green = (ci / 256) % 256;
      blue  = ci % 256;			     			     
      
      tuple = Py_BuildValue("iii", red, green, blue);
      PyList_Append(L, tuple);
    }
  }
	
  return 0;
}


// Interface to Python
PyObject *normalise_and_convert(PyObject *self, PyObject *args) {
  
  //Called 
  //normalise_and_convert(A, L, kmax, rgbmax, exponent)    
  
  int kmax, rgbmax, M, N;
  double exponent;
  PyObject *L;
  PyArrayObject *A;

  // Convert Python arguments to C  
  if (!PyArg_ParseTuple(args, "OOiid", &A, &L, &kmax, &rgbmax, &exponent))
    return NULL;
    
  M = A -> dimensions[0];
  N = A -> dimensions[1];  
  

  // Call underlying routine
  _normalise_and_convert((double *) A -> data, M, N, L, 
			 kmax, rgbmax, exponent); 
  
  // Return None
  return Py_BuildValue("");
}


// Method table for python module
static struct PyMethodDef MethodTable[] = {
  {"normalise_and_convert", normalise_and_convert, METH_VARARGS},  
  {NULL, NULL}
};


// Module initialisation   
void initmandelplot_ext(void){
  Py_InitModule("mandelplot_ext", MethodTable);
  import_array();     //Necessary for handling of numpy structures
}






