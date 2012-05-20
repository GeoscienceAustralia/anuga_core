// Python - C extension for sparse module.
//
// To compile (Python2.3):
//  gcc -c util_ext.c -I/usr/include/python2.3 -o util_ext.o -Wall -O 
//  gcc -shared util_ext.o  -o util_ext.so	
//
// See the module sparse.py
//
//
// Steve Roberts, ANU 2004
	
#include "Python.h"
#include "numpy/arrayobject.h"
#include "math.h"
#include "stdio.h"
#include "numpy_shim.h"

//Matrix-vector routine
int _csr_mv(int M,
	    double* data, 
	    long* colind,
	    long* row_ptr,
	    double* x,
	    double* y) {
  		
  long i, j, ckey;

  for (i=0; i<M; i++ ) 
    for (ckey=row_ptr[i]; ckey<row_ptr[i+1]; ckey++) {
      j = colind[ckey];
      y[i] += data[ckey]*x[j];
    }              
  
  return 0;
}            

//Matrix-matrix routine
int _csr_mm(int M,
	    int columns, 
	    double* data, 
	    long* colind,
	    long* row_ptr,
	    double* x,
	    double* y) {
  		
  long i, j, ckey, c, rowind_i, rowind_j;

  for (i=0; i<M; i++ ) {
    rowind_i = i*columns;
    
    for (ckey=row_ptr[i]; ckey<row_ptr[i+1]; ckey++) {
      j = colind[ckey];
      rowind_j = j*columns;
          
      for (c=0; c<columns; c++) {
        y[rowind_i+c] += data[ckey]*x[rowind_j+c];
      }             
    } 
  }
  
  return 0;
}            

		     
  
/////////////////////////////////////////////////
// Gateways to Python
PyObject *csr_mv(PyObject *self, PyObject *args) {
  
  PyObject *csr_sparse, // input sparse matrix
    *xin, *R;           // output wrapped vector
  
  PyArrayObject 
    *data,            //Non Zeros Data array
    *colind,          //Column indices array
    *row_ptr,         //Row pointers array
    *x,               //Input vector array
    *y;               //Return vector array

  
  int dimensions[2], M, err, columns, rows;
  
  // Convert Python arguments to C  
  if (!PyArg_ParseTuple(args, "OO", &csr_sparse, &xin)) {
    PyErr_SetString(PyExc_RuntimeError, "Csr_mv could not parse input");  
    return NULL;
  }

  x = (PyArrayObject*) PyArray_ContiguousFromObject(xin,PyArray_DOUBLE,1,2);
  if (!x) {
    PyErr_SetString(PyExc_RuntimeError, 
		    "Input array could not be read in csr_mv");    
    return NULL;
  }

/*   printf("x.nd = %i\n",x->nd); */
/*   printf("x.descr->type_num = %i %i\n",x->descr->type_num,PyArray_LONG); */
/*   printf("x.dimensions[0] = %i\n",x->dimensions[0]); */
/*   printf("x.data[0] = %g\n",((double*) x->data)[0]); */
/*   printf("x.data[1] = %g\n",((double*) x->data)[1]); */
/*   printf("x.data[2] = %g\n",((double*) x->data)[2]); */
/*   printf("x.data[3] = %g\n",((double*) x->data)[3]); */
/*   printf("x.data[4] = %g\n",((double*) x->data)[4]); */
/*   printf("x.data[5] = %g\n",((double*) x->data)[5]); */

 


  data = (PyArrayObject*) 
    PyObject_GetAttrString(csr_sparse, "data");     
  if (!data) {
    PyErr_SetString(PyExc_RuntimeError, 
		    "Data array could not be allocated in csr_mv");      
    return NULL;
  }  

  colind = (PyArrayObject*)
    PyObject_GetAttrString(csr_sparse, "colind"); 
  if (!colind) {
    PyErr_SetString(PyExc_RuntimeError, 
		    "Column index array could not be allocated in csr_mv");      
    return NULL;
  }  

  row_ptr = (PyArrayObject*) 
    PyObject_GetAttrString(csr_sparse, "row_ptr");   
  if (!row_ptr) {
    PyErr_SetString(PyExc_RuntimeError, 
		    "Row pointer array could not be allocated in csr_mv"); 
  }
  
  M = (row_ptr -> dimensions[0])-1;
    
  if (x -> nd == 1) {
    // Multiplicant is a vector
  
    //Allocate space for return vectors y (don't DECREF) 
    dimensions[0] = M;
    y = (PyArrayObject *) anuga_FromDims(1, dimensions, PyArray_DOUBLE);
  
    err = _csr_mv(M,
		  (double*) data -> data, 
		  (long*)   colind -> data,
		  (long*)   row_ptr -> data,
		  (double*) x -> data,
		  (double*) y -> data); 

			   
    if (err != 0) {
      PyErr_SetString(PyExc_RuntimeError, "Matrix vector mult could not be calculated");
      return NULL;
    }
  } else if(x -> nd == 2) {
  

    rows = x -> dimensions[0];     //Number of rows in x        
    columns = x -> dimensions[1];  //Number of columns in x        
    
    //Allocate space for return matrix y (don't DECREF) 
    dimensions[0] = M;                   //Number of rows in sparse matrix  
    dimensions[1] = columns;
    y = (PyArrayObject *) anuga_FromDims(2, dimensions, PyArray_DOUBLE);
    
    err = _csr_mm(M, columns,
		  (double*) data -> data, 
		  (long*)   colind -> data,
		  (long*)   row_ptr -> data,
		  (double*) x -> data,
		  (double*) y -> data); 
    
  } else {
    PyErr_SetString(PyExc_RuntimeError, 
		    "Allowed dimensions in sparse_ext.c restricted to 1 or 2");
    return NULL;  
  }
  
  		     
  //Release		     
  Py_DECREF(data);    
  Py_DECREF(colind);    
  Py_DECREF(row_ptr); 
  Py_DECREF(x);           
         
  //Build result, release and return
  R = Py_BuildValue("O", PyArray_Return(y)); 
  Py_DECREF(y);        
  return R;
}




// Method table for python module
static struct PyMethodDef MethodTable[] = {
  {"csr_mv", csr_mv, METH_VARARGS, "Print out"},    
  {NULL, NULL, 0, NULL}   /* sentinel */
};

// Module initialisation   
void initsparse_ext(void){
  Py_InitModule("sparse_ext", MethodTable);
  
  import_array();     //Necessary for handling of NumPY structures  
}

